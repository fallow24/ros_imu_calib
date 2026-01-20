#include "imu_tk/multicalibration.h"
#include "imu_tk/filters.h"
#include "imu_tk/integration.h"
#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include <algorithm>

using namespace imu_tk;
using namespace Eigen;

// Gyro-only rotation alignment residual
// Optimizes rotation R (3 params: axis-angle on SO(3) manifold)
// T = R * M * S where M and S are fixed from intrinsic calibration
// Gyro provides full 3-DOF rotation observability
// The resulting R is then applied to both gyro and acc (same chip)
template <typename _T1>
struct GyroRotationAlignResidual
{
    GyroRotationAlignResidual(
        const TriadData_<_T1> &ref_gyr,        // Reference calibrated gyro
        const TriadData_<_T1> &other_gyr_raw,  // Other raw gyro
        const Eigen::Matrix<_T1, 3, 3> &M_gyr, // Fixed gyro misalignment
        const Eigen::Matrix<_T1, 3, 3> &S_gyr, // Fixed gyro scale
        const Eigen::Matrix<_T1, 3, 1> &B_gyr) // Fixed gyro bias
        : ref_gyr_(ref_gyr), other_gyr_raw_(other_gyr_raw),
          M_gyr_(M_gyr), S_gyr_(S_gyr), B_gyr_(B_gyr)
    {
    }

    template <typename _T2>
    bool operator()(const _T2 *const axis_angle, _T2 *residuals) const
    {
        // Convert axis-angle to rotation matrix using Ceres helper
        _T2 R_data[9];
        ceres::AngleAxisToRotationMatrix(axis_angle, R_data);
        Eigen::Map<const Eigen::Matrix<_T2, 3, 3, Eigen::ColMajor>> R(R_data);

        // T_gyr = R * M_gyr * S_gyr
        Eigen::Matrix<_T2, 3, 3> T_gyr = R * M_gyr_.template cast<_T2>() * S_gyr_.template cast<_T2>();
        Eigen::Matrix<_T2, 3, 1> v_gyr_raw = other_gyr_raw_.data().template cast<_T2>() - B_gyr_.template cast<_T2>();
        Eigen::Matrix<_T2, 3, 1> v_gyr = T_gyr * v_gyr_raw;
        Eigen::Matrix<_T2, 3, 1> diff = v_gyr - ref_gyr_.data().template cast<_T2>();

        residuals[0] = diff(0);
        residuals[1] = diff(1);
        residuals[2] = diff(2);

        return true;
    }

    const TriadData_<_T1> ref_gyr_, other_gyr_raw_;
    const Eigen::Matrix<_T1, 3, 3> M_gyr_, S_gyr_;
    const Eigen::Matrix<_T1, 3, 1> B_gyr_;
};

template <typename _T>
MultiPosRedundantCalibration_<_T>::MultiPosRedundantCalibration_()
    : acc_calibs_(), gyro_calibs_(), extrinsic_vectors_(),
      g_mag_(_T(9.81)),
      init_interval_duration_(_T(30.0)),
      interval_n_samples_(100),
      acc_use_means_(false),
      gyro_dt_(_T(-1.0)),
      lowpass_frequency_cut_(2),
      sigma_(_T(0)),
      win_size_(5),
      verbose_output_(false),
      plot_enabled_(true)
{
}

template <typename _T>
bool MultiPosRedundantCalibration_<_T>::calibrateAcc(
    const std::vector<std::vector<TriadData_<_T>>> &acc_samples_all,
    size_t reference_index,
    int min_recorded_intervals,
    int variance_win_size,
    int th_var_min,
    int th_var_max)
{
    std::cout << "\033[1;32m------------------------------------------------------------------------\033[0m" << std::endl;
    std::cout << "\033[1;32m------------ NOW PERFORMING ACCELEROMETER CALIBRATION ------------------\033[0m" << std::endl;
    std::cout << "\033[1;32m------------------------------------------------------------------------\033[0m" << std::endl;
    if (reference_index >= acc_samples_all.size())
    {
        return false;
    }
    acc_calibs_.clear();
    acc_calibs_.resize(acc_samples_all.size());

    // Configure per-topic calibrators and use cache
    per_topic_calibs_.resize(acc_samples_all.size());
    for (size_t i = 0; i < acc_samples_all.size(); ++i)
    {
        auto &calib_single = per_topic_calibs_[i];
        calib_single.setGravityMagnitude(g_mag_);
        calib_single.setInitStaticIntervalDuration(init_interval_duration_);
        calib_single.setIntarvalsNumSamples(interval_n_samples_);
        calib_single.enableAccUseMeans(acc_use_means_);
        calib_single.enableVerboseOutput(verbose_output_);
        calib_single.enablePlotting(plot_enabled_);

        if (!calib_single.isAccCached())
        {
            if (!calib_single.calibrateAcc(acc_samples_all[i],
                                           min_recorded_intervals,
                                           variance_win_size,
                                           th_var_min,
                                           th_var_max))
                return false;
        }
        acc_calibs_[i] = calib_single.getAccCalib();
    }

    // NOTE: Accelerometer alignment is deferred to calibrateAccGyroExtrinsic
    // where we have gyro data and extrinsic vectors for motion compensation.
    // At this stage we only have intrinsic calibration per sensor.
    return true;
}

template <typename _T>
bool MultiPosRedundantCalibration_<_T>::calibrateAccGyro(
    const std::vector<std::vector<TriadData_<_T>>> &acc_samples_all,
    const std::vector<std::vector<TriadData_<_T>>> &gyro_samples_all,
    size_t reference_index,
    int min_recorded_intervals,
    int variance_win_size,
    int th_var_min,
    int th_var_max)
{
    if (!calibrateAcc(acc_samples_all, reference_index, min_recorded_intervals, variance_win_size, th_var_min, th_var_max))
    {
        return false;
    }

    std::cout << "\033[1;32m------------------------------------------------------------------------\033[0m" << std::endl;
    std::cout << "\033[1;32m----------------- NOW PERFORMING GYRO CALIBRATION ----------------------\033[0m" << std::endl;
    std::cout << "\033[1;32m------------------------------------------------------------------------\033[0m" << std::endl;
    if (reference_index >= gyro_samples_all.size())
    {
        return false;
    }
    gyro_calibs_.clear();
    gyro_calibs_.resize(gyro_samples_all.size());

    // Configure per-topic calibrators and use cache
    per_topic_calibs_.resize(gyro_samples_all.size());
    for (size_t i = 0; i < gyro_samples_all.size(); ++i)
    {
        auto &calib_single = per_topic_calibs_[i];
        calib_single.setGravityMagnitude(g_mag_);
        calib_single.setInitStaticIntervalDuration(init_interval_duration_);
        calib_single.setIntarvalsNumSamples(interval_n_samples_);
        calib_single.enableAccUseMeans(acc_use_means_);
        calib_single.setGyroDataPeriod(gyro_dt_);
        if (lowpass_frequency_cut_ > 0)
            calib_single.setLowpassFreq(lowpass_frequency_cut_);
        calib_single.enableVerboseOutput(verbose_output_);
        calib_single.enablePlotting(plot_enabled_);

        if (!calib_single.isGyroCached())
        {
            if (!calib_single.calibrateAccGyro(acc_samples_all[i], gyro_samples_all[i],
                                               min_recorded_intervals,
                                               variance_win_size,
                                               th_var_min,
                                               th_var_max))
                return false;
        }
        gyro_calibs_[i] = calib_single.getGyroCalib();
    }

    // -------------------------------------------------------------------------
    // Gyro-only rotation alignment on SO(3) manifold (axis-angle)
    // M and S are fixed from intrinsic calibration
    // The resulting R is applied to both gyro and acc (same chip)
    // -------------------------------------------------------------------------
    std::cout << "---------------------- ALIGNING GYRO TO REF (SO(3)) ------------------------" << std::endl;

    // Reference calibrated samples
    const auto &ref_gyr_calib = per_topic_calibs_[reference_index].getCalibGyroSamples();

    for (size_t i = 0; i < gyro_samples_all.size(); ++i)
    {
        if (i == reference_index)
            continue;

        std::cout << "Aligning sensor " << i << " to reference " << reference_index << std::endl;

        // Get intrinsic calibrations (M, S, B) for this sensor
        const auto &gyr_calib = per_topic_calibs_[i].getGyroCalib();
        const auto &acc_calib = per_topic_calibs_[i].getAccCalib();

        Eigen::Matrix<_T, 3, 3> M_gyr = gyr_calib.getMisalignmentMatrix();
        Eigen::Matrix<_T, 3, 3> S_gyr = gyr_calib.getScaleMatrix();
        Eigen::Matrix<_T, 3, 1> B_gyr = gyr_calib.getBiasVector();

        Eigen::Matrix<_T, 3, 3> M_acc = acc_calib.getMisalignmentMatrix();
        Eigen::Matrix<_T, 3, 3> S_acc = acc_calib.getScaleMatrix();
        Eigen::Matrix<_T, 3, 1> B_acc = acc_calib.getBiasVector();

        // Initialize rotation from identity (axis-angle = [0,0,0])
        std::vector<double> axis_angle(3, 0.0);

        ceres::Problem problem;

        // Use gyro samples only for rotation optimization
        size_t n = std::min(ref_gyr_calib.size(), gyro_samples_all[i].size());

        for (size_t k = 0; k < n; ++k)
        {
            ceres::CostFunction *cost_function =
                new ceres::AutoDiffCostFunction<GyroRotationAlignResidual<_T>, 3, 3>(
                    new GyroRotationAlignResidual<_T>(
                        ref_gyr_calib[k],       // Reference calibrated gyro
                        gyro_samples_all[i][k], // Other raw gyro
                        M_gyr, S_gyr, B_gyr));  // Fixed gyro intrinsics
            problem.AddResidualBlock(cost_function, nullptr, axis_angle.data());
        }

        ceres::Solver::Options options;
        options.max_num_iterations = 5000;
        options.function_tolerance = 1e-24;
        options.gradient_tolerance = 1e-28;
        options.parameter_tolerance = 1e-26;
        options.max_linear_solver_iterations = 10000;
        options.linear_solver_type = ceres::DENSE_QR;
        options.minimizer_progress_to_stdout = verbose_output_;

        ceres::Solver::Summary summary;
        ceres::Solve(options, &problem, &summary);

        if (verbose_output_)
            std::cout << summary.BriefReport() << std::endl;

        // Convert optimized axis-angle to rotation matrix
        double R_data[9];
        ceres::AngleAxisToRotationMatrix(axis_angle.data(), R_data);
        Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::ColMajor>> R_opt_d(R_data);
        Eigen::Matrix<_T, 3, 3> R_opt = R_opt_d.template cast<_T>();

        // Build final calibrations: T = R * M * S
        // Store as combined misalignment (R*M) and scale S
        Eigen::Matrix<_T, 3, 3> final_mis_gyr = R_opt * M_gyr;
        Eigen::Matrix<_T, 3, 3> final_mis_acc = R_opt * M_acc;

        gyro_calibs_[i] = CalibratedTriad_<_T>(final_mis_gyr, S_gyr, B_gyr);
        acc_calibs_[i] = CalibratedTriad_<_T>(final_mis_acc, S_acc, B_acc);

        // Print rotation info
        double angle = std::sqrt(axis_angle[0] * axis_angle[0] +
                                 axis_angle[1] * axis_angle[1] +
                                 axis_angle[2] * axis_angle[2]);
        std::cout << "  Optimized rotation: axis-angle=[" << axis_angle[0] << ", "
                  << axis_angle[1] << ", " << axis_angle[2] << "], angle="
                  << (angle * 180.0 / M_PI) << " deg" << std::endl;
    }

    return true;
}

template <typename _T>
bool MultiPosRedundantCalibration_<_T>::calibrateAccGyroExtrinsic(
    const std::vector<std::vector<TriadData_<_T>>> &acc_samples_all,
    const std::vector<std::vector<TriadData_<_T>>> &gyro_samples_all,
    size_t reference_index,
    int min_recorded_intervals,
    int variance_win_size,
    int th_var_min,
    int th_var_max)
{
    if (!calibrateAccGyro(acc_samples_all, gyro_samples_all, reference_index,
                          min_recorded_intervals, variance_win_size, th_var_min, th_var_max))
        return false;

    std::cout << "\033[1;32m------------------------------------------------------------------------\033[0m" << std::endl;
    std::cout << "\033[1;32m---------------------- NOW PERFORMING EXTRINSIC ------------------------\033[0m" << std::endl;
    std::cout << "\033[1;32m------------------------------------------------------------------------\033[0m" << std::endl;
    extrinsic_vectors_.assign(acc_samples_all.size(), std::vector<_T>(3, _T(0)));
    std::vector<_T> init_extr = {_T(0), _T(0), _T(0)};

    // Calibrate extrinsic parameters for each sensor individually using per_topic_calibs_
    // Then transform the extrinsic vectors to the reference frame using the misalignment transformation
    for (size_t i = 0; i < acc_samples_all.size(); ++i)
    {
        // Initialize the extrinsic vector with zero (assuming sensor at origin)
        per_topic_calibs_[i].setInitExtrinsicCalibration(init_extr);
        // Call individual extrinsic calibration on per_topic_calibs_[i]
        if (!per_topic_calibs_[i].calibrateAccGyroExtrinsic(
                acc_samples_all[i], gyro_samples_all[i],
                min_recorded_intervals, variance_win_size, th_var_min, th_var_max))
        {
            std::cerr << "Extrinsic calibration failed for sensor " << i << std::endl;
            return false;
        }

        // Get the extrinsic vector computed in sensor i's local frame (before alignment with reference)
        const auto &extr_local = per_topic_calibs_[i].getExtrinsicCalib();
        Eigen::Matrix<_T, 3, 1> extr_i_local(extr_local[0], extr_local[1], extr_local[2]);

        if (i == reference_index)
        {
            // Reference sensor extrinsic vector is used as-is
            extrinsic_vectors_[i] = {extr_local[0], extr_local[1], extr_local[2]};
            // Print to terminal
            if (verbose_output_)
            {
                std::cout << "Reference sensor " << i << " extrinsic vector (local frame): "
                          << extr_local[0] << " " << extr_local[1] << " " << extr_local[2] << std::endl;
            }
        }
        else
        {
            // use foldRotation to get the rotation part to ref
            const auto &after = gyro_calibs_[i];
            Eigen::Matrix<_T, 3, 3> R, M, S;
            Eigen::Matrix<_T, 3, 3> T_gyro = after.getMisalignmentMatrix() * after.getScaleMatrix();

            imu_tk::foldRotation(T_gyro, R, M, S);

            // Apply rotation to the extrinsic vector
            Eigen::Matrix<_T, 3, 1> extr_i_aligned = R * extr_i_local;

            extrinsic_vectors_[i] = {
                extr_i_aligned(0),
                extr_i_aligned(1),
                extr_i_aligned(2)};

            if (verbose_output_)
            {
                std::cout << "Sensor " << i << " extrinsic vector (local frame): "
                          << extr_local[0] << " " << extr_local[1] << " " << extr_local[2] << std::endl;
                std::cout << "Sensor " << i << " extrinsic vector (aligned with ref): "
                          << extr_i_aligned(0) << " " << extr_i_aligned(1) << " " << extr_i_aligned(2) << std::endl;
            }
        }
    }

    // Accelerometer alignment is now done in calibrateAccGyro() with combined SO(3) optimization
    // No separate accelerometer alignment needed here

    return true;
}

// Explicit instantiation
template class MultiPosRedundantCalibration_<double>;
template class MultiPosRedundantCalibration_<float>;
