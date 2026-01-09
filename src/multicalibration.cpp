#include "imu_tk/multicalibration.h"
#include "imu_tk/filters.h"
#include "imu_tk/integration.h"
#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include <algorithm>

using namespace imu_tk;
using namespace Eigen;

namespace
{

    template <typename _T>
    Eigen::Matrix<_T, 3, 3> nearestRotation(const Eigen::Matrix<_T, 3, 3> &M)
    {
        Eigen::JacobiSVD<Eigen::Matrix<_T, 3, 3>> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Eigen::Matrix<_T, 3, 3> R = svd.matrixU() * svd.matrixV().transpose();
        if (R.determinant() < _T(0))
        {
            Eigen::Matrix<_T, 3, 3> U = svd.matrixU();
            U.col(2) *= _T(-1);
            R = U * svd.matrixV().transpose();
        }
        return R;
    }

    // Align accelerometer of sensor j to reference (already calibrated) on paired samples
    template <typename _T1>
    struct MultiSensorAccAlignResidual
    {
        MultiSensorAccAlignResidual(const TriadData_<_T1> &ref,
                                    const TriadData_<_T1> &other,
                                    const Eigen::Matrix<_T1, 3, 1> &other_bias)
            : ref_(ref), other_(other), other_bias_(other_bias) {}
        template <typename _T2>
        bool operator()(const _T2 *const params, _T2 *residuals) const
        {
            // params: 9 elements of unconstrained 3x3 transformation matrix (row-major)
            Eigen::Map<const Eigen::Matrix<_T2, 3, 3, Eigen::RowMajor>> T(params);

            // Apply transformation: T * (raw - bias)
            Eigen::Matrix<_T2, 3, 1> v_oth0 = other_.data().template cast<_T2>() - other_bias_.template cast<_T2>();
            Eigen::Matrix<_T2, 3, 1> v_other = T * v_oth0;
            // Reference data should be already calibrated
            Eigen::Matrix<_T2, 3, 1> v_ref = ref_.data().template cast<_T2>();
            // Residual is mere difference
            Eigen::Matrix<_T2, 3, 1> diff = v_other - v_ref;
            residuals[0] = diff(0);
            residuals[1] = diff(1);
            residuals[2] = diff(2);
            return true;
        }
        const TriadData_<_T1> ref_, other_;
        const Eigen::Matrix<_T1, 3, 1> other_bias_;
    };

    // Align gyro of sensor j to reference (already calibrated) on paired samples
    template <typename _T1>
    struct MultiSensorGyroAlignResidual
    {
        MultiSensorGyroAlignResidual(const TriadData_<_T1> &ref,
                                     const TriadData_<_T1> &other,
                                     const Eigen::Matrix<_T1, 3, 1> &other_bias)
            : ref_(ref), other_(other), other_bias_(other_bias) {}
        template <typename _T2>
        bool operator()(const _T2 *const params, _T2 *residuals) const
        {
            // params: 9 elements of unconstrained 3x3 transformation matrix (row-major)
            Eigen::Map<const Eigen::Matrix<_T2, 3, 3, Eigen::RowMajor>> T(params);

            // Apply transformation: T * (raw - bias)
            Eigen::Matrix<_T2, 3, 1> v_oth0 = other_.data().template cast<_T2>() - other_bias_.template cast<_T2>();
            Eigen::Matrix<_T2, 3, 1> v_other = T * v_oth0;
            // Reference data should be already calibrated
            Eigen::Matrix<_T2, 3, 1> v_ref = ref_.data().template cast<_T2>();
            // Residual is mere difference
            Eigen::Matrix<_T2, 3, 1> diff = v_other - v_ref;
            residuals[0] = diff(0);
            residuals[1] = diff(1);
            residuals[2] = diff(2);
            return true;
        }
        const TriadData_<_T1> ref_, other_;
        const Eigen::Matrix<_T1, 3, 1> other_bias_;
    };

    // Simple extrinsic residual aligning centrifugal correction between ref and other
    template <typename _T1>
    struct MultiSensorExtrinsicResidual
    {
        MultiSensorExtrinsicResidual(const TriadData_<_T1> &ref_acc,
                                     const TriadData_<_T1> &ref_gyr,
                                     const TriadData_<_T1> &ref_dw,
                                     const TriadData_<_T1> &acc,
                                     const TriadData_<_T1> &gyr,
                                     const TriadData_<_T1> &dw,
                                     _T1 g_mag)
            : ref_acc_(ref_acc), ref_gyr_(ref_gyr), ref_dw_(ref_dw),
              acc_(acc), gyr_(gyr), dw_(dw), g_mag_(g_mag) {}
        template <typename _T2>
        bool operator()(const _T2 *const params, _T2 *residuals) const
        {
            Eigen::Matrix<_T2, 3, 1> extr(params[0], params[1], params[2]);
            auto comp = [&](const TriadData_<_T1> &a, const TriadData_<_T1> &g, const TriadData_<_T1> &d)
            {
                Eigen::Matrix<_T2, 3, 1> acc = a.data().template cast<_T2>();
                Eigen::Matrix<_T2, 3, 1> w = g.data().template cast<_T2>();
                Eigen::Matrix<_T2, 3, 1> dw = d.data().template cast<_T2>();
                return acc - dw.cross(extr) - w.cross(w.cross(extr));
            };
            Eigen::Matrix<_T2, 3, 1> ref_corr = comp(ref_acc_, ref_gyr_, ref_dw_);
            Eigen::Matrix<_T2, 3, 1> oth_corr = comp(acc_, gyr_, dw_);
            Eigen::Matrix<_T2, 3, 1> diff = oth_corr - ref_corr;
            residuals[0] = diff(0);
            residuals[1] = diff(1);
            residuals[2] = diff(2);
            return true;
        }
        const TriadData_<_T1> ref_acc_, ref_gyr_, ref_dw_, acc_, gyr_, dw_;
        const _T1 g_mag_;
    };

} // namespace

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
      verbose_output_(false)
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

    // Align others to reference using residuals
    std::cout << "-------------------- ALLIGNING ACCELEROMETERS TO REF ------------------------" << std::endl;
    const auto &ref_samples = acc_samples_all[reference_index];
    CalibratedTriad_<_T> ref_calib = acc_calibs_[reference_index];
    for (size_t i = 0; i < acc_samples_all.size(); ++i)
    {
        if (i == reference_index)
            continue;

        // Use 9 parameters: unconstrained 3x3 transformation matrix (row-major)
        std::vector<double> params(9);
        // Initialize with current transformation T = M * S
        Eigen::Matrix<_T, 3, 3> T_init = acc_calibs_[i].getMisalignmentMatrix() * acc_calibs_[i].getScaleMatrix();
        // Store in row-major order
        for (int row = 0; row < 3; ++row)
            for (int col = 0; col < 3; ++col)
                params[row * 3 + col] = static_cast<double>(T_init(row, col));

        const auto &before_acc = per_topic_calibs_[i].getAccCalib();
        Eigen::Matrix<_T, 3, 1> other_bias_acc = before_acc.getBiasVector();

        ceres::Problem problem;
        size_t n = std::min(ref_samples.size(), acc_samples_all[i].size());
        for (size_t k = 0; k < n; ++k)
        {
            ceres::CostFunction *cost_function =
                new ceres::AutoDiffCostFunction<MultiSensorAccAlignResidual<_T>, 3, 9>(
                    new MultiSensorAccAlignResidual<_T>(ref_calib.unbiasNormalize(ref_samples[k]),
                                                        acc_samples_all[i][k],
                                                        other_bias_acc));
            problem.AddResidualBlock(cost_function, nullptr, params.data());
        }
        ceres::Solver::Options options;
        // options.linear_solver_type = ceres::DENSE_QR;
        // options.max_num_iterations = 200;
        // options.minimizer_progress_to_stdout = verbose_output_;

        // Rather use overkill
        options.max_num_iterations = 5000;
        options.function_tolerance = 1e-24;
        options.gradient_tolerance = 1e-28;
        options.parameter_tolerance = 1e-26;
        options.max_linear_solver_iterations = 10000;
        options.linear_solver_type = ceres::DENSE_QR;
        options.minimizer_progress_to_stdout = verbose_output_;
        ceres::Solver::Summary summary;
        ceres::Solve(options, &problem, &summary);

        // Rebuild triad from optimized 3x3 matrix, use foldRotation to extract R, M, S
        Eigen::Matrix<_T, 3, 1> B_before_acc = before_acc.getBiasVector();
        Eigen::Matrix<_T, 3, 3> T_opt;
        for (int row = 0; row < 3; ++row)
            for (int col = 0; col < 3; ++col)
                T_opt(row, col) = static_cast<_T>(params[row * 3 + col]);
        
        // Use foldRotation to decompose T = R * M * S (scale applied first)
        Eigen::Matrix<_T, 3, 3> R_opt, S_opt, M_opt;
        imu_tk::foldRotation(T_opt, R_opt, S_opt, M_opt);
        
        // Apply rotation as additional transformation: final_mis = R * M
        // Store as (R*M) * S with proper misalignment and scale
        Eigen::Matrix<_T, 3, 3> final_mis = R_opt * M_opt;
        acc_calibs_[i] = CalibratedTriad_<_T>(final_mis, S_opt, B_before_acc);
    }
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

    // Align gyros to reference
    std::cout << "---------------------- ALIGNING GYROS TO REF ------------------------" << std::endl;
    const auto &ref_samples = gyro_samples_all[reference_index];
    CalibratedTriad_<_T> ref_calib = gyro_calibs_[reference_index];
    for (size_t i = 0; i < gyro_samples_all.size(); ++i)
    {
        if (i == reference_index)
            continue;

        // Use 9 parameters: unconstrained 3x3 transformation matrix (row-major)
        std::vector<double> params(9);
        // Initialize with current transformation T = M * S
        Eigen::Matrix<_T, 3, 3> T_init = gyro_calibs_[i].getMisalignmentMatrix() * gyro_calibs_[i].getScaleMatrix();
        // Store in row-major order
        for (int row = 0; row < 3; ++row)
            for (int col = 0; col < 3; ++col)
                params[row * 3 + col] = static_cast<double>(T_init(row, col));

        const auto &before_gyr = per_topic_calibs_[i].getGyroCalib();
        Eigen::Matrix<_T, 3, 1> other_bias_gyr = before_gyr.getBiasVector();

        ceres::Problem problem;
        size_t n = std::min(ref_samples.size(), gyro_samples_all[i].size());
        for (size_t k = 0; k < n; ++k)
        {
            ceres::CostFunction *cost_function =
                new ceres::AutoDiffCostFunction<MultiSensorGyroAlignResidual<_T>, 3, 9>(
                    new MultiSensorGyroAlignResidual<_T>(ref_calib.unbiasNormalize(ref_samples[k]),
                                                         gyro_samples_all[i][k],
                                                         other_bias_gyr));
            problem.AddResidualBlock(cost_function, nullptr, params.data());
        }
        ceres::Solver::Options options;

        // options.linear_solver_type = ceres::DENSE_QR;
        // options.max_num_iterations = 200;
        // options.minimizer_progress_to_stdout = verbose_output_;

        // Use overkill options instead:
        options.max_num_iterations = 5000;
        options.function_tolerance = 1e-24;
        options.gradient_tolerance = 1e-28;
        options.parameter_tolerance = 1e-26;
        options.max_linear_solver_iterations = 10000;
        options.linear_solver_type = ceres::DENSE_QR;
        options.minimizer_progress_to_stdout = verbose_output_;

        ceres::Solver::Summary summary;
        ceres::Solve(options, &problem, &summary);

        // Rebuild triad from optimized 3x3 matrix, use foldRotation to extract R, M, S
        Eigen::Matrix<_T, 3, 1> B_before_gyr = before_gyr.getBiasVector();
        Eigen::Matrix<_T, 3, 3> T_opt;
        for (int row = 0; row < 3; ++row)
            for (int col = 0; col < 3; ++col)
                T_opt(row, col) = static_cast<_T>(params[row * 3 + col]);
        
        // Use foldRotation to decompose T = R * M * S (scale applied first)
        Eigen::Matrix<_T, 3, 3> R_opt, S_opt, M_opt;
        imu_tk::foldRotation(T_opt, R_opt, S_opt, M_opt);
        
        // Apply rotation as additional transformation: final_mis = R * M
        // Store as (R*M) * S with proper misalignment and scale
        Eigen::Matrix<_T, 3, 3> final_mis = R_opt * M_opt;
        gyro_calibs_[i] = CalibratedTriad_<_T>(final_mis, S_opt, B_before_gyr);
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
            // For non-reference sensors, decompose the after calibration directly
            // T = M_after * S_after = Rotation * M_fold * S_fold 

            const auto &after = gyro_calibs_[i];
            const auto &before = per_topic_calibs_[i].getGyroCalib();

            Eigen::Matrix<_T, 3, 3> M_after = after.getMisalignmentMatrix();
            Eigen::Matrix<_T, 3, 3> S_after = after.getScaleMatrix();
            Eigen::Matrix<_T, 3, 3> M_before = before.getMisalignmentMatrix();
            Eigen::Matrix<_T, 3, 3> S_before = before.getScaleMatrix();

            // Diff transform from before to after (applies to extrinsic vector)
            Eigen::Matrix<_T, 3, 3> T_diff = (M_after * S_after) * (M_before * S_before).inverse();

            // Apply transformation to the extrinsic vector
            Eigen::Matrix<_T, 3, 1> extr_i_aligned = T_diff * extr_i_local;

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

    return true;
}

// Explicit instantiation
template class MultiPosRedundantCalibration_<double>;
template class MultiPosRedundantCalibration_<float>;
