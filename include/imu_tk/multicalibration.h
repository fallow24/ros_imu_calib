#pragma once
#include <vector>
#include <cstddef>
#include <stdexcept> // added
#include <cmath>     // added
#include "imu_tk/calibration.h"

namespace imu_tk
{

    template <typename _T>
    class MultiPosRedundantCalibration_
    {
    public:
        MultiPosRedundantCalibration_();
        ~MultiPosRedundantCalibration_() = default;

        // Configuration API (mirrors MultiPosCalibration)
        void setGravityMagnitude(_T g) { g_mag_ = g; }
        _T setInitStaticIntervalDuration(_T duration)
        {
            init_interval_duration_ = duration;
            return duration;
        }
        int setIntarvalsNumSamples(int num)
        {
            interval_n_samples_ = num;
            return num;
        }
        void enableAccUseMeans(bool enabled) { acc_use_means_ = enabled; }
        void setGyroDataPeriod(_T dt) { gyro_dt_ = dt; }
        void setLowpassFreq(int val)
        {
            lowpass_frequency_cut_ = val;
            sigma_ = _T(1.0) / (_T(2.0) * _T(M_PI) * _T(lowpass_frequency_cut_));
            if (gyro_dt_ <= _T(0))
                throw std::runtime_error("setGyroDataPeriod() must be called before setLowpassFreq() !");
            win_size_ = static_cast<int>(_T(6.0) * sigma_ / gyro_dt_);
            if (win_size_ % 2 == 0)
                win_size_ += 1;
        }
        void enableVerboseOutput(bool enabled) { verbose_output_ = enabled; }

        // Multi-topic calibration
        bool calibrateAcc(const std::vector<std::vector<TriadData_<_T>>> &acc_samples_all,
                          size_t reference_index,
                          int min_recorded_intervals = 9,
                          int variance_win_size = 202,
                          int th_var_min = 2,
                          int th_var_max = 10);

        bool calibrateAccGyro(const std::vector<std::vector<TriadData_<_T>>> &acc_samples_all,
                              const std::vector<std::vector<TriadData_<_T>>> &gyro_samples_all,
                              size_t reference_index,
                              int min_recorded_intervals = 9,
                              int variance_win_size = 202,
                              int th_var_min = 2,
                              int th_var_max = 10);

        bool calibrateAccGyroExtrinsic(const std::vector<std::vector<TriadData_<_T>>> &acc_samples_all,
                                       const std::vector<std::vector<TriadData_<_T>>> &gyro_samples_all,
                                       size_t reference_index,
                                       int min_recorded_intervals = 9,
                                       int variance_win_size = 202,
                                       int th_var_min = 2,
                                       int th_var_max = 10);

        // Accessors returning per-topic results aligned to reference triad
        const std::vector<CalibratedTriad_<_T>> &getAccCalib() const { return acc_calibs_; }
        const std::vector<CalibratedTriad_<_T>> &getGyroCalib() const { return gyro_calibs_; }
        const std::vector<std::vector<_T>> &getExtrinsicCalib() const { return extrinsic_vectors_; }

    private:
        // Per-topic outputs (post allignment)
        std::vector<CalibratedTriad_<_T>> acc_calibs_;
        std::vector<CalibratedTriad_<_T>> gyro_calibs_;
        std::vector<std::vector<_T>> extrinsic_vectors_;

        // Configuration (mirrors MultiPosCalibration)
        _T g_mag_;
        _T init_interval_duration_;
        int interval_n_samples_;
        bool acc_use_means_;
        _T gyro_dt_;
        int lowpass_frequency_cut_;
        _T sigma_;
        int win_size_;
        bool verbose_output_;

        // Persistent per-topic intrinsic calibrators to reuse caches (pre alignment!)
        std::vector<MultiPosCalibration_<_T>> per_topic_calibs_;
    };

    typedef MultiPosRedundantCalibration_<double> MultiPosRedundantCalibration;

} // namespace imu_tk
