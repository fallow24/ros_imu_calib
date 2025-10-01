/* 
 * imu_tk - Inertial Measurement Unit Toolkit
 * 
 *  Copyright (c) 2014, Alberto Pretto <pretto@diag.uniroma1.it>
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 * 
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 * 
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

#pragma once

#include <vector>
#include <ceres/jet.h>
#include "imu_tk/base.h"

namespace imu_tk
{

/**
  * @brief Classify between static and motion intervals checking if for each sample 
  *        of the input signal \f$s\f$ (samples) the local variance magnitude 
  *        is lower or greater then a threshold.
  * 
  * @param samples Input 3D signal (e.g, the acceleremeter readings)
  * @param threshold Threshold used in the classification
  * @param[out] intervals  Ouput detected static intervals
  * @param win_size Size of the sliding window (i.e., number of samples) 
  *                 used to compute the local variance magnitude. It should be equal or
  *                 greater than 11
  * 
  * 
  * The variance magnitude is a scalar computed in a temporal sliding window of size 
  * \f$w_s\f$ (i.e., win_size) as:
  * \f[ \varsigma(t) = 
  *     \sqrt{[var_{w_s}(s^t_x)]^2 + [var_{w_s}(s^t_y)]^2 + [var_{w_s}(s^t_z)]^2} \f] 
  * 
  * Where \f$var_{w_s}(s^t)\f$ is an operator that compute the variance of
  * a general 1D signal in a interval of length \f$w_s\f$ samples
  * centered in \f$t\f$.
  */
template <typename _T> 
  void staticIntervalsDetector ( const std::vector< TriadData_<_T> > &samples,
                                 _T threshold, std::vector< DataInterval > &intervals,
                                 int win_size = 202 );


template<typename _T1>
  std::vector<imu_tk::TriadData_<_T1>> gaussianDerivativeConvolution(
    const std::vector<imu_tk::TriadData_<_T1>>& data,
    _T1 sigma = 1.0)
{
    // Derivative of Gaussian kernel (size 5, normalized, sampled at -2,-1,0,1,2)
    // g'(x) = -x * exp(-x^2/(2*sigma^2)) / (sigma^3 * sqrt(2*pi))
    const int kernel_size = 5;
    const int half_k = kernel_size / 2;
    std::vector<_T1> kernel(kernel_size);
    _T1 sum = _T1(0.0);
    for (int i = -half_k; i <= half_k; ++i) {
        kernel[i + half_k] = _T1(-i) * ceres::exp(-(_T1(i * i)) / (_T1(2.0) * sigma * sigma));
        sum += ceres::abs(kernel[i + half_k]);
    }
    // Normalize the kernel
    for (int i = 0; i < kernel_size; ++i) kernel[i] /= sum;

    int N = data.size();
    std::vector<imu_tk::TriadData_<_T1>> derivative;
    derivative.reserve(N);

    // Convolve kernel with each axis separately
    for (int i = 0; i < N; ++i) {
        Eigen::Matrix<_T1, 3, 1> deriv;
        deriv.setZero();
        _T1 dt_sum = _T1(0.0);
        for (int k = -half_k; k <= half_k; ++k) {
            int idx = i + k;
            if (idx < 0 || idx >= N) continue;
            // Use kernel value and difference in timestamp for scaling
            deriv += kernel[k + half_k] * data[idx].data();
            dt_sum += kernel[k + half_k] * data[idx].timestamp();
        }
        // Estimate dt as average difference between timestamps in kernel window
        _T1 dt = _T1(0.0);
        if (i > 0 && i < N - 1)
            dt = (data[std::min(i + half_k, N - 1)].timestamp() - data[std::max(i - half_k, 0)].timestamp()) / _T1(kernel_size);
        else
            dt = _T1(1.0);
        if (dt != 0)
            deriv /= dt;
        derivative.emplace_back(data[i].timestamp(), deriv);
    }
    return derivative;
}


/*
Helper function that takes a time series and applies smoothing and differentiation. 
*/
// template<typename _T1>
// std::vector<imu_tk::TriadData_<_T1>> movingAverageCentralDiff(
//     const std::vector<imu_tk::TriadData_<_T1>>& data,
//     int window_size = 5) // TODO: fix that magic number
// {
//     int N = data.size();
//     std::vector<imu_tk::TriadData_<_T1>> smoothed;
//     smoothed.reserve(N);

//     // Moving average smoothing
//     for (int i = 0; i < N; ++i) {
//         int start = std::max(0, i - window_size/2);
//         int end = std::min(N - 1, i + window_size/2);
//         int count = end - start + 1;
//         Eigen::Matrix<_T1, 3, 1> avg(0, 0, 0);
//         for (int j = start; j <= end; ++j) {
//             avg += data[j].data();
//         }
//         avg /= static_cast<_T1>(count);
//         smoothed.emplace_back(data[i].timestamp(), avg);
//     }

//     // Central difference derivative
//     std::vector<imu_tk::TriadData_<_T1>> derivative;
//     derivative.reserve(N - 2);
//     for (int i = 1; i < N - 1; ++i) {
//         _T1 dt = smoothed[i + 1].timestamp() - smoothed[i - 1].timestamp();
//         if (dt == 0) continue; // Avoid division by zero
//         Eigen::Matrix<_T1, 3, 1> diff = (smoothed[i + 1].data() - smoothed[i - 1].data()) / dt;
//         derivative.emplace_back(smoothed[i].timestamp(), diff);
//     }
//     return derivative;
// }

}
