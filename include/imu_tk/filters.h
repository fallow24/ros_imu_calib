/*
 * imu_tk - Inertial Measurement Unit Toolkit
 *
 *  Copyright (c) 2025, Fabian Arzberger <fabian.arzberger@uni-wuerzburg.de>
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
#include <Eigen/Dense>
#include <Eigen/LU> 
#include <Eigen/SVD> 

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
    void staticIntervalsDetector(const std::vector<TriadData_<_T>> &samples,
                                 _T threshold, std::vector<DataInterval> &intervals,
                                 int win_size = 202);

    // Decompose transformation T = R * S * M where:
    // Decompose T = R * M * S where:
    // - R is a proper rotation matrix (orthogonal, det = +1)
    // - M is a misalignment matrix with 1's on the diagonal
    // - S is a diagonal scaling matrix (applied first to measurements)
    template <typename _T>
    void foldRotation(const Eigen::Matrix<_T, 3, 3> &T,
                      Eigen::Matrix<_T, 3, 3> &R,
                      Eigen::Matrix<_T, 3, 3> &S,
                      Eigen::Matrix<_T, 3, 3> &M);

    // Compute smoothed derivative of a 3D signal using centered derivative-of-Gaussian convolution.
    // Compatible with Eigen matrices and Ceres Jet types.
    template <typename T>
    std::vector<imu_tk::TriadData_<T>> gaussianDerivativeConvolution(
        const std::vector<imu_tk::TriadData_<T>> &data,
        T sigma = T(0.03),
        int kernel_size = 39)
    {
        if (data.empty())
            return {};

        const int half_k = kernel_size / 2;

        std::vector<T> kernel(kernel_size);

        // --- Derivative-of-Gaussian kernel ---
        const T inv_s2 = T(1.0) / (sigma * sigma);
        const T norm_factor = T(1.0) / (sigma * sigma * sigma * ceres::sqrt(T(2.0) * M_PI));

        for (int i = -half_k; i <= half_k; ++i)
        {
            T x = T(i); // sample position (centered)
            kernel[i + half_k] = norm_factor * -x * ceres::exp(-T(0.5) * x * x * inv_s2);
        }

        // --- Ensure zero mean ---
        T mean = std::accumulate(kernel.begin(), kernel.end(), T(0.0)) / T(kernel_size);
        for (auto &v : kernel)
            v -= mean;

        // --- Normalize kernel energy ---
        T energy = T(0.0);
        for (auto v : kernel)
            energy += v * v;
        energy = ceres::sqrt(energy);
        if (energy > T(0.0))
            for (auto &v : kernel)
                v /= energy;

        // --- Convolution ---
        const int N = static_cast<int>(data.size());
        std::vector<imu_tk::TriadData_<T>> derivative;
        derivative.reserve(N);

        for (int i = 0; i < N; ++i)
        {
            Eigen::Matrix<T, 3, 1> deriv;
            deriv.setZero();

            for (int k = -half_k; k <= half_k; ++k)
            {
                int idx = i + k;
                if (idx < 0 || idx >= N)
                    continue;
                deriv += kernel[k + half_k] * data[idx].data();
            }

            // Store derivative with same timestamp
            derivative.emplace_back(data[i].timestamp(), deriv);
        }

        return derivative;
    }

    // Compute smoothing of a 3D signal using centered Gaussian convolution.
    // Compatible with Eigen matrices and Ceres Jet types.
    template <typename T>
    std::vector<imu_tk::TriadData_<T>> gaussianSmoothingConvolution(
        const std::vector<imu_tk::TriadData_<T>> &data,
        T sigma = T(0.03),
        int kernel_size = 39)
    {
        if (data.empty())
            return {};

        const int half_k = kernel_size / 2;
        std::vector<T> kernel(kernel_size);

        // --- Gaussian kernel ---
        const T inv_s2 = T(1.0) / (sigma * sigma);
        const T norm_factor = T(1.0) / (sigma * ceres::sqrt(T(2.0) * M_PI));

        for (int i = -half_k; i <= half_k; ++i)
        {
            T x = T(i); // sample position (centered)
            kernel[i + half_k] = norm_factor * ceres::exp(-T(0.5) * x * x * inv_s2);
        }

        // --- Normalize kernel area ---
        T sum = std::accumulate(kernel.begin(), kernel.end(), T(0.0));
        if (ceres::abs(sum) > T(0.0))
            for (auto &v : kernel)
                v /= sum;

        // --- Convolution ---
        const int N = static_cast<int>(data.size());
        std::vector<imu_tk::TriadData_<T>> smoothed;
        smoothed.reserve(N);

        for (int i = 0; i < N; ++i)
        {
            Eigen::Matrix<T, 3, 1> smoothed_sample;
            smoothed_sample.setZero();

            for (int k = -half_k; k <= half_k; ++k)
            {
                int idx = i + k;
                if (idx < 0 || idx >= N)
                    continue;
                smoothed_sample += kernel[k + half_k] * data[idx].data();
            }

            smoothed.emplace_back(data[i].timestamp(), smoothed_sample);
        }

        return smoothed;
    }

} // end namespace imutk
