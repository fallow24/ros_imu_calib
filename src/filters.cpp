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

#include "imu_tk/filters.h"

using namespace Eigen;

template <typename _T>
void imu_tk::foldRotation(const Eigen::Matrix<_T, 3, 3> &T,
                          Eigen::Matrix<_T, 3, 3> &R,
                          Eigen::Matrix<_T, 3, 3> &M,
                          Eigen::Matrix<_T, 3, 3> &S)
{
  // Perform polar decomposition: T = U * P where U is orthogonal, P is symmetric positive-definite
  Eigen::JacobiSVD<Eigen::Matrix<_T, 3, 3>> svd(T, Eigen::ComputeFullU | Eigen::ComputeFullV);
  Eigen::Matrix<_T, 3, 3> U = svd.matrixU() * svd.matrixV().transpose();
  Eigen::Matrix<_T, 3, 3> P = svd.matrixV() * svd.singularValues().asDiagonal() * svd.matrixV().transpose();

  // Ensure proper rotation (det = +1, not -1)
  if (U.determinant() < _T(0))
  {
    U.col(2) *= _T(-1);
  }
  R = U;

  // Extract pure diagonal scaling from singular values
  S = Eigen::Matrix<_T, 3, 3>::Zero();
  S(0, 0) = svd.singularValues()(0);
  S(1, 1) = svd.singularValues()(1);
  S(2, 2) = svd.singularValues()(2);

  // We want: P = M_final * S (scale applied first)
  // So: M_final = P * S^{-1}
  Eigen::Matrix<_T, 3, 3> M_temp = P * S.inverse();

  // Extract diagonal from M_temp
  _T diag_m0 = M_temp(0, 0);
  _T diag_m1 = M_temp(1, 1);
  _T diag_m2 = M_temp(2, 2);

  // Create M with 1's on diagonal
  M = M_temp;
  M(0, 0) = _T(1);
  M(1, 1) = _T(1);
  M(2, 2) = _T(1);

  // Fold the diagonal scaling into S
  S(0, 0) *= diag_m0;
  S(1, 1) *= diag_m1;
  S(2, 2) *= diag_m2;


  // Iterative refinement to reduce numerical errors
  // We want to minimize ||T - R * M * S|| while keeping constraints
  for (int iter = 0; iter < 3; ++iter)
  {
    Eigen::Matrix<_T, 3, 3> T_approx = R * M * S;
    Eigen::Matrix<_T, 3, 3> error_matrix = T - T_approx;

    // Correct M: fold off-diagonal errors back
    // M_corrected = M + R^T * error_matrix * S^{-1}
    Eigen::Matrix<_T, 3, 3> correction = R.transpose() * error_matrix * S.inverse();
    M += correction;

    // Re-normalize M diagonal to 1
    _T diag_m0_new = M(0, 0);
    _T diag_m1_new = M(1, 1);
    _T diag_m2_new = M(2, 2);

    M(0, 0) = _T(1);
    M(1, 1) = _T(1);
    M(2, 2) = _T(1);

    // Fold corrected diagonal back into S
    S(0, 0) *= diag_m0_new;
    S(1, 1) *= diag_m1_new;
    S(2, 2) *= diag_m2_new;
  }

  // Check if R is a proper rotation by verifying orthogonality and determinant
  Eigen::Matrix<_T, 3, 3> should_be_identity = R.transpose() * R;
  _T identity_error = (should_be_identity - Eigen::Matrix<_T, 3, 3>::Identity()).norm();
  _T det_R = R.determinant();
  if (identity_error > _T(1e-5) || std::abs(det_R - _T(1)) > _T(1e-5))
  {
    std::cerr << "Warning: Resulting R is not a proper rotation matrix. "
              << "Orthogonality error: " << identity_error
              << ", Determinant: " << det_R << std::endl;
  }

  // Check reconstruction error
  Eigen::Matrix<_T, 3, 3> T_reconstructed = R * M * S;
  _T reconstruction_error = (T - T_reconstructed).norm();
  if (reconstruction_error > _T(1e-5))
  {
    std::cerr << "Warning: High reconstruction error in foldRotation: "
              << reconstruction_error << std::endl;
  }

}

template <typename _T>
void imu_tk::staticIntervalsDetector(const std::vector<imu_tk::TriadData_<_T>> &samples,
                                     _T threshold, std::vector<imu_tk::DataInterval> &intervals,
                                     int win_size)
{
  if (win_size < 22)
    win_size = 22;
  if (!(win_size % 2))
    win_size++;

  int h = win_size / 2;

  if (win_size >= samples.size())
    return;

  intervals.clear();

  bool look_for_start = true;
  imu_tk::DataInterval current_interval;

  for (int i = h; i < samples.size() - h; i++)
  {
    Matrix<_T, 3, 1> variance = dataVariance(samples, DataInterval(i - h, i + h));
    _T norm = variance.norm();

    if (look_for_start)
    {
      if (norm < threshold)
      {
        current_interval.start_idx = i;
        look_for_start = false;
      }
    }
    else
    {
      if (norm >= threshold)
      {
        current_interval.end_idx = i - 1;
        look_for_start = true;
        intervals.push_back(current_interval);
      }
    }
  }

  // If the last interval has not been included in the intervals vector
  if (!look_for_start)
  {
    current_interval.end_idx = samples.size() - h - 1;
    // current_interval.end_ts = samples[current_interval.end_idx].timestamp();
    intervals.push_back(current_interval);
  }
}

template void imu_tk::staticIntervalsDetector<double>(const std::vector<TriadData_<double>> &samples,
                                                      double threshold, std::vector<DataInterval> &intervals,
                                                      int win_size = 101);
template void imu_tk::staticIntervalsDetector<float>(const std::vector<TriadData_<float>> &samples,
                                                     float threshold, std::vector<DataInterval> &intervals,
                                                     int win_size = 101);

template void imu_tk::foldRotation<double>(const Eigen::Matrix<double, 3, 3> &,
                                           Eigen::Matrix<double, 3, 3> &,
                                           Eigen::Matrix<double, 3, 3> &,
                                           Eigen::Matrix<double, 3, 3> &);
template void imu_tk::foldRotation<float>(const Eigen::Matrix<float, 3, 3> &,
                                          Eigen::Matrix<float, 3, 3> &,
                                          Eigen::Matrix<float, 3, 3> &,
                                          Eigen::Matrix<float, 3, 3> &);
