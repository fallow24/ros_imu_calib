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

#include "imu_tk/calibration.h"
#include "imu_tk/filters.h"
#include "imu_tk/integration.h"
#include "imu_tk/visualization.h"

#include <limits>
#include <iostream>
#include "ceres/ceres.h"

using namespace imu_tk;
using namespace Eigen;
using namespace std;

template <typename _T1> struct MultiPosAccResidual
{
  MultiPosAccResidual( const _T1 &g_mag, const Eigen::Matrix< _T1, 3 , 1> &sample ) :
  g_mag_(g_mag),
  sample_(sample){}
  
  template <typename _T2>
    bool operator() ( const _T2* const params, _T2* residuals ) const
  {
    Eigen::Matrix< _T2, 3 , 1> raw_samp( _T2(sample_(0)), _T2(sample_(1)), _T2(sample_(2)) );
    /* Assume body frame same as accelerometer frame, 
     * so bottom left params in the misalignment matris are set to zero */
    CalibratedTriad_<_T2> calib_triad( params[0], params[1], params[2], 
                                     _T2(0), _T2(0), _T2(0),
                                     params[3], params[4], params[5], 
                                     params[6], params[7], params[8] );
    
    Eigen::Matrix< _T2, 3 , 1> calib_samp = calib_triad.unbiasNormalize( raw_samp );
    residuals[0] = _T2 ( g_mag_ ) - calib_samp.norm();
    return true;
  }
  
  static ceres::CostFunction* Create ( const _T1 &g_mag, const Eigen::Matrix< _T1, 3 , 1> &sample )
  {
    return ( new ceres::AutoDiffCostFunction< MultiPosAccResidual, 1, 9 > (
               new MultiPosAccResidual<_T1>( g_mag, sample ) ) );
  }
  
  const _T1 g_mag_;
  const Eigen::Matrix< _T1, 3 , 1> sample_;
};

template <typename _T1> struct MultiPosGyroResidual
{
  MultiPosGyroResidual( const Eigen::Matrix< _T1, 3 , 1> &g_versor_pos0, 
                        const Eigen::Matrix< _T1, 3 , 1> &g_versor_pos1,
                        const std::vector< TriadData_<_T1> > &gyro_samples, 
                        const DataInterval &gyro_interval_pos01, 
                        _T1 dt, bool optimize_bias) :

  g_versor_pos0_(g_versor_pos0), 
  g_versor_pos1_(g_versor_pos1),
  gyro_samples_(gyro_samples),
  interval_pos01_(gyro_interval_pos01),
  dt_(dt), optimize_bias_(optimize_bias){}
  
  template <typename _T2>
    bool operator() ( const _T2* const params, _T2* residuals ) const
  {
    CalibratedTriad_<_T2> calib_triad( params[0], params[1], params[2], 
                                      params[3], params[4], params[5], 
                                      params[6], params[7], params[8],
                                      optimize_bias_?params[9]:_T2(0), 
                                      optimize_bias_?params[10]:_T2(0), 
                                      optimize_bias_?params[11]:_T2(0) );

    std::vector< TriadData_<_T2> > calib_gyro_samples;
    calib_gyro_samples.reserve( interval_pos01_.end_idx - interval_pos01_.start_idx + 1 );
    
    for( int i = interval_pos01_.start_idx; i <= interval_pos01_.end_idx; i++ )
      calib_gyro_samples.push_back( TriadData_<_T2>( calib_triad.unbiasNormalize( gyro_samples_[i] ) ) );
    
    Eigen::Matrix< _T2, 3 , 3> rot_mat;
    integrateGyroInterval( calib_gyro_samples, rot_mat, _T2(dt_) );
    
    Eigen::Matrix< _T2, 3 , 1> diff = rot_mat.transpose()*g_versor_pos0_.template cast<_T2>() -
                                      g_versor_pos1_.template cast<_T2>();
    
    residuals[0] = diff(0);
    residuals[1] = diff(1);
    residuals[2] = diff(2);
    
    return true;
  }
  
  static ceres::CostFunction* Create ( const Eigen::Matrix< _T1, 3 , 1> &g_versor_pos0, 
                                       const Eigen::Matrix< _T1, 3 , 1> &g_versor_pos1,
                                       const std::vector< TriadData_<_T1> > &gyro_samples, 
                                       const DataInterval &gyro_interval_pos01, 
                                       _T1 dt, bool optimize_bias )
  {
    if( optimize_bias )
      return ( new ceres::AutoDiffCostFunction< MultiPosGyroResidual, 3, 12 > (
                new MultiPosGyroResidual( g_versor_pos0, g_versor_pos1, gyro_samples, 
                                          gyro_interval_pos01, dt, optimize_bias ) ) );
    else
      return ( new ceres::AutoDiffCostFunction< MultiPosGyroResidual, 3, 9 > (
                new MultiPosGyroResidual( g_versor_pos0, g_versor_pos1, gyro_samples, 
                                          gyro_interval_pos01, dt, optimize_bias ) ) );
  }
  
  const Eigen::Matrix< _T1, 3 , 1> g_versor_pos0_, g_versor_pos1_;
  const std::vector< TriadData_<_T1> > gyro_samples_;
  const DataInterval interval_pos01_;
  const _T1 dt_;
  const bool optimize_bias_;
};

/*
Implements an extrinsic calibration based on rotation (no translation!)
The robot should be rotated around the base without translation.
During this motion, due to the off-centered placement of the IMU, there are 
centripetal- and tangential-forces being measured by the accelerometer.
This can be compensated for by knowing the offset from the base to the IMU. 
*/
template <typename _T1> struct MultiPosExtrinsicResidual
{
  MultiPosExtrinsicResidual(const _T1 &g_mag,
                            const std::vector< TriadData_<_T1> > &calib_gyro_samples, 
                            const std::vector< TriadData_<_T1> > &angular_acceleration,
                            const std::vector< TriadData_<_T1> > &calib_acc_samples, 
                            const DataInterval &interval_pos01) :
    g_mag_(g_mag), 
    gyro_samples_(calib_gyro_samples), 
    derivative_gyro_(angular_acceleration),
    acc_samples_(calib_acc_samples), 
    interval_pos01_(interval_pos01)
  {}

  template <typename _T2>
  bool operator() ( const _T2* const params, _T2* residuals ) const
  {
    Eigen::Matrix<_T2, 3, 1> extr;
    extr << params[0], params[1], params[2];
    residuals[0] = _T2(0.0);

    for( int i = interval_pos01_.start_idx; i <= interval_pos01_.end_idx; i++ ) {
      Eigen::Matrix<_T2, 3, 1> gyr = gyro_samples_[i].data().template cast<_T2>();
      Eigen::Matrix<_T2, 3, 1> acc = acc_samples_[i].data().template cast<_T2>();
      Eigen::Matrix<_T2, 3, 1> derivative_gyro = derivative_gyro_[i].data().template cast<_T2>();
      residuals[0] += ceres::abs( _T2(g_mag_) - (acc - derivative_gyro.cross(extr) - gyr.cross(gyr.cross(extr)) ).norm() ); 
    }

    return true;
  }

  static ceres::CostFunction* Create( const _T1 &g_mag,
                                      const std::vector< TriadData_<_T1> > &gyro_samples, 
                                      const std::vector< TriadData_<_T1> > &angular_acceleration,
                                      const std::vector< TriadData_<_T1> > &acc_samples, 
                                      const DataInterval &interval_pos01)
  {
    return ( new ceres::AutoDiffCostFunction< MultiPosExtrinsicResidual, 1, 3> (
      new MultiPosExtrinsicResidual( g_mag, gyro_samples, angular_acceleration, acc_samples, interval_pos01)
    ));
  }
  
  const _T1 g_mag_;
  const std::vector< TriadData_<_T1> > gyro_samples_;
  const std::vector< TriadData_<_T1> > derivative_gyro_;
  const std::vector< TriadData_<_T1> > acc_samples_;
  const DataInterval interval_pos01_;
};

template <typename _T>
  MultiPosCalibration_<_T>::MultiPosCalibration_() :
  g_mag_(9.81),
  min_num_intervals_(12),
  init_interval_duration_(_T(30.0)),
  interval_n_samples_(100),
  acc_use_means_(false),
  gyro_dt_(-1.0),
  optimize_gyro_bias_(false),
  verbose_output_(false),
  plot_enabled_(true),
  cached_acc_(false),
  cached_gyr_(false),
  cached_extr_(false)
{}

template <typename _T>
  bool MultiPosCalibration_<_T>::calibrateAcc ( const std::vector< TriadData_<_T> >& acc_samples , 
    int min_recorded_intervals,
    int variance_win_size,
    int th_var_min,
    int th_var_max)
{
  cout<<"Accelerometers calibration: calibrating..."<<endl;
  
  min_cost_static_intervals_.clear();
  min_cost_static_samples_.clear();
  calib_acc_samples_.clear();
  calib_gyro_samples_.clear();
  
  int n_samps = acc_samples.size();
  
  DataInterval init_static_interval = DataInterval::initialInterval( acc_samples, init_interval_duration_ );
  Eigen::Matrix<_T, 3, 1> acc_variance = dataVariance( acc_samples, init_static_interval );
  _T norm_th = acc_variance.norm();

  _T min_cost = std::numeric_limits< _T >::max();
  int min_cost_th = -1;
  std::vector< double > min_cost_calib_params;
  
  for (int th_mult = th_var_min; th_mult <= th_var_max; th_mult++)
  {
    std::vector< imu_tk::DataInterval > static_intervals;
    std::vector< imu_tk::TriadData_<_T> > static_samples;
    std::vector< double > acc_calib_params(9);
    
    acc_calib_params[0] = init_acc_calib_.misYZ();
    acc_calib_params[1] = init_acc_calib_.misZY();
    acc_calib_params[2] = init_acc_calib_.misZX();
    
    acc_calib_params[3] = init_acc_calib_.scaleX();
    acc_calib_params[4] = init_acc_calib_.scaleY();
    acc_calib_params[5] = init_acc_calib_.scaleZ();
    
    acc_calib_params[6] = init_acc_calib_.biasX();
    acc_calib_params[7] = init_acc_calib_.biasY();
    acc_calib_params[8] = init_acc_calib_.biasZ();
    
    std::vector< DataInterval > extracted_intervals;
    staticIntervalsDetector ( acc_samples, th_mult*norm_th, static_intervals, variance_win_size);
    extractIntervalsSamples ( acc_samples, static_intervals, 
                              static_samples, extracted_intervals,
                              interval_n_samples_, acc_use_means_ );
    
    if(verbose_output_)
      cout<<"Accelerometers calibration: extracted "<<extracted_intervals.size()
          <<" intervals using threshold multiplier "<<th_mult<<" -> ";
    
    // TODO Perform here a quality test
    if( extracted_intervals.size() < min_num_intervals_)
    {
      if( verbose_output_) cout<<"Not enough intervals, calibration is not possible"<<endl;
      continue;
    }
    
    if( verbose_output_) cout<<"Trying calibrate... "<<endl;
    
    ceres::Problem problem;
    for( int i = 0; i < static_samples.size(); i++)
    {
      ceres::CostFunction* cost_function =
        MultiPosAccResidual<_T>::Create ( g_mag_, static_samples[i].data() );

      problem.AddResidualBlock ( cost_function, NULL /* squared loss */, acc_calib_params.data() ); 
    }
    
    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = verbose_output_;
    // Below options might be overkill but we dont need to optimize for speed anyway
    options.max_num_iterations = 5000;
    options.function_tolerance = 1e-12;
    options.gradient_tolerance = 1e-16;
    options.parameter_tolerance = 1e-14;
    options.max_linear_solver_iterations = 10000;

    ceres::Solver::Summary summary;
    ceres::Solve ( options, &problem, &summary );

    // TODO: Not only look at the minimum cost, but also the quality of extracted intervals!
    if( summary.final_cost < min_cost && extracted_intervals.size() >= min_recorded_intervals)
    {
      min_cost = summary.final_cost;
      min_cost_th = th_mult;
      min_cost_static_intervals_ = extracted_intervals;
      min_cost_static_samples_ = static_samples;
      min_cost_calib_params = acc_calib_params;
    }
    cout<<"residual "<<summary.final_cost<<endl;
  }
  
  if( min_cost_th < 0 )
  {
    if(verbose_output_) 
      cout<<"Accelerometers calibration: Can't obtain any calibratin with the current dataset"<<endl;
    return false;
  }

  acc_calib_ = CalibratedTriad_<_T>( min_cost_calib_params[0],
                                     min_cost_calib_params[1],
                                     min_cost_calib_params[2],
                                     0,0,0,
                                     min_cost_calib_params[3],
                                     min_cost_calib_params[4],
                                     min_cost_calib_params[5],
                                     min_cost_calib_params[6],
                                     min_cost_calib_params[7],
                                     min_cost_calib_params[8] );
  
  calib_acc_samples_.reserve(n_samps);
  
  // Calibrate the input accelerometer data with the obtained calibration
  for( int i = 0; i < n_samps; i++)
    calib_acc_samples_.push_back( acc_calib_.unbiasNormalize( acc_samples[i]) );
  
  if(verbose_output_) 
  {
    if(plot_enabled_)
    {
      Plot plot;
      plot.plotIntervals( calib_acc_samples_, min_cost_static_intervals_);
      
      cout<<"Accelerometers calibration: Better calibration obtained using threshold multiplier "<<min_cost_th
          <<" with residual "<<min_cost<<endl
          <<acc_calib_<<endl
          <<"Accelerometers calibration: inverse scale factors:"<<endl
          <<1.0/acc_calib_.scaleX()<<endl
          <<1.0/acc_calib_.scaleY()<<endl
          <<1.0/acc_calib_.scaleZ()<<endl;
      
      waitForKey();
    }
    else
    {
      cout<<"Accelerometers calibration: Better calibration obtained using threshold multiplier "<<min_cost_th
          <<" with residual "<<min_cost<<endl
          <<acc_calib_<<endl
          <<"Accelerometers calibration: inverse scale factors:"<<endl
          <<1.0/acc_calib_.scaleX()<<endl
          <<1.0/acc_calib_.scaleY()<<endl
          <<1.0/acc_calib_.scaleZ()<<endl;
    }
  }
  cached_acc_ = true; // mark accelerometer calibration as cached
  return true;
    
}

template <typename _T> 
  bool MultiPosCalibration_<_T>::calibrateAccGyro ( const vector< TriadData_<_T> >& acc_samples, 
                                                   const vector< TriadData_<_T> >& gyro_samples,
                                                   int min_recorded_intervals,
                                                   int variance_win_size,
                                                   int th_var_min,
                                                   int th_var_max)
{
  // Only run accelerometer calibration if not already available
  if( !cached_acc_ )
  {
    if( !calibrateAcc( acc_samples, min_recorded_intervals, 
      variance_win_size, th_var_min, th_var_max ) )
      return false;
  }
  
  cout<<"Gyroscopes calibration: calibrating..."<<endl;
  
  std::vector< TriadData_<_T> > static_acc_means;
  std::vector< DataInterval > extracted_intervals;
  extractIntervalsSamples ( calib_acc_samples_, min_cost_static_intervals_, 
                            static_acc_means, extracted_intervals,
                            interval_n_samples_, true );
  
  int n_static_pos = static_acc_means.size(), n_samps = gyro_samples.size();
  
  // Compute the gyroscopes biases in the (static) initialization interval
  DataInterval init_static_interval = DataInterval::initialInterval( gyro_samples, init_interval_duration_ );
  Eigen::Matrix<_T, 3, 1> gyro_bias = dataMean( gyro_samples, init_static_interval );
  
  gyro_calib_ = CalibratedTriad_<_T>(0, 0, 0, 0, 0, 0, 
                                    1.0, 1.0, 1.0, 
                                    gyro_bias(0), gyro_bias(1), gyro_bias(2) );
  

  // calib_gyro_samples_ already cleared in calibrateAcc()
  calib_gyro_samples_.reserve(n_samps);
  // Remove the bias
  for( int i = 0; i < n_samps; i++ )
    calib_gyro_samples_.push_back(gyro_calib_.unbias(gyro_samples[i]));
  
  std::vector< double > gyro_calib_params(12);

  gyro_calib_params[0] = init_gyro_calib_.misYZ();
  gyro_calib_params[1] = init_gyro_calib_.misZY();
  gyro_calib_params[2] = init_gyro_calib_.misZX();
  gyro_calib_params[3] = init_gyro_calib_.misXZ();
  gyro_calib_params[4] = init_gyro_calib_.misXY();
  gyro_calib_params[5] = init_gyro_calib_.misYX();
  
  gyro_calib_params[6] = init_gyro_calib_.scaleX();
  gyro_calib_params[7] = init_gyro_calib_.scaleY();
  gyro_calib_params[8] = init_gyro_calib_.scaleZ();
  
  // Bias has been estimated and removed in the initialization period
  gyro_calib_params[9] = 0.0;
  gyro_calib_params[10] = 0.0;
  gyro_calib_params[11] = 0.0;
  
  ceres::Problem problem;
  
  // Extract movement (between static intervals)
  // Attention: We count 2 indices, one counts intervals, other counts timestamps.
  // Timestamp index should not be reset in the inner loop, so is initialized here.
  for( int i = 0, t_idx = 0; i < n_static_pos - 1; i++ )
  {
    // Get gravity from static positions i and i+1
    Eigen::Matrix<_T, 3, 1> g_versor_pos0 = static_acc_means[i].data(),
                            g_versor_pos1 = static_acc_means[i + 1].data();
    // Normalize the gravity                   
    g_versor_pos0 /= g_versor_pos0.norm();                           
    g_versor_pos1 /= g_versor_pos1.norm();
    
    // Get gyroscope indices between the static intervals
    int gyro_idx0 = -1, gyro_idx1 = -1; // We have to look at the timestamps
    _T ts0 = calib_acc_samples_[extracted_intervals[i].end_idx].timestamp(), 
       ts1 = calib_acc_samples_[extracted_intervals[i + 1].start_idx].timestamp();
     
    // Assign indices based on the timestamps (assume monotone signal time)
    for( ; t_idx < n_samps; t_idx++ ) {
      if( gyro_idx0 < 0 ) {
        if( calib_gyro_samples_[t_idx].timestamp() >= ts0 )
          gyro_idx0 = t_idx;
      } else {
        if( calib_gyro_samples_[t_idx].timestamp() >= ts1 ){
          gyro_idx1 = t_idx - 1;
          break;
        }
      }
    }
    // gyroscope indices corresponding to movement have been found 
    DataInterval gyro_interval(gyro_idx0, gyro_idx1);
    
    ceres::CostFunction* cost_function =
      MultiPosGyroResidual<_T>::Create ( g_versor_pos0, g_versor_pos1, calib_gyro_samples_,
                                         gyro_interval, gyro_dt_, optimize_gyro_bias_ );

    problem.AddResidualBlock ( cost_function, NULL /* squared loss */, gyro_calib_params.data() ); 
  }
  
  ceres::Solver::Options options;
  options.max_num_iterations = 5000;
  options.function_tolerance = 1e-12;
  options.gradient_tolerance = 1e-16;
  options.parameter_tolerance = 1e-14;
  options.max_linear_solver_iterations = 10000;
  options.linear_solver_type = ceres::DENSE_QR;
  options.minimizer_progress_to_stdout = verbose_output_;

  ceres::Solver::Summary summary;

  ceres::Solve ( options, &problem, &summary );
  gyro_calib_ = CalibratedTriad_<_T>( gyro_calib_params[0],
                                     gyro_calib_params[1],
                                     gyro_calib_params[2],
                                     gyro_calib_params[3],
                                     gyro_calib_params[4],
                                     gyro_calib_params[5],
                                     gyro_calib_params[6],
                                     gyro_calib_params[7],
                                     gyro_calib_params[8],
                                     gyro_bias(0) + gyro_calib_params[9],
                                     gyro_bias(1) + gyro_calib_params[10],
                                     gyro_bias(2) + gyro_calib_params[11]);                            

  // Calibrate the input gyroscopes data with the obtained calibration
  for( int i = 0; i < n_samps; i++)
    calib_gyro_samples_.push_back( gyro_calib_.unbiasNormalize( gyro_samples[i]) );
  
  if(verbose_output_) 
  {
    
    cout<<summary.FullReport()<<endl;
    cout<<"Gyroscopes calibration: residual "<<summary.final_cost<<endl
        <<gyro_calib_<<endl
        <<"Gyroscopes calibration: inverse scale factors:"<<endl
        <<1.0/gyro_calib_.scaleX()<<endl
        <<1.0/gyro_calib_.scaleY()<<endl
        <<1.0/gyro_calib_.scaleZ()<<endl;
  }
  cached_gyr_ = true; // mark gyro calibration as cached
  return true;
}

template <typename _T> 
  bool MultiPosCalibration_<_T>::calibrateAccGyroExtrinsic(const vector< TriadData_<_T> >& acc_samples, 
                                                           const vector< TriadData_<_T> >& gyro_samples,
                                                           int min_recorded_intervals,
                                                           int variance_win_size,
                                                           int th_var_min,
                                                           int th_var_max)
{
  // Perform intrinsic calibration first if needed
  if (!cached_gyr_) {
    if (!calibrateAccGyro(acc_samples, gyro_samples,
      min_recorded_intervals, variance_win_size, th_var_min, th_var_max)) 
    {
      cout << "Could not perform intrinsic calibration!" << endl;
      return false;
    }
  } 
  
  // Extracting static intervals (again... not optimal but ok)
  std::vector< TriadData_<_T> > static_acc_means;
  std::vector< DataInterval > extracted_intervals;
  extractIntervalsSamples ( 
    calib_acc_samples_, min_cost_static_intervals_, 
    static_acc_means, extracted_intervals,
    interval_n_samples_, true 
  );
  int n_static_pos = static_acc_means.size(), n_samps = gyro_samples.size();
  
  // Derive gyro 
  std::vector< TriadData_<_T> > dw = gaussianDerivativeConvolution<_T>(calib_gyro_samples_, sigma_, win_size_); 
  std::vector< TriadData_<_T> > smooth_acc = gaussianSmoothingConvolution<_T>(calib_acc_samples_, sigma_, win_size_);

  // Setup and define the optimization problem
  std::vector< double > extrinsic_calib_params(3);
  extrinsic_calib_params[0] = init_extrinsic_vector_[0]; 
  extrinsic_calib_params[1] = init_extrinsic_vector_[1];
  extrinsic_calib_params[2] = init_extrinsic_vector_[2];
  ceres::Problem problem; 

  // Extract movement (between static intervals)
  // Attention: We count 2 indices, one counts intervals, other counts timestamps.
  // Timestamp index should not be reset in the inner loop, so is initialized here.
  for( int i = 0, t_idx = 0; i < n_static_pos - 1; i++ )
  {
    // Get gyroscope indices between the static intervals
    int gyro_idx0 = -1, gyro_idx1 = -1; // We have to look at the timestamps
    _T ts0 = calib_acc_samples_[extracted_intervals[i].end_idx].timestamp(), 
       ts1 = calib_acc_samples_[extracted_intervals[i + 1].start_idx].timestamp();
     
    // Assign indices based on the timestamps (assume monotone signal time)
    for( ; t_idx < n_samps; t_idx++ ) {
      if( gyro_idx0 < 0 ) {
        if( calib_gyro_samples_[t_idx].timestamp() >= ts0 )
          gyro_idx0 = t_idx;
      } else {
        if( calib_gyro_samples_[t_idx].timestamp() >= ts1 ){
          gyro_idx1 = t_idx - 1;
          break;
        }
      }
    }
    // gyroscope indices corresponding to movement have been found 
    DataInterval gyro_interval(gyro_idx0, gyro_idx1);
    
    ceres::CostFunction* cost_function =
      MultiPosExtrinsicResidual<_T>::Create(g_mag_, gyro_samples, dw, smooth_acc, gyro_interval);

    problem.AddResidualBlock ( cost_function, NULL /* squared loss */, extrinsic_calib_params.data() ); 
  }

  // Solve the optimization problem
  ceres::Solver::Options options;
  options.max_num_iterations = 5000;
  options.function_tolerance = 1e-12;
  options.gradient_tolerance = 1e-16;
  options.parameter_tolerance = 1e-14;
  options.max_linear_solver_iterations = 10000;
  options.linear_solver_type = ceres::DENSE_QR;
  options.minimizer_progress_to_stdout = verbose_output_;

  ceres::Solver::Summary summary;

  ceres::Solve ( options, &problem, &summary );
  extrinsic_vector_ = std::vector<_T>{
    static_cast<_T>(extrinsic_calib_params[0]),
    static_cast<_T>(extrinsic_calib_params[1]),
    static_cast<_T>(extrinsic_calib_params[2])  
  };
  
  if(verbose_output_)
  {
    cout << summary.FullReport() << endl;
    cout << "Extrinsic vector obtained with residual " << summary.final_cost << ":" << endl;
    cout << extrinsic_vector_[0] << endl;
    cout << extrinsic_vector_[1] << endl;
    cout << extrinsic_vector_[2] << endl;
  }
  cached_extr_ = true; // mark extrinsic calibration as cached
  return true;
}

template <typename _T>
  const std::vector< DataInterval >& MultiPosCalibration_<_T>::getStaticIntervals(
    const std::vector< TriadData_<_T> >& acc_samples,
    int variance_win_size,
    int th_mult)
{
  // If static intervals are already computed and cached, return them
  if (!min_cost_static_intervals_.empty())
  {
    return min_cost_static_intervals_;
  }
  
  // Otherwise, compute static intervals now
  DataInterval init_static_interval = DataInterval::initialInterval(acc_samples, init_interval_duration_);
  Eigen::Matrix<_T, 3, 1> acc_variance = dataVariance(acc_samples, init_static_interval);
  _T norm_th = acc_variance.norm();
  
  std::vector< DataInterval > static_intervals;
  staticIntervalsDetector(acc_samples, th_mult * norm_th, static_intervals, variance_win_size);
  extractIntervalsSamples(acc_samples, static_intervals,
                         calib_acc_samples_, min_cost_static_intervals_,
                         interval_n_samples_, acc_use_means_);
  
  return min_cost_static_intervals_;
}

template class MultiPosCalibration_<double>;
template class MultiPosCalibration_<float>;