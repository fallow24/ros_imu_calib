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

#include <limits>
#include <ceres/rotation.h>

#include "imu_tk/visualization.h"

#include <cstdio>
#include <sstream>

#include "ros/ros.h"

#include <boost/utility.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/concept_check.hpp>

using namespace imu_tk;

static int tmp_argc = 1;
static char *tmp_argv[] = { (char *)"" };

class Plot::PlotImpl
{
public:
  
  PlotImpl()
  { 
    gnuplot_pipe_ = popen("gnuplot", "w");
    if( gnuplot_pipe_ == NULL )
      std::cerr<<"WARNING : missing gnuplot application!"<<std::endl;
  };
  
  ~PlotImpl()
  { 
    if( gnuplot_pipe_ != NULL )
      pclose( gnuplot_pipe_ );
  };
  
  bool ready() const { return gnuplot_pipe_ != NULL; };
  void write( const std::string &str )
  {
    if( gnuplot_pipe_ != NULL )
    {
      fprintf( gnuplot_pipe_, "%s\n", str.c_str() );
      fflush( gnuplot_pipe_ );
    }
  };  
  
private:
  
  FILE *gnuplot_pipe_;
};

Plot::Plot()
{
  plot_impl_ptr_ = boost::shared_ptr< PlotImpl > ( new PlotImpl() );
}

template <typename _T> 
  void Plot::plotSamples ( const std::vector< TriadData_<_T> >& samples, 
                                   DataInterval range )
{
  if( !plot_impl_ptr_->ready() )
  {
    std::cerr<<"WARNING : can't plot samples: check that gnuplot is installed in your system!"<<std::endl;
    return;
  }
  
  range = checkInterval( samples, range );
  
  std::stringstream strs;
  strs<<"plot '-' title 'x' with lines, "
      <<"'-' title 'y' with lines, "
      <<"'-' title 'z' with lines"<<std::endl;
   
  _T base_time = samples[0].timestamp();
  for( int i = range.start_idx; i <= range.end_idx; i++)
    strs<<double(samples[i].timestamp() - base_time)<<" "<<double(samples[i].x())<<std::endl;
  strs<<"EOF"<<std::endl;
  for( int i = range.start_idx; i <= range.end_idx; i++)
    strs<<double(samples[i].timestamp() - base_time)<<" "<<double(samples[i].y())<<std::endl;
  strs<<"EOF"<<std::endl;
  for( int i = range.start_idx; i <= range.end_idx; i++)
    strs<<double(samples[i].timestamp() - base_time)<<" "<<double(samples[i].z())<<std::endl;
  strs<<"EOF"<<std::endl;
  plot_impl_ptr_->write( strs.str() );
}

template <typename _T> 
  void Plot::plotIntervals ( const std::vector< TriadData_< _T > >& samples,
                                     const std::vector< DataInterval >& intervals, 
                                     DataInterval range )
{
  if( !plot_impl_ptr_->ready() )
  {
    std::cerr<<"WARNING : can't plot samples: check that gnuplot is installed in your system!"<<std::endl;
    return;
  }
   
  range = checkInterval( samples, range );
  int n_pts = range.end_idx - range.start_idx + 1, 
              n_intervals = intervals.size();
  
  double max = 0, mean = 0;
  for( int i = range.start_idx; i <= range.end_idx; i++)
  {
    if( double(samples[i].x()) > max ) max = double(samples[i].x());
    if( double(samples[i].y()) > max ) max = double(samples[i].y());
    if( double(samples[i].z()) > max ) max = double(samples[i].z());
    
    mean += (double(samples[i].x()) + double(samples[i].y()) + double(samples[i].z()))/3;
  }
  
  mean /= n_pts;
  max -= mean;
  double step_h = mean + max/2, val = 0;
  int interval_idx = 0;
  for( ; interval_idx < n_intervals; interval_idx++ )
  {
    if (intervals[interval_idx].start_idx >= range.start_idx )
      break;
  }
  
  std::stringstream strs;
  strs<<"plot '-' title 'x' with lines, "
      <<"'-' title 'y' with lines, "
      <<"'-' title 'z' with lines, "
      <<"'-' title 'intervals' with lines"<<std::endl;
   
  _T base_time = samples[0].timestamp();
  for( int i = range.start_idx; i <= range.end_idx; i++)
    strs<<double(samples[i].timestamp() - base_time)<<" "<<double(samples[i].x())<<std::endl;
  strs<<"EOF"<<std::endl;
  for( int i = range.start_idx; i <= range.end_idx; i++)
    strs<<double(samples[i].timestamp() - base_time)<<" "<<double(samples[i].y())<<std::endl;
  strs<<"EOF"<<std::endl;
  for( int i = range.start_idx; i <= range.end_idx; i++)
    strs<<double(samples[i].timestamp() - base_time)<<" "<<double(samples[i].z())<<std::endl;
  strs<<"EOF"<<std::endl;
  
  for( int i = range.start_idx; i <= range.end_idx; i++)
  {
    if( interval_idx < n_intervals)
    {
      if( i == intervals[interval_idx].start_idx )
        val = step_h;
      else if( i == intervals[interval_idx].end_idx)
      {
        val = 0;
        interval_idx++;
      }
    }
    strs<<double(samples[i].timestamp() - base_time)<<" "<<val<<std::endl;
  }
  strs<<"EOF"<<std::endl;
  
  plot_impl_ptr_->write( strs.str() );
}

void imu_tk::waitForKey()
{
  do
  {
    std::cout << "\033[38;5;208m"; // Set term orange
    std::cout << "Inspect the static intervals!";
    std::cout<<std::endl<<"Press Enter to continue, Ctrl+C to abort."<<std::endl;
    std::cout << "\033[0m"; // Set term black
    std::cout.flush();
  }
  while ( getchar() != '\n' );
  
}

template void Plot::plotSamples<double> ( const std::vector< TriadData_<double> >& samples, 
                                          DataInterval range );
template void Plot::plotSamples<float> ( const std::vector< TriadData_<float> >& samples, 
                                         DataInterval range );
template void Plot::plotIntervals<double> ( const std::vector< TriadData_<double> >& samples, 
                                            const std::vector< DataInterval >& intervals,
                                            DataInterval range );
template void Plot::plotIntervals<float> ( const std::vector< TriadData_<float> >& samples, 
                                           const std::vector< DataInterval >& intervals,
                                           DataInterval range );