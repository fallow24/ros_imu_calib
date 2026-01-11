#include <iostream>
#include <fstream>
#include <yaml-cpp/yaml.h>

#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#include <rosbag/bag.h>
#include <rosbag/view.h>
#include <sensor_msgs/Imu.h>
#include <boost/foreach.hpp>

#include "imu_tk/io_utils.h"
#include "imu_tk/calibration.h"
#include "imu_tk/filters.h"
#include "imu_tk/integration.h"
#include "imu_tk/visualization.h"

using namespace std;
using namespace imu_tk;
using namespace Eigen;

#define foreach BOOST_FOREACH

int main(int argc, char** argv)
{
  if (argc < 2) {
    cout << "Usage: rosrun ros_imu_calib bag_calib <config.yaml>" << endl;
    cout << "Have a look at \"catkin_ws/src/imu_tk/config/imu_calib.yaml\"" << endl;
    return -1;
  }

  // Load YAML config
  YAML::Node config = YAML::LoadFile(argv[1]);

  // program parameters
  string bag_file      = config["bag_file"].as<string>();
  string imu_topic     = config["imu_topic"].as<string>();
  int num_samples_per_interval = config["num_samples_per_interval"].as<int>();
  int imu_frequency_hz         = config["imu_frequency_hz"].as<int>();
  int min_nr_recorded_intervals = config["min_nr_recorded_intervals"].as<int>();
  int variance_win_size         = config["variance_win_size"].as<int>();
  int th_variance_min          = config["th_variance_min"].as<int>();
  int th_variance_max          = config["th_variance_max"].as<int>();
  double initial_static_interval = config["initial_static_interval"].as<double>();
  double gravity_magnitude      = config["gravity_magnitude"].as<double>();
  bool use_acc_means            = config["use_acc_means"].as<bool>();
  bool plot_inspect_static_intervals = config["plot_inspect_static_intervals"] ? config["plot_inspect_static_intervals"].as<bool>() : true;
  string save_acc_calib         = config["save_acc_calib"].as<string>();
  string save_gyr_calib         = config["save_gyr_calib"].as<string>();

  vector< TriadData > acc_data, gyro_data;
  
  cout << "Reading the bagfile";
  cout.flush();
  rosbag::Bag bag;
  bag.open(bag_file.c_str(), rosbag::bagmode::Read);

  std::vector<std::string> topics;
  topics.push_back(imu_topic);
  
  int i = 0;
  rosbag::View view(bag, rosbag::TopicQuery(topics));
  foreach(rosbag::MessageInstance const m, view) {
    sensor_msgs::Imu::ConstPtr s = m.instantiate<sensor_msgs::Imu>();
    if (s != NULL) {
      if (i++ % 1000 == 0) {
        cout << "."; // Notify user that something's happening
        cout.flush();
      }    
      acc_data.push_back(
        TriadData(
          s->header.stamp.toSec(), 
          s->linear_acceleration.x, 
          s->linear_acceleration.y, 
          s->linear_acceleration.z
        )
      );
      gyro_data.push_back(
        TriadData(
          s->header.stamp.toSec(), 
          s->angular_velocity.x, 
          s->angular_velocity.y, 
          s->angular_velocity.z
        )
      );
    }
  }
  cout << endl << "Read " << i << " messages" << endl;
  bag.close();
  
  std::cout << std::endl << "Now performing calibration" << std::endl << std::endl;
  CalibratedTriad init_acc_calib, init_gyro_calib;
  init_acc_calib.setBias( Vector3d(0, 0, 0) );
  init_gyro_calib.setScale( Vector3d(1, 1, 1) );
  MultiPosCalibration mp_calib;  
  mp_calib.setInitStaticIntervalDuration(initial_static_interval);
  mp_calib.setIntarvalsNumSamples(num_samples_per_interval);
  mp_calib.setInitAccCalibration( init_acc_calib );
  mp_calib.setInitGyroCalibration( init_gyro_calib );
  mp_calib.setGravityMagnitude(gravity_magnitude);
  mp_calib.enableVerboseOutput(true);
  mp_calib.enablePlotting(plot_inspect_static_intervals);
  mp_calib.enableAccUseMeans(use_acc_means);
  mp_calib.setGyroDataPeriod(1.0 / imu_frequency_hz);
  mp_calib.calibrateAccGyro(acc_data, gyro_data, 
    min_nr_recorded_intervals, variance_win_size, 
    th_variance_min, th_variance_max);
  mp_calib.getAccCalib().save(save_acc_calib);
  mp_calib.getGyroCalib().save(save_gyr_calib);
  
  return 0;
}

