#include <iostream>
#include <fstream>
#include <yaml-cpp/yaml.h>
#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#include <rosbag/bag.h>
#include <rosbag/view.h>
#include <sensor_msgs/Imu.h>
#include <boost/foreach.hpp>
#include <algorithm>
#include <unordered_map>

#include "imu_tk/io_utils.h"
#include "imu_tk/calibration.h"
#include "imu_tk/filters.h"
#include "imu_tk/integration.h"
#include "imu_tk/visualization.h"
#include "imu_tk/multicalibration.h"

using namespace std;
using namespace imu_tk;
using namespace Eigen;

#define foreach BOOST_FOREACH

int main(int argc, char **argv)
{
  if (argc < 2)
  {
    cout << "Usage: rosrun ros_imu_calib bag_calib_multi <config.yaml>" << endl;
    cout << "Have a look at \"catkin_ws/src/imu_tk/config/imu_calib_multi.yaml\"" << endl;
    return -1;
  }

  // Load YAML config
  YAML::Node config = YAML::LoadFile(argv[1]);

  // program parameters
  string bag_file = config["bag_file"].as<string>();
  std::vector<std::string> imu_topics = config["imu_topics"].as<std::vector<std::string>>();
  std::string imu_ref_topic = config["imu_ref_topic"].as<std::string>();
  int num_samples_per_interval = config["num_samples_per_interval"].as<int>();
  int imu_frequency_hz = config["imu_frequency_hz"].as<int>();
  int min_nr_recorded_intervals = config["min_nr_recorded_intervals"].as<int>();
  int variance_win_size = config["variance_win_size"].as<int>();
  int th_variance_min = config["th_variance_min"].as<int>();
  int th_variance_max = config["th_variance_max"].as<int>();
  double initial_static_interval = config["initial_static_interval"].as<double>();
  double gravity_magnitude      = config["gravity_magnitude"].as<double>();
  bool use_acc_means            = config["use_acc_means"].as<bool>();
  double lowpass_hz             = config["lowpass_hz"].as<double>();
  std::string save_calib_dir    = config["save_calib_dir"].as<std::string>();
  bool calibrate_extrinsic      = config["calibrate_extrinsic"] ? config["calibrate_extrinsic"].as<bool>() : true;

  // Validate ref topic
  if (std::find(imu_topics.begin(), imu_topics.end(), imu_ref_topic) == imu_topics.end())
  {
    std::cerr << "imu_ref_topic is not listed in imu_topics. Exiting." << std::endl;
    return -1;
  }

  // Per-topic datasets
  std::vector<std::vector<TriadData>> acc_data_all(imu_topics.size()), gyro_data_all(imu_topics.size());
  // Index of the reference to calibrate to
  size_t ref_index = static_cast<size_t>(std::distance(imu_topics.begin(), std::find(imu_topics.begin(), imu_topics.end(), imu_ref_topic)));
  // Looping ROS MessageInstance requires a map from topic to idx in the vector
  std::unordered_map<std::string, size_t> topic_to_idx;
  for (size_t idx = 0; idx < imu_topics.size(); ++idx)
  {
    topic_to_idx[imu_topics[idx]] = idx;
  }

  cout << "Reading the bagfile";
  cout.flush();
  rosbag::Bag bag;
  bag.open(bag_file.c_str(), rosbag::bagmode::Read);

  std::vector<std::string> topics = imu_topics;
  int i = 0;
  rosbag::View view(bag, rosbag::TopicQuery(topics));
  foreach (rosbag::MessageInstance const m, view)
  {
    // Route data to the correct topic dataset
    auto it = topic_to_idx.find(m.getTopic());
    if (it == topic_to_idx.end())
    {
      continue; // Topic not found
    }
    size_t idx = it->second;

    sensor_msgs::Imu::ConstPtr s = m.instantiate<sensor_msgs::Imu>();
    if (s != NULL)
    {
      if (i++ % 1000 == 0)
      {
        cout << ".";
        cout.flush();
      }
      acc_data_all[idx].push_back(
          TriadData(
              s->header.stamp.toSec(),
              s->linear_acceleration.x,
              s->linear_acceleration.y,
              s->linear_acceleration.z));
      gyro_data_all[idx].push_back(
          TriadData(
              s->header.stamp.toSec(),
              s->angular_velocity.x,
              s->angular_velocity.y,
              s->angular_velocity.z));
    }
  }
  cout << endl
       << "Read " << i << " messages across " << topics.size() << " topics. Ref index=" << ref_index << " (" << imu_ref_topic << ")" << endl;
  bag.close();

  std::cout << std::endl
            << "Now performing multi-IMU calibration to align all IMUs to reference" << std::endl
            << std::endl;

  // Configure and run redundant multi-position calibration
  MultiPosRedundantCalibration mpr_calib;
  // Verbose similar to calibration class
  mpr_calib.enableVerboseOutput(true);
  mpr_calib.setGravityMagnitude(gravity_magnitude);
  mpr_calib.setInitStaticIntervalDuration(initial_static_interval);
  mpr_calib.setIntarvalsNumSamples(num_samples_per_interval);
  mpr_calib.enableAccUseMeans(use_acc_means);
  mpr_calib.setGyroDataPeriod(1.0 / imu_frequency_hz);
  mpr_calib.setLowpassFreq(static_cast<int>(lowpass_hz));

  // Calibrate all topics and align to reference triad
  bool ok = false;
  if (calibrate_extrinsic)
  {
    ok = mpr_calib.calibrateAccGyroExtrinsic(acc_data_all, gyro_data_all,
                                             ref_index,
                                             min_nr_recorded_intervals, variance_win_size,
                                             th_variance_min, th_variance_max);
  }
  else
  {
    ok = mpr_calib.calibrateAccGyro(acc_data_all, gyro_data_all,
                                    ref_index,
                                    min_nr_recorded_intervals, variance_win_size,
                                    th_variance_min, th_variance_max);
  }
  if (!ok)
  {
    std::cerr << "Calibration failed." << std::endl;
    return -1;
  }

  // Helper to sanitize topic names (remove '/')
  auto sanitize_topic = [](const std::string &topic)
  {
    std::string out;
    out.reserve(topic.size());
    for (char c : topic)
      if (c != '/')
        out.push_back(c);
    return out;
  };

  const auto &acc_calibs = mpr_calib.getAccCalib();
  const auto &gyro_calibs = mpr_calib.getGyroCalib();

  // Strip tailing "/" from save dir if present
  if (!save_calib_dir.empty() && save_calib_dir.back() == '/')
  {
    save_calib_dir.pop_back();
  }

  // Save per-topic calibration files into save_calib_dir
  for (size_t idx = 0; idx < imu_topics.size(); ++idx)
  {
    std::string base = sanitize_topic(imu_topics[idx]);
    std::string acc_path = save_calib_dir + "/acc_" + base + ".txt";
    std::string gyr_path = save_calib_dir + "/gyr_" + base + ".txt";
    std::string extr_path = save_calib_dir + "/extr_" + base + ".txt";

    acc_calibs[idx].save(acc_path);
    gyro_calibs[idx].save(gyr_path);

    if (calibrate_extrinsic)
    {
      const auto &extr_all = mpr_calib.getExtrinsicCalib();
      std::ofstream ef(extr_path.c_str());
      if (ef.is_open())
      {
        ef << "extrinsic: [" << extr_all[idx][0] << ", "
           << extr_all[idx][1] << ", "
           << extr_all[idx][2] << "]" << std::endl;
      }
      else
      {
        std::cerr << "Failed to open extrinsic output file: " << extr_path << std::endl;
      }
    }
  }

  return 0;
}
