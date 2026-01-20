#include <iostream>
#include <fstream>
#include <yaml-cpp/yaml.h>
#define BOOST_BIND_GLOBAL_PLACEHOLDERS // suppress deprecation warning from bag header
#include <rosbag/bag.h>
#include <rosbag/view.h>
#include <sensor_msgs/Imu.h>
#include <algorithm>
#include <unordered_map>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <cstdlib>
#include <string>

#include "imu_tk/io_utils.h"
#include "imu_tk/calibration.h"
#include "imu_tk/filters.h"
#include "imu_tk/integration.h"
#include "imu_tk/visualization.h"
#include "imu_tk/multicalibration.h"

using namespace std;
using namespace imu_tk;
using namespace Eigen;

// Directory utility functions
static std::string trimTrailingSlashes(const std::string &path) {
  if (path.empty()) return path;
  size_t end = path.size();
  while (end > 1 && path[end - 1] == '/') {
    --end;
  }
  return path.substr(0, end);
}

static std::string expandEnvVars(const std::string &input) {
  std::string out;
  out.reserve(input.size());
  for (size_t i = 0; i < input.size(); ++i) {
    if (input[i] == '$') {
      std::string var;
      size_t j = i + 1;
      if (j < input.size() && input[j] == '{') {
        ++j;
        while (j < input.size() && input[j] != '}') {
          var.push_back(input[j++]);
        }
        if (j < input.size() && input[j] == '}') ++j;
      } else {
        while (j < input.size() && (std::isalnum(static_cast<unsigned char>(input[j])) || input[j] == '_')) {
          var.push_back(input[j++]);
        }
      }
      const char *val = std::getenv(var.c_str());
      if (val) out.append(val);
      i = j - 1;
    } else {
      out.push_back(input[i]);
    }
  }
  return out;
}

static std::string expandUser(const std::string &path) {
  if (path.empty()) return path;
  if (path[0] == '~') {
    const char *home = std::getenv("HOME");
    if (!home) home = std::getenv("USERPROFILE");
    if (home) {
      if (path.size() == 1) return std::string(home);
      if (path.size() >= 2 && path[1] == '/') return std::string(home) + path.substr(1);
    }
  }
  return path;
}

static std::string expandPath(const std::string &path) {
  return trimTrailingSlashes(expandUser(expandEnvVars(path)));
}

bool dirExists(const std::string& path) {
  struct stat info;
  if (stat(path.c_str(), &info) != 0) {
    return false;
  }
  return (info.st_mode & S_IFDIR) != 0;
}

bool createDir(const std::string& path) {
  if (mkdir(path.c_str(), 0755) == 0) {
    return true;
  }
  return errno == EEXIST;
}

bool createDirRecursive(const std::string &path) {
  if (path.empty()) return false;
  if (dirExists(path)) return true;

  std::string current;
  // Handle absolute paths
  if (!path.empty() && path[0] == '/') {
    current = "/";
  }
  size_t start = (current == "/") ? 1 : 0;
  for (size_t i = start; i <= path.size(); ++i) {
    if (i == path.size() || path[i] == '/') {
      std::string part = path.substr(start, i - start);
      if (!part.empty()) {
        if (!current.empty() && current.back() != '/') current.push_back('/');
        current += part;
        if (!dirExists(current)) {
          if (!createDir(current)) {
            return false;
          }
        }
      }
      start = i + 1;
    }
  }
  return true;
}

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
  double gravity_magnitude = config["gravity_magnitude"].as<double>();
  bool use_acc_means = config["use_acc_means"].as<bool>();
  double lowpass_hz = config["lowpass_hz"].as<double>();
  bool plot_inspect_static_intervals = config["plot_inspect_static_intervals"] ? config["plot_inspect_static_intervals"].as<bool>() : true;
  std::string save_calib_dir = config["save_calib_dir"].as<std::string>();
  bool calibrate_extrinsic = config["calibrate_extrinsic"] ? config["calibrate_extrinsic"].as<bool>() : true;

  // Expand paths (handle ~ and env vars) and normalize
  bag_file = expandPath(bag_file);
  save_calib_dir = expandPath(save_calib_dir);

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

  // Check if the save directory exists, if not attempt to create it, if fails exit
  if (!save_calib_dir.empty())
  {
    if (!dirExists(save_calib_dir))
    {
      std::cout << "Save directory " << save_calib_dir << " does not exist. Attempting to create it..." << std::endl;
      if (!createDirRecursive(save_calib_dir))
      {
        std::cerr << "Failed to create directory " << save_calib_dir << ". Exiting." << std::endl;
        return -1;
      }
      else
      {
        std::cout << "Successfully created directory " << save_calib_dir << std::endl;
      }
    }
  }

  cout << "Reading the bagfile";
  cout.flush();
  rosbag::Bag bag;
  bag.open(bag_file.c_str(), rosbag::bagmode::Read);

  std::vector<std::string> topics = imu_topics;
  int i = 0;
  rosbag::View view(bag, rosbag::TopicQuery(topics));
  for (const rosbag::MessageInstance& m : view)
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
  mpr_calib.enablePlotting(plot_inspect_static_intervals);
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

  // save_calib_dir was normalized earlier

  // Save per-topic calibration files into save_calib_dir
  for (size_t idx = 0; idx < imu_topics.size(); ++idx)
  {
    std::string base = sanitize_topic(imu_topics[idx]);
    std::string acc_path = save_calib_dir + "/acc_" + base + ".yaml";
    std::string gyr_path = save_calib_dir + "/gyr_" + base + ".yaml";
    std::string extr_path = save_calib_dir + "/extr_" + base + ".yaml";

    acc_calibs[idx].save(acc_path);
    gyro_calibs[idx].save(gyr_path);

    // Fold accelerometer calibration: T_acc = M_after * S_after
    {
      Eigen::Matrix<double, 3, 3> T_acc =
          acc_calibs[idx].getMisalignmentMatrix() * acc_calibs[idx].getScaleMatrix();
      Eigen::Matrix<double, 3, 3> R_acc, M_acc, S_acc;
      imu_tk::foldRotation(T_acc, R_acc, M_acc, S_acc);

      // Verify reconstruction
      Eigen::Matrix<double, 3, 3> T_acc_rec = R_acc * M_acc * S_acc;

      // Append to accelerometer file as comments
      std::ofstream af(acc_path.c_str(), std::ios::app);
      if (af.is_open())
      {
        af << "\n# ---- Folded calibration (T = R * M * S) ----\n";
        imu_tk::appendMatrixComment(af, "T (misalignment * scale):", T_acc);
        imu_tk::appendMatrixComment(af, "R (rotation):", R_acc);
        imu_tk::appendMatrixComment(af, "S (folded scale):", S_acc);
        imu_tk::appendMatrixComment(af, "M (folded misalignment):", M_acc);
        imu_tk::appendMatrixComment(af, "T reconstructed (R * M * S):", T_acc_rec);
        af << "# Recon error: " << (T_acc - T_acc_rec).norm() << "\n";
        
        // Extract and append Euler angles (ZYX convention)
        Eigen::Vector3d euler_angles_acc = R_acc.eulerAngles(0, 1, 2);
        af << "\n# Euler angles (ZYX convention) in radians [yaw, pitch, roll]:\n";
        af << "# " << euler_angles_acc.reverse().transpose() << "\n";
        af << "# Euler angles in degrees [yaw, pitch, roll]:\n";
        af << "# " << (euler_angles_acc.reverse() * 180.0 / M_PI).transpose() << "\n";
      }
    }

    // Fold gyroscope calibration: T_gyr = M_after * S_after
    {
      Eigen::Matrix<double, 3, 3> T_gyr =
          gyro_calibs[idx].getMisalignmentMatrix() * gyro_calibs[idx].getScaleMatrix();
      Eigen::Matrix<double, 3, 3> R_gyr, M_gyr, S_gyr;
      imu_tk::foldRotation(T_gyr, R_gyr, M_gyr, S_gyr);

      // Verify reconstruction
      Eigen::Matrix<double, 3, 3> T_gyr_rec = R_gyr * M_gyr * S_gyr;

      // Append to gyroscope file as comments
      std::ofstream gf(gyr_path.c_str(), std::ios::app);
      if (gf.is_open())
      {
        gf << "\n# ---- Folded calibration (T = R * M * S) ----\n";
        imu_tk::appendMatrixComment(gf, "T (misalignment * scale):", T_gyr);
        imu_tk::appendMatrixComment(gf, "R (rotation):", R_gyr);
        imu_tk::appendMatrixComment(gf, "S (folded scale):", S_gyr);
        imu_tk::appendMatrixComment(gf, "M (folded misalignment):", M_gyr);
        imu_tk::appendMatrixComment(gf, "T reconstructed (R * M * S):", T_gyr_rec);
        gf << "# Recon error: " << (T_gyr - T_gyr_rec).norm() << "\n";
        
        // Extract and append Euler angles (ZYX convention)
        Eigen::Vector3d euler_angles_gyr = R_gyr.eulerAngles(0, 1, 2);
        gf << "\n# Euler angles (ZYX convention) in radians [yaw, pitch, roll]:\n";
        gf << "# " << euler_angles_gyr.reverse().transpose() << "\n";
        gf << "# Euler angles in degrees [yaw, pitch, roll]:\n";
        gf << "# " << (euler_angles_gyr.reverse() * 180.0 / M_PI).transpose() << "\n";
      }
    }

    if (calibrate_extrinsic)
    {
      const auto &extr_all = mpr_calib.getExtrinsicCalib();
      std::ofstream ef(extr_path.c_str());
      if (ef.is_open())
      {
        ef << "extrinsic: [" << extr_all[idx][0] << ", "
           << extr_all[idx][1] << ", "
           << extr_all[idx][2] << "]" << std::endl;

        std::cout << "Saved calibration for topic " << imu_topics[idx] << " to:\n"
                  << "  Acc: " << acc_path << "\n"
                  << "  Gyr: " << gyr_path << std::endl;
      }
      else
      {
        std::cerr << "Failed to open extrinsic output file: " << extr_path << std::endl;
        std::cerr << "Does the folder exist? If not create it manually." << std::endl;
      }
    }
    if (calibrate_extrinsic)
    {
      std::cout << "  Extrinsic: " << extr_path << std::endl;
    }
  }

  return 0;
}
