#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <yaml-cpp/yaml.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

#include "imu_tk/io_utils.h"
#include "imu_tk/calibration.h"
#include "imu_tk/filters.h"
#include "imu_tk/integration.h"
#include "imu_tk/visualization.h"

using namespace std;
using namespace imu_tk;
using namespace Eigen;

// Directory utility functions
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

bool checkPathAndCreateDir(const std::string& filepath) {
  size_t found = filepath.find_last_of("/\\");
  if (found != string::npos) {
    string dir = filepath.substr(0, found);
    if (!dirExists(dir)) {
      cout << "Creating directory: " << dir << endl;
      if (!createDir(dir)) {
        cerr << "Failed to create directory: " << dir << endl;
        return false;
      }
    }
  }
  return true;
}

int readCSV(string &csv_file, vector< TriadData > &acc_data, vector< TriadData > &gyro_data) 
{
  ifstream csv(csv_file.c_str());
  if (!csv.is_open()) {
    cerr << endl << "Failed to open csv file: " << csv_file << endl;
    return -1;
  }

  // Expected columns: timestamp [us], acc x, y, z, gyro x, y, z
  int i = 0;
  bool have_t0 = false;
  long long t0_us = 0;
  string line;
  while (getline(csv, line)) {
    if (!line.empty() && line[line.size() - 1] == '\r')
      line.erase(line.size() - 1);
    if (line.empty() || line[0] == '#')
      continue;

    vector<string> fields;
    stringstream ss(line);
    string field;
    while (getline(ss, field, ','))
      fields.push_back(field);
    if (fields.size() < 7)
      continue;

    long long t_us;
    double v[6];
    char* end;
    errno = 0;
    t_us = strtoll(fields[0].c_str(), &end, 10);
    bool ok = (end != fields[0].c_str()) && errno == 0;
    for (int k = 0; ok && k < 6; ++k) {
      v[k] = strtod(fields[k + 1].c_str(), &end);
      ok = (end != fields[k + 1].c_str());
    }
    if (!ok)
      continue; // header or malformed line

    // Offset by the first timestamp to keep double precision for microsecond stamps
    if (!have_t0) {
      t0_us = t_us;
      have_t0 = true;
    }
    double t = (t_us - t0_us) * 1e-6;

    if (i++ % 1000 == 0) {
      cout << "."; // Notify user that something's happening
      cout.flush();
    }
    acc_data.push_back(TriadData(t, v[0], v[1], v[2]));
    gyro_data.push_back(TriadData(t, v[3], v[4], v[5]));
  }
  cout << endl << "Read " << i << " samples" << endl;
  csv.close();
  if (i == 0) {
    cerr << "No valid samples found in " << csv_file << endl;
    return -1;
  }
  return 1;
}

int main(int argc, char** argv)
{
  if (argc < 2) {
    cout << "Usage: rosrun ros_imu_calib csv_calib_extr <config.yaml>" << endl;
    cout << "Have a look at \"catkin_ws/src/imu_tk/config/imu_calib_extr.yaml\"" << endl;
    return -1;
  }

  // Load YAML config
  YAML::Node config = YAML::LoadFile(argv[1]);

  // program parameters
  string csv_file      = config["csv_file"].as<string>();
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
  string save_extr_calib        = config["save_extr_calib"].as<string>();
  std::vector<double> extr_init = config["extr_initial_guess"].as<std::vector<double> >();
  double lowpass_hz             = config["lowpass_hz"].as<double>();

  // Check if save directories exist
  checkPathAndCreateDir(save_acc_calib);
  checkPathAndCreateDir(save_gyr_calib);
  checkPathAndCreateDir(save_extr_calib);

  vector< TriadData > acc_data, gyro_data;
  
  cout << "Reading the csv file";
  cout.flush();
  int have_data = readCSV(csv_file, acc_data, gyro_data);
  if (!have_data) return 0;

  // Calculate the mean frequency of the IMU data using the timestamps of the accelerometer data
  double total_time = acc_data.back().timestamp() - acc_data.front().timestamp();
  cout << "Total time of the IMU data: " << total_time << " seconds" << endl;
  double mean_frequency = (acc_data.size() - 1) / total_time;
  cout << "Mean frequency of the IMU data: " << mean_frequency << " Hz" << endl;
  if (std::abs(mean_frequency - imu_frequency_hz) > 0.1 * imu_frequency_hz) {
    std::cerr << "Warning: The mean frequency of the IMU data (" << mean_frequency 
              << " Hz) differs significantly from the expected frequency (" 
              << imu_frequency_hz << " Hz)." << std::endl;
    return -1;
  }
  
  std::cout << std::endl << "Now performing calibration" << std::endl << std::endl;
  CalibratedTriad init_acc_calib, init_gyro_calib;
  init_acc_calib.setBias( Vector3d(0, 0, 0) );
  init_gyro_calib.setScale( Vector3d(1, 1, 1) );
  MultiPosCalibration mp_calib;  
  mp_calib.setInitExtrinsicCalibration(extr_init);
  mp_calib.setInitStaticIntervalDuration(initial_static_interval);
  mp_calib.setIntarvalsNumSamples(num_samples_per_interval);
  mp_calib.setInitAccCalibration( init_acc_calib );
  mp_calib.setInitGyroCalibration( init_gyro_calib );
  mp_calib.setGravityMagnitude(gravity_magnitude);
  mp_calib.enableVerboseOutput(true);
  mp_calib.enablePlotting(plot_inspect_static_intervals);
  mp_calib.enableAccUseMeans(use_acc_means);
  mp_calib.setGyroDataPeriod(1.0 / imu_frequency_hz);
  mp_calib.setLowpassFreq(lowpass_hz);
  mp_calib.calibrateAccGyroExtrinsic(acc_data, gyro_data, 
    min_nr_recorded_intervals, variance_win_size, 
    th_variance_min, th_variance_max);
  mp_calib.getAccCalib().save(save_acc_calib);
  mp_calib.getGyroCalib().save(save_gyr_calib);
  mp_calib.saveExtrinsicCalib(save_extr_calib);
  return 0;
}

