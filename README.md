# IMU-TK

This is a ROS wrapper for imu_tk, learn details from original [README.md](./README_original.md) please.

## Dependencies 

This program depends on [gnuplot](http://www.gnuplot.info/), [yaml-cpp](https://github.com/jbeder/yaml-cpp), [Eigen](https://libeigen.gitlab.io/docs/), and [Ceres](http://ceres-solver.org/index.html).

On Linux (Ubuntu/Debian-based) you can try:
```Bash
sudo apt install libeigen3-dev libyaml-cpp-dev gnuplot gnuplot-qt
```

Install Ceres by doing the following:

```Bash
git clone https://ceres-solver.googlesource.com/ceres-solver
cd ceres-solver
mkdir build
cd build
cmake ..
make
sudo make install 
```

If you have trouble installing Ceres' dependencies refer to [their website](http://ceres-solver.org/installation.html#).  

## Build

Place this package inside your catkin workspace (e.g. `~/catkin_ws/src`), then build it through `catkin_make`.

```Bash
cd ~/catkin_ws
catkin_make
rospack profile
source devel/setup.bash
```

## Test

### Collect IMU Data

Record a bag file with IMU topic (e.g. `rosbag record -O imu /imu`).

Procedure:

1. Left the IMU static for 50 seconds.
2. Rotate the IMU and then lay it in a different attitude.
3. Wait for at least 1 seconds.
4. Have you rotated the IMU 36 ~ 50 times? If not, go back to step 2.
5. Done.

### Calibration

Run `imu_calib_node` to get calibration result(e.g. misalignment, scale and bias).

`imu_calib_node` usage: 

> rosrun imu_tk imu_calib_node [CONFIG]

```Bash
rosrun imu_tk imu_calib_node ~/catkin_ws/src/ros_imu_calib/config/imu_calib.yaml
```

Follow the output of the program on your terminal.
You should be prompted to inspect the static intervals in a plot window.
Refer to the [original paper](http://dx.doi.org/10.1109/ICRA.2014.6907297)) for more details.
This is the corresponding citation:

```
@INPROCEEDINGS{6907297,
  author={Tedaldi, David and Pretto, Alberto and Menegatti, Emanuele},
  booktitle={2014 IEEE International Conference on Robotics and Automation (ICRA)}, 
  title={A robust and easy to implement method for IMU calibration without external equipments}, 
  year={2014},
  volume={},
  number={},
  pages={3042-3049},
  keywords={Calibration;Accelerometers;Gyroscopes;Sensors;Gravity;Accuracy;Vectors},
  doi={10.1109/ICRA.2014.6907297}
}
```

### Correction

Given a raw sensor reading X (e.g., the acceleration ), the calibrated "unbiased" reading X' is obtained.

```Txt
Misalignment matrix:

    [    1     -mis_yz   mis_zy  ]
T = [  mis_xz     1     -mis_zx  ]
    [ -mis_xy   mis_yx     1     ]

Scale matrix:

    [  s_x      0        0  ]
K = [   0      s_y       0  ]
    [   0       0       s_z ]

Bias vector:

    [ b_x ]
B = [ b_y ]
    [ b_z ]

X' = T*K*(X - B)
```

## Future work

I do not plan to maintain this forever.
However, I want to spread this beatiful software and make it easy to use. 
Plans going forward include:

- Add a branch for ROS2
