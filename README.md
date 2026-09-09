# ROS-IMU-CALIB

This originated as a ROS wrapper for imu_tk, learn details from original [README.md](./README_original.md) please.

The program is able to perform intrinsic calibration of an IMU (accelerometer and gyroscope) without external equipment. 

## Updates

- I have added extrinsic calibration with respect to the robot base (the offset vector from the base to the IMU). When performing extrinsic calibration (when collecting the dataset) the operator must ensure that only rotations around the robot base are performed.
  
- I have added multi-IMU calibration (align multiple IMUs to one reference IMU), such that all the IMUs share the same coordinate system as the given reference IMU. For multi-IMU, the program is able to calibrate intrinsic and (optionally) extrinsic parameters.

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

If you desire extrinsic calibration with respect to the robot base, make sure that ONLY rotations around the robot base happen in the bagfile. Additional translation of the IMU will affect the extrinsic results! 

### Calibration

Run either `bag_calib`, `bag_calib_extr`, or `bag_calib_multi` to get calibration result (e.g. misalignment, scale and bias).

`bag_calib` usage: 

> rosrun ros_imu_calib bag_calib [CONFIG]

`bag_calib_extr` usage:

> rosrun ros_imu_calib bag_calib_extr [CONFIG]

`bag_calib_multi` usage:

> rosrun ros_imu_calib bag_calib_multi [CONFIG]

For example:

```Bash
rosrun ros_imu_calib bag_calib ~/catkin_ws/src/ros_imu_calib/config/imu_calib.yaml
```

Please have a look at the example parameters given in the YAML config file.
**Make sure to update the paths accordinly!**

You might use this [example bagfile](https://drive.google.com/file/d/1wRyGn3pAfVio9VZlECjIbqagfqa8B2HJ/view?usp=sharing).


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

The extrinsic calibration is described in detail in [this paper](https://doi.org/10.48550/arXiv.2607.25784).
For now you can cite the preprint from arxiv:

```
@misc{arzberger2026motionaccelerationcalibrationcompensationimus,
      title={Motion-Acceleration Calibration and Compensation in IMUs without External Equipment for Attitude Estimation Filters}, 
      author={Fabian Arzberger and Andreas Nüchter},
      year={2026},
      eprint={2607.25784},
      archivePrefix={arXiv},
      primaryClass={cs.RO},
      url={https://arxiv.org/abs/2607.25784}, 
}
```

### Intrinsic Correction

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

**NOTE:** In the multi-IMU case, the misalignment matrix might not have diagonal entries of 1.
This is intended to represent the axes misalignment to the reference sensor instead of to its own axis.
When using multi-IMU calibration, the resulting yaml files have comments with the decomposed rotation, misalignment, and scale.

### Extrinsic Correction

Given the extrinsic offset vector R one can correct the (calibrated) accelerometer reading A using the (calibrated) gyoroscope reading W (and its derivative dW)

```Txt
A' = A - (dW x R) - (W x W x R) 
```

where `x` denotes the 3d vector cross product

## Future work

I do not plan to maintain this forever.
However, I want to spread this beatiful software and make it easy to use. 
Plans going forward include:

- Add a branch for ROS2
