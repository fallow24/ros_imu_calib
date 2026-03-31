import matplotlib.pyplot as plt
import numpy as np
from imu_simulator.io import bag2pose as bagio

font_increase = 5  # Increase font size by 2 points

# Enable LaTeX rendering for matplotlib
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = 12 +font_increase  # Default font size
plt.rcParams['axes.labelsize'] = 12 +font_increase  # Axis labels
plt.rcParams['axes.titlesize'] = 14 +font_increase  # Title
plt.rcParams['xtick.labelsize'] = 10 +font_increase  # X-axis ticks
plt.rcParams['ytick.labelsize'] = 10 +font_increase  # Y-axis ticks
plt.rcParams['legend.fontsize'] = 11 +font_increase  # Legend

acc, gyr = bagio.load_rosbag_imu("/home/fabi/Documents/Bagfiles/tim_without_lio_bowling_2026-03-27-14-37-10.bag", '/posePub_merged', 1774618640, 1774618648)
acc_raw, gyr_raw = bagio.load_rosbag_imu("/home/fabi/Documents/Bagfiles/tim_without_lio_bowling_2026-03-27-14-37-10.bag", '/posePub_merged_raw', 1774618640, 1774618648)


samples = 10
kernel = np.ones(samples)/samples
acc = np.column_stack([np.convolve(acc[:, i], kernel, mode='valid') for i in range(acc.shape[1])])
acc_raw = np.column_stack([np.convolve(acc_raw[:, i], kernel, mode='valid') for i in range(acc_raw.shape[1])])

acc_time = acc[:, 0] - acc[0, 0]
acc_raw_time = acc_raw[:, 0] - acc_raw[0, 0]
gyr_time = gyr[:, 0] - gyr[0, 0]

# Create two plots above each other, one for the total magnitude of the acceleration and the other for the gyroscope values (x, y, z)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

acc_magnitude = (acc[:, 1]**2 + acc[:, 2]**2 + acc[:, 3]**2)**0.5
acc_raw_magnitude = (acc_raw[:, 1]**2 + acc_raw[:, 2]**2 + acc_raw[:, 3]**2)**0.5
gravity_magnitude = 9.81
ax1.axhline(gravity_magnitude / gravity_magnitude, color='blue', linestyle='--', label='Gravity')
ax1.plot(acc_time, acc_magnitude / gravity_magnitude, label='Compensated', color='green')
ax1.plot(acc_raw_time, acc_raw_magnitude / gravity_magnitude, label='Uncompensated', color='red')
ax1.set_title('Total Magnitude of Acceleration')
ax1.set_ylabel('Acceleration [g]')
ax1.legend(loc='lower left')
ax1.grid()
# Plot the gyroscope values (x, y, z)
ax2.plot(gyr_time, gyr[:, 1], label='X')
ax2.plot(gyr_time, gyr[:, 2], label='Y')
ax2.plot(gyr_time, gyr[:, 3], label='Z')
ax2.set_title('Gyroscope Measurements')
ax2.set_xlabel('Time (s)')
ax2.set_ylabel('Angular Velocity [rad/s]')
ax2.legend(loc='upper left')
ax2.grid()
plt.tight_layout()
print("Saving plot...")
plt.savefig('/home/fabi/catkin_ws/src/ros_imu_calib/pub/img/real_compensation.pgf')
print("Plot saved.")
plt.show()



