# IMU readout: QT Py RP2040 + BNO085

Firmware and host logger that produce the `t_us,ax,ay,az,gx,gy,gz` CSV (m/s², rad/s) read by
`apps/csv_calib_extr.cpp`.

- Board: Adafruit QT Py RP2040 + BNO085 on Wire1 (STEMMA QT), I2C 0x4A, 400 kHz
- `firmware/sensor_read_out_code/`: firmware `imu_stream 2.0`, one ASCII line per sensor event:
  - `A,seq,t_us,x,y,z`: raw acc in counts, sensor clock
  - `G,seq,t_us,x,y,z,temp`: raw gyro in counts, sensor clock
  - `C,...`: calibrated acc in m/s², board clock
  - `U,...`: uncalibrated gyro in rad/s, board clock

  C/U must stay enabled because the BNO08x only emits raw reports while the processed ones run.
  Raw acc runs at about 250 Hz and raw gyro at 200 Hz, with no lost samples. The board waits
  for the host to open the port; there's no timeout.
- `firmware/backup/`: the original sketch. Full flash dumps of the original board state are kept
  outside the repo (24 MB):
  - `flashed_2026-11-11.bin`, sha256 `6e7154e1…dca45231`
  - `flashed_2026-11-11.uf2`, sha256 `360bc1b6…c8844b84`
- LED:
  - orange/black: waiting for host
  - green/white: streaming
  - red+purple: sensor not found
  - red+yellow: setup failed
  - red+blue: sensor stopped responding

## Changes vs. the original firmware
- Events come from an own sh2 callback. `Adafruit_BNO08x::getSensorEvent()` keeps only the last
  event per `sh2_service()` call, which had dropped about 22% of samples.
- Every sample carries the sensor timestamp and a sequence number, so lost samples are visible.
- There's no 5 s startup delay or deadline. Closing and reopening the port works without replugging.

## Build / flash
Uses arduino-pico core `rp2040:rp2040` 6.3.2 and the libraries Adafruit BNO08x 1.2.7 and
Adafruit NeoPixel 1.15.5.
```
arduino-cli compile --fqbn rp2040:rp2040:adafruit_qtpy firmware/sensor_read_out_code
arduino-cli upload  --fqbn rp2040:rp2040:adafruit_qtpy -p /dev/ttyACM0 firmware/sensor_read_out_code
```

## Record
```
sudo cp logger/70-qtpy-imu.rules /etc/udev/rules.d/ && sudo udevadm control --reload-rules && sudo udevadm trigger
python logger/imu_logger.py -d 60 -n NAME     # needs pyserial, numpy
```
This writes two files to `data/` next to `logger/`:
- `NAME.raw.csv`: every board line plus the host arrival time.
- `NAME.csv`: the merged SI CSV from `logger/imu_merge.py`. Gyro timestamps are used and acc is
  interpolated.

The raw-to-SI scale comes from `logger/imu_scale.json`. It's fitted against C/U from a recording
with motion: `python logger/imu_merge.py data/NAME.raw.csv --save-scale`.
