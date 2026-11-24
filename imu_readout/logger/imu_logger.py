#!/usr/bin/env python3
"""Record the QT Py RP2040 + BNO085 stream (firmware imu_stream 2.x) for ros_imu_calib.

Writes two files:
  data/<name>.raw.csv  every line from the board, prefixed with the host arrival time [us]
                       (A/G raw acc/gyro with sensor timestamps, C/U the sensor's own values)
  data/<name>.csv      t_us,ax,ay,az,gx,gy,gz in m/s^2 and rad/s for csv_calib_extr,
                       made by logger/imu_merge.py (rerun that script to reprocess)

The board waits until the port is opened, so start this script whenever; plugging the board
in before or after is both fine.

Run inside the pixi env:
  pixi run python logger/imu_logger.py [-d SECONDS] [-n NAME] [--no-merge]
"""
import argparse
import datetime as dt
import glob
import os
import sys
import time
from pathlib import Path

import serial

sys.path.insert(0, str(Path(__file__).resolve().parent))
import imu_merge  # noqa: E402

BY_ID_GLOB = "/dev/serial/by-id/usb-Adafruit_QT_Py_RP2040_*-if00"
BAUD = 115200  # ignored by USB CDC
NO_DATA_TIMEOUT_S = 8.0
PROJECT_ROOT = Path(__file__).resolve().parent.parent


def log(msg=""):
    print(msg, file=sys.stderr, flush=True)


def find_port(explicit):
    if explicit:
        return explicit if os.path.exists(explicit) else None
    hits = sorted(glob.glob(BY_ID_GLOB))
    return hits[0] if hits else None


def wait_for_port(explicit):
    port = find_port(explicit)
    if not port:
        log("Waiting for the board ... plug it in.")
    while not port:
        time.sleep(0.05)
        port = find_port(explicit)
    return port


class LineReader:
    """Chunked line reader; every line in a chunk gets the chunk's arrival time."""

    def __init__(self, ser):
        self.ser = ser
        self.buf = b""

    def lines(self):
        chunk = self.ser.read(self.ser.in_waiting or 1)
        t_us = time.time_ns() // 1000
        if not chunk:
            return []
        self.buf += chunk
        *complete, self.buf = self.buf.split(b"\n")
        out = []
        for raw in complete:
            text = raw.decode("ascii", errors="replace").strip()
            if text:
                out.append((t_us, text))
        return out


def is_event(text):
    return len(text) > 2 and text[0] in "AGCU" and text[1] == ","


def connect(explicit_port):
    """Open the port and wait for the first sensor event."""
    while True:
        port = wait_for_port(explicit_port)
        try:
            ser = serial.Serial(port, BAUD, timeout=0.2)
        except PermissionError:
            log(f"Permission denied on {port}. Install the udev rule:\n"
                "  sudo cp logger/70-qtpy-imu.rules /etc/udev/rules.d/ "
                "&& sudo udevadm control --reload-rules && sudo udevadm trigger")
            sys.exit(1)
        except serial.SerialException as e:
            log(f"Could not open {port} yet ({e}); retrying.")
            time.sleep(0.2)
            continue

        log(f"Opened {port}.")
        reader = LineReader(ser)
        deadline = time.monotonic() + NO_DATA_TIMEOUT_S
        status = []
        try:
            while time.monotonic() < deadline:
                for t_us, text in reader.lines():
                    if is_event(text):
                        return ser, reader, port, status, (t_us, text)
                    if text.count(",") == 5 and text.replace(",", "").lstrip("-").isdigit():
                        log("The board sends the old 6-column format: it still runs the original "
                            "firmware. Flash firmware/sensor_read_out_code first.")
                        sys.exit(1)
                    log(f"[board] {text}")
                    status.append(text)
                    if text.startswith("[ERR]"):
                        log("The board reports an error (see LED). Replug it to retry.")
                        sys.exit(1)
        except (serial.SerialException, OSError):
            log("Board disconnected while waiting.")
            ser.close()
            continue
        ser.close()
        log(f"No data within {NO_DATA_TIMEOUT_S:.0f} s. Replug the board and try again.")
        sys.exit(1)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-n", "--name", help="recording name (default: imu_YYYYmmdd_HHMMSS)")
    ap.add_argument("-o", "--out-dir", type=Path, default=PROJECT_ROOT / "data")
    ap.add_argument("-p", "--port", help=f"serial port (default: auto, {BY_ID_GLOB})")
    ap.add_argument("-d", "--duration", type=float, help="stop after S seconds (default: Ctrl-C)")
    ap.add_argument("--no-merge", action="store_true", help="only write the .raw.csv")
    args = ap.parse_args()

    name = args.name or f"imu_{dt.datetime.now():%Y%m%d_%H%M%S}"
    args.out_dir.mkdir(parents=True, exist_ok=True)
    raw_path = args.out_dir / f"{name}.raw.csv"

    ser, reader, port, status, first = connect(args.port)
    log(f"Streaming. Writing {raw_path}" + (f" for {args.duration:.0f} s" if args.duration
                                            else " until Ctrl-C"))

    counts = {k: 0 for k in "AGCU"}
    last_seq = {}
    lost = {"A": 0, "G": 0}
    t_start = time.monotonic()
    t_status = t_start
    with open(raw_path, "w") as f:
        f.write(f"# host_us,type,seq,t_us,x,y,z[,temp]  port: {port}  "
                f"start: {dt.datetime.now().isoformat(timespec='seconds')}\n")
        for s in status:
            f.write(f"# [board] {s}\n")
        pending = [first]
        try:
            while True:
                for t_us, text in pending:
                    if not is_event(text):
                        log(f"\n[board] {text}")
                        f.write(f"# [board] {text}\n")
                        continue
                    f.write(f"{t_us},{text}\n")
                    k = text[0]
                    counts[k] += 1
                    if k in lost:
                        seq = int(text.split(",", 2)[1])
                        if k in last_seq:
                            lost[k] += (seq - last_seq[k] - 1) % 256
                        last_seq[k] = seq

                now = time.monotonic()
                if args.duration and now - t_start >= args.duration:
                    break
                if now - t_status >= 1.0:
                    t_status = now
                    el = now - t_start
                    log_line = (f"\r{el:6.1f} s  acc {counts['A']:7d} ({counts['A'] / el:5.1f} Hz, "
                                f"lost {lost['A']})  gyro {counts['G']:7d} "
                                f"({counts['G'] / el:5.1f} Hz, lost {lost['G']})   ")
                    print(log_line, end="", file=sys.stderr, flush=True)
                pending = reader.lines()
        except KeyboardInterrupt:
            pass
        except (serial.SerialException, OSError) as e:
            log(f"\nBoard disconnected ({e}). Keeping what was recorded.")
        finally:
            ser.close()

    log(f"\nSaved {raw_path}\n")
    if args.no_merge:
        return 0
    return imu_merge.main([str(raw_path)])


if __name__ == "__main__":
    sys.exit(main())
