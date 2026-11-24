#!/usr/bin/env python3
"""Turn a raw IMU stream (imu_stream 2.x firmware) into the CSV that csv_calib_extr expects.

Input  (<name>.raw.csv, written by imu_logger.py): one line per sensor event,
    <host_us>,A,<seq>,<t_us>,ax,ay,az          raw acc   [counts], sensor clock
    <host_us>,G,<seq>,<t_us>,gx,gy,gz,<temp>   raw gyro  [counts], sensor clock
    <host_us>,C,<seq>,<board_us>,x,y,z         sensor's calibrated acc   [m/s^2]
    <host_us>,U,<seq>,<board_us>,x,y,z         sensor's uncalibrated gyro [rad/s]
Output (<name>.csv):  t_us,ax,ay,az,gx,gy,gz   in m/s^2 and rad/s
    one row per gyro sample, t_us = gyro sample time on the sensor clock,
    acc linearly interpolated to that time (acc runs at its own rate, ~250 Hz).

Unit conversion: the BNO085 does not document the scale of its raw counts, so the raw values
are fitted against the sensor's own C/U values (C = M*A + b). Only the nominal part of M is
applied: per-axis scale plus the axis permutation/sign, no cross-axis terms and no bias, so the
real calibration is still left to imu_tk. Sources, in order:
    --acc-scale/--gyro-scale  >  scale file (logger/imu_scale.json)  >  fit from this recording.
A fit needs motion: rotate the board through several orientations (acc) and turn it around
all three axes (gyro). Save a good fit with --save-scale; later recordings then reuse it.

Usage:
    pixi run python logger/imu_merge.py data/<name>.raw.csv [-o out.csv] [--save-scale]
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

SCALE_FILE = Path(__file__).resolve().parent / "imu_scale.json"
MIN_FIT_STD = {"acc": 2.0, "gyro": 0.3}  # m/s^2, rad/s of the reference per axis
MIN_FIT_R2 = 0.999


def log(msg=""):
    print(msg, file=sys.stderr, flush=True)


def read_raw(path):
    """Return dict type -> list of rows (host_us, seq, t, x, y, z[, temp]) and comment lines."""
    streams = {k: [] for k in "AGCU"}
    comments = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                comments.append(line)
                continue
            p = line.split(",")
            if len(p) < 7 or p[1] not in streams:
                continue
            try:
                host, seq, t = int(p[0]), int(p[2]), int(p[3])
                vals = [float(v) for v in p[4:]]
            except ValueError:
                continue
            streams[p[1]].append([host, seq, t] + vals)
    return {k: np.array(v, dtype=float) if v else np.zeros((0, 7)) for k, v in streams.items()}, comments


def unwrap_u32(t):
    """Unwrap a uint32 microsecond counter (wraps every ~71.6 min)."""
    t = t.astype(np.int64)
    d = np.diff(t)
    d[d < -(1 << 31)] += 1 << 32
    return np.concatenate([[t[0]], t[0] + np.cumsum(d)]) if len(t) else t


def lost_reports(seq):
    """Reports lost according to the 8-bit sequence numbers."""
    jumps = np.diff(seq.astype(int)) % 256
    return int(np.sum(jumps[jumps != 1] - 1)), int(np.sum(jumps != 1))


def pair_by_seq(raw, ref, window=12):
    """Pair each reference event with the raw event of equal seq that arrived closest to it.

    Both come from the same physical sample and arrive within a few lines of each other.
    Returns index arrays into raw and ref.
    """
    ri, fi = [], []
    order_raw = np.argsort(raw[:, 0], kind="stable")  # by host arrival time
    raw_host = raw[order_raw, 0]
    for j, (host, seq) in enumerate(ref[:, :2]):
        k = np.searchsorted(raw_host, host)
        lo, hi = max(0, k - window), min(len(raw_host), k + window)
        cand = order_raw[lo:hi]
        cand = cand[raw[cand, 1] == seq]
        if len(cand):
            best = cand[np.argmin(np.abs(raw[cand, 0] - host))]
            ri.append(best)
            fi.append(j)
    return np.array(ri, dtype=int), np.array(fi, dtype=int)


def fit_nominal(raw_xyz, ref_xyz):
    """Fit ref = M*raw + b, then reduce M to signed permutation * per-axis scale.

    Returns (N, info) where N is the 3x3 nominal matrix and info describes the fit.
    """
    X = np.hstack([raw_xyz, np.ones((len(raw_xyz), 1))])
    sol, *_ = np.linalg.lstsq(X, ref_xyz, rcond=None)
    M, b = sol[:3].T, sol[3]
    pred = X @ sol
    ss_res = np.sum((ref_xyz - pred) ** 2, axis=0)
    ss_tot = np.sum((ref_xyz - ref_xyz.mean(axis=0)) ** 2, axis=0)
    r2 = 1 - ss_res / np.maximum(ss_tot, 1e-12)
    N = np.zeros((3, 3))
    for out_axis in range(3):
        in_axis = int(np.argmax(np.abs(M[out_axis])))
        N[out_axis, in_axis] = M[out_axis, in_axis]
    is_perm = sorted(int(np.argmax(np.abs(r))) for r in M) == [0, 1, 2]
    return N, {"M": M, "bias": b, "r2": r2, "rms": np.sqrt(ss_res / len(ref_xyz)),
               "ref_std": ref_xyz.std(axis=0), "is_perm": is_perm}


def describe(name, N, info=None):
    lines = []
    for o in range(3):
        i = int(np.argmax(np.abs(N[o])))
        lines.append(f"  {name}_{'xyz'[o]} = {N[o, i]:+.6e} * raw_{'xyz'[i]}")
    if info is not None:
        lines.append(f"  fit R^2 {np.round(info['r2'], 5).tolist()}, residual RMS "
                     f"{np.round(info['rms'], 4).tolist()}, reference std "
                     f"{np.round(info['ref_std'], 3).tolist()}")
    return "\n".join(lines)


def fit_quality_ok(kind, info):
    return (info["is_perm"] and np.all(info["ref_std"] >= MIN_FIT_STD[kind])
            and np.all(info["r2"] >= MIN_FIT_R2))


def get_scales(streams, args):
    """Return {'acc': N, 'gyro': N} and the source of each."""
    saved = {}
    if args.scale_file.exists():
        saved = {k: np.array(v) for k, v in json.loads(args.scale_file.read_text()).items()
                 if k in ("acc", "gyro")}
    out, src, fits = {}, {}, {}
    for kind, raw_key, ref_key, cli in (("acc", "A", "C", args.acc_scale),
                                        ("gyro", "G", "U", args.gyro_scale)):
        raw, ref = streams[raw_key], streams[ref_key]
        if len(raw) > 50 and len(ref) > 50:
            ri, fi = pair_by_seq(raw, ref)
            if len(ri) > 50:
                fits[kind] = fit_nominal(raw[ri, 3:6], ref[fi, 3:6])
        if cli is not None:
            out[kind], src[kind] = np.diag([cli] * 3), f"--{kind}-scale"
        elif kind in saved and not args.refit:
            out[kind], src[kind] = saved[kind], str(args.scale_file)
        elif kind in fits and fit_quality_ok(kind, fits[kind][1]):
            out[kind], src[kind] = fits[kind][0], "fit from this recording"
        else:
            out[kind], src[kind] = None, None
    return out, src, fits


def merge(streams, scales):
    A, G = streams["A"], streams["G"]
    ta, tg = unwrap_u32(A[:, 2]), unwrap_u32(G[:, 2])
    # both raw streams use the same sensor clock; make sure the unwrap agrees
    acc = A[:, 3:6] @ scales["acc"].T
    gyr = G[:, 3:6] @ scales["gyro"].T
    oa = np.argsort(ta, kind="stable")
    ta, acc = ta[oa], acc[oa]
    dt_med = np.median(np.diff(ta))
    keep = (tg >= ta[0]) & (tg <= ta[-1])
    k = np.clip(np.searchsorted(ta, tg), 1, len(ta) - 1)
    gap = ta[k] - ta[k - 1]
    keep &= gap <= 3 * dt_med  # do not interpolate across missing acc data
    w = ((tg - ta[k - 1]) / np.maximum(gap, 1))[:, None]
    acc_i = acc[k - 1] * (1 - w) + acc[k] * w
    return tg[keep], acc_i[keep], gyr[keep], int(np.sum(~keep))


def stream_stats(name, s, unit_clock="sensor"):
    if len(s) < 3:
        return f"  {name}: {len(s)} events"
    t = unwrap_u32(s[:, 2])
    dt = np.diff(t)
    lost, gaps = lost_reports(s[:, 1])
    rate = (len(t) - 1) / ((t[-1] - t[0]) * 1e-6)
    return (f"  {name}: {len(s):6d} events, {rate:6.1f} Hz ({unit_clock} clock), dt [us] "
            f"min/median/max {dt.min():.0f}/{np.median(dt):.0f}/{dt.max():.0f}, "
            f"lost {lost} ({gaps} gaps)")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("raw", type=Path, help="<name>.raw.csv from imu_logger.py")
    ap.add_argument("-o", "--out", type=Path, help="output CSV (default: <name>.csv)")
    ap.add_argument("--acc-scale", type=float, help="force one acc scale [m/s^2 per count]")
    ap.add_argument("--gyro-scale", type=float, help="force one gyro scale [rad/s per count]")
    ap.add_argument("--scale-file", type=Path, default=SCALE_FILE)
    ap.add_argument("--refit", action="store_true", help="ignore the scale file, fit from this recording")
    ap.add_argument("--save-scale", action="store_true",
                    help="store this recording's fit in the scale file (only if the fit is good)")
    args = ap.parse_args(argv)

    streams, comments = read_raw(args.raw)
    log(f"{args.raw}:")
    for name, key, clock in (("raw acc  (A)", "A", "sensor"), ("raw gyro (G)", "G", "sensor"),
                             ("ref acc  (C)", "C", "board"), ("ref gyro (U)", "U", "board")):
        log(stream_stats(name, streams[key], clock))
    if len(streams["A"]) < 2 or len(streams["G"]) < 2:
        log("Not enough A/G events to merge.")
        return 1

    scales, src, fits = get_scales(streams, args)
    for kind in ("acc", "gyro"):
        if kind in fits:
            ok = fit_quality_ok(kind, fits[kind][1])
            log(f"{kind} fit from this recording ({'good' if ok else 'not usable: needs more motion'}):")
            log(describe(kind, *fits[kind]))
    if args.save_scale:
        good = {k: fits[k][0].tolist() for k in fits if fit_quality_ok(k, fits[k][1])}
        if good:
            data = json.loads(args.scale_file.read_text()) if args.scale_file.exists() else {}
            data.update(good)
            data["source"] = str(args.raw)
            args.scale_file.write_text(json.dumps(data, indent=2) + "\n")
            log(f"Saved {', '.join(good)} scale to {args.scale_file}")
        else:
            log("No good fit to save.")

    missing = [k for k in ("acc", "gyro") if scales[k] is None]
    if missing:
        log(f"\nNo scale for {', '.join(missing)}: no scale file and this recording has too little "
            "motion for a fit.\nRecord a scale recording (rotate the board through several "
            "orientations and turn it about all axes), run this script on it with --save-scale, "
            "then rerun on this file.")
        return 2
    for kind in ("acc", "gyro"):
        log(f"{kind} scale from {src[kind]}")

    t, acc, gyr, dropped = merge(streams, scales)
    out = args.out or args.raw.with_name(args.raw.name.replace(".raw.csv", ".csv"))
    if out == args.raw:
        out = args.raw.with_suffix(".merged.csv")
    with open(out, "w") as f:
        f.write(f"# t_us,ax,ay,az,gx,gy,gz  units: m/s^2, rad/s  t: sensor clock  "
                f"source: {args.raw.name}\n")
        for ti, a, g in zip(t, acc, gyr):
            f.write(f"{ti},{a[0]:.6f},{a[1]:.6f},{a[2]:.6f},{g[0]:.7f},{g[1]:.7f},{g[2]:.7f}\n")
    rate = (len(t) - 1) / ((t[-1] - t[0]) * 1e-6) if len(t) > 1 else 0
    log(f"\nWrote {len(t)} rows ({rate:.1f} Hz mean, set imu_frequency_hz to {round(rate)}) "
        f"to {out}" + (f"; skipped {dropped} gyro samples without acc around them" if dropped else ""))
    first = t <= t[0] + 2_000_000
    g_mag = np.linalg.norm(acc[first], axis=1)
    log(f"|acc| over the first 2 s: mean {g_mag.mean():.3f} m/s^2, std {g_mag.std():.3f} "
        "(9.81 if the board was still)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
