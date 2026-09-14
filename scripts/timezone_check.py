"""Empirically determine whether the raw TIME column in the .sec files is
UTC or Taiwan local time (UTC+8). The IAGA-2002 header has no timezone
field, so this is inferred from the data itself using two independent
diurnal signals for the group's nearest available station:

  Signal A: Sq (solar-quiet) diurnal curve of the station's primary
  magnitude signal (H = sqrt(X^2+Y^2) for a vector station, or F itself for
  a scalar-only station). The Sq current system produces a broad extremum
  near local solar noon. Taiwan sits almost exactly on the 120 deg E UTC+8
  standard meridian, so local solar noon is at essentially clock-noon in
  Taiwan local time -- an 8-hour displacement from UTC clock-noon, which is
  large compared to the width of the Sq feature, giving an unambiguous
  discriminator.

  Signal B: diurnal cycle of high-frequency noise power (variance of the
  1-second first difference of the same primary channel). Human/industrial
  activity (power grid, rail, traffic) produces a sharp step-up roughly
  07:00-19:00 in TRUE local clock time and is quiet 01:00-05:00 -- a
  sharper, independent cross-check.

Both signals are computed over a sample of the quietest available days
(lowest daily std, per data/interim/<group>/daily_features.csv) for the
chosen station, excluding any day with excessive missing data
(common.auto_outage_dates).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from common import auto_outage_dates, load_group_config, resolve_day_ref  # noqa: E402
from parser import parse_day_file  # noqa: E402

N_QUIET_DAYS = 18


def _pick_station(cfg) -> str:
    """Prefer the nearest vector station (H = sqrt(X^2+Y^2) is the cleaner
    Sq/noise signal); fall back to the nearest scalar-only station (F) for
    groups with no usable vector pool (e.g. G1)."""
    if cfg.xyz_pool.near:
        return cfg.xyz_pool.near[0]
    if cfg.f_pool.near:
        return cfg.f_pool.near[0]
    raise RuntimeError(f"{cfg.group_id}: no station available at all for timezone check")


def _signals(df: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
    """Return (sq_signal, noise_base_signal) for whichever columns this
    file's header actually reported (see parser.parse_day_file)."""
    if "F" in df.columns:
        return df["F"], df["F"]
    h = np.sqrt(df["X"] ** 2 + df["Y"] ** 2)
    return h, df["X"]


def pick_quiet_days(cfg, station: str) -> list[str]:
    std_col = "F_std" if station in cfg.f_pool.all_stations else "H_std"
    spike_col = "n_spikes"
    daily = pd.read_csv(cfg.interim_dir / "daily_features.csv", dtype={"date": str})
    sub = daily[(daily.station == station) & (daily[spike_col] == 0) & (daily.pct_missing == 0)]
    outage_dates = auto_outage_dates(daily).get(station, set())
    sub = sub[~sub.date.isin(outage_dates)]
    # spread across the whole period: take the lowest-std day from each
    # consecutive ~1-week bucket, up to N_QUIET_DAYS total
    sub = sub.sort_values("date").reset_index(drop=True)
    if len(sub) == 0:
        return []
    sub["bucket"] = np.arange(len(sub)) // max(1, len(sub) // N_QUIET_DAYS)
    picks = sub.loc[sub.groupby("bucket")[std_col].idxmin(), "date"].tolist()
    return sorted(picks)[:N_QUIET_DAYS]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--group", required=True)
    args = ap.parse_args()
    cfg = load_group_config(args.group)
    station = _pick_station(cfg)

    quiet_days = pick_quiet_days(cfg, station)
    print(f"[{args.group}] station={station} quiet days used ({len(quiet_days)}): {quiet_days}", file=sys.stderr)

    sq_by_hour = {h: [] for h in range(24)}
    noise_by_hour = {h: [] for h in range(24)}

    for date_str in quiet_days:
        ref = resolve_day_ref(cfg.gdms_dir, station, date_str)
        df = parse_day_file(ref, station)
        sq_signal, noise_base = _signals(df)
        noise = noise_base.diff().abs()  # 1-second first-difference magnitude
        hour = df.index.hour
        for hr in range(24):
            mask = hour == hr
            sq_by_hour[hr].append(float(sq_signal[mask].mean()))
            noise_by_hour[hr].append(float((noise[mask] ** 2).mean()))  # variance-like proxy

    noise_curve = {hr: float(np.mean(vals)) for hr, vals in noise_by_hour.items()}

    # remove each day's own mean before averaging so day-to-day baseline
    # drift doesn't swamp the diurnal shape
    sq_by_hour2 = {h: [] for h in range(24)}
    for date_str in quiet_days:
        ref = resolve_day_ref(cfg.gdms_dir, station, date_str)
        df = parse_day_file(ref, station)
        sq_signal, _ = _signals(df)
        dev = sq_signal - sq_signal.mean()
        hour = df.index.hour
        for hr in range(24):
            sq_by_hour2[hr].append(float(dev[hour == hr].mean()))
    sq_curve_dev = {hr: float(np.mean(vals)) for hr, vals in sq_by_hour2.items()}

    sq_extreme_hour = max(sq_curve_dev, key=lambda hr: abs(sq_curve_dev[hr])) if quiet_days else None
    noise_peak_hours = sorted(noise_curve, key=lambda hr: -noise_curve[hr])[:6] if quiet_days else []
    noise_quiet_hours = sorted(noise_curve, key=lambda hr: noise_curve[hr])[:6] if quiet_days else []

    # Decision rule: if noise is high in raw-clock hours ~7-19 and low ~1-5,
    # raw clock == local Taiwan time (UTC+8). If that pattern instead sits
    # 8h earlier (peak ~23-11, quiet ~17-21), raw clock == UTC.
    workday_hours = set(range(7, 20))
    night_hours = set(range(1, 6))
    workday_score = sum(1 for hr in noise_peak_hours if hr in workday_hours)
    night_score = sum(1 for hr in noise_quiet_hours if hr in night_hours)

    if not quiet_days:
        conclusion, confidence = "UNCERTAIN", "low"
    elif workday_score >= 4 and night_score >= 4:
        conclusion = "LOCAL (UTC+8)"
        confidence = "high"
    else:
        # check shifted-by-8 hypothesis (UTC labeling)
        workday_hours_utc = {(hr - 8) % 24 for hr in workday_hours}
        night_hours_utc = {(hr - 8) % 24 for hr in night_hours}
        workday_score_utc = sum(1 for hr in noise_peak_hours if hr in workday_hours_utc)
        night_score_utc = sum(1 for hr in noise_quiet_hours if hr in night_hours_utc)
        if workday_score_utc >= 4 and night_score_utc >= 4:
            conclusion = "UTC"
            confidence = "high"
        else:
            conclusion = "UNCERTAIN"
            confidence = "low"

    result = {
        "group": args.group,
        "station": station,
        "quiet_days_used": quiet_days,
        "sq_curve_H_deviation_by_raw_hour": sq_curve_dev,
        "noise_power_by_raw_hour": noise_curve,
        "sq_extreme_raw_hour": sq_extreme_hour,
        "noise_peak_raw_hours": noise_peak_hours,
        "noise_quiet_raw_hours": noise_quiet_hours,
        "workday_score_if_local": workday_score,
        "night_score_if_local": night_score,
        "conclusion_timezone": conclusion,
        "confidence": confidence,
        "note": (
            "conclusion is 'LOCAL (UTC+8)' if noise peaks in raw-clock 07-19 and "
            "is quiet in raw-clock 01-05 (matching real Taiwan work-hours); "
            "'UTC' if that pattern is shifted 8h earlier in the raw clock."
        ),
    }

    out_path = cfg.interim_dir / "timezone_check.json"
    out_path.write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))
    print(f"\nwrote {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
