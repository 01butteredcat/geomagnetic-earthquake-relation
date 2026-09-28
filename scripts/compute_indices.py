"""Basic statistical screening: build a candidate 'local anomaly index' that
attempts to isolate a near-epicenter-only geomagnetic signal from the
common-mode (Sq/storm-driven) variation shared by the whole network.

Runs on whichever of the group's two station pools (see common.py's
GroupConfig) is sufficient: the XYZ pool gives H and Z fields (vector data);
the F pool (scalar-only stations) gives a single F field -- there is no Z/D
analog for scalar data, that's a structural limitation of total-field-only
instruments, not a bug. A group may run one, both, or (if neither pool has
enough stations) no fields at all.

Night window: raw UTC hour in {17,18,19} == Taiwan local time {01,02,03},
confirmed quiet by timezone_check.py for the 2024/G10 data (assumed to hold
across groups, same station/data convention). We deliberately do NOT relabel
this onto a shifted local calendar date -- everything stays keyed by the UTC
file date the samples were read from, which is also how daily_features.csv
and storm_days.csv are keyed.

Method (plan Section 3, generalized to any field/station-pool):
  1. night_mean/std per station/day from the 1-minute series.
  2. MAD z-score of night_mean vs. a trailing TRAILING_WINDOW_DAYS-day
     median (built only from clean -- non-storm, non-outage -- days).
  3. near_index = median dev_z across the pool's near stations; far_index =
     median dev_z across the pool's far stations.
  4. Theil-Sen robust regression near_index ~ far_index on clean days only;
     residual = local_anomaly_index (the candidate local-only signal).
  5. Flag candidate days where |local_anomaly_index| exceeds ~2.5x the MAD
     of its own clean-day distribution.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import theilslopes

sys.path.insert(0, str(Path(__file__).parent))
from common import auto_outage_dates, load_group_config  # noqa: E402

NIGHT_HOURS_UTC = {17, 18, 19}  # == local 01:00-03:59
# A station-night with fewer valid minutes than this (out of 180) contributes
# no index value that day. Before this, a single valid minute was enough, so a
# day whose whole-day pct_missing exceeded OUTAGE_PCT_MISSING_THRESHOLD could
# still be flagged as a candidate from a near-empty night window (e.g. G18
# 2016-02-29: kmn had 68/180 night minutes and hcn 0, and that thin far_index
# produced a candidate). 90 (50% of the window) is a judgment call: 30-60
# leaves that G18 day in, 90-150 all give the same set of removed candidates.
MIN_NIGHT_MINUTES = 90
# 7 days was too short for the G10/2024 case: the 2024-03-21..27 storm +
# 2-day recovery consumed the entire trailing window for dates through early
# April, leaving the most important pre-quake days (03-26..04-04) with no
# baseline at all. 21 days reaches back past that storm into a clean stretch.
# Kept as one dataset-wide constant (not per-group) since it's a property of
# "how long a storm+recovery can plausibly run", not of any one group's data.
TRAILING_WINDOW_DAYS = 21
MIN_CLEAN_POINTS = 5
CANDIDATE_Z_THRESHOLD = 2.5


def night_features_for_station(cfg, station: str, channel: str) -> pd.DataFrame:
    path = cfg.interim_dir / f"minute_series_{station}.parquet"
    df = pd.read_parquet(path)
    night = df[df.index.hour.isin(NIGHT_HOURS_UTC)]
    date_key = night.index.strftime("%Y%m%d")
    g = night.groupby(date_key)
    out = pd.DataFrame(
        {
            f"night_mean_{channel}": g[channel].mean(),
            f"night_std_{channel}": g[channel].std(),
            "night_n": g[channel].count(),
        }
    )
    out.index.name = "date"
    out["station"] = station
    thin = out["night_n"] < MIN_NIGHT_MINUTES
    out.loc[thin, [f"night_mean_{channel}", f"night_std_{channel}"]] = np.nan
    return out.reset_index()


def mad_zscore(series: pd.Series, clean_mask: pd.Series) -> pd.Series:
    """Rolling trailing-window MAD z-score, baseline built from clean days only.
    The window is the TRAILING_WINDOW_DAYS calendar days before each day, not the
    previous TRAILING_WINDOW_DAYS rows: a folder with missing day files (G6_G7_G8
    lacks 2021-12-30 and 2022-01-01) would otherwise reach further back."""
    z = pd.Series(index=series.index, dtype="float64")
    vals = series.values
    clean = clean_mask.values
    days = pd.to_datetime(series.index, format="%Y%m%d")
    starts = np.searchsorted(days, days - pd.Timedelta(days=TRAILING_WINDOW_DAYS))
    for i in range(len(series)):
        lo = int(starts[i])
        window_vals = vals[lo:i][clean[lo:i]] if i > lo else np.array([])
        window_vals = window_vals[~np.isnan(window_vals)]
        if len(window_vals) < MIN_CLEAN_POINTS:
            z.iloc[i] = np.nan
            continue
        med = np.median(window_vals)
        mad = np.median(np.abs(window_vals - med)) * 1.4826
        z.iloc[i] = (vals[i] - med) / mad if mad > 1e-9 else np.nan
    return z


def build_field_index(cfg, field: str, near_stations, far_stations, storm_dates, outage_by_station):
    all_stations = list(near_stations) + list(far_stations)
    frames = [night_features_for_station(cfg, s, field) for s in all_stations]
    night = pd.concat(frames, ignore_index=True)

    all_dates = sorted(night["date"].unique())
    pivot_mean = night.pivot(index="date", columns="station", values=f"night_mean_{field}").reindex(all_dates)

    dev_z = pd.DataFrame(index=all_dates)
    for station in all_stations:
        clean_mask = pd.Series(
            [d not in outage_by_station.get(station, set()) and d not in storm_dates for d in all_dates],
            index=all_dates,
        )
        dev_z[station] = mad_zscore(pivot_mean[station], clean_mask)

    near_index = dev_z[list(near_stations)].median(axis=1)
    far_index = dev_z[list(far_stations)].median(axis=1)

    # A day counts as "clean" (eligible for the regression fit / candidate
    # flagging) if it's not a storm day AND near_index/far_index actually
    # computed a value -- i.e. the station-level median had at least one
    # non-outage contributor, not "every single pool station must
    # individually be outage-free". The stricter all-or-nothing version was
    # fine for G10 (no near/far station there ever had an extended full
    # outage), but a group with e.g. a 3-week total dropout on one near
    # station (seen in G11's zbn right before its anchor event) would
    # otherwise have zero clean days for the whole dropout window even
    # though the median-of-3 near_index remained perfectly computable from
    # the other two stations for that entire stretch.
    clean_overall = pd.Series(
        [d not in storm_dates for d in all_dates], index=all_dates
    ) & near_index.notna() & far_index.notna()
    fit_mask = clean_overall  # already requires near/far to be non-NaN, see clean_overall above
    if fit_mask.sum() >= 10:
        slope, intercept, _, _ = theilslopes(near_index[fit_mask], far_index[fit_mask])
    else:
        slope, intercept = 0.0, 0.0

    predicted = slope * far_index + intercept
    local_anomaly_index = near_index - predicted

    clean_resid = local_anomaly_index[clean_overall.reindex(local_anomaly_index.index, fill_value=False)]
    clean_resid = clean_resid.dropna()
    resid_mad = np.median(np.abs(clean_resid - np.median(clean_resid))) * 1.4826 if len(clean_resid) >= 5 else np.nan

    result = pd.DataFrame(
        {
            "date": all_dates,
            f"near_index_{field}": near_index.values,
            f"far_index_{field}": far_index.values,
            f"predicted_near_{field}": predicted.values,
            f"local_anomaly_index_{field}": local_anomaly_index.values,
        }
    )
    result["is_clean_day"] = clean_overall.values
    if resid_mad and resid_mad > 0:
        result[f"candidate_flag_{field}"] = (
            (result[f"local_anomaly_index_{field}"].abs() > CANDIDATE_Z_THRESHOLD * resid_mad)
            & result["is_clean_day"]
        )
    else:
        result[f"candidate_flag_{field}"] = False

    fit_info = {
        "field": field,
        "near_stations": list(near_stations),
        "far_stations": list(far_stations),
        "slope": float(slope),
        "intercept": float(intercept),
        "n_fit_days": int(fit_mask.sum()),
        "residual_mad": float(resid_mad) if resid_mad == resid_mad else None,
    }
    return result, fit_info


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--group", required=True)
    args = ap.parse_args()
    cfg = load_group_config(args.group)

    storm_days_df = pd.read_csv(cfg.interim_dir / "storm_days.csv", dtype={"date": str})
    storm_dates = set(storm_days_df["date"])
    daily = pd.read_csv(cfg.interim_dir / "daily_features.csv", dtype={"date": str})
    outage_by_station = auto_outage_dates(daily)

    merged = None
    fits: dict = {}
    candidates: dict[str, list[str]] = {}
    methods_run: list[str] = []

    if cfg.xyz_pool.sufficient:
        result_h, fit_h = build_field_index(cfg, "H", cfg.xyz_pool.near, cfg.xyz_pool.far, storm_dates, outage_by_station)
        result_z, fit_z = build_field_index(cfg, "Z", cfg.xyz_pool.near, cfg.xyz_pool.far, storm_dates, outage_by_station)
        merged = result_h.merge(result_z.drop(columns=["is_clean_day"]), on="date", how="outer")
        fits["H"], fits["Z"] = fit_h, fit_z
        candidates["H"] = merged[merged["candidate_flag_H"] == True]["date"].tolist()  # noqa: E712
        candidates["Z"] = merged[merged["candidate_flag_Z"] == True]["date"].tolist()  # noqa: E712
        methods_run += ["H", "Z"]
    else:
        print(f"[{args.group}] XYZ pool insufficient ({len(cfg.xyz_pool.all_stations)} stations) -- skipping H/Z screening", file=sys.stderr)

    if cfg.f_pool.sufficient:
        result_f, fit_f = build_field_index(cfg, "F", cfg.f_pool.near, cfg.f_pool.far, storm_dates, outage_by_station)
        fits["F"] = fit_f
        candidates["F"] = result_f[result_f["candidate_flag_F"] == True]["date"].tolist()  # noqa: E712
        methods_run.append("F")
        if merged is None:
            merged = result_f
        else:
            merged = merged.merge(
                result_f.rename(columns={"is_clean_day": "is_clean_day_F"}), on="date", how="outer"
            )
    else:
        print(f"[{args.group}] F pool insufficient ({len(cfg.f_pool.all_stations)} stations) -- skipping F screening", file=sys.stderr)

    if merged is None:
        print(f"[{args.group}] neither station pool sufficient -- no basic screening possible for this group", file=sys.stderr)
        merged = pd.DataFrame(columns=["date"])
    merged = merged.sort_values("date")
    merged.to_csv(cfg.interim_dir / "local_anomaly_index.csv", index=False)

    both_hz = sorted(set(candidates.get("H", [])) & set(candidates.get("Z", [])))

    out = {
        "group": args.group,
        "methods_run": methods_run,
        "fit_H": fits.get("H"),
        "fit_Z": fits.get("Z"),
        "fit_F": fits.get("F"),
        "candidate_dates_H": candidates.get("H", []),
        "candidate_dates_Z": candidates.get("Z", []),
        "candidate_dates_F": candidates.get("F", []),
        "candidate_dates_both_H_and_Z": both_hz,
        "n_storm_or_recovery_dates_excluded": len(storm_dates),
    }
    (cfg.interim_dir / "candidate_windows.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
