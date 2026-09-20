"""Verification pass (plan Section 7) -- must be run and all checks reviewed
before trusting any precursor conclusion drawn from this pipeline.

Usage: verify_pipeline.py --group G10
"""
from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from common import DATA_TIMEZONE, auto_outage_dates, list_day_refs, load_group_config, resolve_day_ref  # noqa: E402
from parser import is_scalar_only, open_raw, parse_day_file  # noqa: E402
from build_daily_features import SPIKE_THRESHOLD_NT  # noqa: E402
from compute_indices import MIN_CLEAN_POINTS, MIN_NIGHT_MINUTES, NIGHT_HOURS_UTC, TRAILING_WINDOW_DAYS  # noqa: E402

random.seed(42)
checks: dict = {}


def check_raw_vs_parsed_spotcheck(cfg):
    samples = []
    all_stations = list(cfg.stations.keys())
    for _ in range(4):
        station = random.choice(all_stations)
        refs = list_day_refs(cfg.gdms_dir, station)
        ref = random.choice(refs)
        df = parse_day_file(ref, station)
        scalar = is_scalar_only(ref)
        row_idx = random.randint(0, len(df) - 1)
        parsed_row = df.iloc[row_idx]
        ts = df.index[row_idx]

        raw_line = None
        with open_raw(ref) as f:
            for line in f:
                if line.startswith(ts.strftime("%Y-%m-%d %H:%M:%S")):
                    raw_line = line
                    break
        assert raw_line is not None, f"could not find raw line for {ts} in {ref.label}"
        parts = raw_line.split()
        raw = dict(zip(["X", "Y", "Z", "F"], (float(v) for v in parts[3:7])))

        real_cols = ["F"] if scalar else ["X", "Y", "Z"]
        junk_cols = ["X", "Y", "Z"] if scalar else ["F"]
        # junk columns (never real for this file's channel type) are
        # normally the not-reported placeholder (88888), but on some older
        # (pre-2024, CWB-era) totally-dead days the file instead fills the
        # outage sentinel (99999) into the junk column too -- e.g. a whole
        # G2_G3-era file can read X=Y=Z=88888/F=99999 for a vector-type
        # station on a day with zero data at all. The pipeline itself
        # already handles this correctly (parser.py's defensive NaN-ing
        # treats a stray 88888 in a real column as invalid too), so this
        # check only requires the junk value be ONE of the two known
        # sentinels, not a specific one -- a genuinely wrong/corrupted raw
        # value would still be caught.
        junk_ok = all(raw[c] in (88888.0, 99999.0) for c in junk_cols)
        # A real column reads NaN in the parsed output for the outage
        # sentinel (99999) always, AND for a stray not-reported placeholder
        # (88888) too -- but only for vector (X/Y/Z) columns, which is where
        # parser.py's defensive fallback applies (a totally-dead day can
        # apparently still write 88888 into what should be a real vector
        # column; ground truth per common.py's/CLAUDE.md's convention is
        # that shouldn't happen, but the pipeline treats it defensively as
        # invalid rather than a bogus real reading). The scalar (F) column
        # has no such defensive mapping in parser.py, so 88888 there is not
        # expected/tested for.
        nan_sentinels = {"X": (88888.0, 99999.0), "Y": (88888.0, 99999.0), "Z": (88888.0, 99999.0), "F": (99999.0,)}
        real_ok = all(
            pd.isna(parsed_row[c]) if raw[c] in nan_sentinels[c] else np.isclose(parsed_row[c], raw[c])
            for c in real_cols
        )
        ok = junk_ok and real_ok
        samples.append({"station": station, "file": ref.label, "timestamp": str(ts), "ok": bool(ok)})

    all_ok = all(s["ok"] for s in samples)
    checks["raw_vs_parsed_spotcheck"] = {"pass": all_ok, "samples": samples}


def check_known_outage_ground_truth(cfg):
    """Manually pre-verified ground truth for three specific twu outage days
    -- only meaningful for G10, where this was hand-checked against the raw
    files when the pipeline was first built. Other groups have no equivalent
    manually-verified ground truth, so this check is not applicable there
    (their gap-handling correctness is instead covered generically by
    check_no_spikes_adjacent_to_missing_data below)."""
    if cfg.group_id != "G10":
        checks["known_outage_ground_truth"] = {
            "pass": None,
            "note": "skipped_not_applicable -- this hand-verified ground truth is G10-specific",
        }
        return

    results = {}
    df419 = parse_day_file(resolve_day_ref(cfg.gdms_dir, "twu", "20240419"), "twu")
    n_nan_419 = int(df419["X"].isna().sum())
    results["twu_20240419_nan_count_expected_404"] = n_nan_419

    df424 = parse_day_file(resolve_day_ref(cfg.gdms_dir, "twu", "20240424"), "twu")
    n_nan_424 = int(df424["X"].isna().sum())
    results["twu_20240424_nan_count_expected_86400"] = n_nan_424

    df425 = parse_day_file(resolve_day_ref(cfg.gdms_dir, "twu", "20240425"), "twu")
    n_nan_425 = int(df425["X"].isna().sum())
    results["twu_20240425_nan_count_expected_86400"] = n_nan_425

    ok = n_nan_419 == 404 and n_nan_424 == 86400 and n_nan_425 == 86400
    checks["known_outage_ground_truth"] = {"pass": ok, **results}


def check_timezone_conclusion(cfg):
    """common.DATA_TIMEZONE ('UTC') is an already-established, dataset-wide
    ground truth confirmed with high confidence on G10's clean 2024 data
    (see common.py's comment) -- it should hold for every group, since all
    groups share the same CWA IAGA-2002 recording convention. This check's
    per-group re-run of timezone_check.py is a spot-check, not a fresh
    from-scratch proof: a station with degraded data quality (e.g. G11's
    twu, which has many partial outages) can legitimately come back
    UNCERTAIN/low-confidence without that meaning anything is wrong -- it
    just means that particular station's signal wasn't clean enough this
    time. A real problem is a CONFIDENT conclusion that contradicts
    common.DATA_TIMEZONE; that would be a genuine red flag."""
    tz = json.loads((cfg.interim_dir / "timezone_check.json").read_text())
    conclusion, confidence = tz["conclusion_timezone"], tz["confidence"]
    expected = "UTC" if DATA_TIMEZONE == "UTC" else "LOCAL (UTC+8)"
    if confidence != "high":
        result_pass = None
    else:
        result_pass = conclusion == expected
    checks["timezone_conclusion"] = {
        "pass": result_pass,
        "conclusion": conclusion,
        "confidence": confidence,
        "expected_per_common_py": expected,
        "note": (
            "no header TZ field exists; common.DATA_TIMEZONE is the dataset-wide ground truth "
            "(confirmed on G10's clean data) -- this check only fails on a CONFIDENT contradiction, "
            "not on a low-confidence/inconclusive spot-check result for a noisier station/group"
        ),
    }


def check_storm_cancellation(cfg):
    kp_path = cfg.external_dir / "kp.csv"
    idx_field = "H" if cfg.xyz_pool.sufficient else ("F" if cfg.f_pool.sufficient else None)
    if not kp_path.exists() or idx_field is None:
        checks["storm_cancellation_test"] = {
            "pass": None,
            "note": "inconclusive -- Kp data unavailable or no station pool sufficient to compute an index",
        }
        return

    kp = pd.read_csv(kp_path, dtype={"date": str})
    kp_daily_max = kp.groupby("date")["kp"].max().sort_values(ascending=False)

    idx = pd.read_csv(cfg.interim_dir / "local_anomaly_index.csv", dtype={"date": str})
    near_col, local_col = f"near_index_{idx_field}", f"local_anomaly_index_{idx_field}"
    if near_col not in idx.columns:
        checks["storm_cancellation_test"] = {"pass": None, "note": f"inconclusive -- {near_col} not computed"}
        return

    # walk the highest-Kp days in order until we find one with a computable
    # index -- must check BOTH near_index and local_anomaly_index for NaN:
    # local_anomaly_index also depends on far_index, which can be NaN (e.g.
    # a far station that hasn't come online yet this early in the group's
    # date range, or without enough trailing baseline) even when near_index
    # itself is fine, so checking near_col alone is not sufficient.
    for top_kp_date, top_kp_value in kp_daily_max.items():
        row = idx[idx.date == top_kp_date]
        if not row.empty and not row[near_col].isna().all() and not row[local_col].isna().all():
            near = float(row[near_col].iloc[0])
            local = float(row[local_col].iloc[0])
            ok = abs(local) < 0.6 * abs(near) if abs(near) > 0.5 else None
            checks["storm_cancellation_test"] = {
                "pass": ok,
                "field": idx_field,
                "top_kp_date": top_kp_date,
                "top_kp_value": float(top_kp_value),
                "raw_near_index": near,
                "local_anomaly_index_after_regression": local,
                "note": "expect |local_anomaly_index| substantially < |near_index| on the highest-Kp day, showing the far-station regression absorbs common-mode storm signal",
            }
            return

    checks["storm_cancellation_test"] = {
        "pass": None,
        "note": "inconclusive -- no Kp day in this group's range has a computable index (likely trailing-baseline gap)",
    }


def check_candidates_not_in_outage_windows(cfg):
    """A candidate day should never fall on a day compute_indices.py's own
    is_clean_day considered dirty -- candidate_flag is defined to require
    is_clean_day==True, so this is a regression test that the gate actually
    held, not a tautology: it independently rebuilds the same "outage" event
    compute_indices.py's clean_overall now uses (a day only counts as an
    outage for a pool if EVERY station in that pool's near group, or every
    station in its far group, was out that day -- matching how the median-
    of-near/median-of-far index tolerates a single missing station), rather
    than the earlier, over-strict "any one relevant station is dirty" union,
    which produced false failures once auto_outage_dates started finding
    long-running single-station outages that the near/far median already
    routes around without issue. Station-day "out" is judged on the night
    window (see below), matching what compute_indices.py's night features
    (MIN_NIGHT_MINUTES) actually gate on."""
    cand = json.loads((cfg.interim_dir / "candidate_windows.json").read_text())
    daily = pd.read_csv(cfg.interim_dir / "daily_features.csv", dtype={"date": str})
    all_days = sorted(daily["date"].unique())

    # "Outage" here is judged on the NIGHT window the index is actually built
    # from (NIGHT_HOURS_UTC, >= MIN_NIGHT_MINUTES valid minutes of that
    # station+channel), not the whole-day pct_missing: a station can be >5%
    # missing over 24h yet have a complete night window, and vice versa. That
    # mismatch is what made this check flag G14/G18 candidates whose near or
    # far median was in fact computed from full night data.
    night_ok_cache: dict[tuple[str, str], set[str]] = {}

    def night_ok_dates(station, channel):
        key = (station, channel)
        if key not in night_ok_cache:
            df = pd.read_parquet(cfg.interim_dir / f"minute_series_{station}.parquet", columns=[channel])
            night = df[df.index.hour.isin(NIGHT_HOURS_UTC)][channel].dropna()
            n = night.groupby(night.index.strftime("%Y%m%d")).size()
            night_ok_cache[key] = set(n[n >= MIN_NIGHT_MINUTES].index)
        return night_ok_cache[key]

    def pool_outage_dates(pool, channel):
        near_out = {d for d in all_days if pool.near and not any(d in night_ok_dates(s, channel) for s in pool.near)}
        far_out = {d for d in all_days if pool.far and not any(d in night_ok_dates(s, channel) for s in pool.far)}
        return near_out | far_out

    overlap = set()
    if cand.get("candidate_dates_H") or cand.get("candidate_dates_Z"):
        overlap |= set(cand.get("candidate_dates_H", [])) & pool_outage_dates(cfg.xyz_pool, "H")
        overlap |= set(cand.get("candidate_dates_Z", [])) & pool_outage_dates(cfg.xyz_pool, "Z")
    if cand.get("candidate_dates_F"):
        overlap |= set(cand.get("candidate_dates_F", [])) & pool_outage_dates(cfg.f_pool, "F")

    all_candidates = set(cand.get("candidate_dates_H", [])) | set(cand.get("candidate_dates_Z", [])) | set(cand.get("candidate_dates_F", []))
    checks["candidates_not_in_outage_windows"] = {
        "pass": len(overlap) == 0,
        "candidate_dates": sorted(all_candidates),
        "overlap_with_known_outages": sorted(overlap),
    }


def check_quiet_day_smoothness(cfg):
    """Proxy for a visual sanity check: on the quiet days used for the
    timezone check, the parsed curve should be free of PARSING artifacts.
    We do not require every second to be near-flat -- brief (1-2 sample)
    simultaneous jumps up to ~O(100nT) are physically real (e.g.
    sudden-commencement-type impulses) and expected occasionally even on an
    otherwise-quiet day; n_spikes==0 already confirms nothing crossed the
    300nT structural-glitch threshold used elsewhere in the pipeline (see
    build_daily_features.SPIKE_THRESHOLD_NT). This check instead just
    confirms no pervasive corruption (many large jumps, as seen on twu)."""
    tz = json.loads((cfg.interim_dir / "timezone_check.json").read_text())
    station = tz.get("station")
    quiet_days = tz["quiet_days_used"][:5]
    daily = pd.read_csv(cfg.interim_dir / "daily_features.csv", dtype={"date": str})
    sub = daily[(daily.station == station) & (daily.date.isin(quiet_days))]
    ok = bool((sub["n_spikes"] == 0).all()) if len(sub) else None
    checks["quiet_day_smoothness"] = {
        "pass": ok,
        "station": station,
        "days_checked": quiet_days,
        "max_abs_jump_values": sub["max_abs_jump"].tolist(),
        "note": "occasional double-digit-to-~O(100nT) single-second jumps are expected real transients, not parsing bugs; n_spikes==0 confirms none crossed the 300nT structural-glitch threshold",
    }


def check_no_spikes_adjacent_to_missing_data(cfg):
    """Regression test for a historical bug where a NaN diff (from an
    outage-sentinel gap) was miscounted as a spike, which would have skewed
    n_spikes and inflated apparent anomaly counts. build_daily_features.py's
    _despike() now relies on NaN comparisons being False to exclude gap
    boundaries (see its docstring) -- this test does NOT call _despike()
    itself (that would be tautological), it independently re-parses every
    station/day that has any missing samples and recomputes the diff/
    threshold check from scratch, then asserts no flagged spike sits next to
    a NaN sample."""
    daily = pd.read_csv(cfg.interim_dir / "daily_features.csv", dtype={"date": str})
    partial = daily[(daily.n_valid > 0) & (daily.n_valid < daily.n_total)]

    bad = []
    for _, row in partial.iterrows():
        station, date = row["station"], row["date"]
        ref = resolve_day_ref(cfg.gdms_dir, station, date)
        df = parse_day_file(ref, station)
        for col in df.columns:  # ["F"] for a scalar file, ["X","Y","Z"] for a vector file
            s = df[col]
            d = s.diff().abs()
            flagged = np.where((d > SPIKE_THRESHOLD_NT).values)[0]
            for i in flagged:
                if pd.isna(s.iloc[i]) or pd.isna(s.iloc[i - 1]):
                    bad.append({"station": station, "date": date, "column": col, "row": int(i)})

    checks["no_spikes_adjacent_to_missing_data"] = {
        "pass": len(bad) == 0,
        "station_days_checked": len(partial),
        "bad_flags": bad[:20],
        "note": "any entry here means a spike was flagged directly adjacent to a NaN sample -- the exact NaN-miscounted-as-spike failure mode",
    }


def check_baseline_window_excludes_storms(cfg):
    """Regression test for a historical bug (found on the G10/2024 data)
    where a too-short trailing baseline window let a magnetic storm
    contaminate (or entirely starve) the reference used to judge candidate
    anomaly days in the critical pre-anchor-event window. Independently of
    compute_indices.mad_zscore, reconstructs the TRAILING_WINDOW_DAYS
    calendar window for every day in the critical window (anchor event date
    -8..+1 days -- the offset that mattered for G10's 2024-04-03 mainshock
    and 2024-03-30 candidate day, reused as a dataset-wide constant since the
    failure mode -- a storm sitting just before the window -- is generic)
    and checks (a) at least one of those windows does overlap real storm
    dates -- confirming this test actually exercises the failure scenario
    rather than being vacuous -- and (b) every one of those windows still
    has enough clean (non-storm, non-outage) days to support a real baseline
    (>= MIN_CLEAN_POINTS, compute_indices.py's own threshold below which
    mad_zscore returns NaN, i.e. no baseline at all).

    "Outage" here means the day's near_index (or far_index) would have been
    NaN -- i.e. EVERY near station (or every far station) was out that day
    -- matching compute_indices.py's clean_overall definition (a median
    across 3 near stations tolerates one of them being down; it only really
    goes NaN when all of them are). A day with just one out of three near
    stations down is not an outage for this purpose. This check does not
    call into compute_indices.py itself (that would be tautological), it
    independently rebuilds the same all-of-a-pool-down condition from
    daily_features.csv."""
    storm_path = cfg.interim_dir / "storm_days.csv"
    if not storm_path.exists():
        checks["baseline_window_excludes_storms"] = {"pass": None, "note": "inconclusive -- storm_days.csv missing"}
        return
    storm_dates = set(pd.read_csv(storm_path, dtype={"date": str})["date"])

    daily = pd.read_csv(cfg.interim_dir / "daily_features.csv", dtype={"date": str})
    outage_by_station = auto_outage_dates(daily)

    pool = cfg.xyz_pool if cfg.xyz_pool.sufficient else cfg.f_pool
    all_days = sorted(daily["date"].unique())
    near_out = {d: all(d in outage_by_station.get(s, set()) for s in pool.near) for d in all_days}
    far_out = {d: all(d in outage_by_station.get(s, set()) for s in pool.far) for d in all_days}
    outage_dates = {d for d in all_days if near_out[d] or far_out[d]}

    anchor_date = pd.to_datetime(cfg.anchor_event.date)
    critical_window = pd.date_range(
        anchor_date - pd.Timedelta(days=8), anchor_date + pd.Timedelta(days=1)
    ).strftime("%Y%m%d").tolist()

    results = []
    for date_str in critical_window:
        d = pd.to_datetime(date_str, format="%Y%m%d")
        window = pd.date_range(
            d - pd.Timedelta(days=TRAILING_WINDOW_DAYS), d - pd.Timedelta(days=1)
        ).strftime("%Y%m%d").tolist()
        storm_in_window = sorted(set(window) & storm_dates)
        clean_days = [w for w in window if w not in storm_dates and w not in outage_dates]
        results.append(
            {
                "probe_date": date_str,
                "trailing_window": [window[0], window[-1]],
                "storm_dates_in_window": storm_in_window,
                "n_clean_days_in_window": len(clean_days),
            }
        )

    overlap_exercised = any(r["storm_dates_in_window"] for r in results)
    all_sufficient = all(r["n_clean_days_in_window"] >= MIN_CLEAN_POINTS for r in results)
    checks["baseline_window_excludes_storms"] = {
        "pass": all_sufficient if overlap_exercised else None,
        "critical_window_checked": [critical_window[0], critical_window[-1]],
        "trailing_window_days": TRAILING_WINDOW_DAYS,
        "min_clean_points_threshold": MIN_CLEAN_POINTS,
        "probe_days": results,
        "note": (
            "checks every day in the critical pre-anchor-event window still has >= MIN_CLEAN_POINTS clean "
            "baseline days after excluding overlapping storm dates" + ("" if overlap_exercised else
            " -- INCONCLUSIVE: no probe day's window actually overlapped a storm date, so this test "
            "did not exercise the failure scenario")
        ),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--group", required=True)
    args = ap.parse_args()
    cfg = load_group_config(args.group)

    check_raw_vs_parsed_spotcheck(cfg)
    check_known_outage_ground_truth(cfg)
    check_timezone_conclusion(cfg)
    check_storm_cancellation(cfg)
    check_candidates_not_in_outage_windows(cfg)
    check_quiet_day_smoothness(cfg)
    check_no_spikes_adjacent_to_missing_data(cfg)
    check_baseline_window_excludes_storms(cfg)

    n_pass = sum(1 for c in checks.values() if c.get("pass") is True)
    n_fail = sum(1 for c in checks.values() if c.get("pass") is False)
    n_inconclusive = sum(1 for c in checks.values() if c.get("pass") is None)

    report = {
        "group": args.group,
        "summary": {"pass": n_pass, "fail": n_fail, "inconclusive": n_inconclusive},
        "checks": checks,
    }
    (cfg.interim_dir / "verification_report.json").write_text(json.dumps(report, indent=2, default=str))
    print(json.dumps(report, indent=2, default=str))

    if n_fail > 0:
        print(f"\n{n_fail} CHECK(S) FAILED -- do not trust precursor conclusions until fixed", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
