"""Superposed epoch analysis (SEA): the professor's "把一次地震變成很多次"
suggestion. Instead of looking at G10's single 2024-04-03 mainshock, take
every independent M>=5.5 earthquake near the network across the 8 groups
that have ULF near/far data (13 hand-picked M>=6.0 anchor/sub-events from
events.py, plus whatever `fetch_earthquake_catalog.py` additionally found),
align each group's pc3/pc4 near-far polarization z-score series on "days
relative to that earthquake's origin time", and stack (average) across all
events. If the kind of dip seen 4 days before G10's mainshock (2024-03-30)
is a real, reproducible precursor signature rather than a one-off, it should
survive averaging and stand out against a null band built the same way from
random, earthquake-unrelated reference dates.

Per-group z-scoring (via stat_utils.mad_zscore, same formula used
throughout this validation suite) happens BEFORE stacking, since raw
diff_zh magnitudes aren't on a comparable scale across groups/stations with
different noise floors -- stacking raw nT-order differences across groups
would be meaningless.

Usage:
  superposed_epoch_analysis.py --catalog data/external/extended_catalog_m5.5.csv --label m5.5
  superposed_epoch_analysis.py --catalog data/external/extended_catalog_m5.0.csv --label m5.0
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from catalog_utils import load_extended_events  # noqa: E402
from common import PROJECT_DIR, load_group_config  # noqa: E402
from stat_utils import mad_zscore  # noqa: E402

ULF_GROUPS = ("G4", "G5", "G6_G7_G8", "G9", "G10", "G11", "G12", "G13", "G19", "G20", "G23")
BANDS = ("pc3", "pc4")
WINDOW_BEFORE_DAYS = 30
WINDOW_AFTER_DAYS = 10
N_BOOTSTRAP = 2000       # CI on the real stack (resample events with replacement)
N_NULL = 1000            # null realizations (resample the EPOCH DATE at random)
NULL_EXCLUSION_BUFFER_DAYS = 30  # keep fake epochs this far from any real event
SEED = 20260805


def load_group_series(group_id: str) -> dict | None:
    cfg = load_group_config(group_id)
    path = cfg.interim_dir / "ulf_near_far_index.csv"
    if not path.exists():
        return None
    df = pd.read_csv(path, dtype={"date": str}).sort_values("date").reset_index(drop=True)
    out = {"dates": df["date"].tolist()}
    for band in BANDS:
        col = f"{band}_diff_zh"
        if col not in df.columns:
            continue
        vals = df[col].to_numpy(dtype=float)
        z, med, mad = mad_zscore(vals[~np.isnan(vals)]) if np.any(~np.isnan(vals)) else (None, None, None)
        # re-expand z back to full (with-NaN) length, aligned to df rows
        z_full = np.full(len(vals), np.nan)
        if med is not None:
            mask = ~np.isnan(vals)
            z_full[mask] = z
        out[band] = dict(zip(df["date"], z_full))
    return out


def extract_window(series: dict, band: str, event_date: pd.Timestamp) -> np.ndarray:
    lags = np.arange(-WINDOW_BEFORE_DAYS, WINDOW_AFTER_DAYS + 1)
    out = np.full(len(lags), np.nan)
    band_map = series.get(band)
    if band_map is None:
        return out
    for i, d in enumerate(lags):
        date_str = (event_date + pd.Timedelta(days=int(d))).strftime("%Y%m%d")
        v = band_map.get(date_str)
        if v is not None:
            out[i] = v
    return out


def valid_date_range(series: dict) -> tuple[pd.Timestamp, pd.Timestamp]:
    ds = [pd.Timestamp(d) for d in series["dates"]]
    return min(ds), max(ds)


def run_band(band: str, events: list[dict], group_series: dict, rng: np.random.Generator) -> dict:
    lags = np.arange(-WINDOW_BEFORE_DAYS, WINDOW_AFTER_DAYS + 1)
    matrix = []
    used_events = []
    for ev in events:
        series = group_series.get(ev["group"])
        if series is None or band not in series:
            continue
        w = extract_window(series, band, ev["date"])
        if np.all(np.isnan(w)):
            continue
        matrix.append(w)
        used_events.append({"group": ev["group"], "date": ev["date"].strftime("%Y-%m-%d"),
                             "mag": ev["mag"], "source": ev["source"]})
    if not matrix:
        return {"band": band, "error": "no events with usable data"}
    M = np.array(matrix)  # n_events x n_lags

    real_stack_mean = np.nanmean(M, axis=0)
    real_stack_median = np.nanmedian(M, axis=0)
    n_contributing = np.sum(~np.isnan(M), axis=0)

    # bootstrap CI on the real stack (resample events with replacement)
    n_events = M.shape[0]
    boot = np.empty((N_BOOTSTRAP, M.shape[1]))
    for b in range(N_BOOTSTRAP):
        idx = rng.integers(0, n_events, size=n_events)
        boot[b] = np.nanmean(M[idx], axis=0)
    ci_lo = np.nanpercentile(boot, 5, axis=0)
    ci_hi = np.nanpercentile(boot, 95, axis=0)

    # null band: repeat the whole stacking procedure with random, earthquake-
    # unrelated epoch dates (same group membership/event count preserved)
    real_dates_by_group: dict[str, list[pd.Timestamp]] = {}
    for ev in events:
        real_dates_by_group.setdefault(ev["group"], []).append(ev["date"])

    valid_ranges = {g: valid_date_range(s) for g, s in group_series.items()}

    null_stacks = np.empty((N_NULL, M.shape[1]))
    for r in range(N_NULL):
        null_matrix = []
        for ev in events:
            series = group_series.get(ev["group"])
            if series is None or band not in series:
                continue
            lo, hi = valid_ranges[ev["group"]]
            lo_bound = lo + pd.Timedelta(days=WINDOW_BEFORE_DAYS)
            hi_bound = hi - pd.Timedelta(days=WINDOW_AFTER_DAYS)
            if lo_bound >= hi_bound:
                continue
            for _ in range(20):  # retry a few times to avoid the exclusion buffer
                offset_days = int(rng.integers(0, max(1, (hi_bound - lo_bound).days + 1)))
                fake_date = lo_bound + pd.Timedelta(days=offset_days)
                too_close = any(abs((fake_date - rd).days) < NULL_EXCLUSION_BUFFER_DAYS
                                 for rd in real_dates_by_group.get(ev["group"], []))
                if not too_close:
                    break
            w = extract_window(series, band, fake_date)
            if not np.all(np.isnan(w)):
                null_matrix.append(w)
        null_stacks[r] = np.nanmean(np.array(null_matrix), axis=0) if null_matrix else np.nan

    null_lo = np.nanpercentile(null_stacks, 5, axis=0)
    null_hi = np.nanpercentile(null_stacks, 95, axis=0)
    null_median = np.nanpercentile(null_stacks, 50, axis=0)

    return {
        "band": band,
        "n_events_used": len(used_events),
        "events_used": used_events,
        "lags_days": lags.tolist(),
        "stack_mean": [None if np.isnan(v) else round(float(v), 4) for v in real_stack_mean],
        "stack_median": [None if np.isnan(v) else round(float(v), 4) for v in real_stack_median],
        "n_contributing_per_lag": n_contributing.tolist(),
        "bootstrap_ci90_lo": [round(float(v), 4) for v in ci_lo],
        "bootstrap_ci90_hi": [round(float(v), 4) for v in ci_hi],
        "null_band_p5": [round(float(v), 4) for v in null_lo],
        "null_band_p50": [round(float(v), 4) for v in null_median],
        "null_band_p95": [round(float(v), 4) for v in null_hi],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--catalog", type=Path, required=True)
    ap.add_argument("--label", required=True, help="e.g. m5.5 or m5.0, used in output filenames")
    args = ap.parse_args()

    group_series = {}
    for g in ULF_GROUPS:
        s = load_group_series(g)
        if s is not None:
            group_series[g] = s

    events = load_extended_events(args.catalog, ULF_GROUPS)
    print(f"{len(events)} candidate events across {len(group_series)} groups", file=sys.stderr)

    out_dir = PROJECT_DIR / "data" / "interim" / "superposed_epoch"
    out_dir.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(SEED)
    for band in BANDS:
        result = run_band(band, events, group_series, rng)
        out_path = out_dir / f"{band}_stack_{args.label}.json"
        out_path.write_text(json.dumps(result, indent=2))
        if "error" in result:
            print(f"[{band}] {result['error']}", file=sys.stderr)
        else:
            print(f"[{band}] wrote {out_path} ({result['n_events_used']} events used)", file=sys.stderr)


if __name__ == "__main__":
    main()
