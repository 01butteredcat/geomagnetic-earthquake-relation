"""Surrogate-data significance test for the ULF near/far polarization
differential series (`pc3_diff_zh` / `pc4_diff_zh`), per the professor's
first suggestion: since the raw z-score-looks-normal intuition doesn't
rigorously apply to an autocorrelated, non-normal geomagnetic series, shuffle
the series itself (block bootstrap and FFT phase randomization) to get many
surrogate series that look statistically the same but carry no relationship
to any earthquake, then ask how often something as extreme as the observed
dip (or as extreme as the specific z <= -4.1 threshold quoted in
report_template.html for G10's 2024-03-30) shows up by chance.

Runs per (group, band) for every group that has ulf_near_far_index.csv
(the 8 vector-sufficient groups; see stat_utils.py / the plan file for why
G1/G2/G3 are out of scope). Reuses the exact whole-series median/MAD
z-score formula cross_group_analysis.py::ulf_candidate_dates already
established, via stat_utils.mad_zscore, so results here are directly
comparable to (and a rigor upgrade of) that existing analysis.

Usage:
  surrogate_test.py --group G10          # single group, both bands
  surrogate_test.py --all                # every group with ULF data
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import shapiro

sys.path.insert(0, str(Path(__file__).parent))
from common import PROJECT_DIR, load_group_config  # noqa: E402
from stat_utils import (  # noqa: E402
    FIXED_RULE_THRESHOLD,
    block_bootstrap_surrogate,
    estimate_block_length,
    histogram_summary,
    mad_zscore,
    phase_randomize_surrogate,
)

ULF_GROUPS = ("G4", "G5", "G6", "G7", "G8", "G9", "G10", "G11", "G12", "G13", "G19", "G20", "G23", "G24")
BANDS = ("pc3", "pc4")
N_SURROGATES = 2000
SEED = 20260805  # fixed so re-running this script reproduces the same p-values


def run_one(group_id: str, band: str, rng: np.random.Generator) -> dict | None:
    cfg = load_group_config(group_id)
    path = cfg.interim_dir / "ulf_near_far_index.csv"
    if not path.exists():
        return None
    df = pd.read_csv(path, dtype={"date": str}).sort_values("date").reset_index(drop=True)
    col = f"{band}_diff_zh"
    if col not in df.columns:
        return None
    mask = df[col].notna()
    dates = df.loc[mask, "date"].to_numpy()
    vals = df.loc[mask, col].to_numpy(dtype=float)
    n = len(vals)
    if n < 20:
        return {"group": group_id, "band": band, "error": f"only {n} clean days, too few to test"}

    # 1. Normality check -- motivates using median/MAD + surrogates instead
    # of a parametric normal-theory test in the first place.
    shapiro_w, shapiro_p = shapiro(vals)

    # 2. Observed statistic: whole-series median/MAD z-score, most negative day.
    z_obs, med, mad = mad_zscore(vals)
    min_idx = int(np.nanargmin(z_obs))
    t_obs = float(z_obs[min_idx])
    obs_date = str(dates[min_idx])
    n_days_le_threshold = int(np.sum(z_obs <= FIXED_RULE_THRESHOLD))

    block_len = estimate_block_length(vals)

    def null_min_z(surrogate_fn) -> np.ndarray:
        mins = np.empty(N_SURROGATES)
        for i in range(N_SURROGATES):
            surr = surrogate_fn()
            z_surr, _, mad_surr = mad_zscore(surr)
            mins[i] = np.nanmin(z_surr) if mad_surr > 1e-9 else np.nan
        return mins

    block_mins = null_min_z(lambda: block_bootstrap_surrogate(vals, block_len, rng))
    phase_mins = null_min_z(lambda: phase_randomize_surrogate(vals, rng))

    def summarize(mins: np.ndarray, method: str, extra: dict) -> dict:
        valid = mins[~np.isnan(mins)]
        p_vs_observed = float(np.mean(valid <= t_obs)) if len(valid) else float("nan")
        p_vs_threshold = float(np.mean(valid <= FIXED_RULE_THRESHOLD)) if len(valid) else float("nan")
        return {
            "method": method,
            "n_surrogates": int(len(valid)),
            "p_value_vs_observed_extreme": round(p_vs_observed, 5),
            "p_value_vs_fixed_threshold": round(p_vs_threshold, 5),
            "null_min_z_distribution": histogram_summary(valid),
            **extra,
        }

    return {
        "group": group_id,
        "band": band,
        "n_days": n,
        "shapiro": {"W": round(float(shapiro_w), 5), "p": round(float(shapiro_p), 6),
                    "rejects_normality_p05": bool(shapiro_p < 0.05)},
        "observed": {
            "median": round(med, 5), "mad_scaled": round(mad, 5),
            "min_z": round(t_obs, 4), "min_z_date": obs_date,
            "n_days_z_le_fixed_threshold": n_days_le_threshold,
        },
        "fixed_rule_threshold": FIXED_RULE_THRESHOLD,
        "block_bootstrap": summarize(block_mins, "block_bootstrap", {"block_length_days": block_len}),
        "phase_randomization": summarize(phase_mins, "phase_randomization", {}),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--group")
    ap.add_argument("--all", action="store_true")
    args = ap.parse_args()
    if not args.group and not args.all:
        ap.error("pass --group <id> or --all")

    groups = list(ULF_GROUPS) if args.all else [args.group]
    out_dir = PROJECT_DIR / "data" / "interim" / "surrogate_test"
    out_dir.mkdir(parents=True, exist_ok=True)

    for group_id in groups:
        rng = np.random.default_rng(SEED)  # reset per group so band order doesn't matter
        for band in BANDS:
            result = run_one(group_id, band, rng)
            if result is None:
                print(f"[{group_id}/{band}] no data, skipped", file=sys.stderr)
                continue
            out_path = out_dir / f"{group_id}_{band}.json"
            out_path.write_text(json.dumps(result, indent=2))
            if "error" in result:
                print(f"[{group_id}/{band}] {result['error']}", file=sys.stderr)
            else:
                obs = result["observed"]
                bb = result["block_bootstrap"]
                pr = result["phase_randomization"]
                print(f"[{group_id}/{band}] min_z={obs['min_z']} on {obs['min_z_date']} "
                      f"| block-bootstrap p(<=obs)={bb['p_value_vs_observed_extreme']} "
                      f"p(<=-4.1)={bb['p_value_vs_fixed_threshold']} "
                      f"| phase-rand p(<=obs)={pr['p_value_vs_observed_extreme']} "
                      f"p(<=-4.1)={pr['p_value_vs_fixed_threshold']}", file=sys.stderr)


if __name__ == "__main__":
    main()
