"""Dose-response test: does the coseismic geomagnetic anomaly grow with the
ground shaking at the magnetometer?

If what the 1Hz magnetometers record near the origin second is the housing
being shaken, anomaly size should rise with peak ground acceleration (PGA)
at a seismic station near the magnetometer. This needs no noise/signal labels,
unlike coseismic_joint_analysis.py -- and those labels turned out to be mostly
meaningless: only 7 of 93 events have a significant (D7 p < 0.05) anomaly to
time against the shaking in the first place (2026-09-25).

Unit: one row per event that seismometer_comparison.py compared successfully
(comparison_summary.csv, status == "ok", PGA available).

  x  = log10 PGA (gal), max horizontal, at the seismic station nearest the
       magnetometer (<= 12 km, median 7 km) -- seismometer_comparison.py
  y1 = |step30 peak z| within +-180 s (geomag_obs_peak_z; already scaled by
       each station's own off-event noise, so comparable across stations)
  y2 = 1Hz first-difference RMS during shaking / 10 min before
       (geomag_noise_ratio; the shaking-noise signature of the G10 pilot)
  y3 = -log10 of the event's calibrated step30 p-value (geomag_step30_p)

Test: Spearman rho(x, y), one-sided (H1: rho > 0). The p-value comes from
permuting y *within each group* (N_PERM times): G11 alone is ~a quarter of the
events, all on the same magnetometer, and a pooled permutation would let one
group's internal spread drive the result.

Subsets: all; M>=6; seismic station <= 2 km from the magnetometer (the only
cases where the measured PGA really is the shaking at the magnetometer); and
all minus WTP (the station with the largest PGAs, near-field Dapu aftershocks).

Usage:
  coseismic_dose_response.py --self-test
  coseismic_dose_response.py --all
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).parent))
import common  # noqa: E402
from coseismic_step_analysis import SEED  # noqa: E402
from events import GROUPS  # noqa: E402

N_PERM = 2000
COLOCATED_KM = 2.0
PGA_BINS_GAL = (0, 1, 10, 100, np.inf)
MIN_N = 8  # fewer events than this in a subset and no test is attempted

COMPARISON_CSV = common.PROJECT_DIR / "data" / "interim" / "seismometer_comparison" / "comparison_summary.csv"
OUT_DIR = common.PROJECT_DIR / "data" / "interim" / "coseismic_dose_response"

RESPONSES = {
    "peak_abs_z": "地磁 step30 峰值 |z|",
    "noise_ratio": "地磁 1 Hz 雜訊放大倍數",
    "neg_log10_p": "逐事件檢定 −log10 p",
}


def load_events(csv_path: Path = COMPARISON_CSV) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    df = df[(df.status == "ok") & df.pga_gal.notna() & (df.pga_gal > 0)].copy()
    mag = {(gid, e.date): e.magnitude for gid, g in GROUPS.items() for e in g.events}
    df["magnitude"] = [mag.get((g, d)) for g, d in zip(df.group, df.date)]
    df["log10_pga"] = np.log10(df.pga_gal)
    df["peak_abs_z"] = df.geomag_obs_peak_z.abs()
    df["noise_ratio"] = df.geomag_noise_ratio
    df["neg_log10_p"] = -np.log10(df.geomag_step30_p)
    return df


def within_group_permutation_test(x: np.ndarray, y: np.ndarray, groups: np.ndarray,
                                   rng: np.random.Generator, n_perm: int = N_PERM) -> dict:
    """Spearman rho and a one-sided (rho > 0) p-value from shuffling y only
    among events of the same group."""
    rho = float(spearmanr(x, y).statistic)
    idx_by_group = [np.flatnonzero(groups == g) for g in np.unique(groups)]
    null = np.empty(n_perm)
    for i in range(n_perm):
        yp = y.copy()
        for ix in idx_by_group:
            if len(ix) > 1:
                yp[ix] = y[rng.permutation(ix)]
        null[i] = spearmanr(x, yp).statistic
    null = null[~np.isnan(null)]
    p = float((1 + np.sum(null >= rho)) / (len(null) + 1)) if len(null) else float("nan")
    return {"n": int(len(x)), "n_groups": int(len(idx_by_group)), "spearman_rho": round(rho, 4),
            "p_one_sided_within_group": round(p, 5)}


def test_subset(df: pd.DataFrame, rng: np.random.Generator) -> dict:
    out = {}
    for resp in RESPONSES:
        d = df[df[resp].notna() & np.isfinite(df[resp])]
        if len(d) < MIN_N:
            out[resp] = {"n": int(len(d)), "note": f"fewer than {MIN_N} events, not tested"}
            continue
        out[resp] = within_group_permutation_test(d.log10_pga.to_numpy(), d[resp].to_numpy(),
                                                   d.group.to_numpy(), rng)
    return out


def pga_bins(df: pd.DataFrame) -> list[dict]:
    rows = []
    cut = pd.cut(df.pga_gal, PGA_BINS_GAL, right=False)
    for b, d in df.groupby(cut, observed=False):
        rows.append({
            "pga_bin_gal": f"{b.left:g}–{b.right:g}" if np.isfinite(b.right) else f"≥{b.left:g}",
            "n": int(len(d)),
            "median_peak_abs_z": None if d.empty else round(float(d.peak_abs_z.median()), 3),
            "median_noise_ratio": None if d.noise_ratio.dropna().empty else round(float(d.noise_ratio.median()), 3),
            "n_significant_anomaly": int((d.geomag_step30_p < 0.05).sum()),
        })
    return rows


def run_all() -> dict:
    df = load_events()
    rng = np.random.default_rng(SEED)
    subsets = {
        "all": df,
        "m6": df[df.magnitude >= 6],
        "colocated_le_2km": df[df.seismic_colocation_km <= COLOCATED_KM],
        "all_minus_WTP": df[df.seismic_station != "WTP"],
        "all_minus_dropout": df[df.geomag_dropout != True],  # noqa: E712 -- keeps NaN (unknown) rows
    }
    summary = {"seed": SEED, "n_perm": N_PERM, "tail": "one-sided, H1: rho > 0",
               "subsets": {k: test_subset(v, rng) for k, v in subsets.items()},
               "pga_bins_all": pga_bins(df)}
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    cols = ["group", "date", "magnitude", "geomag_station", "seismic_station", "seismic_colocation_km",
            "pga_gal", "pga_component", "peak_abs_z", "noise_ratio", "geomag_step30_p",
            "geomag_dropout", "alignment_verdict", "alignment_verdict_gated"]
    df[cols].to_csv(OUT_DIR / "per_event.csv", index=False)
    (OUT_DIR / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False))
    (OUT_DIR / "summary.md").write_text(render_markdown(summary, df), encoding="utf-8")
    for name, res in summary["subsets"].items():
        for resp, r in res.items():
            print(f"[{name}] {resp}: {r}", file=sys.stderr)
    return summary


def render_markdown(summary: dict, df: pd.DataFrame) -> str:
    label = {"all": "全部事件", "m6": "只用 M≥6",
             "colocated_le_2km": f"地震站距磁力儀 ≤ {COLOCATED_KM:g} km", "all_minus_WTP": "排除 WTP 站",
             "all_minus_dropout": "排除地磁儀器中斷的事件"}
    L = ["# 同震劑量反應檢定：地磁異常 vs 磁力儀附近的地動強度（PGA）", "",
         "H₀：地磁異常大小與 PGA 無關。H₁：地磁異常隨 PGA 增加（單尾）。"
         f"Spearman 等級相關；p 值由組內置換 {summary['n_perm']} 次得到（只在同一組內打亂）。", "",
         "| 子集 | 反應變數 | 事件數 | 組數 | Spearman ρ | p（單尾） |", "|---|---|---|---|---|---|"]
    for name, res in summary["subsets"].items():
        for resp, r in res.items():
            if "spearman_rho" not in r:
                L.append(f"| {label[name]} | {RESPONSES[resp]} | {r['n']} | — | — | 事件數不足 |")
                continue
            L.append(f"| {label[name]} | {RESPONSES[resp]} | {r['n']} | {r['n_groups']} "
                     f"| {r['spearman_rho']} | {r['p_one_sided_within_group']} |")
    L += ["", "## 依 PGA 分級（全部事件）", "",
          "| PGA（gal） | 事件數 | 峰值 \\|z\\| 中位數 | 雜訊放大倍數中位數 | 顯著異常（p<0.05）事件數 |",
          "|---|---|---|---|---|"]
    for b in summary["pga_bins_all"]:
        L.append(f"| {b['pga_bin_gal']} | {b['n']} | {b['median_peak_abs_z']} | {b['median_noise_ratio']} "
                 f"| {b['n_significant_anomaly']} |")
    L += ["", "## 限制", "",
          "- PGA 是在離磁力儀最近的地震站量的（中位數 7 km、最遠約 12 km），不是磁力儀本身的搖晃；只有少數事件的地震站在 2 km 內。",
          "- 3 種反應變數 × 4 個子集同時檢定，未做多重比較校正。",
          "- 地磁儀器中斷：震後 300 秒內缺值超過 20% 的事件（例如 2025-01-21 大埔主震時的 twu，震動一開始讀數就停住、接著資料中斷）。這種事件的「無異常」其實是儀器停擺；主要分析保留它們（對 H₁ 不利、較保守），另附排除後的結果。",
          "- 雜訊放大倍數用 RMS，震前窗口若有突跳（例如 xcg 的幾起 M5 級事件，震前有 2–4 nT 突跳）會低估比值；這只會削弱相關，不會製造假的正相關。",
          "- WTP 站在大埔餘震序列的 PGA 很大（最高約 2000 gal），峰值不是單點突波，視為真實近場紀錄；另附排除 WTP 的結果。", ""]
    return "\n".join(L)


def self_test() -> bool:
    rng = np.random.default_rng(SEED)
    n, groups = 96, np.repeat(np.arange(12), 8)
    x = rng.normal(1.3, 0.8, n)
    y_dose = 0.8 * x + rng.normal(0, 0.6, n)
    y_null = rng.normal(0, 1, n)
    a = within_group_permutation_test(x, y_dose, groups, np.random.default_rng(SEED), 500)
    b = within_group_permutation_test(x, y_null, groups, np.random.default_rng(SEED), 500)
    ok1, ok2 = a["p_one_sided_within_group"] < 0.01, b["p_one_sided_within_group"] > 0.05
    print(f"[self-test] injected dose-response: rho={a['spearman_rho']} p={a['p_one_sided_within_group']}  "
          f"{'PASS' if ok1 else 'FAIL'}", file=sys.stderr)
    print(f"[self-test] no relation:            rho={b['spearman_rho']} p={b['p_one_sided_within_group']}  "
          f"{'PASS' if ok2 else 'FAIL'}", file=sys.stderr)
    return ok1 and ok2


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--all", action="store_true", help="run on real data (default action)")
    args = ap.parse_args()
    if args.self_test:
        sys.exit(0 if self_test() else 1)
    if not self_test():
        print("[main] self-test FAILED -- aborting before touching real data", file=sys.stderr)
        sys.exit(1)
    run_all()
