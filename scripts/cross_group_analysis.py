"""Cross-group reproducibility analysis (plan §7): does a geomagnetic
precursor signal recur across the 13 independent earthquake event-groups
more often than each group's own internal chance rate would predict?

This is explicitly a REPRODUCIBILITY-RATE analysis, not a formal statistical
validation. Per this project's established guidance, a
real statistically-conclusive test needs ~20-30 independent M>=6 events plus
matched quiet-period controls; 13 groups is a meaningful expansion beyond
the original n=1 case study but still falls short of that bar. Findings here
should be read as "how often did this recur across the groups we have",
never as "this proves a precursor exists".

Method, per group:
  1. Candidate dates for each applicable method (H/Z from compute_indices.py
     for vector-sufficient groups, F for scalar-only groups, plus a ULF
     Pc3 Z/H differential candidate flag computed here -- ulf_analysis.py
     only produces the raw near/far/diff series, it does not itself apply a
     MAD z-score threshold the way compute_indices.py does for H/Z/F, so
     that one additional step is done in this script, using the same
     CANDIDATE_Z_THRESHOLD compute_indices.py uses, applied once over the
     group's whole clean-day distribution -- there's no obvious per-day
     causal "trailing window" concept for an already-differenced ULF series
     the way there is for the raw H/Z/F level, so a single whole-series MAD
     z-score is used instead of a rolling one).
  2. For each of three pre-event windows (7/14/30 days immediately before
     the group's anchor event), check whether ANY candidate date falls in
     it. All three window sizes are reported side by side rather than
     picking whichever looks best after the fact (avoiding look-elsewhere
     bias).
  3. Each group supplies its OWN null/base rate for "a random W-day window
     in this group's data contains >=1 candidate", computed by sliding a
     W-day window across every possible start day in the group's own date
     range (excluding the actual pre-event test window) and taking the
     empirical hit fraction. This uses the group's own candidate-count
     characteristics as its control, so no separate quiet-period dataset
     needs to be fetched (that remains a real, larger future step per the
     adaptive-tulip plan, not attempted here).

Cross-group aggregation is reported in two separate confidence tiers (never
pooled together): groups with a sufficient XYZ station pool (can run H/Z and
ULF) vs. groups limited to the F-only method (G1, G2, G3) -- structurally
different methods with different reliability, per common.py's GroupConfig.

Outputs:
  data/interim/cross_group_summary.json
  data/interim/cross_group_summary.md
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest

sys.path.insert(0, str(Path(__file__).parent))
from common import PROJECT_DIR, load_group_config  # noqa: E402
from events import ALL_GROUP_IDS  # noqa: E402
from compute_indices import CANDIDATE_Z_THRESHOLD  # noqa: E402

WINDOWS_DAYS = [7, 14, 30]


def ulf_candidate_dates(cfg) -> list[str]:
    """Whole-series MAD z-score candidate flagging on pc3_diff_zh -- see
    module docstring for why this differs from compute_indices.py's rolling
    trailing-window approach."""
    path = cfg.interim_dir / "ulf_near_far_index.csv"
    if not path.exists():
        return []
    df = pd.read_csv(path, dtype={"date": str})
    vals = df["pc3_diff_zh"].dropna()
    if len(vals) < 10:
        return []
    med = vals.median()
    mad = (vals - med).abs().median() * 1.4826
    if mad <= 1e-9:
        return []
    z = (df["pc3_diff_zh"] - med) / mad
    flagged = df.loc[z.abs() > CANDIDATE_Z_THRESHOLD, "date"]
    return sorted(flagged.tolist())


def pre_event_window(anchor_date: pd.Timestamp, window_days: int) -> list[str]:
    return pd.date_range(anchor_date - pd.Timedelta(days=window_days), anchor_date - pd.Timedelta(days=1)) \
        .strftime("%Y%m%d").tolist()


def sliding_baseline_rate(all_dates: list[str], candidate_dates: set[str], window_days: int, exclude_window: set[str]) -> tuple[float, int]:
    """Fraction of every possible window_days-long sliding window within
    all_dates (excluding every window that shares a day with the actual test
    window -- skipping only the identical window let the tested window's own
    candidates raise its baseline) that contains >= 1 candidate date. Returns
    (rate, n_windows_checked)."""
    if len(all_dates) < window_days + 1:
        return float("nan"), 0
    hits = 0
    n = 0
    for i in range(len(all_dates) - window_days + 1):
        w = all_dates[i:i + window_days]
        if exclude_window.intersection(w):
            continue
        n += 1
        if candidate_dates.intersection(w):
            hits += 1
    return (hits / n if n else float("nan")), n


def analyze_group(group_id: str) -> dict:
    cfg = load_group_config(group_id)
    cand = json.loads((cfg.interim_dir / "candidate_windows.json").read_text())
    local_idx = pd.read_csv(cfg.interim_dir / "local_anomaly_index.csv", dtype={"date": str})
    all_dates = sorted(local_idx["date"].unique())
    anchor = cfg.anchor_event

    methods = {}
    if cfg.xyz_pool.sufficient:
        methods["H"] = set(cand.get("candidate_dates_H", []))
        methods["Z"] = set(cand.get("candidate_dates_Z", []))
        methods["ULF_pc3"] = set(ulf_candidate_dates(cfg))
    if cfg.f_pool.sufficient:
        methods["F"] = set(cand.get("candidate_dates_F", []))

    anchor_ts = pd.to_datetime(anchor.date)
    tier = "vector" if cfg.xyz_pool.sufficient else "scalar_only"

    windows_out = {}
    for w in WINDOWS_DAYS:
        test_window = pre_event_window(anchor_ts, w)
        test_window_set = set(test_window)
        per_method = {}
        for name, dates in methods.items():
            hit = bool(dates.intersection(test_window_set))
            baseline_rate, n_windows = sliding_baseline_rate(all_dates, dates, w, test_window_set)
            per_method[name] = {
                "n_candidates_total": len(dates),
                "hit_in_pre_event_window": hit,
                "baseline_rate": None if baseline_rate != baseline_rate else round(baseline_rate, 4),
                "n_baseline_windows": n_windows,
            }
        windows_out[str(w)] = {
            "test_window": [test_window[0], test_window[-1]] if test_window else None,
            "methods": per_method,
        }

    return {
        "group": group_id,
        "tier": tier,
        "anchor_date": anchor.date,
        "anchor_magnitude": f"{anchor.magnitude_type}{anchor.magnitude}",
        "coord_confidence": anchor.coord_confidence,
        "methods_available": sorted(methods.keys()),
        "n_days_in_range": len(all_dates),
        "windows": windows_out,
    }


def aggregate(per_group: list[dict]) -> dict:
    """Aggregate hit-rate vs baseline-rate across groups, split by tier and
    method, never pooling vector-tier and scalar-only-tier groups together.
    Uses scipy.stats.binomtest against each tier/method/window's own mean
    baseline rate as a simple significance summary -- 3 window sizes x
    several methods means several tests are being run with NO multiple-
    comparison correction; this is reported as an exploratory reproducibility
    signal, not a confirmed detection."""
    out = {}
    tiers = sorted({g["tier"] for g in per_group})
    for tier in tiers:
        tier_groups = [g for g in per_group if g["tier"] == tier]
        all_methods = sorted({m for g in tier_groups for m in g["methods_available"]})
        out[tier] = {"n_groups": len(tier_groups), "groups": [g["group"] for g in tier_groups], "methods": {}}
        for method in all_methods:
            out[tier]["methods"][method] = {}
            for w in WINDOWS_DAYS:
                hits, baselines, groups_with_method = [], [], []
                for g in tier_groups:
                    m = g["windows"][str(w)]["methods"].get(method)
                    if m is None or m["baseline_rate"] is None:
                        continue
                    hits.append(1 if m["hit_in_pre_event_window"] else 0)
                    baselines.append(m["baseline_rate"])
                    groups_with_method.append(g["group"])
                n = len(hits)
                if n == 0:
                    out[tier]["methods"][method][str(w)] = {"n_groups": 0, "note": "no group had a computable baseline for this window"}
                    continue
                observed_hits = sum(hits)
                mean_baseline = float(np.mean(baselines))
                p_value = None
                if 0 < mean_baseline < 1:
                    p_value = binomtest(observed_hits, n, mean_baseline, alternative="greater").pvalue
                out[tier]["methods"][method][str(w)] = {
                    "n_groups": n,
                    "groups": groups_with_method,
                    "observed_hits": observed_hits,
                    "observed_hit_rate": round(observed_hits / n, 4),
                    "mean_group_baseline_rate": round(mean_baseline, 4),
                    "binomtest_p_greater": round(p_value, 4) if p_value is not None else None,
                }
    return out


def render_markdown(per_group: list[dict], agg: dict) -> str:
    # One GROUPS key == one numbered earthquake group (G1 ... G23) since the
    # 2026-09-20 split of the old merged G2_G3 / G6_G7_G8 keys.
    n_groups = len(ALL_GROUP_IDS)
    group_range = f"{ALL_GROUP_IDS[0]}-{ALL_GROUP_IDS[-1]}"
    lines = []
    lines.append(f"# 跨組重現率分析（{group_range}）\n")
    threshold_note = (
        f"首次觸及正式統計檢定文獻常引用的 20-30 起門檻下緣（見 "
        f"`1999-9-21-2010-3-4-2015-2-14-2016-2-6-2-adaptive-tulip.md` 第1節）"
        if n_groups >= 20
        else f"仍遠低於正式統計檢定建議的 20-30 起門檻（見 "
        f"`1999-9-21-2010-3-4-2015-2-14-2016-2-6-2-adaptive-tulip.md` 第1節）"
    )
    lines.append(
        f"**重要限制**：{n_groups} 組獨立序列{threshold_note}。這是樣本數量門檻，不代表資料品質或獨立性"
        "已達標——部分組別震央座標信心較低、部分較舊組別測站覆蓋不完整（見 events.py 的 "
        "`coord_confidence` 與各組 verification_report.json）。"
        "本報告呈現的是「重現率」——地磁異常訊號在幾組中於震前窗口出現——"
        "不是具正式統計顯著性的地震前兆存在證據。三種窗口天數（7/14/30 天）並列呈現，"
        "不是挑選事後看起來最漂亮的單一結果。\n"
    )
    for tier, label in [("vector", "向量法可跑的組別（H/Z 篩選 + ULF 極化）"), ("scalar_only", "僅純量法可跑的組別（F 篩選）")]:
        if tier not in agg:
            continue
        lines.append(f"\n## {label}：{', '.join(agg[tier]['groups'])}\n")
        for method, by_window in agg[tier]["methods"].items():
            lines.append(f"\n### 方法：{method}\n")
            lines.append("| 震前窗口 | 組數 | 出現候選異常的組數 | 觀察重現率 | 各組平均基準率（組內 permutation） | 二項檢定 p (單尾) |")
            lines.append("|---|---|---|---|---|---|")
            for w in WINDOWS_DAYS:
                r = by_window.get(str(w), {})
                if r.get("n_groups", 0) == 0:
                    lines.append(f"| {w} 天 | 0 | - | - | - | - |")
                    continue
                lines.append(
                    f"| {w} 天 | {r['n_groups']} | {r['observed_hits']} | {r['observed_hit_rate']} | "
                    f"{r['mean_group_baseline_rate']} | {r['binomtest_p_greater']} |"
                )
    lines.append("\n## 各組明細\n")
    lines.append("| 組別 | 信心層級 | 錨定事件 | 座標信心 | 可用方法 |")
    lines.append("|---|---|---|---|---|")
    for g in per_group:
        lines.append(f"| {g['group']} | {g['tier']} | {g['anchor_date']} {g['anchor_magnitude']} | {g['coord_confidence']} | {', '.join(g['methods_available'])} |")
    return "\n".join(lines) + "\n"


def main():
    per_group = [analyze_group(g) for g in ALL_GROUP_IDS]
    agg = aggregate(per_group)

    out_dir = PROJECT_DIR / "data" / "interim"
    (out_dir / "cross_group_summary.json").write_text(json.dumps({"per_group": per_group, "aggregate": agg}, indent=2))
    (out_dir / "cross_group_summary.md").write_text(render_markdown(per_group, agg))
    print(f"wrote {out_dir / 'cross_group_summary.json'}", file=sys.stderr)
    print(f"wrote {out_dir / 'cross_group_summary.md'}", file=sys.stderr)
    print(render_markdown(per_group, agg))


if __name__ == "__main__":
    main()
