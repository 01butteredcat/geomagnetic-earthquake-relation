"""跨組重現率分析（計畫 §7）：地磁
前兆訊號在 13 個獨立地震事件組之間重複出現的頻率，
是否高於各組自己的內部機率所預期的？

這明確是**重現率**分析，不是正式的統計
驗證。依照本專案既定的方針，
真正有統計結論的檢定需要約 20-30 起獨立 M>=6 事件，加上
配對的平靜期對照；13 組相對於最初 n=1 的個案研究
是有意義的擴充，但仍未達到那個標準。這裡的發現
應該解讀為「在我們手上的組別中這重複出現了幾次」，
絕不是「這證明了前兆存在」。

方法，逐組：
  1. 每個適用方法的候選日期（向量站足夠的組別用 compute_indices.py 的 H/Z，
     只有純量的組別用 F，再加上在這裡計算的 ULF
     Pc3 Z/H 差值候選旗標——ulf_analysis.py
     只產生原始的近站／遠站／差值序列，不像 compute_indices.py 對 H/Z/F
     那樣自己套用 MAD z-score 門檻，所以
     這一個額外步驟在這支腳本裡做，使用和
     compute_indices.py 相同的 CANDIDATE_Z_THRESHOLD，對
     該組整個乾淨日分布套用一次——對已經做過差分的 ULF 序列，
     並沒有像原始 H/Z/F 水準那樣明顯的逐日
     因果「滑動窗口」概念，所以改用單一的整條序列 MAD
     z-score，而不是滾動的）。
  2. 對三個震前窗口（緊接在該組錨點事件前的 7/14/30 天）
     各自檢查是否有**任何**候選日期落在裡面。三種窗口長度
     並列報告，而不是
     事後挑看起來最好的那個（避免 look-elsewhere
     偏誤）。
  3. 每組提供**自己的**虛無／基準率：「該組資料中隨機一個 W 天窗口
     含有 >=1 個候選日」的比例，計算方式是讓
     W 天窗口滑過該組自己日期範圍內每一個可能的起始日
     （排除實際的震前檢定窗口），取
     經驗命中比例。這用該組自己的候選日數量
     特性當對照，所以不需要另外下載平靜期資料集
     （依 adaptive-tulip 計畫，那仍是未來真正、更大的一步，
     這裡沒有嘗試）。

跨組彙總分成兩個不同的信心層級報告（絕不
合併）：有足夠 XYZ 測站池的組別（可以跑 H/Z 和
ULF）vs. 只能用 F 方法的組別（G1、G2、G3）——依 common.py 的 GroupConfig，
兩者是結構上不同、可靠度也不同的方法。

輸出：
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
    """對 pc3_diff_zh 做整條序列 MAD z-score 候選標記——為什麼這和
    compute_indices.py 的滾動滑動窗口做法不同，見
    模組 docstring。"""
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
    """在 all_dates 內所有可能的 window_days 天滑動窗口中，
    含有 >= 1 個候選日期的比例（排除每一個和實際檢定
    窗口共用任何一天的窗口——只跳過完全相同的窗口，會讓被檢定窗口自己的
    候選日拉高它的基準率）。回傳
    (rate, n_windows_checked)。"""
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
    """跨組彙總命中率 vs 基準率，依層級和
    方法拆開，絕不把向量層和純量層的組別合併。
    用 scipy.stats.binomtest 對每個層級／方法／窗口自己的平均
    基準率做簡單的顯著性摘要——3 種窗口長度 x
    好幾種方法代表同時跑了好幾個檢定，**沒有**做多重
    比較校正；這是當作探索性的重現率
    訊號報告，不是確認的偵測。"""
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
    # 自 2026-09-20 拆分舊的合併鍵 G2_G3 / G6_G7_G8 之後，
    # 一個 GROUPS 鍵 == 一個編號的地震組（G1 ... G23）。
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
