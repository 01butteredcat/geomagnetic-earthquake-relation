"""ULF 近站／遠站極化
差值序列（`pc3_diff_zh` / `pc4_diff_zh`）的替代資料顯著性檢定，依教授的
第一個建議：既然「原始 z-score 看起來像常態」的直覺並不
嚴格適用於有自相關、非常態的地磁序列，就把
序列本身打亂（區塊 bootstrap 和 FFT 相位隨機化），得到許多
統計上看起來一樣、但和任何地震都沒有關係的
替代序列，再問像觀測到的
低谷那麼極端的情況（或像 report_template.html 對 G10 2024-03-30
引用的特定 z <= -4.1 門檻那麼極端）有多常偶然出現。

對每個有 ulf_near_far_index.csv 的組別逐（組, 頻帶）執行
（8 個向量站足夠的組別；為什麼
G1/G2/G3 不在範圍內見 stat_utils.py／計畫檔）。透過 stat_utils.mad_zscore 重用
cross_group_analysis.py::ulf_candidate_dates 已經建立的
整條序列中位數／MAD z-score 公式，所以這裡的結果可以直接
和那個既有分析比較（也是它的嚴謹度升級版）。

## 震前窗口檢定（2026-09-24 加入）

上面的整條序列檢定對區塊 bootstrap 有結構性缺陷：它的
替代序列是用觀測序列自己的區塊重組的，所以被檢定的
極值會漏進它自己的虛無分布。在本專案的序列長度（115-310
天）下，這讓區塊 bootstrap 的 p 值有大約 0.2-0.4 的下限，不管
低谷多大（已用模擬確認；這裡每組的區塊 bootstrap p
都落在 0.23-0.42，例如 G13 的最小 z = -31 仍得到 p = 0.35）。那些欄位
為了 report_template_validation.html 原樣保留，但區塊 bootstrap 在那裡的
「不顯著」不能當成任何事情的證據。

因此每個結果也帶有 `pre_event_window`：該組錨點前
30 天內最負的 z，只用背景（非窗口）日來標準化
和檢定——透過 method_comparison.surrogate_test(tail="lower") 做留一窗區塊 bootstrap 和
相位隨機化——再加上
每 7 天一個、不重疊的假錨點做安慰劑校準。window_summary.{json,md} 中的
跨組二項檢定是和安慰劑比例比較，而不是 5%，因為留一窗虛無分布
在這份資料上本身就偏寬鬆。


用法：
  surrogate_test.py --group G10          # 單一組，兩個頻帶
  surrogate_test.py --all                # 每個有 ULF 資料的組
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
import method_comparison as mc  # noqa: E402
from common import PROJECT_DIR, load_group_config  # noqa: E402
from cross_group_analysis import pre_event_window  # noqa: E402
from scipy.stats import binomtest  # noqa: E402
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
SEED = 20260805  # 固定種子，重跑這支腳本會得到相同的 p 值


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

    # 1. 常態性檢查——這正是一開始改用中位數／MAD + 替代資料，
    # 而不用參數式常態理論檢定的理由。
    shapiro_w, shapiro_p = shapiro(vals)

    # 2. 觀測統計量：整條序列中位數／MAD z-score，最負的那一天。
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

    window = set(pre_event_window(pd.to_datetime(cfg.anchor_event.date), mc.PRE_WINDOW_DAYS))
    # 這裡去掉磁暴／恢復期日（上面的舊檢定沒有）：如果保留它們，
    # G9 和 G13 的震前最小值都落在磁暴日，好幾個震前窗口
    # 的磁暴也比背景期多。乾淨日規則和 method_comparison.py 相同。
    storm = set(pd.read_csv(cfg.interim_dir / "storm_days.csv", dtype={"date": str})["date"])
    keep = ~np.isin(dates, list(storm))
    series = pd.Series(vals[keep], index=dates[keep])
    # 用自己的產生器，所以加上這一節不會改變舊的整條序列抽樣（因而
    # report_template_validation.html 的數字也一個位元都不變）
    rng_w = np.random.default_rng([SEED, sum(map(ord, group_id + band))])
    if len(series) < 20:
        return {"group": group_id, "band": band, "error": f"only {len(series)} non-storm days"}
    win_res = mc.surrogate_test(series.to_numpy(dtype=float), np.flatnonzero(series.index.isin(window)),
                                rng_w, tail="lower")
    pre_event = {k: v for k, v in win_res.items() if "whole" not in k}
    pre_event.update(mc.placebo_counts(series, window, rng_w, tail="lower"))
    pre_event.update(mc.rank_window_test(series, window, tail="lower"))
    # 安慰劑：緊接在真實窗口前的區段假裝成真實窗口——如果
    # 光是位在錨點正前方就足以排到最低，這個也會「顯著」
    prev = mc.rank_window_test(series, set(pre_event_window(pd.to_datetime(cfg.anchor_event.date)
                                                             - pd.Timedelta(days=mc.PRE_WINDOW_DAYS),
                                                             mc.PRE_WINDOW_DAYS)), tail="lower")
    pre_event["placebo_prev_block"] = {k: prev.get(k) for k in ("rank_p", "rank_n_fake_windows", "null_ps")}
    pre_event.update({"anchor": cfg.anchor_event.date, "window_days": mc.PRE_WINDOW_DAYS,
                      "storm_days_excluded": True,
                      "null": "leave-window-out (background days only)", "tail": "lower (most negative z)"})

    floor_note = ("whole-series block bootstrap resamples the observed series itself, so the tested extreme "
                  "leaks into its own null: p has a structural floor (~0.2-0.4 here). Not evidence; see "
                  "pre_event_window.")
    return {
        "group": group_id,
        "band": band,
        "pre_event_window": pre_event,
        "n_days": n,
        "shapiro": {"W": round(float(shapiro_w), 5), "p": round(float(shapiro_p), 6),
                    "rejects_normality_p05": bool(shapiro_p < 0.05)},
        "observed": {
            "median": round(med, 5), "mad_scaled": round(mad, 5),
            "min_z": round(t_obs, 4), "min_z_date": obs_date,
            "n_days_z_le_fixed_threshold": n_days_le_threshold,
        },
        "fixed_rule_threshold": FIXED_RULE_THRESHOLD,
        "block_bootstrap": summarize(block_mins, "block_bootstrap", {"block_length_days": block_len,
                                                                     "structural_floor_note": floor_note}),
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

    results = []
    for group_id in groups:
        rng = np.random.default_rng(SEED)  # 每組重設，讓頻帶的順序不影響結果
        for band in BANDS:
            result = run_one(group_id, band, rng)
            if result is None:
                print(f"[{group_id}/{band}] no data, skipped", file=sys.stderr)
                continue
            out_path = out_dir / f"{group_id}_{band}.json"
            out_path.write_text(json.dumps(result, indent=2))
            results.append(result)
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
                pe = result["pre_event_window"]
                print(f"[{group_id}/{band}] pre-event window: min_z={pe.get('obs_window_min_z')} "
                      f"p(bb)={pe.get('p_window_block_bootstrap')} p(phase)={pe.get('p_window_phase_randomization')} "
                      f"placebo n={pe['placebo_n']}", file=sys.stderr)

    if args.all:
        summ = window_summary(results)
        (out_dir / "window_summary.json").write_text(json.dumps(summ, indent=2, ensure_ascii=False))
        (out_dir / "window_summary.md").write_text(render_window_summary(results, summ), encoding="utf-8")
        print(f"wrote {out_dir / 'window_summary.md'}", file=sys.stderr)


# 2026-09-25 在看到排名檢定數字之前就預先指定：主要的
# 震前檢定是 PRIMARY_BAND + 排名檢定 + 所有組別；其他都是
# 探索性的。敏感度分析（S3）：共用原始資料夾的組別並不
# 獨立，所以每個資料夾家族只保留規模最大的錨點。
PRIMARY_BAND = "pc3"
FAMILY_DROP = ("G6", "G7", "G24")  # G6/G7/G8 保留 G8（ML6.7），G23/G24 保留 G23（ML6.24）


def rank_summary(results: list[dict], rng: np.random.Generator) -> dict:
    out = {}
    for band in BANDS:
        rows = [(r["group"], r["pre_event_window"]) for r in results
                if r.get("band") == band and r.get("pre_event_window", {}).get("rank_p") is not None]
        entry = {}
        for name, keep in (("all_groups", lambda g: True), ("one_per_family", lambda g: g not in FAMILY_DROP)):
            sel = [(g, pe) for g, pe in rows if keep(g)]
            entry[name] = {**mc.fisher_across_groups([pe["rank_p"] for _, pe in sel],
                                                     [pe["null_ps"] for _, pe in sel], rng),
                           "groups": [g for g, _ in sel]}
        entry["sensitivity"] = rank_sensitivity(rows, np.random.default_rng([SEED, sum(map(ord, band)), 1]))
        out[band] = entry
    return out


def rank_sensitivity(rows: list[tuple[str, dict]], rng: np.random.Generator) -> dict:
    """主要 Fisher 合併的檢查（用自己的 rng，所以主要數字
    不會動）：G10 是最早發現的案例（頻帶、方向、窗口和
    -4.1 門檻都來自它 2024-03-30 的低谷），所以拿掉它；同時拿掉 G10 和
    G11；輪流拿掉每一組；以及緊接在
    每個真實窗口前的安慰劑區段，用同樣方式合併。"""
    def fisher(sel):
        return mc.fisher_across_groups([pe["rank_p"] for _, pe in sel], [pe["null_ps"] for _, pe in sel], rng)["p"]
    out = {name: {"groups": [g for g, _ in sel], "p": fisher(sel)}
           for name, sel in (("drop_discovery_G10", [(g, pe) for g, pe in rows if g != "G10"]),
                             ("drop_G10_G11", [(g, pe) for g, pe in rows if g not in ("G10", "G11")]))}
    out["leave_one_out"] = {g: fisher([(h, pe) for h, pe in rows if h != g]) for g, _ in rows}
    prev = [(g, pe["placebo_prev_block"]) for g, pe in rows
            if (pe.get("placebo_prev_block") or {}).get("rank_p") is not None]
    out["placebo_prev_block"] = {"groups": [g for g, _ in prev], "p": fisher(prev) if prev else None}
    return out


def window_summary(results: list[dict]) -> dict:
    """主要：rank_summary（每組做排名檢定，跨組用 Fisher 合併）。
    已被取代但保留供比較：逐頻帶 x 替代資料，震前
    p < 0.05 的組數，和安慰劑比例比較。"""
    out: dict = {"window_days": mc.PRE_WINDOW_DAYS, "tail": "lower",
                 "primary": {"band": PRIMARY_BAND, "test": "rank test, Fisher across all groups"},
                 "rank_test": rank_summary(results, np.random.default_rng(SEED)), "bands": {}}
    for band in BANDS:
        rows = [r["pre_event_window"] for r in results if r.get("band") == band and "pre_event_window" in r]
        entry = {}
        for sur in ("block_bootstrap", "phase_randomization"):
            tested = [r for r in rows if r.get(f"p_window_{sur}") is not None]
            ps = [r[f"p_window_{sur}"] for r in tested]
            k = sum(p < 0.05 for p in ps)
            n_pl = sum(r["placebo_n"] for r in rows)
            k_pl = sum(r[f"placebo_n_lt05_{sur}"] for r in rows)
            rate = k_pl / n_pl if n_pl else None
            # 每組和自己的安慰劑比例比較，再往合併比例收縮
            own = [mc.shrunk_placebo_rate(r[f"placebo_n_lt05_{sur}"], r["placebo_n"], rate or 0.05) for r in tested]
            entry[sur] = {
                "n_groups": len(ps), "n_p_lt_05": k,
                "placebo_n_windows": n_pl, "placebo_rate_p_lt_05": None if rate is None else round(rate, 4),
                "binom_p_vs_placebo_rate": None if not rate or not ps else
                    round(binomtest(k, len(ps), rate, alternative="greater").pvalue, 5),
                "mean_per_group_placebo_rate": round(float(np.mean(own)), 4) if own else None,
                "poisson_binomial_p_vs_own_placebo": round(mc.poisson_binomial_sf(k, own), 5) if own else None,
            }
        out["bands"][band] = entry
    return out


def render_window_summary(results: list[dict], summ: dict) -> str:
    L = ["# ULF 替代資料檢定：震前窗口版", "",
         "## 主要結果：排名檢定（2026-09-25 起）", "",
         "每組把序列切成不重疊的 30 天區塊（對齊真窗口），看真窗口的最低值在所有區塊中排第幾：p =（1 + 比它更低的假區塊數）÷（區塊總數）。"
         "沒有地震效應時，每個區塊一樣可能最低，所以這個 p 是準確的；代價是短序列的單組 p 無法小於 1 ÷ 區塊數。"
         "跨組用 Fisher 合併（−2 Σ ln p），虛無分布由「每組隨機挑一個區塊當真窗口」模擬 20000 次。",
         f"事先指定的主要檢定：{summ['primary']['band']}、全部組別。其餘都是探索性結果。"
         "敏感度：共用原始資料的家族各只留規模最大的一組（G6/G7/G8 留 G8，G23/G24 留 G23）。", "",
         "| 頻帶 | 組別 | 組數 | 單組 p<0.05 | 單組 p<0.10 | Fisher X | 合併 p |", "|---|---|---|---|---|---|---|"]
    for band, e in summ["rank_test"].items():
        for name, v in e.items():
            if name == "sensitivity":
                continue
            label = "全部" if name == "all_groups" else "每家族一組"
            star = "（主要）" if band == summ["primary"]["band"] and name == "all_groups" else ""
            L.append(f"| {band}{star} | {label} | {v['n_groups']} | {v['n_rank_p_lt_05']} | {v['n_rank_p_lt_10']} "
                     f"| {v['fisher_x']} | {v['p']} |")
    L += ["", "### 主要檢定的敏感度（2026-09-28）", "",
          "G10 是發現個案：頻帶、下凹方向、30 天窗口和 −4.1 門檻都來自它 2024-03-30 的下凹，所以另外算不含 G10 的版本。"
          "安慰劑區塊：把每組真窗口「前一個」30 天區塊當成真窗口，用同樣方法合併；若只靠「緊鄰主震前」這個位置就會顯著，這裡也會顯著。", "",
          "| 頻帶 | 不含 G10 | 不含 G10、G11 | 安慰劑區塊 | 逐一拿掉一組時最大的 p |", "|---|---|---|---|---|"]
    for band, e in summ["rank_test"].items():
        sv = e["sensitivity"]
        loo = sv["leave_one_out"]
        worst = max(loo, key=loo.get) if loo else None
        L.append(f"| {band} | {sv['drop_discovery_G10']['p']} | {sv['drop_G10_G11']['p']} | {sv['placebo_prev_block']['p']} "
                 f"| {'—' if worst is None else f'{loo[worst]}（拿掉 {worst}）'} |")
    L += ["", "## 舊版（已被排名檢定取代，保留比較）：留一窗＋安慰劑校準", "",
         "留一窗的區塊拔靴只能重組背景期已有的數值：震前窗口最小值只要低於背景期最小值，p 就固定是下限 0.0005，"
         "不管低多少（例如 G11 Pc3 的 z = −1.2）。下表的 p 值因此偏小，不能當成證據。", "",
         f"統計量：主震前 {summ['window_days']} 天內最負的 z（單尾，檢定「震前下凹」）。"
         "z 以震前窗口以外的背景期中位數／MAD 標準化；虛無分布只由背景期重抽樣（區塊拔靴、相位隨機化各 2000 組）。",
         "只用非磁暴日（磁暴與恢復期日排除，與夜間殘差、日變幅比值相同）。",
         "安慰劑：每 7 天一個、不與真實窗口重疊的假主震日，重跑同一檢定；校準良好時 p<0.05 的比例應約 5%。"
         "合併比例由長序列的 G6/G7/G8 主導，會低估短序列組的偽陽性率，所以主要結果改用每組各自的安慰劑比例（Poisson-二項檢定）。", "",
         "舊版整條序列檢定的區塊拔靴 p 值有結構下限（約 0.2–0.4），不能當成證據；見 surrogate_test.py 說明。", "",
         "## 跨組彙總", "",
         "| 頻帶 | 替代序列 | 震前 p<0.05 組數 | 合併安慰劑比例 | 二項檢定 p（對合併比例） | 各組安慰劑比例平均 | Poisson-二項 p（各組對自己的比例，主要結果） |",
         "|---|---|---|---|---|---|---|"]
    for band, e in summ["bands"].items():
        for sur, v in e.items():
            rate = "—" if v["placebo_rate_p_lt_05"] is None else f"{v['placebo_rate_p_lt_05']:.1%}（{v['placebo_n_windows']} 窗）"
            L.append(f"| {band} | {sur} | {v['n_p_lt_05']}/{v['n_groups']} | {rate} | {v['binom_p_vs_placebo_rate']} "
                     f"| {v['mean_per_group_placebo_rate']} | {v['poisson_binomial_p_vs_own_placebo']} |")
    L += ["", "## 各組結果", "",
          "| 組 | 頻帶 | 震前最低 z | 排名 p（區塊數） | p（區塊拔靴，舊） | p（相位隨機化，舊） | 窗口／背景天數 | 安慰劑窗數 | 整條序列 p（拔靴／相位，舊） |",
          "|---|---|---|---|---|---|---|---|---|"]
    for r in results:
        if "pre_event_window" not in r:
            continue
        pe = r["pre_event_window"]
        n_blocks = pe.get("rank_n_fake_windows", 0) + 1 if pe.get("rank_p") is not None else "—"
        L.append(f"| {r['group']} | {r['band']} | {pe.get('obs_window_min_z')} | {pe.get('rank_p')}（{n_blocks}） "
                 f"| {pe.get('p_window_block_bootstrap')} "
                 f"| {pe.get('p_window_phase_randomization')} | {pe['n_window_days']}／{pe.get('n_background_days')} "
                 f"| {pe['placebo_n']} | {r['block_bootstrap']['p_value_vs_observed_extreme']}／"
                 f"{r['phase_randomization']['p_value_vs_observed_extreme']} |")
    L += ["", "## 限制", "",
          "- 主要檢定只有 1 個（Pc3、全部組別），不需要校正；其餘組合是探索性結果，若要一起看，應以 Bonferroni 調整。",
          "- 排名檢定假設不同 30 天區塊可以互換：若序列有長期漂移，最低值會偏向序列兩端，而真窗口多半在中後段，結果會偏保守。",
          "- G6/G7/G8、G23/G24 共用原始資料，不是完全獨立的樣本。",
          "- 震前窗口與背景期的磁暴比例不同時，排除磁暴日後可用天數也不同，窗口內可能只剩少數幾天。", ""]
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    main()
