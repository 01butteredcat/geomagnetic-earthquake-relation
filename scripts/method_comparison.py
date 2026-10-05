"""跨方法比較：夜間平均值殘差 vs. 日變幅比值。

問的是：在 surrogate_test.py 使用的同一套 2000 組替代資料檢定下，
兩種日尺度異常指標中哪一個給出較小的經驗 p 值：

- **夜間平均值殘差**——compute_indices.py 既有的 `local_anomaly_index_{H,F}`
  （近站池夜間平均 MAD z，遠站池的共模由
  Theil-Sen 迴歸扣除）。向量組用 H，純量組用 F，
  和 compute_indices.py 實際跑的測站池一致。
- **日變幅比值**（Liu et al., 2006）——這裡新加的。總磁場的
  每日變幅，dB_i = 一個 UTC 日內的 max(F_i) - min(F_i)，以及
  R_or = dB_o / dB_r，o 是震央附近的觀測站、
  r 是距離 o 不超過 100 km 的參考站。平靜條件下，
  全球擾動對兩站的影響相同，R 大致維持常數；
  局部地殼導電度的變化會把它往任一方向推離。
  檢定的序列是 log(R)，所以增加和減少是對稱的。

測站配對（每組固定，由錨點震央決定）：o = 測站池中
最近的測站；r = 在距 o 100 km 以內、且比 o 離
震央更遠的測站中，離震央最遠的那一個。如果 o
本身離震央超過 100 km，該組標記為 low_confidence。

向量站的 F 回報的是佔位值，所以用
1 分鐘序列（minute_series_<sta>.parquet）重建為 sqrt(X^2+Y^2+Z^2)。
分鐘平均會略微低估 1Hz 的變幅，但兩站低估的程度
相同，所以比值不受影響。

兩種指標都經過完全相同的檢定，只用乾淨日（不是磁暴
日、兩站都沒有中斷——和 compute_indices.py 的
is_clean_day 標準相同），串接方式和 surrogate_test.py 一樣。統計量是
雙尾最大 |z|——比值可能往任一方向偏離，而 compute_indices.py 的
候選旗標也是用 |index|：

  (B) 震前窗口（主要結果）：錨點前 30 天內乾淨日的
      最大 |z|（cross_group_analysis.pre_event_window），
      對照只用背景日建立的留一窗虛無分布——
      區塊 bootstrap 和相位隨機化。為什麼窗口必須排除在
      自己的虛無分布之外，見 surrogate_test()。
  (A) 整條序列：整條序列 MAD z-score 的最大 |z|，只用相位
      隨機化。和地震時間無關，而且已知相位
      隨機化在厚尾資料上偏寬鬆；只當作
      surrogate_test.py 旁邊的參考點報告。

各 N_SURROGATES 組，p = (1 + #{null >= obs}) / (N + 1)。跨
組：每個指標 p < 0.05 的組數，對 5% 做二項檢定，以及
兩個指標之間 -log10(p) 的配對 Wilcoxon 符號等級檢定。

用法：
  method_comparison.py --self-test
  method_comparison.py --group G10
  method_comparison.py --all
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest, wilcoxon

sys.path.insert(0, str(Path(__file__).parent))
import common  # noqa: E402
from cross_group_analysis import pre_event_window  # noqa: E402
from events import GROUPS  # noqa: E402
from stat_utils import (  # noqa: E402
    block_bootstrap_surrogate,
    estimate_block_length,
    mad_zscore,
    phase_randomize_surrogate,
)

N_SURROGATES = 2000
SEED = 20260805
PAIR_MAX_KM = 100.0          # Liu et al. (2006)：兩站應該靠近
PRE_WINDOW_DAYS = 30
MIN_MINUTES_PER_DAY = 1296   # 1440 的 90%；有大缺口的日子會低估它的變幅
MIN_SERIES_DAYS = 20         # 和 surrogate_test.py 相同的下限
MIN_WINDOW_DAYS = 5          # 震前乾淨日少於這個就不嘗試 (B)
# 安慰劑校準：在每 PLACEBO_STEP_DAYS 天一個、30 天窗口
# 不和真實窗口重疊的假錨點上重跑 (B)。校準良好的檢定在
# 約 5% 的假錨點上給出 p < 0.05；在這份資料上留一窗虛無分布偏寬鬆（區塊
# bootstrap 約 13-15%，相位隨機化約 22-25%，2026-09-24），所以跨組
# 二項檢定是和安慰劑比例比較，而不是 5%。
PLACEBO_STEP_DAYS = 7
N_PLACEBO_SURROGATES = 500
PLACEBO_PRIOR_WINDOWS = 10   # 往合併比例收縮的強度，見 shrunk_placebo_rate()   # 每個安慰劑窗口；只用來算合併比例，所以抽樣次數少一點就夠
INDICES = ("night_residual", "range_ratio")
# 實際跑的（統計量, 替代資料）組合——為什麼整條序列統計量
# 沒有區塊 bootstrap 版本，見 surrogate_test()
TESTS = (("window", "block_bootstrap"), ("window", "phase_randomization"), ("whole", "phase_randomization"))

OUT_DIR = common.PROJECT_DIR / "data" / "interim" / "method_comparison"


# ---------------------------------------------------------------------------
# 日變幅比值
# ---------------------------------------------------------------------------

def pick_station_pair(cfg) -> dict:
    pool = cfg.xyz_pool if cfg.xyz_pool.sufficient else cfg.f_pool
    ev = cfg.anchor_event
    st = cfg.stations
    d_epi = {s: common.haversine_km(st[s]["lat"], st[s]["lon"], ev.lat, ev.lon) for s in pool.all_stations}
    obs = min(d_epi, key=d_epi.get)
    cands = [s for s in d_epi if s != obs and d_epi[s] > d_epi[obs]
             and common.haversine_km(st[s]["lat"], st[s]["lon"], st[obs]["lat"], st[obs]["lon"]) <= PAIR_MAX_KM]
    if not cands:
        return {"error": f"no reference station within {PAIR_MAX_KM:.0f} km of {obs}"}
    ref = max(cands, key=d_epi.get)
    return {
        "field": "H" if pool.channel == "XYZ" else "F",
        "obs_station": obs, "ref_station": ref,
        "obs_epicenter_km": round(d_epi[obs], 1), "ref_epicenter_km": round(d_epi[ref], 1),
        "pair_km": round(common.haversine_km(st[ref]["lat"], st[ref]["lon"], st[obs]["lat"], st[obs]["lon"]), 1),
        "low_confidence": bool(d_epi[obs] > PAIR_MAX_KM),
    }


def daily_f_range(cfg, station: str) -> pd.Series:
    """從 1 分鐘序列算每個 UTC 日的 max(F) - min(F)，有效分鐘數
    少於 MIN_MINUTES_PER_DAY 的日子為 NaN。"""
    m = pd.read_parquet(cfg.interim_dir / f"minute_series_{station}.parquet")
    f = m["F"] if "F" in m.columns else np.sqrt(m["X"] ** 2 + m["Y"] ** 2 + m["Z"] ** 2)
    day = f.index.strftime("%Y%m%d")
    g = f.groupby(day)
    rng = (g.max() - g.min()).where(g.count() >= MIN_MINUTES_PER_DAY)
    return rng.rename(station)


def range_ratio_series(cfg, pair: dict, storm_dates: set[str], outages: dict[str, set[str]]) -> pd.Series:
    """乾淨日上的 log(dB_obs / dB_ref)，以 YYYYMMDD 為索引。"""
    o, r = pair["obs_station"], pair["ref_station"]
    df = pd.concat([daily_f_range(cfg, o), daily_f_range(cfg, r)], axis=1).dropna()
    df = df[(df[o] > 0) & (df[r] > 0)]
    dirty = storm_dates | outages.get(o, set()) | outages.get(r, set())
    df = df[~df.index.isin(dirty)]
    return np.log(df[o] / df[r]).sort_index()


def night_residual_series(cfg, field: str) -> pd.Series:
    df = pd.read_csv(cfg.interim_dir / "local_anomaly_index.csv", dtype={"date": str})
    col = f"local_anomaly_index_{field}"
    df = df[df["is_clean_day"].astype(bool) & df[col].notna()]
    return df.set_index("date")[col].sort_index()


# ---------------------------------------------------------------------------
# 共用的替代資料檢定
# ---------------------------------------------------------------------------

def _extreme(z: np.ndarray, tail: str) -> float:
    """方向調整成越大 = 越極端的檢定統計量：雙尾檢定用 max|z|，
    下尾檢定（低谷）用 -min(z)。"""
    return float(np.nanmax(np.abs(z))) if tail == "two" else float(-np.nanmin(z))


def _p(null: np.ndarray, obs: float) -> float:
    null = null[~np.isnan(null)]
    return float((1 + np.sum(null >= obs)) / (len(null) + 1)) if len(null) else float("nan")


def surrogate_test(values: np.ndarray, window_pos: np.ndarray, rng: np.random.Generator,
                   n_surr: int | None = None, tail: str = "two") -> dict:
    """(B) 震前窗口的 max|z| 對照留一窗替代資料，加上 (A)
    整條序列的 max|z| 對照相位隨機化替代資料。window_pos
    是 `values` 的索引。

    (B) 完全只用背景（非窗口）日建立：它們的
    中位數／MAD 同時用來標準化觀測窗口和每一組替代資料，
    每組替代資料都是從背景重抽出的一段窗口長度序列。
    如果改成重抽整條觀測序列（像 surrogate_test.py 那樣），
    被檢定的極值會漏進它自己的虛無分布——區塊 bootstrap 在
    本專案的序列長度下就會有大約 0.2-0.4 的 p 值下限，不管
    異常多大（2026-09-24 已用模擬確認）。
    (A) 沒有另外的背景可以重抽，所以區塊 bootstrap 在那裡
    結構上無效，只報告相位隨機化。

    tail="two" 檢定 max|z|（這支腳本的預設）；tail="lower" 檢定
    最負的 z，對應 surrogate_test.py 原本的「震前低谷」問題。
    回報的 obs_* 值相應地是 max|z| 或 min z。"""
    sign = 1.0 if tail == "two" else -1.0  # 把 _extreme() 轉回回報用的 min z
    n_surr = n_surr or N_SURROGATES
    out: dict = {"n_days": int(len(values)), "n_window_days": int(len(window_pos))}

    # (A) 整條序列，只用相位隨機化
    z, _, _ = mad_zscore(values)
    obs_a = _extreme(z, tail)
    null_a = np.full(n_surr, np.nan)
    for i in range(n_surr):
        zs, _, mad = mad_zscore(phase_randomize_surrogate(values, rng))
        if mad > 1e-9:
            null_a[i] = _extreme(zs, tail)
    out["obs_whole_max_abs_z" if tail == "two" else "obs_whole_min_z"] = round(sign * obs_a, 4)
    out["p_whole_phase_randomization"] = round(_p(null_a, obs_a), 5)

    # (B) 震前窗口 vs. 只用背景的虛無分布
    in_win = np.zeros(len(values), dtype=bool)
    in_win[window_pos] = True
    bg = values[~in_win]
    w = int(in_win.sum())
    _, med, mad = mad_zscore(bg)
    if w < MIN_WINDOW_DAYS or len(bg) < max(MIN_SERIES_DAYS, w) or mad <= 1e-9:
        out.update({"n_background_days": int(len(bg)),
                    "obs_window_max_abs_z" if tail == "two" else "obs_window_min_z": None,
                    "p_window_block_bootstrap": None, "p_window_phase_randomization": None})
        return out
    obs_b = _extreme((values[in_win] - med) / mad, tail)
    block_len = estimate_block_length(bg)
    out.update({"n_background_days": int(len(bg)), "block_length_days": block_len,
                "obs_window_max_abs_z" if tail == "two" else "obs_window_min_z": round(sign * obs_b, 4)})
    makers = {"block_bootstrap": lambda: block_bootstrap_surrogate(bg, block_len, rng),
              "phase_randomization": lambda: phase_randomize_surrogate(bg, rng)}
    for name, make in makers.items():
        null_b = np.empty(n_surr)
        for i in range(n_surr):
            s = make()
            start = rng.integers(0, len(s) - w + 1)
            null_b[i] = _extreme((s[start:start + w] - med) / mad, tail)
        out[f"p_window_{name}"] = round(_p(null_b, obs_b), 5)
    return out


def run_group(group_id: str) -> dict:
    cfg = common.load_group_config(group_id)
    row: dict = {"group": group_id, "anchor": cfg.anchor_event.date}
    pair = pick_station_pair(cfg)
    if "error" in pair:
        return {**row, "error": pair["error"]}
    row.update(pair)

    daily = pd.read_csv(cfg.interim_dir / "daily_features.csv", dtype={"date": str})
    outages = common.auto_outage_dates(daily)
    storm = pd.read_csv(cfg.interim_dir / "storm_days.csv", dtype={"date": str})
    storm_dates = set(storm["date"])
    window = set(pre_event_window(pd.to_datetime(cfg.anchor_event.date), PRE_WINDOW_DAYS))

    series = {"night_residual": night_residual_series(cfg, pair["field"]),
              "range_ratio": range_ratio_series(cfg, pair, storm_dates, outages)}
    for name, s in series.items():
        rng = np.random.default_rng(SEED)  # 每個指標用相同的抽樣，讓兩者在同樣基礎上檢定
        if len(s) < MIN_SERIES_DAYS:
            row[f"{name}__error"] = f"only {len(s)} clean days"
            continue
        values = s.to_numpy(dtype=float)
        pos = np.flatnonzero(s.index.isin(window))
        for k, v in surrogate_test(values, pos, rng).items():
            row[f"{name}__{k}"] = v
        for k, v in placebo_counts(s, window, rng).items():
            row[f"{name}__{k}"] = v
    return row


def placebo_counts(s: pd.Series, real_window: set[str], rng: np.random.Generator, tail: str = "two") -> dict:
    """在不和真實震前窗口重疊的假錨點上跑 (B)，
    計算每種替代資料給出 p < 0.05 的頻率。"""
    dates = pd.to_datetime(s.index)
    values = s.to_numpy(dtype=float)
    anchors = pd.date_range(dates.min() + pd.Timedelta(days=PRE_WINDOW_DAYS + 1), dates.max(), freq=f"{PLACEBO_STEP_DAYS}D")
    n = 0
    hits = {"block_bootstrap": 0, "phase_randomization": 0}
    for a in anchors:
        win = set(pre_event_window(a, PRE_WINDOW_DAYS))
        if win & real_window:
            continue
        r = surrogate_test(values, np.flatnonzero(s.index.isin(win)), rng, N_PLACEBO_SURROGATES, tail)
        if r["p_window_block_bootstrap"] is None:
            continue
        n += 1
        for sur in hits:
            hits[sur] += int(r[f"p_window_{sur}"] < 0.05)
    return {"placebo_n": n, **{f"placebo_n_lt05_{sur}": k for sur, k in hits.items()}}


def _window_extreme(values: np.ndarray, center: float, tail: str) -> float:
    """越大 = 越極端：下尾用 -min，雙尾用 max|x - center|。"""
    return float(-values.min()) if tail == "lower" else float(np.abs(values - center).max())


def rank_window_test(s: pd.Series, real_window: set[str], tail: str = "lower") -> dict:
    """把真實震前窗口的極值，和同一序列中其他 PRE_WINDOW_DAYS 天
    區段排名，這些區段從真實窗口往前、往後鋪排，讓任兩個
    區段都不重疊（每段需要 >= MIN_WINDOW_DAYS 個乾淨日）。在 H0 下
    各區段可交換，所以 p = (1 + #{fake >= obs}) / (n_fake + 1) 是精確的，
    只差區段邊界之間的逐日相關。

    取代 (B) 的重抽 p 值，成為主要的震前檢定（2026-09-25）：
    區塊 bootstrap 只能重用背景值，所以只要窗口的
    最小值低於背景最小值（G11 pc3：z = -1.2 vs -0.73），它的 p
    就落在下限 1/(N+1)，不管差距多小——而在 H0 下發生這種事的機率
    大約是 w/(w+n)，不是 0.0005。滑動（重疊）的假窗口也修不好：
    一個低值日會同時落在最多 30 個假窗口裡，所以真實窗口
    排第一的頻率遠高於 1/(n+1)（改成鋪排之前，
    rank_self_test 在名目 5% 下是 14%）。誠實的代價：150 天的序列
    只有約 4 個假區段，所以單一組不可能低於 p = 0.2；顯著性
    必須來自跨組合併（fisher_across_groups）。

    `null_ps` 是每個區段對照其他所有區段的 p（包含真實區段）：
    也就是 fisher_across_groups() 抽樣用的精確逐組虛無分布。"""
    dates = pd.to_datetime(s.index)
    values = s.to_numpy(dtype=float)
    idx = s.index.to_numpy()
    center = float(np.median(values))

    def stat(win: set[str]) -> float | None:
        v = values[np.isin(idx, list(win))]
        return _window_extreme(v, center, tail) if len(v) >= MIN_WINDOW_DAYS else None

    obs = stat(real_window)
    if obs is None:
        return {"rank_p": None, "rank_n_fake_windows": 0}
    anchor = pd.to_datetime(max(real_window)) + pd.Timedelta(days=1)
    step = pd.Timedelta(days=PRE_WINDOW_DAYS)
    fake = []
    for direction in (-1, 1):
        k = 1
        while True:
            a = anchor + direction * k * step
            if a - step > dates.max() or a <= dates.min():
                break
            v = stat(set(pre_event_window(a, PRE_WINDOW_DAYS)))
            if v is not None:
                fake.append(v)
            k += 1
    if not fake:
        return {"rank_p": None, "rank_n_fake_windows": 0}
    everything = np.append(np.array(fake), obs)
    all_ps = [(1 + np.sum(np.delete(everything, j) >= everything[j])) / len(everything) for j in range(len(everything))]
    return {"rank_p": round(float(all_ps[-1]), 5),
            "rank_n_fake_windows": len(fake),
            "rank_min_attainable_p": round(1 / len(everything), 5),
            "null_ps": [round(float(p), 5) for p in all_ps]}


def fisher_across_groups(ps: list[float], null_ps: list[list[float]], rng: np.random.Generator,
                         n_sim: int = 20000) -> dict:
    """Fisher 的 X = -2 sum log p，跨組加總，校準方式是每組抽一個假
    窗口的 p（所以短序列的離散性直接放進
    虛無分布，而不是假設 p 是均勻分布）。"""
    obs = float(-2 * np.sum(np.log(ps)))
    sims = np.zeros(n_sim)
    for nps in null_ps:
        sims += -2 * np.log(rng.choice(np.asarray(nps), n_sim))
    return {"n_groups": len(ps), "fisher_x": round(obs, 3),
            "p": round(float((1 + np.sum(sims >= obs)) / (n_sim + 1)), 5),
            "n_rank_p_lt_05": int(sum(p < 0.05 for p in ps)),
            "n_rank_p_lt_10": int(sum(p < 0.10 for p in ps))}


# ---------------------------------------------------------------------------
# 跨組摘要
# ---------------------------------------------------------------------------

def poisson_binomial_sf(k: int, ps: list[float]) -> float:
    """K 為獨立 Bernoulli(p_i) 之和時的 P(K >= k)。"""
    dist = np.array([1.0])
    for p in ps:
        dist = np.convolve(dist, [1 - p, p])
    return float(dist[k:].sum())


def shrunk_placebo_rate(k: int, n: int, pooled: float) -> float:
    """一組自己的安慰劑 p<0.05 比例 k/n，用
    PLACEBO_PRIOR_WINDOWS 個虛擬窗口往合併比例收縮。合併比例本身由
    長序列（G6/G7/G8）主導，可能低估短序列組的偽陽性
    比例；原始的 k/n 在一組只有 0-5 個安慰劑窗口時又太吵。"""
    return (k + PLACEBO_PRIOR_WINDOWS * pooled) / (n + PLACEBO_PRIOR_WINDOWS)


def summarize(rows: list[dict]) -> dict:
    df = pd.DataFrame(rows)
    out: dict = {"n_groups": len(df), "n_surrogates": N_SURROGATES, "pre_window_days": PRE_WINDOW_DAYS, "tests": {}}
    for stat, sur in TESTS:
        if True:
            key = f"{stat}_{sur}"
            entry: dict = {}
            cols = {ix: f"{ix}__p_{stat}_{sur}" for ix in INDICES}
            for ix, col in cols.items():
                p = df[col].dropna() if col in df else pd.Series(dtype=float)
                k = int((p < 0.05).sum())
                entry[ix] = {"n_groups": int(len(p)), "n_p_lt_05": k,
                             "median_p": round(float(p.median()), 4) if len(p) else None}
                pn, pk = f"{ix}__placebo_n", f"{ix}__placebo_n_lt05_{sur}"
                if stat == "window" and pn in df and len(p):
                    n_pl, k_pl = int(df[pn].fillna(0).sum()), int(df[pk].fillna(0).sum())
                    rate = k_pl / n_pl if n_pl else None
                    tested = df[df[col].notna()]
                    own = [shrunk_placebo_rate(int(r[pk]) if r[pk] == r[pk] else 0,
                                               int(r[pn]) if r[pn] == r[pn] else 0, rate or 0.05)
                           for r in tested.to_dict("records")]
                    entry[ix].update({
                        "mean_per_group_placebo_rate": round(float(np.mean(own)), 4) if own else None,
                        "poisson_binomial_p_vs_own_placebo": round(poisson_binomial_sf(k, own), 5) if own else None,
                        "placebo_n_windows": n_pl, "placebo_rate_p_lt_05": None if rate is None else round(rate, 4),
                        "binom_p_vs_placebo_rate": None if not rate else
                            round(binomtest(k, len(p), rate, alternative="greater").pvalue, 5),
                    })
            both = df[[cols["night_residual"], cols["range_ratio"]]].dropna() if all(c in df for c in cols.values()) else pd.DataFrame()
            if len(both) >= 6:
                a = -np.log10(both[cols["night_residual"]].to_numpy())
                b = -np.log10(both[cols["range_ratio"]].to_numpy())
                w = wilcoxon(a, b, zero_method="wilcox") if np.any(a != b) else None
                entry["paired"] = {
                    "n_groups": int(len(both)),
                    "n_night_residual_smaller_p": int(np.sum(a > b)),
                    "n_range_ratio_smaller_p": int(np.sum(b > a)),
                    "wilcoxon_p_two_sided": None if w is None else round(float(w.pvalue), 5),
                }
            out["tests"][key] = entry
    return out


def render_markdown(rows: list[dict], summ: dict) -> str:
    L = ["# 跨方法比較：夜間平均值殘差 vs 全天日變幅比值（Liu et al., 2006）", "",
         f"替代序列各 {N_SURROGATES} 組（區塊拔靴、相位隨機化）；統計量為雙尾最大 |z|；"
         f"(B) 震前窗口為主震前 {PRE_WINDOW_DAYS} 天。只用乾淨日（非磁暴、兩站皆無斷線）。", "",
         "(B) 的虛無分布只用震前窗口以外的背景期資料建構（留一窗），避免被檢定的極值混進自己的虛無分布。",
         f"**校準**：同樣的檢定在不與真實窗口重疊、每 {PLACEBO_STEP_DAYS} 天一個的假主震日上重跑（安慰劑）。"
         "校準良好時安慰劑 p<0.05 的比例應約 5%；實際比例較高，表示虛無分布偏寬鬆。"
         "跨組檢定讓每組對照自己的安慰劑比例（Poisson-二項檢定）；表中「安慰劑比例」為各組比例的平均。", "",
         "## 跨組彙總", "",
         "| 統計量 | 替代序列 | 夜間殘差 p<0.05 組數 | 夜間殘差安慰劑比例 | 日變幅比值 p<0.05 組數 | 日變幅比值安慰劑比例 | 夜間殘差較小 p 的組數 | 比值較小 p 的組數 | Wilcoxon p（雙尾） |",
         "|---|---|---|---|---|---|---|---|---|"]
    label = {"window": "(B) 震前 30 天", "whole": "(A) 整條序列"}
    for key, e in summ["tests"].items():
        stat, sur = key.split("_", 1)
        nr, rr, pr = e["night_residual"], e["range_ratio"], e.get("paired", {})
        def cnt(e):
            b = e.get("poisson_binomial_p_vs_own_placebo")
            return f"{e['n_p_lt_05']}/{e['n_groups']}" + (f"（p={b}）" if b is not None else "")
        def rate(e):
            r = e.get("mean_per_group_placebo_rate")
            return "未校準" if r is None else f"{r:.1%}"
        L.append(f"| {label[stat]} | {sur} | {cnt(nr)} | {rate(nr)} | {cnt(rr)} | {rate(rr)} "
                 f"| {pr.get('n_night_residual_smaller_p', '—')} | {pr.get('n_range_ratio_smaller_p', '—')} "
                 f"| {pr.get('wilcoxon_p_two_sided', '—')} |")
    L += ["", "## 各組結果（(B) 震前窗口）", "",
          "| 組 | 場 | 觀測站 o（距震央） | 參考站 r（距震央） | 兩站距離 | 夜間殘差 p（拔靴／相位） | 日變幅比值 p（拔靴／相位） | 備註 |",
          "|---|---|---|---|---|---|---|---|"]
    for r in rows:
        if "error" in r:
            L.append(f"| {r['group']} | — | — | — | — | — | — | {r['error']} |")
            continue
        def pp(ix):
            if f"{ix}__error" in r:
                return r[f"{ix}__error"]
            a, b = r.get(f"{ix}__p_window_block_bootstrap"), r.get(f"{ix}__p_window_phase_randomization")
            return "窗口內乾淨日不足" if a is None else f"{a} / {b}"
        note = "o 距震央 >100 km，低可信度" if r["low_confidence"] else ""
        L.append(f"| {r['group']} | {r['field']} | {r['obs_station']}（{r['obs_epicenter_km']} km） "
                 f"| {r['ref_station']}（{r['ref_epicenter_km']} km） | {r['pair_km']} km "
                 f"| {pp('night_residual')} | {pp('range_ratio')} | {note} |")
    L += ["", "## 已知限制", "",
          "- 兩個指標用的測站不同：夜間殘差是近站池（3 站中位數）對遠站池回歸的殘差，日變幅比值只用一對站。比較的是兩套完整方法，不是同一批站上的兩種統計量。",
          "- 向量組的夜間殘差用 H，日變幅比值依 Liu (2006) 用總磁場 F，場分量不同。",
          "- 向量站的 F 由每分鐘平均的 X、Y、Z 重建，日變幅略小於 1 秒資料，但兩站一致，比值不受影響。",
          "- (A) 整條序列統計量跟地震時間無關，只檢定極值能否由序列本身結構解釋，而且相位隨機化對厚尾資料偏寬鬆、沒有校準；僅供對照，結論以 (B) 為準。",
          "- (B) 的留一窗虛無分布在安慰劑窗口上偏寬鬆（見上表），可能反映序列的季節漂移等非平穩性；校準後的二項檢定才有意義。",
          "- 共 2 種統計量 × 2 種替代序列 × 2 個指標 × 最多 24 組，未做多重比較校正。", ""]
    return "\n".join(L)


# ---------------------------------------------------------------------------
# 自我測試
# ---------------------------------------------------------------------------

def self_test() -> bool:
    rng = np.random.default_rng(SEED)
    n, w = 150, np.arange(110, 140)
    noise = np.cumsum(rng.normal(0, 0.3, n)) * 0.1 + rng.normal(0, 1, n)  # 輕微的自相關
    shifted = noise.copy()
    shifted[w[10:14]] += 8.0

    global N_SURROGATES
    saved, N_SURROGATES = N_SURROGATES, 500
    try:
        sig = surrogate_test(shifted, w, np.random.default_rng(SEED))
        nul = surrogate_test(noise, w, np.random.default_rng(SEED))
        spike = noise.copy()
        spike[w[15]] += 100.0
        leak = surrogate_test(spike, w, np.random.default_rng(SEED))
    finally:
        N_SURROGATES = saved
    ok1 = sig["p_window_block_bootstrap"] < 0.05 and sig["p_window_phase_randomization"] < 0.05
    ok2 = nul["p_window_block_bootstrap"] > 0.05 and nul["p_window_phase_randomization"] > 0.05
    # 被檢定的極值不能漏進它自己的虛無分布：窗口中 100 個 sigma 的突波
    # 必須達到 p 值下限 1/(N+1)
    ok3 = leak["p_window_block_bootstrap"] <= 2 / 501 and leak["p_window_phase_randomization"] <= 2 / 501
    print(f"[self-test] injected shift: p_window = {sig['p_window_block_bootstrap']} / "
          f"{sig['p_window_phase_randomization']}  {'PASS' if ok1 else 'FAIL'}", file=sys.stderr)
    print(f"[self-test] noise only:     p_window = {nul['p_window_block_bootstrap']} / "
          f"{nul['p_window_phase_randomization']}  {'PASS' if ok2 else 'FAIL'}", file=sys.stderr)
    print(f"[self-test] 100σ spike (leak check): p_window = {leak['p_window_block_bootstrap']} / "
          f"{leak['p_window_phase_randomization']}  {'PASS' if ok3 else 'FAIL'}", file=sys.stderr)
    return ok1 and ok2 and ok3 and rank_self_test()


def rank_self_test(n_rep: int = 100, n_groups: int = 12) -> bool:
    """在 12 個合成的 150 天組別上跑 rank_window_test + fisher_across_groups
    （窗口結束在第 93 天，本專案典型的下載配置）。在 H0 下，
    合併 p 應該在約 5% 的重複中 < 0.05；它所取代的留一窗
    區塊 bootstrap 逐組列出供參考（印出，不做
    斷言）。在一半組別注入的低谷必須被偵測到。單一組別
    無法檢定檢定力：只有約 4 個假區段，它們的 p 不可能低於 0.2。"""
    global N_SURROGATES
    rng = np.random.default_rng(SEED)
    dates = pd.date_range("2024-01-01", periods=150).strftime("%Y%m%d")
    window = set(pre_event_window(pd.Timestamp("2024-04-03"), PRE_WINDOW_DAYS))  # 第 93 天 = 2024-04-03
    w_pos = np.flatnonzero(np.isin(dates, list(window)))

    def series() -> np.ndarray:
        return np.cumsum(rng.normal(0, 0.3, 150)) * 0.1 + rng.normal(0, 1, 150)

    def combined(dip_groups: int) -> float:
        ps, nulls = [], []
        for g in range(n_groups):
            x = series()
            if g < dip_groups:
                x[w_pos[12]] -= 8.0
            r = rank_window_test(pd.Series(x, index=dates), window)
            ps.append(r["rank_p"])
            nulls.append(r["null_ps"])
        return fisher_across_groups(ps, nulls, rng, n_sim=5000)["p"]

    fp = float(np.mean([combined(0) < 0.05 for _ in range(n_rep)]))
    saved, N_SURROGATES = N_SURROGATES, 200
    try:
        bb = float(np.mean([surrogate_test(series(), w_pos, rng, tail="lower")["p_window_block_bootstrap"] < 0.05
                            for _ in range(n_rep)]))
    finally:
        N_SURROGATES = saved
    power = combined(n_groups // 2)
    ok1, ok2 = fp <= 0.09, power < 0.01
    print(f"[self-test] {n_groups} null groups, combined rank p<0.05 in {fp:.1%} of {n_rep}  "
          f"{'PASS' if ok1 else 'FAIL'} (per-group leave-window-out block bootstrap: {bb:.1%}, for reference)",
          file=sys.stderr)
    print(f"[self-test] dip in {n_groups // 2}/{n_groups} groups: combined rank p = {power}  "
          f"{'PASS' if ok2 else 'FAIL'}", file=sys.stderr)
    return ok1 and ok2


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--group", action="append", dest="groups")
    ap.add_argument("--all", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        sys.exit(0 if self_test() else 1)
    if not args.groups and not args.all:
        ap.error("pass --group <id> (repeatable) or --all")
    if not self_test():
        print("[main] self-test FAILED -- aborting before touching real data", file=sys.stderr)
        sys.exit(1)

    group_ids = list(GROUPS) if args.all else args.groups
    rows = []
    for gid in group_ids:
        row = run_group(gid)
        rows.append(row)
        if "error" in row:
            print(f"[{gid}] {row['error']}", file=sys.stderr)
        else:
            print(f"[{gid}] {row['obs_station']}/{row['ref_station']} "
                  f"night p_window(bb)={row.get('night_residual__p_window_block_bootstrap')} "
                  f"ratio p_window(bb)={row.get('range_ratio__p_window_block_bootstrap')}", file=sys.stderr)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summ = summarize(rows)
    pd.DataFrame(rows).to_csv(OUT_DIR / "per_group.csv", index=False)
    (OUT_DIR / "summary.json").write_text(json.dumps(summ, indent=2, ensure_ascii=False))
    (OUT_DIR / "summary.md").write_text(render_markdown(rows, summ), encoding="utf-8")
    print(f"wrote {OUT_DIR}", file=sys.stderr)


if __name__ == "__main__":
    main()
