"""Cross-method comparison: night-mean residual vs. diurnal range ratio.

Asks which of two daily-scale anomaly indices gives the smaller empirical
p-value under the same 2000-surrogate test surrogate_test.py uses:

- **Night-mean residual** -- the existing `local_anomaly_index_{H,F}` from
  compute_indices.py (near-pool night mean MAD z, far-pool common mode
  removed by Theil-Sen regression). H for vector groups, F for scalar ones,
  matching which pool compute_indices.py actually ran.
- **Diurnal range ratio** (Liu et al., 2006) -- new here. Daily range of the
  total field, dB_i = max(F_i) - min(F_i) over one UTC day, and
  R_or = dB_o / dB_r between an observation station o near the epicenter and
  a reference station r no more than 100 km from o. Under quiet conditions
  global disturbances hit both stations alike and R stays near constant; a
  local crustal conductivity change would push it off in either direction.
  The series tested is log(R), so increases and decreases are symmetric.

Station pair (fixed per group, from the anchor epicenter): o = the pool's
nearest station; r = among stations within 100 km of o and farther from the
epicenter than o, the one farthest from the epicenter. Groups whose o is
itself more than 100 km from the epicenter are flagged low_confidence.

F for vector stations is reported as a placeholder, so it's rebuilt from
the 1-minute series (minute_series_<sta>.parquet) as sqrt(X^2+Y^2+Z^2).
Minute means slightly understate the 1Hz range, but identically at both
stations, so the ratio is unaffected.

Both indices go through the identical test, on clean days only (not a storm
day, neither station in an outage -- the same criteria as compute_indices.py's
is_clean_day), concatenated the way surrogate_test.py does. Statistics are
two-sided max |z| -- the ratio can deviate either way and compute_indices.py's
candidate flag also uses |index|:

  (B) pre-event window (the main result): max |z| over the clean days in the
      30 days before the anchor (cross_group_analysis.pre_event_window),
      against a leave-window-out null built from the background days only --
      block-bootstrap and phase-randomized. See surrogate_test() for why the
      window must be left out of its own null.
  (A) whole series: max |z| of the whole-series MAD z-score, phase
      randomization only. Not tied to earthquake timing, and phase
      randomization is known to be liberal on heavy-tailed data; reported
      only as a reference point next to surrogate_test.py.

N_SURROGATES each, p = (1 + #{null >= obs}) / (N + 1). Across
groups: count of p < 0.05 per index with a binomial test against 5%, and a
paired Wilcoxon signed-rank test on -log10(p) between the two indices.

Usage:
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
PAIR_MAX_KM = 100.0          # Liu et al. (2006): the two stations should be close
PRE_WINDOW_DAYS = 30
MIN_MINUTES_PER_DAY = 1296   # 90% of 1440; a day with a big gap would understate its range
MIN_SERIES_DAYS = 20         # same floor as surrogate_test.py
MIN_WINDOW_DAYS = 5          # fewer clean pre-event days than this and (B) isn't attempted
# Placebo calibration: re-run (B) at fake anchors every PLACEBO_STEP_DAYS whose
# 30-day window doesn't overlap the real one. A calibrated test gives p < 0.05 at
# ~5% of them; on this data the leave-window-out null is liberal (~13-15% block
# bootstrap, ~22-25% phase randomization, 2026-09-24), so the cross-group
# binomial tests compare against the placebo rate, not 5%.
PLACEBO_STEP_DAYS = 7
N_PLACEBO_SURROGATES = 500
PLACEBO_PRIOR_WINDOWS = 10   # shrinkage strength toward the pooled rate, see shrunk_placebo_rate()   # per placebo window; only feeds a pooled rate, so fewer draws suffice
INDICES = ("night_residual", "range_ratio")
# (statistic, surrogate) combinations actually run -- see surrogate_test() for why
# the whole-series statistic has no block-bootstrap version
TESTS = (("window", "block_bootstrap"), ("window", "phase_randomization"), ("whole", "phase_randomization"))

OUT_DIR = common.PROJECT_DIR / "data" / "interim" / "method_comparison"


# ---------------------------------------------------------------------------
# Diurnal range ratio
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
    """Per-UTC-day max(F) - min(F) from the 1-minute series, NaN for days
    with fewer than MIN_MINUTES_PER_DAY valid minutes."""
    m = pd.read_parquet(cfg.interim_dir / f"minute_series_{station}.parquet")
    f = m["F"] if "F" in m.columns else np.sqrt(m["X"] ** 2 + m["Y"] ** 2 + m["Z"] ** 2)
    day = f.index.strftime("%Y%m%d")
    g = f.groupby(day)
    rng = (g.max() - g.min()).where(g.count() >= MIN_MINUTES_PER_DAY)
    return rng.rename(station)


def range_ratio_series(cfg, pair: dict, storm_dates: set[str], outages: dict[str, set[str]]) -> pd.Series:
    """log(dB_obs / dB_ref) on clean days, indexed by YYYYMMDD."""
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
# Shared surrogate test
# ---------------------------------------------------------------------------

def _extreme(z: np.ndarray, tail: str) -> float:
    """Test statistic oriented so that larger = more extreme: max|z| for a
    two-sided test, -min(z) for a lower-tail one (a dip)."""
    return float(np.nanmax(np.abs(z))) if tail == "two" else float(-np.nanmin(z))


def _p(null: np.ndarray, obs: float) -> float:
    null = null[~np.isnan(null)]
    return float((1 + np.sum(null >= obs)) / (len(null) + 1)) if len(null) else float("nan")


def surrogate_test(values: np.ndarray, window_pos: np.ndarray, rng: np.random.Generator,
                   n_surr: int | None = None, tail: str = "two") -> dict:
    """(B) pre-window max|z| against leave-window-out surrogates, plus (A)
    whole-series max|z| against phase-randomized surrogates. window_pos
    indexes into `values`.

    (B) builds everything from the background (non-window) days only: their
    median/MAD standardize both the observed window and every surrogate, and
    each surrogate is a window-length stretch resampled from the background.
    Resampling the whole observed series instead (as surrogate_test.py does)
    lets the tested extreme leak into its own null -- block bootstrap then has
    a p-value floor of roughly 0.2-0.4 at this project's series lengths, no
    matter how large the anomaly (checked by simulation 2026-09-24).
    (A) has no separate background to resample, so block bootstrap is
    structurally invalid there and only phase randomization is reported.

    tail="two" tests max|z| (this script's default); tail="lower" tests the
    most negative z, for surrogate_test.py's original "pre-event dip" question.
    Reported obs_* values are max|z| or min z accordingly."""
    sign = 1.0 if tail == "two" else -1.0  # converts _extreme() back to the reported min z
    n_surr = n_surr or N_SURROGATES
    out: dict = {"n_days": int(len(values)), "n_window_days": int(len(window_pos))}

    # (A) whole series, phase randomization only
    z, _, _ = mad_zscore(values)
    obs_a = _extreme(z, tail)
    null_a = np.full(n_surr, np.nan)
    for i in range(n_surr):
        zs, _, mad = mad_zscore(phase_randomize_surrogate(values, rng))
        if mad > 1e-9:
            null_a[i] = _extreme(zs, tail)
    out["obs_whole_max_abs_z" if tail == "two" else "obs_whole_min_z"] = round(sign * obs_a, 4)
    out["p_whole_phase_randomization"] = round(_p(null_a, obs_a), 5)

    # (B) pre-event window vs. background-only null
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
        rng = np.random.default_rng(SEED)  # same draws per index, so the two are tested on equal footing
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
    """Run (B) at fake anchors that don't overlap the real pre-event window and
    count how often each surrogate type gives p < 0.05."""
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
    """Larger = more extreme: -min for the lower tail, max|x - center| for two."""
    return float(-values.min()) if tail == "lower" else float(np.abs(values - center).max())


def rank_window_test(s: pd.Series, real_window: set[str], tail: str = "lower") -> dict:
    """Rank the real pre-event window's extreme among the other PRE_WINDOW_DAYS
    blocks of the same series, tiled back and forward from the real window so no
    two blocks overlap (each needs >= MIN_WINDOW_DAYS clean days). Under H0 the
    blocks are exchangeable, so p = (1 + #{fake >= obs}) / (n_fake + 1) is exact
    up to the day-to-day correlation across block edges.

    Replaces (B)'s resampling p-value as the primary pre-event test (2026-09-25):
    block bootstrap can only reuse background values, so whenever the window's
    minimum is below the background minimum (G11 pc3: z = -1.2 vs -0.73) its p
    sits at the floor 1/(N+1) however small the gap -- the chance of that under
    H0 is about w/(w+n), not 0.0005. Sliding (overlapping) fake windows don't fix
    it either: one low day then sits in up to 30 of them, so the real window
    ranks first far more often than 1/(n+1) (14 % at nominal 5 % in
    rank_self_test before the switch to tiling). Honest cost: a 150-day series
    has ~4 fake blocks, so no single group can go below p = 0.2; significance
    has to come from combining groups (fisher_across_groups).

    `null_ps` is every block's p against all the others (real block included):
    the exact per-group null distribution fisher_across_groups() draws from."""
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
    """Fisher's X = -2 sum log p over groups, calibrated by drawing one fake
    window's p per group (so the discreteness of short series is built into the
    null instead of assuming uniform p)."""
    obs = float(-2 * np.sum(np.log(ps)))
    sims = np.zeros(n_sim)
    for nps in null_ps:
        sims += -2 * np.log(rng.choice(np.asarray(nps), n_sim))
    return {"n_groups": len(ps), "fisher_x": round(obs, 3),
            "p": round(float((1 + np.sum(sims >= obs)) / (n_sim + 1)), 5),
            "n_rank_p_lt_05": int(sum(p < 0.05 for p in ps)),
            "n_rank_p_lt_10": int(sum(p < 0.10 for p in ps))}


# ---------------------------------------------------------------------------
# Cross-group summary
# ---------------------------------------------------------------------------

def poisson_binomial_sf(k: int, ps: list[float]) -> float:
    """P(K >= k) for K a sum of independent Bernoulli(p_i)."""
    dist = np.array([1.0])
    for p in ps:
        dist = np.convolve(dist, [1 - p, p])
    return float(dist[k:].sum())


def shrunk_placebo_rate(k: int, n: int, pooled: float) -> float:
    """A group's own placebo p<0.05 rate k/n, shrunk toward the pooled rate with
    PLACEBO_PRIOR_WINDOWS pseudo-windows. The pooled rate alone is dominated by
    long series (G6/G7/G8) and can understate a short group's false-positive
    rate; the raw k/n is too noisy when a group has only 0-5 placebo windows."""
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
# Self-test
# ---------------------------------------------------------------------------

def self_test() -> bool:
    rng = np.random.default_rng(SEED)
    n, w = 150, np.arange(110, 140)
    noise = np.cumsum(rng.normal(0, 0.3, n)) * 0.1 + rng.normal(0, 1, n)  # mild autocorrelation
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
    # the tested extreme must not leak into its own null: a 100-sigma spike in the
    # window has to reach the p-value floor, 1/(N+1)
    ok3 = leak["p_window_block_bootstrap"] <= 2 / 501 and leak["p_window_phase_randomization"] <= 2 / 501
    print(f"[self-test] injected shift: p_window = {sig['p_window_block_bootstrap']} / "
          f"{sig['p_window_phase_randomization']}  {'PASS' if ok1 else 'FAIL'}", file=sys.stderr)
    print(f"[self-test] noise only:     p_window = {nul['p_window_block_bootstrap']} / "
          f"{nul['p_window_phase_randomization']}  {'PASS' if ok2 else 'FAIL'}", file=sys.stderr)
    print(f"[self-test] 100σ spike (leak check): p_window = {leak['p_window_block_bootstrap']} / "
          f"{leak['p_window_phase_randomization']}  {'PASS' if ok3 else 'FAIL'}", file=sys.stderr)
    return ok1 and ok2 and ok3 and rank_self_test()


def rank_self_test(n_rep: int = 100, n_groups: int = 12) -> bool:
    """rank_window_test + fisher_across_groups on 12 synthetic 150-day groups
    (window ending on day 93, the project's typical fetch layout). Under H0 the
    combined p must be < 0.05 in about 5 % of replicates; the leave-window-out
    block bootstrap it replaces is shown per group for reference (printed, not
    asserted). A dip injected in half the groups must be detected. Single groups
    can't be tested for power: with ~4 fake blocks their p can't go below 0.2."""
    global N_SURROGATES
    rng = np.random.default_rng(SEED)
    dates = pd.date_range("2024-01-01", periods=150).strftime("%Y%m%d")
    window = set(pre_event_window(pd.Timestamp("2024-04-03"), PRE_WINDOW_DAYS))  # day 93 = 2024-04-03
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
