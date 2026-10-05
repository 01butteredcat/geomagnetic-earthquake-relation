"""疊加時間分析（SEA）：教授「把一次地震變成很多次」的
建議。不只看 G10 單一個 2024-04-03 主震，而是拿
8 個有 ULF 近站／遠站資料的組別中，觀測網附近每一起獨立的 M>=5.5 地震
（events.py 中人工挑選的 13 起 M>=6.0 錨點／子事件，
加上 `fetch_earthquake_catalog.py` 另外找到的），
把每組的 pc3/pc4 近站–遠站極化 z-score 序列對齊到「相對於
該地震發震時間的天數」，再對所有事件疊加（平均）。
如果 G10 主震前 4 天（2024-03-30）看到的那種低谷
是真實、可重現的前兆特徵，而不是一次性的，它應該能
撐過平均，並在用同樣方法、以隨機且與地震無關的參考日期
建立的虛無帶中凸顯出來。

逐組的 z-score 化（透過 stat_utils.mad_zscore，和這整套驗證
使用相同公式）在疊加**之前**進行，因為原始
diff_zh 大小在雜訊底不同的組別／測站之間
沒有可比較的尺度——直接跨組疊加 nT 等級的原始差值
沒有意義。

用法：
  superposed_epoch_analysis.py --catalog data/external/extended_catalog_m5.5.csv --label m5.5 --min-mag 5.5
  superposed_epoch_analysis.py --catalog data/external/extended_catalog_m5.0.csv --label m5.0 --min-mag 5.0
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
from events import folder_events  # noqa: E402
from stat_utils import mad_zscore  # noqa: E402

ULF_GROUPS = ("G4", "G5", "G6", "G7", "G8", "G9", "G10", "G11", "G12", "G13", "G19", "G20", "G23", "G24")
BANDS = ("pc3", "pc4")
WINDOW_BEFORE_DAYS = 30
WINDOW_AFTER_DAYS = 10
N_BOOTSTRAP = 2000       # 真實疊加的信賴區間（對事件做有放回重抽）
N_NULL = 1000            # 虛無實現次數（隨機重抽時期日期）
NULL_EXCLUSION_BUFFER_DAYS = 30  # 假時期要離任何真實事件至少這麼遠
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
        # 把 z 展開回完整（含 NaN）長度，對齊 df 的列
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


def _eligible_null_days(band: str, events: list[dict], group_series: dict) -> dict[str, list[pd.Timestamp]]:
    """每組的候選假時期：完整窗口落在資料內**且**距離
    每個真實日期至少 NULL_EXCLUSION_BUFFER_DAYS 的每一天——真實日期包括這個級距的事件，以及
    該組原始資料夾中所有已登錄的事件（兄弟組看到的是同樣的日子）。從
    這個集合抽樣，而不是隨機重試日期、最後退回最後一次嘗試，可以保證沒有
    虛無時期落在真實地震附近。密集的組別（例如 G11 的 32 起事件）可能
    完全沒有符合條件的日子。"""
    real_dates_by_group: dict[str, list[pd.Timestamp]] = {}
    for ev in events:
        real_dates_by_group.setdefault(ev["group"], []).append(ev["date"])
    for g in list(real_dates_by_group):
        real_dates_by_group[g] += [pd.Timestamp(e.time_utc.split(" ")[0]) for e in folder_events(g)]
    eligible: dict[str, list[pd.Timestamp]] = {}
    for g, real in real_dates_by_group.items():
        series = group_series.get(g)
        if series is None or band not in series:
            continue
        lo, hi = valid_date_range(series)
        lo_bound = lo + pd.Timedelta(days=WINDOW_BEFORE_DAYS)
        hi_bound = hi - pd.Timedelta(days=WINDOW_AFTER_DAYS)
        eligible[g] = [d for d in pd.date_range(lo_bound, hi_bound, freq="D")
                       if all(abs((d - rd).days) >= NULL_EXCLUSION_BUFFER_DAYS for rd in real)]
    return eligible


def _event_windows(band: str, events: list[dict], group_series: dict) -> tuple[list, list[dict]]:
    matrix, used = [], []
    for ev in events:
        series = group_series.get(ev["group"])
        if series is None or band not in series:
            continue
        w = extract_window(series, band, ev["date"])
        if np.all(np.isnan(w)):
            continue
        matrix.append(w)
        used.append({"group": ev["group"], "date": ev["date"].strftime("%Y-%m-%d"),
                     "mag": ev["mag"], "source": ev["source"]})
    return matrix, used


def merge_same_day(events: list[dict]) -> list[dict]:
    """每個（組, 日）一筆疊加項：同一組同一天的事件（例如 G8 的三起
    2022-03-22 M6）共用完全相同的每日窗口，所以逐一疊加會讓那個
    窗口被算 2-3 次。保留最大事件那一列，和 backtest_rule.py 的逐日集合一樣。"""
    best: dict[tuple, dict] = {}
    for ev in events:
        key = (ev["group"], ev["date"])
        if key not in best or ev["mag"] > best[key]["mag"]:
            best[key] = ev
    return list(best.values())


def _r4(a) -> list:
    return [None if np.isnan(v) else round(float(v), 4) for v in a]


def run_band(band: str, events: list[dict], group_series: dict, rng: np.random.Generator) -> dict:
    """和虛無帶比較的疊加只用有符合條件虛無日的組別，
    讓真實疊加和虛無帶由相同的組別建立——否則一個密集組
    缺席虛無分布（G11 的疊加比其他組高約 2 個 z）會讓真實疊加
    在每個延遲都跑出虛無帶。這個級距所有事件的疊加另外
    保留為 `*_all_events`，只作描述用。"""
    lags = np.arange(-WINDOW_BEFORE_DAYS, WINDOW_AFTER_DAYS + 1)
    all_matrix, all_used = _event_windows(band, events, group_series)
    if not all_matrix:
        return {"band": band, "error": "no events with usable data"}
    M_all = np.array(all_matrix)

    eligible_days = _eligible_null_days(band, events, group_series)
    groups_without_eligible = sorted(g for g, days in eligible_days.items() if not days)
    comparable = [ev for ev in events if eligible_days.get(ev["group"])]
    matrix, used_events = _event_windows(band, comparable, group_series)
    if not matrix:
        return {"band": band, "error": "no events in groups with eligible null days",
                "null_groups_without_eligible_days": groups_without_eligible}
    M = np.array(matrix)  # n_events x n_lags 矩陣

    real_stack_mean = np.nanmean(M, axis=0)
    real_stack_median = np.nanmedian(M, axis=0)
    n_contributing = np.sum(~np.isnan(M), axis=0)

    # 真實疊加的 bootstrap 信賴區間（對事件做有放回重抽）
    n_events = M.shape[0]
    boot = np.empty((N_BOOTSTRAP, M.shape[1]))
    for b in range(N_BOOTSTRAP):
        idx = rng.integers(0, n_events, size=n_events)
        boot[b] = np.nanmean(M[idx], axis=0)
    ci_lo = np.nanpercentile(boot, 5, axis=0)
    ci_hi = np.nanpercentile(boot, 95, axis=0)

    # 虛無帶：用隨機、與地震無關的時期日期重複整個疊加程序
    # （每個可比較的事件抽一個假時期，從該事件所屬的組抽）
    null_stacks = np.empty((N_NULL, M.shape[1]))
    for r in range(N_NULL):
        null_matrix = []
        for ev in comparable:
            days = eligible_days[ev["group"]]
            w = extract_window(group_series[ev["group"]], band, days[int(rng.integers(0, len(days)))])
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
        "stack_mean": _r4(real_stack_mean),
        "stack_median": _r4(real_stack_median),
        "n_contributing_per_lag": n_contributing.tolist(),
        "bootstrap_ci90_lo": [round(float(v), 4) for v in ci_lo],
        "bootstrap_ci90_hi": [round(float(v), 4) for v in ci_hi],
        "null_band_p5": [round(float(v), 4) for v in null_lo],
        "null_band_p50": [round(float(v), 4) for v in null_median],
        "null_band_p95": [round(float(v), 4) for v in null_hi],
        "null_exclusion_buffer_days": NULL_EXCLUSION_BUFFER_DAYS,
        "null_groups_without_eligible_days": groups_without_eligible,
        "n_events_all": len(all_used),
        "events_all": all_used,
        "stack_mean_all_events": _r4(np.nanmean(M_all, axis=0)),
        "n_contributing_per_lag_all_events": np.sum(~np.isnan(M_all), axis=0).tolist(),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--catalog", type=Path, required=True)
    ap.add_argument("--label", required=True, help="例如 m5.5 或 m5.0，用在輸出檔名中")
    ap.add_argument("--min-mag", type=float, required=True,
                    help="這個級距的規模門檻；必須和 --catalog 檔案的一致")
    args = ap.parse_args()

    group_series = {}
    for g in ULF_GROUPS:
        s = load_group_series(g)
        if s is not None:
            group_series[g] = s

    raw = load_extended_events(args.catalog, ULF_GROUPS, args.min_mag)
    events = merge_same_day(raw)
    print(f"{len(raw)} 起候選事件 -> {len(events)} 個組–日，分布在 {len(group_series)} 組",
          file=sys.stderr)

    out_dir = PROJECT_DIR / "data" / "interim" / "superposed_epoch"
    out_dir.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(SEED)
    for band in BANDS:
        result = run_band(band, events, group_series, rng)
        result["min_mag"] = args.min_mag
        out_path = out_dir / f"{band}_stack_{args.label}.json"
        out_path.write_text(json.dumps(result, indent=2))
        if "error" in result:
            print(f"[{band}] {result['error']}", file=sys.stderr)
        else:
            print(f"[{band}] 已寫入 {out_path}（使用 {result['n_events_used']} 起事件）", file=sys.stderr)


if __name__ == "__main__":
    main()
