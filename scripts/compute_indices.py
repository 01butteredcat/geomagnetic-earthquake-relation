"""基本統計篩檢：建立候選的「局部異常指標」，
嘗試從整個觀測網共有的共模（Sq／磁暴驅動）變化中，
分離出只在震央附近的地磁訊號。

在該組兩個測站池中足夠的那一個上執行（見 common.py 的
GroupConfig）：XYZ 測站池提供 H 和 Z 場（向量資料）；
F 測站池（只有純量的測站）提供單一的 F 場——純量資料沒有 Z/D
的對應量，這是只量總磁場儀器的結構性限制，
不是 bug。一組可能跑其中一種、兩種，或（如果兩個測站池都
測站不足）完全不跑。

夜間窗口：原始 UTC 時 in {17,18,19} == 台灣當地時間 {01,02,03}，
已由 timezone_check.py 在 2024/G10 資料上確認是安靜的（假設在
各組都成立，因為測站／資料慣例相同）。我們刻意**不**把
它重新標到平移後的當地日曆日期——所有東西都以讀取樣本的 UTC
檔案日期為鍵，daily_features.csv
和 storm_days.csv 也是這樣。

方法（計畫第 3 節，推廣到任何場／測站池）：
  1. 從 1 分鐘序列算出每站每日的 night_mean/std。
  2. night_mean 相對於滑動 TRAILING_WINDOW_DAYS 天
     中位數的 MAD z-score（只用乾淨日——非磁暴、非中斷——建立）。
  3. near_index = 測站池近站的 dev_z 中位數；far_index =
     測站池遠站的 dev_z 中位數。
  4. 只在乾淨日上做 Theil-Sen 穩健迴歸 near_index ~ far_index；
     殘差 = local_anomaly_index（候選的局部訊號）。
  5. 標出 |local_anomaly_index| 超過自己乾淨日分布
     MAD 約 2.5 倍的候選日。
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import theilslopes

sys.path.insert(0, str(Path(__file__).parent))
from common import auto_outage_dates, load_group_config  # noqa: E402

NIGHT_HOURS_UTC = {17, 18, 19}  # == 當地 01:00-03:59
# 一個測站–夜的有效分鐘數少於這個（滿分 180）時，當天
# 不提供指標值。在此之前只要有一分鐘有效就夠，所以
# 整天 pct_missing 超過 OUTAGE_PCT_MISSING_THRESHOLD 的日子
# 仍可能因為幾乎全空的夜間窗口被標成候選日（例如 G18
# 2016-02-29：kmn 夜間只有 68/180 分鐘、hcn 是 0，那個很薄的 far_index
# 產生了一個候選日）。90（窗口的 50%）是判斷值：30-60
# 會留下那個 G18 的日子，90-150 去掉的候選日集合都一樣。
MIN_NIGHT_MINUTES = 90
# 7 天對 G10/2024 的情況太短：2024-03-21..27 的磁暴 +
# 2 天恢復期把到四月初為止的整個滑動窗口都吃掉，
# 讓最重要的震前日子（03-26..04-04）完全沒有
# 基準。21 天可以往回越過那次磁暴，延伸到一段乾淨期。
# 2026-09-28 從 21 改成 28：磁暴密集的期間（G6 2021-10、G11 2025-01、G19
# 2024-08）在 21 天下，那些錨點附近的日子仍只有 2-4 個乾淨日（< 5），
# 所以錨點日的指標是 NaN；28 是讓三者
# 都 >= 8 的最短窗口（30/35 也救不回更多）。
# 保留為一個全資料集共用的常數（不逐組設定），因為它是
# 「磁暴＋恢復期合理上可以持續多久」的性質，不是某一組資料的性質。
TRAILING_WINDOW_DAYS = 28
MIN_CLEAN_POINTS = 5
CANDIDATE_Z_THRESHOLD = 2.5


def night_features_for_station(cfg, station: str, channel: str) -> pd.DataFrame:
    path = cfg.interim_dir / f"minute_series_{station}.parquet"
    df = pd.read_parquet(path)
    night = df[df.index.hour.isin(NIGHT_HOURS_UTC)]
    date_key = night.index.strftime("%Y%m%d")
    g = night.groupby(date_key)
    out = pd.DataFrame(
        {
            f"night_mean_{channel}": g[channel].mean(),
            f"night_std_{channel}": g[channel].std(),
            "night_n": g[channel].count(),
        }
    )
    out.index.name = "date"
    out["station"] = station
    thin = out["night_n"] < MIN_NIGHT_MINUTES
    out.loc[thin, [f"night_mean_{channel}", f"night_std_{channel}"]] = np.nan
    return out.reset_index()


def mad_zscore(series: pd.Series, clean_mask: pd.Series) -> pd.Series:
    """滾動滑動窗口 MAD z-score，基準只用乾淨日建立。
    窗口是每一天之前的 TRAILING_WINDOW_DAYS 個日曆天，而不是
    前 TRAILING_WINDOW_DAYS 列：缺日檔的資料夾（G6_G7_G8
    缺 2021-12-30 和 2022-01-01）否則會往回延伸得更遠。"""
    z = pd.Series(index=series.index, dtype="float64")
    vals = series.values
    clean = clean_mask.values
    days = pd.to_datetime(series.index, format="%Y%m%d")
    starts = np.searchsorted(days, days - pd.Timedelta(days=TRAILING_WINDOW_DAYS))
    for i in range(len(series)):
        lo = int(starts[i])
        window_vals = vals[lo:i][clean[lo:i]] if i > lo else np.array([])
        window_vals = window_vals[~np.isnan(window_vals)]
        if len(window_vals) < MIN_CLEAN_POINTS:
            z.iloc[i] = np.nan
            continue
        med = np.median(window_vals)
        mad = np.median(np.abs(window_vals - med)) * 1.4826
        z.iloc[i] = (vals[i] - med) / mad if mad > 1e-9 else np.nan
    return z


def build_field_index(cfg, field: str, near_stations, far_stations, storm_dates, outage_by_station):
    all_stations = list(near_stations) + list(far_stations)
    frames = [night_features_for_station(cfg, s, field) for s in all_stations]
    night = pd.concat(frames, ignore_index=True)

    all_dates = sorted(night["date"].unique())
    pivot_mean = night.pivot(index="date", columns="station", values=f"night_mean_{field}").reindex(all_dates)

    dev_z = pd.DataFrame(index=all_dates)
    for station in all_stations:
        clean_mask = pd.Series(
            [d not in outage_by_station.get(station, set()) and d not in storm_dates for d in all_dates],
            index=all_dates,
        )
        dev_z[station] = mad_zscore(pivot_mean[station], clean_mask)

    near_index = dev_z[list(near_stations)].median(axis=1)
    far_index = dev_z[list(far_stations)].median(axis=1)

    # 一天算「乾淨」（可用於迴歸擬合／候選
    # 標記）的條件是：不是磁暴日，**且** near_index/far_index 真的
    # 算出了值——也就是測站層級的中位數至少有一個
    # 非中斷的貢獻者，而不是「測站池中每一個測站
    # 都必須各自沒有中斷」。比較嚴格的全有或全無版本對
    # G10 沒問題（那裡沒有任何近站／遠站曾長時間完全
    # 中斷），但如果某組的一個近站有例如 3 週的完全中斷
    # （G11 的 zbn 就在錨點事件前出現過），
    # 整段中斷期間就會一個乾淨日都沒有，即使
    # 3 站取中位數的 near_index 在那整段期間都能
    # 由另外兩站完整算出。
    clean_overall = pd.Series(
        [d not in storm_dates for d in all_dates], index=all_dates
    ) & near_index.notna() & far_index.notna()
    fit_mask = clean_overall  # 已經要求 near/far 非 NaN，見上面的 clean_overall
    if fit_mask.sum() >= 10:
        slope, intercept, _, _ = theilslopes(near_index[fit_mask], far_index[fit_mask])
    else:
        slope, intercept = 0.0, 0.0

    predicted = slope * far_index + intercept
    local_anomaly_index = near_index - predicted

    clean_resid = local_anomaly_index[clean_overall.reindex(local_anomaly_index.index, fill_value=False)]
    clean_resid = clean_resid.dropna()
    resid_mad = np.median(np.abs(clean_resid - np.median(clean_resid))) * 1.4826 if len(clean_resid) >= 5 else np.nan

    result = pd.DataFrame(
        {
            "date": all_dates,
            f"near_index_{field}": near_index.values,
            f"far_index_{field}": far_index.values,
            f"predicted_near_{field}": predicted.values,
            f"local_anomaly_index_{field}": local_anomaly_index.values,
        }
    )
    result["is_clean_day"] = clean_overall.values
    if resid_mad and resid_mad > 0:
        result[f"candidate_flag_{field}"] = (
            (result[f"local_anomaly_index_{field}"].abs() > CANDIDATE_Z_THRESHOLD * resid_mad)
            & result["is_clean_day"]
        )
    else:
        result[f"candidate_flag_{field}"] = False

    fit_info = {
        "field": field,
        "near_stations": list(near_stations),
        "far_stations": list(far_stations),
        "slope": float(slope),
        "intercept": float(intercept),
        "n_fit_days": int(fit_mask.sum()),
        "residual_mad": float(resid_mad) if resid_mad == resid_mad else None,
    }
    return result, fit_info


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--group", required=True)
    args = ap.parse_args()
    cfg = load_group_config(args.group)

    storm_days_df = pd.read_csv(cfg.interim_dir / "storm_days.csv", dtype={"date": str})
    storm_dates = set(storm_days_df["date"])
    daily = pd.read_csv(cfg.interim_dir / "daily_features.csv", dtype={"date": str})
    outage_by_station = auto_outage_dates(daily)

    merged = None
    fits: dict = {}
    candidates: dict[str, list[str]] = {}
    methods_run: list[str] = []

    if cfg.xyz_pool.sufficient:
        result_h, fit_h = build_field_index(cfg, "H", cfg.xyz_pool.near, cfg.xyz_pool.far, storm_dates, outage_by_station)
        result_z, fit_z = build_field_index(cfg, "Z", cfg.xyz_pool.near, cfg.xyz_pool.far, storm_dates, outage_by_station)
        merged = result_h.merge(result_z.drop(columns=["is_clean_day"]), on="date", how="outer")
        fits["H"], fits["Z"] = fit_h, fit_z
        candidates["H"] = merged[merged["candidate_flag_H"] == True]["date"].tolist()  # noqa: E712
        candidates["Z"] = merged[merged["candidate_flag_Z"] == True]["date"].tolist()  # noqa: E712
        methods_run += ["H", "Z"]
    else:
        print(f"[{args.group}] XYZ pool insufficient ({len(cfg.xyz_pool.all_stations)} stations) -- skipping H/Z screening", file=sys.stderr)

    if cfg.f_pool.sufficient:
        result_f, fit_f = build_field_index(cfg, "F", cfg.f_pool.near, cfg.f_pool.far, storm_dates, outage_by_station)
        fits["F"] = fit_f
        candidates["F"] = result_f[result_f["candidate_flag_F"] == True]["date"].tolist()  # noqa: E712
        methods_run.append("F")
        if merged is None:
            merged = result_f
        else:
            merged = merged.merge(
                result_f.rename(columns={"is_clean_day": "is_clean_day_F"}), on="date", how="outer"
            )
    else:
        print(f"[{args.group}] F pool insufficient ({len(cfg.f_pool.all_stations)} stations) -- skipping F screening", file=sys.stderr)

    if merged is None:
        print(f"[{args.group}] neither station pool sufficient -- no basic screening possible for this group", file=sys.stderr)
        merged = pd.DataFrame(columns=["date"])
    merged = merged.sort_values("date")
    merged.to_csv(cfg.interim_dir / "local_anomaly_index.csv", index=False)

    both_hz = sorted(set(candidates.get("H", [])) & set(candidates.get("Z", [])))

    out = {
        "group": args.group,
        "methods_run": methods_run,
        "fit_H": fits.get("H"),
        "fit_Z": fits.get("Z"),
        "fit_F": fits.get("F"),
        "candidate_dates_H": candidates.get("H", []),
        "candidate_dates_Z": candidates.get("Z", []),
        "candidate_dates_F": candidates.get("F", []),
        "candidate_dates_both_H_and_Z": both_hz,
        "n_storm_or_recovery_dates_excluded": len(storm_dates),
    }
    (cfg.interim_dir / "candidate_windows.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
