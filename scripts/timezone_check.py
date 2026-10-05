"""以實際資料判斷 .sec 檔原始 TIME 欄位是
UTC 還是台灣當地時間（UTC+8）。IAGA-2002 檔頭沒有時區
欄位，所以這是用資料本身推斷的，使用該組最近可用測站的兩個獨立
日變化訊號：

  訊號 A：測站主要強度訊號的 Sq（太陽靜日）日變化曲線
  （向量站用 H = sqrt(X^2+Y^2)，純量站
  直接用 F）。Sq 電流系統會在當地太陽正午附近產生一個寬的極值。
  台灣幾乎正好位在東經 120 度的 UTC+8
  標準經線上，所以當地太陽正午基本上就是台灣當地時間的
  時鐘正午——和 UTC 時鐘正午差 8 小時，相對於 Sq 特徵的寬度
  這個差距很大，因此可以明確
  區分。

  訊號 B：高頻雜訊功率的日變化（同一主要通道
  1 秒一階差分的變異數）。人為／工業
  活動（電網、鐵路、交通）大約在**真正的**當地時鐘時間
  07:00-19:00 會明顯升高，01:00-05:00 則安靜——這是
  更銳利、獨立的交叉檢查。

兩個訊號都在所選測站最安靜的一批可用日子上計算
（每日標準差最低，依 data/interim/<group>/daily_features.csv），
排除缺值過多的日子
（common.auto_outage_dates）。
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from common import auto_outage_dates, load_group_config, resolve_day_ref  # noqa: E402
from parser import parse_day_file  # noqa: E402

N_QUIET_DAYS = 18


def _pick_station(cfg) -> str:
    """優先用最近的向量站（H = sqrt(X^2+Y^2) 是比較乾淨的
    Sq／雜訊訊號）；沒有可用向量測站池的組別（例如 G1）
    則退而用最近的純量站（F）。"""
    if cfg.xyz_pool.near:
        return cfg.xyz_pool.near[0]
    if cfg.f_pool.near:
        return cfg.f_pool.near[0]
    raise RuntimeError(f"{cfg.group_id}：完全沒有可用來做時區檢查的測站")


def _signals(df: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
    """回傳 (sq_signal, noise_base_signal)，取決於這個
    檔案檔頭實際回報了哪些欄位（見 parser.parse_day_file）。"""
    if "F" in df.columns:
        return df["F"], df["F"]
    h = np.sqrt(df["X"] ** 2 + df["Y"] ** 2)
    return h, df["X"]


def pick_quiet_days(cfg, station: str) -> list[str]:
    std_col = "F_std" if station in cfg.f_pool.all_stations else "H_std"
    spike_col = "n_spikes"
    daily = pd.read_csv(cfg.interim_dir / "daily_features.csv", dtype={"date": str})
    sub = daily[(daily.station == station) & (daily[spike_col] == 0) & (daily.pct_missing == 0)]
    outage_dates = auto_outage_dates(daily).get(station, set())
    sub = sub[~sub.date.isin(outage_dates)]
    # 分散在整段期間：從每個連續約 1 週的區段
    # 各取標準差最低的一天，總共最多 N_QUIET_DAYS 天
    sub = sub.sort_values("date").reset_index(drop=True)
    if len(sub) == 0:
        return []
    sub["bucket"] = np.arange(len(sub)) // max(1, len(sub) // N_QUIET_DAYS)
    picks = sub.loc[sub.groupby("bucket")[std_col].idxmin(), "date"].tolist()
    return sorted(picks)[:N_QUIET_DAYS]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--group", required=True)
    args = ap.parse_args()
    cfg = load_group_config(args.group)
    station = _pick_station(cfg)

    quiet_days = pick_quiet_days(cfg, station)
    print(f"[{args.group}] station={station} 使用的平靜日（{len(quiet_days)} 天）：{quiet_days}", file=sys.stderr)

    sq_by_hour = {h: [] for h in range(24)}
    noise_by_hour = {h: [] for h in range(24)}

    for date_str in quiet_days:
        ref = resolve_day_ref(cfg.gdms_dir, station, date_str)
        df = parse_day_file(ref, station)
        sq_signal, noise_base = _signals(df)
        noise = noise_base.diff().abs()  # 1 秒一階差分的大小
        hour = df.index.hour
        for hr in range(24):
            mask = hour == hr
            sq_by_hour[hr].append(float(sq_signal[mask].mean()))
            noise_by_hour[hr].append(float((noise[mask] ** 2).mean()))  # 類似變異數的代理量

    noise_curve = {hr: float(np.mean(vals)) for hr, vals in noise_by_hour.items()}

    # 平均之前先扣掉每天自己的平均值，以免日與日之間的基線
    # 漂移蓋過日變化的形狀
    sq_by_hour2 = {h: [] for h in range(24)}
    for date_str in quiet_days:
        ref = resolve_day_ref(cfg.gdms_dir, station, date_str)
        df = parse_day_file(ref, station)
        sq_signal, _ = _signals(df)
        dev = sq_signal - sq_signal.mean()
        hour = df.index.hour
        for hr in range(24):
            sq_by_hour2[hr].append(float(dev[hour == hr].mean()))
    sq_curve_dev = {hr: float(np.mean(vals)) for hr, vals in sq_by_hour2.items()}

    sq_extreme_hour = max(sq_curve_dev, key=lambda hr: abs(sq_curve_dev[hr])) if quiet_days else None
    noise_peak_hours = sorted(noise_curve, key=lambda hr: -noise_curve[hr])[:6] if quiet_days else []
    noise_quiet_hours = sorted(noise_curve, key=lambda hr: noise_curve[hr])[:6] if quiet_days else []

    # 判斷規則：如果雜訊在原始時鐘約 7-19 時高、約 1-5 時低，
    # 原始時鐘 == 台灣當地時間（UTC+8）。如果這個型態反而
    # 提早 8 小時（高峰約 23-11、安靜約 17-21），原始時鐘 == UTC。
    workday_hours = set(range(7, 20))
    night_hours = set(range(1, 6))
    workday_score = sum(1 for hr in noise_peak_hours if hr in workday_hours)
    night_score = sum(1 for hr in noise_quiet_hours if hr in night_hours)

    if not quiet_days:
        conclusion, confidence = "UNCERTAIN", "low"
    elif workday_score >= 4 and night_score >= 4:
        conclusion = "LOCAL (UTC+8)"
        confidence = "high"
    else:
        # 檢查平移 8 小時的假設（UTC 標記）
        workday_hours_utc = {(hr - 8) % 24 for hr in workday_hours}
        night_hours_utc = {(hr - 8) % 24 for hr in night_hours}
        workday_score_utc = sum(1 for hr in noise_peak_hours if hr in workday_hours_utc)
        night_score_utc = sum(1 for hr in noise_quiet_hours if hr in night_hours_utc)
        if workday_score_utc >= 4 and night_score_utc >= 4:
            conclusion = "UTC"
            confidence = "high"
        else:
            conclusion = "UNCERTAIN"
            confidence = "low"

    result = {
        "group": args.group,
        "station": station,
        "quiet_days_used": quiet_days,
        "sq_curve_H_deviation_by_raw_hour": sq_curve_dev,
        "noise_power_by_raw_hour": noise_curve,
        "sq_extreme_raw_hour": sq_extreme_hour,
        "noise_peak_raw_hours": noise_peak_hours,
        "noise_quiet_raw_hours": noise_quiet_hours,
        "workday_score_if_local": workday_score,
        "night_score_if_local": night_score,
        "conclusion_timezone": conclusion,
        "confidence": confidence,
        "note": (
            "conclusion is 'LOCAL (UTC+8)' if noise peaks in raw-clock 07-19 and "
            "is quiet in raw-clock 01-05 (matching real Taiwan work-hours); "
            "'UTC' if that pattern is shifted 8h earlier in the raw clock."
        ),
    }

    out_path = cfg.interim_dir / "timezone_check.json"
    out_path.write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))
    print(f"\nwrote {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
