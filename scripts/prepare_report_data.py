"""把所有小型中間產物收集成一個 JSON，給最終的
單檔 HTML 報告用（build_artifact.py 會原封不動地嵌入）。

用法：prepare_report_data.py --group G10

注意：`report_template.html` 的敘述文字（twu 儀器故障
的評論、3/30 候選日的討論等）是專門為
G10/2024 單一事件手寫的，這次重構**沒有**把它推廣到其他組——
這支腳本技術上可以對任何組執行（路徑／近站–
遠站都是動態的），但產生的 report_data.json 中
質性的敘述欄位（突波／故障評論等）只在
report_template.html 原本撰寫的脈絡下才有意義。替
其他 12 組產出精修報告不在範圍內；這支腳本在
13 組批次流程中的角色只限於 run_pipeline.sh 實際用到的部分
（到 verify_pipeline.py 為止的所有步驟）——只有在
產生 G10 原本那份詳細報告時才會呼叫這支腳本。
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from common import G10_KNOWN_OUTAGE_WINDOWS, auto_outage_dates, load_group_config  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--group", required=True)
    args = ap.parse_args()
    cfg = load_group_config(args.group)

    daily = pd.read_csv(cfg.interim_dir / "daily_features.csv", dtype={"date": str})
    local_idx = pd.read_csv(cfg.interim_dir / "local_anomaly_index.csv", dtype={"date": str})
    ulf_nf_path = cfg.interim_dir / "ulf_near_far_index.csv"
    ulf_nf = pd.read_csv(ulf_nf_path, dtype={"date": str}) if ulf_nf_path.exists() else pd.DataFrame()
    storm_days = pd.read_csv(cfg.interim_dir / "storm_days.csv", dtype={"date": str})
    tz = json.loads((cfg.interim_dir / "timezone_check.json").read_text())
    candidates = json.loads((cfg.interim_dir / "candidate_windows.json").read_text())
    verification = json.loads((cfg.interim_dir / "verification_report.json").read_text())
    storm_summary = json.loads((cfg.interim_dir / "storm_days_summary.json").read_text())

    pool = cfg.xyz_pool if cfg.xyz_pool.sufficient else cfg.f_pool
    near_stations, far_stations = pool.near, pool.far

    stations = [
        {"code": code, "name": meta["name"], "lat": meta["lat"], "lon": meta["lon"],
         "distance_km": meta["distance_km"],
         "role": "near" if code in near_stations else ("far" if code in far_stations else "other")}
        for code, meta in sorted(cfg.stations.items(), key=lambda kv: kv[1]["distance_km"])
    ]

    # 突波／故障的敘述資料（twu 儀器故障、ttn 惡化——
    # 針對 G10；其他組則是通用版）
    spikes = daily[daily.n_spikes > 0][["station", "date", "n_spikes", "max_abs_jump", "pct_missing"]]
    spike_by_station = {}
    for station, g in spikes.groupby("station"):
        spike_by_station[station] = g.drop(columns="station").to_dict(orient="records")

    # 整天中斷（例如 G10 的 twu 從 2024-04-24 起）——報告的
    # 品質敘述會引用這些，而不是寫死日期
    full_outage = daily[daily.pct_missing >= 0.99]
    full_outage_by_station = {s: sorted(g["date"].tolist()) for s, g in full_outage.groupby("station")}

    coverage = daily.groupby("station").agg(
        n_days=("date", "count"),
        avg_pct_missing=("pct_missing", "mean"),
        max_pct_missing=("pct_missing", "max"),
    ).reset_index().to_dict(orient="records")

    local_idx_records = local_idx.replace({float("nan"): None}).to_dict(orient="records")
    ulf_nf_records = ulf_nf.replace({float("nan"): None}).to_dict(orient="records") if len(ulf_nf) else []

    # G10 保留它原本人工整理的逐次中斷說明（比
    # 通用的自動偵測摘要更詳細），以相容於
    # report_template.html 既有的表格；其他組則改用
    # 通用的自動偵測摘要（station/date/pct_missing）。
    if args.group == "G10":
        # 人工整理的清單早於資料夾延長到 2024-06-01；
        # 把它的 twu「完全缺漏」那一列延長到實際整天中斷的結束日
        known_outage_windows = [dict(w) for w in G10_KNOWN_OUTAGE_WINDOWS]
        twu_out = full_outage_by_station.get("twu", [])
        for w in known_outage_windows:
            if w["station"] == "twu" and w["note"] == "fully missing" and twu_out:
                last = pd.to_datetime(twu_out[-1], format="%Y%m%d").strftime("%Y-%m-%d")
                w["end"] = f"{last} 23:59:59"
    else:
        auto = auto_outage_dates(daily)
        known_outage_windows = [
            {"station": s, "start": d, "end": d, "note": "auto-detected (pct_missing > threshold)"}
            for s, dates in sorted(auto.items()) for d in sorted(dates)
        ]

    spec_near = spec_far = None
    if near_stations:
        p = cfg.interim_dir / f"ulf_spectrogram_{near_stations[0]}_eq_window.json"
        if p.exists():
            spec_near = json.loads(p.read_text())
    if far_stations:
        p = cfg.interim_dir / f"ulf_spectrogram_{far_stations[-1]}_eq_window.json"
        if p.exists():
            spec_far = json.loads(p.read_text())

    a = cfg.anchor_event
    report_data = {
        "meta": {
            "group": args.group,
            "eq_time_utc": a.time_utc,
            "eq_time_local": a.time_local,
            "eq_lat": a.lat,
            "eq_lon": a.lon,
            "eq_depth_km": a.depth_km,
            "eq_magnitude": a.magnitude,
            "eq_magnitude_type": a.magnitude_type,
            "eq_coord_confidence": a.coord_confidence,
            "data_start": pd.to_datetime(daily["date"].min(), format="%Y%m%d").strftime("%Y-%m-%d"),
            "data_end": pd.to_datetime(daily["date"].max(), format="%Y%m%d").strftime("%Y-%m-%d"),
            "data_timezone": tz["conclusion_timezone"],
            "near_stations": list(near_stations),
            "far_stations": list(far_stations),
        },
        "stations": stations,
        "coverage": coverage,
        "known_outage_windows": known_outage_windows,
        "spike_glitches_by_station": spike_by_station,
        "full_outage_days_by_station": full_outage_by_station,
        "timezone_check": tz,
        "storm_summary": storm_summary,
        "storm_dates": sorted(storm_days["date"].tolist()),
        "storm_onset_dates": sorted(storm_days[storm_days.is_storm_onset]["date"].tolist()),
        "local_anomaly_index": local_idx_records,
        "ulf_near_far_index": ulf_nf_records,
        "candidate_windows": candidates,
        "verification": verification,
        "spectrogram_near": spec_near,
        "spectrogram_far": spec_far,
    }

    out_path = cfg.interim_dir / "report_data.json"
    out_path.write_text(json.dumps(report_data))
    print(f"wrote {out_path} ({out_path.stat().st_size / 1024:.0f} KB)", file=sys.stderr)


if __name__ == "__main__":
    main()
