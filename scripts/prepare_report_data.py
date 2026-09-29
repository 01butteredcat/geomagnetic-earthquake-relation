"""Gather all small interim artifacts into one JSON blob for the final
self-contained HTML report (build_artifact.py embeds this verbatim).

Usage: prepare_report_data.py --group G10

Note: `report_template.html`'s narrative text (twu instrument-failure
commentary, the 3/30 candidate-day discussion, etc.) is hand-written for the
G10/2024 single event specifically and is NOT generalized to other groups by
this refactor -- this script can technically run for any group (paths/near-
far stations are all dynamic), but the resulting report_data.json's
qualitative narrative fields (spike/glitch commentary etc.) are only
meaningful in the context report_template.html was written for. Producing
polished reports for the other 12 groups is out of scope; this script's role
in the 13-group batch flow is limited to what run_pipeline.sh actually uses
(everything up to and including verify_pipeline.py) -- this script is only
invoked when generating G10's original detailed report.
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

    # spike/glitch narrative data (twu instrument failure, ttn escalation --
    # for G10; generic for any other group)
    spikes = daily[daily.n_spikes > 0][["station", "date", "n_spikes", "max_abs_jump", "pct_missing"]]
    spike_by_station = {}
    for station, g in spikes.groupby("station"):
        spike_by_station[station] = g.drop(columns="station").to_dict(orient="records")

    # whole-day outages (e.g. G10's twu from 2024-04-24 on) -- the report's
    # quality narrative quotes these instead of hard-coding dates
    full_outage = daily[daily.pct_missing >= 0.99]
    full_outage_by_station = {s: sorted(g["date"].tolist()) for s, g in full_outage.groupby("station")}

    coverage = daily.groupby("station").agg(
        n_days=("date", "count"),
        avg_pct_missing=("pct_missing", "mean"),
        max_pct_missing=("pct_missing", "max"),
    ).reset_index().to_dict(orient="records")

    local_idx_records = local_idx.replace({float("nan"): None}).to_dict(orient="records")
    ulf_nf_records = ulf_nf.replace({float("nan"): None}).to_dict(orient="records") if len(ulf_nf) else []

    # G10 keeps its original hand-curated per-incident outage notes (richer
    # than the generic auto-detected summary) for backward compatibility
    # with report_template.html's existing table; other groups get a
    # generic auto-detected summary (station/date/pct_missing) instead.
    if args.group == "G10":
        # the hand-curated list predates the folder's extension to 2024-06-01;
        # stretch its twu "fully missing" row to the end of the actual whole-day outage
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
