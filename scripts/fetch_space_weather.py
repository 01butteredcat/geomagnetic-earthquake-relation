"""Fetch historical Kp (GFZ Potsdam) and Dst (Kyoto WDC) space weather
indices for a group's actual date range, and derive a storm-day flag list
used to keep geomagnetic-storm-driven variability from being misread as a
tectonic precursor signal.

Usage: fetch_space_weather.py --group G10

Data sources (plain-text, no auth, fetched once and cached per group):
  Kp:  https://kp.gfz.de/kpdata?startdate=...&enddate=...&format=kp2
       (kp.gfz-potsdam.de 301-redirects here; -L follows it). Kp has been
       recorded since 1932, so any group's date range (2017-2026) is covered.
  Dst: https://wdc.kugi.kyoto-u.ac.jp/dst_<final|provisional|realtime>/<YYYYMM>/dst<YYMM>.for.request
       (requires a browser-like User-Agent or Kyoto returns 403). Which of
       final/provisional/realtime is available depends on how long ago the
       month was -- tried in that order per month since finalized data is
       more likely for older groups and realtime more likely for the most
       recent one (G13, 2026).

If either fetch fails, falls back to an internal proxy: days where the
group's far-reference stations show simultaneous large deviations are
flagged as "globally disturbed", labeled lower confidence.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from common import list_day_refs, load_group_config  # noqa: E402

KP_URL_TMPL = "https://kp.gfz.de/kpdata?startdate={start}&enddate={end}&format=kp2"
DST_URL_TMPL = "https://wdc.kugi.kyoto-u.ac.jp/dst_{kind}/{yyyymm}/dst{yymm}.for.request"
DST_KINDS = ("final", "provisional", "realtime")
UA = "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36"

KP_STORM_THRESHOLD = 5.0  # G1+ geomagnetic storm
DST_STORM_THRESHOLD = -30.0  # nT trough


def _curl(url: str, extra_args: list[str] | None = None) -> str | None:
    args = ["curl", "-sL", "--fail", "-A", UA]
    if extra_args:
        args += extra_args
    args.append(url)
    try:
        out = subprocess.run(args, capture_output=True, timeout=30, check=True)
        return out.stdout.decode("utf-8", errors="replace")
    except Exception as exc:  # noqa: BLE001
        print(f"fetch failed for {url}: {exc}", file=sys.stderr)
        return None


def _group_date_range(cfg) -> tuple[str, str]:
    dates = sorted(r.date_str for r in list_day_refs(cfg.gdms_dir))
    d0, d1 = dates[0], dates[-1]
    return f"{d0[:4]}-{d0[4:6]}-{d0[6:8]}", f"{d1[:4]}-{d1[4:6]}-{d1[6:8]}"


def _month_range(start: str, end: str) -> list[str]:
    months = pd.period_range(start=start, end=end, freq="M")
    return [str(m).replace("-", "") for m in months]


def fetch_kp(start: str, end: str) -> pd.DataFrame | None:
    text = _curl(KP_URL_TMPL.format(start=start, end=end))
    if not text:
        return None
    rows = []
    for line in text.strip().splitlines():
        parts = line.split()
        if len(parts) < 8:
            continue
        year, month, day = int(parts[0]), int(parts[1]), int(parts[2])
        start_hr = float(parts[3])
        kp = float(parts[7])
        rows.append({"date": f"{year:04d}{month:02d}{day:02d}", "start_hour_utc": start_hr, "kp": kp})
    return pd.DataFrame(rows)


def fetch_dst(start: str, end: str) -> pd.DataFrame | None:
    rows = []
    for yyyymm in _month_range(start, end):
        yymm = yyyymm[2:]
        text = None
        for kind in DST_KINDS:
            url = DST_URL_TMPL.format(kind=kind, yyyymm=yyyymm, yymm=yymm)
            text = _curl(url)
            if text:
                break
        if not text:
            return None
        for line in text.strip().splitlines():
            if not line.startswith("DST"):
                continue
            # e.g. "DST2404*01PPX120   0  -4  -5  -9 ... -7"
            head, rest = line[:20], line[20:]
            day = int(head[8:10])
            vals = rest.split()
            if len(vals) < 25:
                continue
            hourly = [float(v) for v in vals[:24]]
            daily_mean = float(vals[24])
            date_str = f"{yyyymm}{day:02d}"
            for hr, v in enumerate(hourly, start=1):
                rows.append({"date": date_str, "hour_utc": hr % 24, "dst": v})
            rows.append({"date": date_str, "hour_utc": -1, "dst": daily_mean})  # -1 = daily mean marker
    return pd.DataFrame(rows)


def internal_proxy_storm_days(cfg) -> pd.DataFrame:
    """Fallback: flag days where far-reference stations show simultaneous
    large H/F deviation, as an internally-inferred (lower-confidence) storm
    proxy. Uses the group's XYZ far pool if available, else its F far pool."""
    far = cfg.xyz_pool.far or cfg.f_pool.far
    range_col = "H_range" if cfg.xyz_pool.far else "F_range"
    daily = pd.read_csv(cfg.interim_dir / "daily_features.csv", dtype={"date": str})
    sub = daily[daily.station.isin(far)].copy()
    sub["z"] = sub.groupby("station")[range_col].transform(
        lambda s: (s - s.median()) / (1.4826 * (s - s.median()).abs().median() + 1e-9)
    )
    agg = sub.groupby("date")["z"].median().reset_index()
    agg["storm_flag"] = agg["z"] > 3.0
    agg["source"] = f"internal_proxy (far-station {range_col} MAD z>3, LOW CONFIDENCE)"
    return agg[agg.storm_flag][["date", "source"]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--group", required=True)
    args = ap.parse_args()
    cfg = load_group_config(args.group)
    start, end = _group_date_range(cfg)
    print(f"[{args.group}] fetching space weather for {start}..{end}", file=sys.stderr)

    kp_df = fetch_kp(start, end)
    dst_df = fetch_dst(start, end)

    source_note = []
    storm_days = set()

    if kp_df is not None and len(kp_df):
        kp_df.to_csv(cfg.external_dir / "kp.csv", index=False)
        kp_daily_max = kp_df.groupby("date")["kp"].max()
        kp_storm_days = set(kp_daily_max[kp_daily_max >= KP_STORM_THRESHOLD].index)
        storm_days |= kp_storm_days
        source_note.append(f"Kp (GFZ): {len(kp_df)} 3h records, {len(kp_storm_days)} storm days (max Kp>={KP_STORM_THRESHOLD})")
    else:
        source_note.append("Kp (GFZ): FETCH FAILED")

    if dst_df is not None and len(dst_df):
        dst_df.to_csv(cfg.external_dir / "dst.csv", index=False)
        dst_hourly = dst_df[dst_df.hour_utc >= 0]
        dst_daily_min = dst_hourly.groupby("date")["dst"].min()
        dst_storm_days = set(dst_daily_min[dst_daily_min <= DST_STORM_THRESHOLD].index)
        storm_days |= dst_storm_days
        source_note.append(f"Dst (Kyoto WDC): {len(dst_df)} records, {len(dst_storm_days)} storm days (min Dst<={DST_STORM_THRESHOLD}nT)")
    else:
        source_note.append("Dst (Kyoto WDC): FETCH FAILED")

    confidence = "high (official Kp/Dst indices)"
    if kp_df is None and dst_df is None:
        proxy = internal_proxy_storm_days(cfg)
        storm_days |= set(proxy["date"])
        source_note.append(f"FALLBACK: internal far-station proxy, {len(proxy)} storm days")
        confidence = "low (internal proxy only, official sources unavailable)"

    # add 1-2 recovery days after each storm day (ring-current decay)
    all_dates = sorted(storm_days)
    extended = set(storm_days)
    for d in all_dates:
        dt = pd.to_datetime(d, format="%Y%m%d")
        for delta in (1, 2):
            extended.add((dt + pd.Timedelta(days=delta)).strftime("%Y%m%d"))

    storm_df = pd.DataFrame(sorted(extended), columns=["date"])
    storm_df["is_storm_or_recovery"] = True
    storm_df["is_storm_onset"] = storm_df["date"].isin(storm_days)
    storm_df.to_csv(cfg.interim_dir / "storm_days.csv", index=False)

    summary = {
        "group": args.group,
        "date_range": [start, end],
        "sources": source_note,
        "confidence": confidence,
        "n_storm_onset_days": len(storm_days),
        "n_storm_or_recovery_days": len(extended),
        "storm_onset_dates": sorted(storm_days),
    }
    (cfg.interim_dir / "storm_days_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
