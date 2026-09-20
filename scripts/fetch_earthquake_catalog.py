"""Fetch an EXTENDED earthquake catalog near the geomagnetic station network,
covering every M>=6.0 event AND the smaller ones in between -- this is the
"把一次地震變成很多次" (turn one earthquake into many) piece of the professor's
suggested methodology: `events.py` only has the 20 hand-picked M>=6.0 events
used to define the 13 groups, which is not enough for a real superposed-epoch
stack or a meaningful backtest population.

Source: USGS FDSN event webservice (public, no auth, well-documented GeoJSON
format) -- https://earthquake.usgs.gov/fdsnws/event/1/query. CWA's own open
data API was not used here since no existing code in this repo already wraps
it and USGS covers the same events (cross-referenced against events.py's
USGS ids/magnitudes already for several entries), so reusing USGS keeps this
consistent with what's already partially in the dataset's own provenance
notes.

Scope, per plan (`~/.claude/plans/block-bootstrap-shiny-pearl.md`):
  - Only the 8 groups that have `ulf_near_far_index.csv` (vector-sufficient
    XYZ pool; G1/G2/G3 are scalar-only and can't run the polarization method
    at all, so they're not part of the superposed-epoch/backtest population).
  - Per-group date window = exactly the window that group's ULF index
    actually covers (read from the CSV itself, not re-derived from
    events.py/docs, so this never drifts from what data is really usable).
  - Bounding box covers the whole station network with margin: lat
    20.5-27 (network spans hcn 21.94N to mtu 26.17N), lon 117.5-123.5
    (kma 118.35E to offshore events out past 122E).
  - Declustering: sort by time; an event is dropped if it falls within
    `decluster_days` days AND `decluster_km` km of an already-kept event of
    >= its own magnitude (i.e. keep the largest event of each tight
    space-time cluster) -- avoids treating an aftershock sequence as many
    independent samples, the same pseudo-replication concern events.py's own
    docstring raises for the 13-group design.
  - Cross-referenced against events.py's 20 curated events (within 6h /
    0.3 magnitude) and flagged `is_known_event` rather than duplicated.

Usage:
  fetch_earthquake_catalog.py --min-mag 5.5
  fetch_earthquake_catalog.py --min-mag 5.0 --output data/external/extended_catalog_m5.0.csv
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import PROJECT_DIR, haversine_km, load_group_config  # noqa: E402
from events import ALL_GROUP_IDS, assign_group_for_time, folder_events, get_group, sibling_group_ids  # noqa: E402

USGS_URL = "https://earthquake.usgs.gov/fdsnws/event/1/query"
BBOX = {"minlatitude": 20.5, "maxlatitude": 27.0, "minlongitude": 117.5, "maxlongitude": 123.5}
ULF_GROUPS = ("G4", "G5", "G6", "G7", "G8", "G9", "G10", "G11", "G12", "G13", "G19", "G20", "G23")

DECLUSTER_DAYS = 3
DECLUSTER_KM = 100
KNOWN_EVENT_HOURS = 6
KNOWN_EVENT_MAG_TOL = 0.3


def group_date_window(group_id: str) -> tuple[str, str] | None:
    cfg = load_group_config(group_id)
    path = cfg.interim_dir / "ulf_near_far_index.csv"
    if not path.exists():
        return None
    dates = []
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            dates.append(row["date"])
    if not dates:
        return None
    dates.sort()
    d0, d1 = dates[0], dates[-1]
    return (f"{d0[:4]}-{d0[4:6]}-{d0[6:8]}", f"{d1[:4]}-{d1[4:6]}-{d1[6:8]}")


def fetch_usgs(start: str, end: str, min_mag: float) -> list[dict]:
    params = {**BBOX, "format": "geojson", "starttime": start, "endtime": end, "minmagnitude": min_mag}
    url = USGS_URL + "?" + "&".join(f"{k}={v}" for k, v in params.items())
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                data = json.load(r)
            break
        except Exception as exc:  # noqa: BLE001
            if attempt == 2:
                raise
            print(f"  retry after error: {exc}", file=sys.stderr)
            time.sleep(2)
    events = []
    for feat in data["features"]:
        p = feat["properties"]
        lon, lat, depth = feat["geometry"]["coordinates"]
        t = datetime.fromtimestamp(p["time"] / 1000, tz=timezone.utc)
        events.append({
            "time_utc": t.strftime("%Y-%m-%d %H:%M:%S"),
            "_dt": t,
            "lat": lat, "lon": lon, "depth_km": depth,
            "mag": p["mag"], "place": p.get("place") or "",
            "usgs_id": feat.get("id", ""),
        })
    return events


def decluster(events: list[dict], days: float = DECLUSTER_DAYS, km: float = DECLUSTER_KM) -> list[dict]:
    kept: list[dict] = []
    for e in sorted(events, key=lambda x: x["_dt"]):
        suppressed = False
        for k in kept:
            dt_days = abs((e["_dt"] - k["_dt"]).total_seconds()) / 86400
            if dt_days > days:
                continue
            if haversine_km(e["lat"], e["lon"], k["lat"], k["lon"]) > km:
                continue
            if k["mag"] >= e["mag"]:
                suppressed = True
                break
        if not suppressed:
            kept.append(e)
    kept_set = {id(k) for k in kept}
    for e in events:
        e["declustered"] = id(e) in kept_set
    return events


def flag_known_events(events: list[dict], group_id: str) -> None:
    # folder_events: a sibling group's registered event is just as "known" as this group's own.
    known = folder_events(group_id)
    for e in events:
        e["is_known_event"] = False
        for ev in known:
            known_dt = datetime.strptime(ev.time_utc, "%Y-%m-%d %H:%M:%S").replace(tzinfo=timezone.utc)
            hours = abs((e["_dt"] - known_dt).total_seconds()) / 3600
            if hours <= KNOWN_EVENT_HOURS and abs(e["mag"] - ev.magnitude) <= KNOWN_EVENT_MAG_TOL:
                e["is_known_event"] = True
                break


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-mag", type=float, default=5.5)
    ap.add_argument("--output", type=Path, default=None)
    ap.add_argument("--groups", nargs="*", default=list(ULF_GROUPS))
    args = ap.parse_args()

    out_path = args.output or (PROJECT_DIR / "data" / "external" / f"extended_catalog_m{args.min_mag}.csv")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    all_rows = []
    for group_id in args.groups:
        window = group_date_window(group_id)
        if window is None:
            print(f"[{group_id}] no ulf_near_far_index.csv -- skipping (scalar-only or not run)", file=sys.stderr)
            continue
        start, end = window
        # USGS endtime is exclusive-ish at day boundary in practice; pad by 1 day
        end_padded = (datetime.strptime(end, "%Y-%m-%d")).strftime("%Y-%m-%d")
        print(f"[{group_id}] querying USGS {start} ~ {end_padded}, M>={args.min_mag}", file=sys.stderr)
        events = fetch_usgs(start, end_padded, args.min_mag)
        events = decluster(events)
        flag_known_events(events, group_id)
        if len(sibling_group_ids(group_id)) > 1:
            # Sibling groups (shared raw-data folder) all query this same date range: keep each
            # catalog event in exactly one of them (decluster/known-flagging above ran on the
            # whole window first, so a foreshock/aftershock pair straddling two groups is still
            # declustered together), otherwise it would be counted once per sibling.
            n_all = len(events)
            events = [e for e in events if assign_group_for_time(group_id, e["_dt"]) == group_id]
            print(f"  folder shared with {sibling_group_ids(group_id)}: kept {len(events)}/{n_all} "
                  f"events nearest {group_id}'s anchor", file=sys.stderr)
        n_kept = sum(1 for e in events if e["declustered"])
        print(f"  {len(events)} raw -> {n_kept} declustered", file=sys.stderr)
        for e in events:
            e["group"] = group_id
            all_rows.append(e)

    fields = ["group", "time_utc", "lat", "lon", "depth_km", "mag", "place", "usgs_id",
              "declustered", "is_known_event"]
    with out_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in all_rows:
            w.writerow({k: r[k] for k in fields})

    n_total = len(all_rows)
    n_decl = sum(1 for r in all_rows if r["declustered"])
    print(f"wrote {out_path} ({n_total} raw rows, {n_decl} declustered/independent events)", file=sys.stderr)


if __name__ == "__main__":
    main()
