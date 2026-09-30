"""Fetch an EXTENDED earthquake catalog near the geomagnetic station network,
covering every M>=6.0 event AND the smaller ones in between -- this is the
"把一次地震變成很多次" (turn one earthquake into many) piece of the professor's
suggested methodology: `events.py` only has the 20 hand-picked M>=6.0 events
used to define the 13 groups, which is not enough for a real superposed-epoch
stack or a meaningful backtest population.

Source: the user-supplied CWA GDMS catalog exports (ML, UTC) -- this project
takes CWA as the authority for every event and magnitude. Two exports tile
2009-01-01 ~ 2026-07-31 (CWA_CATALOGS); a group's window is read from both, so
a window straddling the files (G19) is still pure CWA. The USGS FDSN query is
kept only behind --allow-usgs for windows outside that range: before
2026-09-30 it was the silent fallback for every group outside 2024-09~2026-07,
which put Mw/mb magnitudes into the M>=5.0/5.5 tiers.

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
  - Cross-referenced against events.py's registered events (origin times
    within KNOWN_EVENT_SEC and epicenters within KNOWN_EVENT_KM; magnitude
    is NOT compared, since USGS and CWA ML routinely differ by more than
    0.3 for the same earthquake) and flagged `is_known_event` rather than
    duplicated.

Usage:
  fetch_earthquake_catalog.py --min-mag 5.5
  fetch_earthquake_catalog.py --min-mag 5.0 --output data/external/extended_catalog_m5.0.csv
  fetch_earthquake_catalog.py --min-mag 5.5 --reflag   # recompute is_known_event only, offline
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
import urllib.request
from datetime import datetime, timedelta, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import PROJECT_DIR, haversine_km, load_group_config  # noqa: E402
from events import ALL_GROUP_IDS, assign_group_for_time, folder_events, get_group, sibling_group_ids  # noqa: E402

USGS_URL = "https://earthquake.usgs.gov/fdsnws/event/1/query"
BBOX = {"minlatitude": 20.5, "maxlatitude": 27.0, "minlongitude": 117.5, "maxlongitude": 123.5}
ULF_GROUPS = ("G4", "G5", "G6", "G7", "G8", "G9", "G10", "G11", "G12", "G13", "G19", "G20", "G23", "G24")

DECLUSTER_DAYS = 3
DECLUSTER_KM = 100
# Same earthquake in the catalog and in events.py: origin times within a minute and epicenters
# within 50 km. The old rule (6 h, magnitude within 0.3) missed duplicates whose USGS magnitude
# differs from CWA ML by more than 0.3 (G4 2020-12-10: 6.1 vs 6.64, 0 s apart) and matched
# distinct aftershocks hours away from a registered event.
KNOWN_EVENT_SEC = 60
KNOWN_EVENT_KM = 50

# User-supplied CWA GDMS regional magnitude-report exports (space-delimited, header
# "date time lat lon depth ML nstn dmin gap trms ERH ERZ fixed nph quality"), M>=5.0,
# as (path, first day, last day) of each export's requested range. Both are UTC:
# GDMScatalog.txt by G11's 2025-01-21 anchor (its 2025-01-20 16:17 UTC row only lands on
# 01-21 in Taiwan local time); GDMScatalog_2009-2024.txt by all 62 events.py events before
# 2024-09 matching a row within 2 s. The exports cover a wider area than BBOX (out to lon
# 125.6), so fetch_cwa() applies BBOX itself.
CWA_CATALOGS = (
    (PROJECT_DIR.parent / "GDMScatalog_2009-2024.txt", "2009-01-01", "2024-08-31"),
    (PROJECT_DIR.parent / "GDMScatalog.txt", "2024-09-01", "2026-07-31"),
)
CWA_CATALOG_START = CWA_CATALOGS[0][1]
CWA_CATALOG_END = CWA_CATALOGS[-1][2]


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


def in_bbox(lat: float, lon: float) -> bool:
    return (BBOX["minlatitude"] <= lat <= BBOX["maxlatitude"]
            and BBOX["minlongitude"] <= lon <= BBOX["maxlongitude"])


def fetch_cwa_all(paths: list[Path], start: str, end: str, min_mag: float) -> list[dict]:
    """fetch_cwa() over every export, so a window straddling two of them is read whole."""
    return [e for path in paths for e in fetch_cwa(path, start, end, min_mag)]


def fetch_cwa(path: Path, start: str, end: str, min_mag: float) -> list[dict]:
    """Parse one user-supplied CWA GDMS catalog export, filtered to [start, end]
    (inclusive, YYYY-MM-DD), mag >= min_mag and BBOX (the same box fetch_usgs() queries).
    Returns the same dict shape as fetch_usgs() so decluster()/flag_known_events() work
    unchanged."""
    start_dt = datetime.strptime(start, "%Y-%m-%d").replace(tzinfo=timezone.utc)
    end_dt = datetime.strptime(end, "%Y-%m-%d").replace(tzinfo=timezone.utc) + timedelta(days=1)
    events = []
    with path.open() as f:
        header = f.readline()
        assert header.split()[:6] == ["date", "time", "lat", "lon", "depth", "ML"], \
            f"unexpected GDMS catalog header: {header!r}"
        for line in f:
            parts = line.split()
            if not parts:
                continue
            date_s, time_s, lat_s, lon_s, depth_s, mag_s = parts[:6]
            mag = float(mag_s)
            if mag < min_mag:
                continue
            t = datetime.strptime(f"{date_s} {time_s}", "%Y-%m-%d %H:%M:%S.%f").replace(tzinfo=timezone.utc)
            if not (start_dt <= t < end_dt):
                continue
            if not in_bbox(float(lat_s), float(lon_s)):
                continue
            events.append({
                "time_utc": t.strftime("%Y-%m-%d %H:%M:%S"),
                "_dt": t,
                "lat": float(lat_s), "lon": float(lon_s), "depth_km": float(depth_s),
                "mag": mag, "place": "", "usgs_id": "",
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
            if (abs((e["_dt"] - known_dt).total_seconds()) <= KNOWN_EVENT_SEC
                    and haversine_km(e["lat"], e["lon"], ev.lat, ev.lon) <= KNOWN_EVENT_KM):
                e["is_known_event"] = True
                break


def reflag_existing(path: Path) -> None:
    """Recompute only is_known_event on an already-fetched catalog CSV, offline: rows,
    declustering and group assignment stay exactly as fetched."""
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        fields = reader.fieldnames
        rows = list(reader)
    n_changed = 0
    for group_id in dict.fromkeys(r["group"] for r in rows):
        grp = [r for r in rows if r["group"] == group_id]
        evs = [{"_dt": datetime.strptime(r["time_utc"][:19], "%Y-%m-%d %H:%M:%S").replace(tzinfo=timezone.utc),
                "lat": float(r["lat"]), "lon": float(r["lon"])} for r in grp]
        flag_known_events(evs, group_id)
        for r, e in zip(grp, evs):
            new = str(e["is_known_event"])
            if r["is_known_event"] != new:
                n_changed += 1
                print(f"  {group_id} {r['time_utc']} M{r['mag']} {r['source']}: "
                      f"is_known_event {r['is_known_event']} -> {new}", file=sys.stderr)
                r["is_known_event"] = new
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"reflagged {path}: {n_changed} rows changed", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-mag", type=float, default=5.5)
    ap.add_argument("--output", type=Path, default=None)
    ap.add_argument("--groups", nargs="*", default=list(ULF_GROUPS))
    ap.add_argument("--allow-usgs", action="store_true",
                    help="query USGS for a group whose window falls outside the CWA exports' "
                         f"range ({CWA_CATALOG_START} ~ {CWA_CATALOG_END}); without it such a "
                         "group is an error, never a silent fallback")
    ap.add_argument("--reflag", action="store_true",
                    help="only recompute is_known_event on the existing output CSV (no fetching)")
    args = ap.parse_args()

    out_path = args.output or (PROJECT_DIR / "data" / "external" / f"extended_catalog_m{args.min_mag}.csv")
    if args.reflag:
        reflag_existing(out_path)
        return
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
        if CWA_CATALOG_START <= start and end <= CWA_CATALOG_END:
            paths = [p for p, _, _ in CWA_CATALOGS]
            missing = [str(p) for p in paths if not p.exists()]
            if missing:
                sys.exit(f"missing CWA catalog export(s): {missing}")
            print(f"[{group_id}] reading CWA GDMS catalogs {start} ~ {end_padded}, M>={args.min_mag}", file=sys.stderr)
            events = fetch_cwa_all(paths, start, end_padded, args.min_mag)
            source = "CWA_GDMS"
        elif not args.allow_usgs:
            sys.exit(f"[{group_id}] window {start} ~ {end} is outside the CWA exports "
                     f"({CWA_CATALOG_START} ~ {CWA_CATALOG_END}); add a CWA export or pass --allow-usgs")
        else:
            print(f"[{group_id}] querying USGS {start} ~ {end_padded}, M>={args.min_mag}", file=sys.stderr)
            events = fetch_usgs(start, end_padded, args.min_mag)
            source = "USGS"
        for e in events:
            e["source"] = source
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
              "declustered", "is_known_event", "source"]
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
