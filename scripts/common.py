"""Per-group configuration loader for the multi-event geomagnetic precursor
analysis pipeline.

Originally this module held fixed constants (GDMS_DIR/EQ_*/STATIONS/
SCALAR_ONLY_STATIONS/NEAR_STATIONS/FAR_STATIONS/KNOWN_OUTAGE_WINDOWS) that
only worked for the single G10 (2024-04-03 M7.2) event. `load_group_config`
replaces all of that with a per-group loader driven by `events.py` (event
metadata) and each file's own IAGA-2002 header (station metadata + scalar-
vs-vector status, via `parser.parse_header`), so it works uniformly across
all 13 groups without a hand-maintained per-era station table.
"""
from __future__ import annotations

import json
import math
import re
import sys
import tarfile
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from events import Event, get_group  # noqa: E402
from parser import DayFileRef, parse_header  # noqa: E402

_TGZ_MEMBER_RE = re.compile(r"([a-z]{3})(\d{8})dsec\.sec$")

PROJECT_DIR = Path(__file__).resolve().parent.parent
GX_DATA_ROOT = PROJECT_DIR  # moved 2026-09-14: G1..G23 now live inside geomag_precursor/, not one level up
OUTPUT_DIR = PROJECT_DIR / "output"

# Confirmed empirically by timezone_check.py (Sq diurnal curve + industrial
# noise diurnal curve both corroborate, for the G10/2024 data) that the raw
# TIME column is UTC. Same CWA IAGA-2002 source/format across all groups, so
# treated as a dataset-wide constant rather than re-verified per group.
DATA_TIMEZONE = "UTC"
LOCAL_UTC_OFFSET_HOURS = 8

# Minimum station count in a channel-type pool (F-only or XYZ) to attempt the
# near/far common-mode-regression screening method at all; below this the
# pool is marked insufficient and that method is skipped for the group.
MIN_STATIONS_FOR_METHOD = 5
N_NEAR_STATIONS = 3
N_FAR_STATIONS = 2

# A station/day with more than this fraction of samples missing is excluded
# from "clean" baseline days -- this auto-detected-from-data mechanism
# replaces the old hand-curated KNOWN_OUTAGE_WINDOWS list as the general
# per-group outage-exclusion logic (see auto_outage_dates below), since
# hand-curating a windows list for 13 groups x up to 19 stations doesn't
# scale the way it did for one group.
OUTAGE_PCT_MISSING_THRESHOLD = 0.05

# G10's originally hand-curated outage list (with per-incident notes, e.g.
# "twu ~86% missing on 2024-04-23"). Kept here as a documented historical
# record of what was manually verified for the single-event 2024 pipeline --
# no longer consulted by the general multi-group logic (auto_outage_dates
# below covers G10 the same way it covers every other group), but retained
# since G10/CLAUDE.md references it and the specific incident notes are
# useful context that pct_missing alone doesn't carry.
G10_KNOWN_OUTAGE_WINDOWS = [
    {"station": "ALL", "start": "2024-01-04 23:20:00", "end": "2024-01-04 23:59:59",
     "note": "network-wide brief outage"},
    {"station": "ALL", "start": "2024-03-02 23:20:00", "end": "2024-03-02 23:59:59",
     "note": "network-wide brief outage"},
    {"station": "ALL", "start": "2024-03-05 23:20:00", "end": "2024-03-05 23:59:59",
     "note": "network-wide brief outage"},
    {"station": "ALL", "start": "2024-03-06 23:20:00", "end": "2024-03-06 23:59:59",
     "note": "network-wide brief outage"},
    {"station": "twu", "start": "2024-04-19 00:00:00", "end": "2024-04-19 23:59:59",
     "note": "404s missing"},
    {"station": "twu", "start": "2024-04-23 00:00:00", "end": "2024-04-23 23:59:59",
     "note": "~86% missing"},
    {"station": "twu", "start": "2024-04-24 00:00:00", "end": "2024-04-25 23:59:59",
     "note": "fully missing"},
    {"station": "zbn", "start": "2024-04-13 16:10:00", "end": "2024-04-13 23:59:59",
     "note": "post-quake gap"},
    {"station": "kma", "start": "2024-01-10 00:00:00", "end": "2024-01-10 23:59:59",
     "note": "~6h gap"},
    {"station": "yhg", "start": "2024-01-15 00:00:00", "end": "2024-01-15 23:59:59",
     "note": "~11h gap"},
]


def list_tgz_files(gdms_dir: Path) -> list[Path]:
    """All GDMS batch-download archives directly inside a group folder,
    sorted by filename. This sort order is what "first .tgz wins" means
    when the same station+date turns up in more than one .tgz (their
    filenames don't reliably reflect their date-range contents, so this is
    just a deterministic tie-break, not a claim about which archive is
    "newer")."""
    return sorted(gdms_dir.glob("*.tgz"))


def _tgz_member_index(tgz_path: Path) -> list[dict]:
    """List every <station><YYYYMMDD>dsec.sec member inside a .tgz, without
    extracting any file content -- just tarfile's own member headers.
    Gzip doesn't support random access, so a full pass over a large archive
    (confirmed ~30s for a 276MB/481-file sample) is not cheap; the result is
    cached to a `<name>.tgz.idx.json` sidecar next to the archive, keyed on
    (size, mtime) so replacing a .tgz with a different file of the same name
    invalidates the cache automatically."""
    idx_path = tgz_path.with_name(tgz_path.name + ".idx.json")
    st = tgz_path.stat()
    if idx_path.exists():
        try:
            cached = json.loads(idx_path.read_text())
        except (json.JSONDecodeError, OSError):
            cached = {}
        if cached.get("size") == st.st_size and cached.get("mtime") == st.st_mtime:
            return cached["members"]

    members = []
    with tarfile.open(tgz_path, "r:gz") as tf:
        for m in tf:  # header-only walk; does not decompress member bodies
            if not m.isfile():
                continue
            match = _TGZ_MEMBER_RE.search(m.name)
            if match:
                members.append({"station": match.group(1), "date_str": match.group(2), "member": m.name})
    idx_path.write_text(json.dumps({"size": st.st_size, "mtime": st.st_mtime, "members": members}))
    return members


def list_day_refs(gdms_dir: Path, station: str = "*") -> list[DayFileRef]:
    """List DayFileRefs for a station (or all stations, with the default
    '*') across both loose .sec/.sec.gz files and any *.tgz batch archives
    directly inside gdms_dir -- the one place this precedence is defined, so
    every script that lists day files picks up .tgz support the same way.

    Precedence when the same station+date exists in more than one source:
    loose .sec/.sec.gz always wins over .tgz content (silently -- this is
    the expected steady state, not worth warning about); between multiple
    .tgz archives, the first one encountered in list_tgz_files() order wins,
    and a warning is printed so overlapping batch downloads are noticeable.
    """
    by_key: dict[tuple[str, str], DayFileRef] = {}

    pattern = "*dsec.sec*" if station == "*" else f"{station}*dsec.sec*"
    for p in sorted(gdms_dir.glob(pattern)):
        if p.suffix not in (".sec", ".gz"):
            continue
        stem = p.name[:-3] if p.name.endswith(".sec.gz") else p.name
        st, date_str = stem[:3], stem[3:11]
        key = (st, date_str)
        cur = by_key.get(key)
        # Prefer plain .sec over .sec.gz if both somehow exist.
        if cur is None or (cur.source_path.suffix == ".gz" and p.suffix == ".sec"):
            by_key[key] = DayFileRef(st, date_str, p, None)

    n_loose = len(by_key)
    n_tgz = 0
    n_dup = 0
    for tgz in list_tgz_files(gdms_dir):
        for rec in _tgz_member_index(tgz):
            st, date_str, member = rec["station"], rec["date_str"], rec["member"]
            if station != "*" and st != station:
                continue
            key = (st, date_str)
            existing = by_key.get(key)
            if existing is None:
                by_key[key] = DayFileRef(st, date_str, tgz, member)
                n_tgz += 1
            elif existing.member is not None:
                n_dup += 1
                print(
                    f"[list_day_refs] WARNING: {st}{date_str} duplicate -- "
                    f"{tgz.name}:{member} skipped, already using "
                    f"{existing.source_path.name}:{existing.member}",
                    file=sys.stderr,
                )
            # else: existing is a loose file, which always wins silently.

    if n_tgz or n_dup:
        print(
            f"[list_day_refs] {gdms_dir.name}: {n_loose} loose, {n_tgz} from .tgz, "
            f"{n_dup} duplicate tgz entries skipped",
            file=sys.stderr,
        )
    return sorted(by_key.values(), key=lambda r: (r.station, r.date_str))


def resolve_day_ref(gdms_dir: Path, station: str, date_str: str) -> DayFileRef | None:
    """Return the DayFileRef for one station-day, preferring plain .sec if
    it exists, else .sec.gz, else searching every .tgz in gdms_dir; None if
    not found anywhere."""
    plain = gdms_dir / f"{station}{date_str}dsec.sec"
    if plain.exists():
        return DayFileRef(station, date_str, plain, None)
    gz = gdms_dir / f"{station}{date_str}dsec.sec.gz"
    if gz.exists():
        return DayFileRef(station, date_str, gz, None)
    for tgz in list_tgz_files(gdms_dir):
        for rec in _tgz_member_index(tgz):
            if rec["station"] == station and rec["date_str"] == date_str:
                return DayFileRef(station, date_str, tgz, rec["member"])
    return None


def to_local_hour(utc_index):
    """Given a UTC DatetimeIndex, return the Taiwan local (UTC+8) hour-of-day (0-23)."""
    import pandas as pd

    return (utc_index + pd.Timedelta(hours=LOCAL_UTC_OFFSET_HOURS)).hour


def haversine_km(lat1, lon1, lat2, lon2):
    r = 6371.0
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlambda = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dlambda / 2) ** 2
    return 2 * r * math.asin(min(1.0, a ** 0.5))


def auto_outage_dates(daily_features_df, threshold: float = OUTAGE_PCT_MISSING_THRESHOLD) -> dict[str, set[str]]:
    """Per-station set of dates excluded from 'clean' baseline days, derived
    from daily_features.csv's own pct_missing column rather than a
    hand-curated list -- the general per-group replacement for the old
    G10-only KNOWN_OUTAGE_WINDOWS mechanism. daily_features_df must have
    'station', 'date', 'pct_missing' columns (as written by
    build_daily_features.py)."""
    dirty = daily_features_df[daily_features_df["pct_missing"] > threshold]
    out: dict[str, set[str]] = {}
    for station, g in dirty.groupby("station"):
        out[station] = set(g["date"])
    return out


@dataclass(frozen=True)
class StationPool:
    """Near/far station split for one channel-type pool (F-only or XYZ) within one group."""

    channel: str  # "F" or "XYZ"
    all_stations: tuple[str, ...]  # every station code in this pool, present in the group, nearest-first
    near: tuple[str, ...]
    far: tuple[str, ...]
    sufficient: bool  # True if this pool has enough stations to attempt the near/far method


@dataclass(frozen=True)
class GroupConfig:
    group_id: str
    gdms_dir: Path
    interim_dir: Path
    external_dir: Path
    anchor_event: Event
    all_events: tuple[Event, ...]
    stations: dict  # station code -> {"name","lat","lon","elevation_m","reported","distance_km"}
    f_pool: StationPool
    xyz_pool: StationPool


def _discover_stations(gdms_dir: Path, anchor_event: Event) -> dict:
    """Glob distinct station-code prefixes present in this group's folder and
    read one representative file's header per code for name/lat/lon/
    elevation/Reported status. Scalar-vs-vector status is constant within a
    single group -- confirmed empirically across all 13 groups, the F->XYZF
    upgrade always lands on a group-era boundary, never mid-group -- so one
    representative file per code is sufficient."""
    codes = sorted({r.station for r in list_day_refs(gdms_dir)})
    stations = {}
    for code in codes:
        sample = list_day_refs(gdms_dir, code)[0]
        h = parse_header(sample)
        stations[code] = {
            "name": h["station_name"],
            "lat": h["lat"],
            "lon": h["lon"],
            "elevation_m": h["elevation_m"],
            "reported": h["reported"].strip().upper(),
            "distance_km": round(haversine_km(h["lat"], h["lon"], anchor_event.lat, anchor_event.lon), 1),
        }
    return stations


def _build_pool(stations: dict, channel: str) -> StationPool:
    wanted_reported = "F" if channel == "F" else "XYZF"
    ranked = sorted(
        (c for c, m in stations.items() if m["reported"] == wanted_reported),
        key=lambda c: stations[c]["distance_km"],
    )
    near = tuple(ranked[:N_NEAR_STATIONS])
    far = tuple(c for c in ranked[-N_FAR_STATIONS:] if c not in near)
    return StationPool(
        channel=channel,
        all_stations=tuple(ranked),
        near=near,
        far=far,
        sufficient=len(ranked) >= MIN_STATIONS_FOR_METHOD and len(near) >= 1 and len(far) >= 1,
    )


def load_group_config(group_id: str) -> GroupConfig:
    group = get_group(group_id)
    gdms_dir = GX_DATA_ROOT / group.folder
    interim_dir = PROJECT_DIR / "data" / "interim" / group_id
    external_dir = PROJECT_DIR / "data" / "external" / group_id
    interim_dir.mkdir(parents=True, exist_ok=True)
    external_dir.mkdir(parents=True, exist_ok=True)

    anchor = group.anchor_event
    stations = _discover_stations(gdms_dir, anchor)
    f_pool = _build_pool(stations, "F")
    xyz_pool = _build_pool(stations, "XYZ")

    return GroupConfig(
        group_id=group_id,
        gdms_dir=gdms_dir,
        interim_dir=interim_dir,
        external_dir=external_dir,
        anchor_event=anchor,
        all_events=group.events,
        stations=stations,
        f_pool=f_pool,
        xyz_pool=xyz_pool,
    )


if __name__ == "__main__":
    group_ids = sys.argv[1:] or ["G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9", "G10", "G11", "G12", "G13"]
    for gid in group_ids:
        cfg = load_group_config(gid)
        a = cfg.anchor_event
        print(f"== {gid} == anchor {a.date} {a.magnitude_type}{a.magnitude} ({a.lat},{a.lon}) "
              f"confidence={a.coord_confidence}")
        print(f"  F   pool: {len(cfg.f_pool.all_stations):2d} stations  near={cfg.f_pool.near}  "
              f"far={cfg.f_pool.far}  sufficient={cfg.f_pool.sufficient}")
        print(f"  XYZ pool: {len(cfg.xyz_pool.all_stations):2d} stations  near={cfg.xyz_pool.near}  "
              f"far={cfg.xyz_pool.far}  sufficient={cfg.xyz_pool.sufficient}")
