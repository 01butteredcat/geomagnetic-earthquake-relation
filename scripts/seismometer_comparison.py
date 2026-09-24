"""Co-located/nearby seismometer-vs-geomagnetic comparison for the coseismic
anomalies `coseismic_step_analysis.py` / `coseismic_stacking_analysis.py`
found near each earthquake's origin second.

## Why this exists

Neither of those two scripts can tell "the magnetometer's housing got
physically shaken" apart from "the magnetic field itself changed" -- both
produce the same signature (a step/spike near the origin second) in 1Hz
magnetic data alone. This script adds the one thing that CAN separate them:
independent ground-motion data from a real seismometer/accelerometer, at or
near the same site. The logic (per the user's own framing, which this
script implements directly):

  - If the geomagnetic anomaly's onset and duration line up with the
    seismometer's own strong-motion window, that's consistent with
    (though doesn't prove) shaking-induced instrumental noise.
  - If the geomagnetic anomaly starts measurably BEFORE the ground starts
    moving, or PERSISTS after the ground motion has died down, that favors
    a real geophysical mechanism (piezomagnetic effect, etc.) instead.

## Data used (fetched separately by the user, not by this pipeline)

31 of the 117 events in `events.py` have SAC PoleZero instrument-response
files (`<GROUP_MMDD>/SAC_PZs_TW_<STA>_<CHAN>_...`) and a matching miniSEED
waveform file (`<GROUP_MMDD>/<GROUP_MMDD>_w.mseed`, ~event_utc-60s to
event_utc+600s, 100Hz), both under `seismometer/<GROUP_MMDD>/` -- see
`SEISMIC_DATA_DIRS` below, keyed by (group_id, event.date). (Reorganized
2026-08-17: previously the mseed files sat loose directly under this
project's parent directory, one level up from their matching PZ folder;
both now live together under `seismometer/`, one folder per event, out of
the way of the Gx geomagnetic folders.) Of the 19 uncovered events, 6 are a
permanent, structural gap: G14's 5 (its earliest event, 2009-07-14, predates
the seismic data source's 2012 cutoff) plus G21's 1 (2010-11-21, same reason
-- confirmed 2026-09-19 that the data source's fetchable range only starts
2012-01-01). The remaining 13 (spread across G5, G8, G9, G10, G11, G12,
G15, G17) are simply not yet fetched. (G22's 1 event and G23's 1 event + G24's
1 event -- G23/G24 split 2026-09-22 from a single merged "G23" that had 2
events -- had mseed fetched 2026-09-19 and PZ files added 2026-09-20, and are
now wired in;
the PZ sets lack a few stations that appear in the mseed -- G22: CHK/ELD loc 11
and HEN, G23: HEN and SSH -- which only affects those stations' traces, none
of which were the nearest-station pick.)
(As of 2026-08-13/14, only the 16 anchor events had been fetched; the
remaining 11 non-anchor events plus the previously-missing G9 2022-09-17
foreshock were fetched and verified 2026-08-16 -- see `coverage_summary.json`,
which enumerates all events in `events.py`'s current registry with an honest
`seismic_data_status` so the report never implies more coverage than it has.)
(2026-09-23: registry grew from 49 to 117 events, all 68 new ones non-anchor
and none seismometer-fetched, so the "19 uncovered" breakdown above is stale
-- it describes only the original 19, not the ~86 uncovered now; see
`docs/candidate_events_gdms_2024_2026.md` for the new batch's own provenance.
One of the new events, G11's 2025-01-21b, turned out to already have real
data -- it falls inside the 2025-01-21 anchor's already-fetched mseed window
-- and was wired into `SEISMIC_DATA_DIRS` below; see that key's comment.)

(2026-09-25: the remaining 80 fetchable events were fetched in one batch --
14 M>=6, 66 M5 -- and laid out under one naming rule that
`_register_convention_dirs()` derives rather than lists. Coverage is now 111 of
117; the other 6 predate the data source (G14's 5, G21's 1) and report
`no_data_pre_2012`. The batch carried no PoleZero files; each new folder holds
copies of the existing folders' PZ files whose epoch covers the event, since
same-named PZ files differ only in their CREATED line.)

Empirically confirmed quirks this module works around (see functions below
for where): (1) the PZ-folder-name <-> mseed-filename mapping is NOT a
derivable pattern -- G9's PZ folder is `G9_0918` but its mseed file is
`G09_0918_w.mseed` -- so SEISMIC_DATA_DIRS below is a small hardcoded table,
matching this codebase's own events.py precedent of hardcoding rather than
inferring; (2) miniSEED traces are NOT uniformly -60s/+600s -- some stations
(e.g. G9's `ECS`, otherwise an excellent 1.9km co-location candidate for the
csg persistent-offset case) are short triggered-accelerograph recordings
ending well before +600s, handled via the `triggered_short_trace` flag
rather than assumed away; (3) some (station, channel) pairs have duplicate
traces in the same mseed file (a gap-split recording or two nearby
triggers), handled via `Stream.merge()` rather than a naive `select()[0]`.

## Method

Two independent strong-motion window detectors are computed and both kept
(not one silently preferred), because the confirmed short-triggered-trace
case doesn't have enough pre-event baseline for STA/LTA to work:
  - `detect_window_sta_lta`: classic STA/LTA trigger (needs a pre-event
    baseline; returns `no_pre_event_baseline` rather than a fabricated
    result when the trace is too short).
  - `detect_window_envelope_threshold`: bandpass + envelope, thresholded
    relative to its own peak -- works even on short triggered traces, at
    the cost of being threshold-sensitive.
The geomagnetic side reuses `coseismic_step_analysis.py`'s own detrend +
step30-statistic machinery directly (not `coseismic_stacking_analysis.py`'s
cached output -- the two new scripts have no run-order dependency), z-scored
against that event's own off-event noise floor exactly like the stacking
script.

Usage:
  seismometer_comparison.py --self-test                         # PZ-parse + response-removal sanity check
  seismometer_comparison.py --group G9 --geomag-station csg      # spot check the flagship G9 case
  seismometer_comparison.py --all                                 # every group with seismic data fetched
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import common  # noqa: E402
from events import GROUPS, folder_events, get_group  # noqa: E402
from stat_utils import mad_zscore  # noqa: E402
from coseismic_step_analysis import (  # noqa: E402
    DETREND_WINDOW_SEC,
    EXCLUSION_BUFFER_SEC,
    _build_channels,
    _detrend,
    _load_station_days,
    _pos,
    _rank_stations_for_event,
    _step_statistic,
)

# Reorganized 2026-08-17: each GXX_MMDD folder (with its mseed file moved
# inside it, no longer a loose <GROUP_MMDD>_w.mseed at DATA_ROOT) now lives
# under this subdirectory instead of directly alongside the Gx geomagnetic
# folders.
# Moved 2026-09-14: seismometer/ now lives inside geomag_precursor/ itself
# (alongside Gx, which moved the same way earlier that day), not one level
# up at common.DATA_ROOT anymore -- hence GX_DATA_ROOT (== PROJECT_DIR), not
# DATA_ROOT, below.
SEISMIC_ROOT = common.GX_DATA_ROOT / "seismometer"

# Hardcoded because the folder-name <-> group-id mapping is not reliably
# derivable -- G9's mismatch (pz_dir "G9_0918" vs mseed "G09_0918_w.mseed")
# is the concrete counterexample. G14 intentionally absent: its 2009-12-19
# anchor predates the seismic data source's 2012 cutoff.
#
# Keyed by (group_id, event.date) -- as of 2026-09-20 this covered 30 of the
# 49 events then in the registry (all except G14, G21 and 13 not-yet-fetched; not just the 16
# anchors from the original single-event-per-group fetch), confirmed present on disk with mseed windows correctly
# bracketing each event's origin second. G10's 2024-04-23a/2024-04-23b share
# one PZ folder (same UTC calendar day, station metadata doesn't change
# minute to minute) but have their own separate mseed files -- the same a/b/c-suffixed-`date`
# pattern used below for any other same-group, same-calendar-day events (see events.py's
# 2026-09-23 note on why: `event.date` is used as a same-group unique key here and in
# coseismic_step_analysis.py/coseismic_stacking_analysis.py's per-event output naming, which a
# 2026-09-23 batch of 68 new non-anchor events -- registry now 117 events -- would otherwise
# silently collide on for 7 same-day pairs across G11/G13/G20; all suffixed at the source in
# events.py rather than worked around here). 31 of 117 events have real seismic data as of
# 2026-09-23 (30 pre-existing + 2025-01-21b, recovered from the anchor's already-fetched window
# below -- see that key's own comment).
SEISMIC_DATA_DIRS: dict[str, dict[str, dict]] = {
    "G1": {
        "2018-02-04": {"pz_dir": "G01_0204", "mseed": "G01_0204_w.mseed"},
        "2018-02-06": {"pz_dir": "G01_0206", "mseed": "G01_0206_w.mseed"},
    },
    # G2/G3 and G6/G7/G8 were one merged group each until the 2026-09-20 split; the seismometer
    # folder/file names (fetched under the old merged names) are unchanged.
    "G2": {
        "2019-04-18": {"pz_dir": "G02_G03_0418", "mseed": "G02_G03_0418_w.mseed"},
    },
    "G3": {
        "2019-08-08": {"pz_dir": "G02_G03_0808", "mseed": "G02_G03_0808_w.mseed"},
    },
    "G4": {
        "2020-12-10": {"pz_dir": "G04_1210", "mseed": "G04_1210_w.mseed"},
    },
    "G5": {
        "2021-04-18": {"pz_dir": "G05_0418", "mseed": "G05_0418_w.mseed"},
    },
    "G6": {
        "2021-10-24": {"pz_dir": "G06_G07_G08_1024", "mseed": "G06_G07_G08_1024_w.mseed"},
    },
    "G7": {
        "2022-01-03": {"pz_dir": "G06_G07_G08_0103", "mseed": "G06_G07_G08_0103_w.mseed"},
    },
    "G8": {
        "2022-03-23": {"pz_dir": "G06_G07_G08_0323", "mseed": "G06_G07_G08_0323_w.mseed"},
    },
    "G9": {
        "2022-09-17": {"pz_dir": "G09_0917", "mseed": "G09_0917_w.mseed"},
        "2022-09-18": {"pz_dir": "G9_0918", "mseed": "G09_0918_w.mseed"},
    },
    "G10": {
        "2024-04-03": {"pz_dir": "G10_0403", "mseed": "G10_0403_w.mseed"},
        "2024-04-23a": {"pz_dir": "G10_0422", "mseed": "G10_0422_w.mseed"},
        "2024-04-23b": {"pz_dir": "G10_0422", "mseed": "G10_0422_02_w.mseed"},
        "2024-05-10": {"pz_dir": "G10_0510", "mseed": "G10_0510_w.mseed"},
    },
    "G11": {
        "2025-01-21": {"pz_dir": "G11_0121", "mseed": "G11_0121_w.mseed"},
        # 2025-01-21b (00:26:25 local = 2025-01-20 16:26:25 UTC) falls inside the anchor's
        # already-fetched event_utc-60s~+600s window (16:16:26~16:27:26 UTC) -- same mseed file,
        # genuinely covers this event too. 2025-01-21c (01:42:31 local) does NOT (>1h outside the
        # window) and is correctly left unmapped here (reports not_fetched).
        "2025-01-21b": {"pz_dir": "G11_0121", "mseed": "G11_0121_w.mseed"},
    },
    "G12": {
        "2025-08-27": {"pz_dir": "G12_0827", "mseed": "G12_0827_w.mseed"},
    },
    "G13": {
        "2026-05-01": {"pz_dir": "G13_0501", "mseed": "G13_0501_w.mseed"},
    },
    "G15": {
        "2013-06-02": {"pz_dir": "G15_0602", "mseed": "G15_0602_w.mseed"},
    },
    "G16": {
        "2013-10-31": {"pz_dir": "G16_1031", "mseed": "G16_1031_w.mseed"},
    },
    "G17": {
        "2014-12-11": {"pz_dir": "G17_1211", "mseed": "G17_1211_w.mseed"},
        "2015-02-14": {"pz_dir": "G17_0214", "mseed": "G17_0214_w.mseed"},
    },
    "G18": {
        "2016-02-06": {"pz_dir": "G18_0206", "mseed": "G18_0206_w.mseed"},
        "2016-05-31": {"pz_dir": "G18_0531", "mseed": "G18_0531_w.mseed"},
    },
    "G19": {
        "2024-08-16": {"pz_dir": "G19_0816", "mseed": "G19_0816_w.mseed"},
    },
    "G20": {
        "2025-12-24": {"pz_dir": "G20_1224", "mseed": "G20_1224_w.mseed"},
        "2025-12-27": {"pz_dir": "G20_1227", "mseed": "G20_1227_w.mseed"},
    },
    "G22": {
        "2012-06-10": {"pz_dir": "G22_0610", "mseed": "G22_0610_w.mseed"},
    },
    "G23": {
        "2020-07-26": {"pz_dir": "G23_0726", "mseed": "G23_0726_w.mseed"},
    },
    "G24": {
        "2020-06-14": {"pz_dir": "G23_0614", "mseed": "G23_0614_w.mseed"},
    },
}

# Everything fetched from 2026-09-25 on follows one naming rule, so it is derived
# instead of listed: G<2-digit group>_<MMDD of event.date, i.e. Taiwan local
# date><a/b/c suffix if any>/, holding <that name>_w.mseed plus the SAC PoleZero
# files whose epoch covers the event. (The hand-listed entries above predate the
# rule: mixed zero-padding, some UTC-dated folders, pre-split G06_G07_G08 names.)
# The 2026-09-25 batch added 80 events this way (14 M>=6, 66 M5).
SEISMIC_DATA_SOURCE_START_UTC = pd.Timestamp("2012-01-01")  # nothing earlier is fetchable


def convention_dir_name(group_id: str, event_date: str) -> str:
    m = re.fullmatch(r"\d{4}-(\d{2})-(\d{2})([a-z]?)", event_date)
    return f"G{int(group_id[1:]):02d}_{m.group(1)}{m.group(2)}{m.group(3)}"


def _register_convention_dirs() -> None:
    for group_id, group in GROUPS.items():
        for event in group.events:
            if event.date in SEISMIC_DATA_DIRS.get(group_id, {}):
                continue
            name = convention_dir_name(group_id, event.date)
            if (SEISMIC_ROOT / name / f"{name}_w.mseed").exists():
                SEISMIC_DATA_DIRS.setdefault(group_id, {})[event.date] = {"pz_dir": name, "mseed": f"{name}_w.mseed"}


_register_convention_dirs()

GEOMAG_HALF_SEC = 240        # target: window for the geomagnetic side's z-scored step30 profile
SEARCH_HALF_SEC = 180        # target: peak-search sub-window, same convention as coseismic_step_analysis.py's
                              # SCAN_HALF_SEC (widened 120->180 2026-08-16 -- was pinning G12/G13/G20's
                              # obs_lag at the old +-120s boundary. 2026-08-19: investigated widening
                              # further to 300/360s to check whether G20's 2025-12-24 obs_lag=-178s (2s
                              # from this 180s boundary) was a truncation artifact -- see plan
                              # artifact-wobbly-kettle.md. It was NOT (obs_lag stayed at -178s at 300s),
                              # but the wider window changed several OTHER, non-boundary-pinned events'
                              # reported obs_lag (a wider search finds a larger max-of-N by chance even
                              # under noise alone; see coseismic_step_analysis.py's SCAN_HALF_SEC comment
                              # for the full writeup) -- reverted to 180/240s. Both are per-event TARGETS,
                              # not hard values -- see `_effective_half_sec`: G10's 2024-04-23a/b are
                              # only 357s apart, so their actual window is capped below these targets to
                              # avoid one event's profile/search reaching into the other's real anomaly.)
MIN_OFF_EVENT_SAMPLES = 200
MIN_POST_EVENT_SEC = 60      # a trace with less than this much post-origin data is flagged triggered_short_trace
PERSISTENCE_Z_THRESHOLD = 2.0
PERSISTENCE_TAIL_FRACTION = 0.3
ALIGNMENT_TOLERANCE_SEC = 5

OUT_DIR = common.PROJECT_DIR / "data" / "interim" / "seismometer_comparison"

_MSEED_CACHE: dict[str, object] = {}


# ---------------------------------------------------------------------------
# PZ catalog: parse each PZ file's own dashed-comment header (no separate
# filename-parsing path -- confirmed the header itself carries every field
# needed).
# ---------------------------------------------------------------------------

_PZ_FIELD_RE = {
    "station": re.compile(r"STATION\s*\(KSTNM\):\s*(\S+)"),
    "location": re.compile(r"LOCATION\s*\(KHOLE\):\s*(\S+)"),
    "channel": re.compile(r"CHANNEL\s*\(KCMPNM\):\s*(\S+)"),
    "start": re.compile(r"\bSTART\s*:\s*(\S+)"),
    "end": re.compile(r"\bEND\s*:\s*(\S+)"),
    "lat": re.compile(r"LATITUDE\s*:\s*([\-0-9.]+)"),
    "lon": re.compile(r"LONGITUDE\s*:\s*([\-0-9.]+)"),
    "elevation": re.compile(r"ELEVATION\s*:\s*([\-0-9.]+)"),
}


def _parse_pz_header(text: str) -> dict | None:
    vals = {}
    for key, pat in _PZ_FIELD_RE.items():
        m = pat.search(text)
        if m is None:
            return None
        vals[key] = m.group(1)
    return {
        "station": vals["station"], "location": vals["location"], "channel": vals["channel"],
        "start_utc": pd.Timestamp(vals["start"]), "end_utc": pd.Timestamp(vals["end"]),
        "lat": float(vals["lat"]), "lon": float(vals["lon"]), "elevation_m": float(vals["elevation"]),
    }


def _load_pz_catalog(pz_dir: Path) -> dict[tuple[str, str], list[dict]]:
    """key = (station, channel). value = list of epoch dicts (location,
    start_utc, end_utc, lat, lon, elevation_m, path), one per PZ file --
    a station/channel can have more than one epoch (differentiated by
    `location` and/or `start_utc`)."""
    catalog: dict[tuple[str, str], list[dict]] = {}
    for path in sorted(pz_dir.glob("SAC_PZs_*")):
        text = path.read_text(errors="replace")
        h = _parse_pz_header(text)
        if h is None:
            continue
        h["path"] = path
        catalog.setdefault((h["station"], h["channel"]), []).append(h)
    return catalog


def _select_pz_epoch(catalog: dict, station: str, channel: str, location: str,
                      event_utc: pd.Timestamp) -> dict | None:
    """Match on `location` (read directly off the miniSEED trace's own
    stats.location -- unambiguous per-trace metadata, not a guess), then
    pick the entry in force at event_utc (start_utc <= event_utc <= end_utc)
    with the latest start_utc. Falls back to ignoring the location match only
    if no in-force epoch matches it (shouldn't normally happen)."""
    entries = catalog.get((station, channel), [])
    # the epoch has to still be in force at the event, not just have started --
    # otherwise a station re-instrumented before the event would get its old response
    in_force = [e for e in entries if e["start_utc"] <= event_utc <= e["end_utc"]]
    candidates = [e for e in in_force if e["location"] == location]
    if not candidates:
        candidates = in_force
    if not candidates:
        return None
    return max(candidates, key=lambda e: e["start_utc"])


def find_colocated_seismic_stations(pz_dir: Path, ref_lat: float, ref_lon: float, top_n: int = 3) -> list[dict]:
    """Ranks seismic stations by distance to `ref_lat`/`ref_lon` -- the
    GEOMAGNETIC station's own coordinates, not the epicenter. The goal is
    co-location with the magnetometer being tested, not proximity to the
    earthquake source."""
    catalog = _load_pz_catalog(pz_dir)
    seen: dict[str, float] = {}
    for (station, _channel), entries in catalog.items():
        if station in seen:
            continue
        e = entries[0]
        seen[station] = common.haversine_km(ref_lat, ref_lon, e["lat"], e["lon"])
    ranked = sorted(seen.items(), key=lambda kv: kv[1])
    return [{"station": s, "distance_km": round(d, 2)} for s, d in ranked[:top_n]]


# ---------------------------------------------------------------------------
# miniSEED access
# ---------------------------------------------------------------------------

def _read_mseed_cached(mseed_path: Path):
    from obspy import read
    key = str(mseed_path)
    if key not in _MSEED_CACHE:
        _MSEED_CACHE[key] = read(str(mseed_path))
    return _MSEED_CACHE[key]


def load_trace(mseed_path: Path, pz_catalog: dict, station: str, channel: str,
                event_utc: pd.Timestamp, remove_response: bool = True) -> dict:
    """Reads (cached) the whole mseed file, selects the (station, channel)
    traces, merges any duplicates (confirmed present -- some (station,
    channel) pairs had 2 traces in the same file, a gap-split recording or
    two nearby triggers), removes instrument response via the matching PZ
    epoch if requested, and flags short triggered-accelerograph traces
    rather than assuming a uniform window."""
    from obspy import UTCDateTime

    st = _read_mseed_cached(mseed_path)
    sel = st.select(station=station, channel=channel).copy()
    if len(sel) == 0:
        return {"status": "no_trace"}
    try:
        sel.merge(method=1, fill_value=None)
    except Exception:
        pass

    ev = UTCDateTime(event_utc.isoformat())
    overlapping = [tr for tr in sel if tr.stats.starttime <= ev <= tr.stats.endtime]
    tr = (min(overlapping, key=lambda t: t.stats.endtime - t.stats.starttime)
          if overlapping else min(sel, key=lambda t: abs(t.stats.starttime - ev)))
    tr = tr.copy()

    triggered_short = (tr.stats.endtime - ev) < MIN_POST_EVENT_SEC
    response_removed = False
    if remove_response:
        from obspy.io.sac.sacpz import attach_paz
        epoch = _select_pz_epoch(pz_catalog, station, channel, tr.stats.location, event_utc)
        if epoch is not None:
            try:
                attach_paz(tr, str(epoch["path"]))
                tr.simulate(paz_remove=tr.stats.paz)
                response_removed = True
            except Exception as exc:
                return {"status": "response_removal_failed", "error": str(exc)}
        else:
            return {"status": "no_pz_epoch"}

    return {"status": "ok", "trace": tr, "response_removed": response_removed,
            "triggered_short_trace": bool(triggered_short)}


# ---------------------------------------------------------------------------
# Strong-motion window detectors (both computed, neither silently preferred
# -- see module docstring for why: the confirmed short-triggered-trace case
# doesn't have enough pre-event baseline for STA/LTA).
# ---------------------------------------------------------------------------

def detect_window_sta_lta(tr, event_utc: pd.Timestamp, sta_sec: float = 1.0, lta_sec: float = 10.0,
                           trigger_on: float = 3.5, trigger_off: float = 1.0) -> dict:
    from obspy import UTCDateTime
    from obspy.signal.trigger import classic_sta_lta, trigger_onset

    sr = tr.stats.sampling_rate
    ev = UTCDateTime(event_utc.isoformat())
    pre_event_sec = ev - tr.stats.starttime
    if pre_event_sec < lta_sec + 2:
        return {"status": "no_pre_event_baseline", "pre_event_sec": round(float(pre_event_sec), 2)}

    cft = classic_sta_lta(tr.data.astype(float), max(1, int(sta_sec * sr)), max(2, int(lta_sec * sr)))
    onsets = trigger_onset(cft, trigger_on, trigger_off)
    if len(onsets) == 0:
        return {"status": "no_trigger"}

    ev_idx = int(round((ev - tr.stats.starttime) * sr))
    on_i, off_i = min(onsets, key=lambda w: abs(w[0] - ev_idx))
    return {
        "status": "ok",
        "onset_lag_sec": round(float(on_i / sr - (ev - tr.stats.starttime)), 2),
        "offset_lag_sec": round(float(off_i / sr - (ev - tr.stats.starttime)), 2),
    }


def detect_window_envelope_threshold(tr, event_utc: pd.Timestamp, threshold_fraction: float = 0.1,
                                      band: tuple[float, float] = (1.0, 20.0)) -> dict:
    from obspy import UTCDateTime
    from obspy.signal.filter import envelope

    ev = UTCDateTime(event_utc.isoformat())
    tr2 = tr.copy()
    nyquist = tr.stats.sampling_rate / 2 - 0.5
    try:
        tr2.filter("bandpass", freqmin=band[0], freqmax=min(band[1], nyquist))
    except Exception:
        pass
    env = envelope(tr2.data.astype(float))
    peak = float(np.max(env)) if len(env) else 0.0
    if peak <= 0:
        return {"status": "flat_trace"}
    above = np.where(env > threshold_fraction * peak)[0]
    if len(above) == 0:
        return {"status": "no_signal_above_threshold"}

    sr = tr.stats.sampling_rate
    start_offset = tr.stats.starttime - ev  # seconds, trace start relative to origin
    return {
        "status": "ok",
        "onset_lag_sec": round(float(above[0] / sr + start_offset), 2),
        "offset_lag_sec": round(float(above[-1] / sr + start_offset), 2),
        "peak_amplitude": peak,
    }


def _effective_half_sec(target_half_sec: int, event_utc: pd.Timestamp,
                         exclude_centers: list[pd.Timestamp]) -> int:
    """Local copy of coseismic_step_analysis.py's identically-named helper
    (see that module for the full rationale) -- caps a per-event window at
    half the gap to the nearest *other* real event in the same group, so
    G10's 2024-04-23a/b (357s apart) never let one event's extracted
    profile/search reach into the other's real anomaly. A no-op everywhere
    else (every other group's events are hours-to-years apart)."""
    others = [c for c in exclude_centers if c != event_utc]
    if not others:
        return target_half_sec
    nearest_gap_sec = min(abs((event_utc - c).total_seconds()) for c in others)
    return int(min(target_half_sec, nearest_gap_sec // 2))


# ---------------------------------------------------------------------------
# Geomagnetic side: reuses coseismic_step_analysis.py's detrend/statistic
# machinery directly (independent computation, no dependency on
# coseismic_stacking_analysis.py's cached output).
# ---------------------------------------------------------------------------

def geomag_profile(cfg: "common.GroupConfig", group, event, station: str, channel_type: str,
                    half_sec: int = GEOMAG_HALF_SEC) -> dict | None:
    event_utc = pd.Timestamp(event.time_utc)
    sibling_utcs = [pd.Timestamp(e.time_utc) for e in folder_events(group.group_id)]
    half_sec = _effective_half_sec(half_sec, event_utc, sibling_utcs)
    effective_search_half_sec = min(SEARCH_HALF_SEC, half_sec)
    buffer_sec = half_sec + EXCLUSION_BUFFER_SEC + DETREND_WINDOW_SEC // 2 + 60
    df, missing, wanted = _load_station_days(cfg.gdms_dir, station, event_utc, buffer_sec)
    if df is None:
        return None
    channels = _build_channels(df)
    ch_label = "H" if channel_type == "XYZ" else "F"
    if ch_label not in channels:
        return None
    raw = channels[ch_label]
    idx = raw.index
    if idx.min() + pd.Timedelta(seconds=half_sec) > event_utc or \
       idx.max() - pd.Timedelta(seconds=half_sec) < event_utc:
        return None

    d = _detrend(raw)
    arr = _step_statistic(d, 30)  # step30: same "headline" statistic as coseismic_stacking_analysis.py, used for ONSET timing

    exclude_centers = sibling_utcs  # same list computed above for the half_sec cap; kept as its own name here since
                                      # its role from this point on is null/baseline exclusion, not window sizing

    def _off_event_baseline(series: np.ndarray) -> tuple[float, float] | None:
        mask = np.ones(len(series), dtype=bool)
        for c in exclude_centers:
            p = _pos(idx, c)
            lo, hi = max(0, p - EXCLUSION_BUFFER_SEC), min(len(series), p + EXCLUSION_BUFFER_SEC + 1)
            mask[lo:hi] = False
        vals = series[mask]
        vals = vals[~np.isnan(vals)]
        if len(vals) < MIN_OFF_EVENT_SAMPLES:
            return None
        _, med, mad = mad_zscore(vals)
        if not np.isfinite(mad) or mad <= 1e-9:
            return None
        return med, mad

    def _extract_profile(series: np.ndarray) -> np.ndarray:
        p_center = _pos(idx, event_utc)
        lo_i, hi_i = max(0, p_center - half_sec), min(len(series), p_center + half_sec + 1)
        profile = np.full(2 * half_sec + 1, np.nan)
        dst_lo = lo_i - (p_center - half_sec)
        profile[dst_lo:dst_lo + (hi_i - lo_i)] = series[lo_i:hi_i]
        return profile

    step_base = _off_event_baseline(arr)
    if step_base is None:
        return None
    step_med, step_mad = step_base
    z_profile = (_extract_profile(arr) - step_med) / step_mad
    lags = np.arange(-half_sec, half_sec + 1)

    # `arr` (a moving-window step DIFFERENCE) is well-suited to finding WHEN
    # a transition happens, but is structurally the wrong quantity for
    # asking whether the field STAYS shifted afterward: once both sides of
    # the sliding window sit on the same new plateau, a step-difference
    # statistic returns to ~0 by construction, regardless of whether that
    # plateau is permanent (a real persistent offset) or itself about to
    # revert. Persistence has to be judged on the detrended FIELD LEVEL `d`
    # itself, z-scored the same way -- a separate profile from the one used
    # for onset timing.
    level_base = _off_event_baseline(d)
    level_z_profile = (_extract_profile(d) - level_base[0]) / level_base[1] if level_base is not None else None

    search_mask = np.abs(lags) <= effective_search_half_sec
    sub, sub_lags = z_profile[search_mask], lags[search_mask]
    if np.all(np.isnan(sub)):
        obs_lag_sec, obs_peak_z = None, None
    else:
        j = int(np.nanargmax(np.abs(sub)))
        obs_lag_sec, obs_peak_z = int(sub_lags[j]), round(float(sub[j]), 4)

    return {
        "station": station, "statistic": "step30", "half_sec": half_sec,
        "search_half_sec": effective_search_half_sec,
        "lags_sec": lags.tolist(),
        "z_profile": [None if np.isnan(v) else round(float(v), 4) for v in z_profile],
        "level_z_profile": (None if level_z_profile is None else
                             [None if np.isnan(v) else round(float(v), 4) for v in level_z_profile]),
        "obs_lag_sec": obs_lag_sec, "obs_peak_z": obs_peak_z,
    }


# ---------------------------------------------------------------------------
# Comparison / verdict
# ---------------------------------------------------------------------------

def assess_persistence(lags_sec: list[int], z_profile: list[float | None],
                        shaking_offset_lag_sec: float | None,
                        z_threshold: float = PERSISTENCE_Z_THRESHOLD,
                        tail_frac: float = PERSISTENCE_TAIL_FRACTION) -> bool | None:
    """True if the geomagnetic z-profile stays elevated (|z| > z_threshold)
    for at least `tail_frac` of its samples after the shaking has stopped
    (or after lag=0 if no shaking-offset is known) -- the quantitative form
    of "G9 csg doesn't revert, G10 xcg does"."""
    lags = np.array(lags_sec, dtype=float)
    z = np.array([np.nan if v is None else v for v in z_profile], dtype=float)
    ref = shaking_offset_lag_sec if shaking_offset_lag_sec is not None else 0.0
    tail = z[lags > ref]
    tail = tail[~np.isnan(tail)]
    if len(tail) == 0:
        return None
    return bool(np.mean(np.abs(tail) > z_threshold) >= tail_frac)


def alignment_verdict(geomag_lag_sec: float | None, shaking_onset_lag_sec: float | None,
                       shaking_offset_lag_sec: float | None, persists: bool | None,
                       tolerance_sec: float = ALIGNMENT_TOLERANCE_SEC) -> str:
    """One of: aligned_with_shaking / leads_shaking / persists_after_shaking_ends
    / insufficient_data -- directly implements the user's own framing:
    alignment supports the instrument-noise explanation; leading or
    outlasting the shaking supports a real geophysical mechanism."""
    if geomag_lag_sec is None or shaking_onset_lag_sec is None:
        return "insufficient_data"
    if persists:
        return "persists_after_shaking_ends"
    if geomag_lag_sec < shaking_onset_lag_sec - tolerance_sec:
        return "leads_shaking"
    return "aligned_with_shaking"


def compare_event(cfg: "common.GroupConfig", group, event, geomag_station: str, channel_type: str,
                   seismic_station: str, mseed_path: Path, pz_catalog: dict) -> dict:
    event_utc = pd.Timestamp(event.time_utc)
    geomag = geomag_profile(cfg, group, event, geomag_station, channel_type)
    if geomag is None:
        return {"status": "no_geomag_data"}

    seismic_channels: dict[str, dict] = {}
    best_comp = None
    for comp in ("HLZ", "HLN", "HLE"):
        info = load_trace(mseed_path, pz_catalog, seismic_station, comp, event_utc)
        if info["status"] != "ok":
            seismic_channels[comp] = {"status": info["status"]}
            continue
        sta_lta = detect_window_sta_lta(info["trace"], event_utc)
        env = detect_window_envelope_threshold(info["trace"], event_utc)
        seismic_channels[comp] = {
            "status": "ok",
            "response_removed": info["response_removed"],
            "triggered_short_trace": info["triggered_short_trace"],
            "trace_start_utc": str(info["trace"].stats.starttime),
            "trace_end_utc": str(info["trace"].stats.endtime),
            "sta_lta": sta_lta, "envelope": env,
        }
        if best_comp is None or comp == "HLZ":
            best_comp = comp

    if best_comp is None:
        return {"status": "no_seismic_trace", "geomag": geomag, "seismic_channels": seismic_channels}

    primary = seismic_channels[best_comp]
    shaking_onset = (primary["sta_lta"]["onset_lag_sec"] if primary["sta_lta"]["status"] == "ok"
                      else primary["envelope"].get("onset_lag_sec") if primary["envelope"]["status"] == "ok"
                      else None)
    shaking_offset = (primary["sta_lta"]["offset_lag_sec"] if primary["sta_lta"]["status"] == "ok"
                       else primary["envelope"].get("offset_lag_sec") if primary["envelope"]["status"] == "ok"
                       else None)

    persists = (assess_persistence(geomag["lags_sec"], geomag["level_z_profile"], shaking_offset)
                if geomag.get("level_z_profile") is not None else None)
    verdict = alignment_verdict(geomag["obs_lag_sec"], shaking_onset, shaking_offset, persists)

    return {
        "status": "ok",
        "primary_seismic_channel": best_comp,
        "geomag": geomag,
        "seismic_channels": seismic_channels,
        "shaking_onset_lag_sec": shaking_onset,
        "shaking_offset_lag_sec": shaking_offset,
        "geomag_persists_after_shaking": persists,
        "alignment_verdict": verdict,
    }


# ---------------------------------------------------------------------------
# Case studies: compare_event's result plus a downsampled raw waveform, for
# report plotting. Reuses compare_event rather than recomputing onset
# detection separately.
# ---------------------------------------------------------------------------

def _downsample(data: np.ndarray, factor: int) -> list[float]:
    if factor <= 1:
        return [round(float(v), 4) for v in data]
    n = (len(data) // factor) * factor
    if n == 0:
        return []
    return [round(float(v), 4) for v in data[:n].reshape(-1, factor).mean(axis=1)]


def build_case_study(group_id: str, event_date: str, geomag_station: str) -> dict:
    cfg = common.load_group_config(group_id)
    group = get_group(group_id)
    event = next((e for e in group.events if e.date == event_date), None)
    if event is None:
        return {"status": "unknown_event", "group": group_id, "event_date": event_date}

    if group_id not in SEISMIC_DATA_DIRS or event_date not in SEISMIC_DATA_DIRS[group_id]:
        return {"status": "no_seismic_data", "group": group_id, "event_date": event_date}

    channel_type = "XYZ" if cfg.xyz_pool.all_stations else "F"
    geomag_lat, geomag_lon = cfg.stations[geomag_station]["lat"], cfg.stations[geomag_station]["lon"]
    pz_dir = SEISMIC_ROOT / SEISMIC_DATA_DIRS[group_id][event_date]["pz_dir"]
    mseed_path = pz_dir / SEISMIC_DATA_DIRS[group_id][event_date]["mseed"]
    pz_catalog = _load_pz_catalog(pz_dir)

    colocated = find_colocated_seismic_stations(pz_dir, geomag_lat, geomag_lon, top_n=1)
    if not colocated:
        return {"status": "no_seismic_station", "group": group_id, "event_date": event_date}
    seismic_station = colocated[0]["station"]

    result = compare_event(cfg, group, event, geomag_station, channel_type, seismic_station, mseed_path, pz_catalog)
    result["group"] = group_id
    result["event_date"] = event_date
    result["geomag_station"] = geomag_station
    result["seismic_station"] = seismic_station
    result["seismic_colocation_km"] = colocated[0]["distance_km"]

    if result.get("status") != "ok":
        return result

    event_utc = pd.Timestamp(event.time_utc)
    primary_comp = result["primary_seismic_channel"]
    info = load_trace(mseed_path, pz_catalog, seismic_station, primary_comp, event_utc)
    waveform_payload = None
    if info["status"] == "ok":
        from obspy.signal.filter import envelope as _envelope
        tr = info["trace"]
        sr = tr.stats.sampling_rate
        t0_offset = (tr.stats.starttime.datetime - event_utc.to_pydatetime()).total_seconds()
        factor = max(1, int(round(sr / 20)))  # downsample toward ~20Hz for compact JSON
        wave_ds = _downsample(tr.data.astype(float), factor)
        env_ds = _downsample(_envelope(tr.data.astype(float)), factor)
        lags_ds = [round(t0_offset + i * factor / sr, 3) for i in range(len(wave_ds))]
        waveform_payload = {"sampling_rate_downsampled_hz": round(sr / factor, 3),
                             "lags_sec": lags_ds, "waveform": wave_ds, "envelope": env_ds}
    result["seismic_waveform_downsampled"] = waveform_payload
    return result


# ---------------------------------------------------------------------------
# Self-test: PZ-parse + response-removal round-trip on a known real file
# (not synthetic -- the thing worth sanity-checking here is obspy/PZ
# plumbing, not the statistics, which coseismic_stacking_analysis.py's
# self-test already covers).
# ---------------------------------------------------------------------------

def self_test() -> bool:
    ok = True
    pz_dir = SEISMIC_ROOT / "G10_0403"
    catalog = _load_pz_catalog(pz_dir)
    key = ("ALS", "HLZ")
    status1 = "PASS" if catalog.get(key) else "FAIL"
    print(f"[self-test] PZ catalog parse ALS/HLZ: {len(catalog.get(key, []))} epoch(s)  {status1}")
    ok = ok and status1 == "PASS"

    event_utc = pd.Timestamp("2024-04-02 23:58:11")
    mseed_path = pz_dir / "G10_0403_w.mseed"
    info = load_trace(mseed_path, catalog, "ALS", "HLZ", event_utc, remove_response=True)
    status2 = "PASS" if info.get("status") == "ok" and info.get("response_removed") else "FAIL"
    print(f"[self-test] load_trace ALS/HLZ status={info.get('status')} "
          f"response_removed={info.get('response_removed')}  {status2}")
    ok = ok and status2 == "PASS"

    if info.get("status") == "ok":
        data = info["trace"].data
        finite = bool(np.all(np.isfinite(data)))
        spread = float(np.std(data))
        status3 = "PASS" if finite and spread > 0 else "FAIL"
        print(f"[self-test] response-removed trace sanity: finite={finite} std={spread:.6g}  {status3}")
        ok = ok and status3 == "PASS"

        sta_lta = detect_window_sta_lta(info["trace"], event_utc)
        env = detect_window_envelope_threshold(info["trace"], event_utc)
        status4 = ("PASS" if sta_lta["status"] in ("ok", "no_trigger", "no_pre_event_baseline")
                   and env["status"] in ("ok", "no_signal_above_threshold", "flat_trace") else "FAIL")
        print(f"[self-test] detectors ran without error: sta_lta={sta_lta['status']} "
              f"envelope={env['status']}  {status4}")
        ok = ok and status4 == "PASS"

    colocated = find_colocated_seismic_stations(pz_dir, 24.038, 121.609, top_n=3)  # xcg's coordinates
    status5 = "PASS" if colocated and colocated[0]["station"] == "HWA" else "FAIL"
    print(f"[self-test] find_colocated_seismic_stations(xcg): nearest={colocated[:1]}  {status5}")
    ok = ok and status5 == "PASS"

    return ok


# ---------------------------------------------------------------------------
# Coverage summary + real-data orchestration
# ---------------------------------------------------------------------------

def build_coverage_summary() -> dict:
    items = []
    for group_id, group in GROUPS.items():
        group_dirs = SEISMIC_DATA_DIRS.get(group_id, {})
        for event in group.events:
            if event.date in group_dirs:
                status = "available"
            elif pd.Timestamp(event.time_utc) < SEISMIC_DATA_SOURCE_START_UTC:
                status = "no_data_pre_2012"
            else:
                status = "not_fetched"
            items.append({"group": group_id, "date": event.date, "anchor": event.anchor,
                           "seismic_data_status": status})
    n_available = sum(1 for it in items if it["seismic_data_status"] == "available")
    return {"n_events_total": len(items), "n_available": n_available, "events": items}


def run_available_events(group_ids: tuple[str, ...] | None = None) -> dict:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "events").mkdir(exist_ok=True)
    (OUT_DIR / "case_studies").mkdir(exist_ok=True)

    coverage = build_coverage_summary()
    (OUT_DIR / "coverage_summary.json").write_text(json.dumps(coverage, indent=2))

    targets = group_ids if group_ids else tuple(SEISMIC_DATA_DIRS.keys())
    rows: list[dict] = []
    run_summary = {"events": []}

    for group_id in targets:
        group_dirs = SEISMIC_DATA_DIRS.get(group_id, {})
        if not group_dirs:
            print(f"[{group_id}] no seismic data fetched, skipping", file=sys.stderr)
            continue
        cfg = common.load_group_config(group_id)
        group = get_group(group_id)

        for event in group.events:
            if event.date not in group_dirs:
                continue  # e.g. G14's events (no seismic source coverage at all)
            pz_dir = SEISMIC_ROOT / group_dirs[event.date]["pz_dir"]
            mseed_path = pz_dir / group_dirs[event.date]["mseed"]
            pz_catalog = _load_pz_catalog(pz_dir)

            channel_type = "XYZ" if cfg.xyz_pool.all_stations else "F"
            near = _rank_stations_for_event(cfg, event, channel_type, 1)
            if not near:
                run_summary["events"].append({"group": group_id, "date": event.date, "status": "no_geomag_station"})
                continue
            geomag_station, geomag_distance_km = near[0]
            geomag_lat, geomag_lon = cfg.stations[geomag_station]["lat"], cfg.stations[geomag_station]["lon"]

            colocated = find_colocated_seismic_stations(pz_dir, geomag_lat, geomag_lon, top_n=3)
            if not colocated:
                run_summary["events"].append({"group": group_id, "date": event.date, "status": "no_seismic_station"})
                continue
            seismic_station = colocated[0]["station"]

            result = compare_event(cfg, group, event, geomag_station, channel_type,
                                    seismic_station, mseed_path, pz_catalog)
            result.update({"group": group_id, "geomag_station": geomag_station,
                            "geomag_distance_km": round(geomag_distance_km, 1),
                            "seismic_station": seismic_station,
                            "seismic_colocation_km": colocated[0]["distance_km"]})

            out_path = OUT_DIR / "events" / f"{group_id}__{event.date}__{geomag_station}.json"
            out_path.write_text(json.dumps(result, indent=2, default=str))

            rows.append({
                "group": group_id, "date": event.date, "anchor": event.anchor,
                "geomag_station": geomag_station,
                "geomag_distance_km": round(geomag_distance_km, 1),
                "seismic_station": seismic_station, "seismic_colocation_km": colocated[0]["distance_km"],
                "status": result.get("status"),
                "geomag_obs_lag_sec": (result.get("geomag") or {}).get("obs_lag_sec"),
                "shaking_onset_lag_sec": result.get("shaking_onset_lag_sec"),
                "shaking_offset_lag_sec": result.get("shaking_offset_lag_sec"),
                "geomag_persists_after_shaking": result.get("geomag_persists_after_shaking"),
                "alignment_verdict": result.get("alignment_verdict"),
            })
            run_summary["events"].append({"group": group_id, "date": event.date, "status": result.get("status"),
                                           "alignment_verdict": result.get("alignment_verdict")})
            print(f"[{group_id} {event.date}] geomag={geomag_station} seismic={seismic_station} "
                  f"status={result.get('status')} verdict={result.get('alignment_verdict')}", file=sys.stderr)

    pd.DataFrame(rows).to_csv(OUT_DIR / "comparison_summary.csv", index=False)
    (OUT_DIR / "all_comparisons_run_summary.json").write_text(json.dumps(run_summary, indent=2))
    print(f"[run] {len(rows)} events compared -> {OUT_DIR}", file=sys.stderr)

    if "G9" in targets:
        cs = build_case_study("G9", "2022-09-18", "csg")
        (OUT_DIR / "case_studies" / "G9_csg_vs_ECS.json").write_text(json.dumps(cs, indent=2, default=str))
        print(f"[case-study] G9 csg -> verdict={cs.get('alignment_verdict')}", file=sys.stderr)
    if "G10" in targets:
        cs = build_case_study("G10", "2024-04-03", "xcg")
        (OUT_DIR / "case_studies" / "G10_xcg_vs_HWA.json").write_text(json.dumps(cs, indent=2, default=str))
        print(f"[case-study] G10 xcg -> verdict={cs.get('alignment_verdict')}", file=sys.stderr)

    return run_summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--self-test", action="store_true", help="PZ-parse + response-removal sanity check only")
    ap.add_argument("--group", action="append", dest="groups",
                     help="restrict to this group ID (repeatable); default is every group with seismic data fetched")
    ap.add_argument("--geomag-station", dest="geomag_station", default=None,
                     help="(informational, single-group spot checks) which geomag station to compare against")
    ap.add_argument("--all", action="store_true", help="explicit alias for the default (no --group filter)")
    args = ap.parse_args()

    if args.self_test:
        sys.exit(0 if self_test() else 1)

    if not self_test():
        print("[main] self-test FAILED -- aborting before touching real data", file=sys.stderr)
        sys.exit(1)

    group_ids = tuple(args.groups) if args.groups else None
    run_available_events(group_ids)
