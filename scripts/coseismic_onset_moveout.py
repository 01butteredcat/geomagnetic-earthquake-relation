"""Onset moveout test: does the coseismic geomagnetic anomaly start later at
stations farther from the hypocenter?

This test needs neither seismometer data nor noise/signal labels. Three
candidate sources predict different onset patterns across the network:

  - external disturbance (SSC, substorm onset, Pi2): simultaneous everywhere,
    onset independent of distance (slope ~ 0 s/km);
  - magnetometer shaken by ground motion: onset follows the seismic waves,
    onset ~ R / V with V ~ 6 km/s (P) or ~ 3.5 km/s (S), i.e. a slope of
    ~0.17-0.29 s/km. A 20 -> 200 km spread is ~30-50 s, well resolved at 1 Hz;
  - a lithospheric source at the hypocenter: earliest near the epicenter,
    possibly ahead of the P wave, but no reason to travel at seismic speed.

## Onset detector

On each station's raw 1 Hz first difference (H = sqrt(X^2+Y^2) for vector
stations, F for scalar ones -- same channels as coseismic_step_analysis.py's
`_build_channels`):

  sigma  = 1.4826 * MAD of the first difference over NOISE_PRE_LAGS
           (-660..-60 s, same pre-event reference as
           seismometer_comparison.py::geomag_noise_ratio)
  exceed = |diff| > K_SIGMA * sigma
  onset  = first second t in SEARCH_LAGS (-60..+180 s, capped per event by
           coseismic_step_analysis._effective_half_sec) where exceed[t] holds
           and at least MIN_EXCEED of the MIN_EXCEED_WINDOW seconds t..t+4 exceed.

This is the shaking-noise signature of the G10 pilot (first-difference noise
~11.6x the pre-event level), timed at its first sample rather than at the
step30 peak, which can sit anywhere inside the burst.

## Null calibration

The same detector is run at N_NULL random reference times in the same loaded
data (every real event in the raw-data folder kept NULL_EXCLUSION_SEC away, so
the whole -660..+180 s window stays clear). Since 2026-09-28 the reference
times are **shared by every station of an event** (one draw per event, kept
only where every station tested at the origin also gets a valid result there;
with fewer than MIN_SHARED_NULL such times the event falls back to counting
invalid draws as "no trigger" and is flagged `shared_null_fallback`). Before
that each station drew its own times, which hides the fact that an external
disturbance hits the whole network at once. The fraction that trigger is the
station's false-trigger rate for this detector; `trigger_p` = (1 + #null
triggers at least as early) / (N_NULL + 1) is not used for anything below, but
lets a trigger be read against its own station's noise.

## Moveout statistics

Hypocentral distance R = sqrt(epicentral^2 + depth^2).

  - Per event (>= MIN_STATIONS_FOR_FIT triggered stations): Theil-Sen slope of
    onset lag vs. R, with its 90 % CI, and the apparent velocity 1/slope.
    `simultaneous` if the CI contains 0 and lies below the P slowness;
    `seismic_moveout` if the CI is above 0 and overlaps [1/V_P_KMS, 1/V_S_KMS];
    otherwise `indeterminate`.
  - Pooled over events: Theil-Sen slope of (onset lag vs. R) after removing
    each event's own median lag and median R (so events with different origin
    offsets don't masquerade as moveout), p-value from shuffling onsets among
    stations of the *same event* (one-sided, H1: slope > 0). Reported
    separately for M >= 6 and M < 6 (the 68 backfilled M5 events), and for
    all triggers vs. only stations whose own false-trigger rate is below
    CLEAN_FALSE_RATE.
  - Arrival-window test (primary, event-level null since 2026-09-28): T = the
    number of stations whose onset lands in [R/V_P - 5 s, R/V_S + 30 s]. The
    null evaluates the whole event at one shared random reference time: each
    simulation draws one shared index per event and counts that event's
    stations whose null trigger lands in their own window, summed over events.
    This keeps the correlation between stations of the same event (a substorm
    triggers them all at once), which the older per-trigger Poisson-binomial
    (observed count vs. the sum of each trigger's null in-window fraction, kept
    as `arrival_window_all_triggers` for reference) treats as independent.
    This is the direct "tied to the seismic waves" test; the slope fits are
    easily dragged by a single noisy station.
  - Cross-check against seismometer_comparison.py's co-located seismic onset
    (`shaking_onset_lag_sec`): geomagnetic onset minus seismic onset, for the
    one station per event that script compared.

## Prior-event shaking (sensitivity only)

An event that follows another registered event in the same raw-data folder by
<= PRIOR_EVENT_EXCLUSION_SEC (900 s: the -660 s noise reference plus a few
minutes of the earlier event's shaking) has its noise reference and search
window inside that earlier shaking, so its onsets are not timed against a
quiet baseline. Such events get `prior_event_shaking` = True. The main subsets
(`m6`, `m5`) keep them -- the main test was fixed before this was noticed --
and `m6_no_prior_shaking` / `m5_no_prior_shaking` repeat every statistic
without them, each on its own keyed rng stream so the main numbers do not
move. Only registered events are checked: M5 events are registered only from
2024-09 on and the GDMS json export is M>=6 only, so smaller aftershocks
(e.g. inside G10 2024-04-03c, 2 h after the M7.2) cannot be ruled out.

Usage:
  coseismic_onset_moveout.py --self-test
  coseismic_onset_moveout.py --all [--min-mag 6] [--group G10 ...]
"""
from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm, theilslopes

sys.path.insert(0, str(Path(__file__).parent))
import common  # noqa: E402
from coseismic_step_analysis import (SEED, _build_channels, _effective_half_sec,  # noqa: E402
                                     _load_station_days, keyed_rng)
from events import GROUPS, folder_events  # noqa: E402
from stat_utils import MAD_SCALE  # noqa: E402

NOISE_PRE_LAGS = (-660, -60)
SEARCH_LAGS = (-60, 180)
K_SIGMA = 4.0
MIN_EXCEED_WINDOW = 5
MIN_EXCEED = 3
MIN_PRE_SAMPLES = 300
MAX_SEARCH_MISSING = 0.2
LOAD_BUFFER_SEC = 6 * 3600
N_NULL = 200
MIN_SHARED_NULL = 50  # below this many all-valid shared null times, fall back (see docstring)
NULL_EXCLUSION_SEC = -NOISE_PRE_LAGS[0] + 600
V_P_KMS, V_S_KMS = 6.0, 3.5
MIN_STATIONS_FOR_FIT = 3
CI_ALPHA = 0.90
CLEAN_FALSE_RATE = 0.05
CROSSCHECK_FALSE_RATE = 0.2
PRIOR_EVENT_EXCLUSION_SEC = 900  # -660 s noise reference + the earlier event's shaking
ARRIVAL_PAD_SEC = (-5, 30)  # "arrival window" = [R/V_P - 5, R/V_S + 30] s after origin
N_SIM = 20000
N_PERM = 2000
N_WORKERS = 8

COMPARISON_CSV = common.PROJECT_DIR / "data" / "interim" / "seismometer_comparison" / "comparison_summary.csv"
OUT_DIR = common.PROJECT_DIR / "data" / "interim" / "coseismic_onset_moveout"


# ---------------------------------------------------------------------------
# Detector
# ---------------------------------------------------------------------------

def detect_onset(values: np.ndarray, p0: int, search_hi: int) -> dict:
    """values: gap-free 1 Hz grid (NaN = missing); p0: array position of the
    reference time; search window is SEARCH_LAGS[0]..min(SEARCH_LAGS[1], search_hi)."""
    lo_pre, hi_pre = p0 + NOISE_PRE_LAGS[0], p0 + NOISE_PRE_LAGS[1]
    lo_s, hi_s = p0 + SEARCH_LAGS[0], p0 + min(SEARCH_LAGS[1], search_hi)
    if lo_pre < 1 or hi_s + MIN_EXCEED_WINDOW >= len(values):
        return {"status": "insufficient_window"}
    diff = values[lo_pre:hi_s + MIN_EXCEED_WINDOW + 1] - values[lo_pre - 1:hi_s + MIN_EXCEED_WINDOW]
    pre = diff[:hi_pre - lo_pre]
    pre = pre[~np.isnan(pre)]
    if len(pre) < MIN_PRE_SAMPLES:
        return {"status": "insufficient_pre"}
    sigma = float(np.median(np.abs(pre - np.median(pre)))) * MAD_SCALE
    if sigma <= 1e-9:
        return {"status": "flat_pre"}
    tail = diff[lo_s - lo_pre:]
    search = tail[:hi_s - lo_s + 1]
    if np.mean(np.isnan(search)) > MAX_SEARCH_MISSING:
        return {"status": "gappy_search", "sigma_nt": sigma}
    exceed = np.nan_to_num(np.abs(tail), nan=0.0) > K_SIGMA * sigma
    counts = np.convolve(exceed.astype(int), np.ones(MIN_EXCEED_WINDOW, dtype=int), "full")[
        MIN_EXCEED_WINDOW - 1:]  # counts[i] = exceed[i..i+4]
    hits = np.flatnonzero(exceed[:len(search)] & (counts[:len(search)] >= MIN_EXCEED))
    peak = float(np.nanmax(np.abs(search))) / sigma if not np.all(np.isnan(search)) else None
    if len(hits) == 0:
        return {"status": "no_onset", "sigma_nt": sigma, "peak_ratio": peak}
    i = int(hits[0])
    burst = np.abs(search[i:i + 60])
    return {"status": "onset", "sigma_nt": sigma, "peak_ratio": peak,
            "onset_lag_sec": int(i + SEARCH_LAGS[0]),
            "burst_peak_ratio": float(np.nanmax(burst)) / sigma}


def arrival_window(hypocentral_km: float) -> tuple[float, float]:
    return hypocentral_km / V_P_KMS + ARRIVAL_PAD_SEC[0], hypocentral_km / V_S_KMS + ARRIVAL_PAD_SEC[1]


def prior_event_gap(event_utc: pd.Timestamp, others: list[pd.Timestamp]) -> float | None:
    """Seconds since the nearest earlier registered event (None if there is none)."""
    gaps = [(event_utc - t).total_seconds() for t in others]
    gaps = [g for g in gaps if g > 0]
    return min(gaps) if gaps else None


def _grid(series: pd.Series) -> pd.Series:
    return series.reindex(pd.date_range(series.index[0], series.index[-1], freq="s"))


def shared_null_outcomes(stations: dict[str, tuple[np.ndarray, pd.Timestamp]], event_utc: pd.Timestamp,
                         exclude: list[pd.Timestamp], rng: np.random.Generator,
                         search_hi: int) -> tuple[dict[str, list[int | None]], bool]:
    """Detector outcome of every station at the same N_NULL random reference
    times (None = no trigger). stations: name -> (gridded values, first
    timestamp). A time is kept only if every station gets a valid result there;
    with fewer than MIN_SHARED_NULL such times, the first N_NULL candidates are
    used with invalid results counted as no trigger (fallback = True)."""
    lo = -LOAD_BUFFER_SEC - NOISE_PRE_LAGS[0] + 1
    hi = LOAD_BUFFER_SEC - SEARCH_LAGS[1] - MIN_EXCEED_WINDOW - 2
    ex = np.array([(e - event_utc).total_seconds() for e in exclude])
    kept: dict[str, list] = {st: [] for st in stations}
    loose: dict[str, list] = {st: [] for st in stations}
    n_kept = n_loose = attempts = 0
    while n_kept < N_NULL and attempts < N_NULL * 5:
        attempts += 1
        off = int(rng.integers(lo, hi))
        if len(ex) and np.min(np.abs(ex - off)) < NULL_EXCLUSION_SEC:
            continue
        t = event_utc + pd.Timedelta(seconds=off)
        res = {}
        for st, (values, t0) in stations.items():
            r = detect_onset(values, int(round((t - t0).total_seconds())), search_hi)
            res[st] = r["onset_lag_sec"] if r["status"] == "onset" else (None if r["status"] == "no_onset" else "invalid")
        if n_loose < N_NULL:
            n_loose += 1
            for st, v in res.items():
                loose[st].append(None if v == "invalid" else v)
        if all(v != "invalid" for v in res.values()):
            n_kept += 1
            for st, v in res.items():
                kept[st].append(v)
    return (kept, False) if n_kept >= MIN_SHARED_NULL else (loose, True)


# ---------------------------------------------------------------------------
# Per event / per group
# ---------------------------------------------------------------------------

def process_group(group_id: str, min_mag: float | None = None) -> tuple[list[dict], dict]:
    """Rows (one per event x station) plus, per event, each tested station's
    null in-window indicator at the shared reference times (for the
    event-level arrival-window test)."""
    cfg = common.load_group_config(group_id)
    group = GROUPS[group_id]
    exclude = [pd.Timestamp(e.time_utc) for e in folder_events(group_id)]
    rows: list[dict] = []
    nulls: dict[tuple[str, str], dict[str, np.ndarray]] = {}
    for event in group.events:
        if min_mag is not None and event.magnitude < min_mag:
            continue
        event_utc = pd.Timestamp(event.time_utc)
        search_hi = _effective_half_sec(SEARCH_LAGS[1], event_utc, exclude)
        prior = prior_event_gap(event_utc, exclude)
        ev_rows: list[dict] = []
        loaded: dict[str, tuple[np.ndarray, pd.Timestamp]] = {}
        for station, meta in sorted(cfg.stations.items()):
            if meta["reported"] not in ("F", "XYZF"):
                continue
            epi = common.haversine_km(meta["lat"], meta["lon"], event.lat, event.lon)
            base = {"group": group_id, "event_date": event.date, "event_time_utc": event.time_utc,
                    "magnitude": event.magnitude, "anchor": event.anchor, "depth_km": event.depth_km,
                    "prior_event_sec": prior,
                    "prior_event_shaking": prior is not None and prior <= PRIOR_EVENT_EXCLUSION_SEC,
                    "station": station, "channel": "F" if meta["reported"] == "F" else "H",
                    "epicentral_km": round(epi, 1),
                    "hypocentral_km": round(float(np.hypot(epi, event.depth_km)), 1)}
            df, _, _ = _load_station_days(cfg.gdms_dir, station, event_utc, LOAD_BUFFER_SEC)
            if df is None:
                ev_rows.append({**base, "status": "no_data"})
                continue
            channels = _build_channels(df)
            if base["channel"] not in channels:  # header says F but the day file carries XYZ, or vice versa
                base["channel"] = next(iter(channels)) if len(channels) == 1 else "H"
            series = _grid(channels[base["channel"]])
            if not (series.index[0] <= event_utc <= series.index[-1]):
                ev_rows.append({**base, "status": "no_data"})
                continue
            values = series.to_numpy(dtype=float)
            p0 = int(round((event_utc - series.index[0]).total_seconds()))
            res = detect_onset(values, p0, search_hi)
            if res["status"] in ("onset", "no_onset"):
                loaded[station] = (values, series.index[0])
            ev_rows.append({**base, **res})

        outcomes, fallback = shared_null_outcomes(loaded, event_utc, exclude,
                                                  keyed_rng("onset", group_id, event.date, "shared"), search_hi) \
            if loaded else ({}, False)
        ev_null: dict[str, np.ndarray] = {}
        for row in ev_rows:
            null = outcomes.get(row["station"])
            if null is None:
                continue
            n_trig = sum(v is not None for v in null)
            lo_a, hi_a = arrival_window(row["hypocentral_km"])
            row.update({"n_null": len(null), "shared_null_fallback": fallback,
                        "null_false_rate": round(n_trig / len(null), 4) if null else None})
            ev_null[row["station"]] = np.array([v is not None and lo_a <= v <= hi_a for v in null])
            if row["status"] == "onset":
                null_on = [v for v in null if v is not None]
                row["in_arrival_window"] = bool(lo_a <= row["onset_lag_sec"] <= hi_a)
                # chance of landing in the same window for a trigger of this station's own noise;
                # uniform over the search window if the null never triggered
                row["null_q_arrival"] = round(float(np.mean([lo_a <= v <= hi_a for v in null_on])) if null_on else
                                              max(0.0, min(hi_a, search_hi) - max(lo_a, SEARCH_LAGS[0]) + 1)
                                              / (search_hi - SEARCH_LAGS[0] + 1), 4)
                if null:
                    row["trigger_p"] = round((1 + sum(v is not None and v <= row["onset_lag_sec"] for v in null))
                                             / (len(null) + 1), 4)
        for row in ev_rows:
            for k in ("sigma_nt", "peak_ratio", "burst_peak_ratio"):
                if row.get(k) is not None:
                    row[k] = round(row[k], 4)
        rows.extend(ev_rows)
        nulls[(group_id, event.date)] = ev_null
        print(f"[{group_id}] {event.date} M{event.magnitude}: "
              f"{sum(r.get('status') == 'onset' for r in ev_rows)} onsets"
              f"{' (shared-null fallback)' if fallback else ''}", file=sys.stderr)
    return rows, nulls


# ---------------------------------------------------------------------------
# Moveout statistics
# ---------------------------------------------------------------------------

def classify_slope(lo: float, hi: float) -> str:
    if lo <= 0 <= hi and hi < 1 / V_P_KMS:
        return "simultaneous"
    if lo > 0 and lo <= 1 / V_S_KMS and hi >= 1 / V_P_KMS:
        return "seismic_moveout"
    return "indeterminate"


def per_event_moveout(onsets: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (g, d), ev in onsets.groupby(["group", "event_date"], sort=False):
        trig = ev[ev.status == "onset"]
        row = {"group": g, "event_date": d, "magnitude": ev.magnitude.iloc[0],
               "prior_event_shaking": bool(ev.prior_event_shaking.iloc[0]),
               "n_stations_tested": int(ev.status.isin(["onset", "no_onset"]).sum()),
               "n_triggered": int(len(trig)),
               "nearest_triggered_km": None if trig.empty else float(trig.hypocentral_km.min()),
               "earliest_onset_sec": None if trig.empty else int(trig.onset_lag_sec.min()),
               "median_resid_p_sec": None if trig.empty else
               float(np.median(trig.onset_lag_sec - trig.hypocentral_km / V_P_KMS)),
               "median_resid_s_sec": None if trig.empty else
               float(np.median(trig.onset_lag_sec - trig.hypocentral_km / V_S_KMS))}
        if len(trig) >= MIN_STATIONS_FOR_FIT and trig.hypocentral_km.nunique() > 1:
            s, _, lo, hi = theilslopes(trig.onset_lag_sec, trig.hypocentral_km, alpha=CI_ALPHA)
            row.update({"slope_s_per_km": round(s, 4), "slope_ci_lo": round(lo, 4), "slope_ci_hi": round(hi, 4),
                        "apparent_velocity_kms": round(1 / s, 2) if s > 0 else None,
                        "verdict": classify_slope(lo, hi)})
        else:
            row["verdict"] = "too_few_stations"
        rows.append(row)
    return pd.DataFrame(rows)


def pooled_moveout(trig: pd.DataFrame, rng: np.random.Generator, n_perm: int = N_PERM) -> dict:
    """Theil-Sen on within-event centered (R, onset); p from shuffling onsets
    among stations of the same event, one-sided H1: slope > 0."""
    trig = trig.groupby(["group", "event_date"]).filter(lambda e: len(e) >= 2)
    if len(trig) < 6:
        return {"n_onsets": int(len(trig)), "note": "too few multi-station events"}
    key = trig.group + "|" + trig.event_date
    x = (trig.hypocentral_km - trig.groupby(key).hypocentral_km.transform("median")).to_numpy()
    y = trig.onset_lag_sec.to_numpy(dtype=float)
    ym = trig.groupby(key).onset_lag_sec.transform("median").to_numpy(dtype=float)
    s, _, lo, hi = theilslopes(y - ym, x, alpha=CI_ALPHA)
    idx = [np.flatnonzero((key == k).to_numpy()) for k in key.unique()]
    null = np.empty(n_perm)
    for i in range(n_perm):
        yp = y.copy()
        for ix in idx:
            yp[ix] = y[rng.permutation(ix)]
        null[i] = theilslopes(yp - ym, x)[0]
    p = float((1 + np.sum(null >= s)) / (n_perm + 1))
    return {"n_onsets": int(len(trig)), "n_events": len(idx),
            "slope_s_per_km": round(float(s), 4), "slope_ci90": [round(float(lo), 4), round(float(hi), 4)],
            "apparent_velocity_kms": round(1 / s, 2) if s > 0 else None,
            "p_one_sided_within_event": round(p, 5), "verdict": classify_slope(lo, hi)}


def arrival_window_test(trig: pd.DataFrame, rng: np.random.Generator) -> dict:
    """Do onsets land in the seismic arrival window more often than the same
    stations' own noise triggers do? Observed count vs. a Poisson-binomial
    null with per-trigger probability null_q_arrival (one-sided)."""
    if trig.empty:
        return {"n_onsets": 0}
    q = trig.null_q_arrival.to_numpy(dtype=float)
    obs = int(trig.in_arrival_window.astype(bool).sum())
    sims = (rng.random((N_SIM, len(q))) < q).sum(axis=1)
    return {"n_onsets": int(len(trig)), "n_in_window": obs, "expected_by_chance": round(float(q.sum()), 2),
            "p_one_sided": round(float((1 + np.sum(sims >= obs)) / (N_SIM + 1)), 5),
            "median_resid_s_sec": round(float(np.median(trig.onset_lag_sec - trig.hypocentral_km / V_S_KMS)), 2)}


def arrival_window_event_level(tested: pd.DataFrame, nulls: dict, rng: np.random.Generator) -> dict:
    """Primary arrival-window test. tested: rows with status onset/no_onset.
    T = onsets inside their own arrival window. Null: per simulation, one shared
    reference time per event; count that event's stations whose null trigger
    lands in their own window, sum over events (one-sided)."""
    if tested.empty:
        return {"n_events": 0}
    obs = int((tested.status.eq("onset") & tested.get("in_arrival_window", pd.Series(False, index=tested.index))
               .fillna(False).astype(bool)).sum())
    counts, indep_var = [], 0.0
    for (g, d), ev in tested.groupby(["group", "event_date"], sort=False):
        nl = nulls.get((g, d), {})
        arr = [nl[st] for st in ev.station if st in nl and len(nl[st])]
        if arr and len({len(a) for a in arr}) == 1:
            counts.append(np.sum(arr, axis=0))
            q = np.mean(arr, axis=1)
            indep_var += float(np.sum(q * (1 - q)))
    if not counts:
        return {"n_events": 0}
    sims = np.zeros(N_SIM)
    for c in counts:
        sims += c[rng.integers(0, len(c), N_SIM)]
    event_var = float(sum(c.var() for c in counts))
    return {"n_events": len(counts), "n_stations_tested": int(len(tested)), "n_in_window": obs,
            "expected_by_chance": round(float(sum(c.mean() for c in counts)), 2),
            "null_sd": round(float(sims.std()), 2),
            # null variance of the event sums over what independent stations would give:
            # > 1 means stations of an event fire together (external disturbances)
            "dispersion_ratio": round(event_var / indep_var, 3) if indep_var > 0 else None,
            "p_one_sided": round(float((1 + np.sum(sims >= obs)) / (N_SIM + 1)), 5)}


def timing_test_dispersion_adjusted(trig: pd.DataFrame, dispersion: float | None) -> dict:
    """The per-trigger timing question (given that a station triggered, does
    the onset land in its arrival window more often than its own null
    triggers?) with the Poisson-binomial variance inflated by the event-level
    dispersion ratio, so that stations of one event firing together are not
    counted as independent evidence. Normal approximation, one-sided."""
    if trig.empty or not dispersion:
        return {"n_onsets": int(len(trig))}
    q = trig.null_q_arrival.to_numpy(dtype=float)
    obs = int(trig.in_arrival_window.astype(bool).sum())
    sd = float(np.sqrt(max(dispersion, 1.0) * np.sum(q * (1 - q))))
    z = (obs - q.sum()) / sd if sd > 0 else float("nan")
    return {"n_onsets": int(len(trig)), "n_in_window": obs, "expected_by_chance": round(float(q.sum()), 2),
            "variance_inflation": round(max(dispersion, 1.0), 3), "z": round(z, 2),
            "p_one_sided_normal": float(f"{norm.sf(z):.3g}") if z == z else None}


def trigger_rate_by_distance(onsets: pd.DataFrame) -> list[dict]:
    tested = onsets[onsets.status.isin(["onset", "no_onset"])]
    out = []
    for b, d in tested.groupby(pd.cut(tested.hypocentral_km, [0, 30, 60, 100, 150, 250, 1000]), observed=False):
        if d.empty:
            continue
        out.append({"hypocentral_km": f"{b.left:g}-{b.right:g}", "n": int(len(d)),
                    "trigger_rate": round(float((d.status == "onset").mean()), 4),
                    "median_null_false_rate": round(float(d.null_false_rate.median()), 4)})
    return out


def seismometer_crosscheck(onsets: pd.DataFrame) -> dict:
    if not COMPARISON_CSV.exists():
        return {"note": "comparison_summary.csv missing"}
    cmp_ = pd.read_csv(COMPARISON_CSV)
    cmp_ = cmp_[(cmp_.status == "ok") & cmp_.shaking_onset_lag_sec.notna()]
    m = onsets[onsets.status == "onset"].merge(
        cmp_[["group", "date", "geomag_station", "seismic_station", "seismic_colocation_km", "shaking_onset_lag_sec"]],
        left_on=["group", "event_date", "station"], right_on=["group", "date", "geomag_station"])
    # seismometer traces start at origin-60 s; an onset at the trace start is a truncated pick, not an arrival
    m = m[m.shaking_onset_lag_sec > -55]
    return {"all_triggers": _crosscheck_stats(m),
            # xcg triggers at ~80 % of random times, so its "onset" is mostly the first noise sample of the window
            f"null_false_rate_lt_{CROSSCHECK_FALSE_RATE:g}": _crosscheck_stats(m[m.null_false_rate < CROSSCHECK_FALSE_RATE])}


def _crosscheck_stats(m: pd.DataFrame) -> dict:
    if m.empty:
        return {"n": 0}
    dt = m.onset_lag_sec - m.shaking_onset_lag_sec
    return {"n": int(len(m)), "n_m6": int((m.magnitude >= 6).sum()),
            "median_geomag_minus_seismic_sec": round(float(dt.median()), 2),
            "iqr_sec": [round(float(dt.quantile(0.25)), 2), round(float(dt.quantile(0.75)), 2)],
            "frac_within_5s": round(float((dt.abs() <= 5).mean()), 4),
            "frac_geomag_earlier_by_gt5s": round(float((dt < -5).mean()), 4)}


def summarize(onsets: pd.DataFrame, per_event: pd.DataFrame, nulls: dict) -> dict:
    trig = onsets[onsets.status == "onset"]
    m6, clean_prior = onsets.magnitude >= 6, ~onsets.prior_event_shaking.astype(bool)
    subsets = {"m6": m6, "m5": ~m6,
               "m6_no_prior_shaking": m6 & clean_prior, "m5_no_prior_shaking": ~m6 & clean_prior}
    excluded = onsets[~clean_prior].drop_duplicates(["group", "event_date"])
    out = {"seed": SEED, "k_sigma": K_SIGMA, "search_lags_sec": list(SEARCH_LAGS),
           "v_p_kms": V_P_KMS, "v_s_kms": V_S_KMS, "n_perm": N_PERM,
           "prior_event_exclusion_sec": PRIOR_EVENT_EXCLUSION_SEC,
           "prior_event_shaking_events": [f"{r.group} {r.event_date} M{r.magnitude:g} ({r.prior_event_sec:.0f} s)"
                                          for r in excluded.itertuples()],
           "n_station_tests": int(onsets.status.isin(["onset", "no_onset"]).sum()),
           "n_onsets": int(len(trig)),
           "overall_null_false_rate": round(float(onsets.null_false_rate.mean()), 4),
           "subsets": {}}
    for name, mask in subsets.items():
        t = trig[mask[trig.index]]
        tested = onsets[mask & onsets.status.isin(["onset", "no_onset"])]
        ev_all = arrival_window_event_level(tested, nulls, keyed_rng("onset", name, "window_event_all"))
        ev_clean = arrival_window_event_level(tested[tested.null_false_rate < CLEAN_FALSE_RATE], nulls,
                                              keyed_rng("onset", name, "window_event_clean"))
        pe = per_event[(per_event.magnitude >= 6) if name.startswith("m6") else (per_event.magnitude < 6)]
        if name.endswith("no_prior_shaking"):
            pe = pe[~pe.prior_event_shaking]
        out["subsets"][name] = {
            "n_events": int(len(pe)),
            "event_verdicts": pe.verdict.value_counts().to_dict(),
            # primary: event-level null (stations of an event share reference times)
            "arrival_window_event_level_all": ev_all,
            "arrival_window_event_level_clean": ev_clean,
            # timing only (conditional on triggering), variance inflated by the event-level dispersion
            "timing_dispersion_adjusted_all": timing_test_dispersion_adjusted(t, ev_all.get("dispersion_ratio")),
            "timing_dispersion_adjusted_clean": timing_test_dispersion_adjusted(
                t[t.null_false_rate < CLEAN_FALSE_RATE], ev_clean.get("dispersion_ratio")),
            "n_events_shared_null_fallback": int(onsets[mask & onsets.shared_null_fallback.eq(True)]
                                                 .drop_duplicates(["group", "event_date"]).shape[0]),
            # reference: per-trigger Poisson-binomial (treats stations as independent)
            "arrival_window_all_triggers": arrival_window_test(t, keyed_rng("onset", name, "window_all")),
            "arrival_window_clean_stations": arrival_window_test(t[t.null_false_rate < CLEAN_FALSE_RATE],
                                                                 keyed_rng("onset", name, "window_clean")),
            "pooled_all_triggers": pooled_moveout(t, keyed_rng("onset", name, "pooled_all")),
            "pooled_clean_stations": pooled_moveout(t[t.null_false_rate < CLEAN_FALSE_RATE],
                                                    keyed_rng("onset", name, "pooled_clean")),
            "trigger_rate_by_distance": trigger_rate_by_distance(onsets[mask]),
        }
    out["seismometer_crosscheck"] = seismometer_crosscheck(onsets)
    return out


def plot(onsets: pd.DataFrame, path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    trig = onsets[onsets.status == "onset"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    r = np.linspace(0, max(300, trig.hypocentral_km.max() if len(trig) else 300), 50)
    for ax, (label, mask) in zip(axes, [("M ≥ 6", trig.magnitude >= 6), ("M < 6", trig.magnitude < 6)]):
        d = trig[mask]
        clean = d.null_false_rate < CLEAN_FALSE_RATE
        ax.scatter(d.hypocentral_km[clean], d.onset_lag_sec[clean], s=14, label="station false rate < 5%")
        ax.scatter(d.hypocentral_km[~clean], d.onset_lag_sec[~clean], s=14, marker="x", alpha=0.6,
                   label="noisier station")
        ax.plot(r, r / V_P_KMS, "--", lw=1, label=f"P ({V_P_KMS:g} km/s)")
        ax.plot(r, r / V_S_KMS, ":", lw=1.5, label=f"S ({V_S_KMS:g} km/s)")
        ax.axhline(0, color="gray", lw=0.6)
        ax.set_title(f"{label} (n={len(d)})")
        ax.set_xlabel("hypocentral distance (km)")
    axes[0].set_ylabel("geomagnetic onset lag after origin (s)")
    axes[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def run_all(group_ids: list[str] | None, min_mag: float | None) -> dict:
    group_ids = group_ids or list(GROUPS)
    with ProcessPoolExecutor(max_workers=N_WORKERS) as ex:
        results = list(ex.map(process_group, group_ids, [min_mag] * len(group_ids)))
    onsets = pd.DataFrame([r for rs, _ in results for r in rs])
    nulls = {k: v for _, ns in results for k, v in ns.items()}
    per_event = per_event_moveout(onsets)
    summary = summarize(onsets, per_event, nulls)
    summary.update({"groups": group_ids, "min_mag": min_mag})
    out_dir = OUT_DIR if group_ids == list(GROUPS) and min_mag is None else \
        OUT_DIR / ("subset_" + "_".join(group_ids) + (f"_m{min_mag:g}" if min_mag is not None else ""))
    out_dir.mkdir(parents=True, exist_ok=True)
    onsets.to_csv(out_dir / "onsets.csv", index=False)
    per_event.to_csv(out_dir / "moveout_per_event.csv", index=False)
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False, default=str))
    plot(onsets, out_dir / "travel_time.png")
    print(json.dumps(summary, indent=2, ensure_ascii=False, default=str), file=sys.stderr)
    return summary


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------

def _synthetic(delays: list[float], rng: np.random.Generator) -> list[int | None]:
    out = []
    for dly in delays:
        v = np.cumsum(rng.normal(0, 0.02, 3600))
        p0 = 2000
        s = p0 + int(round(dly))
        v[s:s + 60] += np.cumsum(rng.normal(0, 0.3, 60))
        r = detect_onset(v, p0, SEARCH_LAGS[1])
        out.append(r.get("onset_lag_sec"))
    return out


def _correlated_null_trial(rng: np.random.Generator, n_ev: int = 30, n_st: int = 6, k: int = N_NULL):
    """One H0 dataset where half of all reference times carry a network-wide
    disturbance (every station triggers at the same lag) and otherwise each
    station triggers alone 10 % of the time -- the origin is just another
    reference time. Returns (rows, nulls) shaped like process_group's."""
    rows, nulls = [], {}
    for e in range(n_ev):
        R = rng.uniform(20, 300, n_st)
        lo, hi = R / V_P_KMS + ARRIVAL_PAD_SEC[0], R / V_S_KMS + ARRIVAL_PAD_SEC[1]
        lag = np.full((n_st, k + 1), np.nan)
        common_t = rng.random(k + 1) < 0.5
        lag[:, common_t] = rng.uniform(SEARCH_LAGS[0], SEARCH_LAGS[1], common_t.sum())
        own = (rng.random((n_st, k + 1)) < 0.1) & ~common_t
        lag[own] = rng.uniform(SEARCH_LAGS[0], SEARCH_LAGS[1], own.sum())
        inw = (lag >= lo[:, None]) & (lag <= hi[:, None])
        nulls[("S", str(e))] = {f"s{i}": inw[i, 1:] for i in range(n_st)}
        for i in range(n_st):
            on = lag[i, 1:][~np.isnan(lag[i, 1:])]
            trig = not np.isnan(lag[i, 0])
            rows.append({"group": "S", "event_date": str(e), "station": f"s{i}",
                         "status": "onset" if trig else "no_onset", "hypocentral_km": R[i],
                         "onset_lag_sec": lag[i, 0], "in_arrival_window": bool(inw[i, 0]) if trig else None,
                         "null_q_arrival": float(np.mean((on >= lo[i]) & (on <= hi[i]))) if len(on) else 0.0})
    return pd.DataFrame(rows), nulls


def self_test() -> bool:
    rng = np.random.default_rng(SEED)
    R = np.array([20, 45, 80, 120, 160, 210, 260], dtype=float)
    ok = True
    for name, delays, want in [("simultaneous", np.full(len(R), 12.0), "simultaneous"),
                               ("S-wave moveout", R / V_S_KMS, "seismic_moveout")]:
        on = _synthetic(list(delays), rng)
        if any(o is None for o in on):
            print(f"[self-test] {name}: missed onsets {on}  FAIL", file=sys.stderr)
            ok = False
            continue
        s, _, lo, hi = theilslopes(on, R, alpha=CI_ALPHA)
        verdict = classify_slope(lo, hi)
        err = max(abs(o - d) for o, d in zip(on, delays))
        good = verdict == want and err <= 6  # emergent burst: first samples can sit below K_SIGMA
        print(f"[self-test] {name}: slope={s:.3f} s/km CI=[{lo:.3f},{hi:.3f}] -> {verdict}, "
              f"max onset error {err:.1f}s  {'PASS' if good else 'FAIL'}", file=sys.stderr)
        ok &= good
    quiet = [detect_onset(np.cumsum(rng.normal(0, 0.02, 3600)), 2000, SEARCH_LAGS[1]) for _ in range(200)]
    fr = np.mean([q["status"] == "onset" for q in quiet])
    print(f"[self-test] pure-noise false trigger rate {fr:.3f}  {'PASS' if fr < 0.02 else 'FAIL'}", file=sys.stderr)
    ok &= fr < 0.02
    fp_event = fp_station = 0
    n_rep = 200
    for _ in range(n_rep):
        df, nulls = _correlated_null_trial(rng)
        fp_event += arrival_window_event_level(df, nulls, rng)["p_one_sided"] < 0.05
        fp_station += arrival_window_test(df[df.status == "onset"], rng)["p_one_sided"] < 0.05
    good = fp_event / n_rep <= 0.08 and fp_station > fp_event
    print(f"[self-test] correlated H0 (network-wide disturbances), false positives at 5 %: "
          f"event-level {fp_event / n_rep:.3f}, per-station {fp_station / n_rep:.3f}  "
          f"{'PASS' if good else 'FAIL'}", file=sys.stderr)
    ok &= good
    t0 = pd.Timestamp("2024-04-03 00:00:00")
    gaps = [prior_event_gap(t0 + pd.Timedelta(seconds=d), [t0]) for d in (PRIOR_EVENT_EXCLUSION_SEC,
                                                                        PRIOR_EVENT_EXCLUSION_SEC + 1, 0)]
    flags = [g is not None and g <= PRIOR_EVENT_EXCLUSION_SEC for g in gaps]
    good = flags == [True, False, False]
    print(f"[self-test] prior-event flag at 900/901/0 s: {flags}  {'PASS' if good else 'FAIL'}", file=sys.stderr)
    ok &= good
    return bool(ok)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--all", action="store_true", help="run on real data (default action)")
    ap.add_argument("--group", action="append", help="restrict to this group (repeatable)")
    ap.add_argument("--min-mag", type=float, default=None)
    args = ap.parse_args()
    if args.self_test:
        sys.exit(0 if self_test() else 1)
    if not self_test():
        print("[main] self-test FAILED -- aborting before touching real data", file=sys.stderr)
        sys.exit(1)
    run_all(args.group, args.min_mag)
