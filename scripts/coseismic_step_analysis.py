"""Coseismic (instantaneous, at-origin-time) geomagnetic step/spike detector.

Everything else in this codebase (`build_daily_features.py`, `ulf_analysis.py`,
`compute_indices.py`, `superposed_epoch_analysis.py`, ...) works at daily or
1-minute-resampled granularity and asks "is there a precursor/aftermath
signature over days-to-weeks around an earthquake". This script asks a
narrower question: **at the earthquake's origin second itself, is there a
step-like or spike-like jump in the raw 1Hz X/Y/Z/F field**, the kind of
signal a piezomagnetic-effect or seismic-wave-induction mechanism could in
principle produce (see e.g. Tohoku 2011 Mw9.0 reports of a few-nT coseismic
pulse at near-field stations) -- or, just as plausibly, shaking-induced
sensor-housing noise rather than a genuine field change (see "What the G10
pilot run found" below).

## History / scope

Phase 1 (single-flagship-event pilot, see `~/.claude/plans/recursive-
seeking-nova.md`): only G10's anchor event (2024-04-03 07:58:11 Taiwan time
/ 2024-04-02 23:58:11 UTC, M7.2) against its two nearest XYZ stations.
Phase 2 (current): extended to every event in `events.py`'s 117-event
registry via `run_all()` / `--all` -- see "Extending to all 117 events"
below for exactly what changed and why.

### What the G10 pilot run found (motivates why this needs to be run at
### more than one event before drawing any general conclusion)

Both of G10's nearest stations (xcg 24.8km, cnu 66.3km) showed statistically
significant step/spike anomalies within about a minute of the origin time
(several p <= 0.0025, surviving Bonferroni correction at that scale). But
inspecting the raw waveform directly showed this was *not* a clean DC step:
sample-to-sample noise amplitude jumped to ~11.6-11.7x the pre-event
baseline, starting a few seconds later at the near station than the far
station -- consistent with the arrival of strong ground shaking physically
jostling the magnetometer housing, not a piezomagnetic-type field change.
Whether this "shaking-noise" pattern is the general story across the
dataset, or G10 (the largest, best-instrumented event) was a special case,
is exactly what the 117-event run is for.

## Why this needs new machinery

`parser.parse_day_file` and `common.resolve_day_ref` only address whole
station-days -- there is no sub-range read and no cross-midnight stitching
built anywhere in this codebase (every existing caller works one calendar
day at a time). This script adds that: it loads the full day(s) of raw data
spanning the origin timestamp's required buffer (which, for an origin time
close to UTC midnight like G10's, is two station-days), concatenates them,
and works on the resulting continuous 1Hz series directly -- there was no
computational reason to trim to a sub-window first, since `parse_day_file`
parses an entire 86400-row file regardless of how much of it is used
afterward.

## Statistical design (mirrors `surrogate_test.py` / `superposed_epoch_
analysis.py`'s established pattern: observed statistic + null distribution
built by resampling + p-value + fixed seed, rather than inventing a new
methodology)

1. Detrend each channel with the same 1-hour centered rolling-mean
   subtraction as `ulf_analysis.py::_detrend` (kept here as a local copy,
   documented as such, rather than importing from a script whose own
   purpose is the ULF near/far pipeline -- see plan file for why nothing in
   `ulf_analysis.py` was modified).
2. Two complementary statistics on the detrended series `d`:
   - **step**: `S_M(t) = mean(d[t+1 .. t+M]) - mean(d[t-M+1 .. t])` at
     `M in {10, 30, 90}` seconds, a sliding step/offset detector.
   - **spike**: `|diff(d)[t]|`, a single-sample jump detector (same idea as
     `build_daily_features.py::_despike`'s threshold, reused here as the
     statistic itself rather than a noise filter).
3. **Observed value**: `max(|statistic|)` within +-SCAN_HALF_SEC (180s, capped
   per-event -- see `_effective_half_sec`) of the reported origin second (a
   search window, not a single point -- allows for reporting-time rounding
   and any real propagation delay).
4. **Null distribution**: the *same* "search +-SCAN_HALF_SEC and take the extremum"
   procedure repeated at N_NULL=2000 random reference times drawn from
   elsewhere in the same loaded window (>=600s away from the real origin
   time), so the null is judged on an apples-to-apples footing with the
   observed value rather than a single-point null.
5. p-value = (1 + #{|null| >= |obs|}) / (N_NULL + 1), one-sided since both
   statistics are already |.|.
6. Every result also reports `approx_min_detectable_nt` (the null
   distribution's 95th percentile) -- so a non-significant result can be
   read honestly as either "well-powered, no signal" (this number is small)
   or "this station-day's noise floor was too high to say anything" (this
   number is large), instead of collapsing both into "not significant".

## Known limitations (stated up front, not as an excuse after the fact)

- Physically, coseismic geomagnetic transients (piezomagnetic effect,
  seismic-wave-induced currents) are really only well-documented in the
  literature at near-field stations (single-digit km) for very large
  (M8+) events. G10's M7.2 mainshock and its ~10-30km-distant nearest
  stations put this well outside that regime -- **a null result here is the
  expected outcome**, not evidence the mechanism doesn't exist elsewhere.
- Single event, 2 stations, ~30ish tests total (see
  `multiple_comparisons_context` in the output) -- no multiple-comparisons
  correction is applied (not needed at this scale), but the count is
  reported so it can be compared against a future multi-event expansion.
- No geomagnetic-storm exclusion, by design: storm main/recovery phases
  vary over hours-to-days and are mostly removed by the 1hr detrend, and
  the null is drawn from the same loaded window, so a storm-elevated noise
  floor raises obs and null together. The residual risk is a single
  external transient (SSC, substorm onset, Pi2) inside the search window.
  Each event is flagged instead (`is_storm_day` / `is_storm_onset` from the
  group's `storm_days.csv`, see `_storm_status`) and `--all` reports a
  storm vs. quiet split in `all_events_run_summary.json::storm_sensitivity`.
  Excluding flagged events outright isn't practical: with the daily
  pipeline's Kp>=5 / Dst<=-30nT (+recovery days) definition, about half of
  all events fall on a flagged day.
- The step-statistic's "post minus pre" implementation (see
  `_step_statistic`) has a known ~1-sample edge-alignment approximation from
  the reverse-rolling trick used to compute a forward-looking mean; this is
  irrelevant at the +-SCAN_HALF_SEC (180s) search-window scale used here.

## Extending to all 117 events (done -- `run_all()` / `--all`)

- Uses `events.py`'s `GROUPS` dict and each `Event.time_utc` directly (NOT
  `catalog_utils.load_extended_events()`, which truncates to date-only
  resolution) -- every event in every group, not just the 20 `anchor=True`
  events, since each event's coseismic test is its own independent
  measurement (unlike the daily-scale cross-group statistics, there's no
  foreshock/mainshock/aftershock pseudo-replication concern here).
- Station distance is only pre-computed relative to each group's single
  anchor event (`common.py::_discover_stations`); `_rank_stations_for_event`
  recomputes `common.haversine_km(station_lat, station_lon, ev.lat, ev.lon)`
  fresh for every individual event so non-anchor events get correct
  near-station rankings.
- Every event is tested against *both* the F-only and XYZ station pools
  when the group has stations in both (e.g. G2_G3 has 10 F stations and 2
  XYZ stations -- both get tested, not one at the other's expense); most
  groups only have one non-empty pool (see `common.py`'s per-group
  F/XYZ-pool printout), so this mostly reduces to "whichever pool exists".
- Null-candidate exclusion covers every event in the same group (not just
  the one being tested) within `EXCLUSION_BUFFER_SEC`, since some groups
  have events close enough in time (G9's two events 17h apart, G10's
  2024-04-23a/b six minutes apart) that a naive single-event exclusion
  could let the null distribution get contaminated by a different real
  earthquake landing inside the same loaded window.
- The *observed*-value search window has the same contamination risk as the
  null side, and is capped per-event to guard against it (`_effective_half_sec`):
  G10's 2024-04-23a/b are only 357s apart, so once SCAN_HALF_SEC exceeds
  ~178s an uncapped +-SCAN_HALF_SEC search around one event would reach past
  the midpoint into the other event's own anomaly. `_effective_half_sec`
  caps each event's actual search window at `min(SCAN_HALF_SEC, nearest
  same-group event gap // 2)`, so only this one pair is ever narrowed (every
  other group's events are hours-to-years apart); the capped value is
  recorded per-event in `scan_half_sec` for transparency.
- The O(n) boolean-mask scan in the original `_extremum_in_window` (full
  ~172,800-sample array masked per null draw) was the single-event pilot's
  runtime bottleneck (~20 minutes for 2 stations); replaced with O(window)
  positional-index slicing (`_pos` + the current `_extremum_in_window`),
  since a day's data is guaranteed a complete 86400-row grid and
  concatenated adjacent days have no gap, so timestamp arithmetic against
  the array's own start gives an exact integer array position.
- Positive-control injection testing is NOT repeated across all 117 events
  (see plan file) -- the method was already validated once on G10
  (`injection_power_curve.json`, untouched by `run_all()`); `--all` always
  runs with `run_injection=False`.
- Cross-event stacking (normalizing each event's window via
  `stat_utils.mad_zscore` on its own off-event statistic values, then
  averaging across events the way `superposed_epoch_analysis.py` stacks
  day-scale series, with a null band built from random per-event reference
  times) is the natural next step if the 117-event run shows many
  individually-underpowered-but-suggestive events, but is deliberately not
  implemented in this pass -- see plan file.

Usage:
  coseismic_step_analysis.py                  # G10 anchor only, both pools, full run (smoke test)
  coseismic_step_analysis.py --all             # all 117 events in events.py's registry
  coseismic_step_analysis.py --self-test       # synthetic white-noise sanity check only
  coseismic_step_analysis.py --no-injection    # skip the positive-control injection test (G10-only path)
"""
from __future__ import annotations

import argparse
import json
import sys
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import common  # noqa: E402
import parser as sec_parser  # noqa: E402
from events import GROUPS, folder_events, get_group  # noqa: E402
from stat_utils import histogram_summary  # noqa: E402

SEED = 20260805


def keyed_rng(*keys) -> np.random.Generator:
    """Independent, stable random stream per key tuple, e.g. (group, event,
    station, channel). One shared stream across all events made every null
    draw depend on how many draws earlier events had consumed, so editing
    any one event (or running a single event instead of --all) reshuffled
    every later event's p-value."""
    return np.random.default_rng([SEED, zlib.crc32("|".join(map(str, keys)).encode())])

GROUP_ID = "G10"
N_STATIONS = 2

STEP_WINDOWS_SEC = (10, 30, 90)
SCAN_HALF_SEC = 180          # target: search this far either side of the reported origin second
                              # (widened 120->180 2026-08-16: was pinning G12/G13/G20's obs_lag at the
                              # old +-120s boundary. 2026-08-19: investigated widening further to 300s
                              # to check whether G20's 2025-12-24 obs_lag=-178s (2s from this 180s
                              # boundary) was itself a truncation artifact -- see plan
                              # artifact-wobbly-kettle.md. It was NOT: rerunning at 300s left G20's
                              # obs_lag unchanged at -178s, confirming it's a genuine local extremum,
                              # not a boundary cutoff. But the 300s trial also *changed* the reported
                              # obs_lag for several OTHER, non-boundary-pinned events (e.g. G1
                              # 2018-02-06: -133s at 180s -> 243s at 300s, z=2.84 -> z=3.43) -- a wider
                              # search window finds a larger max-of-N by chance even under noise alone
                              # (the null distribution widens correspondingly, so p-values stay valid,
                              # but *which* lag gets reported as "the" peak becomes window-size-dependent
                              # and less physically interpretable the further it sits from the origin
                              # second). Reverted to 180s given G20 didn't need the wider window after
                              # all. This is a per-event TARGET, not a hard value -- see
                              # `_effective_half_sec`, kept from the 300s trial as a real (if narrow)
                              # fix: G10's 2024-04-23a/b are only 357s apart, so their search window is
                              # capped at ~178s regardless of this target, to avoid one event's search
                              # reaching into the other's real anomaly.)
DETREND_WINDOW_SEC = 3601    # same as ulf_analysis.py::_detrend (1hr centered rolling mean)
REQUIRED_MARGIN_SEC = SCAN_HALF_SEC + max(STEP_WINDOWS_SEC)  # 270s: minimum data support at either edge
                              # (uses the SCAN_HALF_SEC *target*, not the per-event capped value, so this
                              # margin check stays conservative/uniform regardless of any per-event cap)
N_NULL = 2000
EXCLUSION_BUFFER_SEC = 600   # keep null reference times this far from the real origin time
INJECTION_EXCLUSION_BUFFER_SEC = 2000  # wider, for the positive-control run (see run_injection_power_curve)
INJECTION_AMPLITUDES_SIGMA = (1, 2, 3, 5, 10)  # multiples of the local null MAD
MISSING_FRACTION_THRESHOLD = 0.10  # NaN fraction within the +-SCAN_HALF_SEC window that voids a test

OUT_DIR = common.PROJECT_DIR / "data" / "interim" / "coseismic_step_analysis"


# ---------------------------------------------------------------------------
# Signal processing (local copies / adaptations -- see module docstring for
# why these aren't imported from ulf_analysis.py / build_daily_features.py)
# ---------------------------------------------------------------------------

def _detrend(x: pd.Series) -> np.ndarray:
    """Same formula as ulf_analysis.py::_detrend: subtract a 1-hour centered
    rolling mean. Local copy, not an import -- see module docstring."""
    trend = x.rolling(window=DETREND_WINDOW_SEC, center=True, min_periods=1).mean()
    return (x - trend).to_numpy(dtype=float)


def _step_statistic(d: np.ndarray, half_win: int) -> np.ndarray:
    """S_M(t) = mean(d[t+1..t+M]) - mean(d[t-M+1..t]). See module docstring
    for the ~1-sample edge-alignment caveat of the reverse-rolling trick
    used for the forward-looking mean; irrelevant at the +-SCAN_HALF_SEC
    (180s) search-window scale this script uses."""
    s = pd.Series(d)
    min_p = max(3, int(0.6 * half_win))
    pre = s.rolling(half_win, min_periods=min_p).mean()
    post = (
        s[::-1].reset_index(drop=True)
        .rolling(half_win, min_periods=min_p).mean()
        [::-1].reset_index(drop=True).shift(-1)
    )
    post.index = s.index
    return (post - pre).to_numpy(dtype=float)


def _spike_statistic(d: np.ndarray) -> np.ndarray:
    """|diff(d)[t]|, same jump-detection idea as build_daily_features.py::
    _despike, used here as the statistic itself rather than a noise filter."""
    return np.abs(np.diff(d, prepend=d[0]))


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def _load_station_days(gdms_dir: Path, station: str, event_utc: pd.Timestamp,
                        buffer_sec: int) -> tuple[pd.DataFrame | None, list[str], list[str]]:
    """Load and concatenate whichever station-day file(s) cover
    [event_utc - buffer_sec, event_utc + buffer_sec]. There is no sub-range
    read in parser.py, so each covering day is parsed in full (86400 rows)
    and the result sliced afterward by the caller. Returns
    (concatenated_df_or_None, missing_date_strs, wanted_date_strs)."""
    lo = event_utc - pd.Timedelta(seconds=buffer_sec)
    hi = event_utc + pd.Timedelta(seconds=buffer_sec)
    wanted = sorted({d.strftime("%Y%m%d") for d in pd.date_range(lo.normalize(), hi.normalize(), freq="D")})
    frames = []
    missing = []
    for date_str in wanted:
        ref = common.resolve_day_ref(gdms_dir, station, date_str)
        if ref is None:
            missing.append(date_str)
            continue
        frames.append(sec_parser.parse_day_file(ref, station))
    if not frames:
        return None, missing, wanted
    df = pd.concat(frames).sort_index()
    return df, missing, wanted


# ---------------------------------------------------------------------------
# Extremum-in-window helper (shared by observed value and null draws, so
# both sides of the comparison go through identical code). Positional
# (integer-second-offset) slicing, not a boolean mask over the whole array:
# a day's data is a guaranteed-complete 86400-row grid and concatenated
# adjacent days have no gap, so `index[0]` + integer seconds gives an exact
# array position -- O(window) per call instead of O(len(index)). The
# original boolean-mask version was this script's runtime bottleneck at
# 117-event scale (see module docstring's "Extending to all 117 events").
# ---------------------------------------------------------------------------

def _pos(index: pd.DatetimeIndex, t: pd.Timestamp) -> int:
    return int(round((t - index[0]).total_seconds()))


def _extremum_in_window(stat_values: np.ndarray, index: pd.DatetimeIndex,
                         center: pd.Timestamp, half_sec: int) -> tuple[float, int] | None:
    p = _pos(index, center)
    lo, hi = max(0, p - half_sec), min(len(stat_values), p + half_sec + 1)
    sub = stat_values[lo:hi]
    if len(sub) == 0 or np.all(np.isnan(sub)):
        return None
    j = int(np.nanargmax(np.abs(sub)))
    return float(sub[j]), int((lo + j) - p)


def _draw_null_centers(rng: np.random.Generator, valid_lo: pd.Timestamp, valid_hi: pd.Timestamp,
                        exclude_centers: list[pd.Timestamp], exclude_buffer_sec: int, n: int) -> list[pd.Timestamp]:
    """exclude_centers is a list (not a single timestamp) so a candidate is
    rejected if it lands near *any* real event in the group being tested --
    needed because some groups have events close enough in time (see module
    docstring) that a single-event exclusion could let a null draw land on
    top of a different real earthquake."""
    total_sec = (valid_hi - valid_lo).total_seconds()
    centers: list[pd.Timestamp] = []
    attempts = 0
    while len(centers) < n and attempts < n * 5:
        attempts += 1
        cand = valid_lo + pd.Timedelta(seconds=rng.uniform(0, total_sec))
        if any(abs((cand - ec).total_seconds()) < exclude_buffer_sec for ec in exclude_centers):
            continue
        centers.append(cand)
    return centers


def _effective_half_sec(target_half_sec: int, event_utc: pd.Timestamp,
                         exclude_centers: list[pd.Timestamp]) -> int:
    """Cap the search/scan half-window so it never reaches past the midpoint
    to the nearest *other* real event in the same group. Needed because
    SCAN_HALF_SEC is a shared target across all 117 events, but G10's
    2024-04-23a/b are only 357s apart -- once the target exceeds ~178s, an
    uncapped +-SCAN_HALF_SEC search around one of them would reach into the
    other's own real anomaly, contaminating the "observed value" search
    (the null side already handles this via exclude_centers in
    `_draw_null_centers`; this is the analogous guard for the observed side).
    Every other group's events are hours-to-years apart, so this is a no-op
    everywhere except that one pair."""
    others = [c for c in exclude_centers if c != event_utc]
    if not others:
        return target_half_sec
    nearest_gap_sec = min(abs((event_utc - c).total_seconds()) for c in others)
    return int(min(target_half_sec, nearest_gap_sec // 2))


# ---------------------------------------------------------------------------
# Per-(station, channel) test
# ---------------------------------------------------------------------------

def _test_channel(raw: pd.Series, event_utc: pd.Timestamp, exclude_centers: list[pd.Timestamp],
                   rng: np.random.Generator, channel_label: str, scan_half_sec: int) -> dict:
    idx = raw.index
    nan_frac_at_event = float(
        raw[(idx >= event_utc - pd.Timedelta(seconds=scan_half_sec)) &
            (idx <= event_utc + pd.Timedelta(seconds=scan_half_sec))].isna().mean()
    ) if len(raw) else 1.0

    if idx.min() + pd.Timedelta(seconds=REQUIRED_MARGIN_SEC) > event_utc or \
       idx.max() - pd.Timedelta(seconds=REQUIRED_MARGIN_SEC) < event_utc:
        return {"channel": channel_label, "data_status": "insufficient_margin"}
    if nan_frac_at_event > MISSING_FRACTION_THRESHOLD:
        return {"channel": channel_label, "data_status": "insufficient_data",
                "nan_fraction_at_event": round(nan_frac_at_event, 4)}

    d = _detrend(raw)
    valid_lo = idx.min() + pd.Timedelta(seconds=scan_half_sec)
    valid_hi = idx.max() - pd.Timedelta(seconds=scan_half_sec)
    null_centers = _draw_null_centers(rng, valid_lo, valid_hi, exclude_centers, EXCLUSION_BUFFER_SEC, N_NULL)

    out = {"channel": channel_label, "data_status": "ok",
           "nan_fraction_at_event": round(nan_frac_at_event, 4), "scan_half_sec": scan_half_sec, "step": {}}

    for m in STEP_WINDOWS_SEC:
        stat = _step_statistic(d, m)
        obs = _extremum_in_window(stat, idx, event_utc, scan_half_sec)
        null_vals = [r for c in null_centers if (r := _extremum_in_window(stat, idx, c, scan_half_sec)) is not None]
        out["step"][str(m)] = _summarize_test(obs, null_vals)

    spike_stat = _spike_statistic(d)
    obs = _extremum_in_window(spike_stat, idx, event_utc, scan_half_sec)
    null_vals = [r for c in null_centers if (r := _extremum_in_window(spike_stat, idx, c, scan_half_sec)) is not None]
    out["spike"] = _summarize_test(obs, null_vals)

    return out


def _summarize_test(obs: tuple[float, int] | None, null_vals: list[tuple[float, int]]) -> dict:
    if obs is None or not null_vals:
        return {"data_status": "no_observed_or_null_value"}
    obs_val, obs_lag = obs
    null_abs = np.array([abs(v) for v, _ in null_vals])
    n_used = len(null_abs)
    p_value = (1 + int(np.sum(null_abs >= abs(obs_val)))) / (n_used + 1)
    return {
        "obs_value_nt": round(obs_val, 4),
        "obs_lag_sec": obs_lag,
        "p_value": round(p_value, 5),
        "n_null_samples": n_used,
        "approx_min_detectable_nt": round(float(np.percentile(null_abs, 95)), 4),
        "null_summary": histogram_summary(null_abs),
    }


# ---------------------------------------------------------------------------
# Positive control: inject a synthetic pulse into real background noise and
# confirm the detector's p-value drops monotonically with injected amplitude.
# ---------------------------------------------------------------------------

def run_injection_power_curve(raw: pd.Series, exclude_centers: list[pd.Timestamp],
                               rng: np.random.Generator, channel_label: str) -> dict:
    idx = raw.index
    valid_lo = idx.min() + pd.Timedelta(seconds=REQUIRED_MARGIN_SEC + INJECTION_EXCLUSION_BUFFER_SEC)
    valid_hi = idx.max() - pd.Timedelta(seconds=REQUIRED_MARGIN_SEC + INJECTION_EXCLUSION_BUFFER_SEC)
    # Pick an injection point far from the real event(s) and from the array edges.
    candidates = _draw_null_centers(rng, valid_lo, valid_hi, exclude_centers, INJECTION_EXCLUSION_BUFFER_SEC, 1)
    if not candidates:
        return {"channel": channel_label, "status": "no_valid_injection_point"}
    inject_at = candidates[0]
    inject_exclude = exclude_centers + [inject_at]

    # Establish the local noise floor (sigma) from an un-injected null run at this location.
    baseline_d = _detrend(raw)
    baseline_null_centers = _draw_null_centers(rng, valid_lo, valid_hi, inject_exclude,
                                                INJECTION_EXCLUSION_BUFFER_SEC, N_NULL)
    m0 = STEP_WINDOWS_SEC[0]
    baseline_stat = _step_statistic(baseline_d, m0)
    baseline_null = [abs(r[0]) for c in baseline_null_centers
                      if (r := _extremum_in_window(baseline_stat, idx, c, SCAN_HALF_SEC)) is not None]
    sigma = float(np.median(baseline_null)) if baseline_null else 1.0
    if sigma <= 0:
        sigma = 1.0

    duration_sec = 2 * m0  # short plateau: a rise then a fall, both step-like
    curve = []
    for k in INJECTION_AMPLITUDES_SIGMA:
        amplitude = k * sigma
        injected = raw.copy()
        mask = (idx >= inject_at) & (idx < inject_at + pd.Timedelta(seconds=duration_sec))
        injected.loc[mask] = injected.loc[mask] + amplitude
        d = _detrend(injected)
        stat = _step_statistic(d, m0)
        obs = _extremum_in_window(stat, idx, inject_at, SCAN_HALF_SEC)
        null_centers = _draw_null_centers(rng, valid_lo, valid_hi, inject_exclude,
                                           INJECTION_EXCLUSION_BUFFER_SEC, N_NULL)
        null_vals = [r for c in null_centers
                     if (r := _extremum_in_window(stat, idx, c, SCAN_HALF_SEC)) is not None]
        result = _summarize_test(obs, null_vals)
        result["injected_amplitude_nt"] = round(amplitude, 4)
        result["injected_amplitude_sigma_multiple"] = k
        curve.append(result)

    return {
        "channel": channel_label,
        "status": "ok",
        "inject_at_utc": inject_at.strftime("%Y-%m-%d %H:%M:%S"),
        "local_sigma_nt": round(sigma, 4),
        "step_window_sec": m0,
        "duration_sec": duration_sec,
        "curve": curve,
    }


# ---------------------------------------------------------------------------
# Synthetic self-test (no real data): confirms the statistic + null/p-value
# machinery itself is correct, run via --self-test.
# ---------------------------------------------------------------------------

def self_test() -> bool:
    rng = np.random.default_rng(SEED)
    n = 3600
    noise = rng.normal(0, 1.0, size=n)
    step_at = 1800
    amplitude = 20.0
    injected = noise.copy()
    injected[step_at:] += amplitude

    ok = True
    for m in STEP_WINDOWS_SEC:
        stat = _step_statistic(injected, m)
        peak_idx = int(np.nanargmax(np.abs(stat)))
        peak_val = stat[peak_idx]
        near_step = abs(peak_idx - step_at) <= m + 2
        magnitude_ok = abs(peak_val) > amplitude * 0.5  # boxcar smoothing softens the peak; sanity bound only
        status = "PASS" if (near_step and magnitude_ok) else "FAIL"
        if status == "FAIL":
            ok = False
        print(f"[self-test] M={m:3d}s  peak_idx={peak_idx} (expected ~{step_at})  "
              f"peak_val={peak_val:.2f} (injected {amplitude})  {status}")

    spike = _spike_statistic(injected)
    spike_peak = int(np.nanargmax(spike))
    spike_status = "PASS" if abs(spike_peak - step_at) <= 2 else "FAIL"
    if spike_status == "FAIL":
        ok = False
    print(f"[self-test] spike peak_idx={spike_peak} (expected ~{step_at})  {spike_status}")
    return ok


# ---------------------------------------------------------------------------
# Per-event orchestration: tests one Event against every non-empty
# (F/XYZ) station pool, all its stations, all their channels. Shared by the
# single-event `run()` entry point and the 117-event `run_all()` sweep.
# ---------------------------------------------------------------------------

def _rank_stations_for_event(cfg: "common.GroupConfig", event, channel_type: str, n: int) -> list[tuple[str, float]]:
    """Fresh per-event haversine ranking -- cfg.stations[*]['distance_km'] is
    only valid relative to the group's anchor event (see common.py::
    _discover_stations), wrong for any other event in a multi-event group."""
    wanted = "F" if channel_type == "F" else "XYZF"
    ranked = sorted(
        ((code, common.haversine_km(m["lat"], m["lon"], event.lat, event.lon))
         for code, m in cfg.stations.items() if m["reported"] == wanted),
        key=lambda x: x[1],
    )
    return ranked[:n]


def _build_channels(df: pd.DataFrame) -> dict[str, pd.Series]:
    if "F" in df.columns:
        return {"F": df["F"]}
    return {"X": df["X"], "Y": df["Y"], "Z": df["Z"], "H": np.sqrt(df["X"] ** 2 + df["Y"] ** 2)}


def _iter_test_payloads(result: dict):
    for m, payload in result.get("step", {}).items():
        yield "step", int(m), payload
    if "spike" in result:
        yield "spike", None, result["spike"]


def _storm_status(cfg: "common.GroupConfig", event_utc: pd.Timestamp) -> dict:
    """Look up the event's UTC date in the group's daily-scale storm_days.csv
    (fetch_space_weather.py output, keyed by UTC file date). A flag only, for
    sensitivity analysis -- nothing is excluded here. The local in-window null
    already absorbs a uniformly storm-elevated noise floor; what it can't
    absorb is a single transient (SSC, substorm onset) landing inside the
    +-SCAN_HALF_SEC search window, which is what the flag lets you check for.
    Flags are None (unknown) if the file is missing or the event falls outside
    the fetched date range, rather than silently reading as "quiet"."""
    unknown = {"is_storm_day": None, "is_storm_onset": None, "storm_flag_confidence": None}
    csv_path = cfg.interim_dir / "storm_days.csv"
    summary_path = cfg.interim_dir / "storm_days_summary.json"
    if not csv_path.exists() or not summary_path.exists():
        return unknown
    summary = json.loads(summary_path.read_text())
    lo, hi = summary["date_range"]
    if not (pd.Timestamp(lo) <= event_utc.normalize() <= pd.Timestamp(hi)):
        return unknown
    storm_days = pd.read_csv(csv_path, dtype={"date": str})
    row = storm_days[storm_days.date == event_utc.strftime("%Y%m%d")]
    return {
        "is_storm_day": bool(row.is_storm_or_recovery.iloc[0]) if len(row) else False,
        "is_storm_onset": bool(row.is_storm_onset.iloc[0]) if len(row) else False,
        "storm_flag_confidence": summary.get("confidence"),
    }


def process_event(cfg: "common.GroupConfig", group, event,
                   run_injection: bool = False) -> tuple[dict, list[dict], list[dict]]:
    """Run every (channel_type, station, channel) test for one Event.
    exclude_centers covers every event in `group` (not just this one) --
    see module docstring's "Extending to all 117 events" note on why."""
    event_utc = pd.Timestamp(event.time_utc)
    # folder_events, not group.events: G2/G3 and G6/G7/G8 share a raw-data folder, and another
    # group's real event in the same data must stay out of this group's null draws too.
    exclude_centers = [pd.Timestamp(e.time_utc) for e in folder_events(group.group_id)]
    scan_half_sec = _effective_half_sec(SCAN_HALF_SEC, event_utc, exclude_centers)
    storm = _storm_status(cfg, event_utc)

    event_result = {
        "group": cfg.group_id, "event_date": event.date, "event_time_local": event.time_local,
        "event_time_utc": event.time_utc, "magnitude": f"{event.magnitude_type}{event.magnitude}",
        "coord_confidence": event.coord_confidence, "anchor": event.anchor,
        **storm,
        "seed": SEED, "n_null": N_NULL, "step_windows_sec": list(STEP_WINDOWS_SEC),
        "scan_half_sec": scan_half_sec, "scan_half_sec_target": SCAN_HALF_SEC,
        "exclusion_buffer_sec": EXCLUSION_BUFFER_SEC,
        "results": [],
    }
    summary_rows: list[dict] = []
    injection_runs: list[dict] = []
    buffer_sec = REQUIRED_MARGIN_SEC + DETREND_WINDOW_SEC // 2 + 60  # generous margin for detrend edges too

    for channel_type in ("F", "XYZ"):
        stations = _rank_stations_for_event(cfg, event, channel_type, N_STATIONS)
        for station, distance_km in stations:
            df, missing, wanted = _load_station_days(cfg.gdms_dir, station, event_utc, buffer_sec)
            if df is None:
                event_result["results"].append({
                    "station": station, "channel_type": channel_type,
                    "distance_km": round(distance_km, 1),
                    "data_status": "no_data", "missing_dates": missing,
                })
                continue

            for ch_label, series in _build_channels(df).items():
                rng = keyed_rng(cfg.group_id, event.date, channel_type, station, ch_label)
                result = _test_channel(series, event_utc, exclude_centers, rng, ch_label, scan_half_sec)
                result["station"] = station
                result["channel_type"] = channel_type
                result["distance_km"] = round(distance_km, 1)
                result["missing_source_dates"] = missing
                event_result["results"].append(result)

                for stat_type, window, payload in _iter_test_payloads(result):
                    summary_rows.append({
                        "group": cfg.group_id, "event_date": event.date,
                        "magnitude": f"{event.magnitude_type}{event.magnitude}",
                        "coord_confidence": event.coord_confidence,
                        "is_storm_day": storm["is_storm_day"], "is_storm_onset": storm["is_storm_onset"],
                        "channel_type": channel_type, "station": station,
                        "distance_km": round(distance_km, 1), "channel": ch_label,
                        "statistic_type": stat_type, "window_sec": window,
                        "data_status": result.get("data_status"),
                        "obs_value_nt": payload.get("obs_value_nt"),
                        "obs_lag_sec": payload.get("obs_lag_sec"),
                        "p_value": payload.get("p_value"),
                        "n_null_samples": payload.get("n_null_samples"),
                        "approx_min_detectable_nt": payload.get("approx_min_detectable_nt"),
                    })

                if run_injection and result.get("data_status") == "ok":
                    injection_runs.append(
                        run_injection_power_curve(series, exclude_centers,
                                                  keyed_rng(cfg.group_id, event.date, channel_type, station, ch_label, "injection"),
                                                  f"{station}:{ch_label}")
                    )

    n_tests = len(summary_rows)
    event_result["multiple_comparisons_context"] = {
        "n_total_tests": n_tests,
        "bonferroni_alpha_at_p05": round(0.05 / n_tests, 6) if n_tests else None,
        "note": "No correction applied per-event -- see all_events_run_summary.json for the "
                "dataset-wide test count when running --all, and the module docstring.",
    }
    return event_result, summary_rows, injection_runs


# ---------------------------------------------------------------------------
# Entry points
# ---------------------------------------------------------------------------

def run(group_id: str = GROUP_ID, run_injection: bool = True) -> dict:
    """Single-event smoke-test / regression-check entry point: tests just
    `group_id`'s anchor event (default G10) against every non-empty station
    pool. This is what `coseismic_step_analysis.py` (no flags) runs."""
    cfg = common.load_group_config(group_id)
    group = get_group(group_id)
    event = group.anchor_event

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "events").mkdir(exist_ok=True)

    event_result, summary_rows, injection_runs = process_event(cfg, group, event, run_injection=run_injection)

    out_path = OUT_DIR / "events" / f"{group_id}__{event.date}.json"
    out_path.write_text(json.dumps(event_result, indent=2))
    print(f"wrote {out_path} ({len(summary_rows)} tests)", file=sys.stderr)

    if run_injection:
        inj_path = OUT_DIR / "injection_power_curve.json"
        inj_path.write_text(json.dumps({"group": group_id, "event_date": event.date, "runs": injection_runs},
                                        indent=2))
        print(f"wrote {inj_path} ({len(injection_runs)} channel runs)", file=sys.stderr)

    pd.DataFrame(summary_rows).to_csv(OUT_DIR / "summary.csv", index=False)
    print(f"wrote summary.csv ({len(summary_rows)} rows)", file=sys.stderr)
    return event_result


def _storm_sensitivity(items: list[dict]) -> dict:
    """Split run_all's per-event rollup by storm flag and compare the
    fraction of tests with p<0.05 in each subset. If storms were driving the
    detections, the storm subset's rate would sit clearly above the quiet
    subset's; similar rates mean the result isn't storm-sensitive."""
    def _rollup(subset: list[dict]) -> dict:
        n_tests = sum(i["n_tests"] for i in subset)
        n_sig = sum(i["n_significant_p_lt_05"] for i in subset)
        return {
            "n_events": len(subset),
            "n_tests": n_tests,
            "n_significant_p_lt_05": n_sig,
            "frac_tests_p_lt_05": round(n_sig / n_tests, 4) if n_tests else None,
            "n_events_with_any_p_lt_05": sum(1 for i in subset if i["n_significant_p_lt_05"]),
        }

    out = {"all": _rollup(items)}
    for flag in ("is_storm_day", "is_storm_onset"):
        out[f"{flag}=True"] = _rollup([i for i in items if i[flag] is True])
        out[f"{flag}=False"] = _rollup([i for i in items if i[flag] is False])
    out["storm_flag_unknown"] = _rollup([i for i in items if i["is_storm_day"] is None])
    return out


def run_all(run_injection: bool = False) -> dict:
    """Sweep every event in events.py's 117-event registry (all 20 groups,
    all events per group -- not just the 20 anchors). Writes one JSON per
    event under events/, a combined summary.csv across all events, and
    all_events_run_summary.json (the run_all_groups.sh-style per-item
    status rollup, so nothing is silently skipped)."""
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "events").mkdir(exist_ok=True)

    all_summary_rows: list[dict] = []
    all_injection_runs: list[dict] = []
    run_summary = {"seed": SEED, "n_null": N_NULL, "events": []}

    for group_id, group in GROUPS.items():
        cfg = common.load_group_config(group_id)
        for event in group.events:
            event_result, summary_rows, injection_runs = process_event(
                cfg, group, event, run_injection=run_injection)

            out_path = OUT_DIR / "events" / f"{group_id}__{event.date}.json"
            out_path.write_text(json.dumps(event_result, indent=2))
            all_summary_rows.extend(summary_rows)
            all_injection_runs.extend(injection_runs)

            p_values = [r["p_value"] for r in summary_rows if r["p_value"] is not None]
            data_statuses = sorted({r["data_status"] for r in event_result["results"]})
            n_stations = len({(r["station"], r.get("channel_type")) for r in event_result["results"]
                               if "station" in r})
            item = {
                "group": group_id, "event_date": event.date, "anchor": event.anchor,
                "magnitude": f"{event.magnitude_type}{event.magnitude}",
                "coord_confidence": event.coord_confidence,
                "is_storm_day": event_result["is_storm_day"], "is_storm_onset": event_result["is_storm_onset"],
                "n_stations_tested": n_stations, "n_tests": len(summary_rows),
                "n_significant_p_lt_05": sum(1 for p in p_values if p < 0.05),
                "min_p": min(p_values) if p_values else None,
                "data_statuses_seen": data_statuses,
            }
            run_summary["events"].append(item)
            print(f"[{group_id} {event.date}] {item['n_tests']} tests, min_p={item['min_p']}, "
                  f"statuses={data_statuses}", file=sys.stderr)

    run_summary["n_events_total"] = len(run_summary["events"])
    run_summary["n_total_tests"] = len(all_summary_rows)
    run_summary["bonferroni_alpha_at_p05"] = (
        round(0.05 / len(all_summary_rows), 8) if all_summary_rows else None
    )
    run_summary["storm_sensitivity"] = _storm_sensitivity(run_summary["events"])

    (OUT_DIR / "all_events_run_summary.json").write_text(json.dumps(run_summary, indent=2))
    pd.DataFrame(all_summary_rows).to_csv(OUT_DIR / "summary.csv", index=False)
    if run_injection:
        (OUT_DIR / "injection_power_curve.json").write_text(
            json.dumps({"runs": all_injection_runs}, indent=2))

    print(f"[run_all] {run_summary['n_events_total']} events, {run_summary['n_total_tests']} total tests "
          f"-> {OUT_DIR}", file=sys.stderr)
    return run_summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--self-test", action="store_true", help="run only the synthetic sanity check")
    ap.add_argument("--no-injection", action="store_true", help="skip the positive-control injection test")
    ap.add_argument("--all", action="store_true",
                     help="run every event in events.py's 117-event registry instead of just G10's anchor")
    args = ap.parse_args()

    if args.self_test:
        passed = self_test()
        sys.exit(0 if passed else 1)

    if not self_test():
        print("[main] synthetic self-test FAILED -- aborting before touching real data", file=sys.stderr)
        sys.exit(1)

    if args.all:
        if not args.no_injection:
            print("[main] --all always runs with injection testing disabled (already validated once on "
                  "G10, see injection_power_curve.json / plan file) -- ignoring lack of --no-injection",
                  file=sys.stderr)
        run_all(run_injection=False)
    else:
        run(run_injection=not args.no_injection)
