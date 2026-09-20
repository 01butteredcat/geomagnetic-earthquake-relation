"""Cross-event superposed-epoch stacking of the coseismic step/spike
statistics -- the task `coseismic_step_analysis.py`'s own docstring
(see its "Extending to all 49 events" section) explicitly named as the
natural next step but deliberately did not implement:

    "Cross-event stacking (normalizing each event's window via
    stat_utils.mad_zscore on its own off-event statistic values, then
    averaging across events the way superposed_epoch_analysis.py stacks
    day-scale series, with a null band built from random per-event
    reference times) is the natural next step if the 49-event run shows
    many individually-underpowered-but-suggestive events, but is
    deliberately not implemented in this pass."

`coseismic_step_analysis.py` found 19/31 events individually significant
near their origin second, but only ever reports each event's *extremum*
statistic within a +-120s search window -- never the full window shape, so
there was no way to see whether all these events actually look alike (a
consistent step/spike shape that a real shared mechanism would produce) or
just happen to individually cross a p<0.05 threshold in unrelated ways. This
script answers that by extracting each event's full window *profile* (not
just its extremum), z-normalizing it against that event's own off-event
noise floor (so events with very different absolute noise floors -- e.g. a
quiet scalar F station vs a noisy near-field XYZ station -- are on a
comparable footing before averaging), and stacking across events -- exactly
mirroring `superposed_epoch_analysis.py::run_band`'s day-scale bootstrap-CI
+ null-band design, just at 1-second instead of 1-day resolution.

## Why a new file, not an extension of coseismic_step_analysis.py

That script's own output (`data/interim/coseismic_step_analysis/events/*.json`,
`summary.csv`, `all_events_run_summary.json`) is already referenced by a
published report. This script imports its detrend/statistic/null-draw
helpers directly (not a duplicated copy -- unlike that script's own
deliberate duplication of `ulf_analysis.py::_detrend`, which exists because
the two scripts serve genuinely different purposes; this one is a direct,
same-author continuation of coseismic_step_analysis.py itself, reusing
machinery already validated on real data at 49-event scale) but writes its
own new output tree, so nothing about the existing artifacts changes.

## What gets stacked

Four statistics, kept separate (not merged) -- step10/step30/step90/spike,
mirroring how superposed_epoch_analysis.py treats pc3/pc4 as independent
bands rather than combining them. Channel is reduced to one canonical,
cross-station-comparable quantity per pool: H=sqrt(X^2+Y^2) for the XYZ
pool, F for the F-only pool (X/Y individually are station-orientation-
dependent and not directly comparable across stations). Station selection
is per-event (not the group's anchor-relative cfg.stations distances, which
are wrong for non-anchor events -- see coseismic_step_analysis.py's own
_rank_stations_for_event and events.py's G14 note), at two tiers: the
nearest and 2nd-nearest station for that specific event ("near"/"far",
matching common.py's existing near/far vocabulary -- not a claim that "far"
is actually distant, just whichever station ranks 2nd for this event). Both
tiers are stacked so a distance-sensitivity comparison is possible: if the
stacked signal is much weaker at the 2nd-nearest station, that argues for a
spatially localized (near-field shaking-coupled) origin; if similar, that
argues against a cleanly localized signal.

The stacking half-window is STACK_HALF_SEC=300 (wider than coseismic_step_
analysis.py's SCAN_HALF_SEC=120), specifically so the stacked shape can show
whether the average signal *reverts* (G10's "noise burst then recovery"
pattern) or *stays elevated* (G9 csg's "persistent, non-recovering offset"
pattern) well past the detection window -- that reversion-vs-persistence
question is exactly what motivates this script, and 120s doesn't leave much
room to see it.

## Normalization (per event x station x channel x statistic)

1. Detrend + compute the statistic series `d` over the full loaded buffer
   (reuses coseismic_step_analysis.py's _detrend/_step_statistic/
   _spike_statistic unchanged).
2. Off-event baseline: median/MAD (stat_utils.mad_zscore) of `d` restricted
   to samples outside +-EXCLUSION_BUFFER_SEC of *every* real event in that
   event's group (not just the event being tested -- some groups have
   events close enough in time, e.g. G9's two events 17h apart, that a
   single-event exclusion could let a different real earthquake contaminate
   the "off-event" pool).
3. Extract the [-STACK_HALF_SEC, +STACK_HALF_SEC] window profile and
   z-normalize it against that (median, MAD) -- this z-profile is what gets
   stacked.

## Bootstrap CI + null band

Direct 1-second-scale port of superposed_epoch_analysis.py::run_band:
- Bootstrap CI90: N_BOOTSTRAP=2000 resamples of *events* with replacement,
  5th/95th percentile of the resampled stack mean at each lag.
- Null band: N_NULL_STACK=1000 realizations; each draws ONE random
  reference time per contributing event (coseismic_step_analysis.py's own
  `_draw_null_centers`, respecting the same exclusion buffer), re-extracts
  that event's window profile there, re-normalizes with the *same*
  precomputed (median, MAD) (not recomputed -- the null band asks "does a
  random second on this event's own noise floor look like this", so the
  baseline must stay fixed), and stacks. Percentiles 5/50/95 across the
  1000 realizations at each lag give a pointwise band for visualization.
  Significance of the real stack's peak is judged NOT against that pointwise
  band (self_test() showed this produces false positives -- the real peak is
  chosen via argmax over 601 lags, so comparing it to one lag's pointwise
  band is a search-vs-no-search mismatch) but against the distribution of
  each null realization's OWN peak |z| across all lags -- the same
  "identical search procedure for observed and null" design coseismic_step_
  analysis.py itself uses for its own p-values.

Fixed SEED=20260805 (reused from the rest of this codebase) for
reproducibility.

Usage:
  coseismic_stacking_analysis.py --self-test          # synthetic sanity check only
  coseismic_stacking_analysis.py --group G9 --group G10   # smoke test, 2 groups
  coseismic_stacking_analysis.py --all                # every group in events.py
"""
from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
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
    MISSING_FRACTION_THRESHOLD,
    SEED,
    _build_channels,
    _detrend,
    _draw_null_centers,
    _load_station_days,
    _pos,
    _rank_stations_for_event,
    _spike_statistic,
    _step_statistic,
)

STACK_HALF_SEC = 300
N_BOOTSTRAP = 2000
N_NULL_STACK = 1000
STATION_TIERS = ("near", "far")  # rank 0 / rank 1 by per-event haversine distance -- see module docstring
STAT_NAMES = ("step10", "step30", "step90", "spike")
REQUIRED_MARGIN_SEC = STACK_HALF_SEC + EXCLUSION_BUFFER_SEC  # 900s: minimum data support needed at either edge
BUFFER_SEC = REQUIRED_MARGIN_SEC + DETREND_WINDOW_SEC // 2 + 60  # 2760s -- mirrors coseismic_step_analysis.py's own buffer_sec formula
MIN_OFF_EVENT_SAMPLES = 200  # below this, the off-event median/MAD is too noisy an estimate to normalize against

OUT_DIR = common.PROJECT_DIR / "data" / "interim" / "coseismic_stacking_analysis"


# ---------------------------------------------------------------------------
# Statistic / windowing helpers
# ---------------------------------------------------------------------------

def _stat_array(d: np.ndarray, stat_name: str) -> np.ndarray:
    if stat_name == "spike":
        return _spike_statistic(d)
    return _step_statistic(d, int(stat_name.replace("step", "")))


def _window_profile(arr: np.ndarray, idx: pd.DatetimeIndex, center: pd.Timestamp, half_sec: int) -> np.ndarray:
    """Positional slice of `arr` centered on `center` (same O(window)
    positional-index approach as coseismic_step_analysis.py::
    _extremum_in_window, just keeping the whole profile instead of only its
    extremum). NaN-padded if the requested window runs past the array's
    edges -- shouldn't happen given BUFFER_SEC's margin, but handled
    defensively rather than assumed."""
    p = _pos(idx, center)
    n = len(arr)
    out = np.full(2 * half_sec + 1, np.nan)
    lo, hi = p - half_sec, p + half_sec + 1
    src_lo, src_hi = max(0, lo), min(n, hi)
    if src_hi <= src_lo:
        return out
    dst_lo = src_lo - lo
    out[dst_lo:dst_lo + (src_hi - src_lo)] = arr[src_lo:src_hi]
    return out


def _off_event_baseline(arr: np.ndarray, idx: pd.DatetimeIndex,
                         exclude_centers: list[pd.Timestamp], buffer_sec: int) -> tuple[float, float] | None:
    """median/MAD of `arr` excluding every sample within buffer_sec of any
    real event in exclude_centers -- the off-event noise floor each event's
    on-event window gets normalized against. None if too few clean samples
    survive the exclusion to trust the estimate."""
    mask = np.ones(len(arr), dtype=bool)
    for c in exclude_centers:
        p = _pos(idx, c)
        lo, hi = max(0, p - buffer_sec), min(len(arr), p + buffer_sec + 1)
        mask[lo:hi] = False
    vals = arr[mask]
    vals = vals[~np.isnan(vals)]
    if len(vals) < MIN_OFF_EVENT_SAMPLES:
        return None
    _, med, mad = mad_zscore(vals)
    if not np.isfinite(mad) or mad <= 1e-9:
        return None
    return med, mad


# ---------------------------------------------------------------------------
# EventSeries: decouples the stacking core from real-data loading, so
# self_test() can exercise the exact same stack_series() code on synthetic
# data without touching any .sec files.
# ---------------------------------------------------------------------------

@dataclass
class EventSeries:
    group_id: str
    event_date: str
    anchor: bool
    magnitude: str
    station: str
    distance_km: float
    idx: pd.DatetimeIndex
    stat_arrays: dict           # stat_name -> full detrended-statistic series (np.ndarray)
    baselines: dict             # stat_name -> (median, mad) off-event
    event_utc: pd.Timestamp
    valid_lo: pd.Timestamp
    valid_hi: pd.Timestamp
    exclude_centers: list       # list[pd.Timestamp], every real event in this EventSeries' group


def _build_event_series(cfg: "common.GroupConfig", group, event, station: str, distance_km: float,
                         raw: pd.Series) -> EventSeries | None:
    idx = raw.index
    event_utc = pd.Timestamp(event.time_utc)
    exclude_centers = [pd.Timestamp(e.time_utc) for e in folder_events(group.group_id)]

    if idx.min() + pd.Timedelta(seconds=REQUIRED_MARGIN_SEC) > event_utc or \
       idx.max() - pd.Timedelta(seconds=REQUIRED_MARGIN_SEC) < event_utc:
        return None
    nan_frac = float(
        raw[(idx >= event_utc - pd.Timedelta(seconds=STACK_HALF_SEC)) &
            (idx <= event_utc + pd.Timedelta(seconds=STACK_HALF_SEC))].isna().mean()
    )
    if nan_frac > MISSING_FRACTION_THRESHOLD:
        return None

    d = _detrend(raw)
    stat_arrays, baselines = {}, {}
    for stat_name in STAT_NAMES:
        arr = _stat_array(d, stat_name)
        base = _off_event_baseline(arr, idx, exclude_centers, EXCLUSION_BUFFER_SEC)
        if base is None:
            continue
        stat_arrays[stat_name] = arr
        baselines[stat_name] = base
    if not stat_arrays:
        return None

    return EventSeries(
        group_id=cfg.group_id, event_date=event.date, anchor=event.anchor,
        magnitude=f"{event.magnitude_type}{event.magnitude}",
        station=station, distance_km=round(distance_km, 1),
        idx=idx, stat_arrays=stat_arrays, baselines=baselines,
        event_utc=event_utc,
        valid_lo=idx.min() + pd.Timedelta(seconds=STACK_HALF_SEC),
        valid_hi=idx.max() - pd.Timedelta(seconds=STACK_HALF_SEC),
        exclude_centers=exclude_centers,
    )


def load_all_event_series(group_ids: tuple[str, ...]) -> dict[tuple[str, str, str], EventSeries]:
    """Returns {(channel_type_pool, tier, "<group_id>__<event.date>"): EventSeries}.
    Every event in every requested group is attempted at both channel-type
    pools ("F"/"XYZ", whichever the group actually has stations for -- same
    dual-pool loop as coseismic_step_analysis.py::process_event) and both
    station tiers. All 49 events across all 20 groups are eligible -- this
    only needs the existing .sec geomagnetic data, no seismometer dependency."""
    out: dict[tuple[str, str, str], EventSeries] = {}
    for group_id in group_ids:
        cfg = common.load_group_config(group_id)
        group = get_group(group_id)
        for event in group.events:
            event_utc = pd.Timestamp(event.time_utc)
            for channel_type in ("F", "XYZ"):
                stations = _rank_stations_for_event(cfg, event, channel_type, len(STATION_TIERS))
                for rank, tier in enumerate(STATION_TIERS):
                    if rank >= len(stations):
                        continue
                    station, distance_km = stations[rank]
                    df, missing, wanted = _load_station_days(cfg.gdms_dir, station, event_utc, BUFFER_SEC)
                    if df is None:
                        continue
                    channels = _build_channels(df)
                    ch_label = "H" if channel_type == "XYZ" else "F"
                    if ch_label not in channels:
                        continue
                    es = _build_event_series(cfg, group, event, station, distance_km, channels[ch_label])
                    if es is None:
                        continue
                    out[(channel_type, tier, f"{group_id}__{event.date}")] = es
    return out


# ---------------------------------------------------------------------------
# Stacking core (works purely on EventSeries objects -- exercised by both
# real data via run_group_ids() and synthetic data via self_test()).
# ---------------------------------------------------------------------------

def stack_series(events: list[EventSeries], stat_name: str, rng: np.random.Generator) -> dict:
    lags = np.arange(-STACK_HALF_SEC, STACK_HALF_SEC + 1)
    matrix: list[np.ndarray] = []
    usable: list[EventSeries] = []
    used: list[dict] = []
    excluded: list[dict] = []

    for es in events:
        if stat_name not in es.stat_arrays:
            excluded.append({"group": es.group_id, "date": es.event_date, "reason": "no_baseline_for_stat"})
            continue
        med, mad = es.baselines[stat_name]
        profile = _window_profile(es.stat_arrays[stat_name], es.idx, es.event_utc, STACK_HALF_SEC)
        if np.all(np.isnan(profile)):
            excluded.append({"group": es.group_id, "date": es.event_date, "reason": "empty_window"})
            continue
        matrix.append((profile - med) / mad)
        usable.append(es)
        used.append({"group": es.group_id, "date": es.event_date, "anchor": es.anchor,
                      "magnitude": es.magnitude, "station": es.station, "distance_km": es.distance_km})

    if not matrix:
        return {"stat": stat_name, "error": "no events with usable data", "events_excluded": excluded}

    M = np.array(matrix)
    n_events = M.shape[0]
    stack_mean = np.nanmean(M, axis=0)
    stack_median = np.nanmedian(M, axis=0)
    n_contributing = np.sum(~np.isnan(M), axis=0)

    # Bootstrap CI on the real stack (resample events with replacement)
    boot = np.empty((N_BOOTSTRAP, M.shape[1]))
    for b in range(N_BOOTSTRAP):
        sample = rng.integers(0, n_events, size=n_events)
        boot[b] = np.nanmean(M[sample], axis=0)
    ci_lo = np.nanpercentile(boot, 5, axis=0)
    ci_hi = np.nanpercentile(boot, 95, axis=0)

    # Null band: repeat the whole stacking procedure with a random reference
    # time per event instead of its real origin time.
    null_stacks = np.empty((N_NULL_STACK, M.shape[1]))
    for r in range(N_NULL_STACK):
        null_matrix = []
        for es in usable:
            centers = _draw_null_centers(rng, es.valid_lo, es.valid_hi, es.exclude_centers,
                                          EXCLUSION_BUFFER_SEC, 1)
            if not centers:
                continue
            med, mad = es.baselines[stat_name]
            profile = _window_profile(es.stat_arrays[stat_name], es.idx, centers[0], STACK_HALF_SEC)
            if np.all(np.isnan(profile)):
                continue
            null_matrix.append((profile - med) / mad)
        null_stacks[r] = np.nanmean(np.array(null_matrix), axis=0) if null_matrix else np.nan

    null_p5 = np.nanpercentile(null_stacks, 5, axis=0)
    null_p50 = np.nanpercentile(null_stacks, 50, axis=0)
    null_p95 = np.nanpercentile(null_stacks, 95, axis=0)

    peak_i = int(np.nanargmax(np.abs(stack_mean)))
    peak_val = float(stack_mean[peak_i])

    # Significance of the peak: NOT a pointwise comparison against
    # null_p5/p95 at peak_i -- the real peak was chosen via argmax over all
    # 601 lags, so comparing it to a single lag's pointwise band is an
    # apples-to-oranges "search vs no-search" comparison that inflates false
    # positives (confirmed by self_test()'s noise-only check: a pure-noise
    # stack's peak routinely clears the pointwise band just from having 601
    # chances to do so). The correct, matched comparison -- mirroring
    # coseismic_step_analysis.py's own "same search procedure for observed
    # and null" design -- is against the distribution of each NULL
    # realization's *own* peak |z| across all lags.
    null_peak_abs = np.nanmax(np.abs(null_stacks), axis=1)
    null_peak_abs = null_peak_abs[~np.isnan(null_peak_abs)]
    if len(null_peak_abs):
        p_value = (1 + int(np.sum(null_peak_abs >= abs(peak_val)))) / (len(null_peak_abs) + 1)
        null_peak_p95 = float(np.percentile(null_peak_abs, 95))
    else:
        p_value = None
        null_peak_p95 = None
    outside = bool(p_value is not None and p_value < 0.05)

    return {
        "stat": stat_name,
        "half_sec": STACK_HALF_SEC,
        "seed": SEED, "n_bootstrap": N_BOOTSTRAP, "n_null_stacks": N_NULL_STACK,
        "n_events_used": n_events,
        "events_used": used,
        "events_excluded": excluded,
        "lags_sec": lags.tolist(),
        "stack_mean": [None if np.isnan(v) else round(float(v), 4) for v in stack_mean],
        "stack_median": [None if np.isnan(v) else round(float(v), 4) for v in stack_median],
        "n_contributing_per_lag": n_contributing.tolist(),
        "bootstrap_ci90_lo": [round(float(v), 4) for v in ci_lo],
        "bootstrap_ci90_hi": [round(float(v), 4) for v in ci_hi],
        "null_band_p5": [round(float(v), 4) for v in null_p5],
        "null_band_p50": [round(float(v), 4) for v in null_p50],
        "null_band_p95": [round(float(v), 4) for v in null_p95],
        "peak_abs_z": round(abs(peak_val), 4),
        "peak_lag_sec": int(lags[peak_i]),
        "null_peak_abs_z_p95": None if null_peak_p95 is None else round(null_peak_p95, 4),
        "p_value": p_value,
        "outside_null_band_at_peak": outside,
    }


# ---------------------------------------------------------------------------
# Synthetic self-test: exercises stack_series() directly (no real data),
# confirming (a) a stack of events with a known injected step recovers a
# peak near lag=0 that lands outside its own null band, and (b) a stack of
# pure-noise-only events does NOT show a spurious peak outside its null band.
# ---------------------------------------------------------------------------

def self_test() -> bool:
    rng = np.random.default_rng(SEED)
    n_synth = 20
    n_samples = 2 * BUFFER_SEC + 1
    idx = pd.date_range("2024-01-01", periods=n_samples, freq="s")
    amplitude = 15.0
    center_i = BUFFER_SEC

    def make_series(inject: bool) -> EventSeries:
        noise = rng.normal(0, 1.0, size=n_samples)
        if inject:
            jitter = int(rng.integers(-5, 6))
            noise[center_i + jitter:center_i + jitter + 20] += amplitude
        event_utc = idx[center_i]
        arr = _step_statistic(noise, 30)
        base = _off_event_baseline(arr, idx, [event_utc], EXCLUSION_BUFFER_SEC)
        assert base is not None, "self-test off-event baseline computation failed"
        return EventSeries(
            group_id="SYN", event_date=f"synthetic-{int(rng.integers(0, 10**6))}", anchor=True, magnitude="M0",
            station="SYN", distance_km=0.0, idx=idx, stat_arrays={"step30": arr}, baselines={"step30": base},
            event_utc=event_utc,
            valid_lo=idx[0] + pd.Timedelta(seconds=STACK_HALF_SEC),
            valid_hi=idx[-1] - pd.Timedelta(seconds=STACK_HALF_SEC),
            exclude_centers=[event_utc],
        )

    signal_events = [make_series(inject=True) for _ in range(n_synth)]
    noise_events = [make_series(inject=False) for _ in range(n_synth)]

    rng2 = np.random.default_rng(SEED + 1)
    signal_result = stack_series(signal_events, "step30", rng2)
    noise_result = stack_series(noise_events, "step30", rng2)

    ok = True
    near_zero = abs(signal_result["peak_lag_sec"]) <= 10
    outside = signal_result["outside_null_band_at_peak"]
    status1 = "PASS" if (near_zero and outside) else "FAIL"
    print(f"[self-test] injected-signal stack: peak_lag={signal_result['peak_lag_sec']}s "
          f"peak_abs_z={signal_result['peak_abs_z']} outside_null={outside}  {status1}")
    ok = ok and near_zero and outside

    noise_outside = noise_result["outside_null_band_at_peak"]
    status2 = "PASS" if not noise_outside else "FAIL"
    print(f"[self-test] noise-only stack (false-positive check): peak_abs_z={noise_result['peak_abs_z']} "
          f"outside_null={noise_outside}  {status2}")
    ok = ok and not noise_outside

    return ok


# ---------------------------------------------------------------------------
# Real-data orchestration
# ---------------------------------------------------------------------------

def run_group_ids(group_ids: tuple[str, ...]) -> dict:
    all_series = load_all_event_series(group_ids)
    rng = np.random.default_rng(SEED)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "combos").mkdir(exist_ok=True)

    channel_labels = {"XYZ": "H", "F": "F"}
    summary_rows: list[dict] = []
    run_summary = {"seed": SEED, "group_ids": list(group_ids), "combos": []}

    for channel_type, ch_label in channel_labels.items():
        for tier in STATION_TIERS:
            events = [es for (ct, t, _key), es in all_series.items() if ct == channel_type and t == tier]
            if not events:
                continue
            for stat_name in STAT_NAMES:
                result = stack_series(events, stat_name, rng)
                combo_id = f"{tier}__{ch_label}__{stat_name}"
                result["combo_id"] = combo_id
                result["station_tier"] = tier
                result["channel_type_pool"] = channel_type
                result["channel"] = ch_label

                out_path = OUT_DIR / "combos" / f"stack__{combo_id}.json"
                out_path.write_text(json.dumps(result, indent=2))

                if "error" in result:
                    print(f"[{combo_id}] {result['error']}", file=sys.stderr)
                    run_summary["combos"].append({"combo_id": combo_id, "status": "error", "error": result["error"]})
                    continue

                print(f"[{combo_id}] n_events={result['n_events_used']} "
                      f"peak_abs_z={result['peak_abs_z']} lag={result['peak_lag_sec']}s "
                      f"outside_null={result['outside_null_band_at_peak']}", file=sys.stderr)
                run_summary["combos"].append({
                    "combo_id": combo_id, "status": "ok",
                    "n_events_used": result["n_events_used"],
                    "peak_abs_z": result["peak_abs_z"], "peak_lag_sec": result["peak_lag_sec"],
                    "outside_null_band_at_peak": result["outside_null_band_at_peak"],
                })
                summary_rows.append({
                    "combo_id": combo_id, "station_tier": tier, "channel_type_pool": channel_type,
                    "channel": ch_label, "statistic": stat_name,
                    "n_events_used": result["n_events_used"],
                    "peak_abs_z": result["peak_abs_z"], "peak_lag_sec": result["peak_lag_sec"],
                    "outside_null_band_at_peak": result["outside_null_band_at_peak"],
                })

    pd.DataFrame(summary_rows).to_csv(OUT_DIR / "stack_summary.csv", index=False)
    (OUT_DIR / "all_stacks_run_summary.json").write_text(json.dumps(run_summary, indent=2))
    print(f"[run] {len(summary_rows)} combos -> {OUT_DIR}", file=sys.stderr)
    return run_summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--self-test", action="store_true", help="run only the synthetic sanity check")
    ap.add_argument("--group", action="append", dest="groups",
                     help="restrict to this group ID (repeatable); default is every group in events.py")
    ap.add_argument("--all", action="store_true", help="explicit alias for the default (no --group filter)")
    args = ap.parse_args()

    if args.self_test:
        sys.exit(0 if self_test() else 1)

    if not self_test():
        print("[main] synthetic self-test FAILED -- aborting before touching real data", file=sys.stderr)
        sys.exit(1)

    group_ids = tuple(args.groups) if args.groups else tuple(GROUPS.keys())
    run_group_ids(group_ids)
