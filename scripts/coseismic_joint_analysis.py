"""Joint statistical test combining the two, until now independent,
coseismic evidence lines:

  - `coseismic_stacking_analysis.py`: does the geomagnetic 1Hz data show a
    consistent origin-second step/spike signature when stacked across
    events? (says nothing about whether that signature is a real field
    change or shaking-induced instrument noise)
  - `seismometer_comparison.py`: for the subset of events with independent
    seismometer/accelerometer data, does that per-event signature's timing
    line up with real ground shaking (`aligned_with_shaking`, consistent
    with instrument noise) or does it lead/outlast the shaking
    (`leads_shaking` / `persists_after_shaking_ends`, consistent with a
    real mechanism)?

## The question this script asks

Restricted to the events where BOTH lines have data, does the *stacked*
geomagnetic signature actually track the seismometer-derived noise-vs-signal
distinction? If the coseismic signal is mostly shaking noise, events
independently classified `aligned_with_shaking` should stack to a signature
at least as strong as (arguably stronger than, since noise scales with
ground-motion amplitude) the `leads_shaking`/`persists_after_shaking_ends`
events. If a real geophysical mechanism is also present, the
leads/persists group should stack to something visibly stronger and/or
longer-lived than the aligned group.

## Data snapshot (2026-08-16, read directly off
`data/interim/seismometer_comparison/comparison_summary.csv` -- re-derive
rather than hardcode if that file gets regenerated)

27 events have seismic data; 23 also have usable geomag data
(`status == "ok"`; the other 4 -- G2_G3's two events, G4, G16 -- fail with
`no_geomag_data` because their near station's own `.sec` file is missing
for that calendar day, a real historical gap, not a bug). Of those 23:
`aligned_with_shaking`=12, `leads_shaking`=9, `persists_after_shaking_ends`=2.
`persists_after_shaking_ends` alone is too thin (n=2) to stack on its own,
so this script uses a **two-arm** split, not three:
  - `noise_arm`  = alignment_verdict == "aligned_with_shaking"        (12)
  - `signal_arm` = alignment_verdict in ("leads_shaking",
                                          "persists_after_shaking_ends") (11)
This directly reuses `seismometer_comparison.py::alignment_verdict`'s own
framing (its docstring already groups leads/persists together as "favors a
real geophysical mechanism").

(Counts above are from the original 27-event pass. As of 2026-09-24 the arms
are 14 noise / 12 signal = 26 events. The 2026-09-23 registry backfill to 117
events -- 68 M5.0-5.9 non-anchor events from the CWA GDMS export -- added no
arm members, since none of the new events has seismometer data; the only newly
matchable one, G11 2025-01-21b, is `insufficient_data`. It still nudges the
results slightly: the new events widen the off-event exclusion around the 6 arm
events in G11/G12/G13/G19/G20, which shifts their baselines. Rerun 2026-09-24:
far-H deltas moved in the third decimal, min p_tail stayed 0.069, min p_peak
went 0.095 -> 0.093.)

(2026-09-25: after the 80-event seismometer batch, all verdicts give 94 armed
events (48 noise / 46 signal), 38 of them M>=6 (18 / 20). Use --min-mag 6 for
the M>=6-only run; the M5 labels are close to coin flips -- see run_all().
All events: 3 of 32 p-values < 0.05 (far-H spike peak 0.010 and tail 0.048,
both with the NOISE arm stronger; near-H step30 tail 0.024, signal arm
stronger), none surviving Bonferroni. M>=6 only: none below 0.05, min 0.062.)

## Method

1. Build `EventSeries` objects (`coseismic_stacking_analysis.py`'s own
   dataclass + loader machinery, unmodified) for exactly the 23-event
   allowlist via a new loader (`load_event_series_for_events`) that mirrors
   `load_all_event_series` but restricts to explicit (group_id, event_date)
   pairs instead of every event in a group -- added here, not in
   `coseismic_stacking_analysis.py` itself, so that already-published
   script stays untouched.
2. For each (channel_type_pool, station_tier, stat_name) combo (mirrors
   `coseismic_stacking_analysis.py::run_group_ids`'s own combo loop):
   observed Delta = arm difference in TWO statistics (not peak alone,
   because a `persists`-type signal's distinguishing feature is the
   post-shaking TAIL, not necessarily the peak):
     - `delta_peak` = peak_abs_z(signal_arm) - peak_abs_z(noise_arm)
     - `delta_tail` = mean(|stack_mean| for lag>0)(signal_arm) - same(noise_arm)
3. Null distribution: label-permutation test -- shuffle the noise/signal
   arm labels across the same 23 events (keeping arm sizes fixed at 12/11)
   N_PERM=2000 times, recomputing both deltas each time via a lightweight
   stack-mean-only helper (`_lightweight_stack_delta`) that skips
   `stack_series()`'s own bootstrap CI + null band (2000+1000 resamples --
   redoing that on every one of 2000 permutation draws would be off by
   several orders of magnitude too slow). The full `stack_series()` (with
   its bootstrap CI and null band) is still called exactly twice per combo
   -- once per arm -- purely for reporting/plotting the two arms' own stack
   shapes, not for the permutation test itself.
4. p-value = (1 + #{|delta_perm| >= |delta_obs|}) / (N_PERM + 1), same
   one-sided-on-|.| convention as every other resampling test in this
   codebase.

## Honest caveat (stated up front, not after the fact)

N=23 events split 12/11 is an exploratory sample size, not a
high-power one. Results here should be reported as suggestive at most,
never as a confirmatory finding on their own.

Usage:
  coseismic_joint_analysis.py --self-test    # synthetic sanity check only
  coseismic_joint_analysis.py --all          # real data, all 16 combos
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import common  # noqa: E402
from events import GROUPS, get_group  # noqa: E402
from coseismic_step_analysis import (  # noqa: E402
    EXCLUSION_BUFFER_SEC,
    SEED,
    _build_channels,
    _load_station_days,
    _rank_stations_for_event,
    _step_statistic,
)
from coseismic_stacking_analysis import (  # noqa: E402
    BUFFER_SEC,
    STACK_HALF_SEC,
    STAT_NAMES,
    STATION_TIERS,
    EventSeries,
    _build_event_series,
    _off_event_baseline,
    _window_profile,
    stack_series,
)

N_PERM = 2000
NOISE_VERDICTS = ("aligned_with_shaking",)
SIGNAL_VERDICTS = ("leads_shaking", "persists_after_shaking_ends")

OUT_DIR = common.PROJECT_DIR / "data" / "interim" / "coseismic_joint_analysis"
COMPARISON_CSV = common.PROJECT_DIR / "data" / "interim" / "seismometer_comparison" / "comparison_summary.csv"


# ---------------------------------------------------------------------------
# Arm assignment: read seismometer_comparison.py's own output rather than
# recomputing alignment_verdict here -- single source of truth.
# ---------------------------------------------------------------------------

def load_arm_assignment(csv_path: Path = COMPARISON_CSV, verdict_col: str = "alignment_verdict") -> dict[str, str]:
    """Returns {"<group_id>__<event_date>": "noise"|"signal"} for every row
    with status=="ok" and a verdict in one of the two known arms. Rows with
    any other status (e.g. no_geomag_data) or verdict (insufficient_data)
    are simply absent from the returned dict -- not an error, just excluded
    from the joint analysis, exactly as documented in the module docstring."""
    df = pd.read_csv(csv_path)
    arm_of: dict[str, str] = {}
    for _, row in df.iterrows():
        if row.get("status") != "ok":
            continue
        verdict = row.get(verdict_col)
        key = f"{row['group']}__{row['date']}"
        if verdict in NOISE_VERDICTS:
            arm_of[key] = "noise"
        elif verdict in SIGNAL_VERDICTS:
            arm_of[key] = "signal"
    return arm_of


# ---------------------------------------------------------------------------
# EventSeries loader restricted to an explicit event allowlist (new here --
# coseismic_stacking_analysis.py::load_all_event_series is untouched).
# ---------------------------------------------------------------------------

def load_event_series_for_events(event_keys: list[tuple[str, str]]) -> dict[tuple[str, str, str], EventSeries]:
    """Same loading logic as coseismic_stacking_analysis.py::
    load_all_event_series, restricted to an explicit (group_id, event_date)
    allowlist instead of every event in a requested group -- needed because
    the arms are defined per-event (from seismometer_comparison.py's
    per-event alignment_verdict), not per-group."""
    wanted: dict[str, set[str]] = {}
    for group_id, event_date in event_keys:
        wanted.setdefault(group_id, set()).add(event_date)

    out: dict[tuple[str, str, str], EventSeries] = {}
    for group_id, dates in wanted.items():
        cfg = common.load_group_config(group_id)
        group = get_group(group_id)
        for event in group.events:
            if event.date not in dates:
                continue
            event_utc = pd.Timestamp(event.time_utc)
            for channel_type in ("F", "XYZ"):
                stations = _rank_stations_for_event(cfg, event, channel_type, len(STATION_TIERS))
                for rank, tier in enumerate(STATION_TIERS):
                    if rank >= len(stations):
                        continue
                    station, distance_km = stations[rank]
                    df, missing, wanted_dates = _load_station_days(cfg.gdms_dir, station, event_utc, BUFFER_SEC)
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
# Lightweight stack-mean-only path for the permutation loop (skips
# stack_series()'s own bootstrap CI + null band -- see module docstring for
# why: 2000 permutation draws x that machinery would be far too slow).
# ---------------------------------------------------------------------------

def _lightweight_stack_delta(events: list[EventSeries], stat_name: str) -> dict | None:
    lags = np.arange(-STACK_HALF_SEC, STACK_HALF_SEC + 1)
    matrix = []
    for es in events:
        if stat_name not in es.stat_arrays:
            continue
        med, mad = es.baselines[stat_name]
        profile = _window_profile(es.stat_arrays[stat_name], es.idx, es.event_utc, STACK_HALF_SEC)
        if np.all(np.isnan(profile)):
            continue
        matrix.append((profile - med) / mad)
    if not matrix:
        return None
    M = np.array(matrix)
    stack_mean = np.nanmean(M, axis=0)
    peak_abs_z = float(np.nanmax(np.abs(stack_mean)))
    tail_vals = stack_mean[lags > 0]
    tail_vals = tail_vals[~np.isnan(tail_vals)]
    if len(tail_vals) == 0:
        return None
    return {"peak_abs_z": peak_abs_z, "tail_mean_abs_z": float(np.mean(np.abs(tail_vals))), "n_events": len(matrix)}


def _arm_deltas(events: list[EventSeries], labels: np.ndarray, stat_name: str) -> tuple[float, float] | None:
    noise = [es for es, lab in zip(events, labels) if lab == "noise"]
    signal = [es for es, lab in zip(events, labels) if lab == "signal"]
    dn = _lightweight_stack_delta(noise, stat_name)
    ds = _lightweight_stack_delta(signal, stat_name)
    if dn is None or ds is None:
        return None
    return ds["peak_abs_z"] - dn["peak_abs_z"], ds["tail_mean_abs_z"] - dn["tail_mean_abs_z"]


def permutation_test(events: list[EventSeries], arm_labels: np.ndarray, stat_name: str,
                      rng: np.random.Generator, n_perm: int = N_PERM) -> dict:
    obs = _arm_deltas(events, arm_labels, stat_name)
    if obs is None:
        return {"error": "insufficient events in one or both arms for this stat/combo"}
    obs_peak, obs_tail = obs

    perm_peak, perm_tail = [], []
    for _ in range(n_perm):
        shuffled = rng.permutation(arm_labels)
        d = _arm_deltas(events, shuffled, stat_name)
        if d is None:
            continue
        perm_peak.append(d[0])
        perm_tail.append(d[1])
    perm_peak = np.array(perm_peak)
    perm_tail = np.array(perm_tail)

    def _p(obs_val, perm_vals):
        if len(perm_vals) == 0:
            return None
        return (1 + int(np.sum(np.abs(perm_vals) >= abs(obs_val)))) / (len(perm_vals) + 1)

    return {
        "delta_peak_abs_z": round(obs_peak, 4),
        "delta_tail_mean_abs_z": round(obs_tail, 4),
        "p_value_peak": _p(obs_peak, perm_peak),
        "p_value_tail": _p(obs_tail, perm_tail),
        "n_perm_used": len(perm_peak),
    }


# ---------------------------------------------------------------------------
# Synthetic self-test (no real data): confirms the permutation machinery
# itself (a) detects a real injected arm difference and (b) does not
# false-positive when there is none -- same two-sided sanity pattern as
# coseismic_stacking_analysis.py::self_test.
# ---------------------------------------------------------------------------

def _make_synth_series(rng: np.random.Generator, amplitude: float, tag: str) -> EventSeries:
    n_samples = 2 * BUFFER_SEC + 1
    idx = pd.date_range("2024-01-01", periods=n_samples, freq="s")
    center_i = BUFFER_SEC
    noise = rng.normal(0, 1.0, size=n_samples)
    if amplitude > 0:
        jitter = int(rng.integers(-5, 6))
        noise[center_i + jitter:center_i + jitter + 20] += amplitude
    event_utc = idx[center_i]
    arr = _step_statistic(noise, 30)
    base = _off_event_baseline(arr, idx, [event_utc], EXCLUSION_BUFFER_SEC)
    assert base is not None, "self-test off-event baseline computation failed"
    return EventSeries(
        group_id="SYN", event_date=f"{tag}-{int(rng.integers(0, 10**6))}", anchor=True, magnitude="M0",
        station="SYN", distance_km=0.0, idx=idx, stat_arrays={"step30": arr}, baselines={"step30": base},
        event_utc=event_utc,
        valid_lo=idx[0] + pd.Timedelta(seconds=STACK_HALF_SEC),
        valid_hi=idx[-1] - pd.Timedelta(seconds=STACK_HALF_SEC),
        exclude_centers=[event_utc],
    )


def self_test() -> bool:
    rng = np.random.default_rng(SEED)
    n_each = 12
    ok = True

    # Case (a): real arm difference (signal arm gets a stronger injection)
    noise_events = [_make_synth_series(rng, 2.0, "noise") for _ in range(n_each)]
    signal_events = [_make_synth_series(rng, 15.0, "signal") for _ in range(n_each)]
    events_a = noise_events + signal_events
    labels_a = np.array(["noise"] * n_each + ["signal"] * n_each)
    result_a = permutation_test(events_a, labels_a, "step30", rng)
    detected = (result_a.get("p_value_peak") is not None and result_a["p_value_peak"] < 0.05
                and result_a["delta_peak_abs_z"] > 0)
    status_a = "PASS" if detected else "FAIL"
    print(f"[self-test] real-difference case: delta_peak_abs_z={result_a.get('delta_peak_abs_z')} "
          f"p_value_peak={result_a.get('p_value_peak')}  {status_a}")
    ok = ok and detected

    # Case (b): no true arm difference (both arms same amplitude) -- false-positive check
    rng2 = np.random.default_rng(SEED + 1)
    events_b = ([_make_synth_series(rng2, 5.0, "noise") for _ in range(n_each)]
                + [_make_synth_series(rng2, 5.0, "signal") for _ in range(n_each)])
    labels_b = np.array(["noise"] * n_each + ["signal"] * n_each)
    result_b = permutation_test(events_b, labels_b, "step30", rng2)
    no_false_positive = result_b.get("p_value_peak") is not None and result_b["p_value_peak"] >= 0.05
    status_b = "PASS" if no_false_positive else "FAIL"
    print(f"[self-test] no-difference case (false-positive check): "
          f"delta_peak_abs_z={result_b.get('delta_peak_abs_z')} p_value_peak={result_b.get('p_value_peak')}  {status_b}")
    ok = ok and no_false_positive

    return ok


# ---------------------------------------------------------------------------
# Real-data orchestration
# ---------------------------------------------------------------------------

def run_all(min_mag: float | None = None, gated: bool = False) -> dict:
    """min_mag restricts both arms to events.py events of at least that
    magnitude, written to a separate directory (coseismic_joint_analysis_m<min_mag>).
    Worth running alongside the full set since the 2026-09-25 M5 batch: for a
    small event the magnetometer mostly records noise, and a noise peak anywhere
    in the +-180s search window lands before shaking onset about half the time,
    so M5 "leads_shaking" labels are close to coin flips.

    gated uses seismometer_comparison.py's alignment_verdict_gated (a verdict
    only where the event's own step30 anomaly has p < 0.05), written to
    coseismic_joint_analysis_gated[_m<mag>]/. On 2026-09-25 that leaves 7 events
    (5 aligned, 1 leads, 1 persists) -- too few to test; the run exists to make
    that explicit rather than to be read as a result."""
    name = OUT_DIR.name + ("_gated" if gated else "") + ("" if min_mag is None else f"_m{min_mag:g}")
    out_dir = OUT_DIR.with_name(name)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "combos").mkdir(exist_ok=True)

    arm_of = load_arm_assignment(verdict_col="alignment_verdict_gated" if gated else "alignment_verdict")
    if min_mag is not None:
        mag = {f"{gid}__{e.date}": e.magnitude for gid, g in GROUPS.items() for e in g.events}
        arm_of = {k: v for k, v in arm_of.items() if mag.get(k, 0) >= min_mag}
    event_keys = [tuple(k.split("__", 1)) for k in arm_of]
    print(f"[joint] {len(event_keys)} events with a defined arm "
          f"(noise={sum(1 for v in arm_of.values() if v=='noise')}, "
          f"signal={sum(1 for v in arm_of.values() if v=='signal')})", file=sys.stderr)

    all_series = load_event_series_for_events(event_keys)
    rng = np.random.default_rng(SEED)

    channel_labels = {"XYZ": "H", "F": "F"}
    summary_rows: list[dict] = []
    run_summary = {"seed": SEED, "n_perm": N_PERM,
                   "arm_assignment": arm_of, "combos": []}

    for channel_type, ch_label in channel_labels.items():
        for tier in STATION_TIERS:
            events = [es for (ct, t, _key), es in all_series.items() if ct == channel_type and t == tier]
            if not events:
                continue
            labels = np.array([arm_of[f"{es.group_id}__{es.event_date}"] for es in events])
            if "noise" not in labels or "signal" not in labels:
                continue

            for stat_name in STAT_NAMES:
                combo_id = f"{tier}__{ch_label}__{stat_name}"
                perm_result = permutation_test(events, labels, stat_name, rng)

                noise_events = [es for es, lab in zip(events, labels) if lab == "noise"]
                signal_events = [es for es, lab in zip(events, labels) if lab == "signal"]
                noise_stack = stack_series(noise_events, stat_name, rng)
                signal_stack = stack_series(signal_events, stat_name, rng)

                result = {
                    "combo_id": combo_id, "station_tier": tier, "channel_type_pool": channel_type,
                    "channel": ch_label, "stat": stat_name,
                    "n_noise_events": len(noise_events), "n_signal_events": len(signal_events),
                    "permutation_test": perm_result,
                    "noise_arm_stack": noise_stack, "signal_arm_stack": signal_stack,
                }
                (out_dir / "combos" / f"joint__{combo_id}.json").write_text(json.dumps(result, indent=2, default=str))

                if "error" in perm_result:
                    print(f"[{combo_id}] {perm_result['error']}", file=sys.stderr)
                    run_summary["combos"].append({"combo_id": combo_id, "status": "error",
                                                   "error": perm_result["error"]})
                    continue

                print(f"[{combo_id}] n_noise={len(noise_events)} n_signal={len(signal_events)} "
                      f"delta_peak={perm_result['delta_peak_abs_z']} p_peak={perm_result['p_value_peak']} "
                      f"delta_tail={perm_result['delta_tail_mean_abs_z']} p_tail={perm_result['p_value_tail']}",
                      file=sys.stderr)
                row = {"combo_id": combo_id, "station_tier": tier, "channel_type_pool": channel_type,
                       "channel": ch_label, "stat": stat_name,
                       "n_noise_events": len(noise_events), "n_signal_events": len(signal_events), **perm_result}
                summary_rows.append(row)
                run_summary["combos"].append({"combo_id": combo_id, "status": "ok", **perm_result})

    run_summary["min_mag"] = min_mag
    run_summary["gated"] = gated
    pd.DataFrame(summary_rows).to_csv(out_dir / "joint_summary.csv", index=False)
    (out_dir / "all_joint_run_summary.json").write_text(json.dumps(run_summary, indent=2, default=str))
    print(f"[run] {len(summary_rows)} combos -> {out_dir}", file=sys.stderr)
    return run_summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--self-test", action="store_true", help="run only the synthetic sanity check")
    ap.add_argument("--all", action="store_true", help="run the real-data joint analysis (default action)")
    ap.add_argument("--gated", action="store_true",
                    help="arms from alignment_verdict_gated (significant anomalies only), separate _gated directory")
    ap.add_argument("--min-mag", type=float, default=None,
                    help="only events of at least this magnitude, output to a separate _m<mag> directory")
    args = ap.parse_args()

    if args.self_test:
        sys.exit(0 if self_test() else 1)

    if not self_test():
        print("[main] synthetic self-test FAILED -- aborting before touching real data", file=sys.stderr)
        sys.exit(1)

    run_all(min_mag=args.min_mag, gated=args.gated)
