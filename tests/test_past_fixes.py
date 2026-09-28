"""Regression tests for statistical bugs fixed before 2026-09-28."""
import csv
import json
import re
import sys
import types
from collections import Counter
from datetime import datetime, timedelta, timezone

import numpy as np
import pandas as pd
import pytest

import catalog_utils
import fetch_earthquake_catalog as fec
import fetch_space_weather as fsw
import method_comparison as mc
from coseismic_step_analysis import keyed_rng
from cross_group_analysis import pre_event_window
from events import GROUPS, get_group


# --- precursor p-value: rank test instead of a floored resampling p (a5108fe) -----

def _blocks_series(block_minima, anchor="2024-06-30"):
    """One value per day; each 30-day block tiled back from `anchor` has the given
    minimum (last entry = the real pre-event window), every other day is 0."""
    w = mc.PRE_WINDOW_DAYS
    end = pd.Timestamp(anchor) - pd.Timedelta(days=1)
    idx = pd.date_range(end - pd.Timedelta(days=w * len(block_minima) - 1), end)
    vals = np.zeros(len(idx))
    for b, m in enumerate(block_minima):
        vals[b * w + 5] = m
    return pd.Series(vals, index=idx.strftime("%Y%m%d")), set(pre_event_window(pd.Timestamp(anchor), w))


def test_rank_p_is_rank_not_resampling_floor():
    """G11 pc3 case: window min just below the background min. The old leave-window
    bootstrap gave the floor 1/(N+1) = 0.0005; the rank among 4 fake blocks is 1/5."""
    s, real = _blocks_series([-0.5, -0.6, -0.7, -0.73, -1.2])
    out = mc.rank_window_test(s, real, tail="lower")
    assert out["rank_n_fake_windows"] == 4
    assert out["rank_p"] == pytest.approx(0.2)
    assert out["rank_min_attainable_p"] == pytest.approx(0.2)


def test_rank_p_counts_ties_and_worse_blocks():
    s, real = _blocks_series([-2.0, -0.5, -1.0, -1.0])  # real = -1.0; one worse, one tie
    out = mc.rank_window_test(s, real, tail="lower")
    assert out["rank_p"] == pytest.approx((1 + 2) / 4)
    assert out["null_ps"][-1] == out["rank_p"]


def test_fisher_null_is_calibrated_under_h0():
    rng = np.random.default_rng(0)
    null = [[(k + 1) / 5 for k in range(5)] for _ in range(12)]
    rejections = 0
    for _ in range(200):
        ps = [float(rng.choice(n)) for n in null]
        rejections += mc.fisher_across_groups(ps, null, rng, n_sim=2000)["p"] < 0.05
    assert rejections / 200 < 0.09


# --- storm-cache check (a5108fe) ---------------------------------------------------

def _run_check_cache(monkeypatch, tmp_path, summary):
    (tmp_path / "storm_days.csv").write_text("date,is_storm_or_recovery,is_storm_onset\n")
    (tmp_path / "storm_days_summary.json").write_text(json.dumps(summary))
    cfg = types.SimpleNamespace(interim_dir=tmp_path)
    monkeypatch.setattr(fsw, "load_group_config", lambda g: cfg)
    monkeypatch.setattr(fsw, "_group_date_range", lambda c: ("2021-07-23", "2022-05-31"))
    monkeypatch.setattr(sys, "argv", ["fetch_space_weather.py", "--group", "G6", "--check-cache"])
    with pytest.raises(SystemExit) as e:
        fsw.main()
    return e.value.code


def test_cache_built_before_folder_was_extended_is_refetched(monkeypatch, tmp_path):
    stale = {"confidence": "high (official Kp/Dst indices)", "date_range": ["2021-07-23", "2022-04-30"]}
    assert _run_check_cache(monkeypatch, tmp_path, stale) == 1


def test_partial_cache_is_refetched(monkeypatch, tmp_path):
    partial = {"confidence": "medium (Dst unavailable)", "date_range": ["2021-07-23", "2022-05-31"]}
    assert _run_check_cache(monkeypatch, tmp_path, partial) == 1


def test_complete_cache_is_reused(monkeypatch, tmp_path):
    ok = {"confidence": "high (official Kp/Dst indices)", "date_range": ["2021-07-23", "2022-05-31"]}
    assert _run_check_cache(monkeypatch, tmp_path, ok) == 0


# --- same-day event keys (suffixed dates, keyed null streams) -----------------------

def test_event_keys_unique_within_each_group():
    for gid, g in GROUPS.items():
        dup = [d for d, n in Counter(e.date for e in g.events).items() if n > 1]
        assert not dup, f"{gid}: duplicate event.date keys {dup}"


def test_same_local_day_events_are_suffixed_and_consistent():
    for gid, g in GROUPS.items():
        # convention: one event per local day may keep the bare date, the others get a/b/c suffixes
        unsuffixed = Counter(e.date for e in g.events if len(e.date) == 10)
        for e in g.events:
            assert re.fullmatch(r"\d{4}-\d{2}-\d{2}[a-z]?", e.date), (gid, e.date)
            assert e.date[:10] == e.time_local[:10], (gid, e.date, e.time_local)
        assert all(n == 1 for n in unsuffixed.values()), gid


def test_keyed_null_streams_differ_for_same_day_events():
    g10 = get_group("G10")
    same_day = [e for e in g10.events if e.date.startswith("2024-04-23")]
    assert len(same_day) >= 2
    draws = {e.date: keyed_rng("G10", e.date, "XYZ", "xcg", "H").random(5).tobytes() for e in same_day}
    assert len(set(draws.values())) == len(same_day)


def test_keyed_rng_is_order_independent():
    a = keyed_rng("G10", "2024-04-03", "XYZ", "xcg", "H").random(3)
    keyed_rng("G9", "2022-09-18", "F", "csg", "F").random(1000)
    b = keyed_rng("G10", "2024-04-03", "XYZ", "xcg", "H").random(3)
    assert np.array_equal(a, b)


# --- magnitude-based matching and tier filtering (f0a3375, 15b2842) ----------------

def _utc(ev, dt_sec=0):
    return datetime.strptime(ev.time_utc, "%Y-%m-%d %H:%M:%S").replace(tzinfo=timezone.utc) + timedelta(seconds=dt_sec)


def test_catalog_row_with_different_magnitude_scale_is_still_known():
    """G4 2020-12-10: USGS 6.1 vs CWA ML6.64, 0 s apart -- was counted twice."""
    ev = get_group("G4").anchor_event
    rows = [{"_dt": _utc(ev), "lat": ev.lat, "lon": ev.lon, "mag": 6.1}]
    fec.flag_known_events(rows, "G4")
    assert rows[0]["is_known_event"] is True


def test_distinct_aftershock_hours_later_is_not_known():
    """G9 2022-09-18 09:39 M5.9: inside the old +-6 h window but a different event."""
    ev = get_group("G9").anchor_event
    rows = [{"_dt": _utc(ev, 3 * 3600), "lat": ev.lat, "lon": ev.lon, "mag": ev.magnitude}]
    fec.flag_known_events(rows, "G9")
    assert rows[0]["is_known_event"] is False


def test_same_time_far_away_is_not_known():
    ev = get_group("G4").anchor_event
    rows = [{"_dt": _utc(ev), "lat": ev.lat + 1.5, "lon": ev.lon, "mag": ev.magnitude}]
    fec.flag_known_events(rows, "G4")
    assert rows[0]["is_known_event"] is False


def test_tier_filter_drops_lower_magnitudes(tmp_path):
    path = tmp_path / "catalog.csv"
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["group", "time_utc", "mag", "declustered", "is_known_event", "source"])
        w.writeheader()
        w.writerow({"group": "G11", "time_utc": "2024-12-01 00:00:00", "mag": "5.4", "declustered": "True",
                    "is_known_event": "False", "source": "CWA_GDMS"})
        w.writerow({"group": "G11", "time_utc": "2024-12-02 00:00:00", "mag": "6.1", "declustered": "True",
                    "is_known_event": "False", "source": "CWA_GDMS"})
    for tier in (6.0, 5.5, 5.0):
        evs = catalog_utils.load_extended_events(path, ("G11", "G12"), tier)
        assert evs and all(e["mag"] >= tier for e in evs)
    m6 = catalog_utils.load_extended_events(path, ("G11", "G12"), 6.0)
    registered_m6 = sum(e.magnitude >= 6.0 for g in ("G11", "G12") for e in get_group(g).events)
    assert sum(e["source"] == "events.py" for e in m6) == registered_m6
