"""去尖峰改用 INTERMAGNET 技術手冊 5.2.5 的雙邊 100 nT 規則（2026-10-20）。"""
import numpy as np
import pandas as pd

from build_daily_features import SPIKE_THRESHOLD_NT, _despike, _one_sided_jumps


def _flat(n=20):
    return pd.Series(np.full(n, 30000.0))


def test_threshold_is_manual_value():
    assert SPIKE_THRESHOLD_NT == 100.0


def test_single_spike_flagged_but_not_its_neighbours():
    s = _flat(); s[10] += 150.0
    cleaned, mask, _ = _despike(s)
    assert list(np.flatnonzero(mask)) == [10]
    assert np.isnan(cleaned[10]) and cleaned[11] == 30000.0
    assert _one_sided_jumps(s) == 0


def test_threshold_boundary():
    s = _flat(); s[5] += 99.0; s[15] += 101.0
    assert list(np.flatnonzero(_despike(s)[1])) == [15]


def test_step_is_not_a_spike_but_is_counted():
    s = _flat(); s[10:] += 500.0
    assert not _despike(s)[1].any()
    assert _one_sided_jumps(s) == 1


def test_point_next_to_gap_not_flagged():
    s = _flat(); s[9] = np.nan; s[10] += 500.0
    assert not _despike(s)[1].any()


def test_edges_not_flagged():
    s = _flat(); s[0] += 500.0; s[19] += 500.0
    assert not _despike(s)[1].any()
