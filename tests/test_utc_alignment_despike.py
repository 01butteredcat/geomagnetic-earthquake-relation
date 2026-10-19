"""震前窗口用主震 UTC 時刻對齊，以及 ULF 的去尖峰（2026-10-19）。"""
import numpy as np
import pandas as pd
import pytest

from cross_group_analysis import pre_event_anchor, pre_event_window
from events import GROUPS
from ulf_analysis import _has_spike


class _Ev:
    def __init__(self, time_utc):
        self.time_utc = time_utc


@pytest.mark.parametrize("gid, last_night", [
    ("G8", "20220321"),   # UTC 03-22 17:41：當晚夜間時段還沒結束
    ("G11", "20250119"),  # UTC 01-20 16:17：當晚夜間時段還沒開始
    ("G10", "20240402"),  # UTC 04-02 23:58：當晚夜間時段已結束，和台灣日期的舊做法相同
])
def test_last_night_is_the_night_before_the_mainshock(gid, last_night):
    win = pre_event_window(pre_event_anchor(GROUPS[gid].anchor_event), 30)
    assert win[-1] == last_night and len(win) == 30


def test_night_end_boundary():
    assert pre_event_anchor(_Ev("2024-01-10 20:00:00")) == pd.Timestamp("2024-01-11")
    assert pre_event_anchor(_Ev("2024-01-10 19:59:59")) == pd.Timestamp("2024-01-10")


def test_window_ends_before_mainshock_for_every_group():
    for g in GROUPS.values():
        last = pd.Timestamp(pre_event_window(pre_event_anchor(g.anchor_event), 30)[-1])
        assert last + pd.Timedelta(hours=20) <= pd.Timestamp(g.anchor_event.time_utc)


def test_has_spike():
    idx = pd.date_range("2024-01-01", periods=600, freq="s")
    df = pd.DataFrame({c: 30000 + np.sin(np.arange(600) / 50) for c in "XYZ"}, index=idx)
    assert not _has_spike(df)
    df.iloc[300, 2] += 500.0
    assert _has_spike(df)
