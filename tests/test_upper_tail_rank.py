"""上尾排名檢定（2026-11-01，探索性）：文獻多報告震前 Z/H「上升」，
所以另外看震前區塊的最高值排第幾。上尾必須和「把序列取負號後的下尾」完全相同，
而且不能改動現行的下尾結果。"""
import numpy as np
import pandas as pd
import pytest
from hypothesis import given, settings, strategies as st

import method_comparison as mc
from cross_group_analysis import pre_event_window


def _series(values, anchor="2024-04-03"):
    end = pd.Timestamp(anchor) - pd.Timedelta(days=1)
    idx = pd.date_range(end=end, periods=len(values)).strftime("%Y%m%d")
    s = pd.Series(np.asarray(values, dtype=float), index=idx)
    return s, set(pre_event_window(pd.Timestamp(anchor), mc.PRE_WINDOW_DAYS))


def test_upper_peak_in_window_ranks_first():
    w = mc.PRE_WINDOW_DAYS
    vals = np.zeros(5 * w)
    vals[-3] = 9.0  # 震前窗口內的高峰
    for b in range(4):
        vals[b * w + 5] = 1.0 + b  # 其他區塊的峰都比較低
    s, real = _series(vals)
    up = mc.rank_window_test(s, real, tail="upper")
    assert up["rank_n_fake_windows"] == 4
    assert up["rank_p"] == pytest.approx(0.2)
    # 同一筆資料的下尾不受影響：每個區塊的最低都是 0，全部同分
    assert mc.rank_window_test(s, real, tail="lower")["rank_p"] == pytest.approx(1.0)


@settings(max_examples=50, deadline=None)
@given(st.lists(st.floats(-5, 5, allow_nan=False), min_size=4 * mc.PRE_WINDOW_DAYS, max_size=6 * mc.PRE_WINDOW_DAYS))
def test_upper_equals_lower_of_negated(values):
    s, real = _series(values)
    up = mc.rank_window_test(s, real, tail="upper")
    lo = mc.rank_window_test(-s, real, tail="lower")
    assert up["rank_p"] == lo["rank_p"]
    assert up.get("null_ps") == lo.get("null_ps")
