"""cross_group_analysis.sliding_baseline_rate: the tested window's own candidates
must not raise its baseline rate. Until 2026-09-28 only the identical window was
skipped, so every partially overlapping sliding window still counted them."""
import numpy as np
import pandas as pd
from hypothesis import given, settings, strategies as st

from cross_group_analysis import pre_event_window, sliding_baseline_rate


@settings(max_examples=60, deadline=None)
@given(n_days=st.integers(40, 300), w=st.sampled_from([7, 14, 30]), anchor_offset=st.integers(0, 300),
       k=st.integers(1, 5), seed=st.integers(0, 10_000))
def test_baseline_rate_ignores_candidates_inside_tested_window(n_days, w, anchor_offset, k, seed):
    """Candidates that exist only inside the tested window must not raise its own
    baseline rate: no baseline window may share a day with the tested one."""
    days = pd.date_range("2024-01-01", periods=n_days)
    all_dates = list(days.strftime("%Y%m%d"))
    anchor = days[0] + pd.Timedelta(days=w + anchor_offset % max(1, n_days - w))
    test_window = pre_event_window(anchor, w)
    rng = np.random.default_rng(seed)
    cands = set(rng.choice(test_window, size=min(k, len(test_window)), replace=False))
    rate, n = sliding_baseline_rate(all_dates, cands, w, set(test_window))
    assert n == 0 or rate == 0.0
