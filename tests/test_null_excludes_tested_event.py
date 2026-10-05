"""被檢定的窗口／事件絕不能進入它自己的虛無分布。

整條序列區塊 bootstrap 的回歸測試（2026-09-24 修正，32cc9cb）：被檢定的
極值會被重抽進它自己的替代序列，讓 p 有約 0.2–0.4 的下限。
下面的性質測試會檢查每一個應該是
留一窗／留一事件的虛無分布產生器。"""
import numpy as np
import pytest
import pandas as pd
from hypothesis import given, settings, strategies as st

import coseismic_step_analysis as csa
import method_comparison as mc
from cross_group_analysis import pre_event_window

SENTINEL = 1e6  # 任何背景日都不可能出現的值


def _series(n_days, seed, start="2024-01-01"):
    rng = np.random.default_rng(seed)
    idx = pd.date_range(start, periods=n_days).strftime("%Y%m%d")
    return pd.Series(rng.normal(size=n_days), index=idx)


# --- surrogate_test (B)：留一窗 ---------------------------------------

def test_huge_window_anomaly_reaches_min_attainable_p():
    """窗口不放進虛無分布時，50 個 sigma 的下凹會得到 p = 1/(N+1)，
    而不是整條序列 bootstrap 的約 0.2–0.4 下限。"""
    s = _series(150, 0)
    pos = np.arange(100, 130)
    values = s.to_numpy().copy()
    values[115] = -50.0
    n = 200
    out = mc.surrogate_test(values, pos, np.random.default_rng(1), n_surr=n, tail="lower")
    assert out["p_window_block_bootstrap"] == round(1 / (n + 1), 5)
    assert out["p_window_phase_randomization"] == round(1 / (n + 1), 5)


@settings(max_examples=30, deadline=None)
@given(n_days=st.integers(60, 200), win_start=st.integers(0, 200), seed=st.integers(0, 10_000),
       tail=st.sampled_from(["two", "lower"]))
def test_window_values_never_reach_window_surrogates(n_days, win_start, seed, tail):
    w = mc.PRE_WINDOW_DAYS
    win_start = win_start % max(1, n_days - w)
    values = _series(n_days, seed).to_numpy().copy()
    pos = np.arange(win_start, min(win_start + w, n_days))
    values[pos] = SENTINEL + np.arange(len(pos))

    seen = []
    real_bb, real_pr = mc.block_bootstrap_surrogate, mc.phase_randomize_surrogate

    def spy(fn):
        def wrapped(x, *a, **k):
            seen.append(np.asarray(x).copy())
            return fn(x, *a, **k)
        return wrapped

    # (A) 整條序列的相位隨機化本來就可以用全部資料；只有 (B) 是留一窗，
    # 而且 (B) 是唯一會傳入比 `values` 短的序列的呼叫者。
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(mc, "block_bootstrap_surrogate", spy(real_bb))
        mp.setattr(mc, "phase_randomize_surrogate", spy(real_pr))
        mc.surrogate_test(values, pos, np.random.default_rng(seed), n_surr=5, tail=tail)
    window_inputs = [x for x in seen if len(x) == n_days - len(pos)]
    for x in window_inputs:
        assert not np.any(x >= SENTINEL)


# --- rank_window_test：假區段絕不和真實窗口重疊 ----------------

@settings(max_examples=40, deadline=None)
@given(n_days=st.integers(70, 400), anchor_offset=st.integers(31, 400), seed=st.integers(0, 10_000))
def test_rank_fake_blocks_disjoint_from_real_window(n_days, anchor_offset, seed):
    s = _series(n_days, seed)
    anchor = pd.to_datetime(s.index[0]) + pd.Timedelta(days=anchor_offset % n_days)
    real = set(pre_event_window(anchor, mc.PRE_WINDOW_DAYS))

    windows = []

    def spy(a, d):
        w = pre_event_window(a, d)
        windows.append(set(w))
        return w

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(mc, "pre_event_window", spy)
        out = mc.rank_window_test(s, real, tail="lower")
    for w in windows:
        assert not (w & real)
    # 假區段是鋪排的，所以彼此之間也不重疊
    for i in range(len(windows)):
        for j in range(i + 1, len(windows)):
            assert not (windows[i] & windows[j])
    if out["rank_p"] is not None:
        assert out["rank_p"] >= out["rank_min_attainable_p"]


# --- 同震虛無參考時間 -----------------------------------------------

@settings(max_examples=50, deadline=None)
@given(seed=st.integers(0, 2**32 - 1), n_events=st.integers(1, 6), span_h=st.integers(2, 48))
def test_coseismic_null_centers_stay_off_every_real_event(seed, n_events, span_h):
    rng = np.random.default_rng(seed)
    lo = pd.Timestamp("2024-04-01")
    hi = lo + pd.Timedelta(hours=span_h)
    events = [lo + pd.Timedelta(seconds=float(s)) for s in rng.uniform(0, span_h * 3600, n_events)]
    centers = csa._draw_null_centers(rng, lo, hi, events, csa.EXCLUSION_BUFFER_SEC, 200)
    for c in centers:
        assert all(abs((c - e).total_seconds()) >= csa.EXCLUSION_BUFFER_SEC for e in events)


def test_coseismic_null_search_window_cannot_reach_the_origin():
    """距離 ≥ EXCLUSION_BUFFER_SEC 的虛無中心會掃描 ±SCAN_HALF_SEC，而它的
    階躍統計量會再往外讀 max(STEP_WINDOWS_SEC)。"""
    assert csa.EXCLUSION_BUFFER_SEC > csa.SCAN_HALF_SEC + max(csa.STEP_WINDOWS_SEC)


# --- SEA 假時期（15b2842） -------------------------------------------------------

@settings(max_examples=30, deadline=None)
@given(group=st.sampled_from(["G10", "G11", "G13", "G6"]), offsets=st.lists(st.integers(0, 150), max_size=4))
def test_sea_fake_epochs_stay_off_every_real_event(group, offsets):
    import superposed_epoch_analysis as sea
    from events import folder_events, get_group

    anchor = pd.Timestamp(get_group(group).anchor_event.time_utc.split(" ")[0])
    days = pd.date_range(anchor - pd.Timedelta(days=100), anchor + pd.Timedelta(days=60))
    series = {"dates": list(days.strftime("%Y%m%d")), "pc3": {}}
    events = [{"group": group, "date": days[0] + pd.Timedelta(days=o)} for o in offsets]
    events.append({"group": group, "date": anchor})
    eligible = sea._eligible_null_days("pc3", events, {group: series})[group]
    real = [e["date"] for e in events] + [pd.Timestamp(e.time_utc.split(" ")[0]) for e in folder_events(group)]
    for d in eligible:
        assert all(abs((d - r).days) >= sea.NULL_EXCLUSION_BUFFER_DAYS for r in real)
