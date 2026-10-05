"""verify_pipeline 的基準窗口檢查不能把資料範圍以外的日子
算成乾淨的基準日。G12 的資料從 2025-05-26 開始，只比它的
2025-06-11 錨點早 16 天，所以滑動窗口會延伸到資料開始之前——那些
空白日子以前被算成乾淨日，檢查在只有 1 個真正乾淨日的情況下仍會通過。"""
import types

import numpy as np
import pandas as pd

import compute_indices as ci
import verify_pipeline as vp


def _cfg(tmp_path, start, end, anchor, storm_dates):
    days = pd.date_range(start, end).strftime("%Y%m%d")
    rows = [{"station": s, "date": d, "pct_missing": 0.0} for s in ("n1", "f1") for d in days]
    pd.DataFrame(rows).to_csv(tmp_path / "daily_features.csv", index=False)
    pd.DataFrame({"date": storm_dates, "is_storm_or_recovery": True, "is_storm_onset": True}) \
        .to_csv(tmp_path / "storm_days.csv", index=False)
    pool = types.SimpleNamespace(sufficient=True, near=("n1",), far=("f1",))
    return types.SimpleNamespace(interim_dir=tmp_path, xyz_pool=pool, f_pool=pool,
                                 anchor_event=types.SimpleNamespace(date=anchor))


def test_days_before_the_data_are_not_clean(tmp_path):
    cfg = _cfg(tmp_path, "2025-05-26", "2025-07-10", "2025-06-02", ["20250527"])
    vp.checks.clear()
    vp.check_baseline_window_excludes_storms(cfg)
    assert vp.checks["baseline_window_excludes_storms"]["pass"] is False


def test_long_clean_history_still_passes(tmp_path):
    cfg = _cfg(tmp_path, "2025-01-01", "2025-07-10", "2025-06-02", ["20250520"])
    vp.checks.clear()
    vp.check_baseline_window_excludes_storms(cfg)
    assert vp.checks["baseline_window_excludes_storms"]["pass"] is True


def test_compute_indices_has_no_baseline_right_after_data_start():
    idx = pd.date_range("2025-05-26", periods=40).strftime("%Y%m%d")
    s = pd.Series(np.random.default_rng(0).normal(size=40), index=idx)
    z = ci.mad_zscore(s, pd.Series(True, index=idx))
    assert z.iloc[:ci.MIN_CLEAN_POINTS].isna().all()
    assert z.iloc[ci.MIN_CLEAN_POINTS:].notna().all()
