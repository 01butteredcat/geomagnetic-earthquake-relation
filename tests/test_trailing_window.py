"""compute_indices.mad_zscore 的基準期必須是每一天之前的 TRAILING_WINDOW_DAYS 個日曆天
（也就是 verify_pipeline 的基準期檢查所重建的範圍），而不是前 TRAILING_WINDOW_DAYS
列。G6_G7_G8 缺 2021-12-30 和 2022-01-01 的檔案，所以以列數計算的窗口
在 G7（2022-01-03）的震前窗口內會往回延伸到 22–23 個日曆天之前。
"""
import numpy as np
import pandas as pd

import compute_indices as ci


def _series_with_gap(missing):
    days = pd.date_range("2021-12-01", "2022-01-10")
    days = days[~days.isin(pd.to_datetime(missing))]
    idx = days.strftime("%Y%m%d")
    return pd.Series(np.arange(len(idx), dtype=float), index=idx), pd.Series(True, index=idx)


def _expected(series, day):
    d = pd.to_datetime(day)
    lo = (d - pd.Timedelta(days=ci.TRAILING_WINDOW_DAYS)).strftime("%Y%m%d")
    w = series[(series.index >= lo) & (series.index < day)].to_numpy()
    med = np.median(w)
    return (series[day] - med) / (np.median(np.abs(w - med)) * 1.4826)


def test_baseline_is_calendar_window_across_missing_days():
    s, clean = _series_with_gap(["2021-12-30", "2022-01-01"])
    z = ci.mad_zscore(s, clean)
    for day in ("20211231", "20220102", "20220103", "20220110"):
        assert np.isclose(z[day], _expected(s, day)), day


def test_contiguous_series_unchanged():
    s, clean = _series_with_gap([])
    z = ci.mad_zscore(s, clean)
    for day in s.index[ci.MIN_CLEAN_POINTS + 1:]:
        assert np.isclose(z[day], _expected(s, day)), day
