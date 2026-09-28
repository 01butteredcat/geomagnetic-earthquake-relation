"""compute_indices.mad_zscore's baseline must be the TRAILING_WINDOW_DAYS calendar
days before each day (what verify_pipeline's baseline check reconstructs), not
the previous TRAILING_WINDOW_DAYS rows. G6_G7_G8 has no files for 2021-12-30 and
2022-01-01, so a row window reached 22-23 calendar days back inside G7's
(2022-01-03) pre-event window."""
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
