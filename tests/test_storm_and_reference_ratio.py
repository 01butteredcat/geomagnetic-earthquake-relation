"""2026-11-02 兩個敏感度分析的測試：
1. 遠站參考比值（ulf_reference_ratio.py，Han 2014／Warden 2020）：DWT 頻帶、全球性變化被遠站消掉、
   局部增強會讓 P 上升（合成資料的自我測試）。
2. 磁暴日放回去（surrogate_test.storm_kept_summary）：有真實輸出時，「排除磁暴日」版本必須逐位重現
   主要檢定 window_summary.json 的 p 值，而且「同組」比較不能混進排除版本沒有 p 值的組。"""
import json

import numpy as np
import pandas as pd
import pytest

import ulf_reference_ratio as urr
from common import PROJECT_DIR

SURR = PROJECT_DIR / "data" / "interim" / "surrogate_test"


def test_reference_ratio_self_test():
    assert urr.self_test()


def test_dwt_level6_keeps_100s_rejects_10s():
    t = np.arange(urr.N_DWT, dtype=float)
    slow = urr.band_energy(urr.dwt_level_detail(np.sin(2 * np.pi * t / 100)))
    fast = urr.band_energy(urr.dwt_level_detail(np.sin(2 * np.pi * t / 10)))
    assert slow > 0.1 and fast < 1e-4


def test_pair_ratio_cancels_common_scaling():
    rng = np.random.default_rng(0)
    idx = pd.Index([f"2024{m:02d}{d:02d}" for m in (1, 2, 3) for d in range(1, 29)])
    src = pd.Series(rng.lognormal(0, 1, len(idx)), index=idx)  # 共同源場（含磁暴的大起伏）
    near = 2.0 * src * rng.lognormal(0, 0.05, len(idx))
    far = src * rng.lognormal(0, 0.05, len(idx))
    p, _ = urr.pair_ratio(near, far, pd.Series(True, index=idx))
    assert np.nanmedian(np.abs(p - 1)) < 0.1


@pytest.mark.skipif(not (SURR / "storm_kept_summary.json").exists(), reason="沒有 surrogate_test 的輸出")
def test_storm_excluded_reproduces_main_and_same_groups_is_matched():
    main = json.loads((SURR / "window_summary.json").read_text())["rank_test"]
    kept = json.loads((SURR / "storm_kept_summary.json").read_text())["rank_test"]
    for band in ("pc3", "pc4"):
        for name in ("all_groups", "one_per_family"):
            assert kept["storm_excluded"][band][name]["p"] == main[band][name]["p"]
            assert kept["storm_kept_same_groups"][band][name]["groups"] == main[band][name]["groups"]
