"""站間轉換函數差分（ulf_interstation_tf.py，2026-11-13）的測試：
合成資料自我測試（頻帶、磁暴源場被遠站預測掉、局部訊號保留、測站離開不跳動），
以及「轉換函數的訓練日看不到被檢定的窗口」性質測試。"""
import numpy as np
import pandas as pd

import ulf_interstation_tf as ui


def test_self_test():
    assert ui.self_test()


def test_band_bins_partition_kept_frequencies():
    pc3, pc4 = ui.BAND_BINS["pc3"], ui.BAND_BINS["pc4"]
    assert not (pc3 & pc4).any()
    assert (pc3 | pc4).all()
    assert ui.KEPT_FREQS[pc4].min() >= 1 / 150 and ui.KEPT_FREQS[pc3].max() <= 1 / 10


def test_fit_transfer_recovers_known_complex_gain():
    rng = np.random.default_rng(0)
    d, s, k = 20, ui.N_SEG, int(ui.KEEP.sum())
    bf = (rng.normal(size=(d, s, k, 3)) + 1j * rng.normal(size=(d, s, k, 3))).astype(np.complex64)
    true = np.array([0.5 - 0.2j, -0.3 + 0.1j, 1.2 + 0.0j])
    zn = np.einsum("dskc,c->dsk", bf.astype(np.complex128), true)
    t = ui.fit_transfer(zn, bf)
    assert np.allclose(t, true[None, :], atol=1e-4)


def test_training_days_exclude_tested_blocks(monkeypatch):
    """被檢定的窗口（錨點震前窗口、安慰劑區塊）的資料若被改掉，轉換函數必須不變——
    也就是訓練日完全看不到它們。"""
    specs, dates, storm = ui.synthetic_specs(n_days=100, storm_days=range(10, 20), local_days=range(0))
    exclude = set(dates[60:80])
    seen = []
    orig = ui.fit_transfer

    def spy(zn, bf):
        t = orig(zn, bf)
        seen.append(t)
        return t

    monkeypatch.setattr(ui, "fit_transfer", spy)
    ui.index_table(specs, ["n1"], ["f1"], dates, storm, exclude)
    t_before = seen[-1]
    rng = np.random.default_rng(1)
    for d in exclude:  # 把被檢定的日子換成完全不同的資料
        specs["n1"][d] = specs["n1"][d] * 0 + (rng.normal(size=specs["n1"][d].shape) * 50).astype(np.complex64)
    ui.index_table(specs, ["n1"], ["f1"], dates, storm, exclude)
    assert np.allclose(seen[-1], t_before)


def test_storm_days_never_used_for_training(monkeypatch):
    specs, dates, storm = ui.synthetic_specs(n_days=60, storm_days=range(10, 20), local_days=range(0))
    captured = {}
    orig = ui.pair_energies

    def spy(sn, sf, ds, train):
        captured["train"] = set(train)
        return orig(sn, sf, ds, train)

    monkeypatch.setattr(ui, "pair_energies", spy)
    ui.index_table(specs, ["n1"], ["f1"], dates, storm, set(dates[40:50]))
    assert not (captured["train"] & storm)
    assert not (captured["train"] & set(dates[40:50]))


def test_series_diagnostics_window_counts():
    dates = pd.date_range("2024-01-01", periods=60).strftime("%Y%m%d").tolist()
    s = pd.Series(np.exp(np.random.default_rng(0).normal(size=60)), index=dates)
    storm = set(dates[0:10]) | set(dates[50:55])
    window = set(dates[45:60])
    d = ui.series_diagnostics(s, storm, window, set(dates[30:45]), pd.Series(dtype=float), log=True)
    assert d["window_days_storm_kept"] == 15
    assert d["window_days_storm_excluded"] == 10
    assert d["prev_block_days_storm_excluded"] == d["prev_block_days_storm_kept"] == 15
