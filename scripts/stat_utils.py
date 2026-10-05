"""教授建議的驗證方法所共用的統計輔助函式：
穩健的中位數／MAD
z-score 公式（抽出來，讓 `surrogate_test.py`、`superposed_epoch_
analysis.py` 和 `backtest_rule.py` 都使用
`cross_group_analysis.py::ulf_candidate_dates` 已建立的完全相同定義），加上
兩個替代序列產生器，用來建立保留序列
「統計外觀」（自相關／頻譜形狀）、但
和地震時間沒有任何真實關係的虛無分布。

FIXED_RULE_THRESHOLD = -4.1 是 `report_template.html` 對
G10 2024-03-30 Pc3 近站／遠站極化低谷引用的數字（「z ~= -4.1 MAD 單位，
在 113 天可計算窗口中最極端的一天」）。在這個模組之前，它從來不是
程式碼中任何地方的程式化門檻——完整的來源脈絡見
計畫檔。這裡把它升格成具名常數，
就是為了能把它當成真正的規則來回測，而不是只
出現在一次性的敘述文字裡。
"""
from __future__ import annotations

from typing import Callable

import numpy as np

MAD_SCALE = 1.4826  # -> 在真正常態下等同常態的標準差（跨組
                     # 分析自己的選擇；不是宣稱資料**是**常態）
FIXED_RULE_THRESHOLD = -4.1


def mad_zscore(values: np.ndarray) -> tuple[np.ndarray, float, float]:
    """整條序列中位數／MAD z-score，公式和
    cross_group_analysis.py::ulf_candidate_dates 相同。回傳 (z, median, mad)。
    mad 已乘上 MAD_SCALE。由於用中位數／MAD 而非平均／標準差，
    本質上對非常態／厚尾資料穩健。"""
    med = float(np.median(values))
    mad = float(np.median(np.abs(values - med))) * MAD_SCALE
    if mad <= 1e-9:
        return np.full_like(values, np.nan, dtype=float), med, mad
    z = (values - med) / mad
    return z, med, mad


def estimate_block_length(x: np.ndarray, max_lag: int = 30) -> int:
    """區塊 bootstrap 區塊大小的去相關長度經驗法則：
    自相關函數第一次掉到 1/e 以下的延遲，下限
    3 天、上限 n//5，避免單一區塊主導很短的
    序列。這是標準、站得住腳但簡單的選擇（不是像 politis-white
    那樣擬合出的最佳區塊長度）——足以達到「相同的短
    程自相關外觀」，而這正是教授的建議
    所要求的。"""
    n = len(x)
    xc = x - x.mean()
    var = np.dot(xc, xc) / n
    if var <= 1e-12:
        return 3
    acf = []
    for lag in range(1, min(max_lag, n - 1) + 1):
        c = np.dot(xc[:-lag], xc[lag:]) / n / var
        acf.append(c)
        if abs(c) < 1 / np.e:
            return max(3, min(lag, n // 5))
    return max(3, min(max_lag, n // 5))


def block_bootstrap_surrogate(x: np.ndarray, block_len: int, rng: np.random.Generator) -> np.ndarray:
    """循環區塊 bootstrap：用隨機挑選的連續區塊（在尾端
    繞回開頭）重組出等長的替代序列，保留
    區域（區塊內）自相關，同時破壞任何低谷／突波和日曆日期之間
    特定的對應關係——正是教授描述的
    「打亂差值序列，但保留它的統計外觀」。
    """
    n = len(x)
    out = np.empty(n, dtype=float)
    filled = 0
    while filled < n:
        start = rng.integers(0, n)
        take = min(block_len, n - filled)
        idx = (start + np.arange(take)) % n
        out[filled:filled + take] = x[idx]
        filled += take
    return out


def phase_randomize_surrogate(x: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """FFT 相位隨機化替代序列：保留完全相同的振幅頻譜
    （所以替代序列和原序列有相同的自相關／「頻域
    外觀」），但把每個相位換成獨立的
    均勻隨機抽樣，這會破壞任何綁在特定日期上的局部特徵（例如
    持續數天的低谷）。這是標準的非線性時間序列
    替代資料方法（Fourier 轉換替代序列）；強制 Hermitian 對稱，
    讓逆轉換是實數，DC/Nyquist 分量
    保持實數（零相位），因為它們沒有有意義的相位。"""
    n = len(x)
    mean = x.mean()
    xc = x - mean
    X = np.fft.rfft(xc)
    amp = np.abs(X)
    phases = rng.uniform(0, 2 * np.pi, size=len(X))
    phases[0] = 0.0
    if n % 2 == 0:
        phases[-1] = 0.0
    X_new = amp * np.exp(1j * phases)
    x_new = np.fft.irfft(X_new, n=n)
    return x_new + mean


def bootstrap_ci(
    arrays: tuple[np.ndarray, ...],
    statistic: Callable[..., float | np.ndarray],
    n_bootstrap: int,
    rng: np.random.Generator,
    ci: float = 0.90,
) -> dict:
    """statistic(*arrays) 的百分位數 bootstrap 信賴區間，逐列**配對**
    重抽，`arrays` 中每個陣列共用（同一組隨機列
    索引套用到全部陣列），讓配對結構——例如一個磁暴
    日的 (|near_index|, |local_anomaly_index|) 配對——在每次
    重抽中保留，而不是把一條序列和另一條獨立
    重抽的序列錯開。只傳一個陣列時就退化成普通的
    bootstrap，也就是 superposed_epoch_analysis.py 手寫的
    逐事件疊加 bootstrap 已經在做的（n_events x n_lags 矩陣，
    statistic=lambda m: np.nanmean(m, axis=0)）——這個基本函式夠通用，
    未來重構那支腳本（以及同樣模式的 coseismic_stacking_
    analysis.py）時可以直接用 arrays=(M,) 呼叫它，不需要
    更改介面。

    `statistic` 可以回傳純量或固定形狀的 ndarray；百分位數
    會用 nanpercentile 沿著 bootstrap 軸逐元素取，
    所以 `statistic` 回傳 NaN 的重抽（例如分母為零的
    退化重抽）會被容忍，而不會汙染
    整個信賴區間——呼叫者在信任 ci_lo/ci_hi 之前仍應檢查它們是否為 NaN，
    因為全部是 NaN 的 bootstrap 代表這個統計量
    對這份資料沒有定義。

    回傳 {"point_estimate": 在**真實**、未重抽資料上的 statistic(*arrays)
    （不是 bootstrap 值），"ci_lo"、"ci_hi"（在 `ci` 水準，例如
    0.90 -> 第 5／95 百分位數），"ci_level": ci，"n_bootstrap": n_bootstrap}，
    全部轉成普通 float/list，讓呼叫者可以直接放進
    JSON 報告。不回傳原始 bootstrap 陣列——需要
    完整分布的呼叫者應該自己跑迴圈；目前沒有呼叫者
    需要它。"""
    n = arrays[0].shape[0]
    if any(a.shape[0] != n for a in arrays):
        raise ValueError("配對重抽時，所有陣列在 axis 0 的長度必須相同")

    point = statistic(*arrays)
    boot = None
    for b in range(n_bootstrap):
        idx = rng.integers(0, n, size=n)
        val = statistic(*(a[idx] for a in arrays))
        if boot is None:
            boot = np.full((n_bootstrap,) + np.shape(val), np.nan, dtype=float)
        boot[b] = val

    lo_pct, hi_pct = (1 - ci) / 2 * 100, (1 + ci) / 2 * 100
    ci_lo = np.nanpercentile(boot, lo_pct, axis=0)
    ci_hi = np.nanpercentile(boot, hi_pct, axis=0)

    def _cast(x):
        arr = np.asarray(x)
        return float(arr) if arr.ndim == 0 else arr.tolist()

    return {
        "point_estimate": _cast(point),
        "ci_lo": _cast(ci_lo),
        "ci_hi": _cast(ci_hi),
        "ci_level": ci,
        "n_bootstrap": n_bootstrap,
    }


def histogram_summary(values: np.ndarray, bins: int = 40) -> dict:
    """替代統計量分布的精簡、適合 JSON 的摘要，
    給報告畫圖用：分組直方圖 + 重點百分位數，而不是
    把幾千個原始浮點數倒進 report_data.json。"""
    counts, edges = np.histogram(values, bins=bins)
    pct = np.percentile(values, [1, 5, 25, 50, 75, 95, 99])
    return {
        "bin_edges": [round(float(e), 4) for e in edges],
        "counts": [int(c) for c in counts],
        "mean": round(float(np.mean(values)), 4),
        "std": round(float(np.std(values)), 4),
        "p1": round(float(pct[0]), 4), "p5": round(float(pct[1]), 4),
        "p25": round(float(pct[2]), 4), "p50": round(float(pct[3]), 4),
        "p75": round(float(pct[4]), 4), "p95": round(float(pct[5]), 4),
        "p99": round(float(pct[6]), 4),
    }


SELF_TEST_SEED = 20260805


def self_test() -> bool:
    """用已知答案的合成資料健全性檢查 bootstrap_ci()，
    涵蓋這個模組的呼叫者實際需要的兩種形狀：
    兩個等長陣列上的配對比值統計量（verify_pipeline.
    py 的 check_storm_cancellation 用法），以及單一陣列的逐延遲平均
    （superposed_epoch_analysis.py 既有的手寫 bootstrap，
    這個函式可以在不改介面的情況下取代它）。"""
    rng = np.random.default_rng(SELF_TEST_SEED)
    ok = True

    # 情況 A：配對比值統計量，n 很大——信賴區間應該緊緊包住
    # 真實比值。
    def _median_ratio(a, b):
        med_a = np.median(a)
        return float(np.median(b) / med_a) if med_a > 0 else float("nan")

    # 比值本身有 20% 的相對雜訊（不是很小的加法項），讓
    # 統計量的抽樣變異——以及它隨 n 的縮小——
    # 大到能從單次固定種子的抽樣可靠地偵測出來。
    true_ratio = 0.4
    n_large = 200
    near = rng.uniform(1, 5, size=n_large)
    local = near * true_ratio * (1 + rng.normal(0, 0.2, size=n_large))
    res_large = bootstrap_ci((near, local), _median_ratio, 2000, rng, ci=0.90)
    close_to_truth = abs(res_large["point_estimate"] - true_ratio) < 0.1
    status = "PASS" if close_to_truth else "FAIL"
    if status == "FAIL":
        ok = False
    print(f"[self-test] 配對比值 n={n_large}: point={res_large['point_estimate']:.3f} "
          f"ci=[{res_large['ci_lo']:.3f}, {res_large['ci_hi']:.3f}] (true={true_ratio})  {status}")

    # 情況 A 續：同樣的真實比值／雜訊，n 很小——信賴區間應該
    # 明顯比 n 很大的情況寬（這是
    # verify_pipeline.py 中 bootstrap 檢定力不足警告所依賴的性質）。
    n_small = 5
    near_s = rng.uniform(1, 5, size=n_small)
    local_s = near_s * true_ratio * (1 + rng.normal(0, 0.2, size=n_small))
    res_small = bootstrap_ci((near_s, local_s), _median_ratio, 2000, rng, ci=0.90)
    width_large = res_large["ci_hi"] - res_large["ci_lo"]
    width_small = res_small["ci_hi"] - res_small["ci_lo"]
    status = "PASS" if width_small > width_large else "FAIL"
    if status == "FAIL":
        ok = False
    print(f"[self-test] 配對比值 n={n_small}: ci=[{res_small['ci_lo']:.3f}, {res_small['ci_hi']:.3f}] "
          f"width={width_small:.3f} （預期比 n={n_large} 的寬度 {width_large:.3f} 寬）  {status}")

    # 情況 B：單一陣列的逐延遲平均，仿照 superposed_epoch_
    # analysis.py 的 run_band——M 是 (n_events x n_lags)，統計量對
    # 事件做化約，信賴區間應該包住已知注入的逐延遲平均。這裡
    # 刻意用很高（0.999）的信賴水準，讓單次固定種子的
    # 包含檢查穩健（一個校準正確的 90% 信賴區間本來就*預期*
    # 每個延遲約有 10% 的機率沒包住目標——那不是 bug——所以
    # 斷言 90% 信賴區間永遠包住真值本身就是不穩定的
    # 測試；0.999 讓偶然沒包住的機率小到天文數字）。
    n_events, n_lags = 60, 5
    lag_means = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    M = lag_means + rng.normal(0, 0.3, size=(n_events, n_lags))
    res_m = bootstrap_ci((M,), lambda m: np.nanmean(m, axis=0), 2000, rng, ci=0.999)
    ci_lo, ci_hi = np.array(res_m["ci_lo"]), np.array(res_m["ci_hi"])
    right_shape = ci_lo.shape == (n_lags,) and ci_hi.shape == (n_lags,)
    brackets_truth = bool(np.all((ci_lo <= lag_means) & (lag_means <= ci_hi)))
    status = "PASS" if (right_shape and brackets_truth) else "FAIL"
    if status == "FAIL":
        ok = False
    print(f"[self-test] 單一陣列的逐延遲平均：ci_lo={np.round(ci_lo, 2)} "
          f"ci_hi={np.round(ci_hi, 2)} (true={lag_means})  {status}")

    return ok


if __name__ == "__main__":
    import sys

    sys.exit(0 if self_test() else 1)
