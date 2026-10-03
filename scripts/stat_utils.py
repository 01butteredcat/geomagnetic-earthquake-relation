"""Shared statistics helpers for the professor-suggested validation methods:
the robust median/MAD
z-score formula (factored out so `surrogate_test.py`, `superposed_epoch_
analysis.py` and `backtest_rule.py` all use the exact same definition
`cross_group_analysis.py::ulf_candidate_dates` already established), plus
two surrogate-series generators for building a null distribution that
preserves a series' "statistical look" (autocorrelation / spectral shape)
without any real relationship to earthquake timing.

FIXED_RULE_THRESHOLD = -4.1 is the number `report_template.html` quotes for
G10's 2024-03-30 Pc3 near/far polarization dip ("z ~= -4.1 MAD units, the
most extreme day in the 113-day computable window"). It was never a
programmatic threshold anywhere in the codebase before this module -- see the
plan file for the full provenance trail. It's promoted to a named constant
here specifically so it can be backtested as an actual rule instead of only
ever appearing as one-off narrative text.
"""
from __future__ import annotations

from typing import Callable

import numpy as np

MAD_SCALE = 1.4826  # -> normal-equivalent std under a true normal (Cross-group
                     # analysis's own choice; not a claim the data IS normal)
FIXED_RULE_THRESHOLD = -4.1


def mad_zscore(values: np.ndarray) -> tuple[np.ndarray, float, float]:
    """Whole-series median/MAD z-score, same formula as
    cross_group_analysis.py::ulf_candidate_dates. Returns (z, median, mad).
    mad is already scaled by MAD_SCALE. Robust to non-normal/heavy-tailed
    data by construction (median/MAD instead of mean/std)."""
    med = float(np.median(values))
    mad = float(np.median(np.abs(values - med))) * MAD_SCALE
    if mad <= 1e-9:
        return np.full_like(values, np.nan, dtype=float), med, mad
    z = (values - med) / mad
    return z, med, mad


def estimate_block_length(x: np.ndarray, max_lag: int = 30) -> int:
    """Decorrelation-length heuristic for the block bootstrap's block size:
    first lag where the autocorrelation function drops below 1/e, floored at
    3 days and capped at n//5 so a single block can't dominate a short
    series. This is a standard, defensible-but-simple choice (not a fitted
    optimal block length a la politis-white) -- adequate for "same short-
    range autocorrelation look", which is all the professor's suggestion
    asks for."""
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
    """Circular block bootstrap: reassemble a same-length surrogate from
    randomly-chosen contiguous blocks (wrapping around the end), preserving
    local (within-block) autocorrelation while destroying the specific
    alignment between any dip/spike and calendar date -- exactly the
    "shuffle the differential series but keep its statistical look" the
    professor described."""
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
    """FFT phase-randomization surrogate: keeps the exact amplitude spectrum
    (so the surrogate has the same autocorrelation / "frequency-domain
    look" as the original) but replaces every phase with an independent
    uniform random draw, which destroys any localized feature (like a
    multi-day dip) tied to a specific date. Standard nonlinear-time-series
    surrogate-data method (Fourier-transform surrogates); Hermitian symmetry
    is enforced so the inverse transform is real, and DC/Nyquist components
    are kept real (zero phase) since they have no meaningful phase."""
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
    """Percentile-bootstrap CI for statistic(*arrays), row-wise PAIRED
    resampling shared across every array in `arrays` (the same random row
    indices are applied to all of them) so paired structure -- e.g. a storm
    day's (|near_index|, |local_anomaly_index|) pair -- survives each
    replicate instead of shuffling one series against an independently-
    resampled other. A single-array call degenerates to an ordinary
    bootstrap, which is what superposed_epoch_analysis.py's hand-written
    per-event stack bootstrap already does (n_events x n_lags matrix,
    statistic=lambda m: np.nanmean(m, axis=0)) -- this primitive is general
    enough that a future refactor of that script (and coseismic_stacking_
    analysis.py, same pattern) could call it directly with arrays=(M,), no
    interface change needed.

    `statistic` may return a scalar or a fixed-shape ndarray; percentiles
    are then taken elementwise along the bootstrap axis via nanpercentile,
    so a replicate where `statistic` returns NaN (e.g. a degenerate
    resample with a zero denominator) is tolerated rather than poisoning
    the whole CI -- callers should still check ci_lo/ci_hi for NaN before
    trusting them, since an all-NaN bootstrap means the statistic is
    undefined for this data.

    Returns {"point_estimate": statistic(*arrays) on the REAL, un-resampled
    data (not a bootstrap value), "ci_lo", "ci_hi" (at the `ci` level, e.g.
    0.90 -> 5th/95th percentile), "ci_level": ci, "n_bootstrap": n_bootstrap},
    all cast to plain float/list so callers can drop them straight into a
    JSON report. Does not return the raw bootstrap array -- a caller that
    needs the full distribution should run its own loop; no current caller
    needs it."""
    n = arrays[0].shape[0]
    if any(a.shape[0] != n for a in arrays):
        raise ValueError("all arrays must share length along axis 0 for paired resampling")

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
    """Compact JSON-friendly summary of a surrogate statistic's distribution
    for report plotting: binned histogram + headline percentiles, instead of
    dumping thousands of raw floats into report_data.json."""
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
    """Sanity-checks bootstrap_ci() against synthetic data with a known
    answer, covering the two shapes this module's callers actually need:
    a paired-ratio statistic over two same-length arrays (verify_pipeline.
    py's check_storm_cancellation use), and a single-array per-lag mean
    (superposed_epoch_analysis.py's existing hand-written bootstrap, which
    this function could replace without an interface change)."""
    rng = np.random.default_rng(SELF_TEST_SEED)
    ok = True

    # Case A: paired ratio statistic, large n -- CI should tightly bracket
    # the true ratio.
    def _median_ratio(a, b):
        med_a = np.median(a)
        return float(np.median(b) / med_a) if med_a > 0 else float("nan")

    # 20% relative noise on the ratio itself (not a tiny additive term) so
    # the statistic's sampling variance -- and its shrinkage with n -- is
    # large enough to detect reliably from a single fixed-seed draw.
    true_ratio = 0.4
    n_large = 200
    near = rng.uniform(1, 5, size=n_large)
    local = near * true_ratio * (1 + rng.normal(0, 0.2, size=n_large))
    res_large = bootstrap_ci((near, local), _median_ratio, 2000, rng, ci=0.90)
    close_to_truth = abs(res_large["point_estimate"] - true_ratio) < 0.1
    status = "PASS" if close_to_truth else "FAIL"
    if status == "FAIL":
        ok = False
    print(f"[self-test] paired ratio n={n_large}: point={res_large['point_estimate']:.3f} "
          f"ci=[{res_large['ci_lo']:.3f}, {res_large['ci_hi']:.3f}] (true={true_ratio})  {status}")

    # Case A continued: same true ratio/noise, tiny n -- CI should be
    # visibly wider than the large-n case (this is the property the
    # low-bootstrap-power warning in verify_pipeline.py relies on).
    n_small = 5
    near_s = rng.uniform(1, 5, size=n_small)
    local_s = near_s * true_ratio * (1 + rng.normal(0, 0.2, size=n_small))
    res_small = bootstrap_ci((near_s, local_s), _median_ratio, 2000, rng, ci=0.90)
    width_large = res_large["ci_hi"] - res_large["ci_lo"]
    width_small = res_small["ci_hi"] - res_small["ci_lo"]
    status = "PASS" if width_small > width_large else "FAIL"
    if status == "FAIL":
        ok = False
    print(f"[self-test] paired ratio n={n_small}: ci=[{res_small['ci_lo']:.3f}, {res_small['ci_hi']:.3f}] "
          f"width={width_small:.3f} (expected wider than n={n_large}'s width={width_large:.3f})  {status}")

    # Case B: single-array per-lag mean, mirroring superposed_epoch_
    # analysis.py's run_band -- M is (n_events x n_lags), statistic reduces
    # over events, CI should bracket the known injected per-lag mean. Uses a
    # high (0.999) CI level here specifically to make a single fixed-seed
    # bracket check robust (a correctly-calibrated 90% CI is *expected* to
    # miss its target ~10% of the time per lag -- that's not a bug -- so
    # asserting a 90% CI always brackets truth would itself be a flaky
    # test; 0.999 makes a spurious miss astronomically unlikely instead).
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
    print(f"[self-test] single-array per-lag mean: ci_lo={np.round(ci_lo, 2)} "
          f"ci_hi={np.round(ci_hi, 2)} (true={lag_means})  {status}")

    return ok


if __name__ == "__main__":
    import sys

    sys.exit(0 if self_test() else 1)
