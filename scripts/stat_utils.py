"""Shared statistics helpers for the professor-suggested validation methods
(`~/.claude/plans/block-bootstrap-shiny-pearl.md`): the robust median/MAD
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
