"""Numerical primitives: entropy estimation, reference-state entropy, bootstrap
confidence intervals, multi-method outlier detection, HDI, feature direction."""
from __future__ import annotations

from typing import Any, Dict, Iterable, List, Tuple
import re

import numpy as np
import pandas as pd
from scipy.stats import chi2
from sklearn.covariance import LedoitWolf, MinCovDet
from sklearn.ensemble import IsolationForest
from sklearn.neighbors import LocalOutlierFactor

from bioentropy_constants import (
    RNG_SEED,
    ENTROPY_STRONG_OUTLIER_METHODS,
    ENTROPY_OUTLIER_SHRINK,
)

__all__ = [
    "fill_nan_with_col_median",
    "compute_group_entropy",
    "compute_group_entropy_with_diagnostics",
    "prepare_entropy_matrix",
    "compute_group_entropy_robust",
    "compute_distance_state_entropy",
    "bootstrap_distance_state_entropy_ci",
    "compute_hdi",
    "bootstrap_entropy_ci",
    "bootstrap_mean_ci",
    "bootstrap_diff_ci",
    "detect_outliers_methods",
    "benjamini_hochberg",
    "infer_direction",
]


def benjamini_hochberg(pvalues: Iterable[float]) -> np.ndarray:
    """Benjamini-Hochberg FDR-adjusted q-values.

    Accepts NaNs (positions excluded from the rank count and returned as NaN).
    Returns an array aligned to the input order with values in [0, 1].
    """
    p = np.asarray(list(pvalues), dtype=float)
    q = np.full(p.shape, np.nan)
    finite = np.where(np.isfinite(p))[0]
    m = finite.size
    if m == 0:
        return q
    sub = p[finite]
    order = np.argsort(sub)
    ranked = sub[order]
    adj = ranked * m / (np.arange(1, m + 1))
    # enforce monotonic non-decreasing q-values from the largest p downward
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0.0, 1.0)
    out = np.empty(m)
    out[order] = adj
    q[finite] = out
    return q


def fill_nan_with_col_median(X: np.ndarray) -> np.ndarray:
    X = np.asarray(X, dtype=float)
    if X.size == 0:
        return X
    col_median = np.nanmedian(X, axis=0)
    inds = np.where(np.isnan(X))
    if len(inds[0]) > 0:
        X = X.copy()
        X[inds] = np.take(col_median, inds[1])
    return X


def compute_group_entropy(X: np.ndarray) -> float:
    return compute_group_entropy_with_diagnostics(X)[0]


def _entropy_reliability(n_over_p: float) -> str:
    """Qualitative reliability of the log-determinant covariance entropy.

    The estimator needs the sample covariance to be well conditioned. With
    n <= p the sample covariance is singular and the log-det is driven almost
    entirely by the Ledoit-Wolf shrinkage target rather than the data, so the
    cross-group entropy comparison is fragile. Thresholds follow the common
    n/p rules of thumb for stable covariance estimation.
    """
    if not np.isfinite(n_over_p) or n_over_p <= 0:
        return "undefined"
    if n_over_p >= 5.0:
        return "reliable"
    if n_over_p >= 2.0:
        return "limited"
    return "low"


def compute_group_entropy_with_diagnostics(X: np.ndarray) -> Tuple[float, Dict[str, Any]]:
    """Log-determinant Gaussian entropy plus estimator-reliability diagnostics.

    Returns (entropy, diagnostics) where diagnostics reports the sample size n,
    the feature dimension p, the n/p ratio, the Ledoit-Wolf shrinkage intensity
    actually used, and a qualitative reliability label. These let downstream
    tables and the report flag groups where n <= p (covariance singular /
    over-shrunk). The numerical entropy value is identical to before.
    """
    X = fill_nan_with_col_median(X)
    n = int(len(X))
    p = int(X.shape[1]) if X.ndim == 2 else 0
    n_over_p = float(n / p) if p > 0 else float("nan")
    diag: Dict[str, Any] = {
        "entropy_n": n,
        "entropy_p": p,
        "entropy_n_over_p": n_over_p,
        "entropy_ledoitwolf_shrinkage": float("nan"),
        "entropy_cov_reliability": _entropy_reliability(n_over_p),
    }
    if n < 2:
        diag["entropy_cov_reliability"] = "undefined"
        return float("nan"), diag
    lw = LedoitWolf().fit(X)
    diag["entropy_ledoitwolf_shrinkage"] = float(getattr(lw, "shrinkage_", float("nan")))
    cov = lw.covariance_
    sign, logdet = np.linalg.slogdet(cov)
    if sign <= 0:
        eigvals = np.linalg.eigvalsh(cov)
        eigvals = eigvals[eigvals > 1e-10]
        if len(eigvals) == 0:
            return float("nan"), diag
        logdet = np.log(eigvals).sum()
    k = X.shape[1]
    return float(0.5 * logdet + 0.5 * k * np.log(2 * np.pi * np.e)), diag


def prepare_entropy_matrix(
    X: np.ndarray,
    min_methods: int = ENTROPY_STRONG_OUTLIER_METHODS,
    shrink: float = ENTROPY_OUTLIER_SHRINK,
) -> Tuple[np.ndarray, Dict[str, Any]]:
    """Attenuate only high-confidence multivariate outliers before entropy.

    The QC table still reports the broader consensus rule (>=2 methods).  Entropy
    uses a stricter rule so isolated biological extremes are not over-removed.
    """
    X = fill_nan_with_col_median(np.asarray(X, dtype=float))
    n = len(X)
    if n == 0:
        return X, {
            "n_entropy_outliers": 0,
            "pct_entropy_outliers": 0.0,
            "entropy_outlier_rule": f">={min_methods} of 5 methods, shrink={shrink}",
        }

    flags = detect_outliers_methods(X)
    if not flags:
        method_count = np.zeros(n, dtype=int)
    else:
        method_count = np.column_stack([flags[k] for k in flags]).sum(axis=1)
    strong = method_count >= min_methods

    X_adj = X.copy()
    if strong.any():
        non_strong = X_adj[~strong]
        center_source = non_strong if len(non_strong) >= max(3, X_adj.shape[1] // 2) else X_adj
        center = np.nanmedian(center_source, axis=0)
        center = np.where(np.isfinite(center), center, 0.0)
        X_adj[strong] = center.reshape(1, -1) + shrink * (X_adj[strong] - center.reshape(1, -1))

    return X_adj, {
        "n_entropy_outliers": int(strong.sum()),
        "pct_entropy_outliers": float(100 * strong.mean()) if n > 0 else 0.0,
        "entropy_outlier_rule": f">={min_methods} of 5 methods, shrink={shrink}",
    }


def compute_group_entropy_robust(X: np.ndarray) -> Tuple[float, np.ndarray, Dict[str, Any]]:
    X_adj, info = prepare_entropy_matrix(X)
    entropy, diag = compute_group_entropy_with_diagnostics(X_adj)
    info.update(diag)
    return entropy, X_adj, info


def compute_distance_state_entropy(distance_values: np.ndarray) -> float:
    """Entropy-like state disorder from distances to a healthy reference.

    Uses the second moment around zero distance so both group-level displacement
    from Normal and within-group spread contribute to the score.
    """
    d = pd.to_numeric(pd.Series(distance_values), errors="coerce").dropna().to_numpy(dtype=float)
    if len(d) < 2:
        return float("nan")
    d = np.clip(d, 0, None)
    second_moment = float(np.mean(d * d))
    if not np.isfinite(second_moment) or second_moment <= 0:
        return float("nan")
    return float(0.5 * np.log(second_moment) + 0.5 * np.log(2 * np.pi * np.e))


def bootstrap_distance_state_entropy_ci(
    distance_values: np.ndarray, n_boot: int = 200, seed: int = RNG_SEED
) -> Tuple[float, float, float]:
    d = pd.to_numeric(pd.Series(distance_values), errors="coerce").dropna().to_numpy(dtype=float)
    if len(d) < 2:
        return float("nan"), float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    vals = []
    n = len(d)
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        vals.append(compute_distance_state_entropy(d[idx]))
    vals = np.asarray(vals, dtype=float)
    vals = vals[np.isfinite(vals)]
    if len(vals) == 0:
        return float("nan"), float("nan"), float("nan")
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5)), float(vals.std(ddof=1))


def compute_hdi(z_df: pd.DataFrame, features: List[str]) -> pd.Series:
    directions = np.array([infer_direction(f) for f in features], dtype=float)
    Z = z_df[features].to_numpy(dtype=float)
    Z = np.where(np.isnan(Z), 0.0, Z)
    return pd.Series(np.nanmean(Z * directions, axis=1), index=z_df.index, name="HDI")


def bootstrap_entropy_ci(X: np.ndarray, n_boot: int = 200, seed: int = RNG_SEED) -> Tuple[float, float, float]:
    X = np.asarray(X)
    if len(X) < 2:
        return float("nan"), float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    n = len(X)
    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        vals.append(compute_group_entropy(X[idx].copy()))
    vals = np.asarray(vals, dtype=float)
    lo, hi = np.nanquantile(vals, [0.025, 0.975])
    return float(lo), float(hi), float(np.nanstd(vals, ddof=1))


def bootstrap_mean_ci(x: np.ndarray, n_boot: int = 200, seed: int = RNG_SEED) -> Tuple[float, float, float]:
    x = np.asarray(x, dtype=float)
    if len(x) < 2:
        if len(x) == 1:
            return float(x[0]), float(x[0]), 0.0
        return float("nan"), float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    n = len(x)
    idx = rng.integers(0, n, size=(n_boot, n))
    vals = x[idx].mean(axis=1)
    lo, hi = np.quantile(vals, [0.025, 0.975])
    return float(lo), float(hi), float(vals.std(ddof=1))


def bootstrap_diff_ci(x1: np.ndarray, x2: np.ndarray, func, n_boot: int = 400, seed: int = RNG_SEED) -> Tuple[float, float, float]:
    x1 = np.asarray(x1)
    x2 = np.asarray(x2)
    if len(x1) == 0 or len(x2) == 0:
        return float("nan"), float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    n1, n2 = len(x1), len(x2)
    vals = []
    for _ in range(n_boot):
        b1 = x1[rng.integers(0, n1, size=n1)]
        b2 = x2[rng.integers(0, n2, size=n2)]
        vals.append(func(b1) - func(b2))
    vals = np.asarray(vals, dtype=float)
    lo, hi = np.nanquantile(vals, [0.025, 0.975])
    p_boot = 2 * min((vals <= 0).mean(), (vals >= 0).mean())
    return float(lo), float(hi), float(p_boot)


def detect_outliers_methods(X: np.ndarray) -> Dict[str, np.ndarray]:
    X = fill_nan_with_col_median(X)
    n, p = X.shape
    flags: Dict[str, np.ndarray] = {}

    if n == 0:
        zero = np.zeros(0, dtype=bool)
        return {"robust_z": zero, "iqr": zero, "iso": zero, "lof": zero, "rmd": zero}

    med = np.median(X, axis=0)
    mad = np.median(np.abs(X - med), axis=0)
    mad[mad < 1e-9] = np.nan
    rz = 0.6745 * (X - med) / mad
    rz = np.where(np.isnan(rz), 0.0, rz)
    flags["robust_z"] = np.any(np.abs(rz) > 4, axis=1)

    q1 = np.quantile(X, 0.25, axis=0)
    q3 = np.quantile(X, 0.75, axis=0)
    iqr = q3 - q1
    lower = q1 - 1.5 * iqr
    upper = q3 + 1.5 * iqr
    flags["iqr"] = np.any((X < lower) | (X > upper), axis=1)

    contamination = min(0.05, max(0.01, 5 / max(n, 1)))

    try:
        iso = IsolationForest(contamination=contamination, random_state=RNG_SEED)
        flags["iso"] = iso.fit_predict(X) == -1 if n >= 5 else np.zeros(n, dtype=bool)
    except Exception:
        flags["iso"] = np.zeros(n, dtype=bool)

    try:
        k = min(20, max(2, n // 3))
        flags["lof"] = (
            LocalOutlierFactor(n_neighbors=k, contamination=contamination).fit_predict(X) == -1
            if n >= 5 else np.zeros(n, dtype=bool)
        )
    except Exception:
        flags["lof"] = np.zeros(n, dtype=bool)

    try:
        if n >= max(5, p + 1):
            mcd = MinCovDet(random_state=RNG_SEED).fit(X)
            md2 = mcd.mahalanobis(X)
            cutoff = chi2.ppf(0.999, df=p)
            flags["rmd"] = md2 > cutoff
        else:
            flags["rmd"] = np.zeros(n, dtype=bool)
    except Exception:
        flags["rmd"] = np.zeros(n, dtype=bool)

    return flags


def infer_direction(feature_name: str) -> float:
    """
    Direction rule for disease burden aggregation.

    Return +1 when a higher value usually indicates worse status;
    return -1 when a higher value usually indicates better status.

    Important bug fix:
    - '非高密度脂蛋白' / 'non-HDL' must be treated as harmful, not beneficial.
    - Harmful keywords are checked before beneficial HDL keywords to avoid
      false matches caused by substring overlap.
    """
    name = str(feature_name).strip().lower()
    compact = re.sub(r"\s+", "", name)

    harmful_keywords = [
        "非高密度脂蛋白", "non-hdl", "nonhdl",
        "极低密度脂蛋白", "vldl",
        "低密度脂蛋白", "ldl-c", "ldl",
        "总胆固醇", "tc",
        "甘油三酯", "tg",
        "空腹血浆血糖", "空腹血糖", "血糖", "glucose", "glu",
        "糖化血红蛋白", "hba1c", "糖化",
        "c反应蛋白", "hscrp", "crp",
        "alt", "ast", "ggt",
    ]
    beneficial_keywords = [
        "高密度脂蛋白胆固醇", "高密度脂蛋白", "hdl-c",
        "apoa1", "载脂蛋白a1", "脂联素", "adiponectin",
    ]

    if any(kw in compact for kw in harmful_keywords):
        return 1.0
    if any(kw in compact for kw in beneficial_keywords):
        return -1.0
    return 1.0
