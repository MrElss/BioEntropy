import numpy as np
from scipy.spatial.distance import cdist, pdist


def fill_group(x):
    x = np.asarray(x, float)
    keep = np.isfinite(x).any(1)
    a = x[keep]
    if len(a) < 2 or not np.isfinite(a).any(0).all():
        return (None, keep)
    median = np.nanmedian(np.where(np.isfinite(a), a, np.nan), axis=0)
    return (np.where(np.isfinite(a), a, median), keep)


def energy_u(x, reference, scores=False):
    x, reference = (np.asarray(x, float), np.asarray(reference, float))
    if x.ndim != 2 or reference.ndim != 2 or x.shape[1] != reference.shape[1]:
        raise ValueError("能量距离要求两组使用相同的多指标面板。")
    if (
        min(len(x), len(reference)) < 2
        or x.shape[1] == 0
        or (not np.isfinite(x).all())
        or (not np.isfinite(reference).all())
    ):
        return (np.nan, np.full(len(x), np.nan)) if scores else np.nan
    within = 2 * pdist(x).sum() / (len(x) * (len(x) - 1))
    normal = 2 * pdist(reference).sum() / (len(reference) * (len(reference) - 1))
    cross = cdist(x, reference)
    value = float(2 * cross.mean() - within - normal)
    if not scores:
        return value
    contribution = 2 * cross.mean(1) - cdist(x, x).sum(1) / (len(x) - 1) - normal
    return (value, contribution)


def prepare_reference(reference):
    raw = np.asarray(reference, float)
    if raw.ndim != 2 or not raw.shape[1]:
        return (None, None, None, "NA：没有匹配的Normal指标")
    if len(raw) < 2 or (np.isfinite(raw).sum(0) < 3).any():
        return (None, None, None, "NA：Normal至少一个指标有效值不足3个")
    sd = np.nanstd(raw, axis=0, ddof=1)
    if not np.isfinite(sd).all() or (sd <= 0).any():
        return (None, None, None, "NA：Normal标准差无效或为零")
    filled, _ = fill_group(raw)
    if filled is None:
        return (None, None, None, "NA：Normal可用记录不足或存在整列缺失")
    mu = filled.mean(0)
    return ((filled - mu) / sd, mu, sd, "OK")


def labelings(raw, reference, assignments):
    raw = np.asarray(raw, float)
    a = np.asarray(assignments, float)
    output = np.full((a.shape[1], 2), np.nan)
    if len(reference) < 2 or not np.isfinite(reference).all():
        return output
    if np.isfinite(raw).all():
        d = cdist(raw, raw)
        cross = cdist(raw, reference).mean(1)
        q0 = 2 * pdist(reference).sum() / (len(reference) * (len(reference) - 1))
        for start in range(0, a.shape[1], 512):
            for gi, weights in enumerate(
                (a[:, start : start + 512], 1 - a[:, start : start + 512])
            ):
                n = weights.sum(0)
                good = n >= 2
                v = np.full(len(n), np.nan)
                v[good] = (
                    2 * (cross @ weights[:, good]) / n[good]
                    - (weights[:, good] * (d @ weights[:, good])).sum(0) / (n[good] * (n[good] - 1))
                    - q0
                )
                output[start : start + len(n), gi] = v
    else:
        for j in range(a.shape[1]):
            for gi, mask in enumerate((a[:, j].astype(bool), ~a[:, j].astype(bool))):
                x, _ = fill_group(raw[mask])
                if x is not None:
                    output[j, gi] = energy_u(x, reference)
    return output
