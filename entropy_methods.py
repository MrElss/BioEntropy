import numpy as np
from scipy.spatial.distance import cdist
from scipy.special import logsumexp, ndtri
from scipy.stats import qmc


def fixed_bandwidth(n, p):
    if n < 2 or p < 1:
        raise ValueError("固定平滑尺度需要至少2个完整样本及1个指标。")
    return float((4 / (p + 2)) ** (1 / (p + 4)) * n ** (-1 / (p + 4)))


def renyi2(x, h):
    n, p = x.shape
    logs = -cdist(x, x, "sqeuclidean") / (4 * h * h)
    total = logsumexp(logs)
    return (
        float(0.5 * p * np.log(4 * np.pi * h * h) - total + 2 * np.log(n)),
        float(np.exp(np.log(n) - total)),
    )


def shannon_clinical(x, h, power=15, reps=4, seed=20261029):
    n, p = x.shape
    constant = -p * np.log(h * np.sqrt(2 * np.pi)) - np.log(n)
    answers = []
    for r in range(reps):
        u = qmc.Sobol(p + 1, scramble=True, seed=seed + r).random_base2(power)
        component = np.minimum((u[:, 0] * n).astype(int), n - 1)
        samples = x[component] + h * ndtri(np.clip(u[:, 1:], 1e-15, 1 - 1e-15))
        total = 0.0
        for start in range(0, len(samples), 256):
            logs = -cdist(samples[start : start + 256], x, "sqeuclidean") / (2 * h * h)
            total -= (constant + logsumexp(logs, axis=1)).sum()
        answers.append(total / len(samples))
    return (float(np.mean(answers)), float(np.std(answers, ddof=1) / np.sqrt(reps)))


def shannon_animal(x, h, power=16, reps=8, seed=20261030):
    n, p = x.shape
    analytic = 0.5 * p * np.log(2 * np.pi * np.e * h * h) + np.log(n)
    answers = []
    for r in range(reps):
        u = qmc.Sobol(p, scramble=True, seed=seed + r).random_base2(power)
        noise = ndtri(np.clip(u, 1e-15, 1 - 1e-15))
        overlap = 0.0
        for i in range(n):
            delta = x[i] - x
            logits = -(np.sum(delta * delta, axis=1)[None, :] + 2 * h * (noise @ delta.T)) / (
                2 * h * h
            )
            overlap += float(logsumexp(logits, axis=1).mean())
        answers.append(analytic - overlap / n)
    return (float(np.mean(answers)), float(np.std(answers, ddof=1) / np.sqrt(reps)))


def estimate(matrix, h, mode, method="both"):
    x = np.asarray(matrix, dtype=float)
    if method == "none":
        return {}
    if method not in ("Shannon", "Renyi2", "both"):
        raise ValueError("KDE熵方法选择无效。")
    result = dict(Entropy_status="OK")
    if method in ("Shannon", "both"):
        result.update(Entropy_Shannon=np.nan, Shannon_numeric_SE=np.nan, Shannon_upper=np.nan)
    if method in ("Renyi2", "both"):
        result.update(Entropy_Renyi2=np.nan, Renyi2_upper=np.nan, Renyi2_self_fraction=np.nan)
    if x.ndim != 2 or x.shape[1] == 0 or len(x) < 2:
        result["Entropy_status"] = "NA：完整样本不足2个或没有指标"
        return result
    if not np.isfinite(x).all() or not np.isfinite(h) or h <= 0:
        result["Entropy_status"] = "NA：输入或固定尺度无效"
        return result
    if mode not in ("clinical", "animal"):
        raise ValueError("未知实验类型")
    n, p = x.shape
    if method in ("Shannon", "both"):
        hs, se = (shannon_clinical if mode == "clinical" else shannon_animal)(x, h)
        result.update(
            Entropy_Shannon=hs,
            Shannon_numeric_SE=se,
            Shannon_upper=float(0.5 * p * np.log(2 * np.pi * np.e * h * h) + np.log(n)),
        )
    if method in ("Renyi2", "both"):
        hr, fraction = renyi2(x, h)
        result.update(
            Entropy_Renyi2=hr,
            Renyi2_upper=float(0.5 * p * np.log(4 * np.pi * h * h) + np.log(n)),
            Renyi2_self_fraction=fraction,
        )
    if not np.isfinite([v for k, v in result.items() if k != "Entropy_status"]).all():
        raise ValueError("所选熵计算出现非有限值，停止输出。")
    return result
