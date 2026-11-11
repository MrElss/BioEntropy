import numpy as np
import pandas as pd
from scipy.special import logsumexp, ndtri
from scipy.spatial.distance import cdist
from scipy.stats import qmc


def component_scores(x, h, mode, power=None, reps=None):
    x = np.asarray(x, float)
    n, p = x.shape
    if n < 2 or p < 1 or (not np.isfinite(x).all()) or (not np.isfinite(h)) or (h <= 0):
        return (np.full(n, np.nan), np.full(n, np.nan), np.zeros(n, int))
    clinical = mode == "clinical"
    power = power if power is not None else 15 if clinical else 16
    reps = reps if reps is not None else 4 if clinical else 8
    values = np.full((reps, n), np.nan)
    counts = np.zeros(n, int)
    if clinical:
        constant = -p * np.log(h * np.sqrt(2 * np.pi)) - np.log(n)
        for r in range(reps):
            u = qmc.Sobol(p + 1, scramble=True, seed=20261029 + r).random_base2(power)
            owner = np.minimum((u[:, 0] * n).astype(int), n - 1)
            samples = x[owner] + h * ndtri(np.clip(u[:, 1:], 1e-15, 1 - 1e-15))
            information = np.empty(len(samples))
            for start in range(0, len(samples), 256):
                logs = -cdist(samples[start : start + 256], x, "sqeuclidean") / (2 * h * h)
                information[start : start + 256] = -(constant + logsumexp(logs, axis=1))
            nr = np.bincount(owner, minlength=n)
            values[r] = np.divide(
                np.bincount(owner, weights=information, minlength=n),
                nr,
                out=np.full(n, np.nan),
                where=nr > 0,
            )
            counts += nr
    else:
        analytic = 0.5 * p * np.log(2 * np.pi * np.e * h * h) + np.log(n)
        for r in range(reps):
            noise = ndtri(
                np.clip(
                    qmc.Sobol(p, scramble=True, seed=20261030 + r).random_base2(power),
                    1e-15,
                    1 - 1e-15,
                )
            )
            for i in range(n):
                delta = x[i] - x
                logits = -(np.sum(delta * delta, axis=1)[None, :] + 2 * h * (noise @ delta.T)) / (
                    2 * h * h
                )
                values[r, i] = analytic - logsumexp(logits, axis=1).mean()
            counts += len(noise)
    return (values.mean(0), values.std(0, ddof=1) / np.sqrt(reps), counts)


def add_scores(result, progress=print):
    rows = result["observations"]
    z = result["_z"]
    if len(rows) != len(z):
        raise ValueError("样本信息分数与原始记录数不一致。")
    obs = pd.DataFrame(rows)
    weeks = obs.Week.to_numpy() if "Week" in obs else np.zeros(len(obs), int)
    checks = []
    for row in rows:
        row.update(
            Shannon_information_score=np.nan,
            Shannon_weighted_contribution=np.nan,
            Score_numeric_SE=np.nan,
            Score_quadrature_nodes=0,
            Score_status="NA：未纳入Shannon计算或组内完整样本不足2个",
        )
    for group in result["summary"]:
        mask = (
            (obs.Group.to_numpy() == group["Group"])
            & (weeks == group["Week"])
            & obs.Included_entropy.to_numpy(bool)
            & np.isfinite(z).all(1)
        )
        indices = np.flatnonzero(mask)
        if len(indices) != group["N_entropy"]:
            raise ValueError("散点与组熵的纳入样本不一致。")
        if len(indices) < 2:
            continue
        progress(
            f"{result['name']} {group['Group']} {group['Week']}W：计算样本信息分数（{len(indices)}个）"
        )
        scores, se, counts = component_scores(z[mask], group["Bandwidth"], result["mode"])
        for i, score, error, count in zip(indices, scores, se, counts):
            rows[i].update(
                Shannon_information_score=float(score),
                Shannon_weighted_contribution=float(score / len(indices)),
                Score_numeric_SE=float(error),
                Score_quadrature_nodes=int(count),
                Score_status="OK" if np.isfinite(score) else "NA：部分积分重复未覆盖该样本核",
            )
        checks.append(
            dict(
                Experiment=result["name"],
                Group=group["Group"],
                Week=group["Week"],
                N_entropy=len(indices),
                N_scores=int(np.isfinite(scores).sum()),
                Original_Shannon=group["Entropy_Shannon"],
                Mean_information_score=float(scores.mean()),
                Numerical_difference=float(scores.mean() - group["Entropy_Shannon"]),
                Note="数学上均值等于组Shannon；有限数值积分可有小误差；样本分数采用原始计算值",
            )
        )
    result["score_checks"] = checks
