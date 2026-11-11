from report_text import PERMUTATION_NOTE
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.special import ndtri
from scipy.stats import qmc, beta
from threadpoolctl import threadpool_limits
from metric_selection import metrics as result_metrics
from energy_distance import labelings as energy_labelings

METRICS = ("NMD", "Entropy_Shannon", "Entropy_Renyi2")
SEED = 20261030


def assignments(groups, strata, b):
    groups = np.asarray(groups, bool)
    strata = np.asarray(strata, str)
    rng = np.random.default_rng(SEED)
    a = np.empty((len(groups), b + 1), float)
    a[:, 0] = groups
    for stratum in np.unique(strata):
        idx = np.flatnonzero(strata == stratum)
        for j in range(1, b + 1):
            a[idx, j] = rng.permutation(groups[idx])
    return a


def holm(p):
    p = np.asarray(p, float)
    out = np.full(len(p), np.nan)
    ix = np.flatnonzero(np.isfinite(p))
    ix = ix[np.argsort(p[ix])]
    if len(ix):
        out[ix] = np.minimum(1.0, np.maximum.accumulate((len(p) - np.arange(len(ix))) * p[ix]))
    return out


def entropy_labelings(x, h, a, power, selected=None):
    n, p = x.shape
    na = a.sum(0)
    nb = n - na
    selected = METRICS[1:] if selected is None else selected
    result = {key: np.full((a.shape[1], 2), np.nan) for key in selected}
    if n < 4 or np.min(na) < 2 or np.min(nb) < 2:
        return result
    if "Entropy_Shannon" in selected:
        q = 2**power
        noise = ndtri(
            np.clip(
                qmc.Sobol(p, scramble=True, seed=20261101).random_base2(power), 1e-15, 1 - 1e-15
            )
        )
        own = -np.sum(noise * noise, axis=1) / 2
        sums = np.zeros((a.shape[1], 2))
        owner_chunk = max(1, min(32, int(8000000 / (q * n))))
        for start in range(0, n, owner_chunk):
            end = min(n, start + owner_chunk)
            samples = (x[start:end, None, :] + h * noise[None, :, :]).reshape(-1, p)
            logs = -cdist(samples, x, "sqeuclidean") / (2 * h * h)
            rowmax = logs.max(1)
            kernels = np.exp(logs - rowmax[:, None])
            offset = rowmax - np.tile(own, end - start)
            del samples, logs
            for col in range(0, a.shape[1], 256):
                labels = a[:, col : col + 256]
                for gi, weights in enumerate((labels, 1 - labels)):
                    density = kernels @ weights
                    overlap = (
                        (np.log(np.maximum(density, np.finfo(float).tiny)) + offset[:, None])
                        .reshape(end - start, q, -1)
                        .mean(1)
                    )
                    sums[col : col + weights.shape[1], gi] += (overlap * weights[start:end]).sum(0)
        hk = 0.5 * p * np.log(2 * np.pi * np.e * h * h)
        result["Entropy_Shannon"][:, 0] = hk + np.log(na) - sums[:, 0] / na
        result["Entropy_Shannon"][:, 1] = hk + np.log(nb) - sums[:, 1] / nb
    if "Entropy_Renyi2" in selected:
        kernel2 = np.exp(-cdist(x, x, "sqeuclidean") / (4 * h * h))
        for col in range(0, a.shape[1], 512):
            labels = a[:, col : col + 512]
            for gi, weights in enumerate((labels, 1 - labels)):
                total = (weights * (kernel2 @ weights)).sum(0)
                result["Entropy_Renyi2"][col : col + weights.shape[1], gi] = (
                    0.5 * p * np.log(4 * np.pi * h * h) + 2 * np.log(weights.sum(0)) - np.log(total)
                )
    return result


def permutation_p(values):
    values = np.asarray(values, float)
    if len(values) < 2 or not np.isfinite(values).all():
        return (np.nan, np.nan, np.nan, 0)
    b = len(values) - 1
    extreme = int((np.abs(values[1:]) >= abs(values[0]) - 1e-12).sum())
    p = (extreme + 1) / (b + 1)
    lo = 0.0 if extreme == 0 else float(beta.ppf(0.025, extreme, b - extreme + 1))
    hi = 1.0 if extreme == b else float(beta.ppf(0.975, extreme + 1, b - extreme))
    return (p, lo, hi, extreme)


def correction(rows, selected=None):
    for contrast in ("raw",):
        families = (
            [(tuple(selected), "NMD+selected_method")]
            if selected is not None
            else [(METRICS[:2], "NMD+Shannon"), ((METRICS[2],), "Renyi2_supplemental")]
        )
        for metrics, family in families:
            subset = [r for r in rows if r["Contrast"] == contrast and r["Metric"] in metrics]
            for key in ("P_raw", "P_MC_lower", "P_MC_upper"):
                for r, value in zip(subset, holm([r[key] for r in subset])):
                    r[
                        {
                            "P_raw": "P_Holm",
                            "P_MC_lower": "P_Holm_MC_lower",
                            "P_MC_upper": "P_Holm_MC_upper",
                        }[key]
                    ] = float(value)
                    r["Family"] = contrast + " / " + family
                    r["Family_size"] = len(subset)
    return rows


def stars(p):
    return (
        "NA"
        if not np.isfinite(p)
        else "***"
        if p < 0.001
        else "**"
        if p < 0.01
        else "*"
        if p < 0.05
        else "ns"
    )


def compare(result, b=9999, progress=print):
    if b not in (0, 999, 9999):
        raise ValueError("置换次数应为0、999或9999。")
    active = result_metrics(result)
    selected = active[1:]
    has_shannon = "Entropy_Shannon" in selected
    result["inference_settings"] = dict(
        permutations=b,
        seed=SEED,
        quadrature_seed=20261101,
        null="within-trial exchangeability / no group effect on multivariate records",
        clinical_blocks="entire patient trajectories, including missingness, within Trial",
        selected_roles=result.get("roles", {}),
        raw_family="NMD + selected method, across visits of selected pair",
        comparison="raw group difference at each visit, including baseline",
        baseline_alignment=False,
        selected_metrics=active,
        unselected_methods="selected method only",
        refinement="Shannon Holm P in [0.025,0.10]: refine only that visit to 4x quadrature once; the refined estimate is final",
        note=PERMUTATION_NOTE,
    )
    obs = pd.DataFrame(result["observations"]).copy()
    if "Week" not in obs:
        obs["Week"] = 0
    roles = result.get("roles", {})
    treatment, control = (roles.get("treatment"), roles.get("comparator"))
    if treatment not in set(obs.Group) or control not in set(obs.Group) or treatment == control:
        result["inference_settings"]["status"] = "NA：未选择两个不同的治疗组和比较组"
        return []
    selected = obs.Group.isin([treatment, control]).to_numpy()
    if "ID_verified" in obs:
        selected &= obs.ID_verified.to_numpy(bool)
    if result.get("shared_shannon_sample"):
        selected &= obs.Included_entropy.to_numpy(bool)
    selected_rows = selected.copy()
    selected = active[1:]
    z = result["_z"][selected_rows]
    df = obs.loc[selected_rows].reset_index(drop=True)
    if result["mode"] == "clinical":
        if "USUBJID" not in df or df.USUBJID.fillna("").eq("").any():
            result["inference_settings"]["status"] = "NA：缺少可靠患者编号，未进行独立访视置换"
            return []
        trial = df.Trial.astype(str) if "Trial" in df else pd.Series([result["name"]] * len(df))
        keys = list(zip(trial, df.USUBJID.astype(str)))
    else:
        trial = pd.Series(["animal"] * len(df))
        keys = [("animal", str(i)) for i in range(len(df))]
    units = list(dict.fromkeys(keys))
    lookup = {key: i for i, key in enumerate(units)}
    unit_index = np.array([lookup[key] for key in keys], int)
    labels = np.zeros(len(units), bool)
    strata = np.empty(len(units), object)
    for key, sub in df.assign(_unit=unit_index, _stratum=trial).groupby("_unit"):
        if sub.Group.nunique() != 1 or sub.Week.duplicated().any():
            raise ValueError("置换核查：患者分组变化或同次访视重复。")
        labels[key] = sub.Group.iloc[0] == treatment
        strata[key] = sub["_stratum"].iloc[0]
    exchangeable = any(len(set(labels[strata == s])) == 2 for s in np.unique(strata))
    a = assignments(labels, strata, b)
    summary = pd.DataFrame(result["summary"]).set_index(["Group", "Week"])
    weeks = sorted({w for group, w in summary.index if group in (treatment, control)})
    values, checks = ({}, [])
    powers = {}
    h_values = (
        {float(x["Bandwidth"]) for x in result["summary"]}
        if any(k.startswith("Entropy_") for k in selected)
        else {1.0}
    )
    if len(h_values) != 1:
        raise ValueError("组间比较要求相同的固定核带宽。")
    h = h_values.pop()
    no_nmd = roles.get("reference") in (treatment, control)

    def compute_week(w, power=None):
        present = df.Week.eq(w).to_numpy()
        valid = present & np.isfinite(z).all(1) & df.Included_entropy.to_numpy(bool)
        aa = a[unit_index[valid]]
        if power is None:
            power = 12 if result["mode"] == "animal" else 8 if valid.sum() <= 600 else 7
        powers[w] = power
        progress(f"{result['name']} {w}W：{b}次组别置换；仅检验所选方法")
        vv = {k: np.full((b + 1, 2), np.nan) for k in selected}
        if b and exchangeable:
            kde = [k for k in selected if k.startswith("Entropy_")]
            if kde:
                vv.update(entropy_labelings(z[valid], h, aa, power, kde))
            if (
                "Energy_distance" in selected
                and (not no_nmd)
                and (result.get("_energy_reference") is not None)
            ):
                raw = result["_energy_raw"][selected_rows]
                valid_energy = present & np.isfinite(raw).any(1)
                vv["Energy_distance"] = energy_labelings(
                    raw[valid_energy], result["_energy_reference"], a[unit_index[valid_energy]]
                )
        goodn = present & np.isfinite(df.NMD.to_numpy(float))
        an = a[unit_index[goodn]]
        nvalues = df.NMD.to_numpy(float)[goodn]
        counts = an.sum(0)
        other = goodn.sum() - counts
        nv = np.full((b + 1, 2), np.nan)
        if (
            b
            and exchangeable
            and (not no_nmd)
            and len(counts)
            and (min(counts.min(), other.min()) >= 2)
        ):
            nv[:, 0] = nvalues @ an / counts
            nv[:, 1] = nvalues @ (1 - an) / other
        vv["NMD"] = nv
        values[w] = vv
        if has_shannon:
            for group, gi in ((treatment, 0), (control, 1)):
                original = (
                    summary.loc[(group, w), "Entropy_Shannon"]
                    if (group, w) in summary.index
                    else np.nan
                )
                checks.append(
                    dict(
                        Experiment=result["name"],
                        Group=group,
                        Week=w,
                        Nodes_per_kernel=2**power,
                        Original_Shannon=original,
                        Permutation_observed=vv["Entropy_Shannon"][0, gi],
                        Difference=vv["Entropy_Shannon"][0, gi] - original,
                    )
                )

    def rows_from_values():
        rows = []
        for metric in active:
            for w in weeks:
                for contrast in ("raw",):
                    current = values[w][metric][:, 0] - values[w][metric][:, 1]
                    try:
                        original = (
                            summary.loc[(treatment, w), metric] - summary.loc[(control, w), metric]
                        )
                    except KeyError:
                        original = np.nan
                    p, lo, hi, extreme = permutation_p(current)
                    reason = (
                        "用户关闭置换"
                        if not b
                        else "试验内没有可交换的两组患者"
                        if not exchangeable
                        else "正常参照本身属于被比较组；固定参照指标不能用于本置换检验"
                        if metric in ("NMD", "Energy_distance") and no_nmd
                        else "至少一次置换后可用样本不足2个或原始估计不可用"
                    )
                    if not np.isfinite(original):
                        p = lo = hi = np.nan
                    rows.append(
                        dict(
                            Experiment=result["name"],
                            Treatment=treatment,
                            Comparator=control,
                            Metric=metric,
                            Week=int(w),
                            Contrast=contrast,
                            Original_statistic=float(original),
                            Permutation_statistic=float(current[0]),
                            Numeric_difference=float(current[0] - original),
                            P_raw=p,
                            P_MC_lower=lo,
                            P_MC_upper=hi,
                            Permutations=b,
                            Extreme_count=extreme,
                            Nodes_per_kernel=2 ** powers[w] if has_shannon else 0,
                            N_Treatment=int(
                                summary.loc[
                                    (treatment, w), "N_NMD" if metric == "NMD" else "N_entropy"
                                ]
                            )
                            if (treatment, w) in summary.index
                            else 0,
                            N_Comparator=int(
                                summary.loc[
                                    (control, w), "N_NMD" if metric == "NMD" else "N_entropy"
                                ]
                            )
                            if (control, w) in summary.index
                            else 0,
                            Status="探索性：试验内整名患者置换"
                            if np.isfinite(p) and result["mode"] == "clinical"
                            else "探索性：整只动物置换"
                            if np.isfinite(p)
                            else "NA：" + reason,
                        )
                    )
        return correction(rows, active if result.get("selected_method") is not None else None)

    with threadpool_limits(limits=4):
        for w in weeks:
            compute_week(w)
        rows = rows_from_values()
        first = {(r["Metric"], r["Week"], r["Contrast"]): r["P_Holm"] for r in rows}
        refine = {
            r["Week"]
            for r in rows
            if r["Metric"] == "Entropy_Shannon" and 0.025 <= r["P_Holm"] <= 0.1
        }
        for w in sorted(refine):
            compute_week(w, powers[w] + 2)
        if refine:
            rows = rows_from_values()
    for row in rows:
        p = row["P_Holm"]
        key = (row["Metric"], row["Week"], row["Contrast"])
        changed = np.isfinite(p) and stars(first[key]) != stars(p)
        mc_border = np.isfinite(p) and any(
            row["P_Holm_MC_lower"] <= v <= row["P_Holm_MC_upper"] for v in (0.001, 0.01, 0.05)
        )
        numeric = row["Metric"] == "Entropy_Shannon" and abs(row["Numeric_difference"]) > 0.01
        row["Annotation"] = (
            f"p={p:.3g}†" if np.isfinite(p) and (changed or mc_border or numeric) else stars(p)
        )
        row["Caution"] = "; ".join(
            t
            for flag, t in [
                (changed, "积分加密后星号级别变化"),
                (mc_border, "Monte Carlo精度跨越星号阈值"),
                (numeric, "检验积分与原组统计量差异超过0.01 nats"),
            ]
            if flag
        )
        row["Initial_P_Holm"] = first[key]
    result["permutation_numeric_checks"] = checks
    result["inference_settings"]["status"] = "完成" if b else "用户关闭"
    result["inference_settings"]["refined_weeks"] = sorted(map(int, refine))
    return rows
