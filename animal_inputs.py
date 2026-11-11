from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.covariance import LedoitWolf
from clinical_inputs import InputSelectionError


def read_animal(path, reference=None):
    path = Path(path)
    raw = pd.read_excel(path)
    if (
        len(raw.columns) < 3
        or str(raw.columns[0]).strip().lower() != "group"
        or str(raw.columns[1]).strip().lower() not in ("sample id", "sample_id")
    ):
        raise InputSelectionError("动物工作簿前两列应为 Group 和 Sample ID，其后为数值指标。")
    raw = raw.dropna(how="all").copy()
    if raw.iloc[:, :2].isna().any().any():
        raise InputSelectionError("动物分组或编号为空，请核查原表。")
    groups = raw.iloc[:, 0].astype(str).str.strip()
    ids = raw.iloc[:, 1]
    if groups.eq("").any() or ids.astype(str).str.strip().eq("").any():
        raise InputSelectionError("动物分组或编号为空。")
    if pd.DataFrame({"Group": groups, "ID": ids.astype(str)}).duplicated().any():
        raise InputSelectionError("同一组内有重复动物编号，不能自动合并。")
    markers = raw.iloc[:, 2:].drop(
        columns=[
            c
            for c in raw.columns[2:]
            if str(c).strip().lower()
            in ("usubjid", "patient id", "patient_id", "trial", "study", "week")
        ],
        errors="ignore",
    )
    if len({str(c).strip().lower() for c in markers.columns}) != len(markers.columns):
        raise InputSelectionError("存在重复指标列名。")
    numeric = markers.apply(pd.to_numeric, errors="coerce")
    invalid = markers.notna() & numeric.isna()
    if invalid.any().any() or np.isinf(numeric.to_numpy(float)).any():
        raise InputSelectionError("指标中存在非数值或无穷值，请核查后再计算。")
    features = sorted(
        c
        for c in numeric
        if str(c).strip().lower() != "x"
        and numeric[c].notna().mean() >= 0.7
        and (numeric[c].nunique(dropna=True) > 1)
    )
    if len(features) < 2:
        raise InputSelectionError(
            "满足可用性条件的动物指标不足2个（非缺失比例至少70%，且有变异）。"
        )
    excluded_features = [
        {
            "Feature": c,
            "Observed_fraction": float(numeric[c].notna().mean()),
            "Reason": "辅助列x / 非缺失比例低于70% / 无变异",
        }
        for c in numeric
        if c not in features
    ]
    x = numeric[features]
    center = x.mean()
    scale = x.std(ddof=1)
    if not np.isfinite(scale).all() or (scale <= 0).any():
        raise InputSelectionError("动物指标的共同标准差无效。")
    normal_names = [
        g
        for g in groups.unique()
        if "".join(g.lower().split())
        in ("normalcontrol", "normal", "control", "正常组", "正常对照组")
    ]
    normal = (
        reference if reference is not None else normal_names[0] if len(normal_names) == 1 else None
    )
    order = ([normal] if normal in set(groups) else []) + [
        g for g in groups.unique() if g != normal
    ]
    complete = np.isfinite(x).all(axis=1)
    return dict(
        path=path,
        raw=raw,
        x=x,
        features=features,
        groups=groups,
        ids=ids,
        center=center,
        scale=scale,
        z=(x - center) / scale,
        complete=complete,
        group_order=order,
        normal=normal,
        excluded_features=excluded_features,
    )


def nmd_animal(data):
    x = data["x"].to_numpy(float)
    normal = (
        data["reference_x"]
        if "reference_x" in data
        else x[data["groups"].eq(data["normal"]).to_numpy()]
    )
    valid_count = np.isfinite(normal).sum(axis=0)
    if (valid_count < 3).any():
        raise InputSelectionError("正常组至少一个指标有效值不足3个，NMD无法计算。")
    scale = np.nanstd(normal, axis=0, ddof=1)
    if not np.isfinite(scale).all() or (scale <= 0).any():
        raise InputSelectionError("正常组存在无效或零标准差指标，NMD未计算。")
    normal = normal[np.isfinite(normal).any(axis=1)]
    median = np.nanmedian(normal, axis=0)
    norm_filled = np.where(np.isfinite(normal), normal, median)
    mu = norm_filled.mean(axis=0)
    fit = LedoitWolf().fit((norm_filled - mu) / scale)
    if np.linalg.eigvalsh(fit.covariance_).min() <= 0:
        raise InputSelectionError("正常参照协方差不可用。")
    precision = fit.precision_ / np.outer(scale, scale)
    observed = np.isfinite(x).sum(axis=1)
    delta = np.where(np.isfinite(x), x, mu) - mu
    distances = np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", delta, precision, delta), 0))
    distances[observed == 0] = np.nan
    meta = dict(
        Normal_group=data["normal"],
        Normal_n=len(norm_filled),
        features=data["features"],
        mean=mu.tolist(),
        observed_SD_ddof1=scale.tolist(),
        covariance_standardized=fit.covariance_.tolist(),
        shrinkage=float(fit.shrinkage_),
        estimator="LedoitWolf in Normal mean/SD coordinates",
        missing_Normal="exclude entirely missing rows; fill partial values by Normal feature median",
        missing_study="fill partial values by fitted Normal mean; entirely missing row NMD=NA",
        log_transform=False,
        clipping=False,
    )
    status = np.where(
        observed == 0,
        "NA：全部指标缺失",
        np.where(observed < len(data["features"]), "部分缺失按Normal均值填补（原NMD规则）", "完整"),
    )
    return (distances, status, meta)
