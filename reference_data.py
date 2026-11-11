import numpy as np
from sklearn.covariance import LedoitWolf
from clinical_inputs import InputSelectionError


def read_reference(paths, features):
    from clinical_patient_ids import flat_read, numeric

    rows = []
    for path in paths:
        h, rr, _ = flat_read(path)
        cols = {str(x).strip(): j for j, x in enumerate(h) if x is not None}
        missing = [f for f in features if f not in cols]
        if missing:
            raise InputSelectionError(str(path) + "：参照缺少同名同单位指标：" + "、".join(missing))
        for num, r in rr:
            for f in features:
                val = r[cols[f]]
                if (
                    val is not None
                    and str(val).strip() not in ("", "NA", "N/A", "NaN", "nan")
                    and (not np.isfinite(numeric(val)))
                ):
                    raise InputSelectionError(f"{path}第{num}行：参照指标{f}需填写有效数值。")
            rows.append([numeric(r[cols[f]]) for f in features])
    return np.array(rows, dtype=float).reshape((-1, len(features)))


def clinical_nmd(clean, indices, reference, label):
    p = len(indices)
    N = np.full(clean.shape[:2], np.nan)
    A = np.asarray(reference, float).reshape((-1, p))
    good = np.isfinite(A).all(1) & (A >= 0).all(1)
    A = A[good]
    meta = dict(
        n=len(A),
        excluded=int((~good).sum()),
        mean=np.full(p, np.nan),
        sd=np.full(p, np.nan),
        reference=label,
    )
    reason = ""
    if len(A) < 3:
        reason = "完整正常参照不足3行或未指定参照"
    elif np.any(A.std(0, ddof=1) <= 0):
        reason = "正常参照存在零标准差指标"
    if reason:
        meta["status"] = "NA：" + reason
        return (N, meta)
    mu = A.mean(0)
    sd = A.std(0, ddof=1)
    fit = LedoitWolf().fit((A - mu) / sd)
    if np.linalg.eigvalsh(fit.covariance_).min() <= 0:
        meta["status"] = "NA：正常协方差不可用"
        return (N, meta)
    precision = fit.precision_ / np.outer(sd, sd)
    ok = np.isfinite(clean[:, :, indices]).all(2)
    delta = clean[:, :, indices][ok] - mu
    N[ok] = np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", delta, precision, delta), 0))
    meta.update(
        mean=mu,
        sd=sd,
        covariance=fit.covariance_,
        precision_raw=precision,
        shrinkage=fit.shrinkage_,
        status="可计算",
    )
    return (N, meta)
