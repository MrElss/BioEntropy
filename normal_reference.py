from __future__ import annotations
import io
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple, Any
import numpy as np
import pandas as pd
from sklearn.covariance import LedoitWolf
from constants import NMD_PRIMARY_METHOD, SOFTWARE_VERSION
import re as _re


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


GROUP_SYNONYMS = {
    "达格列净组": "达格列净",
    "达格列净 ": "达格列净",
    "安慰剂组": "安慰剂",
    "placebo group": "Placebo",
    "placebo": "Placebo",
    "control group": "Control",
    "vehicle": "Vehicle",
}


def normalize_group_name(x: Any) -> str:
    x = "" if pd.isna(x) else str(x).strip()
    return GROUP_SYNONYMS.get(x, x)


def _read_excel_any(obj: Any) -> Tuple[str, pd.DataFrame]:
    if hasattr(obj, "read") and hasattr(obj, "name"):
        content = obj.read()
        name = obj.name
    elif isinstance(obj, tuple) and len(obj) == 2:
        name, content = obj
    else:
        path = Path(obj)
        name = path.name
        content = path.read_bytes()
    df = pd.read_excel(io.BytesIO(content))
    return (name, df)


def get_weeks(raw: Dict[Tuple[str, int], pd.DataFrame], experiment: str) -> List[int]:
    return sorted([week for exp, week in raw if exp == experiment])


def common_features(
    raw: Dict[Tuple[str, int], pd.DataFrame], experiment: str, missing_thresh: float = 0.3
) -> List[str]:
    feature_sets = []
    baseline_week = anchor_week(raw, experiment)
    for week in get_weeks(raw, experiment):
        df = raw[experiment, week]
        feats = []
        for col in df.columns[2:]:
            if str(col).strip().lower() == "x":
                continue
            s = pd.to_numeric(df[col], errors="coerce")
            if s.notna().mean() >= 1 - missing_thresh and (
                week != baseline_week or s.nunique(dropna=True) > 1
            ):
                feats.append(col)
        feature_sets.append(set(feats))
    common = sorted(set.intersection(*feature_sets)) if feature_sets else []
    return common


def anchor_week(raw: Dict[Tuple[str, int], pd.DataFrame], experiment: str) -> int:
    weeks = get_weeks(raw, experiment)
    return 0 if 0 in weeks else min(weeks)


EXTERNAL_NORMAL_FILE_ALIASES = {
    "normal",
    "norma",
    "normalrange",
    "normal_range",
    "healthy",
    "healthyrange",
    "healthy_range",
    "reference",
    "referencerange",
    "reference_range",
    "ref",
    "refrange",
    "normalref",
    "normalreference",
}


def _canonical_text(x: Any) -> str:
    text = "" if pd.isna(x) else str(x).strip().lower()
    return _re.sub("[^0-9a-zA-Z\\u4e00-\\u9fff]+", "", text)


def _parse_external_normal_range_df(df: pd.DataFrame, file_name: str) -> pd.DataFrame:
    df = df.copy()
    if df.empty or df.shape[1] < 2:
        raise ValueError(f"{file_name} 为空，或列数不足。")
    first_col = df.iloc[:, 0].astype(str).str.strip()
    first_col_norm = first_col.map(_canonical_text)
    has_min_row = first_col_norm.str.contains("normalmin|min|下限").any()
    has_max_row = first_col_norm.str.contains("normalmax|max|上限").any()
    if has_min_row or has_max_row:
        min_idx = first_col_norm[first_col_norm.str.contains("normalmin|min|下限")].index
        max_idx = first_col_norm[first_col_norm.str.contains("normalmax|max|上限")].index
        min_idx = int(min_idx[0]) if len(min_idx) > 0 else None
        max_idx = int(max_idx[0]) if len(max_idx) > 0 else None
        rows = []
        for col in df.columns[2:]:
            lo = pd.to_numeric(
                pd.Series([df.loc[min_idx, col]]) if min_idx is not None else pd.Series([np.nan]),
                errors="coerce",
            ).iloc[0]
            hi = pd.to_numeric(
                pd.Series([df.loc[max_idx, col]]) if max_idx is not None else pd.Series([np.nan]),
                errors="coerce",
            ).iloc[0]
            if pd.isna(lo) and pd.isna(hi):
                continue
            rows.append(
                {
                    "Feature": str(col).strip(),
                    "Normal_min": lo,
                    "Normal_max": hi,
                    "Source_file": file_name,
                    "Format": "wide_min_max_rows",
                }
            )
        if not rows:
            raise ValueError(f"{file_name} 没有识别到任何有效的 Normal_min / Normal_max 指标。")
        return pd.DataFrame(rows)
    col_map = {_canonical_text(c): str(c) for c in df.columns}
    feature_col = None
    min_col = None
    max_col = None
    for key, col in col_map.items():
        if key in {"feature", "biomarker", "indicator", "指标", "特征"}:
            feature_col = col
        elif key in {"normalmin", "min", "下限", "lower", "low"}:
            min_col = col
        elif key in {"normalmax", "max", "上限", "upper", "high"}:
            max_col = col
    if feature_col is not None and (min_col is not None or max_col is not None):
        out = pd.DataFrame(
            {
                "Feature": df[feature_col].astype(str).str.strip(),
                "Normal_min": pd.to_numeric(df[min_col], errors="coerce")
                if min_col is not None
                else np.nan,
                "Normal_max": pd.to_numeric(df[max_col], errors="coerce")
                if max_col is not None
                else np.nan,
            }
        )
        out["Source_file"] = file_name
        out["Format"] = "long_feature_min_max"
        out = out.loc[~(out["Normal_min"].isna() & out["Normal_max"].isna())].reset_index(drop=True)
        if out.empty:
            raise ValueError(f"{file_name} 没有识别到任何有效的 Normal_min / Normal_max 指标。")
        return out
    raise ValueError(
        f"{file_name} 的格式无法识别。请使用两行宽表格式（Normal_min / Normal_max），或长表格式（Feature / Normal_min / Normal_max）。"
    )


TRUE_NORMAL_FILE_ALIASES = {
    "normal",
    "normalsample",
    "normalsamples",
    "normalgroup",
    "normal_sample",
    "healthy",
    "healthysample",
    "healthysamples",
    "healthygroup",
    "healthy_sample",
    "referencenormal",
    "referencesample",
    "truernormal",
    "truenormal",
    "healthreference",
}


@dataclass
class TrueNormalReference:
    experiment: str
    features: List[str]
    mean: np.ndarray
    precision: np.ndarray
    covariance: np.ndarray
    source_feature_map: Dict[str, str]
    unit_factors: Dict[str, float]
    use_log_map: Dict[str, bool]
    per_feature_mean: Dict[str, float]
    per_feature_sd: Dict[str, float]
    n_normal: int
    source_files: List[str]
    clip_bounds: Dict[str, Tuple[float, float]]


def _strip_copy_suffix(text: str) -> str:
    return _re.sub("\\(\\d+\\)$", "", str(text).strip())


def _simplify_text(text: Any) -> str:
    x = "" if pd.isna(text) else str(text)
    x = x.replace("（", "(").replace("）", ")").replace("，", ",").replace("、", "")
    x = x.replace(" ", "").replace("-", "").replace("_", "")
    x = x.lower()
    x = x.replace("γ", "gamma").replace("α", "alpha").replace("β", "beta")
    x = _re.sub("[^0-9a-zA-Z\\u4e00-\\u9fff]+", "", x)
    return x


def feature_semantic_code(feature_name: Any) -> Optional[str]:
    s = _simplify_text(feature_name)
    mapping_rules = [
        (["非高密度脂蛋白胆固醇", "nonhdl"], "non_hdl"),
        (["极低密度脂蛋白胆固醇", "vldl"], "vldl"),
        (["高密度脂蛋白胆固醇", "hdlc", "hdl"], "hdl"),
        (["低密度脂蛋白胆固醇", "ldlc", "ldl"], "ldl"),
        (["空腹血浆血糖", "空腹血糖", "fastingglucose", "fpg"], "glu_fasting"),
        (
            [
                "餐后2h血浆血糖",
                "餐后2h血糖",
                "餐前餐后血糖",
                "postprandialglucose",
                "ppg",
                "2hglucose",
            ],
            "glu_postprandial",
        ),
        (["糖化血红蛋白", "糖化", "hba1c"], "hba1c"),
        (["血糖", "bloodglucose", "plasmaglucose", "glucose"], "glu_blood"),
        (["甘油三酯", "triglyceride", "tg"], "tg"),
        (["总胆固醇", "tc", "cholesteroltotal"], "tc"),
        (["丙氨酸氨基转移酶", "alanineaminotransferase", "alt"], "alt"),
        (["天门冬氨酸氨基转移酶", "aspartateaminotransferase", "ast"], "ast"),
        (["gamma谷氨酰转肽酶", "谷氨酰转肽酶", "ggt", "gammaglutamyltransferase"], "ggt"),
        (["超敏c反应蛋白", "c反应蛋白", "hscrp", "crp"], "crp"),
        (["肌酐", "creatinine"], "creatinine"),
        (["homas", "homair", "homa"], "homa"),
        (["auc"], "auc_like"),
    ]
    for keywords, code in mapping_rules:
        if any(k in s for k in keywords):
            return code
    return None


def parse_input_file_kind(file_name: str) -> Tuple[str, str, Optional[int]]:
    stem = Path(file_name).stem.strip()
    if "-" not in stem:
        canon_stem = _canonical_text(stem)
        if canon_stem in {_canonical_text(x) for x in TRUE_NORMAL_FILE_ALIASES} or canon_stem in {
            _canonical_text(x) for x in EXTERNAL_NORMAL_FILE_ALIASES
        }:
            return ("normal_like", None, None)
        return ("week", stem, 0)
    experiment, suffix = stem.rsplit("-", 1)
    experiment = experiment.strip()
    suffix = _strip_copy_suffix(suffix.strip())
    m = _re.fullmatch("(\\d+)\\s*[Ww]", suffix)
    if m:
        return ("week", experiment, int(m.group(1)))
    canon = _canonical_text(suffix)
    if canon in {_canonical_text(x) for x in EXTERNAL_NORMAL_FILE_ALIASES}:
        return ("normal_like", experiment, None)
    if canon in {_canonical_text(x) for x in TRUE_NORMAL_FILE_ALIASES}:
        return ("normal_like", experiment, None)
    return ("week", stem, 0)


def _looks_like_external_normal_range_df(df: pd.DataFrame) -> bool:
    if df is None or df.empty:
        return False
    first_col = df.iloc[:, 0].astype(str).str.strip()
    first_col_norm = first_col.map(_canonical_text)
    if (
        first_col_norm.str.contains("normalmin|min|下限").any()
        or first_col_norm.str.contains("normalmax|max|上限").any()
    ):
        return True
    col_norm = {_canonical_text(c): str(c) for c in df.columns}
    has_feature_col = any(
        k in {"feature", "biomarker", "indicator", "指标", "特征"} for k in col_norm
    )
    has_min_col = any(k in {"normalmin", "min", "下限", "lower", "low"} for k in col_norm)
    has_max_col = any(k in {"normalmax", "max", "上限", "upper", "high"} for k in col_norm)
    return has_feature_col and (has_min_col or has_max_col)


def _parse_true_normal_sample_df(df: pd.DataFrame, file_name: str) -> pd.DataFrame:
    if df.shape[1] < 3:
        raise ValueError(f"{file_name} 至少需要 3 列：Group、Sample ID、指标列。")
    out = df.copy()
    first_header = _canonical_text(out.columns[0])
    id_tokens = (
        "id",
        "编号",
        "序号",
        "编码",
        "健康人",
        "subject",
        "subjectid",
        "sampleid",
        "sample",
    )
    has_explicit_group = _canonical_text(out.columns[1]) in {
        "sampleid",
        "sample",
        "id",
        "编号",
        "序号",
    }
    if any(tok in first_header for tok in id_tokens) and (not has_explicit_group):
        out.insert(0, "Group", "Normal")
        out = out.rename(columns={out.columns[1]: "Sample ID"})
    else:
        out = out.rename(columns={out.columns[0]: "Group", out.columns[1]: "Sample ID"})
    out["Group"] = out["Group"].ffill().map(normalize_group_name)
    if out["Group"].isna().all():
        out["Group"] = "Normal"
    out["Group"] = out["Group"].fillna("Normal")
    out["Sample ID"] = out["Sample ID"].map(normalize_sample_id)
    numeric_feature_count = 0
    for col in out.columns[2:]:
        s = pd.to_numeric(out[col], errors="coerce")
        if s.notna().sum() >= max(3, int(len(out) * 0.5)):
            out[col] = s
            numeric_feature_count += 1
    if len(out) < 3 or numeric_feature_count < 2:
        raise ValueError(f"{file_name} 看起来不像真实 Normal 样本表。请检查格式。")
    out["__source_file__"] = file_name
    return out


def normalize_sample_id(value):
    if pd.isna(value):
        return np.nan
    if isinstance(value, str):
        return value.strip() or np.nan
    return value


def sample_id_key(value):
    value = normalize_sample_id(value)
    if pd.isna(value):
        return np.nan
    if isinstance(value, (int, float, np.integer, np.floating)):
        return int(value) if float(value).is_integer() else value
    return str(value)


def _validate_week_dataframe(name: str, df: pd.DataFrame) -> None:
    ids = df["Sample ID"].map(sample_id_key).dropna()
    dups = ids[ids.duplicated()].unique()
    if len(dups) > 0:
        shown = ", ".join(str(d) for d in dups[:8])
        raise ValueError(
            f"{name}: 检测到重复的 Sample ID（{shown}{(' …' if len(dups) > 8 else '')}）。同一文件内每个样本的 ID 必须唯一。"
        )
    feature_cols = list(df.columns[2:])
    if not feature_cols:
        raise ValueError(f"{name}: 没有任何指标列（至少需要 Group、Sample ID 和 1 个数值指标列）。")
    bad_cols = []
    for col in feature_cols:
        s = df[col]
        non_empty = s[s.notna() & (s.astype(str).str.strip() != "")]
        if len(non_empty) == 0:
            continue
        if pd.to_numeric(non_empty, errors="coerce").notna().sum() == 0:
            bad_cols.append(str(col))
    if bad_cols:
        raise ValueError(
            f"{name}: 以下指标列没有任何可识别的数值，请检查是否填入了文本、单位或选错了列：{', '.join(bad_cols[:6])}{(' …' if len(bad_cols) > 6 else '')}。"
        )


def load_input_bundle(
    files: Iterable[Any],
) -> Tuple[
    Dict[Tuple[str, int], pd.DataFrame],
    Dict[str, pd.DataFrame],
    Dict[str, pd.DataFrame],
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
]:
    raw: Dict[Tuple[str, int], pd.DataFrame] = {}
    range_refs: Dict[str, pd.DataFrame] = {}
    normal_refs: Dict[str, pd.DataFrame] = {}
    global_normal_frames: List[pd.DataFrame] = []
    global_range_frames: List[pd.DataFrame] = []
    data_rows: List[Dict[str, Any]] = []
    range_rows: List[Dict[str, Any]] = []
    normal_rows: List[Dict[str, Any]] = []
    mixed_rows: List[Dict[str, Any]] = []
    for file_obj in files:
        name, df = _read_excel_any(file_obj)
        kind, exp, week = parse_input_file_kind(name)
        if kind == "week":
            if df.shape[1] < 3:
                raise ValueError(f"{name} 至少需要 3 列：Group、Sample ID、指标列。")
            df = df.copy()
            df = df.rename(columns={df.columns[0]: "Group", df.columns[1]: "Sample ID"})
            df["Group"] = df["Group"].ffill().map(normalize_group_name)
            if df["Group"].isna().all():
                raise ValueError(f"{name} 的 Group 列为空。")
            df["Sample ID"] = df["Sample ID"].map(normalize_sample_id)
            _validate_week_dataframe(name, df)
            groups = sorted(df["Group"].dropna().astype(str).unique().tolist())
            if (exp, week) in raw:
                raise ValueError(f"检测到重复文件：实验 {exp} 的 {week}W 出现了两次。")
            raw[exp, week] = df
            data_rows.append(
                {
                    "文件名": name,
                    "实验号": exp,
                    "周数": week,
                    "样本数": int(len(df)),
                    "分组数": len(groups),
                    "分组列表": " | ".join(groups),
                    "原始列数": int(df.shape[1]),
                    "类型": "实验数据",
                }
            )
            mixed_rows.append(
                {
                    "文件名": name,
                    "类型": "实验数据",
                    "实验号": exp,
                    "周数/说明": f"{week}W",
                    "记录数": int(len(df)),
                    "说明": f"{len(groups)} 个分组",
                }
            )
            continue
        exp_label = exp if exp is not None else "（全部实验）"
        if _looks_like_external_normal_range_df(df):
            parsed = _parse_external_normal_range_df(df, name)
            if exp is None:
                global_range_frames.append(parsed)
            elif exp in range_refs:
                raise ValueError(f"实验 {exp} 的外部 Normal 范围文件重复上传。请仅保留一个。")
            else:
                range_refs[exp] = parsed
            valid_n = int(len(parsed))
            feature_list = " | ".join(parsed["Feature"].astype(str).tolist())
            range_rows.append(
                {
                    "文件名": name,
                    "实验号": exp_label,
                    "参考类型": "外部Normal范围",
                    "有效指标数": valid_n,
                    "格式": parsed["Format"].iloc[0]
                    if "Format" in parsed.columns and len(parsed) > 0
                    else "",
                    "指标列表": feature_list,
                }
            )
            mixed_rows.append(
                {
                    "文件名": name,
                    "类型": "Normal范围",
                    "实验号": exp_label,
                    "周数/说明": "外部区间",
                    "记录数": valid_n,
                    "说明": f"{valid_n} 个有界指标",
                }
            )
        else:
            parsed = _parse_true_normal_sample_df(df, name)
            if exp is None:
                global_normal_frames.append(parsed)
            elif exp in normal_refs:
                normal_refs[exp] = pd.concat([normal_refs[exp], parsed], ignore_index=True)
            else:
                normal_refs[exp] = parsed
            valid_n = int(len(parsed))
            feat_candidates = []
            for col in parsed.columns[2:]:
                if col == "__source_file__":
                    continue
                s = pd.to_numeric(parsed[col], errors="coerce")
                if s.notna().sum() >= 3 and s.nunique(dropna=True) > 1:
                    feat_candidates.append(col)
            normal_rows.append(
                {
                    "文件名": name,
                    "实验号": exp_label,
                    "参考类型": "真实Normal样本",
                    "样本数": valid_n,
                    "可用指标数": len(feat_candidates),
                    "分组": " | ".join(
                        sorted(parsed["Group"].dropna().astype(str).unique().tolist())
                    ),
                }
            )
            mixed_rows.append(
                {
                    "文件名": name,
                    "类型": "Normal样本",
                    "实验号": exp_label,
                    "周数/说明": "独立样本",
                    "记录数": valid_n,
                    "说明": f"{len(feat_candidates)} 个可用指标",
                }
            )
    if not raw:
        raise ValueError("没有读取到任何实验时间点 Excel 文件。")
    experiments_present = sorted({exp for exp, _week in raw.keys()})
    if global_normal_frames:
        shared_normal = pd.concat(global_normal_frames, ignore_index=True)
        for exp in experiments_present:
            if exp not in normal_refs:
                normal_refs[exp] = shared_normal.copy()
    if global_range_frames:
        shared_range = pd.concat(global_range_frames, ignore_index=True)
        for exp in experiments_present:
            if exp not in range_refs:
                range_refs[exp] = shared_range.copy()
    data_overview = (
        pd.DataFrame(data_rows).sort_values(["实验号", "周数", "文件名"]).reset_index(drop=True)
    )
    range_overview = pd.DataFrame(range_rows)
    if not range_overview.empty:
        range_overview = range_overview.sort_values(["实验号", "文件名"]).reset_index(drop=True)
    normal_overview = pd.DataFrame(normal_rows)
    if not normal_overview.empty:
        normal_overview = normal_overview.sort_values(["实验号", "文件名"]).reset_index(drop=True)
    mixed_overview = (
        pd.DataFrame(mixed_rows).sort_values(["实验号", "类型", "文件名"]).reset_index(drop=True)
    )
    return (
        raw,
        range_refs,
        normal_refs,
        mixed_overview,
        data_overview,
        range_overview,
        normal_overview,
    )


def _semantic_match_true_normal_features(
    experiment_features: List[str], normal_df: pd.DataFrame
) -> Dict[str, str]:
    normal_cols = [c for c in normal_df.columns[2:] if c != "__source_file__"]
    normal_by_canon: Dict[str, List[str]] = {}
    for col in normal_cols:
        normal_by_canon.setdefault(_simplify_text(col), []).append(col)
    mapping: Dict[str, str] = {}
    for feat in experiment_features:
        exact_candidates = normal_by_canon.get(_simplify_text(feat), [])
        if exact_candidates:
            mapping[feat] = sorted(exact_candidates, key=lambda c: len(str(c)))[0]
    normal_by_code: Dict[str, List[str]] = {}
    for col in normal_cols:
        code = feature_semantic_code(col)
        if code is not None:
            normal_by_code.setdefault(code, []).append(col)
    for feat in experiment_features:
        if feat in mapping:
            continue
        code = feature_semantic_code(feat)
        if code is None:
            continue
        candidates = normal_by_code.get(code, [])
        if candidates:
            if len(candidates) == 1:
                mapping[feat] = candidates[0]
            else:
                best = sorted(
                    candidates,
                    key=lambda c: (
                        0 if _simplify_text(c) == _simplify_text(feat) else 1,
                        len(str(c)),
                    ),
                )[0]
                mapping[feat] = best
    return mapping


def fit_true_normal_reference(
    raw: Dict[Tuple[str, int], pd.DataFrame], experiment: str, normal_df: pd.DataFrame
) -> Tuple[TrueNormalReference, pd.DataFrame]:
    exp_features = common_features(raw, experiment)
    pooled = pd.concat([raw[experiment, w] for w in get_weeks(raw, experiment)], ignore_index=True)
    mapping = _semantic_match_true_normal_features(exp_features, normal_df)
    if len(mapping) < 2:
        raise ValueError(
            f"实验 {experiment} 的独立 Normal 样本文件与实验共同指标匹配后不足 2 个，无法计算基于 Normal 样本的疾病偏离度。"
        )
    rows = []
    normal_t = pd.DataFrame(index=normal_df.index)
    mean_map: Dict[str, float] = {}
    sd_map: Dict[str, float] = {}
    unit_factors: Dict[str, float] = {}
    use_log_map: Dict[str, bool] = {}
    clip_bounds: Dict[str, Tuple[float, float]] = {}
    source_files = (
        sorted(normal_df["__source_file__"].dropna().astype(str).unique().tolist())
        if "__source_file__" in normal_df.columns
        else []
    )
    for feat, src_feat in mapping.items():
        study_s = pd.to_numeric(pooled[feat], errors="coerce").astype(float)
        norm_s_raw = (
            pd.to_numeric(normal_df[src_feat], errors="coerce")
            .astype(float)
            .replace([np.inf, -np.inf], np.nan)
        )
        norm_s, unit_factor = (norm_s_raw.copy(), 1.0)
        unit_note = "No inferred conversion; study and Normal inputs must use the same units."
        combo = pd.concat([study_s, norm_s], ignore_index=True)
        combo = combo.dropna().astype(float)
        use_log = False
        if use_log:
            t_norm = np.log1p(norm_s)
        else:
            t_norm = norm_s.copy()
        valid = t_norm.dropna().astype(float)
        if valid.shape[0] < 3:
            continue
        lo = hi = float("nan")
        clip_bounds[feat] = (lo, hi)
        t_norm_clip = t_norm.copy()
        normal_t[feat] = t_norm_clip
        mean_map[feat] = float(pd.to_numeric(t_norm_clip, errors="coerce").dropna().mean())
        sd_val = float(pd.to_numeric(t_norm_clip, errors="coerce").dropna().std(ddof=1))
        if not np.isfinite(sd_val) or sd_val <= 0:
            raise ValueError(
                f"Normal reference feature {feat} has no positive SD; NMD is undefined."
            )
        sd_map[feat] = sd_val
        unit_factors[feat] = unit_factor
        use_log_map[feat] = use_log
        rows.append(
            {
                "Experiment": experiment,
                "Feature": feat,
                "Source_feature": src_feat,
                "Semantic_code": feature_semantic_code(feat),
                "Normal_n_nonmissing": int(norm_s.dropna().shape[0]),
                "Use_log1p": use_log,
                "Unit_factor_applied": unit_factor,
                "Unit_note": unit_note,
                "Normal_raw_median_after_unit": float(norm_s.dropna().median())
                if norm_s.dropna().shape[0] > 0
                else np.nan,
                "Study_raw_median": float(study_s.dropna().median())
                if study_s.dropna().shape[0] > 0
                else np.nan,
                "Normal_reference_clipped": False,
                "Clip_lower_transformed": lo,
                "Clip_upper_transformed": hi,
                "Source_files": " | ".join(source_files),
            }
        )
    if normal_t.shape[1] < 2:
        raise ValueError(
            f"实验 {experiment} 的独立 Normal 样本文件经过匹配与预处理后，只剩 {normal_t.shape[1]} 个可用指标，无法稳定计算 Mahalanobis distance。"
        )
    features = list(normal_t.columns)
    normal_observed = normal_t[features].to_numpy(dtype=float)
    normal_evaluable = np.isfinite(normal_observed).any(axis=1)
    normal_retained = (
        normal_observed if normal_evaluable.all() else normal_observed[normal_evaluable]
    )
    X_norm = fill_nan_with_col_median(normal_retained)
    mu = np.asarray(np.mean(X_norm, axis=0), dtype=float)
    scale = np.asarray([sd_map[f] for f in features], dtype=float)
    Z_norm = (X_norm - mu) / scale
    lw = LedoitWolf().fit(Z_norm)
    if np.any(np.linalg.eigvalsh(lw.covariance_) <= 0):
        raise ValueError("Normal covariance is not positive definite; NMD is undefined.")
    precision = np.asarray(lw.precision_) / np.outer(scale, scale)
    covariance = np.asarray(lw.covariance_) * np.outer(scale, scale)
    for row in rows:
        row["Normal_mean_after_imputation"] = float(mu[features.index(row["Feature"])])
        row["Normal_SD_ddof1_observed"] = sd_map[row["Feature"]]
        row["Covariance_fit_scale"] = "Normal mean/SD; dimensionless"
        row["Covariance_estimator"] = "LedoitWolf"
        row["Covariance_shrinkage"] = float(lw.shrinkage_)
        row["Normal_n_input_rows"] = len(normal_observed)
        row["Normal_n_evaluable_rows"] = len(X_norm)
        row["Normal_n_excluded_no_observed_features"] = int((~normal_evaluable).sum())
        row["Missing_Normal"] = "exclude all-missing rows; partial feature median"
        row["Missing_study"] = "Normal fitted mean; all-missing sample NA"
        row["Transformation"] = "identity"
    ref = TrueNormalReference(
        experiment=experiment,
        features=features,
        mean=mu,
        precision=precision,
        covariance=covariance,
        source_feature_map={k: mapping[k] for k in features},
        unit_factors={k: unit_factors[k] for k in features},
        use_log_map={k: use_log_map[k] for k in features},
        per_feature_mean={k: mean_map[k] for k in features},
        per_feature_sd={k: sd_map[k] for k in features},
        n_normal=int(len(X_norm)),
        source_files=source_files,
        clip_bounds={k: clip_bounds[k] for k in features},
    )
    return (ref, pd.DataFrame(rows).sort_values(["Feature"]).reset_index(drop=True))


def apply_true_normal_reference(df: pd.DataFrame, ref: TrueNormalReference) -> pd.DataFrame:
    out = pd.DataFrame(index=df.index)
    feature_count = pd.Series(0, index=df.index, dtype=float)
    cols = []
    for feat in ref.features:
        s = pd.to_numeric(df[feat], errors="coerce").astype(float)
        if ref.use_log_map.get(feat, False):
            s = np.log1p(s)
        cols.append(s)
        feature_count += np.isfinite(s).astype(float)
    X = np.column_stack([c.to_numpy(dtype=float) for c in cols]).astype(float)
    for j in range(X.shape[1]):
        mask = ~np.isfinite(X[:, j])
        if mask.any():
            X[mask, j] = ref.mean[j]
    diff = X - ref.mean.reshape(1, -1)
    d2_cov = np.einsum("ij,jk,ik->i", diff, ref.precision, diff)
    d2_cov = np.clip(d2_cov, 0, None)
    md_cov = np.sqrt(d2_cov)
    sd = np.asarray(
        [max(float(ref.per_feature_sd.get(feat, 1.0)), 1e-09) for feat in ref.features], dtype=float
    )
    z = diff / sd.reshape(1, -1)
    z = np.where(np.isfinite(z), z, 0.0)
    z_capped = z
    d2 = np.sum(z_capped * z_capped, axis=1)
    d2 = np.clip(d2, 0, None)
    md = np.sqrt(d2)
    pval_diagonal = np.full(len(df), np.nan)
    pval_covariance = np.full(len(df), np.nan)
    out["NMD_cov_p_smallsample_F"] = np.full(len(df), np.nan)
    out["NMD"] = md_cov
    out["NMD_sq"] = d2_cov
    out["NMD_p_chi2"] = pval_covariance
    out["NMD_p_smallsample_F"] = out["NMD_cov_p_smallsample_F"]
    out["NMD_covariance_raw"] = md_cov
    out["NMD_sq_covariance_raw"] = d2_cov
    out["NMD_method"] = NMD_PRIMARY_METHOD
    out["NMD_z_cap"] = np.nan
    out["NMD_diagonal_capped"] = md
    out["NMD_sq_diagonal_capped"] = d2
    out["NMD_diagonal_p_chi2"] = pval_diagonal
    out["NMD_diagonal_z_cap"] = np.nan
    out["NMD_feature_count"] = feature_count.to_numpy(dtype=float)
    missing_all = feature_count.eq(0).to_numpy()
    for name in [
        "NMD",
        "NMD_sq",
        "NMD_covariance_raw",
        "NMD_sq_covariance_raw",
        "NMD_diagonal_capped",
        "NMD_sq_diagonal_capped",
    ]:
        out.loc[missing_all, name] = np.nan
    out["NMD_status"] = np.where(
        missing_all,
        "NA_no_observed_features",
        np.where(
            feature_count.to_numpy() < len(ref.features),
            "partial_features_mean_imputed",
            "observed_all_features",
        ),
    )
    out["NMD_tail_probability_status"] = "not_reported_uncalibrated_regularized_reference"
    return out
