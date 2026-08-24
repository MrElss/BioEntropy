
from __future__ import annotations

import hashlib
import io
import tempfile
import zipfile
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple, Any

from matplotlib import font_manager
import matplotlib.pyplot as plt
import numpy as np
import openpyxl
import pandas as pd
from openpyxl.styles import Font, PatternFill
from openpyxl.utils import get_column_letter
from scipy.stats import chi2, f as f_dist, linregress, mannwhitneyu, skew, spearmanr
from sklearn.covariance import LedoitWolf

warnings.filterwarnings("ignore")

from bioentropy_constants import (
    RNG_SEED,
    WINSOR_Q_DEFAULT,
    ENTROPY_STRONG_OUTLIER_METHODS,
    NMD_Z_CAP,
)
from bioentropy_entropy import *  # noqa: F401,F403  (entropy primitives + infer_direction)


def _boot_seed(base: int, week: int, group: object) -> int:
    """Deterministic, process-stable bootstrap seed for a (metric-family, week, group).

    Derives the seed from a stable hash of the group *identity* rather than
    ``len(group_name)``. The old ``RNG_SEED + week + len(group)`` scheme gave two
    groups with equal-length names the same seed (correlated resampling) and made
    results depend on how a group happened to be spelled. Folding in ``RNG_SEED``
    and the per-family ``base`` offset keeps the metric families independent while
    the result is identical across processes (unlike the salted built-in ``hash``).
    """
    key = f"{RNG_SEED}|{int(base)}|{int(week)}|{group}".encode("utf-8")
    return int.from_bytes(hashlib.blake2b(key, digest_size=4).digest(), "big")


GROUP_SYNONYMS = {
    "达格列净组": "达格列净",
    "达格列净 ": "达格列净",
    "安慰剂组": "安慰剂",
    "placebo group": "Placebo",
    "placebo": "Placebo",
    "control group": "Control",
    "vehicle": "Vehicle",
}

PLOT_COLORS = [
    "#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e",
    "#8c564b", "#e377c2", "#17becf", "#7f7f7f", "#bcbd22"
]


CHINESE_FONT_CANDIDATES = [
    "Microsoft YaHei",
    "Microsoft YaHei UI",
    "SimHei",
    "Noto Sans CJK SC",
    "Source Han Sans SC",
    "PingFang SC",
    "WenQuanYi Zen Hei",
    "Arial Unicode MS",
]

GROUP_TRANSLATIONS = {
    "安慰剂": "Placebo",
    "安慰剂组": "Placebo",
    "达格列净": "Dapagliflozin",
    "达格列净组": "Dapagliflozin",
    "对照": "Control",
    "对照组": "Control",
    "模型": "Model",
    "模型组": "Model",
    "正常": "Normal",
    "正常组": "Normal",
    "vehicle": "Vehicle",
    "control": "Control",
    "placebo": "Placebo",
}

def setup_matplotlib_fonts() -> tuple[str, bool]:
    """Configure a Chinese-capable matplotlib font when available."""
    available = {f.name for f in font_manager.fontManager.ttflist}
    for font_name in CHINESE_FONT_CANDIDATES:
        if font_name in available:
            plt.rcParams["font.family"] = "sans-serif"
            plt.rcParams["font.sans-serif"] = [font_name, "DejaVu Sans", "Arial"]
            plt.rcParams["axes.unicode_minus"] = False
            return font_name, True
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Arial"]
    plt.rcParams["axes.unicode_minus"] = False
    return "DejaVu Sans", False

MPL_FONT_NAME, HAS_CJK_FONT = setup_matplotlib_fonts()

def prefer_plot_text(cn: str, en: str) -> str:
    """Prefer Chinese when the runtime font supports it, otherwise use English."""
    return cn if HAS_CJK_FONT else en


@dataclass
class PreprocessParam:
    use_log: bool
    median: float
    mad: float
    lower: float
    upper: float


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
    return name, df


def get_experiments(raw: Dict[Tuple[str, int], pd.DataFrame]) -> List[str]:
    return sorted({exp for exp, _ in raw})


def get_weeks(raw: Dict[Tuple[str, int], pd.DataFrame], experiment: str) -> List[int]:
    return sorted([week for exp, week in raw if exp == experiment])


def get_groups(raw: Dict[Tuple[str, int], pd.DataFrame], experiment: str) -> List[str]:
    groups: set[str] = set()
    for (exp, _week), df in raw.items():
        if exp == experiment:
            groups |= set(df["Group"].dropna().astype(str).tolist())
    return sorted(groups)


def suggest_roles(groups: List[str]) -> Tuple[str, str]:
    """Suggest treatment/comparator using generic control keywords first.

    The defaults are only a starting point and should remain editable in the UI.
    """
    if len(groups) < 2:
        raise ValueError("至少需要 2 个分组。")

    normalized = [(g, g.lower()) for g in groups]
    disease_comparator_keywords = [
        "model", "disease", "diabetic", "db", "kk", "hfd",
        "模型", "造模", "疾病", "糖尿病", "高脂"
    ]
    comparator_keywords = [
        "placebo", "vehicle", "control", "对照", "安慰"
    ]
    healthy_reference_keywords = [
        "normal", "healthy", "sham", "blank", "正常", "健康", "空白", "假手术"
    ]
    treatment_keywords = [
        "treat", "drug", "dose", "therapy", "intervention", "rx",
        "治疗", "给药", "干预", "药物", "剂量", "低剂量", "中剂量", "高剂量"
    ]

    comparator = None
    for raw, low in normalized:
        if any(kw in low for kw in disease_comparator_keywords):
            comparator = raw
            break
    if comparator is None:
        for raw, low in normalized:
            if any(kw in low for kw in comparator_keywords) and not any(kw in low for kw in healthy_reference_keywords):
                comparator = raw
                break
    if comparator is None:
        for raw, low in normalized:
            if any(kw in low for kw in comparator_keywords + healthy_reference_keywords):
                comparator = raw
                break

    treatment = None
    for raw, low in normalized:
        if raw == comparator:
            continue
        if any(kw in low for kw in treatment_keywords):
            treatment = raw
            break

    if treatment is None:
        digit_candidates = [g for g in groups if g != comparator and any(ch.isdigit() for ch in g)]
        if digit_candidates:
            treatment = digit_candidates[0]

    if comparator is not None and treatment is None:
        treatment = next((g for g in groups if g != comparator), groups[0])

    if comparator is None:
        comparator = groups[0]
    if treatment is None:
        treatment = next((g for g in groups if g != comparator), groups[0])

    if treatment == comparator:
        comparator = next((g for g in groups if g != treatment), groups[0])

    return treatment, comparator


def common_features(
    raw: Dict[Tuple[str, int], pd.DataFrame], experiment: str, missing_thresh: float = 0.30
) -> List[str]:
    feature_sets = []
    for week in get_weeks(raw, experiment):
        df = raw[(experiment, week)]
        feats = []
        for col in df.columns[2:]:
            if str(col).strip().lower() == "x":
                continue
            s = pd.to_numeric(df[col], errors="coerce")
            if s.notna().mean() >= (1 - missing_thresh) and s.nunique(dropna=True) > 1:
                feats.append(col)
        feature_sets.append(set(feats))
    common = sorted(set.intersection(*feature_sets)) if feature_sets else []
    return common


def anchor_week(raw: Dict[Tuple[str, int], pd.DataFrame], experiment: str) -> int:
    weeks = get_weeks(raw, experiment)
    return 0 if 0 in weeks else min(weeks)


def fit_preprocessor(
    raw: Dict[Tuple[str, int], pd.DataFrame],
    experiment: str,
    features: List[str],
    winsor_q: float,
    base_week: int,
) -> Dict[str, PreprocessParam]:
    base_df = raw[(experiment, base_week)]
    params: Dict[str, PreprocessParam] = {}

    for feature in features:
        s = pd.to_numeric(base_df[feature], errors="coerce").dropna().astype(float)
        if s.empty:
            raise ValueError(f"实验 {experiment} 的基线周 {base_week}W 中，指标 {feature} 没有有效数值。")
        use_log = bool((s.min() > 0) and (abs(skew(s)) > 1))
        t = np.log1p(s) if use_log else s.copy()

        median = float(np.median(t))
        mad = float(np.median(np.abs(t - median)))
        if mad < 1e-9:
            mad = float((t.quantile(0.75) - t.quantile(0.25)) / 1.349)
        if mad < 1e-9:
            std = float(t.std(ddof=1))
            mad = std if std > 1e-9 else 1.0

        lower = float(t.quantile(winsor_q))
        upper = float(t.quantile(1 - winsor_q))
        params[feature] = PreprocessParam(
            use_log=use_log,
            median=median,
            mad=mad,
            lower=lower,
            upper=upper,
        )
    return params


def apply_preprocessor(
    df: pd.DataFrame, features: List[str], params: Dict[str, PreprocessParam], winsor: bool
) -> pd.DataFrame:
    out = df[["Group", "Sample ID"]].copy()
    for feature in features:
        s = pd.to_numeric(df[feature], errors="coerce").astype(float)
        if params[feature].use_log:
            s = np.log1p(s)
        if winsor:
            s = s.clip(params[feature].lower, params[feature].upper)
        z = 0.6745 * (s - params[feature].median) / params[feature].mad
        out[feature] = z
    return out


def feature_inventory(raw: Dict[Tuple[str, int], pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for exp in get_experiments(raw):
        weeks = get_weeks(raw, exp)
        all_feats = {}
        common = set(common_features(raw, exp))
        for week in weeks:
            df = raw[(exp, week)]
            for col in df.columns[2:]:
                s = pd.to_numeric(df[col], errors="coerce")
                if s.notna().sum() > len(df) * 0.3:
                    all_feats.setdefault(col, []).append(week)
        for feat, feat_weeks in sorted(all_feats.items()):
            rows.append(
                {
                    "Experiment": exp,
                    "Feature": feat,
                    "Present_weeks": ", ".join([f"{w}W" for w in sorted(feat_weeks)]),
                    "Common_to_all_weeks": feat in common,
                    "Direction_higher_is_worse": infer_direction(feat),
                }
            )
    return pd.DataFrame(rows)


def id_consistency_summary(raw: Dict[Tuple[str, int], pd.DataFrame]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    summary_rows = []
    detail_rows = []

    for exp in get_experiments(raw):
        weeks = get_weeks(raw, exp)
        id_sets = {}
        id_group_seen: Dict[int, str] = {}
        inconsistent = []

        for week in weeks:
            df = raw[(exp, week)][["Sample ID", "Group"]].dropna().copy()
            if df.empty:
                continue
            df["Sample ID"] = df["Sample ID"].astype(int)
            id_sets[week] = set(df["Sample ID"].tolist())

            for _, row in df.iterrows():
                sid = int(row["Sample ID"])
                grp = row["Group"]
                if sid in id_group_seen and id_group_seen[sid] != grp:
                    inconsistent.append((sid, id_group_seen[sid], grp, week))
                else:
                    id_group_seen[sid] = grp

        common_ids = set.intersection(*id_sets.values()) if id_sets else set()
        union_ids = set.union(*id_sets.values()) if id_sets else set()

        summary_rows.append(
            {
                "Experiment": exp,
                "Weeks": ", ".join([f"{w}W" for w in weeks]),
                "IDs_common_across_all_weeks": len(common_ids),
                "IDs_union": len(union_ids),
                "Group_inconsistency_count_if_ID_treated_as_longitudinal": len(inconsistent),
                "Interpretation": (
                    "Sample ID 不宜直接当作稳定纵向个体 ID"
                    if len(inconsistent) > 0
                    else "未发现明显组别冲突，可进一步人工核对"
                ),
            }
        )
        for sid, old_grp, new_grp, week in inconsistent:
            detail_rows.append(
                {
                    "Experiment": exp,
                    "Sample ID": sid,
                    "Earlier_group": old_grp,
                    "Later_group": new_grp,
                    "Later_week": week,
                }
            )

    return pd.DataFrame(summary_rows), pd.DataFrame(detail_rows)


def _safe_group_color_map(groups: List[str]) -> Dict[str, str]:
    return {g: PLOT_COLORS[i % len(PLOT_COLORS)] for i, g in enumerate(groups)}


def _format_axis(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(alpha=0.25, linestyle="--", linewidth=0.6)


def save_table_bundle_excel(tables: Dict[str, pd.DataFrame], out_path: Path) -> None:
    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        for sheet_name, df in tables.items():
            safe_name = str(sheet_name)[:31]
            if df is None:
                pd.DataFrame().to_excel(writer, index=False, sheet_name=safe_name)
            else:
                df.to_excel(writer, index=False, sheet_name=safe_name)
        writer.book.save(out_path)

    wb = openpyxl.load_workbook(out_path)
    header_fill = PatternFill("solid", fgColor="1F4E78")
    header_font = Font(color="FFFFFF", bold=True)
    for ws in wb.worksheets:
        ws.freeze_panes = "A2"
        for cell in ws[1]:
            cell.fill = header_fill
            cell.font = header_font
        for col_idx, col in enumerate(ws.columns, start=1):
            max_len = 0
            for cell in col:
                value = "" if cell.value is None else str(cell.value)
                max_len = max(max_len, len(value))
            ws.column_dimensions[get_column_letter(col_idx)].width = min(max(12, max_len + 2), 40)
    wb.save(out_path)


# Result tables exported as separate workbooks (reading order). Only the tables
# needed to read / verify the results are included; verbose diagnostics
# (file/feature inventories, ID-consistency, preprocessing params, outlier stats,
# sample-level dumps, per-indicator mechanism tables, QC/coverage/attrition) are
# intentionally omitted from the exported bundle.
NECESSARY_TABLE_KEYS: List[str] = [
    "实验总览",
    "组时间点汇总",
    "治疗比较",
    "时间趋势",
    "实验内相关性",
    "跨实验合并相关性",
    "TrueNormal终点概览",
    "Display_Normalized_Index",
    "样本级结果",
]


def _sanitize_table_filename(name: str) -> str:
    text = str(name).strip()
    for ch in '\\/:*?"<>|':
        text = text.replace(ch, "_")
    return text or "table"


def save_curated_tables_separately(
    tables: Dict[str, pd.DataFrame],
    out_dir: Path,
    keys: Optional[List[str]] = None,
) -> List[Path]:
    """Write each necessary result table to its own single-sheet workbook.

    Replaces the old ~25-sheet combined workbook: each table lands in its own
    file (numbered by reading order) so the export is easy to browse, and only
    the curated ``NECESSARY_TABLE_KEYS`` are written.
    """
    keys = keys if keys is not None else NECESSARY_TABLE_KEYS
    out_dir.mkdir(parents=True, exist_ok=True)
    written: List[Path] = []
    for idx, key in enumerate(keys, start=1):
        df = tables.get(key)
        if not isinstance(df, pd.DataFrame) or df.empty:
            continue
        path = out_dir / f"{idx:02d}_{_sanitize_table_filename(key)}.xlsx"
        save_table_bundle_excel({key: df}, path)
        written.append(path)
    return written


def run_from_paths(
    input_paths: List[str],
    role_map: Optional[Dict[str, Dict[str, str]]] = None,
    winsor_q: float = WINSOR_Q_DEFAULT,
) -> Dict[str, Any]:
    raw, upload_overview = load_input_files(input_paths)
    results = analyze_all(raw, role_map=role_map, winsor_q=winsor_q)
    results["tables"]["上传文件概览"] = upload_overview
    bundle = generate_output_bundle(results)
    results["bundle"] = bundle
    return results


# --------------------------------------------------------------------------- #
# Publication-style formatting and the harmonized disease index
# --------------------------------------------------------------------------- #

def _format_p_publication(p: float) -> str:
    if pd.isna(p):
        return "p=NA"
    if p < 0.001:
        return "p<0.001"
    if p < 0.01:
        return f"p={p:.3f}"
    if p < 0.1:
        return f"p={p:.3f}"
    return f"p={p:.2f}"


def _sig_star(p: float) -> str:
    if pd.isna(p):
        return ""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def _ordered_groups_for_plot(sub: pd.DataFrame, treatment_group: str | None, comparator_group: str | None) -> list[str]:
    groups = sorted(sub["Group"].dropna().astype(str).unique().tolist())
    ordered = []
    for g in [comparator_group, treatment_group]:
        if g and g in groups and g not in ordered:
            ordered.append(g)
    for g in groups:
        if g not in ordered:
            ordered.append(g)
    return ordered


def _group_color_map(groups: list[str], treatment_group: str | None = None, comparator_group: str | None = None) -> dict[str, str]:
    cmap = {}
    for g in groups:
        if treatment_group and g == treatment_group:
            cmap[g] = "#1f77b4"
        elif comparator_group and g == comparator_group:
            cmap[g] = "#d62728"
    other_colors = ["#2ca02c", "#9467bd", "#ff7f0e", "#8c564b", "#17becf", "#7f7f7f"]
    j = 0
    for g in groups:
        if g not in cmap:
            cmap[g] = other_colors[j % len(other_colors)]
            j += 1
    return cmap


def plot_entropy_trajectory(
    summary_df: pd.DataFrame,
    experiment: str,
    out_path: Path,
    comparison_df: pd.DataFrame | None = None,
    treatment_group: str | None = None,
    comparator_group: str | None = None,
) -> None:
    _draw_publication_trajectory(
        summary_df=summary_df,
        comparison_df=comparison_df,
        experiment=experiment,
        metric="Entropy",
        treatment_group=treatment_group,
        comparator_group=comparator_group,
        out_path=out_path,
    )


def plot_hdi_trajectory(
    summary_df: pd.DataFrame,
    experiment: str,
    out_path: Path,
    comparison_df: pd.DataFrame | None = None,
    treatment_group: str | None = None,
    comparator_group: str | None = None,
) -> None:
    _draw_publication_trajectory(
        summary_df=summary_df,
        comparison_df=comparison_df,
        experiment=experiment,
        metric="HDI_mean",
        treatment_group=treatment_group,
        comparator_group=comparator_group,
        out_path=out_path,
    )


def plot_outlier_counts(outlier_df: pd.DataFrame, experiment: str, out_path: Path) -> None:
    """
    Quality-control figure, not an efficacy figure.
    Shows the number of samples flagged by >=2 methods among:
    robust z / IQR / Isolation Forest / LOF / robust Mahalanobis.
    """
    sub = outlier_df[outlier_df["Experiment"] == experiment].copy()
    if sub.empty:
        return
    groups = sorted(sub["Group"].unique().tolist())
    cmap = _safe_group_color_map(groups)

    fig, ax = plt.subplots(figsize=(8.8, 5.3))
    for idx, group in enumerate(groups):
        g = sub[sub["Group"] == group].sort_values("Week")
        label = display_group_name(group, idx)
        ax.plot(g["Week"], g["n_consensus"], marker="o", linewidth=2.0, label=label, color=cmap[group])
        for x0, y0 in zip(g["Week"], g["n_consensus"]):
            ax.annotate(f"{int(y0)}", (x0, y0), xytext=(0, 7), textcoords="offset points", ha="center", fontsize=8, color=cmap[group])
    ax.set_title(f"Experiment {experiment}: samples flagged by ≥2 outlier methods", fontsize=12.5, pad=12)
    ax.set_xlabel("Week")
    ax.set_ylabel("Flagged sample count")
    _format_axis(ax)
    ax.legend(frameon=False)
    ax.text(
        0.01, 0.01,
        "QC metric only. Consensus outlier = flagged by at least two methods.",
        transform=ax.transAxes, ha="left", va="bottom", fontsize=8.4,
        bbox=dict(boxstyle="round,pad=0.22", facecolor="white", edgecolor="#dddddd", alpha=0.92)
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


# Better fallback labels when Chinese fonts are unavailable.
def display_group_name(name: str, fallback_index: Optional[int] = None) -> str:
    text = "" if pd.isna(name) else str(name)
    if HAS_CJK_FONT:
        return text
    if all(ord(ch) < 128 for ch in text):
        return text
    lowered = text.lower()
    for key, val in GROUP_TRANSLATIONS.items():
        if key.lower() in lowered or key in text:
            return val
    digit_text = "".join(ch for ch in text if ch.isdigit())
    if digit_text:
        return digit_text
    ascii_letters = "".join(ch for ch in text if ord(ch) < 128 and (ch.isalnum() or ch in "-_ "))
    if ascii_letters.strip():
        return ascii_letters.strip()
    return f"Group {fallback_index + 1}" if fallback_index is not None else "Group"


# --------------------------------------------------------------------------- #
# Normal-reference fitting and the mechanism-level summaries
# --------------------------------------------------------------------------- #

NORMAL_REF_NONE = "（无）"

def fit_reference_params(
    raw: Dict[Tuple[str, int], pd.DataFrame],
    experiment: str,
    features: List[str],
    reference_group: str,
) -> Dict[str, PreprocessParam]:
    """Fit robust healthy-reference envelopes from the optional Normal group.

    If the experiment contains a healthy/normal reference group, we pool that group
    across all available weeks to obtain a stable central envelope for each marker.
    This yields an *absolute reference* that is distinct from the anchor-week
    relative score used when no healthy reference exists.
    """
    ref_frames = []
    for week in get_weeks(raw, experiment):
        df = raw[(experiment, week)]
        sub = df[df["Group"] == reference_group]
        if not sub.empty:
            ref_frames.append(sub)
    if not ref_frames:
        raise ValueError(f"实验 {experiment} 中未找到参考组 {reference_group} 的有效数据。")

    ref_df = pd.concat(ref_frames, ignore_index=True)
    params: Dict[str, PreprocessParam] = {}
    for feature in features:
        s = pd.to_numeric(ref_df[feature], errors="coerce").dropna().astype(float)
        if s.empty:
            raise ValueError(f"实验 {experiment} 的参考组 {reference_group} 中，指标 {feature} 没有有效数值。")
        use_log = bool((s.min() > 0) and (abs(skew(s)) > 1))
        t = np.log1p(s) if use_log else s.copy()

        median = float(np.median(t))
        mad = float(np.median(np.abs(t - median)))
        if mad < 1e-9:
            mad = float((t.quantile(0.75) - t.quantile(0.25)) / 1.349)
        if mad < 1e-9:
            std = float(t.std(ddof=1))
            mad = std if std > 1e-9 else 1.0

        if len(t) >= 40:
            lower = float(t.quantile(0.025))
            upper = float(t.quantile(0.975))
        else:
            q1 = float(t.quantile(0.25))
            q3 = float(t.quantile(0.75))
            iqr = q3 - q1
            if iqr < 1e-9:
                iqr = mad * 1.349
            lower = float(q1 - 1.5 * iqr)
            upper = float(q3 + 1.5 * iqr)

        if (not np.isfinite(lower)) or (not np.isfinite(upper)) or lower >= upper:
            lower = float(median - 2.5 * mad)
            upper = float(median + 2.5 * mad)

        params[feature] = PreprocessParam(
            use_log=use_log,
            median=median,
            mad=max(mad, 1e-9),
            lower=lower,
            upper=upper,
        )
    return params


def compute_normal_reference_score(
    df: pd.DataFrame,
    features: List[str],
    ref_params: Dict[str, PreprocessParam],
) -> pd.Series:
    """Compute a healthy-reference proximity score on [0, 1].

    Score = 1 if a sample falls within the robust Normal envelope for every
    available feature; it approaches 0 as the sample moves farther outside the
    reference band. This is a practical, interpretable *absolute* health
    proximity score when a Normal group exists.
    """
    vals = []
    for _, row in df.iterrows():
        dists = []
        for feature in features:
            x = pd.to_numeric(pd.Series([row.get(feature)]), errors="coerce").iloc[0]
            if pd.isna(x):
                continue
            p = ref_params[feature]
            if p.use_log:
                if x <= -1:
                    continue
                x = np.log1p(float(x))
            else:
                x = float(x)

            band = max((p.upper - p.lower) / 2.0, p.mad, 1e-9)
            if x < p.lower:
                dist = (p.lower - x) / band
            elif x > p.upper:
                dist = (x - p.upper) / band
            else:
                dist = 0.0
            dists.append(float(dist))

        if len(dists) == 0:
            vals.append(np.nan)
        else:
            mean_dist = float(np.mean(dists))
            vals.append(float(1.0 / (1.0 + mean_dist)))
    return pd.Series(vals, index=df.index, name="NHPS")


def compute_bri(z_df: pd.DataFrame, features: List[str]) -> pd.Series:
    """Baseline-relative burden index.

    This is the old HDI logic under a clearer name. It measures direction-aware
    deviation from the anchor week, not absolute health.
    """
    return compute_hdi(z_df, features).rename("BRI")


def _feature_mechanism_table(
    z_df: pd.DataFrame,
    experiment: str,
    week: int,
    features: List[str],
) -> pd.DataFrame:
    rows = []
    for group, sub in z_df.groupby("Group"):
        for feature in features:
            direction = infer_direction(feature)
            s = pd.to_numeric(sub[feature], errors="coerce").astype(float)
            valid = s.dropna().to_numpy(dtype=float)
            if len(valid) == 0:
                signed_mean = float("nan")
                var_log = float("nan")
                mean_abs = float("nan")
            else:
                signed_mean = float(np.nanmean(valid * direction))
                mean_abs = float(np.nanmean(np.abs(valid)))
                if len(valid) >= 2:
                    var_log = float(0.5 * np.log(max(np.nanvar(valid, ddof=1), 1e-8)))
                else:
                    var_log = float("nan")
            rows.append(
                {
                    "Experiment": experiment,
                    "Week": week,
                    "Group": group,
                    "Feature": feature,
                    "Direction_higher_is_worse": infer_direction(feature),
                    "Signed_burden_mean": signed_mean,
                    "Mean_abs_z": mean_abs,
                    "Variance_log_component": var_log,
                }
            )
    return pd.DataFrame(rows)


def _feature_comparison_table(
    feature_table: pd.DataFrame,
    experiment: str,
    treatment_group: str,
    comparator_group: str,
) -> pd.DataFrame:
    rows = []
    if feature_table.empty:
        return pd.DataFrame()
    for week in sorted(feature_table["Week"].dropna().unique().tolist()):
        t = feature_table[
            (feature_table["Experiment"] == experiment) &
            (feature_table["Week"] == week) &
            (feature_table["Group"] == treatment_group)
        ][["Feature", "Signed_burden_mean", "Variance_log_component"]].rename(
            columns={
                "Signed_burden_mean": "Signed_burden_treat",
                "Variance_log_component": "VarLog_treat",
            }
        )
        c = feature_table[
            (feature_table["Experiment"] == experiment) &
            (feature_table["Week"] == week) &
            (feature_table["Group"] == comparator_group)
        ][["Feature", "Signed_burden_mean", "Variance_log_component"]].rename(
            columns={
                "Signed_burden_mean": "Signed_burden_comp",
                "Variance_log_component": "VarLog_comp",
            }
        )
        m = t.merge(c, on="Feature", how="inner")
        if m.empty:
            continue
        m["Experiment"] = experiment
        m["Week"] = week
        m["Treatment_group"] = treatment_group
        m["Comparator_group"] = comparator_group
        m["Burden_diff_treat_minus_comp"] = m["Signed_burden_treat"] - m["Signed_burden_comp"]
        m["VarLog_diff_treat_minus_comp"] = m["VarLog_treat"] - m["VarLog_comp"]
        rows.append(m)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def _metric_spec_base(metric: str) -> dict:
    if metric == "Entropy":
        return {
            "value": "Entropy",
            "low": "Entropy_ci_low",
            "high": "Entropy_ci_high",
            "p": "Entropy_p_boot",
            "diff": "Entropy_diff_treat_minus_comp",
            "ylabel": "Entropy",
            "title": "Entropy trajectory",
            "note": "Points are entropy estimates; error bars are 95% bootstrap CIs.",
        }
    if metric == "NHPS_mean":
        return {
            "value": "NHPS_mean",
            "low": "NHPS_mean_ci_low",
            "high": "NHPS_mean_ci_high",
            "p": "NHPS_p_mannwhitney",
            "diff": "NHPS_diff_treat_minus_comp",
            "ylabel": "Normal-reference proximity score",
            "title": "Normal-reference proximity trajectory",
            "note": "Only available when a Normal/healthy reference group is provided. Higher values mean closer to the robust Normal reference envelope.",
        }
    return {
        "value": "BRI_mean",
        "low": "BRI_mean_ci_low",
        "high": "BRI_mean_ci_high",
        "p": "BRI_p_mannwhitney",
        "diff": "BRI_diff_treat_minus_comp",
        "ylabel": "Baseline-relative burden index (BRI)",
        "title": "Baseline-relative burden trajectory",
        "note": "BRI is anchored to the baseline week. Lower values indicate improvement relative to the anchor; it is not an absolute health scale.",
    }


def plot_relative_burden_trajectory(
    summary_df: pd.DataFrame,
    experiment: str,
    out_path: Path,
    comparison_df: pd.DataFrame | None = None,
    treatment_group: str | None = None,
    comparator_group: str | None = None,
) -> None:
    _draw_publication_trajectory(
        summary_df=summary_df,
        comparison_df=comparison_df,
        experiment=experiment,
        metric="BRI_mean",
        treatment_group=treatment_group,
        comparator_group=comparator_group,
        out_path=out_path,
    )


def plot_normal_reference_trajectory(
    summary_df: pd.DataFrame,
    experiment: str,
    out_path: Path,
    comparison_df: pd.DataFrame | None = None,
    treatment_group: str | None = None,
    comparator_group: str | None = None,
) -> None:
    sub = summary_df[(summary_df["Experiment"] == experiment) & (summary_df["NHPS_mean"].notna())].copy()
    if sub.empty:
        return
    _draw_publication_trajectory(
        summary_df=summary_df,
        comparison_df=comparison_df,
        experiment=experiment,
        metric="NHPS_mean",
        treatment_group=treatment_group,
        comparator_group=comparator_group,
        out_path=out_path,
    )


def plot_health_entropy_scatter(summary_df: pd.DataFrame, experiment: str, out_path: Path) -> None:
    sub = summary_df[summary_df["Experiment"] == experiment].copy()
    if sub.empty:
        return
    groups = sorted(sub["Group"].unique().tolist())
    cmap = _safe_group_color_map(groups)

    fig, ax = plt.subplots(figsize=(7.3, 5.6))
    for idx, group in enumerate(groups):
        g = sub[sub["Group"] == group].copy()
        label = display_group_name(group, idx)
        ax.scatter(g["BRI_mean"], g["Entropy"], s=72, alpha=0.92, color=cmap[group], label=label)
        for _, row in g.iterrows():
            ax.annotate(f"{int(row['Week'])}W", (row["BRI_mean"], row["Entropy"]), fontsize=8, xytext=(4, 3), textcoords="offset points")
    if len(sub) >= 3:
        rho, pval = spearmanr(sub["BRI_mean"], sub["Entropy"])
        ax.text(
            0.02, 0.98, f"Spearman rho = {rho:.3f}\np = {pval:.3g}",
            transform=ax.transAxes, va="top", ha="left", fontsize=10,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.9, edgecolor="#cccccc")
        )
    ax.set_title(f"Experiment {experiment}: BRI vs Entropy")
    ax.set_xlabel("Baseline-relative burden index (BRI)")
    ax.set_ylabel("Entropy")
    _format_axis(ax)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Font-safe static figures: feature-name translation and label handling
# --------------------------------------------------------------------------- #
import re as _re_text

FEATURE_TRANSLATIONS = {
    "空腹血糖": "Fasting glucose",
    "血糖": "Glucose",
    "葡萄糖": "Glucose",
    "甘油三酯": "Triglycerides",
    "总胆固醇": "Total cholesterol",
    "胆固醇": "Cholesterol",
    "低密度脂蛋白胆固醇": "LDL-C",
    "高密度脂蛋白胆固醇": "HDL-C",
    "非高密度脂蛋白胆固醇": "non-HDL-C",
    "极低密度脂蛋白胆固醇": "VLDL-C",
    "载脂蛋白A1": "ApoA1",
    "载脂蛋白B": "ApoB",
    "游离脂肪酸": "Free fatty acids",
    "糖化血红蛋白": "HbA1c",
    "胰岛素": "Insulin",
    "脂联素": "Adiponectin",
    "瘦素": "Leptin",
    "谷丙转氨酶": "ALT",
    "丙氨酸氨基转移酶": "ALT",
    "谷草转氨酶": "AST",
    "天门冬氨酸氨基转移酶": "AST",
    "γ-谷氨酰转移酶": "GGT",
    "γ-谷氨酰转肽酶": "GGT",
    "谷氨酰转移酶": "GGT",
    "碱性磷酸酶": "ALP",
    "乳酸脱氢酶": "LDH",
    "肌酸激酶": "CK",
    "尿酸": "Uric acid",
    "尿素": "Urea",
    "尿素氮": "BUN",
    "肌酐": "Creatinine",
    "总蛋白": "Total protein",
    "白蛋白": "Albumin",
    "球蛋白": "Globulin",
    "总胆红素": "Total bilirubin",
    "直接胆红素": "Direct bilirubin",
    "间接胆红素": "Indirect bilirubin",
    "C反应蛋白": "CRP",
    "血红蛋白": "Hemoglobin",
    "红细胞": "RBC",
    "白细胞介素-1β": "IL-1beta",
    "白细胞介素-1": "IL-1",
    "白细胞介素-6": "IL-6",
    "肿瘤坏死因子-α": "TNF-alpha",
    "肿瘤坏死因子": "TNF",
    "脂蛋白a": "Lp(a)",
    "脂蛋白A": "Lp(a)",
    "C肽": "C-peptide",
    "糖化": "HbA1c",
    "白细胞": "WBC",
    "血小板": "Platelets",
    "钠": "Na",
    "钾": "K",
    "氯": "Cl",
    "钙": "Ca",
    "磷": "P",
    "镁": "Mg",
}

_ASCII_REPLACEMENTS = {
    "μ": "u",
    "µ": "u",
    "γ": "gamma",
    "β": "beta",
    "α": "alpha",
    "×": "x",
    "·": ".",
    "－": "-",
    "—": "-",
    "–": "-",
    "％": "%",
    "／": "/",
    "（": "(",
    "）": ")",
}


def _to_ascii_safe(text: Any) -> str:
    s = "" if pd.isna(text) else str(text)
    for old, new in _ASCII_REPLACEMENTS.items():
        s = s.replace(old, new)
    s = s.encode("ascii", "ignore").decode("ascii")
    s = _re_text.sub(r"\s+", " ", s).strip()
    return s


def _extract_units_ascii(feature_name: Any) -> str:
    text = "" if pd.isna(feature_name) else str(feature_name)
    m = _re_text.findall(r"[\(（]([^\)）]+)[\)）]", text)
    if not m:
        return ""
    return _to_ascii_safe(m[-1])


def _strip_units_text(feature_name: Any) -> str:
    text = "" if pd.isna(feature_name) else str(feature_name)
    text = _re_text.sub(r"\s*[\(（][^\)）]*[\)）]\s*", "", text)
    return _re_text.sub(r"\s+", " ", text).strip()


def feature_display_label(feature_name: Any, index: Optional[int] = None) -> tuple[str, str]:
    original = "" if pd.isna(feature_name) else str(feature_name).strip()
    if HAS_CJK_FONT:
        return original, "original"
    if all(ord(ch) < 128 for ch in original):
        return original, "ascii_original"

    base = _strip_units_text(original)
    unit = _extract_units_ascii(original)

    alias = None
    for key, value in sorted(FEATURE_TRANSLATIONS.items(), key=lambda kv: len(kv[0]), reverse=True):
        if key in base or key in original:
            alias = value
            break

    if alias is not None:
        label = alias
        strategy = "english_alias"
    else:
        ascii_part = _to_ascii_safe(base)
        alnum_only = "".join(ch for ch in ascii_part if ch.isalnum())
        if ascii_part and len(alnum_only) >= 3 and ascii_part.lower() not in {"x", "c", "a", "alpha", "beta", "gamma"}:
            label = ascii_part
            strategy = "ascii_cleanup"
        else:
            code = f"F{(index or 0) + 1:02d}" if index is not None else "Feature"
            label = code
            strategy = "feature_code"

    if unit and unit.lower() not in label.lower():
        label = f"{label} ({unit})"
    return label, strategy


def _label_mapping_for_experiment_from_features(experiment: str, features: list[str]) -> pd.DataFrame:
    rows = []
    for idx, feat in enumerate(sorted(dict.fromkeys([str(f) for f in features]))):
        label, strategy = feature_display_label(feat, idx)
        rows.append({
            "Experiment": experiment,
            "Original_feature": feat,
            "Plot_label": label,
            "Label_strategy": strategy,
            "Units_ascii": _extract_units_ascii(feat),
            "Has_non_ascii": not all(ord(ch) < 128 for ch in feat),
            "Matplotlib_font": MPL_FONT_NAME,
            "Has_CJK_font": HAS_CJK_FONT,
        })
    return pd.DataFrame(rows)


def _build_plot_label_mapping_base(tables: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    cmp_df = tables.get("指标比较分解")
    inv_df = tables.get("特征清单")
    experiments = set()
    if isinstance(cmp_df, pd.DataFrame) and not cmp_df.empty and "Experiment" in cmp_df.columns:
        experiments |= set(cmp_df["Experiment"].dropna().astype(str).tolist())
    if isinstance(inv_df, pd.DataFrame) and not inv_df.empty and "Experiment" in inv_df.columns:
        experiments |= set(inv_df["Experiment"].dropna().astype(str).tolist())

    for exp in sorted(experiments):
        feats = []
        if isinstance(cmp_df, pd.DataFrame) and not cmp_df.empty:
            feats.extend(cmp_df.loc[cmp_df["Experiment"].astype(str) == exp, "Feature"].dropna().astype(str).tolist())
        if isinstance(inv_df, pd.DataFrame) and not inv_df.empty:
            feats.extend(inv_df.loc[inv_df["Experiment"].astype(str) == exp, "Feature"].dropna().astype(str).tolist())
        if feats:
            rows.append(_label_mapping_for_experiment_from_features(exp, feats))

    if rows:
        return pd.concat(rows, ignore_index=True)
    return pd.DataFrame(columns=[
        "Experiment", "Original_feature", "Plot_label", "Label_strategy",
        "Units_ascii", "Has_non_ascii", "Matplotlib_font", "Has_CJK_font"
    ])


def _feature_label_map_for_sub(feature_sub: pd.DataFrame, experiment: str) -> tuple[dict[str, str], pd.DataFrame]:
    all_feats = feature_sub.loc[feature_sub["Experiment"].astype(str) == str(experiment), "Feature"].dropna().astype(str).tolist()
    mapping_df = _label_mapping_for_experiment_from_features(experiment, all_feats)
    if mapping_df.empty:
        return {}, mapping_df
    return dict(zip(mapping_df["Original_feature"], mapping_df["Plot_label"])), mapping_df


def _write_plot_label_mapping_files(label_map_df: pd.DataFrame, output_dir: Path) -> None:
    if label_map_df is None or label_map_df.empty:
        return
    csv_path = output_dir / "Plot_Label_Mapping.csv"
    label_map_df.to_csv(csv_path, index=False, encoding="utf-8-sig")

    html_path = output_dir / "Plot_Label_Mapping.html"
    title = "BioEntropy plot label mapping"
    note = (
        "This file records how feature names were rendered in static figures. "
        "When a Chinese-capable font is unavailable, the software falls back to English aliases or feature codes "
        "to avoid garbled labels in exported PNG files."
    )
    html = [
        "<!DOCTYPE html>",
        "<html lang='en'><head><meta charset='utf-8'>",
        f"<title>{title}</title>",
        "<style>body{font-family:Arial,Helvetica,sans-serif;margin:28px;color:#222;line-height:1.6;} h1{color:#1f4e78;} table{border-collapse:collapse;width:100%;} th,td{border:1px solid #d0d7de;padding:8px 10px;text-align:left;} th{background:#eef5fb;} .box{background:#f7fbff;border-left:4px solid #1f78b4;padding:12px 14px;margin:18px 0;}</style>",
        "</head><body>",
        f"<h1>{title}</h1>",
        f"<div class='box'>{note}</div>",
        label_map_df.to_html(index=False, escape=True),
        "</body></html>",
    ]
    html_path.write_text("\n".join(html), encoding="utf-8")


def plot_feature_repair_heatmap(
    feature_cmp_df: pd.DataFrame,
    experiment: str,
    out_path: Path,
    top_n: int = 12,
) -> None:
    sub = feature_cmp_df[feature_cmp_df["Experiment"] == experiment].copy()
    if sub.empty:
        return
    label_map, map_df = _feature_label_map_for_sub(feature_cmp_df, experiment)
    week_max = sub["Week"].max()
    rank = (
        sub[sub["Week"] == week_max]
        .assign(abs_rank=lambda d: d["Burden_diff_treat_minus_comp"].abs())
        .sort_values("abs_rank", ascending=False)
    )
    features = rank["Feature"].head(top_n).tolist()
    if not features:
        return
    pivot = (
        sub[sub["Feature"].isin(features)]
        .pivot_table(index="Feature", columns="Week", values="Burden_diff_treat_minus_comp", aggfunc="mean")
        .reindex(features)
    )
    fig, ax = plt.subplots(figsize=(9.4, max(4.8, 0.46 * len(features) + 2.2)))
    mat = pivot.to_numpy(dtype=float)
    vmax = np.nanmax(np.abs(mat)) if np.isfinite(mat).any() else 1.0
    vmax = max(vmax, 0.25)
    im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ylabels = [label_map.get(feat, feat) for feat in features]
    ax.set_yticks(range(len(features)))
    ax.set_yticklabels(ylabels, fontsize=8.6)
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels([f"{int(x)}W" for x in pivot.columns], fontsize=9)
    ax.set_title(
        f"Experiment {experiment}: feature repair heatmap\n(treatment - comparator; lower = lower burden in treatment)",
        fontsize=12.5, pad=12
    )
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if np.isfinite(mat[i, j]):
                ax.text(j, i, f"{mat[i, j]:+.2f}", ha="center", va="center", fontsize=7.2, color="black")
    cbar = fig.colorbar(im, ax=ax, shrink=0.9)
    cbar.set_label("Signed burden difference", fontsize=9)
    if not HAS_CJK_FONT:
        ax.text(
            0.01, -0.14,
            "CJK font unavailable: exported PNG uses English aliases or feature codes. See Plot_Label_Mapping sheet / CSV.",
            transform=ax.transAxes, ha="left", va="top", fontsize=8.2
        )
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_variance_stabilization_heatmap(
    feature_cmp_df: pd.DataFrame,
    experiment: str,
    out_path: Path,
    top_n: int = 12,
) -> None:
    sub = feature_cmp_df[feature_cmp_df["Experiment"] == experiment].copy()
    if sub.empty:
        return
    label_map, map_df = _feature_label_map_for_sub(feature_cmp_df, experiment)
    week_max = sub["Week"].max()
    rank = (
        sub[sub["Week"] == week_max]
        .assign(abs_rank=lambda d: d["VarLog_diff_treat_minus_comp"].abs())
        .sort_values("abs_rank", ascending=False)
    )
    features = rank["Feature"].head(top_n).tolist()
    if not features:
        return
    pivot = (
        sub[sub["Feature"].isin(features)]
        .pivot_table(index="Feature", columns="Week", values="VarLog_diff_treat_minus_comp", aggfunc="mean")
        .reindex(features)
    )
    fig, ax = plt.subplots(figsize=(9.4, max(4.8, 0.46 * len(features) + 2.3)))
    mat = pivot.to_numpy(dtype=float)
    vmax = np.nanmax(np.abs(mat)) if np.isfinite(mat).any() else 1.0
    vmax = max(vmax, 0.25)
    im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ylabels = [label_map.get(feat, feat) for feat in features]
    ax.set_yticks(range(len(features)))
    ax.set_yticklabels(ylabels, fontsize=8.6)
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels([f"{int(x)}W" for x in pivot.columns], fontsize=9)
    ax.set_title(
        f"Experiment {experiment}: variance-stabilization heatmap\n(treatment - comparator; lower = lower feature variance in treatment)",
        fontsize=12.5, pad=12
    )
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if np.isfinite(mat[i, j]):
                ax.text(j, i, f"{mat[i, j]:+.2f}", ha="center", va="center", fontsize=7.2, color="black")
    cbar = fig.colorbar(im, ax=ax, shrink=0.9)
    cbar.set_label("Approx. log-variance difference", fontsize=9)
    note = (
        "Approximation note: entropy is multivariate and also depends on covariance. "
        "This plot uses per-feature log-variance differences to summarize which markers most likely drive entropy change."
    )
    if not HAS_CJK_FONT:
        note += " CJK font unavailable: see Plot_Label_Mapping for axis labels."
    ax.text(0.01, -0.16, note, transform=ax.transAxes, ha="left", va="top", fontsize=8.0)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _build_markdown_report_base(results: Dict[str, Any]) -> str:
    tables = dict(results["tables"])
    overview = tables["实验总览"]
    compare_all = tables["治疗比较"]
    trends_all = tables["时间趋势"]
    summary_all = tables["组时间点汇总"]
    corr = tables["跨实验合并相关性"]
    label_map_df = build_plot_label_mapping_from_tables(tables)

    lines: List[str] = []
    lines.append("# BioEntropy 自动分析报告")
    lines.append("")
    lines.append("## 1. 指标体系说明")
    lines.append("")
    lines.append("- **Entropy**：组级系统熵，反映多指标系统的离散度/无序度。")
    lines.append("- **BRI（Baseline-relative burden index）**：锚定基线相对负担指数。它由方向校正后的稳健 z-score 聚合得到，反映样本相对锚定周的改善或恶化；**它不是绝对健康分数**。")
    lines.append("- **NHPS（Normal-reference health proximity score）**：只有在指定 Normal/健康参考组时才计算。它表示样本与稳健 Normal 参考包络的贴近程度，数值越高越接近参考健康状态。")
    lines.append("- 当前软件默认不会自动删除样本；异常值模块用于 QC 标记。winsor_q 仅对极端尾部做稳健截尾，降低单个极端值对结果的支配。")
    lines.append("")
    lines.append("## 2. 图形标签与字体策略")
    lines.append("")
    if HAS_CJK_FONT:
        lines.append(f"- 当前运行环境检测到可用中文字体：`{MPL_FONT_NAME}`。静态图优先使用原始中文标签。\n")
    else:
        lines.append(f"- 当前运行环境未检测到可用中文字体，matplotlib 回退为 `{MPL_FONT_NAME}`。为避免 PNG 图出现乱码或方框字，静态图会自动改用英文别名或特征代码。")
        lines.append("- 对应关系已写入 `图标签映射` 工作表，以及 `Plot_Label_Mapping.csv` / `Plot_Label_Mapping.html`。")
    lines.append("")
    if not overview.empty:
        lines.append("## 3. 各实验模式")
        lines.append("")
        for _, row in overview.iterrows():
            lines.append(
                f"- 实验 **{row['Experiment']}**：时间点 {row['Weeks']}；共同指标 {int(row['Common_features'])} 个；"
                f"锚定周 {int(row['Anchor_week'])}W；治疗组 `{row['Treatment_group']}`；比较组 `{row['Comparator_group']}`；"
                f"健康参考组 `{row['Reference_group']}`；模式 `{row['Score_mode']}`。"
            )
        lines.append("")
    if not corr.empty:
        r = corr.iloc[0]
        lines.append("## 4. 跨实验总体相关")
        lines.append("")

        def _assoc_line(label: str, rho_col: str, p_col: str) -> Optional[str]:
            if rho_col not in r.index or pd.isna(r.get(rho_col, np.nan)):
                return None
            q = r.get(p_col + "_fdr_bh", np.nan)
            q_txt = f"，FDR q = **{q:.3g}**" if pd.notna(q) else ""
            return f"- {label} 与 entropy：rho = **{r[rho_col]:.3f}**, p = **{r[p_col]:.3g}**{q_txt}。"

        for ln in (
            _assoc_line("BRI", "Spearman_rho_BRI", "BRI_p_value"),
            _assoc_line("NMD", "Spearman_rho_NMD", "NMD_p_value"),
            _assoc_line("NRBS", "Spearman_rho_NRBS", "NRBS_p_value"),
            _assoc_line("NHPS", "Spearman_rho_NHPS", "NHPS_p_value"),
        ):
            if ln:
                lines.append(ln)
        lines.append(
            "- 注：p 值已用 Benjamini–Hochberg FDR 校正（同一族内的疾病偏离–熵关联检验），"
            "请优先参考 q 值。相关性在组级别计算，组数很少，属探索性结果，不能替代受试者层面的预设统计分析。"
        )
        lines.append("")

    # Entropy reliability caveat: the log-determinant covariance entropy is
    # fragile when the per-group sample size n is not much larger than the
    # feature dimension p. Surface how many groups fall in each reliability bin.
    if not summary_all.empty and "Entropy_cov_reliability" in summary_all.columns:
        rel = summary_all["Entropy_cov_reliability"].value_counts().to_dict()
        n_low = int(rel.get("low", 0))
        n_limited = int(rel.get("limited", 0))
        n_reliable = int(rel.get("reliable", 0))
        total = int(summary_all.shape[0])
        if (n_low + n_limited) > 0:
            lines.append("## 4b. 熵估计可靠性（n/p）提示")
            lines.append("")
            lines.append(
                f"- 协方差对数行列式熵需要每组样本量 n 显著大于指标维度 p。本次共 {total} 个组-时间点中，"
                f"可靠（n/p≥5）{n_reliable} 个、有限（2≤n/p<5）{n_limited} 个、偏弱（n/p<2）{n_low} 个。"
            )
            lines.append(
                "- 对于“偏弱/有限”的组，协方差熵主要受 Ledoit–Wolf 收缩先验支配，组间熵比较应谨慎；"
                "在有真实 Normal 样本时建议优先参考 normal-reference state entropy 与 NMD。"
                "每组的 n、p、n/p、收缩强度与可靠性标签见“组时间点汇总”表。"
            )
            lines.append("")
    for exp in sorted(overview["Experiment"].tolist() if not overview.empty else []):
        exp_cmp = compare_all[compare_all["Experiment"] == exp].copy() if not compare_all.empty else pd.DataFrame()
        exp_trend = trends_all[trends_all["Experiment"] == exp].copy() if not trends_all.empty else pd.DataFrame()
        exp_info = overview[overview["Experiment"] == exp].iloc[0]
        lines.append(f"## 5. 实验 {exp} 摘要")
        lines.append("")
        if exp_info["Reference_group"] == NORMAL_REF_NONE:
            lines.append("- 本实验未指定 Normal/健康参考组，因此**不输出绝对健康程度判断**；请将 BRI 解释为“相对锚定周的负担变化”，而不是“离健康有多远”。")
        else:
            lines.append(f"- 本实验已指定健康参考组 `{exp_info['Reference_group']}`，可同时参考 NHPS。")
        if not exp_cmp.empty:
            for _, row in exp_cmp.sort_values("Week").iterrows():
                lines.append(
                    f"- {int(row['Week'])}W：BRI差值(治疗-比较) = {row['BRI_diff_treat_minus_comp']:+.3f} "
                    f"[{row['BRI_diff_ci_low']:+.3f}, {row['BRI_diff_ci_high']:+.3f}], "
                    f"p = {row['BRI_p_mannwhitney']:.3g}; "
                    f"Entropy差值 = {row['Entropy_diff_treat_minus_comp']:+.3f} "
                    f"[{row['Entropy_diff_ci_low']:+.3f}, {row['Entropy_diff_ci_high']:+.3f}], "
                    f"p = {row['Entropy_p_boot']:.3g}."
                )
                if pd.notna(row.get("NHPS_diff_treat_minus_comp", np.nan)):
                    lines.append(
                        f"  - NHPS差值(治疗-比较) = {row['NHPS_diff_treat_minus_comp']:+.3f} "
                        f"[{row['NHPS_diff_ci_low']:+.3f}, {row['NHPS_diff_ci_high']:+.3f}], "
                        f"p = {row['NHPS_p_mannwhitney']:.3g}."
                    )
        if not exp_trend.empty:
            lines.append("- 时间趋势：")
            for _, row in exp_trend.iterrows():
                msg = (
                    f"  - {row['Group']}：Entropy slope = {row['Entropy_linear_slope_per_week']:.4f}/week；"
                    f"BRI slope = {row['BRI_linear_slope_per_week']:.4f}/week。"
                )
                if pd.notna(row.get("NHPS_linear_slope_per_week", np.nan)):
                    msg += f" NHPS slope = {row['NHPS_linear_slope_per_week']:.4f}/week。"
                lines.append(msg)
        lines.append("")
    lines.append("## 6. 图形与机理说明")
    lines.append("")
    lines.append("- **relative_burden_trajectory**：显示锚定基线相对负担指数 BRI 的时间变化。")
    lines.append("- **normal_proximity_trajectory**：仅在存在健康参考组时生成，显示 NHPS 的时间变化。")
    lines.append("- **feature_repair_heatmap**：基于方向校正后的组均值差异，显示治疗组相对比较组在哪些指标上负担更低。")
    lines.append("- **variance_stabilization_heatmap**：基于每个指标的近似 log-variance 差异，概括哪些指标更可能驱动降熵。该图用于机理理解，是近似分解，不是严格熵加和分解。")
    lines.append("- **consensus_outlier_counts**：显示被 ≥2 种异常值算法同时标记的样本数量，仅用于质量控制。")
    if label_map_df is not None and not label_map_df.empty:
        lines.append("- 若静态图启用了英文别名或特征代码，请结合 `图标签映射` 工作表或 `Plot_Label_Mapping.csv/html` 查看原始指标名称。")
    lines.append("")
    return "\n".join(lines)


def _generate_output_bundle_base(results: Dict[str, Any], package_name: str = "BioEntropy_Results") -> Dict[str, Any]:
    temp_root = Path(tempfile.mkdtemp(prefix="bioentropy_run_"))
    output_dir = temp_root / package_name
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    tables = dict(results["tables"])
    label_map_df = build_plot_label_mapping_from_tables(tables)
    if label_map_df is not None and not label_map_df.empty:
        tables["图标签映射"] = label_map_df

    tables_dir = output_dir / "数据表"
    save_curated_tables_separately(tables, tables_dir)
    xlsx_path = tables_dir

    figure_paths: Dict[str, str] = {}
    for exp, exp_res in results["experiments"].items():
        summary_df = exp_res["summary"]
        outlier_df = exp_res["outliers"]
        comparison_df = exp_res["comparisons"]
        feature_cmp_df = exp_res.get("feature_comparisons", pd.DataFrame())
        treatment_group = exp_res.get("treatment_group")
        comparator_group = exp_res.get("comparator_group")

        p1 = figures_dir / f"{exp}_entropy_trajectory.png"
        p2 = figures_dir / f"{exp}_relative_burden_trajectory.png"
        p3 = figures_dir / f"{exp}_bri_entropy_scatter.png"
        p4 = figures_dir / f"{exp}_consensus_outlier_counts.png"
        p5 = figures_dir / f"{exp}_feature_repair_heatmap.png"
        p6 = figures_dir / f"{exp}_variance_stabilization_heatmap.png"

        plot_entropy_trajectory(summary_df, exp, p1, comparison_df, treatment_group, comparator_group)
        plot_relative_burden_trajectory(summary_df, exp, p2, comparison_df, treatment_group, comparator_group)
        plot_health_entropy_scatter(summary_df, exp, p3)
        plot_outlier_counts(outlier_df, exp, p4)
        plot_feature_repair_heatmap(feature_cmp_df, exp, p5)
        plot_variance_stabilization_heatmap(feature_cmp_df, exp, p6)

        figure_paths[f"{exp}_entropy"] = str(p1)
        figure_paths[f"{exp}_relative_burden"] = str(p2)
        figure_paths[f"{exp}_scatter"] = str(p3)
        figure_paths[f"{exp}_outliers"] = str(p4)
        figure_paths[f"{exp}_feature_repair"] = str(p5)
        figure_paths[f"{exp}_variance_stabilization"] = str(p6)

        if exp_res.get("reference_group") not in [None, NORMAL_REF_NONE]:
            p7 = figures_dir / f"{exp}_normal_proximity_trajectory.png"
            plot_normal_reference_trajectory(summary_df, exp, p7, comparison_df, treatment_group, comparator_group)
            if p7.exists():
                figure_paths[f"{exp}_normal_proximity"] = str(p7)

    zip_path = temp_root / f"{package_name}.zip"
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for path in output_dir.rglob("*"):
            if path.is_file():
                zf.write(path, arcname=str(path.relative_to(output_dir.parent)))

    return {
        "temp_root": str(temp_root),
        "output_dir": str(output_dir),
        "xlsx_path": str(xlsx_path),
        "report_path": "",
        "figure_paths": figure_paths,
        "zip_path": str(zip_path),
    }


# --------------------------------------------------------------------------- #
# External Normal-range support: NRBS / NRPS on top of the base pipeline
# --------------------------------------------------------------------------- #
import re as _re


EXTERNAL_NORMAL_FILE_ALIASES = {
    "normal", "norma", "normalrange", "normal_range",
    "healthy", "healthyrange", "healthy_range",
    "reference", "referencerange", "reference_range",
    "ref", "refrange", "normalref", "normalreference"
}


@dataclass
class ExternalNormalRangeParam:
    feature: str
    lower: float
    upper: float
    source_feature: str
    source_file: str
    unit_factor: float = 1.0
    bound_type: str = "two_sided"


def _canon_text_v5(x: Any) -> str:
    text = "" if pd.isna(x) else str(x).strip().lower()
    return _re.sub(r"[^0-9a-zA-Z\u4e00-\u9fff]+", "", text)


def _parse_external_normal_range_df(df: pd.DataFrame, file_name: str) -> pd.DataFrame:
    """Parse an external Normal-range reference file.

    Supported formats:
    1) Wide 2-row format:
       row labels in column 1 contain Normal_min / Normal_max, columns 3+ are features.
    2) Long format:
       columns contain Feature, Normal_min, Normal_max.
    """
    df = df.copy()
    if df.empty or df.shape[1] < 2:
        raise ValueError(f"{file_name} 为空，或列数不足。")

    # Wide format like the user's uploaded 106-Norma.xlsx
    first_col = df.iloc[:, 0].astype(str).str.strip()
    first_col_norm = first_col.map(_canon_text_v5)
    has_min_row = first_col_norm.str.contains("normalmin|min|下限").any()
    has_max_row = first_col_norm.str.contains("normalmax|max|上限").any()

    if has_min_row or has_max_row:
        min_idx = first_col_norm[first_col_norm.str.contains("normalmin|min|下限")].index
        max_idx = first_col_norm[first_col_norm.str.contains("normalmax|max|上限")].index
        min_idx = int(min_idx[0]) if len(min_idx) > 0 else None
        max_idx = int(max_idx[0]) if len(max_idx) > 0 else None

        rows = []
        for col in df.columns[2:]:
            lo = pd.to_numeric(pd.Series([df.loc[min_idx, col]]) if min_idx is not None else pd.Series([np.nan]), errors="coerce").iloc[0]
            hi = pd.to_numeric(pd.Series([df.loc[max_idx, col]]) if max_idx is not None else pd.Series([np.nan]), errors="coerce").iloc[0]
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

    # Long format
    col_map = {_canon_text_v5(c): str(c) for c in df.columns}
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
                "Normal_min": pd.to_numeric(df[min_col], errors="coerce") if min_col is not None else np.nan,
                "Normal_max": pd.to_numeric(df[max_col], errors="coerce") if max_col is not None else np.nan,
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


# Built-in published adult clinical-lab reference ranges. Only the pathological
# bound is set for each marker so the burden penalty (NRBS) is direction-aware:
# for higher-worse labs we set an upper limit and leave the lower open; for
# HDL-C (higher is better) we set only a lower limit. Values are standard adult
# reference limits (e.g. CDS/NCEP/ATP-III/lab norms) and provide an ABSOLUTE
# health anchor when a study has no enrolled Normal/healthy arm. This is a
# transparent, published-reference anchor — it is not tuned to any ranking.
#   "tokens":  core marker names matched as substrings of the study column name
#   "lower"/"upper":  reference bound in the study's native unit (NaN = open side)
PUBLISHED_LAB_REFERENCE_RANGES: List[Dict[str, Any]] = [
    {"tokens": ["空腹血浆血糖", "空腹血糖", "fastingplasmaglucose", "fpg"], "lower": np.nan, "upper": 6.1},
    {"tokens": ["血糖", "glucose"], "lower": np.nan, "upper": 6.1},
    {"tokens": ["糖化血红蛋白", "hba1c"], "lower": np.nan, "upper": 6.0},
    {"tokens": ["甘油三酯", "triglyceride", "tg"], "lower": np.nan, "upper": 1.7},
    {"tokens": ["非高密度脂蛋白胆固醇", "nonhdl", "非高密度"], "lower": np.nan, "upper": 4.1},
    {"tokens": ["低密度脂蛋白胆固醇", "ldl"], "lower": np.nan, "upper": 3.4},
    {"tokens": ["高密度脂蛋白胆固醇", "hdl"], "lower": 1.0, "upper": np.nan},
    {"tokens": ["总胆固醇", "totalcholesterol"], "lower": np.nan, "upper": 5.2},
    {"tokens": ["丙氨酸氨基转移酶", "alt", "谷丙转氨酶"], "lower": np.nan, "upper": 40.0},
    {"tokens": ["天门冬氨酸氨基转移酶", "ast", "谷草转氨酶"], "lower": np.nan, "upper": 40.0},
    {"tokens": ["谷氨酰转肽酶", "ggt", "gammagt"], "lower": np.nan, "upper": 60.0},
    {"tokens": ["超敏c反应蛋白", "c反应蛋白", "hscrp", "crp"], "lower": np.nan, "upper": 3.0},
]


def build_default_reference_range_df(
    raw: Dict[Tuple[str, int], pd.DataFrame],
    experiment: str,
) -> pd.DataFrame:
    """Build a Normal-range reference table from built-in published lab limits.

    Each study feature (column name intersected across all weeks) is matched to
    the most specific published reference marker whose canonical token appears in
    the canonicalized feature name. The longest matching token wins, so e.g.
    「非高密度脂蛋白胆固醇」 binds to the non-HDL limit rather than the HDL limit.
    Returns an empty frame when no feature matches (e.g. animal panels), leaving
    the pipeline to fall back to the baseline-relative BRI.
    """
    common = common_features(raw, experiment)
    rows = []
    for feat in common:
        canon = _canon_text_v5(feat)
        best_tok = None
        best_entry = None
        for entry in PUBLISHED_LAB_REFERENCE_RANGES:
            for tok in entry["tokens"]:
                tk = _canon_text_v5(tok)
                if tk and tk in canon and (best_tok is None or len(tk) > len(best_tok)):
                    best_tok = tk
                    best_entry = entry
        if best_entry is None:
            continue
        rows.append(
            {
                "Feature": feat,
                "Normal_min": best_entry["lower"],
                "Normal_max": best_entry["upper"],
                "Source_file": "built_in_published_reference_ranges",
                "Format": "built_in_published",
            }
        )
    return pd.DataFrame(rows)


def load_input_files(files: Iterable[Any]) -> Tuple[Dict[Tuple[str, int], pd.DataFrame], pd.DataFrame]:
    # load_input_bundle returns a 7-tuple:
    #   raw, range_refs, normal_refs, mixed_overview, data_overview, range_overview, normal_overview
    raw, _range_refs, _normal_refs, mixed_overview, _data_overview, _range_overview, _normal_overview = load_input_bundle(files)
    return raw, mixed_overview


def _harmonize_external_range_bounds(feature: str, lower: float, upper: float, study_series: pd.Series) -> Tuple[float, float, float, str]:
    """Harmonize obvious unit mismatches between uploaded range file and study data."""
    factor = 1.0
    note = ""
    name = str(feature).lower()
    med = float(pd.to_numeric(study_series, errors="coerce").dropna().median()) if study_series is not None and pd.to_numeric(study_series, errors="coerce").dropna().shape[0] > 0 else np.nan

    finite_bounds = [float(v) for v in [lower, upper] if pd.notna(v)]
    if finite_bounds:
        max_abs = max(abs(v) for v in finite_bounds)
    else:
        max_abs = np.nan

    if finite_bounds and max_abs <= 1 and pd.notna(med) and med > 2:
        if ("糖化血红蛋白" in name) or ("hba1c" in name):
            factor = 100.0
            note = "auto_scaled_x100_hba1c_fraction_to_percent"

    if pd.notna(lower):
        lower = float(lower) * factor
    if pd.notna(upper):
        upper = float(upper) * factor

    return lower, upper, factor, note


def fit_external_normal_range_params(
    raw: Dict[Tuple[str, int], pd.DataFrame],
    experiment: str,
    range_df: pd.DataFrame,
) -> Tuple[Dict[str, ExternalNormalRangeParam], pd.DataFrame]:
    """Fit external Normal-range parameters using only features common across all weeks."""
    common = common_features(raw, experiment)
    pooled = pd.concat([raw[(experiment, w)] for w in get_weeks(raw, experiment)], ignore_index=True)
    exp_map = {_canon_text_v5(f): f for f in common}

    rows = []
    params: Dict[str, ExternalNormalRangeParam] = {}
    for _, row in range_df.iterrows():
        src_feat = str(row["Feature"]).strip()
        key = _canon_text_v5(src_feat)
        if key not in exp_map:
            rows.append(
                {
                    "Experiment": experiment,
                    "Feature": src_feat,
                    "Matched_feature": np.nan,
                    "Used_in_primary_NRBS": False,
                    "Normal_min": row.get("Normal_min", np.nan),
                    "Normal_max": row.get("Normal_max", np.nan),
                    "Source_file": row.get("Source_file", ""),
                    "Unit_factor_applied": 1.0,
                    "Unit_note": "not_matched_to_common_features",
                }
            )
            continue

        feat = exp_map[key]
        lower = pd.to_numeric(pd.Series([row.get("Normal_min", np.nan)]), errors="coerce").iloc[0]
        upper = pd.to_numeric(pd.Series([row.get("Normal_max", np.nan)]), errors="coerce").iloc[0]
        lower, upper, factor, note = _harmonize_external_range_bounds(feat, lower, upper, pooled[feat])

        if pd.isna(lower) and pd.isna(upper):
            continue
        if pd.notna(lower) and pd.notna(upper) and lower > upper:
            lower, upper = upper, lower
            note = (note + "; swapped_min_max").strip("; ")

        if pd.notna(lower) and pd.notna(upper):
            bound_type = "two_sided"
        elif pd.notna(upper):
            bound_type = "upper_only"
        else:
            bound_type = "lower_only"

        params[feat] = ExternalNormalRangeParam(
            feature=feat,
            lower=lower,
            upper=upper,
            source_feature=src_feat,
            source_file=str(row.get("Source_file", "")),
            unit_factor=factor,
            bound_type=bound_type,
        )
        rows.append(
            {
                "Experiment": experiment,
                "Feature": feat,
                "Matched_feature": feat,
                "Used_in_primary_NRBS": True,
                "Normal_min": lower,
                "Normal_max": upper,
                "Source_feature": src_feat,
                "Source_file": row.get("Source_file", ""),
                "Unit_factor_applied": factor,
                "Unit_note": note,
                "Bound_type": bound_type,
            }
        )

    params_df = pd.DataFrame(rows)
    if not params_df.empty:
        params_df = params_df.sort_values(["Used_in_primary_NRBS", "Matched_feature", "Feature"], ascending=[False, True, True]).reset_index(drop=True)
    return params, params_df


def _external_range_penalty_series(values: pd.Series, lower: float, upper: float) -> pd.Series:
    s = pd.to_numeric(values, errors="coerce").astype(float).copy()
    out = pd.Series(np.nan, index=s.index, dtype=float)
    mask = s.notna()
    if mask.sum() == 0:
        return out

    x = s.loc[mask].astype(float).to_numpy()
    if pd.notna(lower) and pd.notna(upper):
        scale = max(float(upper - lower), abs(float(upper)) * 0.1, abs(float(lower)) * 0.1, 1e-9)
        excess = np.where(
            x < lower,
            (lower - x) / scale,
            np.where(x > upper, (x - upper) / scale, 0.0),
        )
    elif pd.notna(upper):
        scale = max(abs(float(upper)), 1e-9)
        excess = np.maximum(0.0, x - upper) / scale
    elif pd.notna(lower):
        scale = max(abs(float(lower)), 1e-9)
        excess = np.maximum(0.0, lower - x) / scale
    else:
        excess = np.full_like(x, np.nan, dtype=float)

    penalty = np.log1p(excess)
    out.loc[mask] = penalty
    return out


def compute_external_normal_range_scores(
    df: pd.DataFrame,
    range_params: Dict[str, ExternalNormalRangeParam],
) -> Tuple[pd.DataFrame, pd.Series, pd.Series, pd.Series]:
    """Compute external Normal-range burden / proximity scores.

    NRBS = mean(log1p(relative excess outside uploaded normal interval))
    NRPS = 1 / (1 + NRBS)
    """
    penalty_df = pd.DataFrame(index=df.index)
    for feat, p in range_params.items():
        if feat not in df.columns:
            continue
        penalty_df[feat] = _external_range_penalty_series(df[feat], p.lower, p.upper)

    if penalty_df.empty:
        nrbs = pd.Series(np.nan, index=df.index, name="NRBS")
        nrps = pd.Series(np.nan, index=df.index, name="NRPS")
        coverage = pd.Series(0, index=df.index, name="NR_feature_count")
        return penalty_df, nrbs, nrps, coverage

    nrbs = penalty_df.mean(axis=1, skipna=True).rename("NRBS")
    nrps = (1.0 / (1.0 + nrbs)).rename("NRPS")
    coverage = penalty_df.notna().sum(axis=1).astype(int).rename("NR_feature_count")
    nrps = nrps.where(coverage > 0, np.nan)
    nrbs = nrbs.where(coverage > 0, np.nan)
    return penalty_df, nrbs, nrps, coverage


def _external_range_feature_table(
    penalty_df: pd.DataFrame,
    groups: pd.Series,
    experiment: str,
    week: int,
    range_params: Dict[str, ExternalNormalRangeParam],
) -> pd.DataFrame:
    rows = []
    if penalty_df.empty:
        return pd.DataFrame()
    tmp = penalty_df.copy()
    tmp["Group"] = groups.values
    for group, sub in tmp.groupby("Group"):
        for feat in penalty_df.columns:
            s = pd.to_numeric(sub[feat], errors="coerce").dropna().astype(float)
            if s.empty:
                continue
            if len(s) >= 2:
                var_log = float(0.5 * np.log(max(s.var(ddof=1), 1e-8)))
            else:
                var_log = float("nan")
            p = range_params.get(feat)
            rows.append(
                {
                    "Experiment": experiment,
                    "Week": week,
                    "Group": group,
                    "Feature": feat,
                    "Signed_burden_mean": float(s.mean()),
                    "Mean_abs_z": float(s.mean()),
                    "Variance_log_component": var_log,
                    "Pct_in_normal_range": float((s <= 1e-12).mean() * 100.0),
                    "Normal_min": p.lower if p is not None else np.nan,
                    "Normal_max": p.upper if p is not None else np.nan,
                    "Bound_type": p.bound_type if p is not None else "",
                }
            )
    return pd.DataFrame(rows)


def _range_scatter(summary_df: pd.DataFrame, experiment: str, out_path: Path, metric_col: str, xlabel: str, title: str) -> None:
    sub = summary_df[(summary_df["Experiment"] == experiment) & (summary_df[metric_col].notna())].copy()
    if sub.empty:
        return
    groups = list(dict.fromkeys(sub["Group"].tolist()))
    cmap = _group_color_map(groups)
    fig, ax = plt.subplots(figsize=(7.2, 5.8))
    for idx, group in enumerate(groups):
        g = sub[sub["Group"] == group].sort_values("Week")
        label = display_group_name(group, idx)
        ax.scatter(g[metric_col], g["Entropy"], s=72, alpha=0.92, color=cmap[group], label=label)
        for _, row in g.iterrows():
            ax.annotate(f"{int(row['Week'])}W", (row[metric_col], row["Entropy"]), fontsize=8, xytext=(4, 3), textcoords="offset points")
    if len(sub) >= 3:
        rho, pval = spearmanr(sub[metric_col], sub["Entropy"])
        ax.text(
            0.02, 0.98, f"Spearman rho = {rho:.3f}\np = {pval:.3g}",
            transform=ax.transAxes, ha="left", va="top", fontsize=9.2,
            bbox=dict(boxstyle="round,pad=0.24", facecolor="white", edgecolor="#dddddd", alpha=0.9)
        )
    ax.set_title(f"Experiment {experiment}: {title}")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Entropy")
    _format_axis(ax)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_normal_range_burden_trajectory(
    summary_df: pd.DataFrame,
    experiment: str,
    out_path: Path,
    comparison_df: pd.DataFrame | None = None,
    treatment_group: str | None = None,
    comparator_group: str | None = None,
) -> None:
    sub = summary_df[(summary_df["Experiment"] == experiment) & (summary_df["NRBS_mean"].notna())].copy()
    if sub.empty:
        return
    _draw_publication_trajectory(
        summary_df=summary_df,
        comparison_df=comparison_df,
        experiment=experiment,
        metric="NRBS_mean",
        treatment_group=treatment_group,
        comparator_group=comparator_group,
        out_path=out_path,
    )


def plot_normal_range_burden_entropy_scatter(summary_df: pd.DataFrame, experiment: str, out_path: Path) -> None:
    _range_scatter(
        summary_df=summary_df,
        experiment=experiment,
        out_path=out_path,
        metric_col="NRBS_mean",
        xlabel="Normal-range burden score (NRBS)",
        title="NRBS vs Entropy",
    )


def plot_normal_range_feature_repair_heatmap(
    feature_cmp_df: pd.DataFrame,
    experiment: str,
    out_path: Path,
    top_n: int = 12,
) -> None:
    sub = feature_cmp_df[feature_cmp_df["Experiment"] == experiment].copy()
    if sub.empty:
        return
    label_map, _map_df = _feature_label_map_for_sub(feature_cmp_df, experiment)
    week_max = sub["Week"].max()
    rank = (
        sub[sub["Week"] == week_max]
        .assign(abs_rank=lambda d: d["Burden_diff_treat_minus_comp"].abs())
        .sort_values("abs_rank", ascending=False)
    )
    features = rank["Feature"].head(top_n).tolist()
    if not features:
        return
    pivot = (
        sub[sub["Feature"].isin(features)]
        .pivot_table(index="Feature", columns="Week", values="Burden_diff_treat_minus_comp", aggfunc="mean")
        .reindex(features)
    )
    fig, ax = plt.subplots(figsize=(9.4, max(4.8, 0.46 * len(features) + 2.2)))
    mat = pivot.to_numpy(dtype=float)
    vmax = np.nanmax(np.abs(mat)) if np.isfinite(mat).any() else 1.0
    vmax = max(vmax, 0.15)
    im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ylabels = [label_map.get(feat, feat) for feat in features]
    ax.set_yticks(range(len(features)))
    ax.set_yticklabels(ylabels, fontsize=8.6)
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels([f"{int(x)}W" for x in pivot.columns], fontsize=9)
    ax.set_title(
        f"Experiment {experiment}: normal-range repair heatmap\n(treatment - comparator; lower = closer to uploaded normal interval in treatment)",
        fontsize=12.5, pad=12
    )
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if np.isfinite(mat[i, j]):
                ax.text(j, i, f"{mat[i, j]:+.2f}", ha="center", va="center", fontsize=7.2, color="black")
    cbar = fig.colorbar(im, ax=ax, shrink=0.9)
    cbar.set_label("Signed NRBS difference", fontsize=9)
    if not HAS_CJK_FONT:
        ax.text(
            0.01, -0.14,
            "CJK font unavailable: exported PNG uses English aliases or feature codes. See Plot_Label_Mapping sheet / CSV.",
            transform=ax.transAxes, ha="left", va="top", fontsize=8.2
        )
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def build_normal_range_template(out_path: Path) -> Path:
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "NormalRange"

    headers = ["Group", "Sample ID", "空腹血浆血糖 (mmol/L)", "低密度脂蛋白胆固醇（LDL-C, mmol/L）", "高密度脂蛋白胆固醇（HDL-C, mmol/L）", "糖化血红蛋白", "丙氨酸氨基转移酶（ALT, U/L）"]
    ws.append(headers)
    ws.append(["Normal_min", 1, 3.9, None, 1.03, None, None])
    ws.append(["Normal_max", 2, 5.6, 2.6, None, 5.7, 55])

    fill = PatternFill("solid", fgColor="1F4E78")
    font = Font(color="FFFFFF", bold=True)
    for cell in ws[1]:
        cell.fill = fill
        cell.font = font
    for idx, width in enumerate([16, 12, 22, 26, 26, 16, 20], start=1):
        ws.column_dimensions[get_column_letter(idx)].width = width

    note = wb.create_sheet("Instructions")
    note["A1"] = "BioEntropy 外部 Normal 范围模板说明"
    note["A1"].font = Font(bold=True, size=14)
    note["A3"] = "1. 文件命名请使用：实验号-Normal.xlsx，例如 106-Normal.xlsx。"
    note["A4"] = "2. 第一行表头，第1列为 Group，第2列为 Sample ID（仅占位，可忽略），第3列起为指标名。"
    note["A5"] = "3. 第二行 Group=Normal_min：填写正常下限；没有下限可留空。"
    note["A6"] = "4. 第三行 Group=Normal_max：填写正常上限；没有上限可留空。"
    note["A7"] = "5. 本文件用于构建区间型疾病负担分数（NRBS/NRPS），不是实验样本组。"
    note["A8"] = "6. 真正的 Normal/健康样本应作为一个真实分组放在实验时间点文件中，不要填写在本模板内。"
    note.column_dimensions["A"].width = 110

    wb.save(out_path)
    return out_path


# ======================= V6 TRUE NORMAL SAMPLE FILE EXTENSION =======================

TRUE_NORMAL_FILE_ALIASES = {
    "normal", "normalsample", "normalsamples", "normalgroup", "normal_sample",
    "healthy", "healthysample", "healthysamples", "healthygroup", "healthy_sample",
    "referencenormal", "referencesample", "truernormal", "truenormal", "healthreference"
}
UPLOADED_TRUE_NORMAL_LABEL = "（使用上传的Normal样本文件）"


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


def _strip_copy_suffix_v6(text: str) -> str:
    return _re.sub(r"\(\d+\)$", "", str(text).strip())


def _simplify_text_v6(text: Any) -> str:
    x = "" if pd.isna(text) else str(text)
    x = x.replace("（", "(").replace("）", ")").replace("，", ",").replace("、", "")
    x = x.replace(" ", "").replace("-", "").replace("_", "")
    x = x.lower()
    x = x.replace("γ", "gamma").replace("α", "alpha").replace("β", "beta")
    x = _re.sub(r"[^0-9a-zA-Z\u4e00-\u9fff]+", "", x)
    return x


def feature_semantic_code(feature_name: Any) -> Optional[str]:
    s = _simplify_text_v6(feature_name)

    mapping_rules = [
        (["非高密度脂蛋白胆固醇", "nonhdl"], "non_hdl"),
        (["极低密度脂蛋白胆固醇", "vldl"], "vldl"),
        (["高密度脂蛋白胆固醇", "hdlc", "hdl"], "hdl"),
        (["低密度脂蛋白胆固醇", "ldlc", "ldl"], "ldl"),
        (["空腹血浆血糖", "空腹血糖", "fastingglucose", "fpg"], "glu_fasting"),
        (["餐后2h血浆血糖", "餐后2h血糖", "餐前餐后血糖", "postprandialglucose", "ppg", "2hglucose"], "glu_postprandial"),
        (["糖化血红蛋白", "糖化", "hba1c"], "hba1c"),
        # Bare blood glucose (no fasting/postprandial/HbA1c qualifier). Placed
        # after the qualified-glucose and HbA1c rules so it only catches a plain
        # 血糖 / glucose column.
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
        # A bare Normal/健康 workbook (e.g. Normal.xlsx) carries no experiment
        # prefix; treat it as a GLOBAL Normal reference applied to every
        # experiment in the batch (experiment=None). Content still decides
        # sample-vs-range downstream.
        canon_stem = _canon_text_v5(stem)
        if canon_stem in {_canon_text_v5(x) for x in TRUE_NORMAL_FILE_ALIASES} or \
           canon_stem in {_canon_text_v5(x) for x in EXTERNAL_NORMAL_FILE_ALIASES}:
            return "normal_like", None, None
        # Animal studies are commonly delivered as one wide workbook such as
        # DB鼠.xlsx / KK鼠.xlsx.  Treat such files as a single 0W experiment;
        # the sheet layout validation below will still reject non-data files.
        return "week", stem, 0
    experiment, suffix = stem.rsplit("-", 1)
    experiment = experiment.strip()
    suffix = _strip_copy_suffix_v6(suffix.strip())

    m = _re.fullmatch(r"(\d+)\s*[Ww]", suffix)
    if m:
        return "week", experiment, int(m.group(1))

    canon = _canon_text_v5(suffix)
    if canon in {_canon_text_v5(x) for x in EXTERNAL_NORMAL_FILE_ALIASES}:
        return "normal_like", experiment, None
    if canon in {_canon_text_v5(x) for x in TRUE_NORMAL_FILE_ALIASES}:
        return "normal_like", experiment, None

    # Fallback for animal/non-longitudinal files whose names contain hyphens
    # but do not encode a week, e.g. 动物实验-DB-DB鼠.xlsx.
    return "week", stem, 0


def parse_uploaded_filename(file_name: str) -> Tuple[str, int]:
    kind, experiment, week = parse_input_file_kind(file_name)
    if kind != "week" or week is None:
        raise ValueError(f"文件名 {file_name} 不是实验时间点文件。请使用“实验号-周数W.xlsx”。")
    return experiment, week


def _looks_like_external_normal_range_df(df: pd.DataFrame) -> bool:
    if df is None or df.empty:
        return False
    first_col = df.iloc[:, 0].astype(str).str.strip()
    first_col_norm = first_col.map(_canon_text_v5)
    if first_col_norm.str.contains("normalmin|min|下限").any() or first_col_norm.str.contains("normalmax|max|上限").any():
        return True
    col_norm = {_canon_text_v5(c): str(c) for c in df.columns}
    has_feature_col = any(k in {"feature", "biomarker", "indicator", "指标", "特征"} for k in col_norm)
    has_min_col = any(k in {"normalmin", "min", "下限", "lower", "low"} for k in col_norm)
    has_max_col = any(k in {"normalmax", "max", "上限", "upper", "high"} for k in col_norm)
    return has_feature_col and (has_min_col or has_max_col)


def _parse_true_normal_sample_df(df: pd.DataFrame, file_name: str) -> pd.DataFrame:
    if df.shape[1] < 3:
        raise ValueError(f"{file_name} 至少需要 3 列：Group、Sample ID、指标列。")
    out = df.copy()
    # A Normal-sample workbook often carries just an ID column followed by marker
    # columns (no Group column), e.g. 健康人ID | 血糖 | 总胆固醇 | ...  In that case
    # the whole file is the Normal cloud: synthesise Group="Normal" and keep the
    # ID column as Sample ID so no marker column is consumed as a group/ID.
    first_header = _canon_text_v5(out.columns[0])
    id_tokens = ("id", "编号", "序号", "编码", "健康人", "subject", "subjectid", "sampleid", "sample")
    has_explicit_group = _canon_text_v5(out.columns[1]) in {"sampleid", "sample", "id", "编号", "序号"}
    if (any(tok in first_header for tok in id_tokens)) and not has_explicit_group:
        out.insert(0, "Group", "Normal")
        out = out.rename(columns={out.columns[1]: "Sample ID"})
    else:
        out = out.rename(columns={out.columns[0]: "Group", out.columns[1]: "Sample ID"})
    out["Group"] = out["Group"].ffill().map(normalize_group_name)
    if out["Group"].isna().all():
        out["Group"] = "Normal"
    out["Group"] = out["Group"].fillna("Normal")
    out["Sample ID"] = pd.to_numeric(out["Sample ID"], errors="coerce")
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


def _validate_week_dataframe(name: str, df: pd.DataFrame) -> None:
    """Fail fast with a specific, actionable message on malformed experiment files.

    Runs before the (slow) analysis so users see exactly which file/column/ID is
    wrong instead of a mid-run traceback. Only raises on genuinely broken input;
    valid numeric data passes untouched.
    """
    # 1. Sample IDs (when numeric) must be unique within a file. Non-numeric IDs
    #    are coerced to NaN upstream and tolerated (e.g. clinical patient codes),
    #    so we only check duplicates among the usable numeric IDs.
    ids = df["Sample ID"].dropna()
    dups = ids[ids.duplicated()].unique()
    if len(dups) > 0:
        shown = ", ".join(str(int(d)) if float(d).is_integer() else str(d) for d in dups[:8])
        raise ValueError(
            f"{name}: 检测到重复的 Sample ID（{shown}{' …' if len(dups) > 8 else ''}）。"
            f"同一文件内每个样本的 ID 必须唯一。"
        )

    # 3. Every feature column that has data must contain at least one number;
    #    an all-text column usually means a shifted header or the wrong column.
    feature_cols = list(df.columns[2:])
    if not feature_cols:
        raise ValueError(f"{name}: 没有任何指标列（至少需要 Group、Sample ID 和 1 个数值指标列）。")
    bad_cols = []
    for col in feature_cols:
        s = df[col]
        non_empty = s[s.notna() & (s.astype(str).str.strip() != "")]
        if len(non_empty) == 0:
            continue  # entirely empty column is tolerated
        if pd.to_numeric(non_empty, errors="coerce").notna().sum() == 0:
            bad_cols.append(str(col))
    if bad_cols:
        raise ValueError(
            f"{name}: 以下指标列没有任何可识别的数值，请检查是否填入了文本、单位或选错了列："
            f"{', '.join(bad_cols[:6])}{' …' if len(bad_cols) > 6 else ''}。"
        )


def load_input_bundle(files: Iterable[Any]) -> Tuple[
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
            df["Sample ID"] = pd.to_numeric(df["Sample ID"], errors="coerce")
            _validate_week_dataframe(name, df)
            groups = sorted(df["Group"].dropna().astype(str).unique().tolist())

            if (exp, week) in raw:
                raise ValueError(f"检测到重复文件：实验 {exp} 的 {week}W 出现了两次。")
            raw[(exp, week)] = df

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

        # normal-like files are auto-detected by content
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
                    "格式": parsed["Format"].iloc[0] if "Format" in parsed.columns and len(parsed) > 0 else "",
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
                    "分组": " | ".join(sorted(parsed["Group"].dropna().astype(str).unique().tolist())),
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

    # Distribute global (prefix-less) Normal references to every experiment that
    # did not receive its own experiment-specific Normal file.
    experiments_present = sorted({exp for (exp, _week) in raw.keys()})
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

    data_overview = pd.DataFrame(data_rows).sort_values(["实验号", "周数", "文件名"]).reset_index(drop=True)
    range_overview = pd.DataFrame(range_rows)
    if not range_overview.empty:
        range_overview = range_overview.sort_values(["实验号", "文件名"]).reset_index(drop=True)
    normal_overview = pd.DataFrame(normal_rows)
    if not normal_overview.empty:
        normal_overview = normal_overview.sort_values(["实验号", "文件名"]).reset_index(drop=True)
    mixed_overview = pd.DataFrame(mixed_rows).sort_values(["实验号", "类型", "文件名"]).reset_index(drop=True)
    return raw, range_refs, normal_refs, mixed_overview, data_overview, range_overview, normal_overview


def _semantic_match_true_normal_features(
    experiment_features: List[str],
    normal_df: pd.DataFrame,
) -> Dict[str, str]:
    normal_cols = [c for c in normal_df.columns[2:] if c != "__source_file__"]
    normal_by_canon: Dict[str, List[str]] = {}
    for col in normal_cols:
        normal_by_canon.setdefault(_simplify_text_v6(col), []).append(col)

    mapping: Dict[str, str] = {}
    for feat in experiment_features:
        exact_candidates = normal_by_canon.get(_simplify_text_v6(feat), [])
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
            # Prefer same units / same keywords when multiple candidates exist
            if len(candidates) == 1:
                mapping[feat] = candidates[0]
            else:
                best = sorted(candidates, key=lambda c: (0 if _simplify_text_v6(c) == _simplify_text_v6(feat) else 1, len(str(c))))[0]
                mapping[feat] = best
    return mapping


def _harmonize_true_normal_units(feature: str, normal_series: pd.Series, study_series: pd.Series) -> Tuple[pd.Series, float, str]:
    s_norm = pd.to_numeric(normal_series, errors="coerce").astype(float).copy()
    s_study = pd.to_numeric(study_series, errors="coerce").astype(float)
    factor = 1.0
    note = ""

    name = str(feature).lower()
    med_norm = float(s_norm.dropna().median()) if s_norm.dropna().shape[0] > 0 else np.nan
    med_study = float(s_study.dropna().median()) if s_study.dropna().shape[0] > 0 else np.nan

    # Conventional (US, mg/dL) → SI (mmol/L) conversion constants for the analytes
    # a Normal-sample file commonly ships in US units while the study is in SI.
    code = feature_semantic_code(feature)
    conv = None
    if code in {"glu_fasting", "glu_postprandial"} or "血糖" in str(feature) or "glucose" in name:
        conv = 0.0555084   # glucose mg/dL → mmol/L
    elif code in {"tc", "hdl", "ldl", "non_hdl", "vldl"}:
        conv = 0.0258600   # cholesterol mg/dL → mmol/L
    elif code == "tg":
        conv = 0.0112990   # triglyceride mg/dL → mmol/L

    def _plausible(f: float) -> bool:
        # Same analyte in matched units: healthy vs diseased medians differ by at
        # most a few-fold, so the ratio should sit in a wide biological band.
        if not (np.isfinite(med_norm) and np.isfinite(med_study)) or med_study <= 0 or med_norm <= 0:
            return False
        r = (med_norm * f) / med_study
        return 0.2 <= r <= 5.0

    if pd.notna(med_norm) and pd.notna(med_study):
        if (("糖化" in name) or ("hba1c" in name)) and med_norm < 1 and med_study > 2:
            factor = 100.0
            note = "auto_scaled_x100_hba1c_fraction_to_percent"
        elif conv is not None and not _plausible(1.0):
            # Same-unit is implausible; adopt the unit conversion that brings the
            # Normal median into a biologically plausible ratio of the study median.
            if _plausible(conv):
                factor = conv
                note = "auto_converted_mgdl_to_mmoll"
            elif _plausible(1.0 / conv):
                factor = 1.0 / conv
                note = "auto_converted_mmoll_to_mgdl"
        elif factor == 1.0 and med_norm > 0 and med_study > 0:
            ratio = med_study / med_norm if med_norm != 0 else np.nan
            if ("肌酐" in name or "creatinine" in name) and ratio > 10 and ratio < 1000:
                # do not auto-rescale because creatinine units can genuinely differ; leave review note instead
                note = "review_creatinine_unit_if_needed"
    if factor != 1.0:
        s_norm = s_norm * factor
    return s_norm, factor, note


def fit_true_normal_reference(
    raw: Dict[Tuple[str, int], pd.DataFrame],
    experiment: str,
    normal_df: pd.DataFrame,
) -> Tuple[TrueNormalReference, pd.DataFrame]:
    exp_features = common_features(raw, experiment)
    pooled = pd.concat([raw[(experiment, w)] for w in get_weeks(raw, experiment)], ignore_index=True)

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
    source_files = sorted(normal_df["__source_file__"].dropna().astype(str).unique().tolist()) if "__source_file__" in normal_df.columns else []

    for feat, src_feat in mapping.items():
        study_s = pd.to_numeric(pooled[feat], errors="coerce").astype(float)
        norm_s_raw = pd.to_numeric(normal_df[src_feat], errors="coerce").astype(float)
        norm_s, unit_factor, unit_note = _harmonize_true_normal_units(feat, norm_s_raw, study_s)

        combo = pd.concat([study_s, norm_s], ignore_index=True)
        combo = combo.dropna().astype(float)
        use_log = bool((combo.shape[0] > 0) and (combo.min() > 0) and (abs(skew(combo)) > 1))
        if use_log:
            t_norm = np.log1p(norm_s)
        else:
            t_norm = norm_s.copy()

        valid = t_norm.dropna().astype(float)
        if valid.shape[0] < 3:
            continue

        median = float(np.median(valid))
        mad = float(np.median(np.abs(valid - median)))
        if mad < 1e-9:
            mad = float((valid.quantile(0.75) - valid.quantile(0.25)) / 1.349)
        if mad < 1e-9:
            std_tmp = float(valid.std(ddof=1))
            mad = std_tmp if std_tmp > 1e-9 else 1.0

        lo = float(median - 4.0 * mad)
        hi = float(median + 4.0 * mad)
        if not np.isfinite(lo) or not np.isfinite(hi) or lo >= hi:
            q1 = float(valid.quantile(0.25))
            q3 = float(valid.quantile(0.75))
            iqr = q3 - q1
            if iqr < 1e-9:
                iqr = mad * 1.349
            lo = float(q1 - 1.5 * iqr)
            hi = float(q3 + 1.5 * iqr)
        clip_bounds[feat] = (lo, hi)
        t_norm_clip = t_norm.clip(lo, hi)

        normal_t[feat] = t_norm_clip
        mean_map[feat] = float(pd.to_numeric(t_norm_clip, errors="coerce").dropna().mean())
        sd_val = float(pd.to_numeric(t_norm_clip, errors="coerce").dropna().std(ddof=1))
        sd_map[feat] = sd_val if np.isfinite(sd_val) and sd_val > 1e-9 else 1.0
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
                "Normal_raw_median_after_unit": float(norm_s.dropna().median()) if norm_s.dropna().shape[0] > 0 else np.nan,
                "Study_raw_median": float(study_s.dropna().median()) if study_s.dropna().shape[0] > 0 else np.nan,
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
    X_norm = fill_nan_with_col_median(normal_t[features].to_numpy(dtype=float))
    lw = LedoitWolf().fit(X_norm)
    mu = np.asarray(np.nanmean(X_norm, axis=0), dtype=float)
    precision = np.asarray(lw.precision_, dtype=float)
    covariance = np.asarray(lw.covariance_, dtype=float)

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
        n_normal=int(len(normal_df)),
        source_files=source_files,
        clip_bounds={k: clip_bounds[k] for k in features},
    )
    return ref, pd.DataFrame(rows).sort_values(["Feature"]).reset_index(drop=True)


def apply_true_normal_reference(
    df: pd.DataFrame,
    ref: TrueNormalReference,
) -> pd.DataFrame:
    out = pd.DataFrame(index=df.index)
    feature_count = pd.Series(0, index=df.index, dtype=float)

    cols = []
    for feat in ref.features:
        s = pd.to_numeric(df[feat], errors="coerce").astype(float)
        if ref.use_log_map.get(feat, False):
            s = np.log1p(s)
        cols.append(s)
        feature_count += s.notna().astype(float)

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
        [max(float(ref.per_feature_sd.get(feat, 1.0)), 1e-9) for feat in ref.features],
        dtype=float,
    )
    z = diff / sd.reshape(1, -1)
    z = np.where(np.isfinite(z), z, 0.0)
    z_capped = np.clip(z, -NMD_Z_CAP, NMD_Z_CAP)
    d2 = np.sum(z_capped * z_capped, axis=1)
    d2 = np.clip(d2, 0, None)
    md = np.sqrt(d2)
    # Diagonal chi-square p-value: assumes the per-feature scales are KNOWN and
    # the features independent-normal. It ignores that the scales/mean are
    # estimated from a finite Normal reference, so it is optimistic at small m.
    pval = 1.0 - chi2.cdf(d2, df=len(ref.features))

    # Small-sample-correct p-value for the FULL-covariance Mahalanobis distance.
    # For a new observation vs a Normal reference of m samples and p features,
    #   ((m - p) / (p (m - 1))) * (m / (m + 1)) * D^2  ~  F(p, m - p)
    # (out-of-sample Hotelling T^2). Defined only when m > p + 1.
    m_ref = int(getattr(ref, "n_normal", 0))
    p_ref = len(ref.features)
    if m_ref > p_ref + 1 and p_ref >= 1:
        f_stat = ((m_ref - p_ref) / (p_ref * (m_ref - 1))) * (m_ref / (m_ref + 1)) * d2_cov
        out["NMD_cov_p_smallsample_F"] = 1.0 - f_dist.cdf(f_stat, p_ref, m_ref - p_ref)
    else:
        out["NMD_cov_p_smallsample_F"] = np.full(len(df), np.nan)

    out["NMD"] = md
    out["NMD_sq"] = d2
    out["NMD_p_chi2"] = pval
    out["NMD_covariance_raw"] = md_cov
    out["NMD_sq_covariance_raw"] = d2_cov
    out["NMD_method"] = f"robust_diagonal_zcap_{NMD_Z_CAP:g}"
    out["NMD_z_cap"] = NMD_Z_CAP
    out["NMD_feature_count"] = feature_count.to_numpy(dtype=float)
    return out


def compute_normal_reference_nmd_band(
    ref: "TrueNormalReference",
    normal_df: pd.DataFrame,
) -> Optional[Dict[str, float]]:
    """NMD distribution of the Normal reference subjects against their own cloud.

    Applying the same robust-capped Mahalanobis transform to the Normal samples
    yields the distance-from-Normal that *healthy* individuals themselves carry.
    This is the reference band a treatment arm would have to reach to be
    indistinguishable from Normal. Returns per-subject central range + the group
    mean with a bootstrap CI, or None if it cannot be computed.
    """
    if ref is None or normal_df is None or normal_df.empty:
        return None
    nd = pd.DataFrame(index=normal_df.index)
    for feat in ref.features:
        src = ref.source_feature_map.get(feat, feat)
        col = src if src in normal_df.columns else (feat if feat in normal_df.columns else None)
        if col is None:
            return None
        nd[feat] = pd.to_numeric(normal_df[col], errors="coerce").astype(float) * float(ref.unit_factors.get(feat, 1.0))
    md = apply_true_normal_reference(nd, ref)
    vals = pd.to_numeric(md["NMD"], errors="coerce").dropna().to_numpy(dtype=float)
    if vals.shape[0] < 3:
        return None
    lo, hi, _se = bootstrap_mean_ci(vals, n_boot=500, seed=RNG_SEED + 4200)
    p2_5, p25, p50, p75, p97_5 = (float(x) for x in np.percentile(vals, [2.5, 25, 50, 75, 97.5]))
    return {
        "metric": "NMD",
        "n": int(vals.shape[0]),
        "mean": float(np.mean(vals)),
        "mean_ci_low": float(lo),
        "mean_ci_high": float(hi),
        "median": p50,
        "p2_5": p2_5,
        "p25": p25,
        "p75": p75,
        "p97_5": p97_5,
    }


def _true_normal_feature_table(
    df: pd.DataFrame,
    experiment: str,
    week: int,
    ref: TrueNormalReference,
) -> pd.DataFrame:
    rows = []
    for group, sub in df.groupby("Group"):
        for feat in ref.features:
            s = pd.to_numeric(sub[feat], errors="coerce").astype(float)
            if ref.use_log_map.get(feat, False):
                s_t = np.log1p(s)
            else:
                s_t = s.copy()
            group_mean_t = float(s_t.mean()) if s_t.notna().sum() > 0 else np.nan
            z_abs = abs((group_mean_t - ref.per_feature_mean[feat]) / ref.per_feature_sd[feat]) if pd.notna(group_mean_t) else np.nan
            rows.append(
                {
                    "Experiment": experiment,
                    "Week": week,
                    "Group": group,
                    "Feature": feat,
                    "Source_feature": ref.source_feature_map.get(feat, ""),
                    "Semantic_code": feature_semantic_code(feat),
                    "Group_mean_transformed": group_mean_t,
                    "Normal_mean_transformed": ref.per_feature_mean[feat],
                    "Normal_sd_transformed": ref.per_feature_sd[feat],
                    "Abs_z_to_normal": z_abs,
                    "Direction_for_health": "closer_to_normal_is_better",
                }
            )
    return pd.DataFrame(rows)


def plot_true_normal_md_trajectory(
    summary_df: pd.DataFrame,
    experiment: str,
    out_path: Path,
    comparison_df: pd.DataFrame | None = None,
    treatment_group: str | None = None,
    comparator_group: str | None = None,
) -> None:
    _draw_publication_trajectory(
        summary_df=summary_df,
        comparison_df=comparison_df,
        experiment=experiment,
        metric="NMD_mean",
        treatment_group=treatment_group,
        comparator_group=comparator_group,
        out_path=out_path,
    )


def plot_true_normal_feature_repair_heatmap(feature_cmp_df: pd.DataFrame, experiment: str, out_path: Path) -> None:
    sub = feature_cmp_df[feature_cmp_df["Experiment"] == experiment].copy()
    if sub.empty:
        return
    pivot = sub.pivot(index="Feature", columns="Week", values="Abs_z_to_normal_diff_treat_minus_comp")
    if pivot.empty:
        return
    pivot = pivot.sort_index()
    pivot = pivot.loc[pivot.abs().mean(axis=1).sort_values(ascending=False).index]
    if len(pivot) > 12:
        pivot = pivot.iloc[:12]

    labels = [feature_display_label(f, idx)[0] for idx, f in enumerate(pivot.index.tolist())]
    fig_w = max(7.5, 1.25 * pivot.shape[1] + 3.6)
    fig_h = max(5.0, 0.46 * pivot.shape[0] + 2.1)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    vmax = np.nanmax(np.abs(pivot.values)) if np.isfinite(pivot.values).any() else 1.0
    vmax = max(vmax, 0.2)
    im = ax.imshow(pivot.values, aspect="auto", cmap="coolwarm", vmin=-vmax, vmax=vmax)

    ax.set_xticks(range(pivot.shape[1]))
    ax.set_xticklabels([f"{int(x)}W" for x in pivot.columns], fontsize=10)
    ax.set_yticks(range(pivot.shape[0]))
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_title(
        f"Experiment {experiment}: deviation-from-Normal repair heatmap\n"
        "(treatment - comparator; lower = treatment closer to Normal)",
        fontsize=12,
        pad=14,
    )
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Difference in |z distance to Normal|", fontsize=10)

    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            val = pivot.values[i, j]
            if np.isfinite(val):
                ax.text(j, i, f"{val:+.2f}", ha="center", va="center", fontsize=8)

    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _metric_spec(metric: str) -> dict:
    if metric == "NMD_mean":
        return {
            "value": "NMD_mean",
            "low": "NMD_mean_ci_low",
            "high": "NMD_mean_ci_high",
            "p": "NMD_p_mannwhitney",
            "diff": "NMD_diff_treat_minus_comp",
            "ylabel": "True-Normal robust distance (NMD)",
            "title": "Disease-burden trajectory vs true Normal",
            "note": "NMD is the robust capped Normal-reference distance; raw covariance distance is retained as NMD_covariance_raw. Lower values mean closer to the Normal reference cloud.",
        }
    if metric == "NRBS_mean":
        return {
            "value": "NRBS_mean",
            "low": "NRBS_mean_ci_low",
            "high": "NRBS_mean_ci_high",
            "p": "NRBS_p_mannwhitney",
            "diff": "NRBS_diff_treat_minus_comp",
            "ylabel": "Normal-range burden score (NRBS)",
            "title": "Disease-burden trajectory vs published Normal ranges",
            "note": "NRBS is the mean log-excess outside published adult reference ranges (an absolute health anchor, not baseline-relative). 0 means every marker is within its Normal range; higher means farther from Normal.",
        }
    return _metric_spec_base(metric)


def suggest_reference_group(groups: List[str]) -> Optional[str]:
    """Keep original behavior for in-file healthy groups only."""
    keywords = ["normal", "healthy", "sham", "blank", "naive", "正常", "健康", "空白", "假手术"]
    normalized = [(g, str(g).strip().lower()) for g in groups]
    for raw, low in normalized:
        if any(kw in low for kw in keywords):
            return raw
    return None


def analyze_experiment(
    raw: Dict[Tuple[str, int], pd.DataFrame],
    experiment: str,
    treatment_group: Optional[str] = None,
    comparator_group: Optional[str] = None,
    reference_group: Optional[str] = None,
    winsor_q: float = WINSOR_Q_DEFAULT,
    external_range_df: Optional[pd.DataFrame] = None,
    normal_sample_df: Optional[pd.DataFrame] = None,
) -> Dict[str, Any]:
    weeks = get_weeks(raw, experiment)
    if not weeks:
        raise ValueError(f"实验 {experiment} 没有找到文件。")
    groups = get_groups(raw, experiment)
    if len(groups) < 2:
        raise ValueError(f"实验 {experiment} 至少需要两个分组。")

    if treatment_group is None or comparator_group is None:
        sug_t, sug_c = suggest_roles(groups)
        treatment_group = treatment_group or sug_t
        comparator_group = comparator_group or sug_c

    if treatment_group not in groups or comparator_group not in groups:
        raise ValueError(f"实验 {experiment} 的治疗组或比较组不在分组列表里。")
    if treatment_group == comparator_group:
        raise ValueError(f"实验 {experiment} 的治疗组和比较组不能相同。")

    if reference_group in [NORMAL_REF_NONE, UPLOADED_TRUE_NORMAL_LABEL, "", None]:
        reference_group = None
    if reference_group is None:
        auto_reference = suggest_reference_group(groups)
        if auto_reference is not None and auto_reference not in {treatment_group, comparator_group}:
            reference_group = auto_reference
    if reference_group is not None and reference_group not in groups:
        raise ValueError(f"实验 {experiment} 的健康参考组 {reference_group} 不在分组列表里。")

    features = common_features(raw, experiment)
    if len(features) < 2:
        raise ValueError(f"实验 {experiment} 共同指标太少（{len(features)} 个）。至少需要 2 个共同指标。")

    base_week = anchor_week(raw, experiment)
    params = fit_preprocessor(raw, experiment, features, winsor_q=winsor_q, base_week=base_week)

    ref_params = None
    in_file_true_normal_available = False
    if reference_group is not None:
        ref_params = fit_reference_params(raw, experiment, features, reference_group=reference_group)
        in_file_true_normal_available = True

    external_range_params: Dict[str, ExternalNormalRangeParam] = {}
    external_range_param_df = pd.DataFrame()
    if external_range_df is not None and not external_range_df.empty:
        external_range_params, external_range_param_df = fit_external_normal_range_params(raw, experiment, external_range_df)

    true_normal_ref = None
    true_normal_param_df = pd.DataFrame()
    effective_normal_df = normal_sample_df
    if (effective_normal_df is None or effective_normal_df.empty) and reference_group is not None:
        ref_frames = []
        for w in weeks:
            sub_ref = raw[(experiment, w)][raw[(experiment, w)]["Group"] == reference_group].copy()
            if not sub_ref.empty:
                ref_frames.append(sub_ref)
        if ref_frames:
            effective_normal_df = pd.concat(ref_frames, ignore_index=True)
            effective_normal_df["__source_file__"] = f"in-file reference group: {reference_group}"
    normal_reference_band: Optional[Dict[str, float]] = None
    if effective_normal_df is not None and not effective_normal_df.empty:
        try:
            true_normal_ref, true_normal_param_df = fit_true_normal_reference(raw, experiment, effective_normal_df)
        except ValueError:
            # A Normal file whose markers do not overlap this experiment's panel
            # (e.g. a shared clinical Normal.xlsx reaching an animal study) must
            # not abort the run — fall back to the range/baseline anchors instead.
            true_normal_ref = None
            true_normal_param_df = pd.DataFrame()
        if true_normal_ref is not None:
            normal_reference_band = compute_normal_reference_nmd_band(true_normal_ref, effective_normal_df)

    # Absolute anchor fallback: when no enrolled Normal/healthy arm and no
    # uploaded Normal-range file exist, anchor clinical-lab markers to built-in
    # published reference ranges so an absolute disease-deviation score (NRBS)
    # can replace the purely baseline-relative BRI. No-op for panels whose
    # markers do not match any published lab (e.g. animal studies).
    used_builtin_reference_ranges = False
    if not external_range_params and true_normal_ref is None and reference_group is None:
        default_range_df = build_default_reference_range_df(raw, experiment)
        if not default_range_df.empty:
            external_range_params, external_range_param_df = fit_external_normal_range_params(
                raw, experiment, default_range_df
            )
            used_builtin_reference_ranges = bool(external_range_params)

    if true_normal_ref is not None:
        if normal_sample_df is not None and not normal_sample_df.empty:
            score_mode = "true_normal_file"
        else:
            score_mode = "in_file_true_normal_group"
    elif in_file_true_normal_available and external_range_params:
        score_mode = "true_normal_group + external_normal_range"
    elif in_file_true_normal_available:
        score_mode = "true_normal_group"
    elif external_range_params:
        score_mode = "builtin_reference_range" if used_builtin_reference_ranges else "external_normal_range"
    else:
        score_mode = "anchor_relative_only"

    preproc_rows = []
    summary_rows = []
    comparison_rows = []
    trend_rows = []
    outlier_rows = []
    mechanism_rows = []
    range_feature_rows = []
    true_normal_feature_rows = []
    sample_tables: Dict[int, Tuple[pd.DataFrame, List[str]]] = {}

    for feat, p in params.items():
        preproc_rows.append(
            {
                "Experiment": experiment,
                "Feature": feat,
                "Anchor_week": base_week,
                "Use_log1p_anchor_relative": p.use_log,
                "Anchor_median": p.median,
                "Anchor_MAD": p.mad,
                "Winsor_lower": p.lower,
                "Winsor_upper": p.upper,
                "Reference_group": reference_group if reference_group is not None else NORMAL_REF_NONE,
                "Analysis_mode": score_mode,
            }
        )

    for week in weeks:
        raw_df = raw[(experiment, week)].copy()
        z_df = apply_preprocessor(raw_df, features, params, winsor=True)
        z_df["Experiment"] = experiment
        z_df["Week"] = week
        z_df["BRI"] = compute_bri(z_df, features)
        z_df["HDI"] = z_df["BRI"]

        if ref_params is not None:
            z_df["NHPS"] = compute_normal_reference_score(raw_df, features, ref_params)
        else:
            z_df["NHPS"] = np.nan

        if true_normal_ref is not None:
            md_df = apply_true_normal_reference(raw_df, true_normal_ref)
            for col in md_df.columns:
                z_df[col] = md_df[col].values
            tnf = _true_normal_feature_table(raw_df.assign(Group=raw_df["Group"]), experiment, week, true_normal_ref)
            if not tnf.empty:
                true_normal_feature_rows.append(tnf)
        else:
            z_df["NMD"] = np.nan
            z_df["NMD_sq"] = np.nan
            z_df["NMD_p_chi2"] = np.nan
            z_df["NMD_cov_p_smallsample_F"] = np.nan
            z_df["NMD_covariance_raw"] = np.nan
            z_df["NMD_sq_covariance_raw"] = np.nan
            z_df["NMD_method"] = ""
            z_df["NMD_z_cap"] = np.nan
            z_df["NMD_feature_count"] = 0

        if external_range_params:
            penalty_df, nrbs, nrps, nr_cov = compute_external_normal_range_scores(raw_df, external_range_params)
            z_df["NRBS"] = nrbs
            z_df["NRPS"] = nrps
            z_df["NR_feature_count"] = nr_cov
            rf = _external_range_feature_table(penalty_df, raw_df["Group"], experiment, week, external_range_params)
            if not rf.empty:
                range_feature_rows.append(rf)
        else:
            z_df["NRBS"] = np.nan
            z_df["NRPS"] = np.nan
            z_df["NR_feature_count"] = 0

        z_df["Score_mode"] = score_mode
        z_df["Reference_group"] = reference_group if reference_group is not None else NORMAL_REF_NONE
        z_df["External_range_source"] = external_range_df["Source_file"].iloc[0] if external_range_df is not None and not external_range_df.empty and "Source_file" in external_range_df.columns else ""
        z_df["True_normal_source_files"] = " | ".join(true_normal_ref.source_files) if true_normal_ref is not None else ""
        sample_tables[week] = (z_df.copy(), features)

        mechanism_rows.append(_feature_mechanism_table(z_df, experiment, week, features))

        for group, sub in z_df.groupby("Group"):
            X = sub[features].to_numpy(dtype=float)
            ent_raw = compute_group_entropy(X.copy())
            ent, X_entropy, ent_info = compute_group_entropy_robust(X.copy())
            ent_lo, ent_hi, ent_se = bootstrap_entropy_ci(X_entropy.copy(), n_boot=200, seed=_boot_seed(0, week, group))
            b_lo, b_hi, b_se = bootstrap_mean_ci(sub["BRI"].values, n_boot=200, seed=_boot_seed(1000, week, group))
            n_lo, n_hi, n_se = bootstrap_mean_ci(sub["NHPS"].dropna().values, n_boot=200, seed=_boot_seed(1500, week, group)) if sub["NHPS"].notna().sum() >= 2 else (float("nan"), float("nan"), float("nan"))
            r_lo, r_hi, r_se = bootstrap_mean_ci(sub["NRBS"].dropna().values, n_boot=200, seed=_boot_seed(1750, week, group)) if sub["NRBS"].notna().sum() >= 2 else (float("nan"), float("nan"), float("nan"))
            rp_lo, rp_hi, rp_se = bootstrap_mean_ci(sub["NRPS"].dropna().values, n_boot=200, seed=_boot_seed(1800, week, group)) if sub["NRPS"].notna().sum() >= 2 else (float("nan"), float("nan"), float("nan"))
            md_lo, md_hi, md_se = bootstrap_mean_ci(sub["NMD"].dropna().values, n_boot=200, seed=_boot_seed(1850, week, group)) if sub["NMD"].notna().sum() >= 2 else (float("nan"), float("nan"), float("nan"))
            ent_cov = ent
            ent_ref = float("nan")
            ent_method = "robust_covariance_entropy"
            if sub["NMD"].notna().sum() >= 2:
                nmd_values = sub["NMD"].dropna().to_numpy(dtype=float)
                ent_ref = compute_distance_state_entropy(nmd_values)
                if np.isfinite(ent_ref):
                    ent = ent_ref
                    ent_lo, ent_hi, ent_se = bootstrap_distance_state_entropy_ci(
                        nmd_values, n_boot=200, seed=_boot_seed(1950, week, group)
                    )
                    ent_method = "normal_reference_state_entropy"

            flags = detect_outliers_methods(X.copy())
            if len(X) > 0:
                stacked = np.column_stack([flags[k] for k in flags])
                method_counts = stacked.sum(axis=1)
                consensus = method_counts >= 2
                strong_consensus = method_counts >= ENTROPY_STRONG_OUTLIER_METHODS
            else:
                consensus = np.zeros(0, dtype=bool)
                strong_consensus = np.zeros(0, dtype=bool)

            outlier_rows.append(
                {
                    "Experiment": experiment,
                    "Week": week,
                    "Group": group,
                    "n": len(sub),
                    "n_robust_z": int(flags["robust_z"].sum()),
                    "n_iqr": int(flags["iqr"].sum()),
                    "n_iso": int(flags["iso"].sum()),
                    "n_lof": int(flags["lof"].sum()),
                    "n_rmd": int(flags["rmd"].sum()),
                    "n_consensus": int(consensus.sum()),
                    "pct_consensus": float(100 * consensus.mean()) if len(consensus) > 0 else 0.0,
                    "Consensus_rule": "flagged by >=2 methods",
                    "n_entropy_strong_outlier": int(strong_consensus.sum()),
                    "pct_entropy_strong_outlier": float(100 * strong_consensus.mean()) if len(strong_consensus) > 0 else 0.0,
                    "Entropy_outlier_rule": ent_info["entropy_outlier_rule"],
                }
            )

            summary_rows.append(
                {
                    "Experiment": experiment,
                    "Week": week,
                    "Group": group,
                    "n": len(sub),
                    "Entropy_raw": ent_raw,
                    "Entropy_robust_covariance": ent_cov,
                    "Entropy_reference_state": ent_ref,
                    "Entropy_method": ent_method,
                    "Entropy": ent,
                    "Entropy_ci_low": ent_lo,
                    "Entropy_ci_high": ent_hi,
                    "Entropy_boot_se": ent_se,
                    "Entropy_n_high_confidence_outliers": ent_info["n_entropy_outliers"],
                    "Entropy_pct_high_confidence_outliers": ent_info["pct_entropy_outliers"],
                    "Entropy_outlier_policy": ent_info["entropy_outlier_rule"],
                    "Entropy_cov_n_samples": ent_info.get("entropy_n", len(sub)),
                    "Entropy_cov_feature_count": ent_info.get("entropy_p", len(features)),
                    "Entropy_cov_n_over_p": ent_info.get("entropy_n_over_p", float("nan")),
                    "Entropy_cov_ledoitwolf_shrinkage": ent_info.get("entropy_ledoitwolf_shrinkage", float("nan")),
                    "Entropy_cov_reliability": ent_info.get("entropy_cov_reliability", "undefined"),
                    "BRI_mean": float(sub["BRI"].mean()),
                    "BRI_median": float(sub["BRI"].median()),
                    "BRI_sd": float(sub["BRI"].std(ddof=1)),
                    "BRI_mean_ci_low": b_lo,
                    "BRI_mean_ci_high": b_hi,
                    "BRI_boot_se": b_se,
                    "HDI_mean": float(sub["BRI"].mean()),
                    "HDI_median": float(sub["BRI"].median()),
                    "HDI_sd": float(sub["BRI"].std(ddof=1)),
                    "HDI_mean_ci_low": b_lo,
                    "HDI_mean_ci_high": b_hi,
                    "HDI_boot_se": b_se,
                    "NHPS_mean": float(sub["NHPS"].mean()) if sub["NHPS"].notna().sum() > 0 else float("nan"),
                    "NHPS_median": float(sub["NHPS"].median()) if sub["NHPS"].notna().sum() > 0 else float("nan"),
                    "NHPS_sd": float(sub["NHPS"].std(ddof=1)) if sub["NHPS"].notna().sum() > 1 else float("nan"),
                    "NHPS_mean_ci_low": n_lo,
                    "NHPS_mean_ci_high": n_hi,
                    "NHPS_boot_se": n_se,
                    "NMD_mean": float(sub["NMD"].mean()) if sub["NMD"].notna().sum() > 0 else float("nan"),
                    "NMD_median": float(sub["NMD"].median()) if sub["NMD"].notna().sum() > 0 else float("nan"),
                    "NMD_sd": float(sub["NMD"].std(ddof=1)) if sub["NMD"].notna().sum() > 1 else float("nan"),
                    "NMD_mean_ci_low": md_lo,
                    "NMD_mean_ci_high": md_hi,
                    "NMD_boot_se": md_se,
                    "NMD_feature_count_mean": float(sub["NMD_feature_count"].mean()) if len(sub) > 0 else 0.0,
                    "NMD_method": sub["NMD_method"].dropna().astype(str).iloc[0] if "NMD_method" in sub.columns and sub["NMD_method"].dropna().shape[0] > 0 else "",
                    "NMD_z_cap": float(sub["NMD_z_cap"].dropna().iloc[0]) if "NMD_z_cap" in sub.columns and sub["NMD_z_cap"].dropna().shape[0] > 0 else float("nan"),
                    "NMD_covariance_raw_mean": float(sub["NMD_covariance_raw"].mean()) if "NMD_covariance_raw" in sub.columns and sub["NMD_covariance_raw"].notna().sum() > 0 else float("nan"),
                    "NRBS_mean": float(sub["NRBS"].mean()) if sub["NRBS"].notna().sum() > 0 else float("nan"),
                    "NRBS_median": float(sub["NRBS"].median()) if sub["NRBS"].notna().sum() > 0 else float("nan"),
                    "NRBS_sd": float(sub["NRBS"].std(ddof=1)) if sub["NRBS"].notna().sum() > 1 else float("nan"),
                    "NRBS_mean_ci_low": r_lo,
                    "NRBS_mean_ci_high": r_hi,
                    "NRBS_boot_se": r_se,
                    "NRPS_mean": float(sub["NRPS"].mean()) if sub["NRPS"].notna().sum() > 0 else float("nan"),
                    "NRPS_median": float(sub["NRPS"].median()) if sub["NRPS"].notna().sum() > 0 else float("nan"),
                    "NRPS_sd": float(sub["NRPS"].std(ddof=1)) if sub["NRPS"].notna().sum() > 1 else float("nan"),
                    "NRPS_mean_ci_low": rp_lo,
                    "NRPS_mean_ci_high": rp_hi,
                    "NRPS_boot_se": rp_se,
                    "NR_feature_count_mean": float(sub["NR_feature_count"].mean()) if len(sub) > 0 else 0.0,
                    "Score_mode": score_mode,
                    "Reference_group": reference_group if reference_group is not None else NORMAL_REF_NONE,
                    "True_normal_available": bool(true_normal_ref is not None),
                    "True_normal_features_common": len(true_normal_ref.features) if true_normal_ref is not None else 0,
                    "True_normal_n_samples": true_normal_ref.n_normal if true_normal_ref is not None else 0,
                    "External_range_available": bool(external_range_params),
                    "External_range_features_common": len(external_range_params),
                }
            )

        if treatment_group in z_df["Group"].unique() and comparator_group in z_df["Group"].unique():
            s1 = z_df[z_df["Group"] == treatment_group]
            s2 = z_df[z_df["Group"] == comparator_group]

            bri_diff = float(s1["BRI"].mean() - s2["BRI"].mean())
            try:
                _, p_bri = mannwhitneyu(s1["BRI"].values, s2["BRI"].values, alternative="two-sided")
            except Exception:
                p_bri = float("nan")
            bri_lo, bri_hi, _ = bootstrap_diff_ci(s1["BRI"].values, s2["BRI"].values, np.mean, n_boot=400, seed=RNG_SEED + 2000 + week)

            if s1["NHPS"].notna().sum() >= 2 and s2["NHPS"].notna().sum() >= 2:
                nhps_diff = float(s1["NHPS"].mean() - s2["NHPS"].mean())
                try:
                    _, p_nhps = mannwhitneyu(s1["NHPS"].dropna().values, s2["NHPS"].dropna().values, alternative="two-sided")
                except Exception:
                    p_nhps = float("nan")
                nhps_lo, nhps_hi, _ = bootstrap_diff_ci(s1["NHPS"].dropna().values, s2["NHPS"].dropna().values, np.mean, n_boot=400, seed=RNG_SEED + 2500 + week)
            else:
                nhps_diff = nhps_lo = nhps_hi = p_nhps = float("nan")

            if s1["NMD"].notna().sum() >= 2 and s2["NMD"].notna().sum() >= 2:
                nmd_diff = float(s1["NMD"].mean() - s2["NMD"].mean())
                try:
                    _, p_nmd = mannwhitneyu(s1["NMD"].dropna().values, s2["NMD"].dropna().values, alternative="two-sided")
                except Exception:
                    p_nmd = float("nan")
                nmd_lo, nmd_hi, _ = bootstrap_diff_ci(s1["NMD"].dropna().values, s2["NMD"].dropna().values, np.mean, n_boot=400, seed=RNG_SEED + 2550 + week)
            else:
                nmd_diff = nmd_lo = nmd_hi = p_nmd = float("nan")

            if s1["NRBS"].notna().sum() >= 2 and s2["NRBS"].notna().sum() >= 2:
                nrbs_diff = float(s1["NRBS"].mean() - s2["NRBS"].mean())
                try:
                    _, p_nrbs = mannwhitneyu(s1["NRBS"].dropna().values, s2["NRBS"].dropna().values, alternative="two-sided")
                except Exception:
                    p_nrbs = float("nan")
                nrbs_lo, nrbs_hi, _ = bootstrap_diff_ci(s1["NRBS"].dropna().values, s2["NRBS"].dropna().values, np.mean, n_boot=400, seed=RNG_SEED + 2600 + week)
                nrps_diff = float(s1["NRPS"].mean() - s2["NRPS"].mean())
                try:
                    _, p_nrps = mannwhitneyu(s1["NRPS"].dropna().values, s2["NRPS"].dropna().values, alternative="two-sided")
                except Exception:
                    p_nrps = float("nan")
                nrps_lo, nrps_hi, _ = bootstrap_diff_ci(s1["NRPS"].dropna().values, s2["NRPS"].dropna().values, np.mean, n_boot=400, seed=RNG_SEED + 2650 + week)
            else:
                nrbs_diff = nrbs_lo = nrbs_hi = p_nrbs = float("nan")
                nrps_diff = nrps_lo = nrps_hi = p_nrps = float("nan")

            X1 = s1[features].to_numpy(dtype=float)
            X2 = s2[features].to_numpy(dtype=float)
            e1_raw = compute_group_entropy(X1.copy())
            e2_raw = compute_group_entropy(X2.copy())
            e1, X1_entropy, _ = compute_group_entropy_robust(X1.copy())
            e2, X2_entropy, _ = compute_group_entropy_robust(X2.copy())
            ent_lo, ent_hi, p_ent = bootstrap_diff_ci(
                X1_entropy, X2_entropy, lambda a: compute_group_entropy(a.copy()), n_boot=300, seed=RNG_SEED + 3000 + week
            )
            entropy_method = "robust_covariance_entropy"
            if s1["NMD"].notna().sum() >= 2 and s2["NMD"].notna().sum() >= 2:
                nmd1 = s1["NMD"].dropna().to_numpy(dtype=float)
                nmd2 = s2["NMD"].dropna().to_numpy(dtype=float)
                e1_state = compute_distance_state_entropy(nmd1)
                e2_state = compute_distance_state_entropy(nmd2)
                if np.isfinite(e1_state) and np.isfinite(e2_state):
                    e1 = e1_state
                    e2 = e2_state
                    ent_lo, ent_hi, p_ent = bootstrap_diff_ci(
                        nmd1, nmd2, lambda a: compute_distance_state_entropy(a), n_boot=300, seed=RNG_SEED + 3050 + week
                    )
                    entropy_method = "normal_reference_state_entropy"

            comparison_rows.append(
                {
                    "Experiment": experiment,
                    "Week": week,
                    "Treatment_group": treatment_group,
                    "Comparator_group": comparator_group,
                    "Reference_group": reference_group if reference_group is not None else NORMAL_REF_NONE,
                    "BRI_diff_treat_minus_comp": bri_diff,
                    "BRI_diff_ci_low": bri_lo,
                    "BRI_diff_ci_high": bri_hi,
                    "BRI_p_mannwhitney": p_bri,
                    "HDI_diff_treat_minus_comp": bri_diff,
                    "HDI_diff_ci_low": bri_lo,
                    "HDI_diff_ci_high": bri_hi,
                    "HDI_p_mannwhitney": p_bri,
                    "NHPS_diff_treat_minus_comp": nhps_diff,
                    "NHPS_diff_ci_low": nhps_lo,
                    "NHPS_diff_ci_high": nhps_hi,
                    "NHPS_p_mannwhitney": p_nhps,
                    "NMD_diff_treat_minus_comp": nmd_diff,
                    "NMD_diff_ci_low": nmd_lo,
                    "NMD_diff_ci_high": nmd_hi,
                    "NMD_p_mannwhitney": p_nmd,
                    "NRBS_diff_treat_minus_comp": nrbs_diff,
                    "NRBS_diff_ci_low": nrbs_lo,
                    "NRBS_diff_ci_high": nrbs_hi,
                    "NRBS_p_mannwhitney": p_nrbs,
                    "NRPS_diff_treat_minus_comp": nrps_diff,
                    "NRPS_diff_ci_low": nrps_lo,
                    "NRPS_diff_ci_high": nrps_hi,
                    "NRPS_p_mannwhitney": p_nrps,
                    "Entropy_raw_diff_treat_minus_comp": float(e1_raw - e2_raw),
                    "Entropy_method": entropy_method,
                    "Entropy_diff_treat_minus_comp": float(e1 - e2),
                    "Entropy_diff_ci_low": ent_lo,
                    "Entropy_diff_ci_high": ent_hi,
                    "Entropy_p_boot": p_ent,
                }
            )

    feature_df = pd.concat(mechanism_rows, ignore_index=True) if mechanism_rows else pd.DataFrame()
    feature_cmp_df = _feature_comparison_table(feature_df, experiment, treatment_group, comparator_group)
    nr_feature_df = pd.concat(range_feature_rows, ignore_index=True) if range_feature_rows else pd.DataFrame()
    nr_feature_cmp_df = _feature_comparison_table(nr_feature_df, experiment, treatment_group, comparator_group) if not nr_feature_df.empty else pd.DataFrame()
    true_normal_feature_df = pd.concat(true_normal_feature_rows, ignore_index=True) if true_normal_feature_rows else pd.DataFrame()
    true_normal_feature_cmp_df = _true_normal_feature_comparison_table(true_normal_feature_df, experiment, treatment_group, comparator_group) if not true_normal_feature_df.empty else pd.DataFrame()

    summary_df = pd.DataFrame(summary_rows).sort_values(["Week", "Group"]).reset_index(drop=True)
    comparison_df = pd.DataFrame(comparison_rows).sort_values(["Week"]).reset_index(drop=True)
    preproc_df = pd.DataFrame(preproc_rows).sort_values(["Feature"]).reset_index(drop=True)
    outlier_df = pd.DataFrame(outlier_rows).sort_values(["Week", "Group"]).reset_index(drop=True)

    for group in summary_df["Group"].unique():
        sub = summary_df[summary_df["Group"] == group].sort_values("Week")
        if len(sub) >= 2:
            rho_e, p_e = spearmanr(sub["Week"], sub["Entropy"])
            rho_b, p_b = spearmanr(sub["Week"], sub["BRI_mean"])
            lin_e = linregress(sub["Week"], sub["Entropy"])
            lin_b = linregress(sub["Week"], sub["BRI_mean"])
            if sub["NHPS_mean"].notna().sum() >= 2:
                rho_n, p_n = spearmanr(sub["Week"], sub["NHPS_mean"])
                lin_n = linregress(sub.loc[sub["NHPS_mean"].notna(), "Week"], sub.loc[sub["NHPS_mean"].notna(), "NHPS_mean"])
            else:
                rho_n = p_n = float("nan")
                lin_n = type("obj", (object,), {"slope": float("nan")})()
            if sub["NMD_mean"].notna().sum() >= 2:
                rho_m, p_m = spearmanr(sub["Week"], sub["NMD_mean"])
                lin_m = linregress(sub.loc[sub["NMD_mean"].notna(), "Week"], sub.loc[sub["NMD_mean"].notna(), "NMD_mean"])
            else:
                rho_m = p_m = float("nan")
                lin_m = type("obj", (object,), {"slope": float("nan")})()
            if sub["NRBS_mean"].notna().sum() >= 2:
                rho_r, p_r = spearmanr(sub["Week"], sub["NRBS_mean"])
                lin_r = linregress(sub.loc[sub["NRBS_mean"].notna(), "Week"], sub.loc[sub["NRBS_mean"].notna(), "NRBS_mean"])
            else:
                rho_r = p_r = float("nan")
                lin_r = type("obj", (object,), {"slope": float("nan")})()
        else:
            rho_e = p_e = rho_b = p_b = rho_n = p_n = rho_m = p_m = rho_r = p_r = float("nan")
            lin_e = lin_b = lin_n = lin_m = lin_r = type("obj", (object,), {"slope": float("nan")})()
        trend_rows.append(
            {
                "Experiment": experiment,
                "Group": group,
                "Entropy_spearman_rho": rho_e,
                "Entropy_spearman_p": p_e,
                "Entropy_linear_slope_per_week": lin_e.slope,
                "BRI_spearman_rho": rho_b,
                "BRI_spearman_p": p_b,
                "BRI_linear_slope_per_week": lin_b.slope,
                "HDI_spearman_rho": rho_b,
                "HDI_spearman_p": p_b,
                "HDI_linear_slope_per_week": lin_b.slope,
                "NHPS_spearman_rho": rho_n,
                "NHPS_spearman_p": p_n,
                "NHPS_linear_slope_per_week": lin_n.slope,
                "NMD_spearman_rho": rho_m,
                "NMD_spearman_p": p_m,
                "NMD_linear_slope_per_week": lin_m.slope,
                "NRBS_spearman_rho": rho_r,
                "NRBS_spearman_p": p_r,
                "NRBS_linear_slope_per_week": lin_r.slope,
            }
        )

    trend_df = pd.DataFrame(trend_rows).sort_values(["Group"]).reset_index(drop=True)

    if len(summary_df) >= 3:
        corr_df = summary_df.copy()
        rho_b, pval_b = spearmanr(corr_df["BRI_mean"], corr_df["Entropy"])
        if corr_df["NHPS_mean"].notna().sum() >= 3:
            tmp = corr_df.dropna(subset=["NHPS_mean", "Entropy"])
            rho_n, pval_n = spearmanr(tmp["NHPS_mean"], tmp["Entropy"]) if len(tmp) >= 3 else (float("nan"), float("nan"))
        else:
            rho_n, pval_n = float("nan"), float("nan")
        if corr_df["NMD_mean"].notna().sum() >= 3:
            tmpm = corr_df.dropna(subset=["NMD_mean", "Entropy"])
            rho_m, pval_m = spearmanr(tmpm["NMD_mean"], tmpm["Entropy"]) if len(tmpm) >= 3 else (float("nan"), float("nan"))
        else:
            rho_m, pval_m = float("nan"), float("nan")
        if corr_df["NRBS_mean"].notna().sum() >= 3:
            tmp2 = corr_df.dropna(subset=["NRBS_mean", "Entropy"])
            rho_r, pval_r = spearmanr(tmp2["NRBS_mean"], tmp2["Entropy"]) if len(tmp2) >= 3 else (float("nan"), float("nan"))
        else:
            rho_r, pval_r = float("nan"), float("nan")
    else:
        rho_b = pval_b = rho_n = pval_n = rho_m = pval_m = rho_r = pval_r = float("nan")

    corr_meta = pd.DataFrame(
        [
            {
                "Experiment": experiment,
                "Scope": "Group-time summaries within experiment",
                "Spearman_rho_BRI_vs_Entropy": rho_b,
                "BRI_p_value": pval_b,
                "Spearman_rho_NHPS_vs_Entropy": rho_n,
                "NHPS_p_value": pval_n,
                "Spearman_rho_NMD_vs_Entropy": rho_m,
                "NMD_p_value": pval_m,
                "Spearman_rho_NRBS_vs_Entropy": rho_r,
                "NRBS_p_value": pval_r,
            }
        ]
    )

    overview = pd.DataFrame(
        [
            {
                "Experiment": experiment,
                "Weeks": ", ".join([f"{w}W" for w in weeks]),
                "Anchor_week": base_week,
                "n_files": len(weeks),
                "n_groups": len(groups),
                "Groups": " | ".join(groups),
                "Common_features": len(features),
                "Treatment_group": treatment_group,
                "Comparator_group": comparator_group,
                "Reference_group": reference_group if reference_group is not None else NORMAL_REF_NONE,
                "Score_mode": score_mode,
                "True_normal_source_files": " | ".join(true_normal_ref.source_files) if true_normal_ref is not None else "",
                "True_normal_features_common": len(true_normal_ref.features) if true_normal_ref is not None else 0,
                "True_normal_n_samples": true_normal_ref.n_normal if true_normal_ref is not None else 0,
                "External_range_file": external_range_df["Source_file"].iloc[0] if external_range_df is not None and not external_range_df.empty and "Source_file" in external_range_df.columns else "",
                "External_range_features_common": len(external_range_params),
            }
        ]
    )

    sample_level = pd.concat([sample_tables[w][0] for w in sorted(sample_tables)], ignore_index=True)

    return {
        "experiment": experiment,
        "weeks": weeks,
        "groups": groups,
        "treatment_group": treatment_group,
        "comparator_group": comparator_group,
        "reference_group": reference_group,
        "score_mode": score_mode,
        "anchor_week": base_week,
        "features": features,
        "true_normal_available": bool(true_normal_ref is not None),
        "normal_reference_band": normal_reference_band,
        "true_normal_params": true_normal_param_df,
        "true_normal_feature_table": true_normal_feature_df,
        "true_normal_feature_comparisons": true_normal_feature_cmp_df,
        "external_range_available": bool(external_range_params),
        "builtin_reference_range_used": bool(used_builtin_reference_ranges),
        "external_range_params": external_range_param_df,
        "external_range_feature_table": nr_feature_df,
        "external_range_feature_comparisons": nr_feature_cmp_df,
        "overview": overview,
        "preprocess": preproc_df,
        "summary": summary_df,
        "comparisons": comparison_df,
        "trends": trend_df,
        "outliers": outlier_df,
        "sample_level": sample_level,
        "corr_meta": corr_meta,
        "feature_mechanism": feature_df,
        "feature_comparisons": feature_cmp_df,
    }


def _analyze_all_base(
    raw: Dict[Tuple[str, int], pd.DataFrame],
    role_map: Optional[Dict[str, Dict[str, str]]] = None,
    winsor_q: float = WINSOR_Q_DEFAULT,
    external_range_map: Optional[Dict[str, pd.DataFrame]] = None,
    true_normal_map: Optional[Dict[str, pd.DataFrame]] = None,
) -> Dict[str, Any]:
    role_map = role_map or {}
    external_range_map = external_range_map or {}
    true_normal_map = true_normal_map or {}

    file_overview_rows = []
    for (exp, week), df in raw.items():
        file_overview_rows.append(
            {
                "Experiment": exp,
                "Week": week,
                "Rows": len(df),
                "Groups": " | ".join(sorted(df["Group"].dropna().astype(str).unique().tolist())),
            }
        )
    file_overview = pd.DataFrame(file_overview_rows).sort_values(["Experiment", "Week"]).reset_index(drop=True)

    exp_results = {}
    for exp in get_experiments(raw):
        roles = role_map.get(exp, {})
        exp_results[exp] = analyze_experiment(
            raw=raw,
            experiment=exp,
            treatment_group=roles.get("treatment"),
            comparator_group=roles.get("comparator"),
            reference_group=roles.get("reference"),
            winsor_q=winsor_q,
            external_range_df=external_range_map.get(exp),
            normal_sample_df=true_normal_map.get(exp),
        )

    id_summary, id_detail = id_consistency_summary(raw)
    feat_inventory = feature_inventory(raw)

    summary_all = pd.concat([exp_results[exp]["summary"] for exp in exp_results], ignore_index=True)
    comparisons_all = pd.concat([exp_results[exp]["comparisons"] for exp in exp_results], ignore_index=True) if any(
        len(exp_results[exp]["comparisons"]) > 0 for exp in exp_results
    ) else pd.DataFrame()
    trends_all = pd.concat([exp_results[exp]["trends"] for exp in exp_results], ignore_index=True)
    preprocess_all = pd.concat([exp_results[exp]["preprocess"] for exp in exp_results], ignore_index=True)
    outliers_all = pd.concat([exp_results[exp]["outliers"] for exp in exp_results], ignore_index=True)
    sample_level_all = pd.concat([exp_results[exp]["sample_level"] for exp in exp_results], ignore_index=True)
    overview_all = pd.concat([exp_results[exp]["overview"] for exp in exp_results], ignore_index=True)
    corr_all = pd.concat([exp_results[exp]["corr_meta"] for exp in exp_results], ignore_index=True)
    feature_mech_all = pd.concat([exp_results[exp]["feature_mechanism"] for exp in exp_results], ignore_index=True) if any(
        len(exp_results[exp]["feature_mechanism"]) > 0 for exp in exp_results
    ) else pd.DataFrame()
    feature_cmp_all = pd.concat([exp_results[exp]["feature_comparisons"] for exp in exp_results], ignore_index=True) if any(
        len(exp_results[exp]["feature_comparisons"]) > 0 for exp in exp_results
    ) else pd.DataFrame()
    tn_param_all = pd.concat([exp_results[exp]["true_normal_params"] for exp in exp_results], ignore_index=True) if any(
        isinstance(exp_results[exp]["true_normal_params"], pd.DataFrame) and len(exp_results[exp]["true_normal_params"]) > 0 for exp in exp_results
    ) else pd.DataFrame()
    tn_feature_all = pd.concat([exp_results[exp]["true_normal_feature_table"] for exp in exp_results], ignore_index=True) if any(
        isinstance(exp_results[exp]["true_normal_feature_table"], pd.DataFrame) and len(exp_results[exp]["true_normal_feature_table"]) > 0 for exp in exp_results
    ) else pd.DataFrame()
    tn_feature_cmp_all = pd.concat([exp_results[exp]["true_normal_feature_comparisons"] for exp in exp_results], ignore_index=True) if any(
        isinstance(exp_results[exp]["true_normal_feature_comparisons"], pd.DataFrame) and len(exp_results[exp]["true_normal_feature_comparisons"]) > 0 for exp in exp_results
    ) else pd.DataFrame()
    nr_param_all = pd.concat([exp_results[exp]["external_range_params"] for exp in exp_results], ignore_index=True) if any(
        isinstance(exp_results[exp]["external_range_params"], pd.DataFrame) and len(exp_results[exp]["external_range_params"]) > 0 for exp in exp_results
    ) else pd.DataFrame()
    nr_feature_all = pd.concat([exp_results[exp]["external_range_feature_table"] for exp in exp_results], ignore_index=True) if any(
        isinstance(exp_results[exp]["external_range_feature_table"], pd.DataFrame) and len(exp_results[exp]["external_range_feature_table"]) > 0 for exp in exp_results
    ) else pd.DataFrame()
    nr_feature_cmp_all = pd.concat([exp_results[exp]["external_range_feature_comparisons"] for exp in exp_results], ignore_index=True) if any(
        isinstance(exp_results[exp]["external_range_feature_comparisons"], pd.DataFrame) and len(exp_results[exp]["external_range_feature_comparisons"]) > 0 for exp in exp_results
    ) else pd.DataFrame()

    if len(summary_all) >= 4 and summary_all["Experiment"].nunique() >= 1:
        pooled_corr = summary_all.copy()
        pooled_corr["BRI_z_within_exp"] = pooled_corr.groupby("Experiment")["BRI_mean"].transform(
            lambda s: (s - s.mean()) / (s.std(ddof=0) if s.std(ddof=0) > 0 else 1.0)
        )
        pooled_corr["Entropy_z_within_exp"] = pooled_corr.groupby("Experiment")["Entropy"].transform(
            lambda s: (s - s.mean()) / (s.std(ddof=0) if s.std(ddof=0) > 0 else 1.0)
        )
        rho_b, pval_b = spearmanr(pooled_corr["BRI_z_within_exp"], pooled_corr["Entropy_z_within_exp"])

        if pooled_corr["NMD_mean"].notna().sum() >= 4:
            tmpm = pooled_corr.dropna(subset=["NMD_mean", "Entropy"]).copy()
            tmpm["NMD_z_within_exp"] = tmpm.groupby("Experiment")["NMD_mean"].transform(
                lambda s: (s - s.mean()) / (s.std(ddof=0) if s.std(ddof=0) > 0 else 1.0)
            )
            tmpm["Entropy_z_within_exp"] = tmpm.groupby("Experiment")["Entropy"].transform(
                lambda s: (s - s.mean()) / (s.std(ddof=0) if s.std(ddof=0) > 0 else 1.0)
            )
            rho_m, pval_m = spearmanr(tmpm["NMD_z_within_exp"], tmpm["Entropy_z_within_exp"]) if len(tmpm) >= 4 else (float("nan"), float("nan"))
        else:
            rho_m, pval_m = float("nan"), float("nan")

        if pooled_corr["NRBS_mean"].notna().sum() >= 4:
            tmp3 = pooled_corr.dropna(subset=["NRBS_mean", "Entropy"]).copy()
            tmp3["NRBS_z_within_exp"] = tmp3.groupby("Experiment")["NRBS_mean"].transform(
                lambda s: (s - s.mean()) / (s.std(ddof=0) if s.std(ddof=0) > 0 else 1.0)
            )
            tmp3["Entropy_z_within_exp"] = tmp3.groupby("Experiment")["Entropy"].transform(
                lambda s: (s - s.mean()) / (s.std(ddof=0) if s.std(ddof=0) > 0 else 1.0)
            )
            rho_r, pval_r = spearmanr(tmp3["NRBS_z_within_exp"], tmp3["Entropy_z_within_exp"]) if len(tmp3) >= 4 else (float("nan"), float("nan"))
        else:
            rho_r, pval_r = float("nan"), float("nan")

        pooled_corr_meta = pd.DataFrame(
            [{
                "Scope": "Pooled after within-experiment z-normalization",
                "Spearman_rho_BRI": rho_b,
                "BRI_p_value": pval_b,
                "Spearman_rho_NMD": rho_m,
                "NMD_p_value": pval_m,
                "Spearman_rho_NRBS": rho_r,
                "NRBS_p_value": pval_r,
            }]
        )
    else:
        pooled_corr_meta = pd.DataFrame(
            [{"Scope": "Pooled after within-experiment z-normalization", "Spearman_rho_BRI": np.nan, "BRI_p_value": np.nan, "Spearman_rho_NMD": np.nan, "NMD_p_value": np.nan, "Spearman_rho_NRBS": np.nan, "NRBS_p_value": np.nan}]
        )

    # --- Benjamini-Hochberg FDR across the reported association tests ---
    # Two pre-specified families: (1) deviation-vs-Entropy associations
    # (within-experiment + pooled), (2) metric-vs-time trends. Each p-value gets
    # a "<col>_fdr_bh" q-value column; HDI mirrors BRI so it is left out of the
    # family to avoid double-counting the same test.
    def _apply_bh_fdr(frames_and_cols):
        entries = []  # (frame_index, pcol, row_idx, value)
        for fi, (df, pcols) in enumerate(frames_and_cols):
            if not isinstance(df, pd.DataFrame):
                continue
            for pcol in pcols:
                if pcol in df.columns:
                    df[pcol + "_fdr_bh"] = np.nan
                    for idx, val in df[pcol].items():
                        if pd.notna(val):
                            entries.append((fi, pcol, idx, float(val)))
        if not entries:
            return
        q = benjamini_hochberg([e[3] for e in entries])
        for (fi, pcol, idx, _), qi in zip(entries, q):
            frames_and_cols[fi][0].loc[idx, pcol + "_fdr_bh"] = qi

    _apply_bh_fdr([
        (corr_all, ["BRI_p_value", "NHPS_p_value", "NMD_p_value", "NRBS_p_value"]),
        (pooled_corr_meta, ["BRI_p_value", "NMD_p_value", "NRBS_p_value"]),
    ])
    _apply_bh_fdr([
        (trends_all, ["Entropy_spearman_p", "BRI_spearman_p", "NHPS_spearman_p",
                      "NMD_spearman_p", "NRBS_spearman_p"]),
    ])

    tables = {
        "实验总览": overview_all,
        "文件总览": file_overview,
        "特征清单": feat_inventory,
        "ID一致性摘要": id_summary,
        "ID一致性明细": id_detail,
        "组时间点汇总": summary_all,
        "治疗比较": comparisons_all,
        "时间趋势": trends_all,
        "预处理参数": preprocess_all,
        "异常值统计": outliers_all,
        "样本级结果": sample_level_all,
        "实验内相关性": corr_all,
        "跨实验合并相关性": pooled_corr_meta,
        "指标机制分解": feature_mech_all,
        "指标比较分解": feature_cmp_all,
    }
    if not tn_param_all.empty:
        tables["TrueNormal参数"] = tn_param_all
    if not tn_feature_all.empty:
        tables["TrueNormal指标偏离"] = tn_feature_all
    if not tn_feature_cmp_all.empty:
        tables["TrueNormal指标比较"] = tn_feature_cmp_all
    if not nr_param_all.empty:
        tables["Normal范围参数"] = nr_param_all
    if not nr_feature_all.empty:
        tables["Normal范围指标负担"] = nr_feature_all
    if not nr_feature_cmp_all.empty:
        tables["Normal范围指标比较"] = nr_feature_cmp_all

    return {"experiments": exp_results, "tables": tables}


def build_plot_label_mapping_from_tables(tables: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    base = _build_plot_label_mapping_base(tables)
    rows = []
    if isinstance(base, pd.DataFrame) and not base.empty:
        rows.append(base)
    for key in ["TrueNormal指标比较", "TrueNormal参数", "Normal范围指标比较", "Normal范围参数"]:
        df = tables.get(key)
        if isinstance(df, pd.DataFrame) and not df.empty and "Experiment" in df.columns:
            for exp in sorted(df["Experiment"].dropna().astype(str).unique().tolist()):
                feats = []
                if "Feature" in df.columns:
                    feats.extend(df.loc[df["Experiment"].astype(str) == exp, "Feature"].dropna().astype(str).tolist())
                if feats:
                    rows.append(_label_mapping_for_experiment_from_features(exp, feats))
    if rows:
        out = pd.concat(rows, ignore_index=True)
        out = out.drop_duplicates(subset=["Experiment", "Original_feature", "Plot_label"]).reset_index(drop=True)
        return out
    return pd.DataFrame(columns=["Experiment", "Original_feature", "Plot_label"])


def _build_markdown_report_ranged(results: Dict[str, Any]) -> str:
    base = _build_markdown_report_base(results)
    lines = [base, "\n\n## 独立 Normal 样本文件模块补充说明\n"]
    lines.append("- 本版本支持上传 **真实的独立 Normal/健康样本文件**，例如 `106-Normal.xlsx`、`109-Normal.xlsx`。")
    lines.append("- 当检测到这类文件时，软件会在与各时间点共同存在的指标上构建真正的 Normal 多变量参考云；主 NMD 使用 **稳健封顶 Normal-reference distance**，原始 LedoitWolf 协方差距离保留在 `NMD_covariance_raw`。")
    lines.append("- 在此基础上新增：")
    lines.append("  - **NMD（True-Normal robust distance）**：样本到真实 Normal 参考云的稳健封顶距离，**越低越健康**；")
    lines.append("  - **deviation-from-Normal heatmap**：显示治疗组相对比较组，哪些指标更接近真实 Normal。")
    lines.append("- NMD 定义为样本相对真实 Normal 参考云的多变量偏离度；越低表示越接近 Normal 参考状态。为避免小样本 Normal 方差过小导致单个指标支配结果，主 NMD 对单指标 z 偏离做封顶处理。")

    tables = results.get("tables", {})
    tn_params = tables.get("TrueNormal参数", pd.DataFrame())
    corr_df = tables.get("实验内相关性", pd.DataFrame())
    overview_df = tables.get("实验总览", pd.DataFrame())

    if isinstance(tn_params, pd.DataFrame) and not tn_params.empty:
        lines.append("\n### 真实 Normal 样本文件概况")
        for exp in sorted(tn_params["Experiment"].dropna().astype(str).unique().tolist()):
            sub = tn_params[tn_params["Experiment"].astype(str) == exp].copy()
            src = sub["Source_files"].dropna().astype(str).unique().tolist() if "Source_files" in sub.columns else []
            unit_notes = [x for x in sub.get("Unit_note", pd.Series(dtype=str)).dropna().astype(str).tolist() if x]
            n_feat = len(sub["Feature"].dropna().astype(str).unique().tolist()) if "Feature" in sub.columns else len(sub)
            lines.append(f"- 实验 **{exp}**：真实 Normal 参考文件 = `{src[0] if src else ''}`；进入主 NMD 的共同特征数 = **{n_feat}**。")
            if unit_notes:
                lines.append(f"  - 单位/语义匹配提示：{'; '.join(sorted(set(unit_notes)))}。")

    if isinstance(corr_df, pd.DataFrame) and not corr_df.empty:
        lines.append("\n### Entropy 与疾病偏离度的关系")
        for _, row in corr_df.iterrows():
            exp = row.get("Experiment", "")
            msg = [f"- 实验 **{exp}**："]
            if pd.notna(row.get("Spearman_rho_BRI_vs_Entropy", np.nan)):
                msg.append(f"BRI vs Entropy rho = {row['Spearman_rho_BRI_vs_Entropy']:.3f} (p = {row['BRI_p_value']:.3g})")
            if pd.notna(row.get("Spearman_rho_NMD_vs_Entropy", np.nan)):
                msg.append(f"NMD vs Entropy rho = {row['Spearman_rho_NMD_vs_Entropy']:.3f} (p = {row['NMD_p_value']:.3g})")
            if pd.notna(row.get("Spearman_rho_NRBS_vs_Entropy", np.nan)):
                msg.append(f"NRBS vs Entropy rho = {row['Spearman_rho_NRBS_vs_Entropy']:.3f} (p = {row['NRBS_p_value']:.3g})")
            lines.append("  " + "；".join(msg[1:]) + "。")

    if isinstance(overview_df, pd.DataFrame) and not overview_df.empty:
        lines.append("\n### 方法学解释边界")
        for _, row in overview_df.iterrows():
            exp = row.get("Experiment", "")
            mode = str(row.get("Score_mode", ""))
            if "true_normal_file" in mode:
                lines.append(f"- 实验 **{exp}**：已启用真实 Normal 样本模式。此时 NMD 可以被解释为“相对真实 Normal 参考云的疾病偏离度”。")
            elif "external_normal_range" in mode:
                lines.append(f"- 实验 **{exp}**：已启用外部 Normal 区间模式。请将 NRBS/NRPS 解释为“相对外部正常区间的贴近程度”；不要把它误写成基于真实 Normal 样本协方差的 Mahalanobis distance。")

    return "\n".join(lines)


def _generate_output_bundle_ranged(results: Dict[str, Any], package_name: str = "BioEntropy_Results") -> Dict[str, Any]:
    bundle = _generate_output_bundle_base(results, package_name=package_name)
    output_dir = Path(bundle["output_dir"])
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    figure_paths = dict(bundle.get("figure_paths", {}))

    for exp, exp_res in results["experiments"].items():
        summary_df = exp_res["summary"]
        comparison_df = exp_res["comparisons"]
        treatment_group = exp_res.get("treatment_group")
        comparator_group = exp_res.get("comparator_group")

        if exp_res.get("true_normal_available"):
            p1 = figures_dir / f"{exp}_true_normal_md_trajectory.png"
            p2 = figures_dir / f"{exp}_true_normal_md_entropy_scatter.png"
            p3 = figures_dir / f"{exp}_true_normal_feature_repair_heatmap.png"
            plot_true_normal_md_trajectory(summary_df, exp, p1, comparison_df, treatment_group, comparator_group)
            plot_true_normal_md_entropy_scatter(summary_df, exp, p2)
            tn_cmp = exp_res.get("true_normal_feature_comparisons", pd.DataFrame())
            if isinstance(tn_cmp, pd.DataFrame) and not tn_cmp.empty:
                plot_true_normal_feature_repair_heatmap(tn_cmp, exp, p3)

            if p1.exists():
                figure_paths[f"{exp}_true_normal_md"] = str(p1)
            if p2.exists():
                figure_paths[f"{exp}_true_normal_md_scatter"] = str(p2)
            if p3.exists():
                figure_paths[f"{exp}_true_normal_feature_repair"] = str(p3)

        if exp_res.get("external_range_available"):
            p4 = figures_dir / f"{exp}_normal_range_burden_trajectory.png"
            p5 = figures_dir / f"{exp}_normal_range_burden_entropy_scatter.png"
            p6 = figures_dir / f"{exp}_normal_range_feature_repair_heatmap.png"
            plot_normal_range_burden_trajectory(summary_df, exp, p4, comparison_df, treatment_group, comparator_group)
            plot_normal_range_burden_entropy_scatter(summary_df, exp, p5)
            nr_cmp = exp_res.get("external_range_feature_comparisons", pd.DataFrame())
            if isinstance(nr_cmp, pd.DataFrame) and not nr_cmp.empty:
                plot_normal_range_feature_repair_heatmap(nr_cmp, exp, p6)
            if p4.exists():
                figure_paths[f"{exp}_normal_range_burden"] = str(p4)
            if p5.exists():
                figure_paths[f"{exp}_normal_range_scatter"] = str(p5)
            if p6.exists():
                figure_paths[f"{exp}_normal_range_feature_repair"] = str(p6)

    zip_path = Path(bundle["zip_path"])
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for path in output_dir.rglob("*"):
            if path.is_file():
                zf.write(path, arcname=str(path.relative_to(output_dir.parent)))

    bundle["figure_paths"] = figure_paths
    return bundle


def build_input_guide_html(out_path: Path) -> Path:
    html = """
<!DOCTYPE html>
<html lang="zh-CN"><head><meta charset="utf-8"><title>BioEntropy 输入说明</title>
<style>
body{font-family:Arial,'Microsoft YaHei UI',sans-serif;max-width:1000px;margin:24px auto;line-height:1.65;color:#222;padding:0 18px}
code{background:#f4f4f4;padding:2px 5px;border-radius:4px}
table{border-collapse:collapse;width:100%;margin:14px 0}
th,td{border:1px solid #ddd;padding:8px 10px;text-align:left}
th{background:#f5f7fb}
.small{color:#666;font-size:13px}
h1,h2{color:#1f4e79}
</style></head><body>
<h1>BioEntropy 输入文件命名与填写格式说明（V6）</h1>

<h2>1. 实验时间点文件命名规则</h2>
<p>实验数据文件请命名为：<code>实验号-周数W.xlsx</code></p>
<table><tr><th>示例</th><th>含义</th></tr>
<tr><td>106-0W.xlsx</td><td>106 实验，0 周</td></tr>
<tr><td>106-24W.xlsx</td><td>106 实验，24 周</td></tr>
<tr><td>109-18W.xlsx</td><td>109 实验，18 周</td></tr></table>
<p>表格格式：第 1 列 <code>Group</code>，第 2 列 <code>Sample ID</code>，后面为各指标列。</p>
<p class="small">动物实验或非纵向单文件数据也可直接上传，例如 <code>DB鼠.xlsx</code>、<code>KK鼠.xlsx</code>、<code>地鼠.xlsx</code>；软件会自动按 0W 单时间点实验处理。</p>
<p class="small">单时间点实验没有真实 Week 维度，结果图会自动使用不同给药组作为横坐标。</p>

<h2>2. 真实 Normal 样本文件命名规则（可选）</h2>
<p>如果你手头有真实的 Normal/健康样本，请额外上传：<code>实验号-Normal.xlsx</code></p>
<table><tr><th>示例</th><th>含义</th></tr>
<tr><td>106-Normal.xlsx</td><td>为 106 实验提供真实 Normal 样本</td></tr>
<tr><td>109-Normal.xlsx</td><td>为 109 实验提供真实 Normal 样本</td></tr></table>
<p class="small">软件允许同一实验上传多个 Normal 样本文件；会自动合并。</p>
<p class="small">当提供真实 Normal 样本时，软件可计算 NMD（True-Normal robust distance）：越低越接近 Normal；原始协方差距离会保留为追溯列。</p>

<h2>3. 外部 Normal 区间文件命名规则（可选）</h2>
<p>如果你没有真实 Normal 样本，但有外部参考区间，请上传：<code>实验号-NormalRange.xlsx</code>（也兼容 <code>实验号-Normal.xlsx</code> 的区间格式）。</p>
<p class="small">这类文件只提供 <code>Normal_min</code> / <code>Normal_max</code>，不会被当作真正的 Normal 样本去估计协方差；软件会输出 NRBS/NRPS。</p>

<h2>4. 真实 Normal 样本文件填写格式</h2>
<table><tr><th>列名</th><th>说明</th></tr>
<tr><td>Group</td><td>建议填写 Normal</td></tr>
<tr><td>Sample ID</td><td>样本编号</td></tr>
<tr><td>第 3 列起</td><td>每个真实 Normal 样本的指标值；一行代表 1 个 Normal 样本</td></tr></table>

<h2>5. 外部区间文件填写格式</h2>
<p>支持两种格式：</p>
<ol>
<li>宽表：两行分别填写 <code>Normal_min</code> 和 <code>Normal_max</code>；</li>
<li>长表：列名包含 <code>Feature</code>、<code>Normal_min</code>、<code>Normal_max</code>。</li>
</ol>

<h2>6. 分析解释边界</h2>
<ul>
<li>没有任何 Normal 参考时：主要输出 Entropy 和 BRI（相对锚定周变化）。</li>
<li>实验数据内含 Normal/健康组时：软件会自动把该组作为 in-file Normal 参考云，输出稳健封顶 NMD；动物实验通常走这个模式。</li>
<li>有独立真实 Normal 样本文件时：输出稳健封顶 NMD，适合解释“离真实 Normal 有多远”。</li>
<li>只有外部区间时：输出 NRBS/NRPS，适合解释“离上传的正常区间有多远”。</li>
<li>主 Entropy 会根据数据条件选择方法：没有 Normal 参考时使用高置信异常值收缩后的稳健协方差熵；有 NMD 时使用 Normal-reference state entropy。原始协方差熵保留在 <code>Entropy_raw</code>。</li>
</ul>
<h2>7. 输出图说明</h2>
<ul>
<li><code>figures</code> 文件夹只保留 7 个主图/主文件：原始马氏距离/疾病偏离度、归一化疾病偏离度展示图、Entropy、二维柱形图、可旋转 3D HTML、二维状态图、药物机制驱动指标图。</li>
<li><code>figures_auxiliary</code> 文件夹会额外输出辅助图：相图、治疗反应排序图、交互式 response landscape。</li>
<li>动物实验或单时间点实验的横坐标为给药组/处理组；临床纵向实验的横坐标为 Week。</li>
<li>归一化疾病偏离度仅用于展示：有 Normal/健康组和 model/模型组时，0–100 尺度采用生物学锚点归一化，即 normal control = 0、model control = 100；真实判断仍看原始数值、置信区间和统计检验。</li>
<li>组别显示顺序可以自定义。桌面端“自定义组排序”可填写全局顺序，也可按实验号分别填写，例如：<code>exp1: normal control, treatment, comparator, model control; exp2: normal control, low-dose, mid-dose, high-dose, comparator, model control</code>。</li>
<li>如果没有真实 Normal 参考，第一张主图会使用 BRI/HDI 作为基线参考疾病偏离度；只有存在 Normal 参考时才严格解释为 True-Normal NMD。</li>
</ul>
</body></html>
"""
    out_path.write_text(html, encoding="utf-8")
    return out_path


def build_excel_template(out_path: Path) -> Path:
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "StudyData"
    headers = [
        "Group", "Sample ID",
        "空腹血浆血糖 (mmol/L)",
        "低密度脂蛋白胆固醇（LDL-C, mmol/L）",
        "高密度脂蛋白胆固醇（HDL-C, mmol/L）",
        "甘油三酯（TG,mmol/L）",
        "总胆固醇（TC, mmol/L）",
        "非高密度脂蛋白胆固醇（mmol/L）",
        "丙氨酸氨基转移酶（ALT, U/L）",
        "天门冬氨酸氨基转移酶（AST,U/L）",
        "γ-谷氨酰转肽酶（GGT, U/L）",
        "超敏C反应蛋白（hsCRP, mg/L）",
        "糖化血红蛋白",
    ]
    ws.append(headers)
    demo_rows = [
        ["Placebo", 1, 8.9, 2.5, 1.0, 1.7, 5.4, 4.4, 48, 55, 33, 2.1, 8.4],
        ["Placebo", 2, 9.3, 2.6, 0.9, 1.9, 5.7, 4.8, 45, 50, 36, 1.8, 8.8],
        ["DrugA", 1, 8.0, 2.2, 1.1, 1.4, 4.9, 3.8, 39, 47, 29, 1.3, 7.6],
        ["DrugA", 2, 7.7, 2.1, 1.2, 1.3, 4.7, 3.5, 37, 45, 27, 1.1, 7.3],
    ]
    for row in demo_rows:
        ws.append(row)
    for cell in ws[1]:
        cell.font = Font(bold=True, color="FFFFFF")
        cell.fill = PatternFill("solid", fgColor="4F81BD")
    for idx, col in enumerate(headers, start=1):
        ws.column_dimensions[get_column_letter(idx)].width = max(14, min(28, len(str(col)) * 1.2))
    note = wb.create_sheet("README")
    note["A1"] = "BioEntropy 实验数据模板说明"
    note["A1"].font = Font(bold=True, size=14)
    note["A3"] = "1. 文件命名请使用：实验号-周数W.xlsx，例如 106-12W.xlsx。"
    note["A4"] = "2. 第一列 Group，第二列 Sample ID。"
    note["A5"] = "3. 第三列起为指标值；不同时间点尽量保持相同指标名称。"
    note["A6"] = "4. 真实 Normal 样本请另存为：实验号-Normal.xlsx。"
    note["A7"] = "5. 外部 min/max 区间请另存为：实验号-NormalRange.xlsx。"
    note.column_dimensions["A"].width = 96
    wb.save(out_path)
    return out_path


def _true_normal_feature_comparison_table(
    feature_table: pd.DataFrame,
    experiment: str,
    treatment_group: str,
    comparator_group: str,
) -> pd.DataFrame:
    if feature_table is None or feature_table.empty:
        return pd.DataFrame()
    rows = []
    for week in sorted(feature_table["Week"].dropna().unique().tolist()):
        t = feature_table[
            (feature_table["Experiment"] == experiment) &
            (feature_table["Week"] == week) &
            (feature_table["Group"] == treatment_group)
        ][["Feature", "Abs_z_to_normal"]].rename(columns={"Abs_z_to_normal": "Abs_z_to_normal_treat"})
        c = feature_table[
            (feature_table["Experiment"] == experiment) &
            (feature_table["Week"] == week) &
            (feature_table["Group"] == comparator_group)
        ][["Feature", "Abs_z_to_normal"]].rename(columns={"Abs_z_to_normal": "Abs_z_to_normal_comp"})
        m = t.merge(c, on="Feature", how="inner")
        if m.empty:
            continue
        m["Abs_z_to_normal_diff_treat_minus_comp"] = m["Abs_z_to_normal_treat"] - m["Abs_z_to_normal_comp"]
        m["Experiment"] = experiment
        m["Week"] = week
        rows.append(m[["Experiment", "Week", "Feature", "Abs_z_to_normal_treat", "Abs_z_to_normal_comp", "Abs_z_to_normal_diff_treat_minus_comp"]])
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


# --------------------------------------------------------------------------- #
# Normal-coverage and attrition tables, and the assembled report
# --------------------------------------------------------------------------- #
from scipy.stats import fisher_exact

def _format_p_value_v62(p: float) -> str:
    if p is None or pd.isna(p):
        return "NA"
    if p < 0.001:
        return "<0.001"
    if p < 0.01:
        return f"{p:.3f}"
    if p < 0.1:
        return f"{p:.3f}"
    return f"{p:.2f}"

def _sig_stars_v62(p: float) -> str:
    if p is None or pd.isna(p):
        return "NA"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"

def _ambiguous_ids_for_experiment_v62(raw: Dict[Tuple[str, int], pd.DataFrame], experiment: str) -> set:
    group_map: Dict[int, set] = {}
    for week in get_weeks(raw, experiment):
        df = raw[(experiment, week)][["Sample ID", "Group"]].dropna().copy()
        if df.empty:
            continue
        df["Sample ID"] = pd.to_numeric(df["Sample ID"], errors="coerce")
        df = df.dropna(subset=["Sample ID"])
        if df.empty:
            continue
        df["Sample ID"] = df["Sample ID"].astype(int)
        for _, row in df.iterrows():
            group_map.setdefault(int(row["Sample ID"]), set()).add(str(row["Group"]))
    return {sid for sid, gs in group_map.items() if len(gs) > 1}

def build_true_normal_coverage_tables(exp_res: Dict[str, Any]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    sample = exp_res.get("sample_level", pd.DataFrame())
    if not isinstance(sample, pd.DataFrame) or sample.empty or "NMD_sq" not in sample.columns:
        return pd.DataFrame(), pd.DataFrame()

    df = sample.copy()
    df["NMD_sq"] = pd.to_numeric(df["NMD_sq"], errors="coerce")
    df["NMD_feature_count"] = pd.to_numeric(df.get("NMD_feature_count"), errors="coerce")
    df["nmd_evaluable"] = df["NMD_sq"].notna() & df["NMD_feature_count"].fillna(0).gt(0)

    def _cut(q, k):
        if pd.isna(k) or k <= 0:
            return np.nan
        return float(chi2.ppf(q, int(k)))

    df["chi2_cut95"] = df["NMD_feature_count"].map(lambda k: _cut(0.95, k))
    df["chi2_cut99"] = df["NMD_feature_count"].map(lambda k: _cut(0.99, k))
    df["within_normal95"] = df["nmd_evaluable"] & (df["NMD_sq"] <= df["chi2_cut95"])
    df["within_normal99"] = df["nmd_evaluable"] & (df["NMD_sq"] <= df["chi2_cut99"])

    summary_rows = []
    for (week, group), sub in df.groupby(["Week", "Group"]):
        n_total = int(len(sub))
        eval_mask = sub["nmd_evaluable"].fillna(False).astype(bool)
        n_eval = int(eval_mask.sum())
        n95 = int(sub.loc[eval_mask, "within_normal95"].sum()) if n_eval > 0 else 0
        n99 = int(sub.loc[eval_mask, "within_normal99"].sum()) if n_eval > 0 else 0
        summary_rows.append({
            "Experiment": exp_res.get("experiment"),
            "Week": int(week),
            "Group": str(group),
            "n_total": n_total,
            "n_evaluable": n_eval,
            "n_within_normal95": n95,
            "pct_within_normal95": float(100.0 * n95 / n_eval) if n_eval > 0 else np.nan,
            "n_within_normal99": n99,
            "pct_within_normal99": float(100.0 * n99 / n_eval) if n_eval > 0 else np.nan,
            "Mean_NMD": float(sub["NMD"].mean()) if sub["NMD"].notna().sum() > 0 else np.nan,
            "Median_NMD": float(sub["NMD"].median()) if sub["NMD"].notna().sum() > 0 else np.nan,
            "Theoretical_cut95_sq_median": float(sub["chi2_cut95"].median()) if sub["chi2_cut95"].notna().sum() > 0 else np.nan,
            "Theoretical_cut99_sq_median": float(sub["chi2_cut99"].median()) if sub["chi2_cut99"].notna().sum() > 0 else np.nan,
        })
    summary_df = pd.DataFrame(summary_rows).sort_values(["Week", "Group"]).reset_index(drop=True)

    comparison_rows = []
    treat = exp_res.get("treatment_group")
    comp = exp_res.get("comparator_group")
    if not summary_df.empty and treat is not None and comp is not None:
        for week in sorted(summary_df["Week"].dropna().unique().tolist()):
            t = summary_df[(summary_df["Week"] == week) & (summary_df["Group"] == treat)]
            c = summary_df[(summary_df["Week"] == week) & (summary_df["Group"] == comp)]
            if t.empty or c.empty:
                continue
            t = t.iloc[0]
            c = c.iloc[0]
            row = {
                "Experiment": exp_res.get("experiment"),
                "Week": int(week),
                "Treatment_group": treat,
                "Comparator_group": comp,
                "pct95_diff_treat_minus_comp": float(t["pct_within_normal95"] - c["pct_within_normal95"]) if pd.notna(t["pct_within_normal95"]) and pd.notna(c["pct_within_normal95"]) else np.nan,
                "pct99_diff_treat_minus_comp": float(t["pct_within_normal99"] - c["pct_within_normal99"]) if pd.notna(t["pct_within_normal99"]) and pd.notna(c["pct_within_normal99"]) else np.nan,
                "n_evaluable_treat": int(t["n_evaluable"]),
                "n_evaluable_comp": int(c["n_evaluable"]),
                "n_within95_treat": int(t["n_within_normal95"]),
                "n_within95_comp": int(c["n_within_normal95"]),
                "n_within99_treat": int(t["n_within_normal99"]),
                "n_within99_comp": int(c["n_within_normal99"]),
            }
            if t["n_evaluable"] > 0 and c["n_evaluable"] > 0:
                a = int(t["n_within_normal95"]); b = int(t["n_evaluable"] - t["n_within_normal95"])
                c1 = int(c["n_within_normal95"]); d = int(c["n_evaluable"] - c["n_within_normal95"])
                try:
                    _, p95 = fisher_exact([[a, b], [c1, d]], alternative="two-sided")
                except Exception:
                    p95 = np.nan
                a = int(t["n_within_normal99"]); b = int(t["n_evaluable"] - t["n_within_normal99"])
                c1 = int(c["n_within_normal99"]); d = int(c["n_evaluable"] - c["n_within_normal99"])
                try:
                    _, p99 = fisher_exact([[a, b], [c1, d]], alternative="two-sided")
                except Exception:
                    p99 = np.nan
            else:
                p95 = p99 = np.nan
            row["p_fisher_95"] = p95
            row["p_fisher_99"] = p99
            comparison_rows.append(row)
    comparison_df = pd.DataFrame(comparison_rows).sort_values(["Week"]).reset_index(drop=True) if comparison_rows else pd.DataFrame()
    return summary_df, comparison_df

def build_true_normal_terminal_overview(exp_res: Dict[str, Any], coverage_df: pd.DataFrame) -> pd.DataFrame:
    summary_df = exp_res.get("summary", pd.DataFrame())
    if not isinstance(summary_df, pd.DataFrame) or summary_df.empty or not exp_res.get("true_normal_available"):
        return pd.DataFrame()
    last_week = int(summary_df["Week"].max())
    sub = summary_df[summary_df["Week"] == last_week].copy()
    keep_cols = [
        "Experiment", "Week", "Group", "n", "NMD_mean", "NMD_mean_ci_low", "NMD_mean_ci_high",
        "Entropy", "Entropy_ci_low", "Entropy_ci_high"
    ]
    keep_cols = [c for c in keep_cols if c in sub.columns]
    out = sub[keep_cols].copy()
    if isinstance(coverage_df, pd.DataFrame) and not coverage_df.empty:
        cov = coverage_df[coverage_df["Week"] == last_week][["Group", "n_evaluable", "n_within_normal95", "pct_within_normal95", "n_within_normal99", "pct_within_normal99"]].copy()
        out = out.merge(cov, on="Group", how="left")
    out["Interpretation"] = np.where(
        out["pct_within_normal95"].fillna(0) >= 50,
        "多数样本进入稳健 Normal 参考包络",
        np.where(out["pct_within_normal95"].fillna(0) > 0, "仅少部分样本进入稳健 Normal 参考包络", "几乎没有样本进入稳健 Normal 参考包络")
    )
    return out.sort_values("NMD_mean").reset_index(drop=True)

def build_attrition_bias_tables(raw: Dict[Tuple[str, int], pd.DataFrame], sample_level_df: pd.DataFrame, experiment: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if not isinstance(sample_level_df, pd.DataFrame) or sample_level_df.empty:
        return pd.DataFrame(), pd.DataFrame()
    weeks = get_weeks(raw, experiment)
    if len(weeks) < 2:
        return pd.DataFrame(), pd.DataFrame()

    base_week = min(weeks)
    terminal_week = max(weeks)
    base = sample_level_df[(sample_level_df["Experiment"].astype(str) == str(experiment)) & (sample_level_df["Week"] == base_week)].copy()
    terminal = sample_level_df[(sample_level_df["Experiment"].astype(str) == str(experiment)) & (sample_level_df["Week"] == terminal_week)].copy()

    if base.empty or terminal.empty or "Sample ID" not in base.columns:
        return pd.DataFrame(), pd.DataFrame()

    for df in (base, terminal):
        df["Sample ID"] = pd.to_numeric(df["Sample ID"], errors="coerce")
        df.dropna(subset=["Sample ID"], inplace=True)
        df["Sample ID"] = df["Sample ID"].astype(int)

    ambiguous_ids = _ambiguous_ids_for_experiment_v62(raw, experiment)
    terminal_group_map = terminal.groupby("Sample ID")["Group"].agg(lambda s: " | ".join(sorted(set(s.astype(str))))).to_dict()

    summary_rows = []
    detail_rows = []

    mode_specs = [
        ("raw_by_id", False),
        ("strict_excluding_ambiguous_ids", True),
    ]
    for mode_name, drop_amb in mode_specs:
        base_mode = base.copy()
        terminal_mode = terminal.copy()
        if drop_amb:
            base_mode = base_mode[~base_mode["Sample ID"].isin(ambiguous_ids)].copy()
            terminal_mode = terminal_mode[~terminal_mode["Sample ID"].isin(ambiguous_ids)].copy()

        for group in sorted(base_mode["Group"].astype(str).unique().tolist()):
            g_base = base_mode[base_mode["Group"].astype(str) == str(group)].copy()
            g_term = terminal_mode[terminal_mode["Group"].astype(str) == str(group)].copy()
            term_ids = set(g_term["Sample ID"].tolist())
            if g_base.empty:
                continue
            retained_ids = set(g_base["Sample ID"].tolist()) & term_ids
            dropped_ids = set(g_base["Sample ID"].tolist()) - retained_ids

            retained = g_base[g_base["Sample ID"].isin(retained_ids)].copy()
            dropped = g_base[g_base["Sample ID"].isin(dropped_ids)].copy()

            p = np.nan
            if retained["NMD"].notna().sum() >= 2 and dropped["NMD"].notna().sum() >= 2:
                try:
                    _, p = mannwhitneyu(retained["NMD"].dropna().values, dropped["NMD"].dropna().values, alternative="two-sided")
                except Exception:
                    p = np.nan

            retained_mean = float(retained["NMD"].mean()) if retained["NMD"].notna().sum() > 0 else np.nan
            dropped_mean = float(dropped["NMD"].mean()) if dropped["NMD"].notna().sum() > 0 else np.nan
            delta = float(dropped_mean - retained_mean) if pd.notna(dropped_mean) and pd.notna(retained_mean) else np.nan

            if len(dropped) == 0:
                interp = "基线样本全部保留到终点。"
            elif pd.notna(delta) and delta > 0.5:
                interp = "未保留样本的基线 NMD 更高，提示可能存在组成漂移。"
            elif pd.notna(delta) and delta < -0.5:
                interp = "保留样本的基线 NMD 更高，需核对样本构成。"
            else:
                interp = "保留与未保留样本的基线 NMD 接近。"

            summary_rows.append({
                "Experiment": experiment,
                "Mode": mode_name,
                "Baseline_week": base_week,
                "Terminal_week": terminal_week,
                "Group": str(group),
                "n_baseline": int(len(g_base)),
                "n_retained_to_terminal": int(len(retained)),
                "n_not_retained_to_terminal": int(len(dropped)),
                "retention_pct": float(100.0 * len(retained) / len(g_base)) if len(g_base) > 0 else np.nan,
                "Baseline_NMD_retained_mean": retained_mean,
                "Baseline_NMD_not_retained_mean": dropped_mean,
                "Delta_not_retained_minus_retained": delta,
                "p_mannwhitney": p,
                "Ambiguous_ids_excluded_n": int(len(ambiguous_ids)) if drop_amb else 0,
                "Interpretation": interp,
            })

            for _, row in g_base.iterrows():
                sid = int(row["Sample ID"])
                retained_flag = sid in retained_ids
                detail_rows.append({
                    "Experiment": experiment,
                    "Mode": mode_name,
                    "Baseline_week": base_week,
                    "Terminal_week": terminal_week,
                    "Group_at_baseline": str(group),
                    "Sample ID": sid,
                    "Retained_to_terminal_same_group": bool(retained_flag),
                    "Terminal_group_seen": terminal_group_map.get(sid, ""),
                    "Baseline_NMD": float(row["NMD"]) if pd.notna(row["NMD"]) else np.nan,
                    "Baseline_BRI": float(row["BRI"]) if pd.notna(row["BRI"]) else np.nan,
                    "Is_ambiguous_id_any_week": bool(sid in ambiguous_ids),
                })

    summary_df = pd.DataFrame(summary_rows).sort_values(["Mode", "Group"]).reset_index(drop=True) if summary_rows else pd.DataFrame()
    detail_df = pd.DataFrame(detail_rows).sort_values(["Mode", "Group_at_baseline", "Sample ID"]).reset_index(drop=True) if detail_rows else pd.DataFrame()
    return summary_df, detail_df

def analyze_all(
    raw: Dict[Tuple[str, int], pd.DataFrame],
    role_map: Optional[Dict[str, Dict[str, str]]] = None,
    winsor_q: float = WINSOR_Q_DEFAULT,
    external_range_map: Optional[Dict[str, pd.DataFrame]] = None,
    true_normal_map: Optional[Dict[str, pd.DataFrame]] = None,
) -> Dict[str, Any]:
    results = _analyze_all_base(
        raw=raw,
        role_map=role_map,
        winsor_q=winsor_q,
        external_range_map=external_range_map,
        true_normal_map=true_normal_map,
    )
    coverage_tables = []
    coverage_cmp_tables = []
    attr_summary_tables = []
    attr_detail_tables = []
    terminal_tables = []

    for exp, exp_res in results.get("experiments", {}).items():
        if exp_res.get("true_normal_available"):
            cov_sum, cov_cmp = build_true_normal_coverage_tables(exp_res)
            exp_res["true_normal_coverage_summary"] = cov_sum
            exp_res["true_normal_coverage_comparison"] = cov_cmp
            if isinstance(cov_sum, pd.DataFrame) and not cov_sum.empty:
                coverage_tables.append(cov_sum)
            if isinstance(cov_cmp, pd.DataFrame) and not cov_cmp.empty:
                coverage_cmp_tables.append(cov_cmp)

            attr_sum, attr_detail = build_attrition_bias_tables(raw, exp_res.get("sample_level", pd.DataFrame()), exp)
            exp_res["attrition_bias_summary"] = attr_sum
            exp_res["attrition_bias_detail"] = attr_detail
            if isinstance(attr_sum, pd.DataFrame) and not attr_sum.empty:
                attr_summary_tables.append(attr_sum)
            if isinstance(attr_detail, pd.DataFrame) and not attr_detail.empty:
                attr_detail_tables.append(attr_detail)

            terminal_df = build_true_normal_terminal_overview(exp_res, cov_sum)
            exp_res["true_normal_terminal_overview"] = terminal_df
            if isinstance(terminal_df, pd.DataFrame) and not terminal_df.empty:
                terminal_tables.append(terminal_df)
        else:
            exp_res["true_normal_coverage_summary"] = pd.DataFrame()
            exp_res["true_normal_coverage_comparison"] = pd.DataFrame()
            exp_res["attrition_bias_summary"] = pd.DataFrame()
            exp_res["attrition_bias_detail"] = pd.DataFrame()
            exp_res["true_normal_terminal_overview"] = pd.DataFrame()

    if coverage_tables:
        results["tables"]["TrueNormal覆盖率"] = pd.concat(coverage_tables, ignore_index=True)
    if coverage_cmp_tables:
        results["tables"]["TrueNormal覆盖率比较"] = pd.concat(coverage_cmp_tables, ignore_index=True)
    if attr_summary_tables:
        results["tables"]["基线留存偏倚摘要"] = pd.concat(attr_summary_tables, ignore_index=True)
    if attr_detail_tables:
        results["tables"]["基线留存偏倚明细"] = pd.concat(attr_detail_tables, ignore_index=True)
    if terminal_tables:
        results["tables"]["TrueNormal终点概览"] = pd.concat(terminal_tables, ignore_index=True)
    return results

def build_normal_sample_template(out_path: Path) -> Path:
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "NormalSamples"
    headers = [
        "Group", "Sample ID",
        "空腹血浆血糖 (mmol/L)", "餐后2h血浆血糖 (mmol/L)",
        "低密度脂蛋白胆固醇（LDL-C, mmol/L）", "高密度脂蛋白胆固醇（HDL-C, mmol/L）",
        "甘油三酯（TG,mmol/L）", "总胆固醇（TC, mmol/L）",
        "非高密度脂蛋白胆固醇（mmol/L）", "丙氨酸氨基转移酶（ALT, U/L）",
        "天门冬氨酸氨基转移酶（AST,U/L）", "γ-谷氨酰转肽酶（GGT, U/L）",
        "超敏C反应蛋白（hsCRP, mg/L）", "糖化血红蛋白", "肌酐 (umol/L)",
    ]
    ws.append(headers)
    demo_rows = [
        ["Normal", 1, 4.8, 6.2, 1.5, 1.3, 0.7, 3.9, 2.6, 35, 42, 18, 0.8, 5.2, 60],
        ["Normal", 2, 5.0, 6.8, 1.6, 1.2, 0.8, 4.1, 2.9, 28, 38, 22, 0.9, 5.4, 66],
        ["Normal", 3, 4.6, 6.0, 1.4, 1.4, 0.6, 3.8, 2.4, 31, 40, 17, 0.7, 5.1, 58],
    ]
    for row in demo_rows:
        ws.append(row)
    for cell in ws[1]:
        cell.font = Font(bold=True, color="FFFFFF")
        cell.fill = PatternFill("solid", fgColor="4F81BD")
    for idx, col in enumerate(headers, start=1):
        ws.column_dimensions[get_column_letter(idx)].width = max(14, min(28, len(str(col)) * 1.2))
    note = wb.create_sheet("README")
    note["A1"] = "BioEntropy 真实 Normal 样本模板说明"
    note["A1"].font = Font(bold=True, size=14)
    note["A3"] = "1. 文件命名请使用：实验号-Normal.xlsx，例如 106-Normal.xlsx。"
    note["A4"] = "2. 第一列是 Group，建议全部填写 Normal；第二列是 Sample ID。"
    note["A5"] = "3. 第三列起为真实 Normal 样本的指标值，每行代表 1 个真实样本。"
    note["A6"] = "4. 如果同一实验有多个 Normal 样本文件，软件会自动合并。"
    note["A7"] = "5. 如果你手里只有 Normal_min / Normal_max 区间，而不是真实样本，请改用实验号-NormalRange.xlsx。"
    note.column_dimensions["A"].width = 100
    wb.save(out_path)
    return out_path

def generate_output_bundle(results: Dict[str, Any], package_name: str = "BioEntropy_Results") -> Dict[str, Any]:
    bundle = _generate_output_bundle_ranged(results, package_name=package_name)
    output_dir = Path(bundle["output_dir"])
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    figure_paths = dict(bundle.get("figure_paths", {}))

    for exp, exp_res in results.get("experiments", {}).items():
        cov_sum = exp_res.get("true_normal_coverage_summary", pd.DataFrame())
        cov_cmp = exp_res.get("true_normal_coverage_comparison", pd.DataFrame())
        if isinstance(cov_sum, pd.DataFrame) and not cov_sum.empty:
            p_cov = figures_dir / f"{exp}_true_normal_coverage_trajectory.png"
            plot_true_normal_coverage_trajectory(
                cov_sum, exp, p_cov, cov_cmp,
                exp_res.get("treatment_group"), exp_res.get("comparator_group")
            )
            if p_cov.exists():
                figure_paths[f"{exp}_true_normal_coverage"] = str(p_cov)

        attr_sum = exp_res.get("attrition_bias_summary", pd.DataFrame())
        if isinstance(attr_sum, pd.DataFrame) and not attr_sum.empty:
            p_attr = figures_dir / f"{exp}_attrition_bias_diagnostic.png"
            plot_attrition_bias_diagnostic(attr_sum, exp, p_attr)
            if p_attr.exists():
                figure_paths[f"{exp}_attrition_bias"] = str(p_attr)

    # refresh zip
    zip_path = Path(bundle["zip_path"])
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for path in output_dir.rglob("*"):
            if path.is_file():
                zf.write(path, arcname=str(path.relative_to(output_dir.parent)))
    bundle["figure_paths"] = figure_paths
    return bundle


# --------------------------------------------------------------------------- #
# Report wording clean-up applied to the assembled report
# --------------------------------------------------------------------------- #
def build_markdown_report(results: Dict[str, Any]) -> str:
    base = _build_markdown_report_ranged(results)
    cleanup_pairs = [
        ("- NMD 定义为样本相对真实 Normal 参考云的多变量偏离度；越低表示越接近 Normal 参考状态。\n", ""),
        ("；", "。"),
    ]
    for old, new in cleanup_pairs:
        base = base.replace(old, new)

    lines = [base, "\n\n## 新增诊断：Normal coverage 与 ID-based attrition bias\n"]
    tables = results.get("tables", {})
    coverage_df = tables.get("TrueNormal覆盖率", pd.DataFrame())
    coverage_cmp = tables.get("TrueNormal覆盖率比较", pd.DataFrame())
    attr_df = tables.get("基线留存偏倚摘要", pd.DataFrame())
    terminal_df = tables.get("TrueNormal终点概览", pd.DataFrame())

    if isinstance(coverage_df, pd.DataFrame) and not coverage_df.empty:
        lines.append("\n### Normal coverage 解释")
        lines.append("- **Normal coverage** 表示：在给定阈值下，落在稳健 Normal 参考包络内的样本比例。")
        lines.append("- 95% coverage 使用 χ²(0.95, df=k) 阈值；99% coverage 使用 χ²(0.99, df=k) 阈值，其中 k 为进入 NMD 的特征数；主 NMD 已做单指标封顶，因此 coverage 是近似参考。")
        lines.append("- 它回答的是“有多少样本进入 Normal 参考包络附近”，不是治愈率，也不能单独替代疾病严重程度判断。")
        for exp in sorted(coverage_df["Experiment"].astype(str).unique().tolist()):
            sub = coverage_df[coverage_df["Experiment"].astype(str) == exp]
            last_week = int(sub["Week"].max())
            last = sub[sub["Week"] == last_week].copy()
            if not last.empty:
                parts = []
                for _, row in last.sort_values("pct_within_normal95", ascending=False).iterrows():
                    parts.append(f"{row['Group']}: 95% coverage = {row['pct_within_normal95']:.2f}%")
                lines.append(f"- 实验 **{exp}** 终点（{last_week}W）：{'；'.join(parts)}。")

    if isinstance(coverage_cmp, pd.DataFrame) and not coverage_cmp.empty:
        lines.append("\n### 治疗组与比较组的 coverage 比较")
        for exp in sorted(coverage_cmp["Experiment"].astype(str).unique().tolist()):
            sub = coverage_cmp[coverage_cmp["Experiment"].astype(str) == exp].sort_values("Week")
            parts = []
            for _, row in sub.iterrows():
                parts.append(
                    f"{int(row['Week'])}W：Δ95 = {row['pct95_diff_treat_minus_comp']:+.2f} pct-pts, p = {_format_p_value_v62(row['p_fisher_95'])}"
                )
            lines.append(f"- 实验 **{exp}**：{'；'.join(parts)}。")

    if isinstance(attr_df, pd.DataFrame) and not attr_df.empty:
        lines.append("\n### ID-based attrition / composition 诊断")
        lines.append("- 该诊断基于 Sample ID 在基线与终点是否仍出现在同一分组中。")
        lines.append("- 若实验存在跨组重复或切换的 Sample ID，这部分应解释为**组成漂移诊断**，而不是严格的个体级失访因果结论。")
        for exp in sorted(attr_df["Experiment"].astype(str).unique().tolist()):
            strict = attr_df[(attr_df["Experiment"].astype(str) == exp) & (attr_df["Mode"] == "strict_excluding_ambiguous_ids")].copy()
            if strict.empty:
                strict = attr_df[attr_df["Experiment"].astype(str) == exp].copy()
            parts = []
            for _, row in strict.iterrows():
                delta = row["Delta_not_retained_minus_retained"]
                delta_txt = f"{delta:+.2f}" if pd.notna(delta) else "NA"
                parts.append(
                    f"{row['Group']}: retention = {row['retention_pct']:.1f}%, Δ基线NMD(未保留-保留) = {delta_txt}, p = {_format_p_value_v62(row['p_mannwhitney'])}"
                )
            lines.append(f"- 实验 **{exp}**：{'；'.join(parts)}。")

    if isinstance(terminal_df, pd.DataFrame) and not terminal_df.empty:
        lines.append("\n### True-Normal 终点概览")
        for exp in sorted(terminal_df["Experiment"].astype(str).unique().tolist()):
            sub = terminal_df[terminal_df["Experiment"].astype(str) == exp].sort_values("NMD_mean")
            parts = []
            for _, row in sub.iterrows():
                pct = row["pct_within_normal95"]
                pct_text = f"{pct:.2f}%" if pd.notna(pct) else "NA"
                parts.append(f"{row['Group']}: NMD = {row['NMD_mean']:.2f}, 95% coverage = {pct_text}")
            lines.append(f"- 实验 **{exp}**：{'；'.join(parts)}。")

    return "\n".join(lines)


# ======================= V6.4 plotting/layout overrides =======================
from matplotlib.lines import Line2D as _Line2D
from matplotlib.transforms import Bbox as _Bbox
from textwrap import fill as _tw_fill

def _add_top_title_and_note(
    fig,
    title: str,
    note: Optional[str] = None,
    legend_handles: Optional[list] = None,
    legend_labels: Optional[list] = None,
    legend_ncol: int = 2,
    legend_y: float = 0.885,
    title_y: float = 0.985,
    note_y: float = 0.944,
    note_width: int = 115,
    note_fontsize: float = 8.8,
) -> None:
    fig.suptitle(title, fontsize=13.5, y=title_y)
    if note:
        fig.text(
            0.5,
            note_y,
            _tw_fill(str(note), note_width),
            ha="center",
            va="top",
            fontsize=note_fontsize,
            color="#444444",
        )
    if legend_handles and legend_labels:
        fig.legend(
            legend_handles,
            legend_labels,
            loc="upper center",
            bbox_to_anchor=(0.5, legend_y),
            frameon=False,
            ncol=max(1, legend_ncol),
            fontsize=9.8,
            columnspacing=1.8,
            handlelength=2.3,
        )

def _bbox_overlaps_any(candidate_bbox, existing_bboxes, pad: float = 2.0) -> bool:
    if candidate_bbox is None:
        return False
    padded = _Bbox.from_extents(
        candidate_bbox.x0 - pad,
        candidate_bbox.y0 - pad,
        candidate_bbox.x1 + pad,
        candidate_bbox.y1 + pad,
    )
    for bbox in existing_bboxes:
        other = _Bbox.from_extents(bbox.x0 - pad, bbox.y0 - pad, bbox.x1 + pad, bbox.y1 + pad)
        if padded.overlaps(other):
            return True
    return False

def _draw_publication_trajectory(
    summary_df: pd.DataFrame,
    comparison_df: pd.DataFrame | None,
    experiment: str,
    metric: str,
    treatment_group: str | None,
    comparator_group: str | None,
    out_path: Path,
) -> None:
    sub = summary_df[summary_df["Experiment"] == experiment].copy()
    if sub.empty:
        return

    spec = _metric_spec(metric)
    groups = _ordered_groups_for_plot(sub, treatment_group, comparator_group)
    cmap = _group_color_map(groups, treatment_group, comparator_group)
    weeks_unique = pd.to_numeric(sub["Week"], errors="coerce").dropna().unique()

    if len(weeks_unique) <= 1:
        rows = []
        value_scores = {}
        for group in sub["Group"].dropna().astype(str).unique().tolist():
            vals = pd.to_numeric(sub.loc[sub["Group"].astype(str) == str(group), spec["value"]], errors="coerce")
            value_scores[group] = float(vals.mean()) if vals.notna().any() else float("inf")
        groups = sorted(value_scores, key=lambda g: (value_scores.get(g, float("inf")), str(g)))
        for group in groups:
            rows.append(sub[sub["Group"].astype(str) == str(group)].iloc[0])
        if not rows:
            return
        plot_df = pd.DataFrame(rows).reset_index(drop=True)
        x = np.arange(len(groups), dtype=float)
        y = pd.to_numeric(plot_df[spec["value"]], errors="coerce").to_numpy(dtype=float)
        fig, ax = plt.subplots(figsize=(max(8.4, 1.05 * len(groups) + 3.0), 6.3))
        colors = [cmap.get(g, PLOT_COLORS[i % len(PLOT_COLORS)]) for i, g in enumerate(groups)]
        bars = ax.bar(x, y, color=colors, alpha=0.9, width=0.62)
        if spec["low"] in plot_df.columns and spec["high"] in plot_df.columns:
            lows = pd.to_numeric(plot_df[spec["low"]], errors="coerce").to_numpy(dtype=float)
            highs = pd.to_numeric(plot_df[spec["high"]], errors="coerce").to_numpy(dtype=float)
            if np.isfinite(lows).any() and np.isfinite(highs).any():
                yerr = np.vstack([np.maximum(0, y - lows), np.maximum(0, highs - y)])
                ax.errorbar(x, y, yerr=yerr, fmt="none", ecolor="#333333", elinewidth=1.1, capsize=4, zorder=4)
        for bar, raw in zip(bars, y):
            if np.isfinite(raw):
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f"{raw:.2f}",
                        ha="center", va="bottom", fontsize=8.5)
        ax.set_xticks(x)
        ax.set_xticklabels([display_group_name(g, i) for i, g in enumerate(groups)], rotation=20, ha="right")
        ax.set_xlabel("Treatment group")
        ax.set_ylabel(spec["ylabel"])
        ax.grid(axis="y", linestyle="--", alpha=0.28)
        _format_axis(ax)
        finite_y = y[np.isfinite(y)]
        if finite_y.size:
            ymin = min(float(finite_y.min()), 0.0)
            ymax = float(finite_y.max())
            yrange = max(ymax - ymin, 1.0)
            ax.set_ylim(ymin - 0.04 * yrange, ymax + 0.18 * yrange)
        _add_top_title_and_note(
            fig,
            title=f"Experiment {experiment}: {spec['title'].replace('trajectory', 'by treatment group')}",
            note="Single-timepoint experiment: the x-axis shows treatment groups, not weeks. " + spec["note"],
            legend_handles=None,
            legend_labels=None,
            title_y=0.985,
            note_y=0.945,
            note_width=120,
            note_fontsize=8.8,
        )
        fig.tight_layout(rect=[0.06, 0.12, 0.98, 0.82])
        fig.savefig(out_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        return

    fig, ax = plt.subplots(figsize=(9.8, 6.6))
    all_y = []

    for idx, group in enumerate(groups):
        g = sub[sub["Group"] == group].sort_values("Week")
        if g.empty:
            continue
        label = display_group_name(group, idx)
        y = g[spec["value"]].to_numpy(dtype=float)
        low = g[spec["low"]].to_numpy(dtype=float)
        high = g[spec["high"]].to_numpy(dtype=float)
        weeks = g["Week"].to_numpy(dtype=float)
        yerr = np.vstack([np.maximum(0, y - low), np.maximum(0, high - y)])
        ax.errorbar(
            weeks,
            y,
            yerr=yerr,
            marker="o",
            markersize=6.2,
            linewidth=2.15,
            elinewidth=1.1,
            capsize=3.5,
            color=cmap[group],
            label=label,
            alpha=0.96,
            zorder=3,
        )
        all_y.extend(list(y))
        all_y.extend(list(low))
        all_y.extend(list(high))

        y_offset = 7 if group == treatment_group else -12
        for x0, y0 in zip(weeks, y):
            ax.annotate(
                f"{y0:.2f}",
                (x0, y0),
                xytext=(0, y_offset),
                textcoords="offset points",
                ha="center",
                fontsize=8.2,
                color=cmap[group],
                bbox=dict(boxstyle="round,pad=0.14", facecolor="white", edgecolor="none", alpha=0.75),
            )

    if metric == "HDI_mean":
        ax.axhline(0, color="#8c8c8c", linestyle="--", linewidth=1)

    if comparison_df is not None and not comparison_df.empty and treatment_group and comparator_group:
        cmp = comparison_df[comparison_df["Experiment"] == experiment].copy()
        if not cmp.empty:
            sub_map = {}
            for group in [treatment_group, comparator_group]:
                temp = sub[sub["Group"] == group].copy().set_index("Week")
                if not temp.empty:
                    sub_map[group] = temp

            ymin = np.nanmin(all_y) if all_y else 0.0
            ymax = np.nanmax(all_y) if all_y else 1.0
            yrange = max(ymax - ymin, 1e-6)

            for _, row in cmp.sort_values("Week").iterrows():
                week = row["Week"]
                if treatment_group not in sub_map or comparator_group not in sub_map:
                    continue
                if week not in sub_map[treatment_group].index or week not in sub_map[comparator_group].index:
                    continue
                y1 = float(sub_map[treatment_group].loc[week, spec["high"]])
                y2 = float(sub_map[comparator_group].loc[week, spec["high"]])
                y_top = max(y1, y2) + 0.06 * yrange
                label = f"{_sig_star(row[spec['p']])}\nΔ={row[spec['diff']]:+.2f}\n{_format_p_publication(row[spec['p']])}"
                ax.text(
                    week,
                    y_top,
                    label,
                    ha="center",
                    va="bottom",
                    fontsize=7.5,
                    bbox=dict(boxstyle="round,pad=0.20", facecolor="white", edgecolor="#cccccc", alpha=0.94),
                    zorder=4,
                )
                all_y.append(y_top)

    ymin = np.nanmin(all_y) if all_y else 0.0
    ymax = np.nanmax(all_y) if all_y else 1.0
    yrange = max(ymax - ymin, 1e-6)
    ax.set_ylim(ymin - 0.08 * yrange, ymax + 0.24 * yrange)
    ax.set_xlabel("Week")
    ax.set_ylabel(spec["ylabel"])
    ax.margins(x=0.05)
    _format_axis(ax)

    handles, labels = ax.get_legend_handles_labels()
    _add_top_title_and_note(
        fig,
        title=f"Experiment {experiment}: {spec['title']}",
        note=spec["note"],
        legend_handles=handles,
        legend_labels=labels,
        legend_ncol=min(len(labels), 3),
        legend_y=0.888,
        title_y=0.985,
        note_y=0.948,
        note_width=116,
        note_fontsize=8.7,
    )
    fig.tight_layout(rect=[0.04, 0.06, 0.98, 0.80])
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

def plot_true_normal_md_entropy_scatter(summary_df: pd.DataFrame, experiment: str, out_path: Path) -> None:
    sub = summary_df[summary_df["Experiment"] == experiment].copy()
    sub = sub.dropna(subset=["NMD_mean", "Entropy"])
    if sub.empty:
        return

    groups = _ordered_groups_for_plot(sub, None, None)
    cmap = _group_color_map(groups, None, None)

    fig, ax = plt.subplots(figsize=(7.6, 5.8))
    items = []
    for idx, (_, row) in enumerate(sub.iterrows()):
        label = f"{display_group_name(row['Group'], idx)} {int(row['Week'])}W"
        color = cmap.get(row["Group"], PLOT_COLORS[idx % len(PLOT_COLORS)])
        ax.scatter(row["NMD_mean"], row["Entropy"], s=82, color=color, alpha=0.9, zorder=3)
        items.append((float(row["NMD_mean"]), float(row["Entropy"]), label, color))

    ax.margins(x=0.10, y=0.14)
    _annotate_points_nonoverlap(ax, items, fontsize=8.0)

    if len(sub) >= 3:
        rho, p = spearmanr(sub["NMD_mean"], sub["Entropy"])
        ax.text(
            0.02,
            0.98,
            f"Spearman rho={rho:.3f}\np={p:.3g}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            bbox=dict(boxstyle="round,pad=0.22", facecolor="white", edgecolor="#dddddd", alpha=0.92),
            fontsize=9,
        )
    ax.set_xlabel("True-Normal robust distance (NMD)")
    ax.set_ylabel("Entropy")
    ax.set_title(f"Experiment {experiment}: NMD vs Entropy", fontsize=13, pad=12)
    _format_axis(ax)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

def plot_true_normal_coverage_trajectory(
    coverage_df: pd.DataFrame,
    experiment: str,
    out_path: Path,
    comparison_df: Optional[pd.DataFrame] = None,
    treatment_group: Optional[str] = None,
    comparator_group: Optional[str] = None,
) -> None:
    if coverage_df is None or coverage_df.empty:
        return
    sub = coverage_df.copy()
    if "Experiment" in sub.columns:
        sub = sub[sub["Experiment"].astype(str) == str(experiment)].copy()
    if sub.empty:
        return

    groups = _ordered_groups_for_plot(sub, treatment_group, comparator_group)
    cmap = _safe_group_color_map(groups)
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.1), sharey=True)
    specs = [("pct_within_normal95", "95%"), ("pct_within_normal99", "99%")]
    single_timepoint = len(pd.to_numeric(sub["Week"], errors="coerce").dropna().unique()) <= 1
    overall_max = float(np.nanmax(sub[["pct_within_normal95", "pct_within_normal99"]].to_numpy(dtype=float))) if sub.shape[0] > 0 else np.nan
    if not np.isfinite(overall_max):
        overall_max = 5.0
    if overall_max <= 5:
        y_top = 6.0
        ann_fmt = "{:.2f}"
    elif overall_max <= 20:
        y_top = min(28.0, overall_max + 6)
        ann_fmt = "{:.1f}"
    else:
        y_top = min(100.0, overall_max + 12)
        ann_fmt = "{:.0f}"

    legend_handles = []
    legend_labels = []
    for idx, group in enumerate(groups):
        legend_handles.append(_Line2D([0], [0], color=cmap[group], marker="o", linewidth=2.2))
        legend_labels.append(display_group_name(group, idx))

    for ax, (col, label) in zip(axes, specs):
        if single_timepoint:
            rows = []
            plot_groups = []
            for group in groups:
                g = sub[sub["Group"].astype(str) == group]
                if not g.empty:
                    rows.append(g.iloc[0])
                    plot_groups.append(group)
            x = np.arange(len(plot_groups), dtype=float)
            vals = np.asarray([float(pd.to_numeric(row[col], errors="coerce")) for row in rows], dtype=float)
            colors = [cmap.get(g, PLOT_COLORS[i % len(PLOT_COLORS)]) for i, g in enumerate(plot_groups)]
            bars = ax.bar(x, vals, color=colors, alpha=0.9, width=0.62, zorder=3)
            for bar, val in zip(bars, vals):
                if np.isfinite(val):
                    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), ann_fmt.format(val),
                            ha="center", va="bottom", fontsize=8)
            ax.set_xticks(x)
            ax.set_xticklabels([display_group_name(g, i) for i, g in enumerate(plot_groups)], rotation=20, ha="right")
        else:
            for idx, group in enumerate(groups):
                g = sub[sub["Group"].astype(str) == group].sort_values("Week")
                if g.empty:
                    continue
                ax.plot(g["Week"], g[col], marker="o", linewidth=2.2, color=cmap[group], zorder=3)
                dy = 6 + (idx % 2) * 7
                for _, row in g.iterrows():
                    if pd.notna(row[col]):
                        ax.annotate(
                            ann_fmt.format(row[col]),
                            (row["Week"], row[col]),
                            xytext=(0, dy),
                            textcoords="offset points",
                            ha="center",
                            fontsize=8,
                            color=cmap[group],
                            bbox=dict(boxstyle="round,pad=0.12", facecolor="white", edgecolor="none", alpha=0.7),
                        )
        ax.set_title(f"{label} coverage", fontsize=11.5, pad=10)
        ax.set_xlabel("Treatment group" if single_timepoint else "Week")
        ax.set_ylim(-0.2, y_top)
        ax.grid(True, axis="y", linestyle="--", alpha=0.3)
        _format_axis(ax)
        if ax is axes[0]:
            ax.set_ylabel("Normal coverage (%)")
        if comparison_df is not None and not comparison_df.empty and treatment_group and comparator_group:
            comp = comparison_df[comparison_df["Experiment"].astype(str) == str(experiment)].sort_values("Week")
            txt_y = y_top * 0.90
            for _, row in comp.iterrows():
                if single_timepoint:
                    continue
                pcol = "p_fisher_95" if "95" in label else "p_fisher_99"
                dcol = "pct95_diff_treat_minus_comp" if "95" in label else "pct99_diff_treat_minus_comp"
                if pd.notna(row.get(dcol, np.nan)):
                    txt = f"Δ={row[dcol]:+.2f}\np={_format_p_value_v62(row.get(pcol, np.nan))}"
                    ax.text(
                        row["Week"],
                        txt_y,
                        txt,
                        ha="center",
                        va="top",
                        fontsize=7.2,
                        bbox=dict(boxstyle="round,pad=0.22", facecolor="white", edgecolor="#dddddd", alpha=0.92),
                        zorder=4,
                    )

    subtitle = (
        "Coverage = percentage of samples inside the robust Normal reference envelope. "
        + (
            "Single-timepoint experiment: bars compare treatment groups; the two panels correspond to 95% and 99% thresholds."
            if single_timepoint else
            "Solid markers show per-group trajectories; the two panels correspond to 95% and 99% thresholds."
        )
    )
    _add_top_title_and_note(
        fig,
        title=f"Experiment {experiment}: True-Normal coverage by treatment group" if single_timepoint else f"Experiment {experiment}: True-Normal coverage trajectory",
        note=subtitle,
        legend_handles=legend_handles,
        legend_labels=legend_labels,
        legend_ncol=min(len(legend_labels), 3),
        legend_y=0.887,
        title_y=0.985,
        note_y=0.948,
        note_width=120,
        note_fontsize=8.7,
    )
    fig.tight_layout(rect=[0.03, 0.08, 0.98, 0.80])
    fig.savefig(out_path, dpi=220, bbox_inches="tight")
    plt.close(fig)

def plot_attrition_bias_diagnostic(
    attrition_summary: pd.DataFrame,
    experiment: str,
    out_path: Path,
) -> None:
    if attrition_summary is None or attrition_summary.empty:
        return
    sub = attrition_summary.copy()
    if "Experiment" in sub.columns:
        sub = sub[sub["Experiment"].astype(str) == str(experiment)].copy()
    if sub.empty:
        return
    modes = [m for m in ["raw_by_id", "strict_excluding_ambiguous_ids"] if m in sub["Mode"].astype(str).unique().tolist()]
    if not modes:
        modes = sorted(sub["Mode"].astype(str).unique().tolist())
    mode_labels = {
        "raw_by_id": "raw_by_id",
        "strict_excluding_ambiguous_ids": "strict_excluding_ambiguous_ids",
    }
    ncols = len(modes)
    fig, axes = plt.subplots(1, ncols, figsize=(6.5 * ncols, 5.0), squeeze=False)
    axes = axes[0]

    legend_handles = [
        _Line2D([0], [0], color="#1f77b4", marker="s", linewidth=0, markersize=10),
        _Line2D([0], [0], color="#d62728", marker="s", linewidth=0, markersize=10),
    ]
    legend_labels = ["基线后续保留", "基线后续未保留"]

    for ax, mode in zip(axes, modes):
        s = sub[sub["Mode"].astype(str) == mode].copy()
        groups = s["Group"].astype(str).tolist()
        x = np.arange(len(groups))
        retained = s["Baseline_NMD_retained_mean"].to_numpy(dtype=float)
        dropped = s["Baseline_NMD_not_retained_mean"].to_numpy(dtype=float)
        w = 0.36
        ax.bar(x - w/2, retained, width=w, color="#1f77b4", zorder=3)
        ax.bar(x + w/2, dropped, width=w, color="#d62728", zorder=3)
        ymax_global = np.nanmax(np.vstack([retained, dropped])) if len(groups) else np.nan
        if not np.isfinite(ymax_global):
            ymax_global = 1.0
        ax.set_ylim(0, ymax_global * 1.35 + 0.6)

        for i, row in s.reset_index(drop=True).iterrows():
            ymax = np.nanmax([retained[i] if i < len(retained) else np.nan, dropped[i] if i < len(dropped) else np.nan])
            if np.isfinite(ymax):
                delta = row['Delta_not_retained_minus_retained']
                delta_txt = f"{delta:+.2f}" if pd.notna(delta) else "NA"
                txt = f"ret={int(row['n_retained_to_terminal'])}/{int(row['n_baseline'])}\nΔ={delta_txt}\np={_format_p_value_v62(row['p_mannwhitney'])}"
                ax.text(i, ymax + ax.get_ylim()[1]*0.035, txt, ha="center", va="bottom", fontsize=7.3, zorder=4)
        ax.set_xticks(x)
        ax.set_xticklabels([display_group_name(g, i) for i, g in enumerate(groups)])
        ax.set_ylabel("Baseline NMD")
        ax.set_title(mode_labels.get(mode, mode), fontsize=11.2, pad=10)
        ax.grid(True, axis="y", linestyle="--", alpha=0.3, zorder=0)
        _format_axis(ax)

    subtitle = (
        "Bars compare baseline NMD between samples later retained to the terminal week and samples not retained. "
        "Interpret this as a sample-composition diagnostic rather than a strict subject-level dropout claim."
    )
    _add_top_title_and_note(
        fig,
        title=f"Experiment {experiment}: ID-based attrition / composition diagnostic",
        note=subtitle,
        legend_handles=legend_handles,
        legend_labels=legend_labels,
        legend_ncol=2,
        legend_y=0.888,
        title_y=0.985,
        note_y=0.948,
        note_width=122,
        note_fontsize=8.6,
    )
    fig.tight_layout(rect=[0.03, 0.08, 0.98, 0.80])
    fig.savefig(out_path, dpi=220, bbox_inches="tight")
    plt.close(fig)

def plot_attrition_bias_nmd(
    attrition_summary: pd.DataFrame,
    experiment: str,
    out_path: Path,
) -> None:
    if attrition_summary is None or attrition_summary.empty:
        return
    sub = attrition_summary.copy()
    if "Experiment" in sub.columns:
        sub = sub[sub["Experiment"].astype(str) == str(experiment)].copy()
    if sub.empty:
        return

    groups = sub["Baseline_group"].astype(str).tolist()
    x = np.arange(len(groups))
    width = 0.34
    retained = pd.to_numeric(sub["Baseline_NMD_retained_mean"], errors="coerce").fillna(np.nan).to_numpy()
    other = pd.to_numeric(sub["Baseline_NMD_other_mean"], errors="coerce").fillna(np.nan).to_numpy()

    fig, ax = plt.subplots(figsize=(8.9, 5.4))
    ax.bar(x - width/2, retained, width=width, label="Retained same group", color="#4CAF50", zorder=3)
    ax.bar(x + width/2, other, width=width, label="Other IDs", color="#9E9E9E", zorder=3)

    ymax = np.nanmax(np.r_[retained, other]) if np.isfinite(np.r_[retained, other]).any() else 1.0
    ymax = max(1.0, ymax * 1.32 + 0.6)
    ax.set_ylim(0, ymax)
    for i, (rv, ov) in enumerate(zip(retained, other)):
        if np.isfinite(rv):
            ax.text(x[i] - width/2, rv + ymax*0.02, f"{rv:.2f}", ha="center", va="bottom", fontsize=9)
        if np.isfinite(ov):
            ax.text(x[i] + width/2, ov + ymax*0.02, f"{ov:.2f}", ha="center", va="bottom", fontsize=9)
        p = sub.iloc[i].get("Baseline_NMD_p_retained_vs_other", np.nan)
        if pd.notna(p):
            txt = f"p={p:.3g}"
            ax.text(x[i], ymax*0.89, txt, ha="center", va="top", fontsize=9,
                    bbox=dict(boxstyle="round,pad=0.18", facecolor="white", edgecolor="#dddddd", alpha=0.92))

    ax.set_xticks(x)
    ax.set_xticklabels([display_group_name(g, i) for i, g in enumerate(groups)])
    ax.set_ylabel("Baseline NMD mean")
    ax.grid(True, axis="y", alpha=0.25, linestyle="--", zorder=0)
    _format_axis(ax)

    handles, labels = ax.get_legend_handles_labels()
    subtitle = (
        "Bars compare baseline NMD of terminal-week retained same-group IDs versus all other baseline IDs. "
        "If retained IDs start closer to Normal, late-stage convergence may be partly influenced by sample-composition drift."
    )
    _add_top_title_and_note(
        fig,
        title=f"Experiment {experiment}: baseline NMD by final-week status",
        note=subtitle,
        legend_handles=handles,
        legend_labels=labels,
        legend_ncol=2,
        legend_y=0.888,
        title_y=0.985,
        note_y=0.948,
        note_width=122,
        note_fontsize=8.6,
    )
    fig.tight_layout(rect=[0.04, 0.08, 0.98, 0.80])
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ======================= V6.4.1 label-placement speed fix =======================
def _annotate_points_nonoverlap(ax, items, fontsize: float = 8.0) -> None:
    """
    Fast, deterministic label placement without repeated canvas draws.
    Items: iterable of (x, y, label, color)
    """
    items = list(items)
    if not items:
        return
    xs = np.array([x for x, _, _, _ in items], dtype=float)
    ys = np.array([y for _, y, _, _ in items], dtype=float)
    xspan = max(float(np.nanmax(xs) - np.nanmin(xs)), 1e-6)
    yspan = max(float(np.nanmax(ys) - np.nanmin(ys)), 1e-6)
    near_dx = xspan * 0.12
    near_dy = yspan * 0.12
    offsets = [(6, 6), (6, -10), (-6, 6), (-6, -10), (12, 0), (-12, 0), (0, 12), (0, -12)]

    placed = []
    for x, y, label, color in sorted(items, key=lambda t: (t[0], t[1])):
        local_count = 0
        for px, py in placed:
            if abs(x - px) <= near_dx and abs(y - py) <= near_dy:
                local_count += 1
        dx, dy = offsets[local_count % len(offsets)]
        ax.annotate(
            label,
            (x, y),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=fontsize,
            color=color,
            ha="left" if dx >= 0 else "right",
            va="bottom" if dy >= 0 else "top",
            bbox=dict(boxstyle="round,pad=0.16", facecolor="white", edgecolor="none", alpha=0.78),
            arrowprops=dict(arrowstyle="-", lw=0.5, color=color, alpha=0.35, shrinkA=2, shrinkB=2) if local_count > 0 else None,
            annotation_clip=False,
        )
        placed.append((x, y))
