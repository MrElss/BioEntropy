"""Figure generation, group ordering, and display-normalization helpers."""
from __future__ import annotations

import math
import json
from pathlib import Path
from typing import Any, Dict, List, Optional

import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgb
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Patch
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

try:
    import plotly.graph_objects as go
    PLOTLY_AVAILABLE = True
except Exception:
    go = None
    PLOTLY_AVAILABLE = False

import bioentropy_core_legacy as _legacy
from bioentropy_core_legacy import *  # noqa: F401,F403
from bioentropy_constants import (
    DISPLAY_CONTRAST_GAMMA,
    _FALLBACK_COLORS,
    _NMD_COLOR,
    _ENTROPY_COLOR,
    _BRI_COLOR,
    _MODEL_COLOR,
    _NORMAL_COLOR,
    _GRID_COLOR,
    _TEXT_COLOR,
    _LIGHT_BG,
)
import bioentropy_display_state as _dstate

_LEGACY_FORMAT_AXIS = _legacy._format_axis

__all__ = [
    "_remove_horizontal_grid",
    "_format_axis_no_horizontal_grid",
    "_semantic_group_color",
    "_ensure_group_color_map",
    "_apply_academic_axis",
    "_lighten",
    "_normalize_0_100",
    "_is_single_timepoint",
    "_group_key",
    "_control_group_rank",
    "default_group_display_order",
    "parse_group_order_text",
    "_experiment_from_frame",
    "_custom_order_for_frame",
    "_apply_custom_group_order",
    "_ordered_groups_for_plot",
    "_metric_group_rows",
    "_ordered_single_timepoint_groups",
    "_plot_single_metric_by_group",
    "_plot_burden_entropy_single_group_bar",
    "_plot_burden_entropy_single_group_heatmap",
    "_plot_burden_entropy_single_group_interactive_html",
    "build_nmd_entropy_joint_table",
    "_get_weeks_and_groups",
    "_comparison_lookup",
    "_signif_label",
    "plot_entropy_trajectory",
    "plot_entropy_baseline_aligned",
    "plot_true_normal_nmd_trajectory",
    "plot_true_normal_nmd_baseline_aligned",
    "plot_true_normal_nmd_entropy_3d",
    "plot_true_normal_nmd_entropy_grouped_bar",
    "plot_true_normal_nmd_entropy_joint_trajectory",
    "plot_true_normal_nmd_entropy_dual_heatmap",
    "_burden_metadata",
    "plot_burden_entropy_joint_trajectory",
    "_annotate_labels_avoid_overlap",
    "plot_burden_entropy_comparison_arrow",
    "plot_clickable_burden_entropy_state_html",
    "plot_state_time_small_multiples",
    "plot_single_timepoint_nmd_entropy_state",
    "plot_burden_entropy_dual_heatmap",
    "plot_burden_entropy_grouped_bar",
    "_front_sorted_groups_by_low_values",
    "plot_burden_entropy_3d",
    "_hex_to_rgb_tuple",
    "_rgba",
    "_plotly_add_cuboid",
    "plot_burden_entropy_3d_interactive_html",
    "plot_burden_entropy_3d_perspective_html",
    "_preferred_burden_kind",
    "_normalize_within_experiment",
    "_is_healthy_reference_group",
    "_is_model_control_group",
    "_normalize_anchor_linear_0_100",
    "_contrast_enhance_0_100",
    "_normalize_display_window",
    "build_display_normalized_metric_table",
    "plot_normalized_metric_index",
    "plot_normalized_burden_index",
    "plot_normalized_entropy_index",
    "_copy_or_regenerate_main_metric_plot",
    "plot_burden_entropy_main_3d",
    "plot_mechanism_diagram",
    "_plot_feature_mechanism_unavailable",
    "plot_feature_mechanism_drivers",
]


def _stable_html_div_id(out_path) -> str:
    """Deterministic Plotly div id derived from the output filename.

    Plotly otherwise assigns a random uuid to each interactive HTML div, which
    makes the exported files non-reproducible run-to-run (and breaks the
    byte-level regression guard). Basing the id on the file stem keeps it stable
    and unique per figure.
    """
    stem = "".join(ch if ch.isalnum() else "_" for ch in Path(out_path).stem)
    return f"bioentropy_{stem}"


def _remove_horizontal_grid(ax) -> None:
    try:
        ax.yaxis.grid(False)
    except Exception:
        try:
            ax.grid(False, axis="y")
        except Exception:
            pass


def _format_axis_no_horizontal_grid(ax):
    _LEGACY_FORMAT_AXIS(ax)
    _remove_horizontal_grid(ax)


_legacy._format_axis = _format_axis_no_horizontal_grid


def _dose_level(name: str) -> Optional[str]:
    """Dose tier of an arm name, or None when the name carries no dose.

    The single-letter forms are only accepted at the end of the name: a bare
    " h" also occurs inside names such as "100mg/kg HTD1801", where it would
    otherwise mark every arm as high dose.
    """
    for level, words, suffixes in (
        ("high", ("300", "high", "高"), ("-h", "_h", " h")),
        ("mid", ("200", "mid", "中"), ("-m", "_m", " m")),
        ("low", ("100", "low", "低"), ("-l", "_l", " l")),
    ):
        if any(w in name for w in words) or any(name.endswith(s) for s in suffixes):
            return level
    return None


def _named_group_color(group: str) -> Optional[str]:
    """Colour implied by the group name itself, or None if the name says nothing."""
    name = str(group).lower()
    if any(tok in name for tok in ["normal", "healthy", "正常", "control normal"]):
        return _NORMAL_COLOR
    if any(tok in name for tok in ["model", "disease", "db", "kk", "hfd", "模型"]):
        return _MODEL_COLOR
    if any(tok in name for tok in ["placebo", "安慰"]):
        return "#8A3B3B"
    if "metformin" in name or "二甲" in name:
        return "#00A087"
    if "pm" in name or "physical" in name or "物理" in name:
        return "#F39B7F"
    if "bbr" in name:
        return "#7E6148"
    if "udca" in name:
        return "#B09C85"
    if "htd" in name or "1801" in name:
        dose = _dose_level(name)
        if dose is not None:
            return {"high": "#3C5488", "mid": "#4DBBD5", "low": "#91D1C2"}[dose]
        return _NMD_COLOR
    return None


def _semantic_group_color(group: str, idx: int) -> str:
    named = _named_group_color(group)
    return named if named is not None else _FALLBACK_COLORS[idx % len(_FALLBACK_COLORS)]


def _ensure_group_color_map(groups: List[str], treatment_group: Optional[str], comparator_group: Optional[str]) -> Dict[str, str]:
    cmap = {g: _semantic_group_color(g, i) for i, g in enumerate(groups)}
    if treatment_group in cmap:
        cmap[treatment_group] = _NMD_COLOR
    if comparator_group in cmap:
        cmap[comparator_group] = _MODEL_COLOR

    # The role overrides above can hand a group a colour the positional palette
    # already gave to another group, which draws two series in the same colour.
    # Treatment and comparator keep their colour; the remaining groups move to
    # the first unused palette entry, name-derived colours taking precedence.
    roles = [g for g in groups if g in (treatment_group, comparator_group)]
    named = [g for g in groups if g not in roles and _named_group_color(g) is not None]
    rest = [g for g in groups if g not in roles and g not in named]
    used = {cmap[g] for g in roles}
    for g in named + rest:
        if cmap[g] in used:
            free = next((c for c in _FALLBACK_COLORS if c not in used), None)
            if free is not None:
                cmap[g] = free
        used.add(cmap[g])
    return cmap


def _apply_academic_axis(ax, grid_axis: str = "y") -> None:
    ax.set_facecolor(_LIGHT_BG)
    ax.tick_params(colors=_TEXT_COLOR, labelsize=9)
    ax.xaxis.label.set_color(_TEXT_COLOR)
    ax.yaxis.label.set_color(_TEXT_COLOR)
    ax.title.set_color(_TEXT_COLOR)
    ax.grid(False)
    if grid_axis and "x" in str(grid_axis).lower():
        ax.xaxis.grid(True, color=_GRID_COLOR, linestyle="-", linewidth=0.65, alpha=0.78)
    _remove_horizontal_grid(ax)
    for side in ["top", "right"]:
        ax.spines[side].set_visible(False)
    for side in ["left", "bottom"]:
        ax.spines[side].set_color("#A8B3C1")
        ax.spines[side].set_linewidth(0.8)


def _lighten(color: str, amount: float = 0.35) -> tuple[float, float, float]:
    r, g, b = to_rgb(color)
    return (min(1.0, r + (1 - r) * amount), min(1.0, g + (1 - g) * amount), min(1.0, b + (1 - b) * amount))


def _normalize_0_100(series: pd.Series) -> pd.Series:
    s = pd.to_numeric(series, errors="coerce").astype(float)
    finite = s[np.isfinite(s)]
    if finite.empty:
        return pd.Series(np.nan, index=s.index, dtype=float)
    vmin = float(finite.min())
    vmax = float(finite.max())
    if math.isclose(vmin, vmax):
        return pd.Series(50.0, index=s.index, dtype=float)
    return (s - vmin) / (vmax - vmin) * 100.0


def _is_single_timepoint(sub: pd.DataFrame) -> bool:
    if "Week" not in sub.columns:
        return True
    weeks = pd.to_numeric(sub["Week"], errors="coerce").dropna().unique()
    return len(weeks) <= 1


def _group_key(value: Any) -> str:
    return " ".join(str(value).strip().lower().replace("，", ",").split())


def _control_group_rank(group: Any) -> int:
    if _is_healthy_reference_group(group):
        return 0
    if _is_model_control_group(group):
        return 1
    return 2


def default_group_display_order(groups: List[str]) -> List[str]:
    indexed = list(enumerate(groups))
    ordered = sorted(indexed, key=lambda item: (_control_group_rank(item[1]), item[0]))
    return [group for _, group in ordered]


def parse_group_order_text(text: str) -> Dict[str, List[str]]:
    """Parse optional user-defined group ordering text.

    Supported formats:
    - Global: normal control, treatment, comparator, model control
    - Per experiment: exp1: normal control, treatment, comparator; exp2: ...
    """
    raw = str(text or "").strip()
    if not raw:
        return {}
    normalized = raw.replace("；", ";").replace("\n", ";").replace("＞", ">").replace(">", ",").replace("，", ",")
    out: Dict[str, List[str]] = {}
    global_items: List[str] = []
    for segment in [s.strip() for s in normalized.split(";") if s.strip()]:
        exp_name: Optional[str] = None
        body = segment
        for sep in (":", "："):
            if sep in segment:
                left, right = segment.split(sep, 1)
                if left.strip() and right.strip():
                    exp_name = left.strip()
                    body = right.strip()
                break
        items = [x.strip() for x in body.split(",") if x.strip()]
        if not items:
            continue
        if exp_name:
            out[exp_name] = items
        else:
            global_items.extend(items)
    if global_items:
        out["__all__"] = global_items
    return out


def _experiment_from_frame(sub: pd.DataFrame) -> Optional[str]:
    if isinstance(sub, pd.DataFrame) and "Experiment" in sub.columns:
        vals = sub["Experiment"].dropna().astype(str).unique().tolist()
        if len(vals) == 1:
            return vals[0]
    return None


def _custom_order_for_frame(sub: pd.DataFrame) -> List[str]:
    if not _dstate.DISPLAY_GROUP_ORDER_MAP:
        return []
    exp = _experiment_from_frame(sub) or _dstate.CURRENT_DISPLAY_EXPERIMENT
    if exp:
        exp_key = _group_key(exp)
        for key, value in _dstate.DISPLAY_GROUP_ORDER_MAP.items():
            if _group_key(key) == exp_key:
                return value
    return _dstate.DISPLAY_GROUP_ORDER_MAP.get("__all__", [])


def _apply_custom_group_order(groups: List[str], sub: Optional[pd.DataFrame] = None) -> List[str]:
    order = _custom_order_for_frame(sub) if sub is not None else _dstate.DISPLAY_GROUP_ORDER_MAP.get("__all__", [])
    if not order:
        return default_group_display_order(groups)
    remaining = list(groups)
    ordered: List[str] = []
    for wanted in order:
        wanted_key = _group_key(wanted)
        match = next((g for g in remaining if _group_key(g) == wanted_key), None)
        if match is not None:
            ordered.append(match)
            remaining.remove(match)
    ordered.extend(default_group_display_order(remaining))
    return ordered


def _ordered_groups_for_plot(
    sub: pd.DataFrame,
    treatment_group: Optional[str] = None,
    comparator_group: Optional[str] = None,
) -> List[str]:
    groups = _legacy._ordered_groups_for_plot(sub, treatment_group, comparator_group)
    return _apply_custom_group_order(groups, sub)


def _metric_group_rows(sub: pd.DataFrame, groups: List[str]) -> pd.DataFrame:
    rows = []
    for group in groups:
        g = sub[sub["Group"].astype(str) == str(group)]
        if not g.empty:
            rows.append(g.iloc[0])
    return pd.DataFrame(rows).reset_index(drop=True) if rows else pd.DataFrame()


def _ordered_single_timepoint_groups(sub: pd.DataFrame, value_col: str) -> List[str]:
    groups = sub["Group"].dropna().astype(str).unique().tolist()
    scores: Dict[str, float] = {}
    for group in groups:
        vals = pd.to_numeric(sub.loc[sub["Group"].astype(str) == group, value_col], errors="coerce")
        scores[group] = float(vals.mean()) if vals.notna().any() else float("inf")
    return _apply_custom_group_order(sorted(groups, key=lambda g: (scores.get(g, float("inf")), str(g))), sub)


def _per_sample_points_for_metric(
    experiment: str, metric_col: str, groups: List[str], normalize: Optional[str] = None,
) -> Dict[str, np.ndarray]:
    """Per-subject values to overlay on a single-timepoint group chart.

    Reads the render-pass sample-level frame from display state. NMD / NRBS / BRI
    map to their own per-sample column directly; Entropy has no per-sample value
    (entropy is a group-level covariance quantity), so each subject's
    reference-state entropy contribution 0.5·ln(2πe·NMD²) is used — the same
    formula the bar aggregates over its NMD distribution, so the dots sit on the
    bar's scale. ``normalize`` puts the values on the chart's 0–100 scale:
    ``"display"`` = linear anchor + gamma (the normalized bars), ``"linear"`` =
    linear anchor only (the state plots); ``None`` keeps raw values.
    """
    sl = getattr(_dstate, "CURRENT_SAMPLE_LEVEL", None)
    if not isinstance(sl, pd.DataFrame) or sl.empty or "Group" not in sl.columns:
        return {}
    sub = sl
    if "Experiment" in sub.columns:
        sub = sub[sub["Experiment"].astype(str) == str(experiment)]
    if sub.empty:
        return {}
    base = str(metric_col)[:-5] if str(metric_col).endswith("_mean") else str(metric_col)
    if base == "Entropy":
        if "NMD" not in sub.columns:
            return {}
        src, transform = "NMD", True
    elif base in sub.columns:
        src, transform = base, False
    else:
        return {}
    out: Dict[str, np.ndarray] = {}
    for g in groups:
        v = pd.to_numeric(sub.loc[sub["Group"].astype(str) == str(g), src], errors="coerce").to_numpy(dtype=float)
        v = v[np.isfinite(v)]
        if transform:
            v = v[v > 0]
            if v.size:
                v = 0.5 * np.log(2 * np.pi * np.e * v * v)
        if v.size:
            out[str(g)] = v
    if normalize and out:
        flat_vals: List[float] = []
        flat_grp: List[str] = []
        for g, arr in out.items():
            flat_vals.extend(float(x) for x in arr)
            flat_grp.extend([g] * len(arr))
        series, grp = pd.Series(flat_vals, dtype=float), pd.Series(flat_grp)
        if normalize == "linear":
            norm = _normalize_anchor_linear_0_100(series, grp).to_numpy(dtype=float)
        else:
            norm = _normalize_display_window(series, grp).to_numpy(dtype=float)
        out2: Dict[str, np.ndarray] = {}
        i = 0
        for g, arr in out.items():
            out2[g] = norm[i:i + len(arr)]
            i += len(arr)
        return out2
    return out


def _overlay_group_points(
    ax, positions, groups: List[str], points_by_group: Dict[str, np.ndarray],
    jitter_width: float = 0.42, size: float = 16, seed: int = 7,
) -> List[float]:
    """Overlay each group's per-subject values as jittered points at the given x
    positions. Purely additive — the underlying chart is untouched. Returns the
    flat list of plotted values so callers can pad the y-limit to keep them
    visible."""
    if not points_by_group:
        return []
    rnd = np.random.RandomState(seed)
    plotted: List[float] = []
    for xi, g in zip(positions, groups):
        pts = points_by_group.get(str(g))
        if pts is None or len(pts) == 0:
            continue
        pts = np.asarray(pts, dtype=float)
        pts = pts[np.isfinite(pts)]
        if pts.size == 0:
            continue
        jit = (rnd.rand(pts.size) - 0.5) * jitter_width if pts.size > 1 else np.zeros(pts.size)
        ax.scatter(np.full(pts.size, float(xi)) + jit, pts, s=size, facecolors="#1f2933",
                   edgecolors="white", linewidths=0.5, alpha=0.55, zorder=6)
        plotted.extend(float(v) for v in pts)
    return plotted


def _overlay_trajectory_points(
    ax, experiment: str, metric_col: str, weeks: List[int], groups: List[str],
    cmap: Dict[str, str], size: float = 7, alpha: float = 0.16, jitter: float = 1.4,
    normalize: Optional[str] = None,
) -> List[float]:
    """Overlay each subject's per-week value (jittered in x) under a longitudinal
    trajectory, coloured by group. Reads the render-pass sample-level frame from
    display state. Points sit behind the lines (low zorder / alpha) so the group
    trajectories stay legible. Returns the flat list of plotted values so callers
    can pad the y-limit. NMD / NRBS / BRI use their own per-sample column;
    Entropy uses each subject's reference-state contribution 0.5·ln(2πe·NMD²).
    ``normalize="display"`` puts the values on the 0–100 display scale (linear
    anchor + gamma) so they share the normalized charts' y-axis."""
    sl = getattr(_dstate, "CURRENT_SAMPLE_LEVEL", None)
    if not isinstance(sl, pd.DataFrame) or sl.empty:
        return []
    sub = sl
    if "Experiment" in sub.columns:
        sub = sub[sub["Experiment"].astype(str) == str(experiment)]
    if sub.empty or "Group" not in sub.columns or "Week" not in sub.columns:
        return []
    base = str(metric_col)[:-5] if str(metric_col).endswith("_mean") else str(metric_col)
    if base == "Entropy":
        if "NMD" not in sub.columns:
            return []
        src, transform = "NMD", True
    elif base in sub.columns:
        src, transform = base, False
    else:
        return []
    wk = pd.to_numeric(sub["Week"], errors="coerce")
    # Pass 1: collect each (group, week) transformed array.
    cells: Dict[Tuple[int, int], np.ndarray] = {}
    flat_vals: List[float] = []
    flat_grp: List[str] = []
    for gi, group in enumerate(groups):
        gmask = sub["Group"].astype(str) == str(group)
        for w in weeks:
            v = pd.to_numeric(sub.loc[gmask & (wk == int(w)), src], errors="coerce").to_numpy(dtype=float)
            v = v[np.isfinite(v)]
            if transform:
                v = v[v > 0]
                if v.size:
                    v = 0.5 * np.log(2 * np.pi * np.e * v * v)
            if v.size == 0:
                continue
            cells[(gi, int(w))] = v
            flat_vals.extend(float(x) for x in v)
            flat_grp.extend([str(group)] * v.size)
    if not cells:
        return []
    if normalize:
        norm = _normalize_display_window(pd.Series(flat_vals, dtype=float), pd.Series(flat_grp)).to_numpy(dtype=float)
        i = 0
        for key in list(cells):
            n = cells[key].size
            cells[key] = norm[i:i + n]
            i += n
    # Pass 2: plot with jitter.
    rnd = np.random.RandomState(7)
    plotted: List[float] = []
    for gi, group in enumerate(groups):
        color = cmap.get(group, _FALLBACK_COLORS[gi % len(_FALLBACK_COLORS)])
        for w in weeks:
            v = cells.get((gi, int(w)))
            if v is None or v.size == 0:
                continue
            jx = float(w) + (rnd.rand(v.size) - 0.5) * jitter
            ax.scatter(jx, v, s=size, c=[color], edgecolors="none", alpha=alpha, zorder=1)
            plotted.extend(float(x) for x in v)
    return plotted


def _plot_single_metric_by_group(
    summary_df: pd.DataFrame,
    experiment: str,
    out_path: Path,
    metric_col: str,
    ylabel: str,
    title: str,
    note: str,
    treatment_group: Optional[str] = None,
    comparator_group: Optional[str] = None,
    low_col: Optional[str] = None,
    high_col: Optional[str] = None,
) -> bool:
    sub = summary_df.copy()
    if "Experiment" in sub.columns:
        sub = sub[sub["Experiment"].astype(str) == str(experiment)].copy()
    keep = [c for c in ["Week", "Group", metric_col, low_col, high_col] if c and c in sub.columns]
    sub = sub[keep].dropna(subset=[metric_col]).copy()
    if sub.empty or not _is_single_timepoint(sub):
        return False
    groups = _ordered_single_timepoint_groups(sub, metric_col)
    groups = [g for g in groups if not sub[sub["Group"].astype(str) == str(g)].empty]
    plot_df = _metric_group_rows(sub, groups)
    if plot_df.empty:
        return False

    cmap = _ensure_group_color_map(groups, treatment_group, comparator_group)
    x = np.arange(len(groups), dtype=float)
    y = pd.to_numeric(plot_df[metric_col], errors="coerce").to_numpy(dtype=float)
    colors = [cmap.get(g, _FALLBACK_COLORS[i % len(_FALLBACK_COLORS)]) for i, g in enumerate(groups)]
    yerr = None
    if low_col and high_col and low_col in plot_df.columns and high_col in plot_df.columns:
        lows = pd.to_numeric(plot_df[low_col], errors="coerce").to_numpy(dtype=float)
        highs = pd.to_numeric(plot_df[high_col], errors="coerce").to_numpy(dtype=float)
        if np.isfinite(lows).any() and np.isfinite(highs).any():
            yerr = np.vstack([np.maximum(0, y - lows), np.maximum(0, highs - y)])

    fig, ax = plt.subplots(figsize=(max(8.2, 1.05 * len(groups) + 3.0), 6.2))
    bars = ax.bar(x, y, color=colors, alpha=0.9, width=0.62)
    if yerr is not None:
        ax.errorbar(x, y, yerr=yerr, fmt="none", ecolor="#333333", elinewidth=1.1, capsize=4, zorder=4)
    for bar, raw in zip(bars, y):
        if np.isfinite(raw):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f"{raw:.2f}",
                    ha="center", va="bottom", fontsize=8.5)
    # Overlay each subject's raw value as a jittered point on top of its bar.
    # The bar/CI/labels are untouched; the dots are purely additive.
    sample_points = _per_sample_points_for_metric(experiment, metric_col, groups)
    point_vals = _overlay_group_points(ax, x, groups, sample_points)
    ax.set_xticks(x)
    ax.set_xticklabels([display_group_name(g, i) for i, g in enumerate(groups)], rotation=20, ha="right")
    ax.set_xlabel("Treatment group")
    ax.set_ylabel(ylabel)
    _legacy._format_axis(ax)
    _remove_horizontal_grid(ax)
    finite_y = y[np.isfinite(y)]
    if finite_y.size:
        hi_candidates = [float(finite_y.max())] + ([max(point_vals)] if point_vals else [])
        lo_candidates = [float(finite_y.min()), 0.0] + ([min(point_vals)] if point_vals else [])
        ymin = min(lo_candidates)
        ymax = max(hi_candidates)
        span = max(ymax - ymin, 1.0)
        ax.set_ylim(ymin - 0.04 * span, ymax + 0.18 * span)
    _legacy._add_top_title_and_note(
        fig,
        title=title,
        note=note,
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
    return True


def _plot_burden_entropy_single_group_bar(
    sub: pd.DataFrame,
    experiment: str,
    out_path: Path,
    burden_col: str,
    meta: Dict[str, Any],
    treatment_group: Optional[str],
    comparator_group: Optional[str],
) -> bool:
    if sub.empty or not _is_single_timepoint(sub):
        return False
    groups = _ordered_single_timepoint_groups(sub, burden_col)
    groups = [g for g in groups if not sub[sub["Group"].astype(str) == str(g)].empty]
    plot_df = _metric_group_rows(sub, groups)
    if plot_df.empty:
        return False
    plot_df["Burden_display"] = _normalize_display_window(pd.to_numeric(plot_df[burden_col], errors="coerce"), plot_df.get("Group"))
    plot_df["Entropy_display"] = _normalize_display_window(pd.to_numeric(plot_df["Entropy"], errors="coerce"), plot_df.get("Group"))
    cmap = _ensure_group_color_map(groups, treatment_group, comparator_group)
    x = np.arange(len(groups), dtype=float)
    width = 0.34
    fig, ax = plt.subplots(figsize=(max(9.0, 1.15 * len(groups) + 3.2), 6.4))
    burden_color = _NMD_COLOR if str(meta.get("short", "")).upper() == "NMD" else _BRI_COLOR
    group_edges = [cmap.get(g, _FALLBACK_COLORS[i % len(_FALLBACK_COLORS)]) for i, g in enumerate(groups)]
    b1 = ax.bar(x - width / 2, plot_df["Burden_display"], width=width, color=burden_color, alpha=0.92,
                edgecolor=group_edges, linewidth=0.9)
    b2 = ax.bar(x + width / 2, plot_df["Entropy_display"], width=width, color=_ENTROPY_COLOR, alpha=0.92,
                edgecolor=group_edges, linewidth=0.9)
    for bi, raw in zip(b1, plot_df[burden_col].tolist()):
        if pd.notna(raw):
            ax.text(bi.get_x() + bi.get_width() / 2, bi.get_height() + 2.1, f"{float(raw):.2f}",
                    ha="center", va="bottom", fontsize=8.0)
    for bi, raw in zip(b2, plot_df["Entropy"].tolist()):
        if pd.notna(raw):
            ax.text(bi.get_x() + bi.get_width() / 2, bi.get_height() + 2.1, f"{float(raw):.2f}",
                    ha="center", va="bottom", fontsize=8.0)
    # Overlay each subject's value on the two bars (0–100 display scale).
    _overlay_group_points(ax, x - width / 2, groups,
                          _per_sample_points_for_metric(experiment, burden_col, groups, normalize="display"),
                          jitter_width=width * 0.7, size=12)
    _overlay_group_points(ax, x + width / 2, groups,
                          _per_sample_points_for_metric(experiment, "Entropy", groups, normalize="display"),
                          jitter_width=width * 0.7, size=12)
    ax.set_xticks(x)
    ax.set_xticklabels([display_group_name(g, i) for i, g in enumerate(groups)], rotation=20, ha="right")
    ax.set_xlabel("Treatment group")
    ax.set_ylabel("Normalized group mean (0–100)")
    ax.set_ylim(0, 115)
    _legacy._format_axis(ax)
    _apply_academic_axis(ax, "y")
    note = (
        f"Single-timepoint experiment: the x-axis shows treatment groups, not weeks. "
        f"Within each group, the left / right bars show {meta['short']} and Entropy. "
        f"Bar heights use anchored display normalization when Normal and model controls are present (Normal=0, model=100), then gamma={DISPLAY_CONTRAST_GAMMA:g} display contrast; labels above bars are raw values."
    )
    _legacy._add_top_title_and_note(
        fig,
        title=f"Experiment {experiment}: group-wise {meta['short']} and Entropy",
        note=note,
        legend_handles=[Patch(facecolor=burden_color, alpha=0.92), Patch(facecolor=_ENTROPY_COLOR, alpha=0.92)],
        legend_labels=[meta["axis"], "Entropy"],
        legend_ncol=2,
        legend_y=0.885,
        title_y=0.985,
        note_y=0.948,
        note_width=122,
        note_fontsize=8.8,
    )
    fig.tight_layout(rect=[0.05, 0.14, 0.98, 0.80])
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return True


def _plot_burden_entropy_single_group_heatmap(
    sub: pd.DataFrame,
    experiment: str,
    out_path: Path,
    burden_col: str,
    meta: Dict[str, Any],
    treatment_group: Optional[str],
    comparator_group: Optional[str],
) -> bool:
    if sub.empty or not _is_single_timepoint(sub):
        return False
    groups = _ordered_single_timepoint_groups(sub, burden_col)
    groups = [g for g in groups if not sub[sub["Group"].astype(str) == str(g)].empty]
    plot_df = _metric_group_rows(sub, groups)
    if plot_df.empty:
        return False
    raw = np.vstack([
        pd.to_numeric(plot_df[burden_col], errors="coerce").to_numpy(dtype=float),
        pd.to_numeric(plot_df["Entropy"], errors="coerce").to_numpy(dtype=float),
    ])
    norm = np.vstack([
        _normalize_0_100(pd.Series(raw[0])).to_numpy(dtype=float),
        _normalize_0_100(pd.Series(raw[1])).to_numpy(dtype=float),
    ])
    fig, ax = plt.subplots(figsize=(max(8.8, 1.05 * len(groups) + 3.0), 4.6))
    im = ax.imshow(norm, aspect="auto", cmap="RdYlBu_r", vmin=0, vmax=100)
    ax.set_yticks([0, 1])
    ax.set_yticklabels([meta["short"], "Entropy"])
    ax.set_xticks(np.arange(len(groups)))
    ax.set_xticklabels([display_group_name(g, i) for i, g in enumerate(groups)], rotation=20, ha="right")
    ax.set_xlabel("Treatment group")
    for i in range(raw.shape[0]):
        for j in range(raw.shape[1]):
            if np.isfinite(raw[i, j]):
                ax.text(j, i, f"{raw[i, j]:.2f}", ha="center", va="center", fontsize=8.0, color="#222")
    cbar = fig.colorbar(im, ax=ax, fraction=0.036, pad=0.02)
    cbar.ax.tick_params(labelsize=8)
    _legacy._add_top_title_and_note(
        fig,
        title=f"Experiment {experiment}: group-wise {meta['short']} / Entropy heatmap",
        note="Single-timepoint experiment: columns are treatment groups. Colors are normalized within each metric; numbers are raw values.",
        legend_handles=None,
        legend_labels=None,
        title_y=0.985,
        note_y=0.940,
        note_width=116,
        note_fontsize=8.6,
    )
    fig.tight_layout(rect=[0.05, 0.16, 0.98, 0.78])
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return True


def _plot_burden_entropy_single_group_interactive_html(
    sub: pd.DataFrame,
    experiment: str,
    out_path: Path,
    burden_col: str,
    meta: Dict[str, Any],
    treatment_group: Optional[str],
    comparator_group: Optional[str],
) -> bool:
    if not PLOTLY_AVAILABLE or sub.empty or not _is_single_timepoint(sub):
        return False
    groups = _ordered_single_timepoint_groups(sub, burden_col)
    groups = [g for g in groups if not sub[sub["Group"].astype(str) == str(g)].empty]
    plot_df = _metric_group_rows(sub, groups)
    if plot_df.empty:
        return False
    labels = [display_group_name(g, i) for i, g in enumerate(groups)]
    burden_display = _normalize_0_100(pd.to_numeric(plot_df[burden_col], errors="coerce"))
    entropy_display = _normalize_0_100(pd.to_numeric(plot_df["Entropy"], errors="coerce"))
    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=labels, y=burden_display, name=meta["short"],
        text=[f"{float(v):.2f}" if pd.notna(v) else "" for v in plot_df[burden_col]],
        hovertemplate="Group: %{x}<br>Metric: " + meta["short"] + "<br>Raw: %{text}<br>Display height: %{y:.1f}<extra></extra>",
        marker_color=_NMD_COLOR if str(meta.get("short", "")).upper() == "NMD" else _BRI_COLOR,
    ))
    fig.add_trace(go.Bar(
        x=labels, y=entropy_display, name="Entropy",
        text=[f"{float(v):.2f}" if pd.notna(v) else "" for v in plot_df["Entropy"]],
        hovertemplate="Group: %{x}<br>Metric: Entropy<br>Raw: %{text}<br>Display height: %{y:.1f}<extra></extra>",
        marker_color=_ENTROPY_COLOR,
    ))
    fig.update_layout(
        title=f"Experiment {experiment}: group-wise {meta['short']} and Entropy (interactive)",
        barmode="group",
        xaxis_title="Treatment group",
        yaxis_title="Normalized group mean (0–100)",
        legend=dict(orientation="h", y=1.02, x=0.5, xanchor="center"),
        margin=dict(l=70, r=30, t=105, b=120),
        height=650,
        annotations=[dict(
            text="Single-timepoint experiment: the x-axis shows treatment groups, not weeks. Labels and hover tooltips show raw values.",
            x=0.5, y=0.965, xref="paper", yref="paper", showarrow=False,
            font=dict(size=13, color="#555"), align="center"
        )],
    )
    fig.write_html(str(out_path), include_plotlyjs="cdn", full_html=True, div_id=_stable_html_div_id(out_path))
    return True


def build_nmd_entropy_joint_table(results: Dict[str, Any]) -> pd.DataFrame:
    rows: List[pd.DataFrame] = []
    for exp, exp_res in results.get("experiments", {}).items():
        if not bool(exp_res.get("true_normal_available")):
            continue
        summary_df = exp_res.get("summary", pd.DataFrame())
        if not isinstance(summary_df, pd.DataFrame) or summary_df.empty:
            continue
        keep = [c for c in ["Experiment", "Week", "Group", "n", "NMD_mean", "Entropy"] if c in summary_df.columns]
        sub = summary_df[keep].dropna(subset=["NMD_mean", "Entropy"]).copy()
        if sub.empty:
            continue
        nmd_linear = _normalize_anchor_linear_0_100(sub["NMD_mean"], sub.get("Group"))
        entropy_linear = _normalize_anchor_linear_0_100(sub["Entropy"], sub.get("Group"))
        sub["NMD_linear_0_100"] = nmd_linear.round(4)
        sub["Entropy_linear_0_100"] = entropy_linear.round(4)
        sub["NMD_display_contrast_0_100"] = _contrast_enhance_0_100(nmd_linear).round(4)
        sub["Entropy_display_contrast_0_100"] = _contrast_enhance_0_100(entropy_linear).round(4)
        sub["NMD_display_0_100"] = sub["NMD_display_contrast_0_100"]
        sub["Entropy_display_0_100"] = sub["Entropy_display_contrast_0_100"]
        sub["Joint_burden_0_100"] = ((sub["NMD_display_0_100"] + sub["Entropy_display_0_100"]) / 2.0).round(4)
        sub["Interpretation"] = (
            f"数值越低表示越接近 Normal 且熵越低；linear_0_100 保留 Normal=0/model=100 的线性锚定值，"
            f"display_contrast_0_100 使用 gamma={DISPLAY_CONTRAST_GAMMA:g} 增强展示差异；正式解释不应替代原始 NMD 或 Entropy。"
        )
        if _custom_order_for_frame(sub):
            ordered_groups = _apply_custom_group_order(sub["Group"].dropna().astype(str).unique().tolist(), sub)
            rank_map = {str(g): i for i, g in enumerate(ordered_groups)}
            sub["_Group_order_rank"] = sub["Group"].astype(str).map(rank_map).fillna(9999).astype(int)
        rows.append(sub)
    if not rows:
        return pd.DataFrame()
    out = pd.concat(rows, ignore_index=True)
    sort_cols = [c for c in ["Experiment", "Week", "_Group_order_rank", "Group"] if c in out.columns]
    return out.sort_values(sort_cols).drop(columns=[c for c in ["_Group_order_rank"] if c in out.columns]).reset_index(drop=True)


def _get_weeks_and_groups(summary_df: pd.DataFrame, treatment_group: Optional[str] = None, comparator_group: Optional[str] = None):
    sub = summary_df.copy()
    sub["Week"] = pd.to_numeric(sub["Week"], errors="coerce").astype(int)
    weeks = sorted(sub["Week"].dropna().astype(int).unique().tolist())
    groups = _ordered_groups_for_plot(sub, treatment_group, comparator_group)
    return weeks, groups


def _comparison_lookup(comparison_df: Optional[pd.DataFrame]) -> pd.DataFrame:
    if comparison_df is None or not isinstance(comparison_df, pd.DataFrame) or comparison_df.empty:
        return pd.DataFrame()
    out = comparison_df.copy()
    out["Week"] = pd.to_numeric(out["Week"], errors="coerce").astype(int)
    return out


def _signif_label(p):
    try:
        p = float(p)
    except Exception:
        return "ns"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def plot_entropy_trajectory(summary_df: pd.DataFrame, experiment: str, out_path: Path,
                            comparison_df: Optional[pd.DataFrame] = None,
                            treatment_group: Optional[str] = None,
                            comparator_group: Optional[str] = None) -> None:
    sub = summary_df.copy()
    if "Experiment" in sub.columns:
        sub = sub[sub["Experiment"].astype(str) == str(experiment)].copy()
    keep = [c for c in ["Week", "Group", "Entropy", "Entropy_ci_low", "Entropy_ci_high"] if c in sub.columns]
    sub = sub[keep].dropna(subset=["Entropy"]).copy()
    if sub.empty:
        return
    sub["Week"] = pd.to_numeric(sub["Week"], errors="coerce").astype(int)
    if _is_single_timepoint(sub):
        _plot_single_metric_by_group(
            sub, experiment, out_path,
            metric_col="Entropy",
            ylabel="Entropy",
            title=f"Experiment {experiment}: Entropy by treatment group",
            note="Single-timepoint experiment: the x-axis shows treatment groups, not weeks. Bars are entropy estimates; error bars are 95% bootstrap CIs when available.",
            treatment_group=treatment_group,
            comparator_group=comparator_group,
            low_col="Entropy_ci_low",
            high_col="Entropy_ci_high",
        )
        return
    weeks, groups = _get_weeks_and_groups(sub, treatment_group, comparator_group)
    cmap = _ensure_group_color_map(groups, treatment_group, comparator_group)
    comp = _comparison_lookup(comparison_df)
    
    fig, ax = plt.subplots(figsize=(10.2, 6.6))
    legend_handles, legend_labels = [], []
    for idx, group in enumerate(groups):
        g = sub[sub["Group"].astype(str) == group].sort_values("Week")
        if g.empty:
            continue
        color = cmap.get(group, _FALLBACK_COLORS[idx % len(_FALLBACK_COLORS)])
        y = g["Entropy"].to_numpy(dtype=float)
        if {"Entropy_ci_low", "Entropy_ci_high"}.issubset(g.columns):
            yerr = np.vstack([
                y - g["Entropy_ci_low"].to_numpy(dtype=float),
                g["Entropy_ci_high"].to_numpy(dtype=float) - y,
            ])
        else:
            yerr = None
        ax.errorbar(g["Week"], y, yerr=yerr, marker="o", linewidth=2.0, capsize=4, color=color)
        legend_handles.append(_legacy._Line2D([0], [0], color=color, marker="o", linewidth=2.0))
        legend_labels.append(display_group_name(group, idx))
        # Fixed pixel offset above/below by group so labels never overlap.
        for _, row in g.iterrows():
            ax.annotate(f"{row['Entropy']:.2f}", (float(row["Week"]), float(row["Entropy"])),
                        xytext=(0, 8 if idx == 0 else -8), textcoords="offset points",
                        ha="center", va="bottom" if idx == 0 else "top", fontsize=8.5, color=color,
                        zorder=7)
    
    # Overlay each subject's per-week entropy contribution as a faint cloud.
    pv = _overlay_trajectory_points(ax, experiment, "Entropy", weeks, groups, cmap)
    ax.set_xticks(weeks)
    ax.set_xticklabels([f"{w}W" for w in weeks])
    ax.set_xlabel("Week")
    ax.set_ylabel("Entropy")
    _legacy._format_axis(ax)

    if not comp.empty and treatment_group and comparator_group and len(groups) >= 2:
        y_min = min(float(sub["Entropy"].min()), min(pv) if pv else float(sub["Entropy"].min()))
        y_max = max(float(sub["Entropy"].max()), max(pv) if pv else float(sub["Entropy"].max()))
        y_span = max(y_max - y_min, 1.0)
        ann_y = y_max + 0.05 * y_span
        for i, week in enumerate(weeks):
            c = comp[comp["Week"] == week]
            if c.empty:
                continue
            row = c.iloc[0]
            delta = row.get("Entropy_diff_treat_minus_comp", row.get("Entropy_diff", np.nan))
            p = row.get("Entropy_p_boot", row.get("Entropy_p_mannwhitney", np.nan))
            txt = f"{_signif_label(p)}\nΔ={delta:+.2f}\np={p:.2g}" if pd.notna(delta) and pd.notna(p) else ""
            if txt:
                ax.text(week, ann_y, txt, ha="center", va="bottom", fontsize=9,
                        bbox=dict(boxstyle="round,pad=0.22", facecolor="white", edgecolor="#c0c0c0", alpha=0.88))
        ax.set_ylim(y_min - 0.06 * y_span, y_max + 0.30 * y_span)

    _legacy._add_top_title_and_note(
        fig,
        title=f"Experiment {experiment}: Entropy trajectory",
        note="Points are entropy estimates; error bars are 95% bootstrap CIs. X-axis follows uploaded weeks (e.g. 0W, 6W, 12W, 18W, 24W).",
        legend_handles=legend_handles,
        legend_labels=legend_labels,
        legend_ncol=min(len(legend_labels), 4),
        legend_y=0.885,
        title_y=0.985,
        note_y=0.948,
        note_width=120,
        note_fontsize=8.8,
    )
    fig.tight_layout(rect=[0.05, 0.06, 0.98, 0.80])
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_entropy_baseline_aligned(summary_df: pd.DataFrame, experiment: str, out_path: Path,
                                  treatment_group: Optional[str] = None,
                                  comparator_group: Optional[str] = None) -> None:
    sub = summary_df.copy()
    if "Experiment" in sub.columns:
        sub = sub[sub["Experiment"].astype(str) == str(experiment)].copy()
    sub = sub[[c for c in ["Week", "Group", "Entropy"] if c in sub.columns]].dropna(subset=["Entropy"]).copy()
    if sub.empty:
        return
    sub["Week"] = pd.to_numeric(sub["Week"], errors="coerce").astype(int)
    if _is_single_timepoint(sub):
        _plot_single_metric_by_group(
            sub, experiment, out_path,
            metric_col="Entropy",
            ylabel="Entropy",
            title=f"Experiment {experiment}: Entropy by treatment group",
            note="Single-timepoint experiment: baseline-aligned change is not applicable, so groups are compared directly on the x-axis.",
            treatment_group=treatment_group,
            comparator_group=comparator_group,
        )
        return
    weeks, groups = _get_weeks_and_groups(sub, treatment_group, comparator_group)
    cmap = _ensure_group_color_map(groups, treatment_group, comparator_group)
    
    rows = []
    for g in groups:
        gg = sub[sub["Group"].astype(str) == g].sort_values("Week").copy()
        if gg.empty:
            continue
        baseline_week = 0 if 0 in gg["Week"].values else int(gg["Week"].min())
        baseline = float(gg.loc[gg["Week"] == baseline_week, "Entropy"].iloc[0])
        gg["DeltaEntropy"] = gg["Entropy"] - baseline
        rows.append(gg)
    if not rows:
        return
    al = pd.concat(rows, ignore_index=True)
    
    fig, ax = plt.subplots(figsize=(10.2, 6.2))
    legend_handles, legend_labels = [], []
    for idx, group in enumerate(groups):
        g = al[al["Group"].astype(str) == group].sort_values("Week")
        color = cmap.get(group, _FALLBACK_COLORS[idx % len(_FALLBACK_COLORS)])
        ax.plot(g["Week"], g["DeltaEntropy"], marker="o", linewidth=2.0, color=color)
        legend_handles.append(_legacy._Line2D([0], [0], color=color, marker="o", linewidth=2.0))
        legend_labels.append(display_group_name(group, idx))
        for _, row in g.iterrows():
            ax.text(float(row["Week"]), float(row["DeltaEntropy"]) + 0.06, f"{row['DeltaEntropy']:+.2f}",
                    ha="center", va="bottom", fontsize=8.4, color=color)
    ax.axhline(0, color="#666666", linestyle="--", linewidth=1.0)
    ax.set_xticks(weeks)
    ax.set_xticklabels([f"{w}W" for w in weeks])
    ax.set_xlabel("Week")
    ax.set_ylabel("ΔEntropy vs baseline")
    _legacy._format_axis(ax)
    _legacy._add_top_title_and_note(
        fig,
        title=f"Experiment {experiment}: Baseline-aligned entropy change",
        note="Each group is aligned to its own baseline (0W if available; otherwise the earliest uploaded week). This answers the request of placing 0W at the same starting point.",
        legend_handles=legend_handles,
        legend_labels=legend_labels,
        legend_ncol=min(len(legend_labels), 4),
        legend_y=0.885,
        title_y=0.985,
        note_y=0.948,
        note_width=120,
        note_fontsize=8.7,
    )
    fig.tight_layout(rect=[0.05, 0.06, 0.98, 0.82])
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_true_normal_nmd_trajectory(summary_df: pd.DataFrame, experiment: str, out_path: Path,
                                    comparison_df: Optional[pd.DataFrame] = None,
                                    treatment_group: Optional[str] = None,
                                    comparator_group: Optional[str] = None) -> None:
    sub = summary_df.copy()
    if "Experiment" in sub.columns:
        sub = sub[sub["Experiment"].astype(str) == str(experiment)].copy()
    keep = [c for c in ["Week", "Group", "NMD_mean", "NMD_ci_low", "NMD_ci_high"] if c in sub.columns]
    sub = sub[keep].dropna(subset=["NMD_mean"]).copy()
    if sub.empty:
        return
    sub["Week"] = pd.to_numeric(sub["Week"], errors="coerce").astype(int)
    if _is_single_timepoint(sub):
        _plot_single_metric_by_group(
            sub, experiment, out_path,
            metric_col="NMD_mean",
            ylabel="True-Normal robust distance (NMD)",
            title=f"Experiment {experiment}: NMD by treatment group",
            note="Single-timepoint experiment: the x-axis shows treatment groups, not weeks. Lower NMD indicates closer to the Normal reference cloud.",
            treatment_group=treatment_group,
            comparator_group=comparator_group,
            low_col="NMD_ci_low",
            high_col="NMD_ci_high",
        )
        return
    weeks, groups = _get_weeks_and_groups(sub, treatment_group, comparator_group)
    cmap = _ensure_group_color_map(groups, treatment_group, comparator_group)
    comp = _comparison_lookup(comparison_df)
    fig, ax = plt.subplots(figsize=(10.2, 6.6))
    legend_handles, legend_labels = [], []
    for idx, group in enumerate(groups):
        g = sub[sub["Group"].astype(str) == group].sort_values("Week")
        color = cmap.get(group, _FALLBACK_COLORS[idx % len(_FALLBACK_COLORS)])
        y = g["NMD_mean"].to_numpy(dtype=float)
        if {"NMD_ci_low", "NMD_ci_high"}.issubset(g.columns):
            yerr = np.vstack([y - g["NMD_ci_low"].to_numpy(dtype=float), g["NMD_ci_high"].to_numpy(dtype=float) - y])
        else:
            yerr = None
        ax.errorbar(g["Week"], y, yerr=yerr, marker="o", linewidth=2.0, capsize=4, color=color)
        legend_handles.append(_legacy._Line2D([0], [0], color=color, marker="o", linewidth=2.0))
        legend_labels.append(display_group_name(group, idx))
        # Value labels: alternate above/below by group using a fixed pixel offset
        # so the two groups' labels never overlap even when the means coincide.
        for _, row in g.iterrows():
            ax.annotate(f"{row['NMD_mean']:.2f}", (float(row["Week"]), float(row["NMD_mean"])),
                        xytext=(0, 8 if idx == 0 else -8), textcoords="offset points",
                        ha="center", va="bottom" if idx == 0 else "top", fontsize=8.5, color=color,
                        zorder=7)
    # Overlay each subject's per-week NMD as a faint jittered cloud behind the
    # group lines (purely additive; the lines/labels are untouched).
    pv = _overlay_trajectory_points(ax, experiment, "NMD_mean", weeks, groups, cmap)
    ax.set_xticks(weeks)
    ax.set_xticklabels([f"{w}W" for w in weeks])
    ax.set_xlabel("Week")
    ax.set_ylabel("True-Normal robust distance (NMD)")
    _legacy._format_axis(ax)

    if not comp.empty and treatment_group and comparator_group and len(groups) >= 2:
        y_min = min(float(sub["NMD_mean"].min()), min(pv) if pv else float(sub["NMD_mean"].min()))
        y_max = max(float(sub["NMD_mean"].max()), max(pv) if pv else float(sub["NMD_mean"].max()))
        y_span = max(y_max - y_min, 1.0)
        ann_y = y_max + 0.05 * y_span  # annotation row above the point cloud
        for week in weeks:
            c = comp[comp["Week"] == week]
            if c.empty:
                continue
            row = c.iloc[0]
            delta = row.get("NMD_diff_treat_minus_comp", row.get("NMD_diff", np.nan))
            p = row.get("NMD_p_mannwhitney", row.get("NMD_p_boot", np.nan))
            txt = f"{_signif_label(p)}\nΔ={delta:+.2f}\np={p:.2g}" if pd.notna(delta) and pd.notna(p) else ""
            if txt:
                ax.text(week, ann_y, txt, ha="center", va="bottom", fontsize=9,
                        bbox=dict(boxstyle="round,pad=0.22", facecolor="white", edgecolor="#c0c0c0", alpha=0.88))
        ax.set_ylim(y_min - 0.06 * y_span, y_max + 0.30 * y_span)

    _legacy._add_top_title_and_note(
        fig,
        title=f"Experiment {experiment}: Disease-burden trajectory vs true Normal",
        note="NMD is the robust capped Normal-reference distance; raw covariance distance is retained as NMD_covariance_raw. Lower values indicate samples closer to the Normal cloud. X-axis follows uploaded weeks.",
        legend_handles=legend_handles,
        legend_labels=legend_labels,
        legend_ncol=min(len(legend_labels), 4),
        legend_y=0.885,
        title_y=0.985,
        note_y=0.948,
        note_width=118,
        note_fontsize=8.5,
    )
    fig.tight_layout(rect=[0.05, 0.06, 0.98, 0.80])
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_true_normal_nmd_baseline_aligned(summary_df: pd.DataFrame, experiment: str, out_path: Path,
                                          treatment_group: Optional[str] = None,
                                          comparator_group: Optional[str] = None) -> None:
    sub = summary_df.copy()
    if "Experiment" in sub.columns:
        sub = sub[sub["Experiment"].astype(str) == str(experiment)].copy()
    sub = sub[[c for c in ["Week", "Group", "NMD_mean"] if c in sub.columns]].dropna(subset=["NMD_mean"]).copy()
    if sub.empty:
        return
    sub["Week"] = pd.to_numeric(sub["Week"], errors="coerce").astype(int)
    if _is_single_timepoint(sub):
        _plot_single_metric_by_group(
            sub, experiment, out_path,
            metric_col="NMD_mean",
            ylabel="True-Normal robust distance (NMD)",
            title=f"Experiment {experiment}: NMD by treatment group",
            note="Single-timepoint experiment: baseline-aligned change is not applicable, so groups are compared directly on the x-axis.",
            treatment_group=treatment_group,
            comparator_group=comparator_group,
        )
        return
    weeks, groups = _get_weeks_and_groups(sub, treatment_group, comparator_group)
    cmap = _ensure_group_color_map(groups, treatment_group, comparator_group)
    rows = []
    for g in groups:
        gg = sub[sub["Group"].astype(str) == g].sort_values("Week").copy()
        baseline_week = 0 if 0 in gg["Week"].values else int(gg["Week"].min())
        baseline = float(gg.loc[gg["Week"] == baseline_week, "NMD_mean"].iloc[0])
        gg["DeltaNMD"] = gg["NMD_mean"] - baseline
        rows.append(gg)
    if not rows:
        return
    al = pd.concat(rows, ignore_index=True)
    fig, ax = plt.subplots(figsize=(10.2, 6.2))
    legend_handles, legend_labels = [], []
    for idx, group in enumerate(groups):
        g = al[al["Group"].astype(str) == group].sort_values("Week")
        color = cmap.get(group, _FALLBACK_COLORS[idx % len(_FALLBACK_COLORS)])
        ax.plot(g["Week"], g["DeltaNMD"], marker="o", linewidth=2.0, color=color)
        legend_handles.append(_legacy._Line2D([0], [0], color=color, marker="o", linewidth=2.0))
        legend_labels.append(display_group_name(group, idx))
        for _, row in g.iterrows():
            ax.text(float(row["Week"]), float(row["DeltaNMD"]) + 0.08, f"{row['DeltaNMD']:+.2f}",
                    ha="center", va="bottom", fontsize=8.4, color=color)
    ax.axhline(0, color="#666666", linestyle="--", linewidth=1.0)
    ax.set_xticks(weeks)
    ax.set_xticklabels([f"{w}W" for w in weeks])
    ax.set_xlabel("Week")
    ax.set_ylabel("ΔNMD vs baseline")
    _legacy._format_axis(ax)
    _legacy._add_top_title_and_note(
        fig,
        title=f"Experiment {experiment}: Baseline-aligned NMD change",
        note="Each group is aligned to its own baseline (0W if available; otherwise the earliest uploaded week). Values below 0 indicate reduced disease burden relative to baseline.",
        legend_handles=legend_handles,
        legend_labels=legend_labels,
        legend_ncol=min(len(legend_labels), 4),
        legend_y=0.885,
        title_y=0.985,
        note_y=0.948,
        note_width=118,
        note_fontsize=8.5,
    )
    fig.tight_layout(rect=[0.05, 0.06, 0.98, 0.82])
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_true_normal_nmd_entropy_3d(
    summary_df: pd.DataFrame,
    experiment: str,
    out_path: Path,
    treatment_group: Optional[str] = None,
    comparator_group: Optional[str] = None,
) -> None:
    sub = summary_df.copy()
    if "Experiment" in sub.columns:
        sub = sub[sub["Experiment"].astype(str) == str(experiment)].copy()
    sub = sub[[c for c in ["Week", "Group", "NMD_mean", "Entropy"] if c in sub.columns]].dropna(subset=["NMD_mean", "Entropy"]).copy()
    if sub.empty:
        return

    sub["Week"] = pd.to_numeric(sub["Week"], errors="coerce").astype(int)
    weeks = sorted(sub["Week"].dropna().astype(int).unique().tolist())
    groups = _ordered_groups_for_plot(sub, treatment_group, comparator_group)
    cmap = _ensure_group_color_map(groups, treatment_group, comparator_group)
    
    # normalize within metric for joint 3D display
    sub["NMD_display"] = _normalize_0_100(sub["NMD_mean"])
    sub["Entropy_display"] = _normalize_0_100(sub["Entropy"])

    fig = plt.figure(figsize=(12.4, 7.6))
    ax = fig.add_subplot(111, projection="3d")
    
    # Each week occupies a block; within each block, groups are separated, and for each group two bars are shown (NMD / Entropy)
    week_positions = {w: i * (len(groups) + 1.1) for i, w in enumerate(weeks)}
    dx = 0.55
    dy = 0.38
    metric_y = {"NMD": 0.0, "Entropy": 0.65}
    show_labels = len(weeks) * len(groups) <= 20
    max_height = 0
    
    for g_idx, group in enumerate(groups):
        base_color = cmap.get(group, _FALLBACK_COLORS[g_idx % len(_FALLBACK_COLORS)])
        ent_color = _lighten(base_color, 0.45)
        gsub = sub[sub["Group"].astype(str) == group].sort_values("Week")
        for _, row in gsub.iterrows():
            xbase = week_positions[int(row["Week"])] + g_idx * 0.72
            nmd_h = float(row["NMD_display"])
            ent_h = float(row["Entropy_display"])
            ax.bar3d(xbase, metric_y["NMD"], 0, dx, dy, nmd_h, color=base_color, shade=True, alpha=0.95, edgecolor="white", linewidth=0.35)
            ax.bar3d(xbase, metric_y["Entropy"], 0, dx, dy, ent_h, color=ent_color, shade=True, alpha=0.97, edgecolor="white", linewidth=0.35)
            if show_labels:
                ax.text(xbase + dx/2, metric_y["NMD"] + dy/2, nmd_h + 2, f"{row['NMD_mean']:.2f}", ha="center", va="bottom", fontsize=6.5)
                ax.text(xbase + dx/2, metric_y["Entropy"] + dy/2, ent_h + 2, f"{row['Entropy']:.2f}", ha="center", va="bottom", fontsize=6.5)
            max_height = max(max_height, nmd_h, ent_h)
    
    xticks = [week_positions[w] + (len(groups)-1)*0.36 + dx/2 for w in weeks]
    ax.set_xticks(xticks)
    ax.set_xticklabels([f"{w}W" for w in weeks], fontsize=9)
    ax.set_yticks([metric_y["NMD"] + dy/2, metric_y["Entropy"] + dy/2])
    ax.set_yticklabels(["NMD", "Entropy"], fontsize=10)
    ax.set_xlim(-0.5, max(xticks)+1.0)
    ax.set_ylim(-0.1, 1.2)
    ax.set_zlim(0, max(105, max_height + 10))
    ax.set_xlabel("Week", labelpad=8, fontsize=10)
    ax.set_ylabel("Metric", labelpad=10, fontsize=10)
    ax.set_zlabel("Normalized display height (0–100 within metric)", labelpad=10, fontsize=10)
    ax.view_init(elev=24, azim=-58)
    ax.grid(False)
    ax.xaxis.pane.set_facecolor((0.97, 0.97, 0.97, 1.0))
    ax.yaxis.pane.set_facecolor((0.97, 0.97, 0.97, 1.0))
    ax.zaxis.pane.set_facecolor((1.0, 1.0, 1.0, 0.0))
    
    handles, labels = [], []
    for g_idx, group in enumerate(groups):
        handles.append(plt.Line2D([0], [0], marker="s", color="w", markerfacecolor=cmap.get(group, _FALLBACK_COLORS[g_idx % len(_FALLBACK_COLORS)]), markersize=8.5))
        labels.append(display_group_name(group, g_idx))
    note = (
        "3D grouped bars are kept as a supplementary display. Within each week, groups are separated and NMD / Entropy are shown as two depth levels. "
        "Heights are normalized within NMD and Entropy separately (0–100) to allow joint display. For formal interpretation, prefer the 2D joint trajectory, the dual heatmap, or the grouped bar chart."
    )
    _legacy._add_top_title_and_note(
        fig,
        title=f"Experiment {experiment}: NMD × Entropy 3D supplementary chart",
        note=note,
        legend_handles=handles,
        legend_labels=labels,
        legend_ncol=min(len(labels), 4),
        legend_y=0.885,
        title_y=0.985,
        note_y=0.948,
        note_width=118,
        note_fontsize=8.5,
    )
    fig.tight_layout(rect=[0.02, 0.04, 0.98, 0.81])
    fig.savefig(out_path, dpi=260, bbox_inches="tight")
    plt.close(fig)


def plot_true_normal_nmd_entropy_grouped_bar(
    summary_df: pd.DataFrame,
    experiment: str,
    out_path: Path,
    treatment_group: Optional[str] = None,
    comparator_group: Optional[str] = None,
) -> None:
    sub = summary_df.copy()
    if "Experiment" in sub.columns:
        sub = sub[sub["Experiment"].astype(str) == str(experiment)].copy()
    sub = sub[[c for c in ["Week", "Group", "NMD_mean", "Entropy"] if c in sub.columns]].dropna(subset=["NMD_mean", "Entropy"]).copy()
    if sub.empty:
        return
    sub["Week"] = pd.to_numeric(sub["Week"], errors="coerce").astype(int)
    weeks, groups = _get_weeks_and_groups(sub, treatment_group, comparator_group)
    cmap = _ensure_group_color_map(groups, treatment_group, comparator_group)
    
    # Normalize within metric to 0-100 to allow grouped display similar to Fig4
    sub["NMD_display"] = _normalize_0_100(sub["NMD_mean"])
    sub["Entropy_display"] = _normalize_0_100(sub["Entropy"])

    n = len(groups)
    fig, axes = plt.subplots(1, n, figsize=(max(7.0, 4.0 * n), 4.8), squeeze=False)
    axes = axes.flatten()
    width = 0.36
    x = np.arange(len(weeks))
    
    for idx, group in enumerate(groups):
        ax = axes[idx]
        g = sub[sub["Group"].astype(str) == group].sort_values("Week")
        if g.empty:
            ax.axis('off')
            continue
        nmd = g["NMD_display"].to_numpy(dtype=float)
        ent = g["Entropy_display"].to_numpy(dtype=float)
        base_color = cmap.get(group, _FALLBACK_COLORS[idx % len(_FALLBACK_COLORS)])
        ent_color = _lighten(base_color, 0.35)
        b1 = ax.bar(x - width/2, nmd, width=width, color=base_color, alpha=0.88, label="Disease burden (NMD)")
        b2 = ax.bar(x + width/2, ent, width=width, color=ent_color, alpha=0.96, label="Entropy")
        ax.set_title(display_group_name(group, idx), fontsize=11)
        ax.set_xticks(x)
        ax.set_xticklabels([f"{w}W" for w in weeks])
        ax.set_ylim(0, 105)
        _remove_horizontal_grid(ax)
        if idx == 0:
            ax.set_ylabel("Normalized group mean (0–100)")
        for bi, raw in zip(b1, g["NMD_mean"].tolist()):
            ax.text(bi.get_x() + bi.get_width()/2, bi.get_height()+2, f"{raw:.2f}", ha="center", va="bottom", fontsize=7)
        for bi, raw in zip(b2, g["Entropy"].tolist()):
            ax.text(bi.get_x() + bi.get_width()/2, bi.get_height()+2, f"{raw:.2f}", ha="center", va="bottom", fontsize=7)
    
    handles = [plt.Rectangle((0,0),1,1,color="#4e79a7"), plt.Rectangle((0,0),1,1,color="#9ecae1")]
    labels = ["Disease burden (NMD)", "Entropy"]
    note = (
        "This grouped bar chart is analogous to the example 'group-wise disease score and bio-entropy' figure. "
        "Within each group panel, bar heights are normalized within metric (0–100) to allow NMD and Entropy to be shown together; labels above bars are raw values. "
        "This avoids the occlusion of 3D bars and is recommended for intuitive side-by-side comparison."
    )
    _legacy._add_top_title_and_note(
        fig,
        title=f"Experiment {experiment}: Group-wise disease burden and entropy",
        note=note,
        legend_handles=handles,
        legend_labels=labels,
        legend_ncol=2,
        legend_y=0.885,
        title_y=0.985,
        note_y=0.948,
        note_width=118,
        note_fontsize=8.5,
    )
    fig.tight_layout(rect=[0.05, 0.06, 0.98, 0.81])
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_true_normal_nmd_entropy_joint_trajectory(
    summary_df: pd.DataFrame,
    experiment: str,
    out_path: Path,
    treatment_group: Optional[str] = None,
    comparator_group: Optional[str] = None,
) -> None:
    sub = summary_df.copy()
    if "Experiment" in sub.columns:
        sub = sub[sub["Experiment"].astype(str) == str(experiment)].copy()
    sub = sub[[c for c in ["Week", "Group", "NMD_mean", "Entropy"] if c in sub.columns]].dropna(subset=["NMD_mean", "Entropy"]).copy()
    if sub.empty:
        return
    sub["Week"] = pd.to_numeric(sub["Week"], errors="coerce").astype(int)
    groups = _ordered_groups_for_plot(sub, treatment_group, comparator_group)
    cmap = _ensure_group_color_map(groups, treatment_group, comparator_group)
    single_timepoint = _is_single_timepoint(sub)

    # Use faceting when group count is large to avoid overlap.
    if len(groups) <= 4:
        fig, ax = plt.subplots(figsize=(8.0, 6.2))
        legend_handles, legend_labels = [], []
        text_items = []
        for idx, group in enumerate(groups):
            g = sub[sub["Group"].astype(str) == group].sort_values("Week")
            if g.empty:
                continue
            color = cmap.get(group, _FALLBACK_COLORS[idx % len(_FALLBACK_COLORS)])
            ax.plot(g["NMD_mean"], g["Entropy"], marker="o", linewidth=2.2, color=color, zorder=3)
            # arrows along time
            xvals = g["NMD_mean"].to_numpy(dtype=float)
            yvals = g["Entropy"].to_numpy(dtype=float)
            for i in range(len(g) - 1):
                ax.annotate("", xy=(xvals[i+1], yvals[i+1]), xytext=(xvals[i], yvals[i]),
                            arrowprops=dict(arrowstyle="->", lw=1.0, color=color, alpha=0.70))
            legend_handles.append(_legacy._Line2D([0], [0], color=color, marker="o", linewidth=2.2))
            legend_labels.append(display_group_name(group, idx))
            for _, row in g.iterrows():
                label = display_group_name(group, idx) if single_timepoint else f"{int(row['Week'])}W"
                text_items.append((float(row["NMD_mean"]), float(row["Entropy"]), label, color))
        _legacy._annotate_points_nonoverlap(ax, text_items, fontsize=8.0)
        ax.set_xlabel("True-Normal robust distance (NMD)")
        ax.set_ylabel("Entropy")
        _legacy._format_axis(ax)
        note = (
            "Single-timepoint experiment: each point is one treatment group; there is no week trajectory. "
            "Lower-left indicates lower NMD and lower entropy."
            if single_timepoint else
            "Arrows show the time direction. Movement toward the lower-left indicates that samples become closer to the "
            "true Normal reference cloud while system entropy also decreases. This is the recommended primary figure "
            "for joint interpretation of NMD and Entropy."
        )
        _legacy._add_top_title_and_note(
            fig,
            title=f"Experiment {experiment}: NMD–Entropy group comparison" if single_timepoint else f"Experiment {experiment}: NMD–Entropy joint trajectory",
            note=note,
            legend_handles=legend_handles,
            legend_labels=legend_labels,
            legend_ncol=min(len(legend_labels), 3),
            legend_y=0.885,
            title_y=0.985,
            note_y=0.948,
            note_width=118,
            note_fontsize=8.7,
        )
        fig.tight_layout(rect=[0.05, 0.06, 0.98, 0.80])
        fig.savefig(out_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        return

    # Faceted layout for many groups
    n = len(groups)
    ncols = 2 if n <= 6 else 3
    nrows = int(math.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.2 * ncols, 3.9 * nrows), squeeze=False)
    axes_flat = axes.flatten()
    x_min = float(sub["NMD_mean"].min())
    x_max = float(sub["NMD_mean"].max())
    y_min = float(sub["Entropy"].min())
    y_max = float(sub["Entropy"].max())
    x_pad = max((x_max - x_min) * 0.08, 0.2)
    y_pad = max((y_max - y_min) * 0.08, 0.15)
    for idx, group in enumerate(groups):
        ax = axes_flat[idx]
        g = sub[sub["Group"].astype(str) == group].sort_values("Week")
        color = cmap.get(group, _FALLBACK_COLORS[idx % len(_FALLBACK_COLORS)])
        ax.plot(g["NMD_mean"], g["Entropy"], marker="o", linewidth=2.1, color=color, zorder=3)
        xvals = g["NMD_mean"].to_numpy(dtype=float)
        yvals = g["Entropy"].to_numpy(dtype=float)
        for i in range(len(g) - 1):
            ax.annotate("", xy=(xvals[i+1], yvals[i+1]), xytext=(xvals[i], yvals[i]),
                        arrowprops=dict(arrowstyle="->", lw=1.0, color=color, alpha=0.70))
        items = [
            (float(r["NMD_mean"]), float(r["Entropy"]), display_group_name(group, idx) if single_timepoint else f"{int(r['Week'])}W", color)
            for _, r in g.iterrows()
        ]
        _legacy._annotate_points_nonoverlap(ax, items, fontsize=7.6)
        ax.set_title(display_group_name(group, idx), fontsize=11.0, pad=8)
        ax.set_xlim(x_min - x_pad, x_max + x_pad)
        ax.set_ylim(y_min - y_pad, y_max + y_pad)
        _legacy._format_axis(ax)
        _remove_horizontal_grid(ax)
        if idx % ncols == 0:
            ax.set_ylabel("Entropy")
        else:
            ax.set_ylabel("")
        if idx >= (nrows - 1) * ncols:
            ax.set_xlabel("NMD")
        else:
            ax.set_xlabel("")
    for j in range(len(groups), len(axes_flat)):
        axes_flat[j].axis("off")
    note = (
        "Single-timepoint experiment: each panel corresponds to one treatment group; there is no week trajectory. "
        "The lower-left direction represents simultaneous lower NMD and lower entropy."
        if single_timepoint else
        "Faceted layout is used automatically when the number of experimental groups is large. "
        "Each panel keeps the same axes so the lower-left direction still represents simultaneous reduction "
        "of disease burden and entropy."
    )
    _legacy._add_top_title_and_note(
        fig,
        title=f"Experiment {experiment}: NMD–Entropy group comparison (faceted)" if single_timepoint else f"Experiment {experiment}: NMD–Entropy joint trajectory (faceted)",
        note=note,
        legend_handles=None,
        legend_labels=None,
        title_y=0.985,
        note_y=0.948,
        note_width=122,
        note_fontsize=8.6,
    )
    fig.tight_layout(rect=[0.04, 0.05, 0.98, 0.87])
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_true_normal_nmd_entropy_dual_heatmap(
    summary_df: pd.DataFrame,
    experiment: str,
    out_path: Path,
    treatment_group: Optional[str] = None,
    comparator_group: Optional[str] = None,
) -> None:
    sub = summary_df.copy()
    if "Experiment" in sub.columns:
        sub = sub[sub["Experiment"].astype(str) == str(experiment)].copy()
    sub = sub[[c for c in ["Week", "Group", "NMD_mean", "Entropy"] if c in sub.columns]].dropna(subset=["NMD_mean", "Entropy"]).copy()
    if sub.empty:
        return
    sub["Week"] = pd.to_numeric(sub["Week"], errors="coerce").astype(int)
    groups = _ordered_groups_for_plot(sub, treatment_group, comparator_group)
    weeks = sorted(sub["Week"].dropna().astype(int).unique().tolist())
    if _plot_burden_entropy_single_group_heatmap(
        sub, experiment, out_path, "NMD_mean", _burden_metadata("NMD"), treatment_group, comparator_group
    ):
        return

    nmd_wide = sub.pivot(index="Group", columns="Week", values="NMD_mean").reindex(index=groups, columns=weeks)
    ent_wide = sub.pivot(index="Group", columns="Week", values="Entropy").reindex(index=groups, columns=weeks)
    # robust dataframe-wide normalization
    nmd_vals = nmd_wide.to_numpy(dtype=float)
    ent_vals = ent_wide.to_numpy(dtype=float)
    def _norm_arr(arr):
        finite = arr[np.isfinite(arr)]
        if finite.size == 0:
            return np.full_like(arr, np.nan, dtype=float)
        vmin = float(finite.min()); vmax = float(finite.max())
        if math.isclose(vmin, vmax):
            return np.full_like(arr, 50.0, dtype=float)
        return (arr - vmin) / (vmax - vmin) * 100.0
    nmd_norm_vals = _norm_arr(nmd_vals)
    ent_norm_vals = _norm_arr(ent_vals)

    fig, axes = plt.subplots(2, 1, figsize=(1.22 * max(6, len(weeks)) + 3.2, 1.0 * max(3, len(groups)) + 3.8), sharex=True)
    panels = [
        (axes[0], nmd_norm_vals, nmd_wide, "NMD burden (0–100 within experiment)"),
        (axes[1], ent_norm_vals, ent_wide, "Entropy burden (0–100 within experiment)"),
    ]
    group_labels = [display_group_name(g, i) for i, g in enumerate(groups)]
    for ax, zvals, raw_df, title in panels:
        im = ax.imshow(zvals, aspect="auto", cmap="RdYlBu_r", vmin=0, vmax=100)
        ax.set_title(title, fontsize=11.0, pad=8)
        ax.set_yticks(np.arange(len(groups)))
        ax.set_yticklabels(group_labels, fontsize=9)
        ax.set_xticks(np.arange(len(weeks)))
        ax.set_xticklabels([f"{w}W" for w in weeks], fontsize=9)
        for i in range(len(groups)):
            for j in range(len(weeks)):
                raw_val = raw_df.iloc[i, j]
                if pd.notna(raw_val) and len(groups) <= 6 and len(weeks) <= 9:
                    txt = f"{raw_val:.2f}"
                    ax.text(j, i, txt, ha="center", va="center", fontsize=7.2, color="#222222")
        cbar = fig.colorbar(im, ax=ax, fraction=0.024, pad=0.02)
        cbar.ax.tick_params(labelsize=8)

    axes[1].set_xlabel("Week")
    note = (
        "This figure avoids 3D occlusion and is recommended when many experimental groups are present. "
        "Both panels are normalized within the experiment to the same 0–100 display range; lower values indicate "
        "lower burden. If NMD and Entropy both cool down together across the same cells, the two measures are changing synergistically."
    )
    _legacy._add_top_title_and_note(
        fig,
        title=f"Experiment {experiment}: NMD / Entropy dual heatmap",
        note=note,
        legend_handles=None,
        legend_labels=None,
        title_y=0.985,
        note_y=0.948,
        note_width=118,
        note_fontsize=8.7,
    )
    fig.tight_layout(rect=[0.05, 0.05, 0.98, 0.84])
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _burden_metadata(burden_kind: str) -> Dict[str, str]:
    key = str(burden_kind).upper()
    if key == "NRBS":
        return {
            "mean": "NRBS_mean",
            "low": "NRBS_mean_ci_low",
            "high": "NRBS_mean_ci_high",
            "diff": "NRBS_diff_treat_minus_comp",
            "p": "NRBS_p_mannwhitney",
            "axis": "Normal-range burden score (NRBS)",
            "short": "NRBS",
            "prefix": "normal_range_nrbs",
            "joint_title": "NRBS–Entropy joint trajectory",
            "heat_title": "NRBS / Entropy dual heatmap",
            "grouped_title": "Group-wise disease burden and entropy",
            "3d_title": "NRBS × Entropy 3D grouped chart",
            "note": "NRBS is anchored to published adult Normal ranges (an absolute health scale). 0 means all markers within Normal range; higher indicates greater deviation from Normal.",
        }
    if key == "BRI":
        return {
            "mean": "BRI_mean",
            "low": "BRI_mean_ci_low",
            "high": "BRI_mean_ci_high",
            "diff": "BRI_diff_treat_minus_comp",
            "p": "BRI_p_mannwhitney",
            "axis": "Baseline-referenced burden index (BRI)",
            "short": "BRI",
            "prefix": "baseline_bri",
            "joint_title": "BRI–Entropy joint trajectory",
            "heat_title": "BRI / Entropy dual heatmap",
            "grouped_title": "Group-wise BRI and entropy",
            "3d_title": "BRI × Entropy 3D grouped chart",
            "note": "BRI is referenced to the baseline / 0W cohort. Lower values indicate lower disease burden relative to baseline.",
        }
    else:
        return {
            "mean": "NMD_mean",
            "low": "NMD_ci_low",
            "high": "NMD_ci_high",
            "diff": "NMD_diff",
            "p": "NMD_p_boot",
            "axis": "True-Normal robust distance (NMD)",
            "short": "NMD",
            "prefix": "true_normal_nmd",
            "joint_title": "NMD–Entropy joint trajectory",
            "heat_title": "NMD / Entropy dual heatmap",
            "grouped_title": "Group-wise disease burden and entropy",
            "3d_title": "NMD × Entropy 3D grouped chart",
            "note": "NMD is the robust capped Normal-reference distance; raw covariance distance is retained as NMD_covariance_raw. Lower values indicate samples closer to the Normal reference.",
        }


def plot_burden_entropy_joint_trajectory(
    summary_df: pd.DataFrame,
    experiment: str,
    out_path: Path,
    burden_kind: str = "BRI",
    treatment_group: Optional[str] = None,
    comparator_group: Optional[str] = None,
) -> None:
    meta = _burden_metadata(burden_kind)
    burden_col = meta["mean"]

    sub = summary_df.copy()
    if "Experiment" in sub.columns:
        sub = sub[sub["Experiment"].astype(str) == str(experiment)].copy()
    keep = [c for c in ["Week", "Group", burden_col, "Entropy"] if c in sub.columns]
    sub = sub[keep].dropna(subset=[burden_col, "Entropy"]).copy()
    if sub.empty:
        return

    sub["Week"] = pd.to_numeric(sub["Week"], errors="coerce").astype(int)
    groups = _ordered_groups_for_plot(sub, treatment_group, comparator_group)
    cmap = _ensure_group_color_map(groups, treatment_group, comparator_group)
    single_timepoint = _is_single_timepoint(sub)

    # Use faceting when group count is large to avoid overlap.
    if len(groups) <= 4:
        fig, ax = plt.subplots(figsize=(8.0, 6.2))
        legend_handles, legend_labels = [], []
        text_items = []
        for idx, group in enumerate(groups):
            g = sub[sub["Group"].astype(str) == group].sort_values("Week")
            if g.empty:
                continue
            color = cmap.get(group, _FALLBACK_COLORS[idx % len(_FALLBACK_COLORS)])
            ax.plot(g[burden_col], g["Entropy"], marker="o", linewidth=2.2, color=color, zorder=3)
            xvals = g[burden_col].to_numpy(dtype=float)
            yvals = g["Entropy"].to_numpy(dtype=float)
            for i in range(len(g) - 1):
                ax.annotate("", xy=(xvals[i+1], yvals[i+1]), xytext=(xvals[i], yvals[i]),
                            arrowprops=dict(arrowstyle="->", lw=1.0, color=color, alpha=0.72))
            legend_handles.append(_legacy._Line2D([0], [0], color=color, marker="o", linewidth=2.2))
            legend_labels.append(display_group_name(group, idx))
            for _, row in g.iterrows():
                label = display_group_name(group, idx) if single_timepoint else f"{int(row['Week'])}W"
                text_items.append((float(row[burden_col]), float(row["Entropy"]), label, color))
        _legacy._annotate_points_nonoverlap(ax, text_items, fontsize=8.0)
        ax.set_xlabel(meta["axis"])
        ax.set_ylabel("Entropy")
        _legacy._format_axis(ax)
        note = (
            f"Single-timepoint experiment: each point is one treatment group; there is no week trajectory. "
            f"Lower-left indicates lower {meta['short']} and lower Entropy."
            if single_timepoint else
            f"Arrows show the time direction. Movement toward the lower-left indicates that {meta['short']} and Entropy decrease together. "
            f"This figure is recommended for visually judging whether burden reduction and entropy reduction are positively correlated across time."
        )
        _legacy._add_top_title_and_note(
            fig,
            title=f"Experiment {experiment}: {meta['short']}–Entropy group comparison" if single_timepoint else f"Experiment {experiment}: {meta['joint_title']}",
            note=note,
            legend_handles=legend_handles,
            legend_labels=legend_labels,
            legend_ncol=min(len(legend_labels), 3),
            legend_y=0.885,
            title_y=0.985,
            note_y=0.948,
            note_width=118,
            note_fontsize=8.7,
        )
        fig.tight_layout(rect=[0.05, 0.06, 0.98, 0.80])
        fig.savefig(out_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        return

    n = len(groups)
    ncols = 2 if n <= 6 else 3
    nrows = int(math.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.2 * ncols, 3.9 * nrows), squeeze=False)
    axes_flat = axes.flatten()
    x_min = float(sub[burden_col].min())
    x_max = float(sub[burden_col].max())
    y_min = float(sub["Entropy"].min())
    y_max = float(sub["Entropy"].max())
    x_pad = max((x_max - x_min) * 0.08, 0.2)
    y_pad = max((y_max - y_min) * 0.08, 0.15)
    for idx, group in enumerate(groups):
        ax = axes_flat[idx]
        g = sub[sub["Group"].astype(str) == group].sort_values("Week")
        color = cmap.get(group, _FALLBACK_COLORS[idx % len(_FALLBACK_COLORS)])
        ax.plot(g[burden_col], g["Entropy"], marker="o", linewidth=2.1, color=color, zorder=3)
        xvals = g[burden_col].to_numpy(dtype=float)
        yvals = g["Entropy"].to_numpy(dtype=float)
        for i in range(len(g) - 1):
            ax.annotate("", xy=(xvals[i+1], yvals[i+1]), xytext=(xvals[i], yvals[i]),
                        arrowprops=dict(arrowstyle="->", lw=1.0, color=color, alpha=0.70))
        items = [
            (float(r[burden_col]), float(r["Entropy"]), display_group_name(group, idx) if single_timepoint else f"{int(r['Week'])}W", color)
            for _, r in g.iterrows()
        ]
        _legacy._annotate_points_nonoverlap(ax, items, fontsize=7.6)
        ax.set_title(display_group_name(group, idx), fontsize=11.0, pad=8)
        ax.set_xlim(x_min - x_pad, x_max + x_pad)
        ax.set_ylim(y_min - y_pad, y_max + y_pad)
        _legacy._format_axis(ax)
        _remove_horizontal_grid(ax)
        if idx % ncols == 0:
            ax.set_ylabel("Entropy")
        else:
            ax.set_ylabel("")
        if idx >= (nrows - 1) * ncols:
            ax.set_xlabel(meta["short"])
        else:
            ax.set_xlabel("")
    for j in range(len(groups), len(axes_flat)):
        axes_flat[j].axis("off")
    note = (
        f"Single-timepoint experiment: each panel corresponds to one treatment group; there is no week trajectory. "
        f"The lower-left direction represents simultaneous lower {meta['short']} and lower Entropy."
        if single_timepoint else
        "Faceted layout is used automatically when the number of experimental groups is large. "
        "Each panel keeps the same axes so the lower-left direction still represents simultaneous reduction "
        f"of {meta['short']} and Entropy."
    )
    _legacy._add_top_title_and_note(
        fig,
        title=f"Experiment {experiment}: {meta['short']}–Entropy group comparison (faceted)" if single_timepoint else f"Experiment {experiment}: {meta['joint_title']} (faceted)",
        note=note,
        legend_handles=None,
        legend_labels=None,
        title_y=0.985,
        note_y=0.948,
        note_width=122,
        note_fontsize=8.6,
    )
    fig.tight_layout(rect=[0.04, 0.05, 0.98, 0.87])
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _annotate_labels_avoid_overlap(ax, points: List[tuple[float, float, str, str]], fontsize: float = 8.0) -> None:
    if not points:
        return
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    xspan = max(xmax - xmin, 1e-9)
    yspan = max(ymax - ymin, 1e-9)
    y_gap = yspan * 0.045
    x_mid = xmin + xspan * 0.52
    grouped = {"left": [], "right": []}
    for x, y, label, color in points:
        side = "right" if x <= x_mid else "left"
        grouped[side].append((x, y, label, color))

    for side, items in grouped.items():
        if not items:
            continue
        items = sorted(items, key=lambda p: p[1])
        placed = []
        last_y = ymin - yspan
        for x, y, label, color in items:
            ty = max(float(y), last_y + y_gap)
            placed.append([x, y, label, color, ty])
            last_y = ty
        overflow = placed[-1][4] - (ymax - yspan * 0.04)
        if overflow > 0:
            for item in placed:
                item[4] -= overflow
        for x, y, label, color, ty in placed:
            if side == "right":
                tx = min(x + xspan * 0.045, xmax - xspan * 0.02)
                ha = "left"
            else:
                tx = max(x - xspan * 0.045, xmin + xspan * 0.02)
                ha = "right"
            ax.annotate(
                label,
                xy=(x, y),
                xytext=(tx, ty),
                textcoords="data",
                ha=ha,
                va="center",
                fontsize=fontsize,
                color=color,
                bbox=dict(boxstyle="round,pad=0.18", facecolor="white", edgecolor="none", alpha=0.82),
                arrowprops=dict(arrowstyle="-", color=color, lw=0.8, alpha=0.55),
                zorder=5,
            )


def plot_burden_entropy_comparison_arrow(
    summary_df: pd.DataFrame,
    experiment: str,
    out_path: Path,
    burden_kind: str = "NMD",
    treatment_group: Optional[str] = None,
    comparator_group: Optional[str] = None,
) -> None:
    meta = _burden_metadata(burden_kind)
    burden_col = meta["mean"]
    if burden_col not in summary_df.columns or "Entropy" not in summary_df.columns:
        return
    sub = summary_df.copy()
    if "Experiment" in sub.columns:
        sub = sub[sub["Experiment"].astype(str) == str(experiment)].copy()
    sub = sub[[c for c in ["Week", "Group", burden_col, "Entropy"] if c in sub.columns]].dropna(subset=[burden_col, "Entropy"]).copy()
    if sub.empty:
        return
    sub["Week"] = pd.to_numeric(sub["Week"], errors="coerce").astype(int)
    single_timepoint = _is_single_timepoint(sub)
    groups = _ordered_single_timepoint_groups(sub, burden_col) if single_timepoint else _ordered_groups_for_plot(sub, treatment_group, comparator_group)
    cmap = _ensure_group_color_map(groups, treatment_group, comparator_group)

    x_min, x_max = float(sub[burden_col].min()), float(sub[burden_col].max())
    y_min, y_max = float(sub["Entropy"].min()), float(sub["Entropy"].max())
    x_pad = max((x_max - x_min) * 0.22, 0.45)
    y_pad = max((y_max - y_min) * 0.22, 0.18)

    fig, ax = plt.subplots(figsize=(9.0, 6.8))
    ax.set_xlim(x_min - x_pad, x_max + x_pad)
    ax.set_ylim(y_min - y_pad, y_max + y_pad)

    # Direction guide: conceptual severity decreases from upper-right toward lower-left.
    start = (x_max + x_pad * 0.55, y_max + y_pad * 0.55)
    end = (x_min - x_pad * 0.35, y_min - y_pad * 0.35)
    ax.annotate(
        "",
        xy=end,
        xytext=start,
        arrowprops=dict(arrowstyle="-|>", lw=2.2, color="#6B7280", alpha=0.35, mutation_scale=18),
        zorder=1,
    )
    ax.text(start[0], start[1], "more severe", ha="right", va="bottom", fontsize=9.2, color="#6B7280")
    ax.text(end[0], end[1], "closer to Normal", ha="left", va="top", fontsize=9.2, color="#6B7280")

    legend_handles, legend_labels, label_points = [], [], []
    for idx, group in enumerate(groups):
        g = sub[sub["Group"].astype(str) == str(group)].sort_values("Week")
        if g.empty:
            continue
        color = cmap.get(group, _FALLBACK_COLORS[idx % len(_FALLBACK_COLORS)])
        ax.scatter(
            g[burden_col],
            g["Entropy"],
            s=92 if single_timepoint else 68,
            color=color,
            edgecolor="white",
            linewidth=1.0,
            zorder=4,
        )
        if not single_timepoint and len(g) > 1:
            ax.plot(g[burden_col], g["Entropy"], color=color, linewidth=1.5, alpha=0.66, zorder=3)
            xs = g[burden_col].to_numpy(dtype=float)
            ys = g["Entropy"].to_numpy(dtype=float)
            for i in range(len(g) - 1):
                ax.annotate(
                    "",
                    xy=(xs[i + 1], ys[i + 1]),
                    xytext=(xs[i], ys[i]),
                    arrowprops=dict(arrowstyle="->", lw=1.0, color=color, alpha=0.62, mutation_scale=11),
                    zorder=3,
                )
        legend_handles.append(_legacy._Line2D([0], [0], color=color, marker="o", linewidth=1.8))
        legend_labels.append(display_group_name(group, idx))
        if single_timepoint:
            for _, row in g.iterrows():
                label_points.append((float(row[burden_col]), float(row["Entropy"]), display_group_name(group, idx), color))
        else:
            first = g.iloc[0]
            last = g.iloc[-1]
            label_points.append((float(first[burden_col]), float(first["Entropy"]), f"{display_group_name(group, idx)} {int(first['Week'])}W", color))
            if int(last["Week"]) != int(first["Week"]):
                label_points.append((float(last[burden_col]), float(last["Entropy"]), f"{display_group_name(group, idx)} {int(last['Week'])}W", color))

    _annotate_labels_avoid_overlap(ax, label_points, fontsize=8.0 if len(label_points) <= 10 else 7.4)
    ax.set_xlabel(meta["axis"])
    ax.set_ylabel("Entropy")
    _legacy._format_axis(ax)
    _apply_academic_axis(ax, "")
    note = (
        f"Each point is one treatment group; the grey arrow marks the desired interpretation direction from higher {meta['short']} / higher entropy toward lower values."
        if single_timepoint else
        "Each point is a group-week state. Colored arrows connect time within each group; the grey arrow marks the conceptual direction from severe to Normal-like states."
    )
    _legacy._add_top_title_and_note(
        fig,
        title=f"Experiment {experiment}: {meta['short']}–Entropy state map",
        note=note,
        legend_handles=legend_handles,
        legend_labels=legend_labels,
        legend_ncol=min(4, max(1, len(legend_labels))),
        legend_y=0.875,
        title_y=0.985,
        note_y=0.948,
        note_width=124,
        note_fontsize=8.7,
    )
    fig.tight_layout(rect=[0.07, 0.07, 0.98, 0.79])
    fig.savefig(out_path, dpi=330, bbox_inches="tight")
    plt.close(fig)


_CLICKABLE_STATE_JS = r"""
(function() {
  var gd = document.getElementById(__DIV_ID__);
  if (!gd) return;
  var LABELS = __LABELS__;
  var CURVES = __CURVES__;
  var CONNECTOR = __CONNECTOR__;
  var LINES = __LINES__;
  var SINGLE = __SINGLE__;
  var selBox = document.getElementById(__SEL_ID__);
  var linesBox = document.getElementById(__LINES_ID__);
  var listBox = document.getElementById(__LIST_ID__);
  var toolbar = document.getElementById(__TOOLBAR_ID__);
  var MATERIALS = { flat:{type:'none',edge:null}, glossy:{type:'radial',edge:'#ffffff'},
                    metallic:{type:'radial',edge:'#5b6270'}, pearl:{type:'horizontal',edge:'#eef1f6'} };
  var MATERIAL_OPTIONS = [['flat','实心'],['glossy','3D 光泽'],['metallic','金属'],['pearl','柔光/珍珠']];
  var curveByNumber = {};
  CURVES.forEach(function(c){ curveByNumber[c.curveNumber] = c; });
  var selected = null;

  function cloneArray(v){ return Array.prototype.slice.call(v || []); }
  function isPointTrace(tr){ return tr && tr.meta && tr.meta.bioentropyRole === 'point'; }
  function nonNullXY(cn){
    var tr = gd.data[cn], out = [];
    if (!tr) return out;
    var xs = tr.x || [], ys = tr.y || [];
    for (var i = 0; i < xs.length; i++) {
      if (xs[i] != null && ys[i] != null) out.push([Number(xs[i]), Number(ys[i])]);
    }
    return out;
  }
  function firstXY(cn){ var p = nonNullXY(cn); return p.length ? p[0] : null; }

  // Re-thread the single-timepoint connector spline through the surviving points
  // after a delete/undo (the multi-week group traces re-spline on their own via
  // connectgaps). No arrows are drawn on the map.
  function rebuildConnector(){
    if (!CONNECTOR) return;
    var pts = [];
    CONNECTOR.members.forEach(function(m){ var p = firstXY(m); if (p) pts.push(p); });
    var cx = pts.map(function(p){ return p[0]; }), cy = pts.map(function(p){ return p[1]; });
    Plotly.restyle(gd, {x: [cx], y: [cy]}, [CONNECTOR.curveNumber]);
  }
  function updateTrajectory(cn){
    if (SINGLE) rebuildConnector();
  }

  function relayoutLabel(rec, fields){
    if (rec.annotationIndex == null) return;
    var u = {};
    Object.keys(fields).forEach(function(k){ u['annotations[' + rec.annotationIndex + '].' + k] = fields[k]; });
    Plotly.relayout(gd, u);
  }
  function renderLabel(rec){
    relayoutLabel(rec, {
      text: String(rec.text || ''), visible: !!rec.visible,
      xshift: Number(rec.xshift || 0), yshift: Number(rec.yshift || 0),
      'font.size': Number(rec.fontSize || 9)
    });
  }
  function setMarkerAttr(rec, attr, val){
    var tr = gd.data[rec.curveNumber];
    if (!tr || !tr.marker) return;
    var arr = cloneArray(tr.marker[attr]);
    while (arr.length <= rec.pointNumber) arr.push(null);
    arr[rec.pointNumber] = val;
    var u = {}; u['marker.' + attr] = [arr];
    Plotly.restyle(gd, u, [rec.curveNumber]);
  }
  function setMarkerNested(rec, path, val){
    var tr = gd.data[rec.curveNumber];
    if (!tr || !tr.marker) return;
    var parts = path.split('.'), obj = tr.marker;
    for (var k = 0; k < parts.length - 1 && obj; k++) obj = obj[parts[k]];
    var arr = obj ? cloneArray(obj[parts[parts.length - 1]]) : [];
    while (arr.length <= rec.pointNumber) arr.push(null);
    arr[rec.pointNumber] = val;
    var u = {}; u['marker.' + path] = [arr];
    Plotly.restyle(gd, u, [rec.curveNumber]);
  }
  function setPointSize(rec, s){ s = Number(s); if (!isFinite(s) || s <= 0) return; rec.markerSize = s; setMarkerAttr(rec, 'size', s); }
  function setPointColor(rec, c){ rec.markerColor = c; setMarkerAttr(rec, 'color', c); }
  function setPointMaterial(rec, m){
    var spec = MATERIALS[m] || MATERIALS.flat;
    rec.markerMaterial = m;
    setMarkerNested(rec, 'gradient.type', spec.type);
    setMarkerNested(rec, 'gradient.color', spec.edge || rec.markerColor || '#3C5488');
  }
  function setLabelText(rec, t){ rec.text = t; renderLabel(rec); }
  function setLabelFont(rec, s){ s = Number(s); if (!isFinite(s) || s <= 0) return; rec.fontSize = s; relayoutLabel(rec, {'font.size': s}); }
  function setLabelVisible(rec, v){ rec.visible = !!v; renderLabel(rec); }
  function nudgeScale(ev){ if (ev && ev.shiftKey) return 3.0; if (ev && ev.altKey) return 0.35; return 1.0; }
  function nudgeLabel(rec, dx, dy, ev){
    var s = nudgeScale(ev);
    rec.xshift = Number(rec.xshift || 0) + dx * 8 * s;
    rec.yshift = Number(rec.yshift || 0) + dy * 8 * s;
    renderLabel(rec);
  }
  function resetLabelPosition(rec){ rec.xshift = Number(rec.defaultXshift || 0); rec.yshift = Number(rec.defaultYshift || 0); renderLabel(rec); }
  function resetAllLabelPositions(){ LABELS.forEach(resetLabelPosition); }
  function pointHidden(rec){
    var tr = gd.data[rec.curveNumber];
    if (!tr || !isPointTrace(tr)) return false;
    return tr.x[rec.pointNumber] == null || tr.y[rec.pointNumber] == null;
  }
  var undoStack = [];
  function deletePoint(rec){
    var tr = gd.data[rec.curveNumber];
    if (!tr || !isPointTrace(tr)) return;
    var i = rec.pointNumber, xs = cloneArray(tr.x), ys = cloneArray(tr.y), cd = cloneArray(tr.customdata);
    if (i < 0 || i >= xs.length || xs[i] == null) return;
    // remember what we removed so it can be restored (undo)
    undoStack.push({labelIndex: rec.labelIndex, curveNumber: rec.curveNumber, pointNumber: i,
                    x: xs[i], y: ys[i], customdata: (i < cd.length ? cd[i] : null), wasVisible: !!rec.visible});
    xs[i] = null; ys[i] = null; if (i < cd.length) cd[i] = null;
    Plotly.restyle(gd, {x: [xs], y: [ys], customdata: [cd]}, [rec.curveNumber]);
    rec.deleted = true; rec.visible = false; renderLabel(rec);
    updateTrajectory(rec.curveNumber);
    if (selected === rec.labelIndex) { selected = null; renderSelected(); }
    renderList(); updateUndoBtn();
  }
  function undoDelete(){
    var last = undoStack.pop();
    if (!last) return;
    var tr = gd.data[last.curveNumber];
    if (tr) {
      var xs = cloneArray(tr.x), ys = cloneArray(tr.y), cd = cloneArray(tr.customdata);
      while (xs.length <= last.pointNumber) { xs.push(null); ys.push(null); }
      xs[last.pointNumber] = last.x; ys[last.pointNumber] = last.y;
      if (last.customdata != null) { while (cd.length <= last.pointNumber) cd.push(null); cd[last.pointNumber] = last.customdata; }
      Plotly.restyle(gd, {x: [xs], y: [ys], customdata: [cd]}, [last.curveNumber]);
    }
    var rec = LABELS[last.labelIndex];
    if (rec) { rec.deleted = false; rec.visible = last.wasVisible; renderLabel(rec); updateTrajectory(rec.curveNumber); selectPoint(rec.labelIndex); }
    renderList(); updateUndoBtn();
  }
  function updateUndoBtn(){
    var b = document.getElementById(__UNDO_ID__);
    if (b) { b.disabled = undoStack.length === 0; b.textContent = '撤销删除' + (undoStack.length ? ' (' + undoStack.length + ')' : ''); }
  }
  function setLineColor(idx, color){
    var ln = LINES[idx]; if (!ln) return;
    ln.color = color;
    Plotly.restyle(gd, {'line.color': color}, [ln.curveNumber]);
    if (ln.arrowCurve != null) Plotly.restyle(gd, {'marker.color': color}, [ln.arrowCurve]);
  }
  function setAllLabels(mode){
    LABELS.forEach(function(rec){
      if (rec.deleted) { rec.visible = false; renderLabel(rec); return; }
      rec.text = String(rec.defaultText || '');
      if (pointHidden(rec) || mode === 'hide') { rec.visible = false; }
      else if (mode === 'zero') { var t = String(rec.defaultText || ''); var wk = /\d+W$/.test(t); rec.visible = (!wk || /0W$/.test(t)); }
      else { rec.visible = true; }
      renderLabel(rec);
    });
    renderList();
  }

  function el(tag, cls, txt){ var e = document.createElement(tag); if (cls) e.className = cls; if (txt != null) e.textContent = txt; return e; }
  function row(labelText){ var r = el('div', 'be-row'); if (labelText != null) r.appendChild(el('label', null, labelText)); return r; }

  function selectPoint(labelIndex){
    if (labelIndex == null || !LABELS[labelIndex]) return;
    selected = labelIndex; renderSelected(); renderList();
    if (selBox && selBox.closest) { var acc = selBox.closest('details'); if (acc) acc.open = true; }
  }
  function renderSelected(){
    if (!selBox) return;
    selBox.innerHTML = '';
    if (selected == null || !LABELS[selected]) { selBox.appendChild(el('div', 'be-empty', '点击图中的点或标签进行编辑')); return; }
    var rec = LABELS[selected];
    selBox.appendChild(el('div', 'be-id', String(rec.groupName || '') + ' / ' + String(rec.pointName || '')));

    var rT = row('文字'); var iT = el('input'); iT.type = 'text'; iT.value = String(rec.text || '');
    iT.addEventListener('input', function(){ setLabelText(rec, iT.value); }); rT.appendChild(iT); selBox.appendChild(rT);

    var rF = row('字号'); var iF = el('input'); iF.type = 'number'; iF.min = '6'; iF.max = '28'; iF.value = String(rec.fontSize || 9);
    iF.addEventListener('input', function(){ setLabelFont(rec, iF.value); }); rF.appendChild(iF); selBox.appendChild(rF);

    var rS = row('点大小'); var iS = el('input'); iS.type = 'number'; iS.min = '3'; iS.max = '48'; iS.value = String(rec.markerSize || 11);
    iS.addEventListener('input', function(){ setPointSize(rec, iS.value); }); rS.appendChild(iS); selBox.appendChild(rS);

    var rC = row('点颜色'); var iC = el('input'); iC.type = 'color'; iC.value = String(rec.markerColor || '#3C5488');
    iC.addEventListener('input', function(){ setPointColor(rec, iC.value); }); rC.appendChild(iC); selBox.appendChild(rC);

    var rM = row('材质'); var sM = el('select'); sM.className = 'point-material';
    MATERIAL_OPTIONS.forEach(function(o){ var op = el('option', null, o[1]); op.value = o[0]; if (o[0] === (rec.markerMaterial || 'flat')) op.selected = true; sM.appendChild(op); });
    sM.addEventListener('change', function(){ setPointMaterial(rec, sM.value); }); rM.appendChild(sM); selBox.appendChild(rM);

    var rP = row('位置'); var pad = el('div', 'dpad');
    [['u','↑',0,1],['d','↓',0,-1],['l','←',-1,0],['r','→',1,0]].forEach(function(sp){
      var b = el('button', sp[0], sp[1]); b.type = 'button'; b.title = '移动标签；Shift=大步，Alt=小步';
      b.addEventListener('click', function(ev){ nudgeLabel(rec, sp[2], sp[3], ev); }); pad.appendChild(b);
    });
    var cbtn = el('button', 'c', '复位'); cbtn.type = 'button'; cbtn.addEventListener('click', function(){ resetLabelPosition(rec); }); pad.appendChild(cbtn);
    rP.appendChild(pad); selBox.appendChild(rP);

    var rD = el('div', 'be-row'); var bD = el('button', 'be-btn be-del', '删除该点'); bD.type = 'button';
    bD.addEventListener('click', function(){ deletePoint(rec); }); rD.appendChild(bD); selBox.appendChild(rD);
  }
  function renderLines(){
    if (!linesBox) return;
    linesBox.innerHTML = '';
    if (!LINES.length) { linesBox.appendChild(el('div', 'be-empty', '此图暂无可调线条')); return; }
    LINES.forEach(function(ln, idx){
      var r = row(ln.label); var c = el('input'); c.type = 'color'; c.value = String(ln.color || '#8A9099');
      c.addEventListener('input', function(){ setLineColor(idx, c.value); }); r.appendChild(c); linesBox.appendChild(r);
    });
  }
  function renderList(){
    if (!listBox) return;
    listBox.innerHTML = '';
    LABELS.forEach(function(rec){
      var it = el('div', 'be-item' + (selected === rec.labelIndex ? ' sel' : ''));
      var cbx = el('input'); cbx.type = 'checkbox'; cbx.checked = !!rec.visible; cbx.disabled = !!rec.deleted;
      cbx.addEventListener('change', function(){ setLabelVisible(rec, cbx.checked); });
      var nm = el('button', null, (rec.deleted ? '（已删除）' : '') + String(rec.groupName || '') + ' / ' + String(rec.pointName || ''));
      nm.type = 'button'; nm.addEventListener('click', function(){ selectPoint(rec.labelIndex); });
      it.appendChild(cbx); it.appendChild(nm); listBox.appendChild(it);
    });
  }

  renderSelected(); renderLines(); renderList(); updateUndoBtn();

  var undoBtn = document.getElementById(__UNDO_ID__);
  if (undoBtn) undoBtn.addEventListener('click', undoDelete);

  // Force Plotly to size into its flex column right away; otherwise the plot can
  // render at full window width on load (pushing the control panel off-screen)
  // until a manual resize triggers a re-layout.
  function fit(){ try { Plotly.Plots.resize(gd); } catch (e) {} }
  fit();
  if (window.requestAnimationFrame) window.requestAnimationFrame(fit);
  setTimeout(fit, 60); setTimeout(fit, 300);
  window.addEventListener('resize', fit);

  gd.on('plotly_click', function(ev){
    if (!ev || !ev.points || !ev.points.length) return;
    var pt = ev.points[0], tr = gd.data[pt.curveNumber];
    if (!tr || !isPointTrace(tr)) return;
    var i = pt.pointNumber; if (i == null) return;
    var rec = null;
    LABELS.forEach(function(r){ if (r.curveNumber === pt.curveNumber && Number(r.pointNumber) === Number(i)) rec = r; });
    if (rec) selectPoint(rec.labelIndex);
  });
  gd.on('plotly_clickannotation', function(ev){
    if (!ev || ev.index == null) return;
    var rec = null;
    LABELS.forEach(function(r){ if (Number(r.annotationIndex) === Number(ev.index)) rec = r; });
    if (rec) selectPoint(rec.labelIndex);
  });
  if (toolbar) {
    toolbar.addEventListener('click', function(ev){
      var t = ev.target; if (!t || !t.getAttribute) return;
      var a = t.getAttribute('data-action'); if (!a) return;
      if (a === 'reset_positions') resetAllLabelPositions(); else setAllLabels(a);
    });
  }
})();
"""


def plot_clickable_burden_entropy_state_html(
    summary_df: pd.DataFrame,
    experiment: str,
    out_path: Path,
    burden_kind: str = "NMD",
    treatment_group: Optional[str] = None,
    comparator_group: Optional[str] = None,
    normalized: bool = False,
) -> None:
    if not PLOTLY_AVAILABLE:
        out_path.with_suffix(".plotly_missing.txt").write_text(
            "Plotly is not installed. Rebuild with requirements_desktop.txt.", encoding="utf-8"
        )
        return
    meta = _burden_metadata(burden_kind)
    burden_col = meta["mean"]
    if burden_col not in summary_df.columns or "Entropy" not in summary_df.columns:
        return
    sub = summary_df.copy()
    if "Experiment" in sub.columns:
        sub = sub[sub["Experiment"].astype(str) == str(experiment)].copy()
    sub = sub[[c for c in ["Week", "Group", burden_col, "Entropy"] if c in sub.columns]].dropna(subset=[burden_col, "Entropy"]).copy()
    if sub.empty:
        return
    if "Week" in sub.columns:
        sub["Week"] = pd.to_numeric(sub["Week"], errors="coerce").fillna(0).astype(int)
    else:
        sub["Week"] = 0
    if "Group" not in sub.columns:
        sub["Group"] = "Group"
    sub["Group"] = sub["Group"].astype(str)

    if normalized:
        sub["Burden_plot"] = _normalize_display_window(sub[burden_col], sub.get("Group"))
        sub["Entropy_plot"] = _normalize_display_window(sub["Entropy"], sub.get("Group"))
        x_axis_title = f"{meta['short']} normalized display index (0-100)"
        y_axis_title = "Entropy normalized display index (0-100)"
        x_range = [-5.0, 105.0]
        y_range = [-5.0, 105.0]
        title_mode = "normalized"
        hover_value_lines = (
            f"{meta['short']} display: " + "%{customdata[4]:.2f}<br>"
            "Entropy display: %{customdata[5]:.2f}<br>"
            f"{meta['short']} raw: " + "%{customdata[2]:.4g}<br>"
            "Entropy raw: %{customdata[3]:.4g}<br>"
        )
        instruction = (
            "Hover shows normalized display values plus raw values. "
            "Use the label editor above the plot to edit or hide individual labels."
        )
    else:
        sub["Burden_plot"] = pd.to_numeric(sub[burden_col], errors="coerce").astype(float)
        sub["Entropy_plot"] = pd.to_numeric(sub["Entropy"], errors="coerce").astype(float)
        x_min, x_max = float(sub["Burden_plot"].min()), float(sub["Burden_plot"].max())
        y_min, y_max = float(sub["Entropy_plot"].min()), float(sub["Entropy_plot"].max())
        x_pad = max((x_max - x_min) * 0.20, 0.45)
        y_pad = max((y_max - y_min) * 0.20, 0.18)
        x_axis_title = meta["axis"]
        y_axis_title = "Entropy"
        x_range = [x_min - x_pad, x_max + x_pad]
        y_range = [y_min - y_pad, y_max + y_pad]
        title_mode = "raw-value"
        hover_value_lines = (
            f"{meta['short']}: " + "%{customdata[2]:.4g}<br>"
            "Entropy: %{customdata[3]:.4g}<br>"
        )
        instruction = (
            "Hover shows raw values. Use the label editor above the plot to edit or hide individual labels."
        )
    sub = sub.dropna(subset=["Burden_plot", "Entropy_plot"]).copy()
    if sub.empty:
        return
    single_timepoint = _is_single_timepoint(sub)
    groups = _ordered_single_timepoint_groups(sub, burden_col) if single_timepoint else _ordered_groups_for_plot(sub, treatment_group, comparator_group)
    cmap = _ensure_group_color_map(groups, treatment_group, comparator_group)

    if not single_timepoint:
        # Time-series redesign for multi-week (clinical) experiments: x = Week and
        # y = a normalized 0-100 display index, plotting BOTH the burden index and
        # Entropy as one line each per group. The old BRI×Entropy phase plot made
        # weeks pile up and could not be read across time. The raw variant uses the
        # linear anchored 0-100; the normalized variant uses the contrast-enhanced
        # window (both keep raw values in the hover; inference stays on raw metrics).
        if normalized:
            sub["BRI_disp"] = _normalize_display_window(sub[burden_col], sub.get("Group"))
            sub["Ent_disp"] = _normalize_display_window(sub["Entropy"], sub.get("Group"))
            norm_kind = "contrast-enhanced"
        else:
            sub["BRI_disp"] = _normalize_anchor_linear_0_100(sub[burden_col], sub.get("Group"))
            sub["Ent_disp"] = _normalize_anchor_linear_0_100(sub["Entropy"], sub.get("Group"))
            norm_kind = "linear-anchored"
        _wk = pd.to_numeric(sub["Week"], errors="coerce")
        wk_min, wk_max = float(_wk.min()), float(_wk.max())
        wk_pad = max((wk_max - wk_min) * 0.06, 1.0)
        x_axis_title = "Week (时间点)"
        y_axis_title = "Normalized display index (0-100)"
        x_range = [wk_min - wk_pad, wk_max + wk_pad]
        y_range = [-6.0, 106.0]
        title_mode = f"time-series · {norm_kind}"
        instruction = (
            f"x = 时间 (Week)；y = 归一化展示指数 (0-100)。每组显示两条线："
            f"{meta['short']}(实心圆/实线) 与 Entropy(空心圆/虚线)。悬停查看原始值。"
        )

    def _annotation_shift(position: str) -> tuple[int, int]:
        return {
            "top center": (0, 14),
            "bottom center": (0, -14),
            "middle right": (20, 0),
            "middle left": (-20, 0),
            "top right": (16, 12),
            "top left": (-16, 12),
            "bottom right": (16, -12),
            "bottom left": (-16, -12),
        }.get(str(position), (0, 14))

    fig = go.Figure()
    label_records: List[Dict[str, Any]] = []
    curve_specs: List[Dict[str, Any]] = []
    base_marker_size = 14 if single_timepoint else 10
    base_font_size = 10 if single_timepoint else 9

    # Build the list of plotted series. Single-timepoint (animal) = one phase-plot
    # marker series per group. Multi-week (clinical) = a time series per group AND
    # per metric (burden + Entropy), x=Week, y=normalized 0-100.
    series_list: List[Dict[str, Any]] = []
    for idx, group in enumerate(groups):
        g = sub[sub["Group"].astype(str) == str(group)].sort_values("Week")
        if g.empty:
            continue
        color = cmap.get(group, _FALLBACK_COLORS[idx % len(_FALLBACK_COLORS)])
        display_name = display_group_name(group, idx)
        if single_timepoint:
            xs = [float(v) for v in g["Burden_plot"]]
            ys = [float(v) for v in g["Entropy_plot"]]
            custom = [[display_name, display_name, float(r[burden_col]), float(r["Entropy"]),
                       float(r["Burden_plot"]), float(r["Entropy_plot"])] for _, r in g.iterrows()]
            hover = ("Group: %{customdata[0]}<br>Point: %{customdata[1]}<br>"
                     + hover_value_lines + "<b>Click point to delete point</b><extra></extra>")
            series_list.append(dict(
                name=display_name, color=color, xs=xs, ys=ys, custom=custom,
                labels=[display_name for _ in xs], first_visible=False,
                mode="markers", symbol="circle", dash=None, hover=hover, hasLine=False,
            ))
        else:
            weeks = [int(r["Week"]) for _, r in g.iterrows()]
            for mlabel, dispcol, rawcol, dash, symbol in (
                (meta["short"], "BRI_disp", burden_col, "solid", "circle"),
                ("Entropy", "Ent_disp", "Entropy", "dash", "circle-open"),
            ):
                xs = [float(w) for w in weeks]
                ys = [float(r[dispcol]) for _, r in g.iterrows()]
                sname = f"{display_name} · {mlabel}"
                custom = [[display_name, f"{w}W", mlabel, float(r[rawcol]), float(r[dispcol])]
                          for w, (_, r) in zip(weeks, g.iterrows())]
                labels = [sname if i == 0 else f"{w}W" for i, w in enumerate(weeks)]
                hover = ("Group: %{customdata[0]}<br>Metric: %{customdata[2]}<br>"
                         "Week: %{customdata[1]}<br>Normalized: %{customdata[4]:.1f}<br>"
                         "Raw: %{customdata[3]:.4g}<br><b>Click point to delete point</b><extra></extra>")
                series_list.append(dict(
                    name=sname, color=color, xs=xs, ys=ys, custom=custom, labels=labels,
                    first_visible=True, mode="lines+markers", symbol=symbol, dash=dash,
                    hover=hover, hasLine=True,
                ))

    for s in series_list:
        n_pts = len(s["xs"])
        if n_pts == 0:
            continue
        color = s["color"]
        point_curve_number = len(fig.data)
        line_kw = dict(color=color, width=2.1, shape="spline")
        if s.get("dash"):
            line_kw["dash"] = s["dash"]
        fig.add_trace(go.Scatter(
            x=s["xs"], y=s["ys"], mode=s["mode"], name=s["name"], customdata=s["custom"],
            meta=dict(bioentropyRole="point", bioentropyPair=s["name"]),
            # Per-point marker arrays so each dot can be individually resized,
            # recolored, deleted, or given a different rendering "material" (gradient)
            # from the editor. Multi-week Entropy uses open circles + dashed line to
            # tell it apart from the burden line within the same (colored) group.
            marker=dict(
                size=[base_marker_size] * n_pts,
                color=[color] * n_pts,
                symbol=[s["symbol"]] * n_pts,
                gradient=dict(type=["none"] * n_pts, color=[color] * n_pts),
                line=dict(color="white", width=1.2),
            ),
            line=line_kw,
            connectgaps=True,
            hovertemplate=s["hover"],
        ))
        curve_specs.append({
            "curveNumber": point_curve_number, "name": s["name"], "color": color,
            "baseSize": base_marker_size, "hasLine": s["hasLine"], "arrowCurve": None,
        })
        text_positions = [["top center", "bottom center", "middle right", "middle left"][i % 4]
                          for i in range(n_pts)]
        for point_number in range(n_pts):
            xshift, yshift = _annotation_shift(text_positions[point_number])
            # Multi-week time series: the 4 lines are identified by the legend, so
            # hide the per-point week labels by default (avoids heavy clutter); they
            # can be revealed from the editor. Animal phase plot keeps group labels.
            vis = bool(single_timepoint)
            label_records.append({
                "labelIndex": len(label_records),
                "annotationIndex": None,
                "curveNumber": point_curve_number,
                "pointNumber": point_number,
                "pairName": s["name"],
                "groupName": s["custom"][point_number][0],
                "pointName": s["labels"][point_number],
                "defaultText": s["labels"][point_number],
                "text": s["labels"][point_number],
                "x": float(s["xs"][point_number]),
                "y": float(s["ys"][point_number]),
                "xshift": xshift,
                "yshift": yshift,
                "defaultXshift": xshift,
                "defaultYshift": yshift,
                "visible": vis,
                "fontSize": base_font_size,
                "defaultFontSize": base_font_size,
                "markerSize": base_marker_size,
                "markerColor": color,
                "markerMaterial": "flat",
                "defaultMarkerColor": color,
            })

    # Smooth connector for the single-timepoint (animal) case: one spline threading
    # through the group points. The trajectory is ordered by DESCENDING disease
    # deviation, so it originates at the Model group (peak entropy / farthest from
    # Normal) and progresses down toward the Normal group; the head-arrow (drawn on
    # the last segment) therefore points at the Normal end.
    connector_color = "#8A9099"
    connector_spec = None
    if single_timepoint and len(curve_specs) >= 2:
        member_points = [
            (s["curveNumber"], float(fig.data[s["curveNumber"]].x[0]), float(fig.data[s["curveNumber"]].y[0]))
            for s in curve_specs if len(fig.data[s["curveNumber"]].x)
        ]
        member_points.sort(key=lambda t: t[1], reverse=True)  # Model (high burden) -> Normal (low)
        members = [m[0] for m in member_points]
        cx = [m[1] for m in member_points]
        cy = [m[2] for m in member_points]
        connector_curve = len(fig.data)
        fig.add_trace(go.Scatter(
            x=cx,
            y=cy,
            mode="lines",
            line=dict(color=connector_color, width=2.4, shape="spline"),
            connectgaps=True,
            meta=dict(bioentropyRole="connector"),
            hoverinfo="skip",
            showlegend=False,
        ))
        connector_spec = {"curveNumber": connector_curve, "members": members, "arrowCurve": None, "color": connector_color}

    # Labels first, so each label's annotationIndex == its labelIndex (keeps the
    # JS label<->annotation mapping and plotly_clickannotation trivial).
    label_annotations: List[Dict[str, Any]] = []
    for record in label_records:
        record["annotationIndex"] = record["labelIndex"]
        label_annotations.append(dict(
            x=record["x"],
            y=record["y"],
            xref="x",
            yref="y",
            text=record["text"],
            showarrow=False,
            xshift=record["xshift"],
            yshift=record["yshift"],
            xanchor="center",
            yanchor="middle",
            font=dict(size=record["fontSize"], color="#2D3748"),
            bgcolor="rgba(255,255,255,0.0)",
            captureevents=True,  # clickable -> select the point for editing
            visible=bool(record.get("visible", True)),
        ))

    # No directional arrows on the trajectory map (per request): the smooth spline
    # alone conveys the path. `arrowCurve` stays None so the JS arrow hooks no-op.
    instruction_annotation = dict(
        text=instruction,
        x=0.5, y=1.02, xref="paper", yref="paper", showarrow=False,
        font=dict(size=12, color="#4B5563"), align="center",
    )
    fig.update_layout(
        title=f"Experiment {experiment}: {title_mode} {meta['short']}–Entropy trajectory map",
        template="plotly_white",
        xaxis=dict(title=x_axis_title, range=x_range, zeroline=False),
        yaxis=dict(title=y_axis_title, range=y_range, zeroline=False, showgrid=False),
        legend=dict(orientation="h", y=1.08, x=0.5, xanchor="center"),
        margin=dict(l=78, r=34, t=120, b=68),
        annotations=label_annotations + [instruction_annotation],
    )
    safe_exp = "".join(ch if ch.isalnum() else "_" for ch in str(experiment))
    div_id = f"bioentropy_clickable_{safe_exp}_{'norm' if normalized else 'raw'}"
    panel_id = f"{div_id}_panel"
    sel_id = f"{div_id}_sel"
    lines_id = f"{div_id}_lines"
    list_id = f"{div_id}_list"
    toolbar_id = f"{div_id}_toolbar"
    undo_id = f"{div_id}_undo"
    plot_html = fig.to_html(
        include_plotlyjs="cdn", full_html=False, div_id=div_id,
        config={"scrollZoom": True, "displaylogo": False, "responsive": True},
    )
    # Lines whose colour the user can customize: the animal connector, or each
    # clinical group trajectory. Setting a colour recolours the line and its arrow.
    if single_timepoint and connector_spec is not None:
        lines_spec = [{
            "label": "连接线 (Model→Normal)",
            "curveNumber": connector_spec["curveNumber"],
            "arrowCurve": connector_spec.get("arrowCurve"),
            "color": connector_color,
        }]
    else:
        lines_spec = [{
            "label": s["name"],
            "curveNumber": s["curveNumber"],
            "arrowCurve": s.get("arrowCurve"),
            "color": s["color"],
        } for s in curve_specs if s.get("hasLine")]

    script = (
        "<script>\n"
        + _CLICKABLE_STATE_JS
        .replace("__DIV_ID__", json.dumps(div_id))
        .replace("__SEL_ID__", json.dumps(sel_id))
        .replace("__LINES_ID__", json.dumps(lines_id))
        .replace("__LIST_ID__", json.dumps(list_id))
        .replace("__TOOLBAR_ID__", json.dumps(toolbar_id))
        .replace("__UNDO_ID__", json.dumps(undo_id))
        .replace("__LABELS__", json.dumps(label_records, ensure_ascii=False))
        .replace("__CURVES__", json.dumps(curve_specs, ensure_ascii=False))
        .replace("__CONNECTOR__", json.dumps(connector_spec, ensure_ascii=False))
        .replace("__LINES__", json.dumps(lines_spec, ensure_ascii=False))
        .replace("__SINGLE__", "true" if single_timepoint else "false")
        + "\n</script>"
    )
    css = (
        "<style>"
        "body{margin:0;font-family:Arial,'Microsoft YaHei UI',sans-serif;background:#fff;color:#243B53;}"
        ".be-wrap{display:flex;gap:14px;align-items:stretch;padding:10px 12px;box-sizing:border-box;}"
        ".be-plot{flex:1 1 auto;min-width:0;height:82vh;}"
        ".be-plot .plotly-graph-div{width:100%!important;height:100%!important;}"
        ".be-panel{flex:0 0 340px;max-height:94vh;overflow:auto;border:1px solid #D8E0EA;border-radius:10px;background:#FBFDFF;padding:10px 12px;box-sizing:border-box;}"
        ".be-title{font-size:14px;font-weight:700;margin-bottom:4px;}"
        ".be-note{font-size:11px;color:#7A8699;margin-bottom:8px;line-height:1.5;}"
        ".be-acc{border:1px solid #E6EDF5;border-radius:8px;margin-bottom:8px;background:#fff;}"
        ".be-acc>summary{cursor:pointer;padding:8px 10px;font-size:13px;font-weight:600;color:#243B53;}"
        ".be-acc>div{padding:8px 10px;border-top:1px solid #EEF2F7;}"
        ".be-empty{color:#8593A8;font-size:12px;}"
        ".be-id{font-size:12px;font-weight:600;margin-bottom:6px;color:#334E68;}"
        ".be-row{display:flex;align-items:center;gap:8px;margin:6px 0;font-size:12px;}"
        ".be-row>label{flex:0 0 58px;color:#4A5568;}"
        ".be-row input[type='text']{flex:1 1 auto;min-width:0;border:1px solid #CBD5E1;border-radius:5px;padding:4px 6px;font-size:12px;}"
        ".be-row input[type='number']{width:60px;border:1px solid #CBD5E1;border-radius:5px;padding:3px 5px;font-size:12px;}"
        ".be-row input[type='color']{width:44px;height:26px;border:1px solid #CBD5E1;border-radius:5px;padding:0;background:#fff;cursor:pointer;}"
        ".be-row select{border:1px solid #CBD5E1;border-radius:5px;padding:3px 5px;font-size:12px;background:#fff;}"
        ".be-btn{border:1px solid #CBD5E1;background:#F8FAFC;border-radius:6px;padding:5px 9px;font-size:12px;cursor:pointer;color:#243B53;}"
        ".be-btn:hover{background:#EEF2F7;}"
        ".be-del{border:1px solid #F0B4B4;background:#FEF2F2;color:#B4232A;}"
        ".be-del:hover{background:#FEE2E2;}"
        ".be-undo{border:1px solid #B7C7DE;background:#EEF3FA;color:#274472;font-weight:600;}"
        ".be-undo:hover:not(:disabled){background:#E1EAF6;}"
        ".be-btn:disabled{opacity:0.5;cursor:not-allowed;}"
        ".be-toolbar{display:flex;flex-wrap:wrap;gap:6px;margin-bottom:8px;}"
        ".be-list{max-height:240px;overflow:auto;}"
        ".be-item{display:flex;align-items:center;gap:6px;padding:2px;border-bottom:1px solid #F0F4F9;font-size:12px;}"
        ".be-item>button{flex:1 1 auto;text-align:left;border:none;background:none;cursor:pointer;color:#243B53;padding:2px 4px;border-radius:4px;}"
        ".be-item>button:hover{background:#EEF2F7;}"
        ".be-item.sel>button{background:#E3ECF7;font-weight:600;}"
        ".dpad{position:relative;width:92px;height:92px;}"
        ".dpad button{position:absolute;border:1px solid #CBD5E1;background:#fff;border-radius:8px;width:28px;height:28px;padding:0;font-size:14px;line-height:1;cursor:pointer;}"
        ".dpad button:hover{background:#EEF2F7;}"
        ".dpad .u{top:0;left:32px;}.dpad .d{bottom:0;left:32px;}.dpad .l{left:0;top:32px;}.dpad .r{right:0;top:32px;}"
        ".dpad .c{top:32px;left:32px;width:28px;height:28px;border-radius:50%;font-size:10px;background:#F3F7FB;}"
        "@media(max-width:900px){.be-wrap{flex-direction:column;}.be-panel{flex:1 1 auto;}.be-plot{height:60vh;}}"
        "</style>"
    )
    body = (
        "<div class='be-wrap'>"
        f"<div class='be-plot'>{plot_html}</div>"
        f"<aside class='be-panel' id='{panel_id}'>"
        "<div class='be-title'>交互控制面板</div>"
        "<div class='be-note'>点击左侧图中的点或标签进行编辑。此交互只影响本页显示，不改变结果文件；刷新页面即可复位。</div>"
        f"<div class='be-toolbar'><button type='button' id='{undo_id}' class='be-btn be-undo' disabled>撤销删除</button></div>"
        "<details class='be-acc' open><summary>① 选中点编辑（文字/字号/大小/颜色/材质/位置/删除）</summary>"
        f"<div id='{sel_id}'></div></details>"
        "<details class='be-acc' open><summary>② 线条颜色</summary>"
        f"<div id='{lines_id}'></div></details>"
        "<details class='be-acc'><summary>③ 全部标签 / 批量操作</summary>"
        f"<div id='{toolbar_id}' class='be-toolbar'>"
        "<button type='button' class='be-btn' data-action='hide'>隐藏全部文字</button>"
        "<button type='button' class='be-btn' data-action='zero'>只显示0W/组名</button>"
        "<button type='button' class='be-btn' data-action='show'>恢复全部文字</button>"
        "<button type='button' class='be-btn' data-action='reset_positions'>复位全部位置</button>"
        "</div>"
        f"<div id='{list_id}' class='be-list'></div>"
        "</details>"
        "</aside>"
        "</div>"
    )
    out_path.write_text(
        "<!DOCTYPE html><html><head><meta charset='utf-8'>"
        "<title>BioEntropy trajectory map</title>"
        + css
        + "</head><body>"
        + body + script
        + "</body></html>",
        encoding="utf-8",
    )


def plot_state_time_small_multiples(
    summary_df: pd.DataFrame,
    experiment: str,
    out_path: Path,
    burden_kind: str = "NMD",
    treatment_group: Optional[str] = None,
    comparator_group: Optional[str] = None,
    comparison_df: Optional[pd.DataFrame] = None,
    normal_band: Optional[Dict[str, float]] = None,
) -> None:
    """Small-multiples phase portrait: one panel per week, groups placed by state.

    Each week is its own facet; within it the two arms are bubbles positioned in
    the (entropy, burden) state plane on axes shared across all panels, so the
    week-to-week march is directly comparable. The diseased corner is top-right,
    the healthy corner bottom-left; a thin connector + significance star shows the
    per-week gap. Multi-week experiments only.
    """
    meta = _burden_metadata(burden_kind)
    bcol = meta["mean"]
    short = str(meta["short"])
    if bcol not in summary_df.columns or "Entropy" not in summary_df.columns:
        return
    sub = summary_df.copy()
    if "Experiment" in sub.columns:
        sub = sub[sub["Experiment"].astype(str) == str(experiment)].copy()
    sub = sub[[c for c in ["Week", "Group", bcol, "Entropy"] if c in sub.columns]].dropna(subset=[bcol, "Entropy"]).copy()
    if sub.empty or _is_single_timepoint(sub):
        return
    sub["Week"] = pd.to_numeric(sub["Week"], errors="coerce").astype(int)
    sub["Group"] = sub["Group"].astype(str)
    weeks = sorted(sub["Week"].unique().tolist())
    groups = _ordered_groups_for_plot(sub, treatment_group, comparator_group)
    cmap = _ensure_group_color_map(groups, treatment_group, comparator_group)
    comp = _comparison_lookup(comparison_df)

    # Per-subject burden values (for the point overlay), keyed by (group, week).
    psrc = str(bcol)[:-5] if str(bcol).endswith("_mean") else str(bcol)
    sl = getattr(_dstate, "CURRENT_SAMPLE_LEVEL", None)
    sl_exp = None
    if isinstance(sl, pd.DataFrame) and not sl.empty and {psrc, "Week", "Group"}.issubset(sl.columns):
        sl_exp = sl
        if "Experiment" in sl_exp.columns:
            sl_exp = sl_exp[sl_exp["Experiment"].astype(str) == str(experiment)]
        if sl_exp.empty:
            sl_exp = None

    # Shared, tight state-plane limits so inter-group separation is legible;
    # widen the burden axis to fit the per-subject spread (robust 0.5–99.5%).
    ex = pd.to_numeric(sub["Entropy"], errors="coerce")
    by = pd.to_numeric(sub[bcol], errors="coerce")
    by_lo, by_hi = float(by.min()), float(by.max())
    # Adaptive label precision from the magnitude of the plotted burden values.
    # NMD-scale figures (values >= 1) keep one decimal so their pixels stay
    # identical to the regression baseline; small NRBS values (~0.07-0.12)
    # would all collapse to "0.1" at one decimal, so give them enough decimals
    # to stay distinguishable.
    _byabs = np.abs(by.to_numpy(dtype=float))
    _byabs = _byabs[np.isfinite(_byabs) & (_byabs > 0)]
    _absmax = float(np.max(_byabs)) if _byabs.size else 1.0
    if _absmax >= 1.0:
        _lab_dec = 1
    elif _absmax >= 0.1:
        _lab_dec = 3
    elif _absmax >= 0.01:
        _lab_dec = 4
    else:
        _lab_dec = 5
    if sl_exp is not None:
        psv = pd.to_numeric(sl_exp[psrc], errors="coerce")
        psv = psv[np.isfinite(psv)]
        if len(psv):
            by_lo = min(by_lo, float(psv.quantile(0.005)))
            by_hi = max(by_hi, float(psv.quantile(0.995)))
    exs = max(float(ex.max() - ex.min()), 1e-6)
    bys = max(by_hi - by_lo, 1e-6)
    xlim = (float(ex.min()) - 0.18 * exs, float(ex.max()) + 0.18 * exs)
    ylim = (by_lo - 0.10 * bys, by_hi + 0.14 * bys)

    ncols = len(weeks)
    fig, axes = plt.subplots(1, ncols, figsize=(max(1.55 * ncols + 1.4, 8.0), 4.7),
                             sharex=True, sharey=True)
    if ncols == 1:
        axes = [axes]
    p_col = meta["p"]
    _rng = np.random.RandomState(7)

    for i, wk in enumerate(weeks):
        ax = axes[i]
        ax.set_facecolor(_LIGHT_BG)
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        wk_df = sub[sub["Week"] == wk]
        pts = {}
        for gi, group in enumerate(groups):
            row = wk_df[wk_df["Group"] == str(group)]
            if row.empty:
                continue
            xv = float(pd.to_numeric(row["Entropy"], errors="coerce").iloc[0])
            yv = float(pd.to_numeric(row[bcol], errors="coerce").iloc[0])
            pts[group] = (xv, yv)
        # connector between the two arms
        if len(pts) == 2:
            (x0, y0), (x1, y1) = list(pts.values())
            ax.plot([x0, x1], [y0, y1], color="#9aa0a6", linewidth=1.1, alpha=0.8, zorder=2)
        for gi, group in enumerate(groups):
            if group not in pts:
                continue
            xv, yv = pts[group]
            color = cmap.get(group, _FALLBACK_COLORS[gi % len(_FALLBACK_COLORS)])
            # Per-subject burden spread, jittered around the group's entropy x.
            if sl_exp is not None:
                gv = pd.to_numeric(
                    sl_exp.loc[(sl_exp["Group"].astype(str) == str(group))
                               & (pd.to_numeric(sl_exp["Week"], errors="coerce") == wk), psrc],
                    errors="coerce").to_numpy(dtype=float)
                gv = gv[np.isfinite(gv)]
                if gv.size:
                    jx = xv + (_rng.rand(gv.size) - 0.5) * (exs * 0.14)
                    ax.scatter(jx, gv, s=5, color=color, alpha=0.14, edgecolors="none", zorder=1)
            ax.scatter([xv], [yv], s=360, color=color, alpha=0.9, edgecolor="white",
                       linewidth=1.4, zorder=4)
            # Value label offset above (first group) / below (second) the bubble
            # with a white halo, so the two labels never overlap when bubbles do.
            ax.annotate(f"{yv:.{_lab_dec}f}", (xv, yv), xytext=(0, 11 if gi == 0 else -11),
                        textcoords="offset points", fontsize=7.0, color=color,
                        ha="center", va="bottom" if gi == 0 else "top", zorder=6, fontweight="bold",
                        path_effects=[pe.withStroke(linewidth=2.2, foreground="white")])
        # per-week significance of the burden gap
        c = comp[comp["Week"] == wk] if not comp.empty else pd.DataFrame()
        if not c.empty:
            crow = c.iloc[0]
            p = crow.get(f"{short}_p_mannwhitney", crow.get(p_col, crow.get("NMD_p_boot", np.nan)))
            star = _signif_label(p) if pd.notna(p) else ""
            if star:
                ax.text(0.5, 0.965, star, transform=ax.transAxes, ha="center", va="top",
                        fontsize=10.5, fontweight="bold", color=_TEXT_COLOR)
        ax.set_title(f"{wk}W", fontsize=11, fontweight="bold", pad=3)
        ax.tick_params(labelbottom=False, labelleft=(i == 0), labelsize=8, length=0)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.grid(True, color=_GRID_COLOR, linewidth=0.5, alpha=0.5)

    axes[0].set_ylabel(f"{short}  ↓ closer to Normal", fontsize=11)
    handles = [plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=cmap.get(g, _FALLBACK_COLORS[i % len(_FALLBACK_COLORS)]),
                          markeredgecolor="white", markersize=12) for i, g in enumerate(groups)]
    labels = [display_group_name(g, i) for i, g in enumerate(groups)]
    fig.legend(handles, labels, loc="upper right", ncol=len(groups), frameon=False,
               bbox_to_anchor=(0.995, 0.995), fontsize=10)
    fig.suptitle(f"Experiment {experiment}: weekly {short}–entropy state (shared axes; diseased top-right → healthy bottom-left)",
                 fontsize=12.5, fontweight="bold", x=0.5, y=1.02)
    band_note = ""
    if normal_band and np.isfinite(normal_band.get("median", np.nan)):
        band_note = (f"  Normal {short} ≈ {normal_band['median']:.1f} "
                     f"(healthy 95% {normal_band.get('p2_5', float('nan')):.1f}–{normal_band.get('p97_5', float('nan')):.1f}), below the plotted range.")
    fig.text(0.5, -0.02,
             f"x = system entropy, y = {short} (distance from Normal). Stars: burden-gap Mann–Whitney p (*<0.05, **<0.01, ***<0.001)." + band_note,
             ha="center", va="top", fontsize=8.6, color="#555555")
    fig.text(0.5, -0.075, "System entropy →", ha="center", va="top", fontsize=10.5, color=_TEXT_COLOR)
    fig.tight_layout(rect=(0.02, 0.02, 1, 0.94))
    fig.savefig(str(out_path), dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_single_timepoint_nmd_entropy_state(
    summary_df: pd.DataFrame,
    experiment: str,
    out_path: Path,
    burden_kind: str = "NMD",
    treatment_group: Optional[str] = None,
    comparator_group: Optional[str] = None,
    normalized: bool = False,
) -> None:
    """Single-timepoint burden-vs-entropy state plot.

    One marker per group placed along the entropy ordering (healthy → diseased,
    evenly spaced so labels never overlap), y = burden, joined by a smooth PCHIP
    curve. With ``normalized=True`` the y-axis is the 0–100 display index (Normal
    anchored to 0, model to 100). No-op for multi-week data.
    """
    meta = _burden_metadata(burden_kind)
    bcol = meta["mean"]
    short = str(meta["short"])
    if bcol not in summary_df.columns or "Entropy" not in summary_df.columns:
        return
    sub = summary_df.copy()
    if "Experiment" in sub.columns:
        sub = sub[sub["Experiment"].astype(str) == str(experiment)].copy()
    sub = sub[[c for c in ["Week", "Group", bcol, "Entropy"] if c in sub.columns]].dropna(subset=[bcol, "Entropy"]).copy()
    if sub.empty or not _is_single_timepoint(sub):
        return
    sub["Group"] = sub["Group"].astype(str)
    # y = normalized display index (0-100) or the raw burden value.
    if normalized:
        sub["_yv"] = pd.to_numeric(_normalize_anchor_linear_0_100(sub[bcol], sub["Group"]), errors="coerce")
        ylabel = f"{short} display index (0–100)"
        title_extra = " · normalized display"
    else:
        sub["_yv"] = pd.to_numeric(sub[bcol], errors="coerce")
        ylabel = f"{short}  (disease deviation)"
        title_extra = ""
    sub = sub.dropna(subset=["_yv"])
    # Sort groups left→right strictly by increasing raw burden (NMD): healthy →
    # diseased. This is a direct value sort (NOT the custom/default display group
    # order), so the curve reads as a clean monotone rise. Both the raw and
    # normalized versions use the same raw-NMD order (the normalized index is
    # monotone in it), so the two figures line up group-for-group.
    _scores: Dict[str, float] = {}
    for _g in sub["Group"].unique():
        _v = pd.to_numeric(sub.loc[sub["Group"] == _g, bcol], errors="coerce")
        _scores[str(_g)] = float(_v.mean()) if _v.notna().any() else float("inf")
    groups = sorted(_scores, key=lambda g: (_scores[g], g))
    rows = _metric_group_rows(sub, groups)
    if rows.empty or len(rows) < 2:
        return
    xpos = np.arange(len(groups), dtype=float)
    y = pd.to_numeric(rows["_yv"], errors="coerce").to_numpy(dtype=float)
    labels = [display_group_name(g, i) for i, g in enumerate(groups)]
    cmap = _ensure_group_color_map(groups, treatment_group, comparator_group)
    colors = [cmap.get(g, _FALLBACK_COLORS[i % len(_FALLBACK_COLORS)]) for i, g in enumerate(groups)]

    fig, ax = plt.subplots(figsize=(8.0, 6.4))
    # Smooth curve through the points (PCHIP avoids the overshoot a cubic spline
    # would add); fall back to straight segments if scipy is unavailable.
    if len(xpos) >= 3:
        try:
            from scipy.interpolate import PchipInterpolator
            xi = np.linspace(float(xpos.min()), float(xpos.max()), 200)
            ax.plot(xi, PchipInterpolator(xpos, y)(xi), color="#1F3B57", linewidth=1.8, zorder=2)
        except Exception:
            ax.plot(xpos, y, color="#1F3B57", linewidth=1.8, zorder=2)
    else:
        ax.plot(xpos, y, color="#1F3B57", linewidth=1.8, zorder=2)

    ax.scatter(xpos, y, s=300, c=colors, edgecolor="white", linewidth=1.6, zorder=4)
    # Overlay each subject's burden value (jittered) around its group marker.
    pv = _overlay_group_points(
        ax, xpos, groups,
        _per_sample_points_for_metric(experiment, bcol, groups, normalize=("linear" if normalized else None)),
        jitter_width=0.34, size=14)
    ax.set_xticks(list(xpos))
    rot = 20 if len(groups) > 4 else 0
    ax.set_xticklabels(labels, fontsize=10.5, rotation=rot, ha=("right" if rot else "center"))
    ax.set_xlabel(f"Groups ordered by {short}  (healthy → diseased)", fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(f"Experiment {experiment}: {short}–entropy state (healthy → diseased){title_extra}",
                 fontsize=12.5, fontweight="bold")
    ax.set_xlim(xpos.min() - 0.4, xpos.max() + 0.4)
    y_lo = min(float(np.nanmin(y)), min(pv) if pv else float(np.nanmin(y)))
    y_hi = max(float(np.nanmax(y)), max(pv) if pv else float(np.nanmax(y)))
    ypad = max((y_hi - y_lo) * 0.14, 1e-6)
    ax.set_ylim(y_lo - ypad, y_hi + ypad)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(str(out_path), dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_burden_entropy_dual_heatmap(
    summary_df: pd.DataFrame,
    experiment: str,
    out_path: Path,
    burden_kind: str = "BRI",
    treatment_group: Optional[str] = None,
    comparator_group: Optional[str] = None,
) -> None:
    meta = _burden_metadata(burden_kind)
    burden_col = meta["mean"]

    sub = summary_df.copy()
    if "Experiment" in sub.columns:
        sub = sub[sub["Experiment"].astype(str) == str(experiment)].copy()
    keep = [c for c in ["Week", "Group", burden_col, "Entropy"] if c in sub.columns]
    sub = sub[keep].dropna(subset=[burden_col, "Entropy"]).copy()
    if sub.empty:
        return
    sub["Week"] = pd.to_numeric(sub["Week"], errors="coerce").astype(int)
    groups = _ordered_groups_for_plot(sub, treatment_group, comparator_group)
    weeks = sorted(sub["Week"].dropna().astype(int).unique().tolist())

    burden_wide = sub.pivot(index="Group", columns="Week", values=burden_col).reindex(index=groups, columns=weeks)
    ent_wide = sub.pivot(index="Group", columns="Week", values="Entropy").reindex(index=groups, columns=weeks)

    def _norm_arr(arr):
        finite = arr[np.isfinite(arr)]
        if finite.size == 0:
            return np.full_like(arr, np.nan, dtype=float)
        vmin = float(finite.min()); vmax = float(finite.max())
        if math.isclose(vmin, vmax):
            return np.full_like(arr, 50.0, dtype=float)
        return (arr - vmin) / (vmax - vmin) * 100.0

    burden_norm = _norm_arr(burden_wide.to_numpy(dtype=float))
    ent_norm = _norm_arr(ent_wide.to_numpy(dtype=float))

    fig, axes = plt.subplots(2, 1, figsize=(1.22 * max(6, len(weeks)) + 3.2, 1.0 * max(3, len(groups)) + 3.8), sharex=True)
    panels = [
        (axes[0], burden_norm, burden_wide, f"{meta['short']} burden (0–100 within experiment)"),
        (axes[1], ent_norm, ent_wide, "Entropy burden (0–100 within experiment)"),
    ]
    group_labels = [display_group_name(g, i) for i, g in enumerate(groups)]

    for ax, zvals, raw_df, title in panels:
        im = ax.imshow(zvals, aspect="auto", cmap="RdYlBu_r", vmin=0, vmax=100)
        ax.set_title(title, fontsize=11.0, pad=8)
        ax.set_yticks(np.arange(len(groups)))
        ax.set_yticklabels(group_labels, fontsize=9)
        ax.set_xticks(np.arange(len(weeks)))
        ax.set_xticklabels([f"{w}W" for w in weeks], fontsize=9)
        for i in range(len(groups)):
            for j in range(len(weeks)):
                raw_val = raw_df.iloc[i, j]
                if pd.notna(raw_val) and len(groups) <= 6 and len(weeks) <= 9:
                    ax.text(j, i, f"{raw_val:.2f}", ha="center", va="center", fontsize=7.2, color="#222222")
        cbar = fig.colorbar(im, ax=ax, fraction=0.024, pad=0.02)
        cbar.ax.tick_params(labelsize=8)

    axes[1].set_xlabel("Week")
    note = (
        "This figure avoids 3D occlusion and is recommended when many experimental groups are present. "
        "Both panels are normalized to the same 0–100 display range within the experiment. "
        f"If {meta['short']} and Entropy cool down together across the same cells, the two measures are changing synergistically."
    )
    _legacy._add_top_title_and_note(
        fig,
        title=f"Experiment {experiment}: {meta['heat_title']}",
        note=note,
        legend_handles=None,
        legend_labels=None,
        title_y=0.985,
        note_y=0.948,
        note_width=118,
        note_fontsize=8.7,
    )
    fig.tight_layout(rect=[0.05, 0.05, 0.98, 0.84])
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_burden_entropy_grouped_bar(
    summary_df: pd.DataFrame,
    experiment: str,
    out_path: Path,
    burden_kind: str = "BRI",
    treatment_group: Optional[str] = None,
    comparator_group: Optional[str] = None,
) -> None:
    meta = _burden_metadata(burden_kind)
    burden_col = meta["mean"]

    sub = summary_df.copy()
    if "Experiment" in sub.columns:
        sub = sub[sub["Experiment"].astype(str) == str(experiment)].copy()
    keep = [c for c in ["Week", "Group", burden_col, "Entropy"] if c in sub.columns]
    sub = sub[keep].dropna(subset=[burden_col, "Entropy"]).copy()
    if sub.empty:
        return
    sub["Week"] = pd.to_numeric(sub["Week"], errors="coerce").astype(int)
    groups = _ordered_groups_for_plot(sub, treatment_group, comparator_group)
    weeks = sorted(sub["Week"].dropna().astype(int).unique().tolist())
    if _plot_burden_entropy_single_group_bar(sub, experiment, out_path, burden_col, meta, treatment_group, comparator_group):
        return
    sub["Burden_display"] = _normalize_display_window(sub[burden_col], sub.get("Group"))
    sub["Entropy_display"] = _normalize_display_window(sub["Entropy"], sub.get("Group"))

    n = len(groups)
    fig, axes = plt.subplots(1, n, figsize=(1.20 * max(6, len(weeks)) * max(1.0, n * 0.58) + 3.0, 6.6), sharey=True, squeeze=False)
    axes = axes.ravel()

    metric_handles = [
        Patch(facecolor=_NMD_COLOR if str(meta.get("short", "")).upper() == "NMD" else _BRI_COLOR, edgecolor="none", alpha=0.95),
        Patch(facecolor=_ENTROPY_COLOR, edgecolor="none", alpha=0.95),
    ]
    metric_labels = [meta["axis"], "Entropy"]

    for idx_g, group in enumerate(groups):
        ax = axes[idx_g]
        g = sub[sub["Group"].astype(str) == group].set_index("Week").reindex(weeks)
        burden_color = _NMD_COLOR if str(meta.get("short", "")).upper() == "NMD" else _BRI_COLOR
        entropy_color = _ENTROPY_COLOR
        x = np.arange(len(weeks), dtype=float)
        width = 0.36

        b1 = ax.bar(x - width/2, g["Burden_display"], width=width, color=burden_color, alpha=0.92,
                    edgecolor="#FFFFFF", linewidth=0.55)
        b2 = ax.bar(x + width/2, g["Entropy_display"], width=width, color=entropy_color, alpha=0.92,
                    edgecolor="#FFFFFF", linewidth=0.55)

        for bi, raw in zip(b1, g[burden_col].tolist()):
            if pd.notna(raw):
                ax.text(bi.get_x() + bi.get_width()/2, bi.get_height() + 2.2, f"{raw:.2f}",
                        ha="center", va="bottom", fontsize=8.0)
        for bi, raw in zip(b2, g["Entropy"].tolist()):
            if pd.notna(raw):
                ax.text(bi.get_x() + bi.get_width()/2, bi.get_height() + 2.2, f"{raw:.2f}",
                        ha="center", va="bottom", fontsize=8.0)

        ax.set_title(display_group_name(group, idx_g), fontsize=15, pad=10)
        ax.set_xticks(x)
        ax.set_xticklabels([f"{w}W" for w in weeks], fontsize=10)
        ax.set_xlabel("Week")
        _legacy._format_axis(ax)
        _apply_academic_axis(ax, "y")

    axes[0].set_ylabel("Normalized group mean (0–100)")
    ymax = max(110, np.nanmax(sub[["Burden_display", "Entropy_display"]].to_numpy(dtype=float)) + 15)
    for ax in axes:
        ax.set_ylim(0, ymax)

    note = (
        f"Each panel represents one group. Within each week, the left / right bars show {meta['short']} and Entropy, respectively. "
        f"Bar heights use anchored display normalization when Normal and model controls are present (Normal=0, model=100), then gamma={DISPLAY_CONTRAST_GAMMA:g} display contrast; labels above bars are raw values. "
        "This 2D grouped-bar view avoids 3D occlusion and is recommended for intuitive side-by-side comparison across weeks."
    )
    _legacy._add_top_title_and_note(
        fig,
        title=f"Experiment {experiment}: {meta['grouped_title']}",
        note=note,
        legend_handles=metric_handles,
        legend_labels=metric_labels,
        legend_ncol=2,
        legend_y=0.885,
        title_y=0.985,
        note_y=0.948,
        note_width=128,
        note_fontsize=8.8,
    )
    fig.tight_layout(rect=[0.05, 0.05, 0.98, 0.80])
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _front_sorted_groups_by_low_values(sub: pd.DataFrame, groups: List[str], burden_col: str) -> List[str]:
    """Sort groups so lower combined display height is placed closer to the default camera/front."""
    if not groups:
        return groups
    tmp = sub.copy()
    if "Burden_display" not in tmp.columns:
        tmp["Burden_display"] = _normalize_0_100(tmp[burden_col])
    if "Entropy_display" not in tmp.columns:
        tmp["Entropy_display"] = _normalize_0_100(tmp["Entropy"])
    scores = {}
    for g in groups:
        gg = tmp[tmp["Group"].astype(str) == str(g)]
        if gg.empty:
            scores[g] = float("inf")
        else:
            vals = pd.concat([
                pd.to_numeric(gg["Burden_display"], errors="coerce"),
                pd.to_numeric(gg["Entropy_display"], errors="coerce")
            ], ignore_index=True)
            scores[g] = float(vals.mean()) if vals.notna().any() else float("inf")
    return sorted(groups, key=lambda g: (scores.get(g, float("inf")), str(g)))


def plot_burden_entropy_3d(
    summary_df: pd.DataFrame,
    experiment: str,
    out_path: Path,
    burden_kind: str = "BRI",
    treatment_group: Optional[str] = None,
    comparator_group: Optional[str] = None,
) -> None:
    meta = _burden_metadata(burden_kind)
    burden_col = meta["mean"]

    sub = summary_df.copy()
    if "Experiment" in sub.columns:
        sub = sub[sub["Experiment"].astype(str) == str(experiment)].copy()
    keep = [c for c in ["Week", "Group", burden_col, "Entropy"] if c in sub.columns]
    sub = sub[keep].dropna(subset=[burden_col, "Entropy"]).copy()
    if sub.empty:
        return
    sub["Week"] = pd.to_numeric(sub["Week"], errors="coerce").astype(int)
    groups = _ordered_groups_for_plot(sub, treatment_group, comparator_group)
    if _plot_burden_entropy_single_group_heatmap(sub, experiment, out_path, burden_col, meta, treatment_group, comparator_group):
        return
    weeks = sorted(sub["Week"].dropna().astype(int).unique().tolist())
    cmap = _ensure_group_color_map(groups, treatment_group, comparator_group)

    sub["Burden_display"] = _normalize_0_100(sub[burden_col])
    sub["Entropy_display"] = _normalize_0_100(sub["Entropy"])
    # For the default 3D view, put lower combined values in the front row to reduce visual occlusion.
    groups_for_3d = _front_sorted_groups_by_low_values(sub, groups, burden_col)

    fig = plt.figure(figsize=(1.18 * max(6, len(weeks)) + 4.0, 7.4))
    ax = fig.add_subplot(111, projection="3d")

    # Arrange bars like the example: x=Week, y=Group (front / back rows), and within each (week, group)
    # place two neighboring bars for burden and Entropy.
    xgap = 1.20
    ygap = 1.05
    dx = 0.30
    dy = 0.32
    metric_offsets = [-0.20, 0.20]
    y_positions = {g: i * ygap for i, g in enumerate(groups_for_3d)}
    x_positions = {w: i * xgap for i, w in enumerate(weeks)}

    legend_handles, legend_labels = [], []
    for idx_g, group in enumerate(groups):
        c = cmap.get(group, _FALLBACK_COLORS[idx_g % len(_FALLBACK_COLORS)])
        legend_handles.append(_legacy._Line2D([0], [0], marker="s", color="w", markerfacecolor=c,
                                              markeredgecolor="w", markersize=10))
        legend_labels.append(display_group_name(group, idx_g))

    metric_handles = [
        _legacy._Line2D([0], [0], marker="s", color="w", markerfacecolor="#666666", markeredgecolor="w", markersize=10),
        _legacy._Line2D([0], [0], marker="s", color="w", markerfacecolor="#BBBBBB", markeredgecolor="w", markersize=10),
    ]

    for idx_g, group in enumerate(groups_for_3d):
        group_df = sub[sub["Group"].astype(str) == group].set_index("Week").reindex(weeks)
        y0 = y_positions[group]
        base_color = cmap.get(group, _FALLBACK_COLORS[idx_g % len(_FALLBACK_COLORS)])
        ent_color = _lighten(base_color, 0.52)

        for week in weeks:
            row = group_df.loc[week] if week in group_df.index else None
            if row is None or pd.isna(row[burden_col]) or pd.isna(row["Entropy"]):
                continue
            x0 = x_positions[week]
            bh = float(row["Burden_display"])
            eh = float(row["Entropy_display"])

            # Two bars side-by-side along x within the same (week, group) cell
            xb = x0 + metric_offsets[0]
            xe = x0 + metric_offsets[1]
            ax.bar3d(xb, y0, 0, dx, dy, bh, color=base_color, shade=True, alpha=0.62,
                     edgecolor="white", linewidth=0.4)
            ax.bar3d(xe, y0, 0, dx, dy, eh, color=ent_color, shade=True, alpha=0.62,
                     edgecolor="white", linewidth=0.4)

            if len(weeks) <= 10 and len(groups) <= 4:
                ax.text(xb + dx/2, y0 + dy/2, bh + 2.0, f"{row[burden_col]:.2f}", ha="center", va="bottom", fontsize=7.4)
                ax.text(xe + dx/2, y0 + dy/2, eh + 2.0, f"{row['Entropy']:.2f}", ha="center", va="bottom", fontsize=7.4)

    ax.set_xticks([x_positions[w] for w in weeks])
    ax.set_xticklabels([f"{w}W" for w in weeks], fontsize=10)
    ax.set_yticks([y_positions[g] + dy/2 for g in groups_for_3d])
    ax.set_yticklabels([display_group_name(g, groups.index(g) if g in groups else i) for i, g in enumerate(groups_for_3d)], fontsize=10)
    ax.set_xlabel("Week", labelpad=10)
    ax.set_ylabel("Group", labelpad=12)
    ax.set_zlabel("Normalized display height (0–100 within metric)", labelpad=12)
    ax.set_zlim(0, 110)
    ax.view_init(elev=22, azim=-56)

    note = (
        f"Within each week, the front / back rows distinguish groups and the left / right bars represent {meta['short']} and Entropy, respectively. "
        "Bar heights are normalized within metric (0–100) to allow joint display; labels above bars are raw values. "
        "This 3D grouped-bar view is a supplementary chart and can be interactively rotated after export if desired."
    )
    group_legend = ax.legend(legend_handles, legend_labels, loc="upper left", bbox_to_anchor=(0.14, 0.98),
                             ncol=min(max(1, len(legend_labels)), 4), frameon=False)
    ax.add_artist(group_legend)
    ax.legend(handles=metric_handles, labels=[meta["short"], "Entropy"], loc="upper left",
              bbox_to_anchor=(0.72, 0.98), ncol=2, frameon=False)

    _legacy._add_top_title_and_note(
        fig,
        title=f"Experiment {experiment}: {meta['3d_title']}",
        note=note,
        legend_handles=None,
        legend_labels=None,
        title_y=0.985,
        note_y=0.948,
        note_width=126,
        note_fontsize=8.8,
    )
    fig.tight_layout(rect=[0.02, 0.04, 0.98, 0.87])
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _hex_to_rgb_tuple(color: str) -> tuple[int, int, int]:
    try:
        color = str(color).strip()
        if color.startswith('#'):
            color = color[1:]
        if len(color) == 3:
            color = ''.join(ch * 2 for ch in color)
        return tuple(int(color[i:i+2], 16) for i in (0, 2, 4))
    except Exception:
        return (100, 100, 100)


def _rgba(color: str, alpha: float) -> str:
    r, g, b = _hex_to_rgb_tuple(color)
    return f"rgba({r},{g},{b},{float(alpha):.3f})"


def _plotly_add_cuboid(fig, x0, y0, z0, dx, dy, dz, color, opacity, name=None, showlegend=False, hovertext=""):
    # vertices: bottom and top cuboid corners
    x = [x0, x0+dx, x0+dx, x0, x0, x0+dx, x0+dx, x0]
    y = [y0, y0, y0+dy, y0+dy, y0, y0, y0+dy, y0+dy]
    z = [z0, z0, z0, z0, z0+dz, z0+dz, z0+dz, z0+dz]
    # 12 triangles
    i = [0,0,1,1,2,2,3,3,4,4,4,5]
    j = [1,2,5,6,6,7,7,4,5,6,7,6]
    k = [2,3,6,2,7,3,4,0,6,7,0,1]
    fig.add_trace(go.Mesh3d(
        x=x, y=y, z=z, i=i, j=j, k=k,
        color=color, opacity=opacity, flatshading=True,
        lighting=dict(ambient=0.72, diffuse=0.60, specular=0.12, roughness=0.62),
        name=name or "", showlegend=showlegend,
        hovertemplate=hovertext + "<extra></extra>" if hovertext else None,
    ))


def plot_burden_entropy_3d_interactive_html(
    summary_df: pd.DataFrame,
    experiment: str,
    out_path: Path,
    burden_kind: str = "BRI",
    treatment_group: Optional[str] = None,
    comparator_group: Optional[str] = None,
) -> None:
    """Export an interactive Plotly 3D grouped-bar chart.

    The PNG 3D figure is useful for a quick static preview, but the HTML version is more reliable:
    users can rotate/zoom it and inspect each bar interactively.
    """
    if not PLOTLY_AVAILABLE:
        out_path.with_suffix('.plotly_missing.txt').write_text(
            'Plotly is not installed. Install plotly or rebuild with requirements_desktop.txt.', encoding='utf-8'
        )
        return
    meta = _burden_metadata(burden_kind)
    burden_col = meta["mean"]
    sub = summary_df.copy()
    if "Experiment" in sub.columns:
        sub = sub[sub["Experiment"].astype(str) == str(experiment)].copy()
    keep = [c for c in ["Week", "Group", burden_col, "Entropy"] if c in sub.columns]
    sub = sub[keep].dropna(subset=[burden_col, "Entropy"]).copy()
    if sub.empty:
        return
    sub["Week"] = pd.to_numeric(sub["Week"], errors="coerce").astype(int)
    groups = _ordered_groups_for_plot(sub, treatment_group, comparator_group)
    weeks = sorted(sub["Week"].dropna().astype(int).unique().tolist())
    cmap = _ensure_group_color_map(groups, treatment_group, comparator_group)
    if _plot_burden_entropy_single_group_interactive_html(sub, experiment, out_path, burden_col, meta, treatment_group, comparator_group):
        return
    sub["Burden_display"] = _normalize_0_100(sub[burden_col])
    sub["Entropy_display"] = _normalize_0_100(sub["Entropy"])
    # For the default camera, put lower combined values closer to the viewer/front.
    groups_for_3d = _front_sorted_groups_by_low_values(sub, groups, burden_col)

    fig = go.Figure()
    xgap = 1.55
    ygap = 1.10
    dx, dy = 0.34, 0.36
    metric_offsets = [-(dx + 0.12) / 2, +(0.12) / 2]
    x_positions = {w: i * xgap for i, w in enumerate(weeks)}
    y_positions = {g: i * ygap for i, g in enumerate(groups_for_3d)}

    # add dummy legend traces for groups and metrics
    for idx_g, group in enumerate(groups):
        c = cmap.get(group, _FALLBACK_COLORS[idx_g % len(_FALLBACK_COLORS)])
        fig.add_trace(go.Scatter3d(
            x=[None], y=[None], z=[None], mode='markers',
            marker=dict(size=9, color=c, symbol='square'),
            name=display_group_name(group, idx_g), showlegend=True
        ))
    for idx_g, group in enumerate(groups_for_3d):
        gdf = sub[sub["Group"].astype(str) == group].set_index('Week').reindex(weeks)
        base_color = cmap.get(group, _FALLBACK_COLORS[idx_g % len(_FALLBACK_COLORS)])
        ent_color = _lighten(base_color, 0.55)
        y0 = y_positions[group]
        for week in weeks:
            if week not in gdf.index:
                continue
            row = gdf.loc[week]
            if pd.isna(row[burden_col]) or pd.isna(row['Entropy']):
                continue
            x0 = x_positions[week]
            bh = float(row['Burden_display'])
            eh = float(row['Entropy_display'])
            raw_b = float(row[burden_col])
            raw_e = float(row['Entropy'])
            xb = x0 + metric_offsets[0]
            xe = x0 + metric_offsets[1]
            hover_b = f"Experiment: {experiment}<br>Week: {week}W<br>Group: {display_group_name(group, groups.index(group) if group in groups else idx_g)}<br>Metric: {meta['short']}<br>Raw value: {raw_b:.4g}<br>Display height: {bh:.1f}"
            hover_e = f"Experiment: {experiment}<br>Week: {week}W<br>Group: {display_group_name(group, groups.index(group) if group in groups else idx_g)}<br>Metric: Entropy<br>Raw value: {raw_e:.4g}<br>Display height: {eh:.1f}"
            _plotly_add_cuboid(fig, xb, y0, 0, dx, dy, bh, _rgba(base_color, 0.72), 0.72, hovertext=hover_b)
            _plotly_add_cuboid(fig, xe, y0, 0, dx, dy, eh, _rgba(ent_color, 0.55), 0.55, hovertext=hover_e)
            fig.add_trace(go.Scatter3d(x=[xb+dx/2], y=[y0+dy/2], z=[bh+3], text=[f"{raw_b:.2f}"], mode='text', showlegend=False, textfont=dict(size=10, color='#222')))
            fig.add_trace(go.Scatter3d(x=[xe+dx/2], y=[y0+dy/2], z=[eh+3], text=[f"{raw_e:.2f}"], mode='text', showlegend=False, textfont=dict(size=10, color='#333')))

    fig.update_layout(
        title=f"Experiment {experiment}: {meta['3d_title']} (interactive)",
        scene=dict(
            xaxis=dict(title='Week', tickmode='array', tickvals=[x_positions[w] for w in weeks], ticktext=[f"{w}W" for w in weeks]),
            yaxis=dict(title='Group', tickmode='array', tickvals=[y_positions[g]+dy/2 for g in groups_for_3d], ticktext=[display_group_name(g, groups.index(g) if g in groups else i) for i, g in enumerate(groups_for_3d)]),
            zaxis=dict(title='Normalized display height (0–100 within metric)', range=[0, 112]),
            camera=dict(eye=dict(x=1.55, y=-2.05, z=1.35)),
            aspectmode='manual', aspectratio=dict(x=max(1.8, len(weeks)*0.45), y=max(0.8, len(groups_for_3d)*0.36), z=0.9),
        ),
        legend=dict(orientation='h', y=1.02, x=0.5, xanchor='center', title=dict(text='Group; dark=BRI/NMD, light=Entropy')),
        margin=dict(l=0, r=0, t=95, b=0),
        height=760,
        annotations=[dict(
            text=f"Each (Week × Group) cell has two bars: dark group color = {meta['short']}, light group color = Entropy. The default view places lower combined values in the front row to reduce occlusion. Drag to rotate; scroll to zoom. Heights are normalized separately within each metric; labels and hover tooltips show raw values.",
            x=0.5, y=0.965, xref='paper', yref='paper', showarrow=False,
            font=dict(size=13, color='#555'), align='center'
        )]
    )
    fig.write_html(str(out_path), include_plotlyjs='cdn', full_html=True, div_id=_stable_html_div_id(out_path))


def plot_burden_entropy_3d_perspective_html(
    summary_df: pd.DataFrame,
    experiment: str,
    out_path: Path,
    burden_kind: str = "NMD",
    treatment_group: Optional[str] = None,
    comparator_group: Optional[str] = None,
) -> None:
    if not PLOTLY_AVAILABLE:
        out_path.with_suffix(".plotly_missing.txt").write_text(
            "Plotly is not installed. Rebuild with requirements_desktop.txt.", encoding="utf-8"
        )
        return
    meta = _burden_metadata(burden_kind)
    burden_col = meta["mean"]
    if burden_col not in summary_df.columns or "Entropy" not in summary_df.columns:
        return
    sub = summary_df.copy()
    if "Experiment" in sub.columns:
        sub = sub[sub["Experiment"].astype(str) == str(experiment)].copy()
    sub = sub[[c for c in ["Week", "Group", burden_col, "Entropy"] if c in sub.columns]].dropna(subset=[burden_col, "Entropy"]).copy()
    if sub.empty:
        return
    sub["Week"] = pd.to_numeric(sub["Week"], errors="coerce").astype(int)
    single_timepoint = _is_single_timepoint(sub)
    groups = _ordered_single_timepoint_groups(sub, burden_col) if single_timepoint else _ordered_groups_for_plot(sub, treatment_group, comparator_group)
    weeks = sorted(sub["Week"].dropna().astype(int).unique().tolist())
    sub["Burden_display"] = _normalize_display_window(sub[burden_col], sub.get("Group"))
    sub["Entropy_display"] = _normalize_display_window(sub["Entropy"], sub.get("Group"))

    fig = go.Figure()
    burden_color = _NMD_COLOR if str(meta.get("short", "")).upper() == "NMD" else _BRI_COLOR
    entropy_color = _ENTROPY_COLOR
    dx, dy = 0.36, 0.38

    # Legend proxies.
    fig.add_trace(go.Scatter3d(
        x=[None], y=[None], z=[None], mode="markers",
        marker=dict(size=10, color=burden_color, symbol="square"),
        name=meta["short"], showlegend=True,
    ))
    fig.add_trace(go.Scatter3d(
        x=[None], y=[None], z=[None], mode="markers",
        marker=dict(size=10, color=entropy_color, symbol="square"),
        name="Entropy", showlegend=True,
    ))

    if single_timepoint:
        x_positions = {g: i * 1.25 for i, g in enumerate(groups)}
        y_positions = {meta["short"]: 0.0, "Entropy": 0.92}
        plot_df = _metric_group_rows(sub, groups)
        for idx, group in enumerate(groups):
            row_df = plot_df[plot_df["Group"].astype(str) == str(group)]
            if row_df.empty:
                continue
            row = row_df.iloc[0]
            x0 = x_positions[group]
            for metric_name, disp_col, raw_col, y0, color in [
                (meta["short"], "Burden_display", burden_col, y_positions[meta["short"]], burden_color),
                ("Entropy", "Entropy_display", "Entropy", y_positions["Entropy"], entropy_color),
            ]:
                h = float(row[disp_col])
                raw = float(row[raw_col])
                hover = (
                    f"Experiment: {experiment}<br>Group: {display_group_name(group, idx)}"
                    f"<br>Metric: {metric_name}<br>Raw value: {raw:.4g}<br>Display height: {h:.1f}"
                )
                _plotly_add_cuboid(fig, x0, y0, 0, dx, dy, h, _rgba(color, 0.82), 0.82, hovertext=hover)
                fig.add_trace(go.Scatter3d(
                    x=[x0 + dx / 2], y=[y0 + dy / 2], z=[h + 3],
                    text=[f"{raw:.2f}"], mode="text", showlegend=False,
                    textfont=dict(size=10, color="#1f2933"),
                    hoverinfo="skip",
                ))
        scene = dict(
            xaxis=dict(title="Treatment group", tickmode="array",
                       tickvals=[x_positions[g] + dx / 2 for g in groups],
                       ticktext=[display_group_name(g, i) for i, g in enumerate(groups)]),
            yaxis=dict(title="Metric", tickmode="array",
                       tickvals=[y_positions[meta["short"]] + dy / 2, y_positions["Entropy"] + dy / 2],
                       ticktext=[meta["short"], "Entropy"]),
            zaxis=dict(title="Display height (0-100)", range=[0, 112]),
            camera=dict(eye=dict(x=1.55, y=-2.05, z=1.35), projection=dict(type="perspective")),
            aspectmode="manual",
            aspectratio=dict(x=max(1.8, len(groups) * 0.33), y=0.9, z=0.85),
        )
        subtitle = "Single-timepoint experiment: x-axis is disease/treatment group. Drag to rotate; scroll to zoom."
    else:
        x_positions = {w: i * 1.35 for i, w in enumerate(weeks)}
        y_positions = {g: i * 1.02 for i, g in enumerate(groups)}
        metric_offsets = [-0.20, 0.22]
        for idx, group in enumerate(groups):
            gdf = sub[sub["Group"].astype(str) == str(group)].set_index("Week").reindex(weeks)
            for week in weeks:
                if week not in gdf.index:
                    continue
                row = gdf.loc[week]
                if pd.isna(row[burden_col]) or pd.isna(row["Entropy"]):
                    continue
                x_base = x_positions[week]
                y_base = y_positions[group]
                for metric_name, disp_col, raw_col, offset, color in [
                    (meta["short"], "Burden_display", burden_col, metric_offsets[0], burden_color),
                    ("Entropy", "Entropy_display", "Entropy", metric_offsets[1], entropy_color),
                ]:
                    h = float(row[disp_col])
                    raw = float(row[raw_col])
                    hover = (
                        f"Experiment: {experiment}<br>Week: {week}W<br>Group: {display_group_name(group, idx)}"
                        f"<br>Metric: {metric_name}<br>Raw value: {raw:.4g}<br>Display height: {h:.1f}"
                    )
                    _plotly_add_cuboid(fig, x_base + offset, y_base, 0, dx, dy, h, _rgba(color, 0.80), 0.80, hovertext=hover)
        scene = dict(
            xaxis=dict(title="Week", tickmode="array",
                       tickvals=[x_positions[w] for w in weeks],
                       ticktext=[f"{w}W" for w in weeks]),
            yaxis=dict(title="Group", tickmode="array",
                       tickvals=[y_positions[g] + dy / 2 for g in groups],
                       ticktext=[display_group_name(g, i) for i, g in enumerate(groups)]),
            zaxis=dict(title="Display height (0-100)", range=[0, 112]),
            camera=dict(eye=dict(x=1.55, y=-2.20, z=1.35), projection=dict(type="perspective")),
            aspectmode="manual",
            aspectratio=dict(x=max(1.8, len(weeks) * 0.42), y=max(0.9, len(groups) * 0.34), z=0.85),
        )
        subtitle = "Each Week x Group cell has two bars. Drag to rotate; scroll to zoom."

    fig.update_layout(
        title=f"Experiment {experiment}: interactive {meta['short']} x Entropy 3D bars",
        template="plotly_white",
        scene=scene,
        legend=dict(orientation="h", y=1.02, x=0.5, xanchor="center"),
        margin=dict(l=10, r=10, t=108, b=15),
        height=760,
        annotations=[dict(
            text=(
                f"{subtitle} Heights use anchored 0-100 display normalization when Normal and model controls are present; "
                f"gamma={DISPLAY_CONTRAST_GAMMA:g} display contrast spreads high disease-burden values; hover and labels show raw values."
            ),
            x=0.5, y=0.965, xref="paper", yref="paper", showarrow=False,
            font=dict(size=13, color="#4B5563"), align="center",
        )],
    )
    fig.write_html(str(out_path), include_plotlyjs="cdn", full_html=True, div_id=_stable_html_div_id(out_path))


def _preferred_burden_kind(exp_res: Dict[str, Any]) -> str:
    summary_df = exp_res.get("summary", pd.DataFrame())
    has_summary = isinstance(summary_df, pd.DataFrame) and not summary_df.empty
    # 1. True enrolled Normal/healthy reference cloud → Mahalanobis NMD.
    if bool(exp_res.get("true_normal_available")) and has_summary and "NMD_mean" in summary_df.columns:
        return "NMD"
    # 2. Published / uploaded Normal-range anchor → absolute NRBS (preferred over
    #    the purely baseline-relative BRI for an absolute health scale).
    if has_summary and "NRBS_mean" in summary_df.columns and pd.to_numeric(summary_df["NRBS_mean"], errors="coerce").notna().any():
        return "NRBS"
    # 3. No absolute anchor available → baseline-relative BRI.
    return "BRI"


def _normalize_within_experiment(values: pd.Series) -> pd.Series:
    s = pd.to_numeric(values, errors="coerce").astype(float)
    finite = s[np.isfinite(s)]
    if finite.empty:
        return pd.Series(np.nan, index=s.index, dtype=float)
    vmin = float(finite.min())
    vmax = float(finite.max())
    if math.isclose(vmin, vmax):
        return pd.Series(50.0, index=s.index, dtype=float)
    return (s - vmin) / (vmax - vmin) * 100.0


def _is_healthy_reference_group(group: Any) -> bool:
    name = str(group).strip().lower()
    return any(tok in name for tok in ["normal", "healthy", "sham", "blank", "naive", "正常", "健康", "空白", "假手术"])


def _is_model_control_group(group: Any) -> bool:
    name = str(group).strip().lower()
    return any(tok in name for tok in ["model control", "model", "disease control", "disease", "vehicle", "模型", "疾病", "模型对照"])


def _normalize_anchor_linear_0_100(values: pd.Series, groups: Optional[pd.Series] = None) -> pd.Series:
    """Linear anchored normalization for traceable table values.

    Prefer biological anchors when both are available:
    Normal/healthy reference = 0 and model/disease control = 100.
    This does not change raw values or statistical inference.
    """
    s = pd.to_numeric(values, errors="coerce").astype(float)
    finite_mask = np.isfinite(s)
    ref = s[finite_mask]
    if groups is not None:
        g = pd.Series(groups, index=s.index).astype(str)
        normal_vals = s[finite_mask & g.map(_is_healthy_reference_group)]
        model_vals = s[finite_mask & g.map(_is_model_control_group)]
        if normal_vals.notna().any() and model_vals.notna().any():
            normal_anchor = float(normal_vals.mean())
            model_anchor = float(model_vals.mean())
            if not math.isclose(normal_anchor, model_anchor):
                return ((s - normal_anchor) / (model_anchor - normal_anchor) * 100.0).clip(lower=0.0, upper=100.0)
    if ref.empty:
        return pd.Series(np.nan, index=s.index, dtype=float)
    vmin = float(ref.min())
    vmax = float(ref.max())
    if math.isclose(vmin, vmax):
        return pd.Series(50.0, index=s.index, dtype=float)
    return ((s - vmin) / (vmax - vmin) * 100.0).clip(lower=0.0, upper=100.0)


def _contrast_enhance_0_100(linear_values: pd.Series, gamma: float = DISPLAY_CONTRAST_GAMMA) -> pd.Series:
    """Nonlinear display transform that preserves 0 and 100 anchors.

    Gamma > 1 stretches the upper disease-range visually. It is display-only; raw values and
    linear anchored values are retained in tables and labels.
    """
    s = pd.to_numeric(linear_values, errors="coerce").astype(float).clip(lower=0.0, upper=100.0)
    return ((s / 100.0) ** float(gamma) * 100.0).clip(lower=0.0, upper=100.0)


def _normalize_display_window(values: pd.Series, groups: Optional[pd.Series] = None) -> pd.Series:
    return _contrast_enhance_0_100(_normalize_anchor_linear_0_100(values, groups))


def build_display_normalized_metric_table(results: Dict[str, Any]) -> pd.DataFrame:
    rows: List[pd.DataFrame] = []
    for exp, exp_res in results.get("experiments", {}).items():
        summary_df = exp_res.get("summary", pd.DataFrame())
        if not isinstance(summary_df, pd.DataFrame) or summary_df.empty:
            continue
        burden_kind = _preferred_burden_kind(exp_res)
        meta = _burden_metadata(burden_kind)
        burden_col = meta["mean"]
        if burden_col not in summary_df.columns:
            continue
        keep = [c for c in ["Experiment", "Week", "Group", "n", burden_col, "Entropy"] if c in summary_df.columns]
        sub = summary_df[keep].dropna(subset=[burden_col]).copy()
        if sub.empty:
            continue
        sub["Display_metric"] = meta["short"]
        sub["Raw_disease_deviation"] = pd.to_numeric(sub[burden_col], errors="coerce")
        burden_linear = _normalize_anchor_linear_0_100(sub[burden_col], sub.get("Group"))
        sub["Disease_deviation_linear_0_100"] = burden_linear.round(4)
        sub["Disease_deviation_display_contrast_0_100"] = _contrast_enhance_0_100(burden_linear).round(4)
        sub["Disease_deviation_normalized_0_100"] = sub["Disease_deviation_display_contrast_0_100"]
        if "Entropy" in sub.columns:
            entropy_linear = _normalize_anchor_linear_0_100(sub["Entropy"], sub.get("Group"))
            sub["Entropy_linear_0_100"] = entropy_linear.round(4)
            sub["Entropy_display_contrast_0_100"] = _contrast_enhance_0_100(entropy_linear).round(4)
            sub["Entropy_normalized_0_100"] = sub["Entropy_display_contrast_0_100"]
        sub["Normalization_note"] = (
            "Linear anchored values set Normal/healthy control = 0 and model/disease control = 100 when available. "
            f"Display_contrast values apply gamma={DISPLAY_CONTRAST_GAMMA:g} to spread high disease-deviation and entropy values visually. "
            "Use raw values, linear anchored values, and statistical tests for formal inference."
        )
        if _custom_order_for_frame(sub):
            ordered_groups = _apply_custom_group_order(sub["Group"].dropna().astype(str).unique().tolist(), sub)
            rank_map = {str(g): i for i, g in enumerate(ordered_groups)}
            sub["_Group_order_rank"] = sub["Group"].astype(str).map(rank_map).fillna(9999).astype(int)
        rows.append(sub)
    if not rows:
        return pd.DataFrame()
    out = pd.concat(rows, ignore_index=True)
    sort_cols = [c for c in ["Experiment", "Week", "_Group_order_rank", "Disease_deviation_normalized_0_100", "Group"] if c in out.columns]
    out = out.sort_values(sort_cols).reset_index(drop=True)
    return out.drop(columns=[c for c in ["_Group_order_rank"] if c in out.columns])


def plot_normalized_metric_index(
    summary_df: pd.DataFrame,
    experiment: str,
    out_path: Path,
    metric_col: str,
    metric_label: str,
    title_label: str,
    raw_label: str,
    treatment_group: Optional[str] = None,
    comparator_group: Optional[str] = None,
) -> None:
    if metric_col not in summary_df.columns:
        return
    sub = summary_df.copy()
    if "Experiment" in sub.columns:
        sub = sub[sub["Experiment"].astype(str) == str(experiment)].copy()
    keep = [c for c in ["Week", "Group", metric_col] if c in sub.columns]
    sub = sub[keep].dropna(subset=[metric_col]).copy()
    if sub.empty:
        return
    sub["Week"] = pd.to_numeric(sub["Week"], errors="coerce").astype(int)
    sub["Metric_linear"] = _normalize_anchor_linear_0_100(sub[metric_col], sub.get("Group"))
    sub["Metric_norm"] = _contrast_enhance_0_100(sub["Metric_linear"])
    single_timepoint = _is_single_timepoint(sub)
    cmap = _ensure_group_color_map(sub["Group"].dropna().astype(str).unique().tolist(), treatment_group, comparator_group)

    if single_timepoint:
        groups = _ordered_single_timepoint_groups(sub, metric_col)
        plot_df = _metric_group_rows(sub, groups)
        if plot_df.empty:
            return
        x = np.arange(len(groups), dtype=float)
        fig, ax = plt.subplots(figsize=(max(8.8, 1.05 * len(groups) + 3.8), 5.9))
        colors = [cmap.get(g, _FALLBACK_COLORS[i % len(_FALLBACK_COLORS)]) for i, g in enumerate(groups)]
        bars = ax.bar(x, plot_df["Metric_norm"], color=colors, alpha=0.92, width=0.62,
                      edgecolor="#FFFFFF", linewidth=0.8)
        for bar, raw, norm in zip(bars, plot_df[metric_col], plot_df["Metric_norm"]):
            if pd.notna(raw) and pd.notna(norm):
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 2.0,
                        f"{float(norm):.1f}\nraw {float(raw):.2f}",
                        ha="center", va="bottom", fontsize=7.0, color=_TEXT_COLOR)
        # Overlay each subject's value on the same 0–100 display scale.
        _overlay_group_points(ax, x, groups,
                              _per_sample_points_for_metric(experiment, metric_col, groups, normalize="display"))
        ax.set_xticks(x)
        ax.set_xticklabels([display_group_name(g, i) for i, g in enumerate(groups)], rotation=24, ha="right")
        ax.set_xlabel("Treatment group")
        note = (
            f"Single-timepoint experiment: {metric_label} is normalized for display. Linear anchors set Normal=0 and model=100 when available, "
            f"otherwise fall back to experiment min-max; plotted bars use gamma={DISPLAY_CONTRAST_GAMMA:g} therapeutic-window contrast. Labels show display and raw values; linear values remain in Excel."
        )
    else:
        weeks, groups = _get_weeks_and_groups(sub, treatment_group, comparator_group)
        fig, ax = plt.subplots(figsize=(9.4, 5.6))
        legend_handles, legend_labels = [], []
        for idx, group in enumerate(groups):
            g = sub[sub["Group"].astype(str) == str(group)].sort_values("Week")
            if g.empty:
                continue
            color = cmap.get(group, _FALLBACK_COLORS[idx % len(_FALLBACK_COLORS)])
            ax.plot(g["Week"], g["Metric_norm"], marker="o", linewidth=2.25, color=color,
                    markeredgecolor="white", markeredgewidth=0.8, zorder=4)
            legend_handles.append(_legacy._Line2D([0], [0], color=color, marker="o", linewidth=2.2))
            legend_labels.append(display_group_name(group, idx))
            for _, row in g.iterrows():
                ax.annotate(f"{float(row['Metric_norm']):.0f}", (row["Week"], row["Metric_norm"]),
                            xytext=(0, 7 if idx == 0 else -7), textcoords="offset points",
                            ha="center", va="bottom" if idx == 0 else "top", fontsize=7.2, color=color, zorder=7)
        # Overlay each subject's value on the same 0–100 display scale.
        _overlay_trajectory_points(ax, experiment, metric_col, weeks, groups, cmap, normalize="display")
        ax.set_xticks(weeks)
        ax.set_xticklabels([f"{w}W" for w in weeks])
        ax.set_xlabel("Week")
        note = (
            f"{metric_label} is normalized for display. Linear anchors set Normal=0 and model=100 when available, otherwise fall back to experiment min-max; "
            f"plotted values use gamma={DISPLAY_CONTRAST_GAMMA:g} therapeutic-window contrast. Linear values are retained in Excel."
        )
    ax.set_ylabel(f"{raw_label}\ndisplay index (0-100)")
    ax.set_ylim(-4, 112)
    _legacy._format_axis(ax)
    _apply_academic_axis(ax, "y")
    _legacy._add_top_title_and_note(
        fig,
        title=f"Experiment {experiment}: normalized {title_label}",
        note=note,
        legend_handles=legend_handles if not single_timepoint else None,
        legend_labels=legend_labels if not single_timepoint else None,
        legend_ncol=min(len(legend_labels or []), 4) if not single_timepoint and legend_labels else 1,
        legend_y=0.875,
        title_y=0.985,
        note_y=0.946,
        note_width=122,
        note_fontsize=8.7,
    )
    fig.tight_layout(rect=[0.06, 0.10, 0.98, 0.80 if single_timepoint else 0.82])
    fig.savefig(out_path, dpi=330, bbox_inches="tight")
    plt.close(fig)


def plot_normalized_burden_index(
    summary_df: pd.DataFrame,
    experiment: str,
    out_path: Path,
    burden_kind: str = "NMD",
    treatment_group: Optional[str] = None,
    comparator_group: Optional[str] = None,
) -> None:
    meta = _burden_metadata(burden_kind)
    plot_normalized_metric_index(
        summary_df, experiment, out_path,
        metric_col=meta["mean"],
        metric_label=f"{meta['short']} / disease deviation",
        title_label=f"{meta['short']} disease-deviation display",
        raw_label=meta["short"],
        treatment_group=treatment_group,
        comparator_group=comparator_group,
    )


def plot_normalized_entropy_index(
    summary_df: pd.DataFrame,
    experiment: str,
    out_path: Path,
    treatment_group: Optional[str] = None,
    comparator_group: Optional[str] = None,
) -> None:
    plot_normalized_metric_index(
        summary_df, experiment, out_path,
        metric_col="Entropy",
        metric_label="Entropy",
        title_label="Entropy display",
        raw_label="Entropy",
        treatment_group=treatment_group,
        comparator_group=comparator_group,
    )


def _copy_or_regenerate_main_metric_plot(
    summary_df: pd.DataFrame,
    experiment: str,
    out_path: Path,
    burden_kind: str,
    comparison_df: Optional[pd.DataFrame],
    treatment_group: Optional[str],
    comparator_group: Optional[str],
) -> None:
    if str(burden_kind).upper() == "NMD" and "NMD_mean" in summary_df.columns:
        plot_true_normal_nmd_trajectory(summary_df, experiment, out_path, comparison_df, treatment_group, comparator_group)
    elif str(burden_kind).upper() == "NRBS" and "NRBS_mean" in summary_df.columns:
        plot_normal_range_burden_trajectory(summary_df, experiment, out_path, comparison_df, treatment_group, comparator_group)
    else:
        plot_relative_burden_trajectory(summary_df, experiment, out_path, comparison_df, treatment_group, comparator_group)


def plot_burden_entropy_main_3d(
    summary_df: pd.DataFrame,
    experiment: str,
    out_path: Path,
    burden_kind: str = "NMD",
    treatment_group: Optional[str] = None,
    comparator_group: Optional[str] = None,
) -> None:
    meta = _burden_metadata(burden_kind)
    burden_col = meta["mean"]
    if burden_col not in summary_df.columns or "Entropy" not in summary_df.columns:
        return
    sub = summary_df.copy()
    if "Experiment" in sub.columns:
        sub = sub[sub["Experiment"].astype(str) == str(experiment)].copy()
    keep = [c for c in ["Week", "Group", burden_col, "Entropy"] if c in sub.columns]
    sub = sub[keep].dropna(subset=[burden_col, "Entropy"]).copy()
    if sub.empty:
        return
    sub["Week"] = pd.to_numeric(sub["Week"], errors="coerce").astype(int)
    single_timepoint = _is_single_timepoint(sub)
    groups = (
        _ordered_single_timepoint_groups(sub, burden_col)
        if single_timepoint else
        _ordered_groups_for_plot(sub, treatment_group, comparator_group)
    )
    weeks = sorted(sub["Week"].dropna().astype(int).unique().tolist())
    sub["Burden_display"] = _normalize_0_100(sub[burden_col])
    sub["Entropy_display"] = _normalize_0_100(sub["Entropy"])
    burden_color = _NMD_COLOR if str(meta.get("short", "")).upper() == "NMD" else _BRI_COLOR

    fig = plt.figure(figsize=(max(10.4, 1.08 * max(len(groups), len(weeks)) + 7.0), 7.4))
    ax = fig.add_subplot(111, projection="3d")
    dx, dy = 0.34, 0.34

    if single_timepoint:
        plot_df = _metric_group_rows(sub, groups)
        x_positions = {g: i * 1.05 for i, g in enumerate(groups)}
        metric_y = {meta["short"]: 0.0, "Entropy": 0.72}
        for idx, group in enumerate(groups):
            row = plot_df[plot_df["Group"].astype(str) == str(group)]
            if row.empty:
                continue
            row = row.iloc[0]
            x0 = x_positions[group]
            for label, col, color in [
                (meta["short"], "Burden_display", burden_color),
                ("Entropy", "Entropy_display", _ENTROPY_COLOR),
            ]:
                h = float(row[col])
                raw_col = burden_col if label == meta["short"] else "Entropy"
                ax.bar3d(x0, metric_y[label], 0, dx, dy, h, color=color, alpha=0.88,
                         edgecolor="white", linewidth=0.45, shade=True)
                if len(groups) <= 9:
                    ax.text(x0 + dx / 2, metric_y[label] + dy / 2, h + 2.0, f"{float(row[raw_col]):.2f}",
                            ha="center", va="bottom", fontsize=7.2, color=_TEXT_COLOR)
        ax.set_xticks([x_positions[g] + dx / 2 for g in groups])
        ax.set_xticklabels([display_group_name(g, i) for i, g in enumerate(groups)], rotation=18, ha="right", fontsize=9)
        ax.set_xlabel("Treatment group", labelpad=26)
        ax.set_yticks([metric_y[meta["short"]] + dy / 2, metric_y["Entropy"] + dy / 2])
        ax.set_yticklabels([meta["short"], "Entropy"], fontsize=10)
        note = (
            f"Single-timepoint experiment: x-axis is treatment group. Bar heights are normalized separately for {meta['short']} "
            "and Entropy (0-100); labels show raw values."
        )
    else:
        x_positions = {w: i * 1.25 for i, w in enumerate(weeks)}
        y_positions = {g: i * 0.96 for i, g in enumerate(groups)}
        metric_offsets = [-0.19, 0.20]
        for idx_g, group in enumerate(groups):
            gdf = sub[sub["Group"].astype(str) == str(group)].set_index("Week").reindex(weeks)
            for week in weeks:
                if week not in gdf.index:
                    continue
                row = gdf.loc[week]
                if pd.isna(row[burden_col]) or pd.isna(row["Entropy"]):
                    continue
                x0 = x_positions[week]
                y0 = y_positions[group]
                for offset, label, raw_col, disp_col, color in [
                    (metric_offsets[0], meta["short"], burden_col, "Burden_display", burden_color),
                    (metric_offsets[1], "Entropy", "Entropy", "Entropy_display", _ENTROPY_COLOR),
                ]:
                    h = float(row[disp_col])
                    ax.bar3d(x0 + offset, y0, 0, dx, dy, h, color=color, alpha=0.86,
                             edgecolor="white", linewidth=0.42, shade=True)
                    if len(groups) * len(weeks) <= 24:
                        ax.text(x0 + offset + dx / 2, y0 + dy / 2, h + 2.0, f"{float(row[raw_col]):.2f}",
                                ha="center", va="bottom", fontsize=6.8, color=_TEXT_COLOR)
        ax.set_xticks([x_positions[w] for w in weeks])
        ax.set_xticklabels([f"{w}W" for w in weeks], fontsize=9)
        ax.set_xlabel("Week", labelpad=10)
        ax.set_yticks([y_positions[g] + dy / 2 for g in groups])
        ax.set_yticklabels([display_group_name(g, i) for i, g in enumerate(groups)], fontsize=9)
        ax.set_ylabel("Group", labelpad=12)
        note = (
            f"Each Week x Group cell contains two bars: {meta['short']} and Entropy. Heights are normalized separately "
            "within each metric (0-100); labels show raw values when space allows."
        )

    ax.set_zlim(0, 112)
    ax.set_zlabel("Normalized display height (0-100)", labelpad=10)
    ax.view_init(elev=24, azim=-55)
    for pane in [ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane]:
        pane.set_edgecolor("#EDF2F7")
        pane.set_facecolor((0.98, 0.99, 1.0, 0.82))
    ax.tick_params(colors=_TEXT_COLOR)
    legend_handles = [
        Patch(facecolor=burden_color, edgecolor="none", label=meta["short"], alpha=0.88),
        Patch(facecolor=_ENTROPY_COLOR, edgecolor="none", label="Entropy", alpha=0.88),
    ]
    _legacy._add_top_title_and_note(
        fig,
        title=f"Experiment {experiment}: {meta['short']} x Entropy 3D bar chart",
        note=note,
        legend_handles=legend_handles,
        legend_labels=[meta["short"], "Entropy"],
        legend_ncol=2,
        legend_y=0.885,
        title_y=0.985,
        note_y=0.948,
        note_width=118,
        note_fontsize=8.8,
    )
    fig.tight_layout(rect=[0.02, 0.04, 0.98, 0.84])
    fig.savefig(out_path, dpi=320, bbox_inches="tight")
    plt.close(fig)


def plot_mechanism_diagram(
    summary_df: pd.DataFrame,
    experiment: str,
    out_path: Path,
    burden_kind: str = "NMD",
) -> None:
    meta = _burden_metadata(burden_kind)
    fig, ax = plt.subplots(figsize=(12.2, 6.8))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    boxes = [
        (0.05, 0.58, 0.22, 0.22, "Disease state", "Biochemical markers drift\naway from the Normal cloud", _MODEL_COLOR),
        (0.39, 0.58, 0.22, 0.22, f"{meta['short']} rises", "Distance from healthy\nreference increases", _NMD_COLOR if meta["short"] == "NMD" else _BRI_COLOR),
        (0.73, 0.58, 0.22, 0.22, "Entropy rises", "Multi-marker covariance\nbecomes more dispersed", _ENTROPY_COLOR),
        (0.39, 0.17, 0.22, 0.22, "Effective treatment", "State shifts toward\nlower burden", "#2F855A"),
        (0.73, 0.17, 0.22, 0.22, "System re-ordering", f"{meta['short']} and Entropy\nmove downward together", "#2C7A7B"),
    ]
    for x, y, w, h, title, body, color in boxes:
        rect = FancyBboxPatch(
            (x, y), w, h,
            boxstyle="round,pad=0.018,rounding_size=0.018",
            linewidth=1.35,
            edgecolor=color,
            facecolor=_lighten(color, 0.86),
        )
        ax.add_patch(rect)
        ax.text(x + w / 2, y + h - 0.052, title, ha="center", va="center",
                fontsize=13.0, fontweight="bold", color=color)
        ax.text(x + w / 2, y + h / 2 - 0.025, body, ha="center", va="center",
                fontsize=10.2, color=_TEXT_COLOR, linespacing=1.25)

    arrows = [
        ((0.27, 0.69), (0.39, 0.69), _MODEL_COLOR),
        ((0.61, 0.69), (0.73, 0.69), _NMD_COLOR),
        ((0.50, 0.58), (0.50, 0.39), "#2F855A"),
        ((0.61, 0.28), (0.73, 0.28), "#2C7A7B"),
        ((0.84, 0.58), (0.84, 0.39), _ENTROPY_COLOR),
    ]
    for start, end, color in arrows:
        ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=18,
                                     linewidth=1.8, color=color, shrinkA=4, shrinkB=4))

    ax.text(0.50, 0.925, f"Experiment {experiment}: BioEntropy working hypothesis",
            ha="center", va="center", fontsize=16.2, fontweight="bold", color=_TEXT_COLOR)
    ax.text(
        0.50, 0.875,
        f"Therapeutic improvement should reduce both {meta['short']} (disease deviation) and Entropy (multi-marker disorder).",
        ha="center", va="center", fontsize=10.5, color="#4A5568"
    )
    ax.text(0.055, 0.065,
            "Interpretation boundary: this diagram is a mechanistic hypothesis map, not an independent statistical test.",
            ha="left", va="center", fontsize=9.0, color="#697386")
    fig.savefig(out_path, dpi=330, bbox_inches="tight")
    plt.close(fig)


def _plot_feature_mechanism_unavailable(experiment: str, out_path: Path, reason: str) -> None:
    fig, ax = plt.subplots(figsize=(10.5, 5.6))
    ax.axis("off")
    ax.set_facecolor(_LIGHT_BG)
    message = (
        "Feature-level mechanism driver plot could not be generated for this experiment.\n\n"
        f"Reason: {reason}\n\n"
        "This figure is intended to identify biochemical indicators whose treatment-vs-comparator changes "
        "may jointly explain lower disease deviation and lower entropy. Generate or retain the feature comparison "
        "table to support this mechanism-level visualization."
    )
    ax.text(
        0.5, 0.52, message,
        ha="center", va="center", transform=ax.transAxes,
        fontsize=11.0, color=_TEXT_COLOR, linespacing=1.55,
        bbox=dict(boxstyle="round,pad=0.65", facecolor="white", edgecolor=_GRID_COLOR, linewidth=1.0),
    )
    _legacy._add_top_title_and_note(
        fig,
        title=f"Experiment {experiment}: feature mechanism driver check",
        note="This is a mechanism-data availability notice, not a software-mechanism diagram.",
        legend_handles=None,
        legend_labels=None,
        title_y=0.96,
        note_y=0.905,
        note_width=110,
        note_fontsize=9.0,
    )
    fig.tight_layout(rect=[0.04, 0.04, 0.96, 0.84])
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_feature_mechanism_drivers(
    exp_res: Dict[str, Any],
    experiment: str,
    out_path: Path,
    burden_kind: str = "NMD",
) -> None:
    feature_cmp = exp_res.get("feature_comparisons", pd.DataFrame())
    if not isinstance(feature_cmp, pd.DataFrame) or feature_cmp.empty:
        _plot_feature_mechanism_unavailable(experiment, out_path, "feature_comparisons table is empty or unavailable")
        return
    cmp_df = feature_cmp.copy()
    if "Experiment" in cmp_df.columns:
        cmp_df = cmp_df[cmp_df["Experiment"].astype(str) == str(experiment)].copy()
    if cmp_df.empty:
        _plot_feature_mechanism_unavailable(experiment, out_path, "no feature comparison rows match this experiment")
        return
    if "Week" in cmp_df.columns and cmp_df["Week"].notna().any():
        terminal_week = int(pd.to_numeric(cmp_df["Week"], errors="coerce").max())
        cmp_df = cmp_df[pd.to_numeric(cmp_df["Week"], errors="coerce") == terminal_week].copy()
    else:
        terminal_week = None
    needed = ["Feature", "Burden_diff_treat_minus_comp", "VarLog_diff_treat_minus_comp"]
    if not all(c in cmp_df.columns for c in needed):
        missing = [c for c in needed if c not in cmp_df.columns]
        _plot_feature_mechanism_unavailable(experiment, out_path, "missing required columns: " + ", ".join(missing))
        return
    cmp_df["Burden_diff_treat_minus_comp"] = pd.to_numeric(cmp_df["Burden_diff_treat_minus_comp"], errors="coerce")
    cmp_df["VarLog_diff_treat_minus_comp"] = pd.to_numeric(cmp_df["VarLog_diff_treat_minus_comp"], errors="coerce")
    cmp_df = cmp_df.dropna(subset=["Burden_diff_treat_minus_comp", "VarLog_diff_treat_minus_comp"]).copy()
    if cmp_df.empty:
        _plot_feature_mechanism_unavailable(experiment, out_path, "feature comparison rows have no numeric burden/variance differences")
        return
    # Negative treatment-minus-comparator values mean the treatment group is lower than comparator.
    cmp_df["Burden_reduction"] = -cmp_df["Burden_diff_treat_minus_comp"]
    cmp_df["Entropy_component_reduction"] = -cmp_df["VarLog_diff_treat_minus_comp"]
    cmp_df["Burden_reduction_norm"] = _normalize_within_experiment(cmp_df["Burden_reduction"])
    cmp_df["Entropy_component_reduction_norm"] = _normalize_within_experiment(cmp_df["Entropy_component_reduction"])
    cmp_df["Mechanism_score"] = (cmp_df["Burden_reduction_norm"] + cmp_df["Entropy_component_reduction_norm"]) / 2.0
    top = cmp_df.sort_values("Mechanism_score", ascending=False).head(10).copy()
    if top.empty:
        _plot_feature_mechanism_unavailable(experiment, out_path, "no feature driver rows remain after ranking")
        return
    top = top.sort_values("Mechanism_score", ascending=True).reset_index(drop=True)
    y = np.arange(len(top))
    treatment = str(top["Treatment_group"].dropna().iloc[0]) if "Treatment_group" in top.columns and top["Treatment_group"].notna().any() else str(exp_res.get("treatment_group") or "Treatment")
    comparator = str(top["Comparator_group"].dropna().iloc[0]) if "Comparator_group" in top.columns and top["Comparator_group"].notna().any() else str(exp_res.get("comparator_group") or "Comparator")
    week_text = f"{terminal_week}W" if terminal_week is not None else "selected endpoint"

    fig, axes = plt.subplots(1, 2, figsize=(13.8, max(5.0, 0.48 * len(top) + 2.2)), sharey=True)
    panels = [
        (axes[0], "Burden_reduction", "Lower pathological burden", _BRI_COLOR),
        (axes[1], "Entropy_component_reduction", "Lower entropy-related variance", _ENTROPY_COLOR),
    ]
    for ax, col, title, color in panels:
        vals = pd.to_numeric(top[col], errors="coerce").to_numpy(dtype=float)
        colors = [color if v >= 0 else "#A0AEC0" for v in vals]
        ax.barh(y, vals, color=colors, alpha=0.92, edgecolor="white", linewidth=0.7)
        ax.axvline(0, color="#4A5568", linewidth=0.9)
        for yi, v in zip(y, vals):
            if np.isfinite(v):
                ha = "left" if v >= 0 else "right"
                dx = 0.02 * max(np.nanmax(np.abs(vals)), 1.0)
                ax.text(v + (dx if v >= 0 else -dx), yi, f"{v:+.2f}", va="center", ha=ha, fontsize=8.0, color=_TEXT_COLOR)
        ax.set_title(title, fontsize=11.0, color=_TEXT_COLOR, pad=10)
        _apply_academic_axis(ax, "x")
        ax.set_xlabel(f"{treatment} improvement vs {comparator}")
    axes[0].set_yticks(y)
    axes[0].set_yticklabels(top["Feature"].astype(str).tolist(), fontsize=9)
    axes[1].tick_params(axis="y", left=False, labelleft=False)

    note = (
        f"Endpoint {week_text}. Bars show treatment-vs-comparator reductions after direction-aware preprocessing. "
        "Positive bars suggest the treatment group is lower than comparator for that feature component. "
        "Left panel reflects disease-burden direction; right panel reflects the feature's entropy/variance component."
    )
    _legacy._add_top_title_and_note(
        fig,
        title=f"Experiment {experiment}: feature drivers potentially linking {treatment} to lower entropy",
        note=note,
        legend_handles=None,
        legend_labels=None,
        title_y=0.985,
        note_y=0.945,
        note_width=130,
        note_fontsize=8.7,
    )
    fig.tight_layout(rect=[0.13, 0.06, 0.98, 0.82])
    fig.savefig(out_path, dpi=330, bbox_inches="tight")
    plt.close(fig)
