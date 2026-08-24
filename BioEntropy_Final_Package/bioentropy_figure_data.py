"""Per-figure numerical data export in a GraphPad Prism-ready layout.

Writes the numbers each main figure draws: bar/marker heights, every scatter
point, and the delta / p-value / significance annotations printed on the figure.
Values come from the same helpers the plotting code calls, so a table cannot
drift from its figure.

The bar is not always the mean of the dots, so every ``*_bars.csv`` carries
``Mean_of_plotted_points`` and a ``Bar_equals_mean_of_points`` flag. For raw
disease deviation (Fig 01) the bar is the mean of its dots. For the normalized
figures (Fig 02/04/05/08) the bar is normalize(group mean) while the dots are
normalize(subject value); the 0-100 transform applies a gamma exponent and is
therefore non-linear, so the two differ. For entropy (Fig 03) the bar is a
group-level covariance quantity and the dots are per-subject reference-state
contributions.

Also writes 图表-数据对照表.md/.csv mapping each visual element to the file and
column holding its number, with the formula used to derive it.
"""
from __future__ import annotations

from functools import partial
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

import bioentropy_display_state as _dstate  # noqa: F401  (render-pass state)
import bioentropy_plots as _plots

__all__ = ["export_figure_data_tables", "FIGURE_DATA_DIRNAME"]

FIGURE_DATA_DIRNAME = "图形数据_GraphPad"


# --------------------------------------------------------------------------- #
# small helpers
# --------------------------------------------------------------------------- #
def _base_metric(col: str) -> str:
    return col[:-5] if str(col).endswith("_mean") else str(col)


def _write(df: pd.DataFrame, path: Path) -> Optional[Path]:
    """Write a CSV (utf-8-sig so Excel/Prism read Chinese headers correctly)."""
    if not isinstance(df, pd.DataFrame) or df.empty:
        return None
    df.to_csv(path, index=False, encoding="utf-8-sig")
    return path


def _first_col(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    for c in candidates:
        if c and c in df.columns:
            return c
    return None


def _points_long(points: Dict[tuple, np.ndarray], multi_week: bool) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    for (group, week), arr in points.items():
        for i, v in enumerate(np.asarray(arr, dtype=float), start=1):
            row: Dict[str, Any] = {"Group": group}
            if multi_week:
                row["Week"] = week
            row["Replicate"] = i
            row["Value"] = float(v)
            rows.append(row)
    return pd.DataFrame(rows)


def _points_prism(points: Dict[tuple, np.ndarray], multi_week: bool) -> pd.DataFrame:
    cols: Dict[str, List[float]] = {}
    for (group, week), arr in points.items():
        name = f"{group}_{week}W" if multi_week else str(group)
        cols[name] = [float(x) for x in np.asarray(arr, dtype=float)]
    if not cols:
        return pd.DataFrame()
    h = max(len(v) for v in cols.values())
    return pd.DataFrame({k: v + [np.nan] * (h - len(v)) for k, v in cols.items()})


def _collect_points(sample_level: pd.DataFrame, experiment: str, base: str,
                    groups: List[str], weeks: List[int], multi_week: bool,
                    entropy_contribution: bool = False,
                    normalize: Optional[str] = None) -> Dict[tuple, np.ndarray]:
    """Per-subject plotted values keyed by (group, week), mirroring the overlays."""
    if not isinstance(sample_level, pd.DataFrame) or sample_level.empty:
        return {}
    sub = sample_level
    if "Experiment" in sub.columns:
        sub = sub[sub["Experiment"].astype(str) == str(experiment)]
    if sub.empty or "Group" not in sub.columns:
        return {}
    src = "NMD" if entropy_contribution else base
    if src not in sub.columns:
        return {}
    wk = pd.to_numeric(sub.get("Week"), errors="coerce") if "Week" in sub.columns else None

    cells: Dict[tuple, np.ndarray] = {}
    flat_vals: List[float] = []
    flat_grp: List[str] = []
    for group in groups:
        gmask = sub["Group"].astype(str) == str(group)
        for w in (weeks if multi_week else [weeks[0] if weeks else 0]):
            mask = gmask if (wk is None or not multi_week) else (gmask & (wk == int(w)))
            v = pd.to_numeric(sub.loc[mask, src], errors="coerce").to_numpy(dtype=float)
            v = v[np.isfinite(v)]
            if entropy_contribution:
                v = v[v > 0]
                if v.size:
                    v = 0.5 * np.log(2 * np.pi * np.e * v * v)
            if v.size == 0:
                continue
            cells[(str(group), int(w))] = v
            flat_vals.extend(float(x) for x in v)
            flat_grp.extend([str(group)] * v.size)

    if normalize and cells:
        series, grp = pd.Series(flat_vals, dtype=float), pd.Series(flat_grp)
        if normalize == "linear":
            norm = _plots._normalize_anchor_linear_0_100(series, grp).to_numpy(dtype=float)
        else:
            norm = _plots._normalize_display_window(series, grp).to_numpy(dtype=float)
        i = 0
        for key in list(cells):
            n = cells[key].size
            cells[key] = norm[i:i + n]
            i += n
    return cells


def _emit(out_dir: Path, stem: str, bars: Optional[pd.DataFrame],
          points: Dict[tuple, np.ndarray], multi_week: bool,
          written: List[str], bar_value_col: Optional[str] = None) -> None:
    """Write the bars / points / prism trio, annotating bar-vs-points agreement."""
    if bars is not None and not bars.empty:
        bars = bars.copy()
        if points and bar_value_col and bar_value_col in bars.columns:
            means = {k: float(np.mean(v)) for k, v in points.items() if len(v)}
            if multi_week and "Week" in bars.columns:
                mv = [means.get((str(g), int(w)), np.nan)
                      for g, w in zip(bars["Group"], bars["Week"])]
            else:
                only_week = next(iter({k[1] for k in means}), 0)
                mv = [means.get((str(g), only_week), np.nan) for g in bars["Group"]]
            bars["Mean_of_plotted_points"] = mv
            d = (bars["Mean_of_plotted_points"] - bars[bar_value_col]).abs()
            bars["Bar_equals_mean_of_points"] = np.where(
                d.isna(), "", np.where(d < 1e-6, "YES", "NO"))
        if _write(bars, out_dir / f"{stem}_bars.csv"):
            written.append(f"{stem}_bars.csv")
    if points:
        if _write(_points_long(points, multi_week), out_dir / f"{stem}_points.csv"):
            written.append(f"{stem}_points.csv")
        if _write(_points_prism(points, multi_week), out_dir / f"{stem}_prism.csv"):
            written.append(f"{stem}_prism.csv")


def _annotation_frame(comparisons: Any, weeks: List[int], diff_cols: List[str],
                      p_cols: List[str]) -> pd.DataFrame:
    """Δ / p-value / significance label printed on a figure (from comparisons)."""
    if not isinstance(comparisons, pd.DataFrame) or comparisons.empty:
        return pd.DataFrame()
    c = comparisons.copy()
    if "Week" not in c.columns:
        return pd.DataFrame()
    c["Week"] = pd.to_numeric(c["Week"], errors="coerce")
    dcol = _first_col(c, diff_cols)
    pcol = _first_col(c, p_cols)
    if not dcol or not pcol:
        return pd.DataFrame()
    rows: List[Dict[str, Any]] = []
    for w in weeks:
        sel = c[c["Week"] == w]
        if sel.empty:
            continue
        r = sel.iloc[0]
        d, p = r.get(dcol), r.get(pcol)
        if pd.isna(d) or pd.isna(p):
            continue
        rows.append({
            "Week": int(w),
            "Delta_treatment_minus_comparator": float(d),
            "p_value": float(p),
            "Significance_label": _plots._signif_label(p),
            "Delta_source_column": dcol,
            "p_source_column": pcol,
            "Source_table": "03_治疗比较.xlsx",
        })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# main entry point
# --------------------------------------------------------------------------- #

def _bar_frame(sub: pd.DataFrame, multi_week: bool, groups: List[str],
               cols: Dict[str, Any], order: Optional[List[str]] = None) -> pd.DataFrame:
    """One bar/marker table for a figure, in the figure's own group order."""
    bc: Dict[str, Any] = {"Group": sub["Group"]}
    if multi_week:
        bc["Week"] = sub["Week"]
    if "n" in sub.columns:
        bc["n"] = sub["n"]
    bc.update(cols)
    df = pd.DataFrame(bc)
    seq = order if order is not None else groups
    df["_o"] = df["Group"].map({g: i for i, g in enumerate(seq)}).fillna(999)
    df = df.sort_values(["Week", "_o"] if multi_week else ["_o"]).drop(columns=["_o"])
    return df.reset_index(drop=True)


def export_figure_data_tables(results: Dict[str, Any], output_dir: Path) -> List[Path]:
    """Write the plotted numbers for every value-bearing main figure."""
    out_root = Path(output_dir) / FIGURE_DATA_DIRNAME
    out_root.mkdir(parents=True, exist_ok=True)
    files: List[Path] = []
    manifest: List[Dict[str, Any]] = []
    mapping_rows: List[Dict[str, Any]] = []

    def _map(exp: str, fig: str, element: str, file: str, column: str, how: str) -> None:
        mapping_rows.append({"实验": exp, "图": fig, "图中元素": element,
                             "文件": file, "列名": column, "计算方式": how})

    for exp, exp_res in (results.get("experiments") or {}).items():
        summary = exp_res.get("summary", pd.DataFrame())
        if not isinstance(summary, pd.DataFrame) or summary.empty:
            continue
        sample_level = exp_res.get("sample_level", pd.DataFrame())
        comparisons = exp_res.get("comparisons", exp_res.get("comparison", pd.DataFrame()))

        burden_kind = _plots._preferred_burden_kind(exp_res)
        meta = _plots._burden_metadata(burden_kind)
        bcol, short = meta["mean"], str(meta["short"])
        if bcol not in summary.columns or "Entropy" not in summary.columns:
            continue
        base = _base_metric(bcol)

        sub = summary.copy()
        if "Experiment" in sub.columns:
            sub = sub[sub["Experiment"].astype(str) == str(exp)].copy()
        if sub.empty:
            continue
        sub["Week"] = (pd.to_numeric(sub["Week"], errors="coerce").fillna(0).astype(int)
                       if "Week" in sub.columns else 0)
        sub["Group"] = sub["Group"].astype(str)
        weeks = sorted(sub["Week"].unique().tolist())
        multi_week = len(weeks) > 1
        groups = _plots._ordered_groups_for_plot(
            sub, exp_res.get("treatment_group"), exp_res.get("comparator_group"))

        exp_dir = out_root / str(exp)
        exp_dir.mkdir(parents=True, exist_ok=True)
        written: List[str] = []

        disp_b_gamma = _plots._normalize_display_window(sub[bcol], sub["Group"])
        disp_b_lin = _plots._normalize_anchor_linear_0_100(sub[bcol], sub["Group"])
        disp_e_gamma = _plots._normalize_display_window(sub["Entropy"], sub["Group"])
        disp_e_lin = _plots._normalize_anchor_linear_0_100(sub["Entropy"], sub["Group"])

        bar_frame = partial(_bar_frame, sub, multi_week, groups)

        # ================= Figure 01 : raw disease deviation ================ #
        burden_pts = _collect_points(sample_level, exp, base, groups, weeks, multi_week)
        f1: Dict[str, Any] = {f"{short}_mean_plotted": sub[bcol]}
        lo = _first_col(sub, [meta.get("low"), f"{base}_mean_ci_low", f"{base}_ci_low"])
        hi = _first_col(sub, [meta.get("high"), f"{base}_mean_ci_high", f"{base}_ci_high"])
        sd = _first_col(sub, [f"{base}_sd"])
        se = _first_col(sub, [f"{base}_boot_se"])
        if lo:
            f1[f"{short}_CI_low"] = sub[lo]
        if hi:
            f1[f"{short}_CI_high"] = sub[hi]
        if sd:
            f1[f"{short}_SD"] = sub[sd]
        if se:
            f1[f"{short}_bootstrap_SE"] = sub[se]
        s1 = f"Fig01_{short}_raw"
        _emit(exp_dir, s1, bar_frame(f1), burden_pts, multi_week, written,
              bar_value_col=f"{short}_mean_plotted")
        _map(exp, "图01", f"柱高/折线点（{short} 组均值）", f"{s1}_bars.csv",
             f"{short}_mean_plotted", f"该组所有样本 {short} 的算术平均")
        _map(exp, "图01", "散点（每个受试对象）", f"{s1}_points.csv / _prism.csv", "Value",
             f"样本级 {short}（见 09_样本级结果.xlsx 的 {base} 列）")
        if lo:
            _map(exp, "图01", "误差棒（若显示）", f"{s1}_bars.csv",
                 f"{short}_CI_low / {short}_CI_high", "自助法 200 次重抽样、seed=42 的 95% CI")
        ann1 = _annotation_frame(comparisons, weeks,
                                 [f"{base}_diff_treat_minus_comp", f"{base}_diff"],
                                 [f"{base}_p_mannwhitney", f"{base}_p_boot"])
        if _write(ann1, exp_dir / f"{s1}_annotations.csv"):
            written.append(f"{s1}_annotations.csv")
            _map(exp, "图01", "显著性标注（ns/*/**/***、Δ、p）", f"{s1}_annotations.csv",
                 "Significance_label / Delta_treatment_minus_comparator / p_value",
                 "来自 03_治疗比较.xlsx：Δ=治疗组−对照组均值差；p 为 Mann–Whitney U 检验")

        # ================= Figure 02 : normalized deviation ================= #
        burden_pts_disp = _collect_points(sample_level, exp, base, groups, weeks,
                                          multi_week, normalize="display")
        s2 = f"Fig02_{short}_normalized"
        _emit(exp_dir, s2,
              bar_frame({f"{short}_raw": sub[bcol],
                          "Linear_anchored_0_100": disp_b_lin,
                          "Display_plotted_0_100_gamma3": disp_b_gamma}),
              burden_pts_disp, multi_week, written,
              bar_value_col="Display_plotted_0_100_gamma3")
        _map(exp, "图02", "柱高/折线点（归一化显示值）", f"{s2}_bars.csv",
             "Display_plotted_0_100_gamma3",
             "先线性锚定 Normal=0、model=100，再施加 gamma=3 显示增强")
        _map(exp, "图02", "散点", f"{s2}_points.csv / _prism.csv", "Value",
             "对每个样本值施加与柱相同的归一化变换（注意：非线性，故散点均值≠柱高）")

        # ================= Figure 03 : raw entropy ========================== #
        ent_pts = _collect_points(sample_level, exp, base, groups, weeks, multi_week,
                                  entropy_contribution=True)
        f3: Dict[str, Any] = {"Entropy_plotted": sub["Entropy"]}
        for c, name in (("Entropy_ci_low", "Entropy_CI_low"),
                        ("Entropy_ci_high", "Entropy_CI_high"),
                        ("Entropy_boot_se", "Entropy_bootstrap_SE"),
                        ("Entropy_method", "Entropy_method")):
            if c in sub.columns:
                f3[name] = sub[c]
        s3 = "Fig03_Entropy_raw"
        _emit(exp_dir, s3, bar_frame(f3), ent_pts, multi_week, written,
              bar_value_col="Entropy_plotted")
        _map(exp, "图03", "柱高/折线点（系统熵）", f"{s3}_bars.csv", "Entropy_plotted",
             "组内多指标分布的对数行列式熵 H=½·log[(2πe)^p·|Σ|]，Σ 为 Ledoit–Wolf 收缩估计")
        _map(exp, "图03", "误差棒", f"{s3}_bars.csv", "Entropy_CI_low / Entropy_CI_high",
             "自助法 200 次重抽样、seed=42 的 95% CI")
        _map(exp, "图03", "散点", f"{s3}_points.csv / _prism.csv", "Value",
             "每个样本的参考态熵贡献 0.5·ln(2πe·NMD²)（与柱是不同的量）")
        ann3 = _annotation_frame(comparisons, weeks,
                                 ["Entropy_diff_treat_minus_comp", "Entropy_diff"],
                                 ["Entropy_p_boot", "Entropy_p_mannwhitney"])
        if _write(ann3, exp_dir / f"{s3}_annotations.csv"):
            written.append(f"{s3}_annotations.csv")
            _map(exp, "图03", "显著性标注（ns/*/**/***、Δ、p）", f"{s3}_annotations.csv",
                 "Significance_label / Delta_treatment_minus_comparator / p_value",
                 "来自 03_治疗比较.xlsx：Δ=治疗组−对照组熵差；p 为自助法/Mann–Whitney")

        # ================= Figure 04 : normalized entropy =================== #
        ent_pts_disp = _collect_points(sample_level, exp, base, groups, weeks,
                                       multi_week, entropy_contribution=True,
                                       normalize="display")
        s4 = "Fig04_Entropy_normalized"
        _emit(exp_dir, s4,
              bar_frame({"Entropy_raw": sub["Entropy"],
                          "Linear_anchored_0_100": disp_e_lin,
                          "Display_plotted_0_100_gamma3": disp_e_gamma}),
              ent_pts_disp, multi_week, written,
              bar_value_col="Display_plotted_0_100_gamma3")
        _map(exp, "图04", "柱高/折线点（熵归一化显示值）", f"{s4}_bars.csv",
             "Display_plotted_0_100_gamma3", "同图02，对熵做线性锚定 + gamma=3")
        _map(exp, "图04", "散点", f"{s4}_points.csv / _prism.csv", "Value",
             "对每个样本的熵贡献施加同一归一化变换")

        # ================= Figure 05 : grouped bar ========================== #
        s5 = f"Fig05_{short}_Entropy_grouped_bar"
        _emit(exp_dir, s5,
              bar_frame({f"{short}_raw_label": sub[bcol],
                          f"{short}_bar_height_0_100": disp_b_gamma,
                          "Entropy_raw_label": sub["Entropy"],
                          "Entropy_bar_height_0_100": disp_e_gamma}),
              {}, multi_week, written)
        _emit(exp_dir, f"{s5}_{short}Dots", None, burden_pts_disp, multi_week, written)
        _emit(exp_dir, f"{s5}_EntropyDots", None, ent_pts_disp, multi_week, written)
        _map(exp, "图05", f"左柱高（{short}）", f"{s5}_bars.csv",
             f"{short}_bar_height_0_100", "归一化显示值（与图02柱高相同）")
        _map(exp, "图05", "右柱高（Entropy）", f"{s5}_bars.csv",
             "Entropy_bar_height_0_100", "归一化显示值（与图04柱高相同）")
        _map(exp, "图05", "柱上标注的数字", f"{s5}_bars.csv",
             f"{short}_raw_label / Entropy_raw_label", "对应的原始值（非柱高）")
        _map(exp, "图05", f"左柱散点（{short}）", f"{s5}_{short}Dots_prism.csv", "各组一列",
             "与图02散点相同")
        _map(exp, "图05", "右柱散点（Entropy）", f"{s5}_EntropyDots_prism.csv", "各组一列",
             "与图04散点相同")

        if not multi_week:
            # ============ Figures 07 / 08 : state plots (single week) ======= #
            order = sub.groupby("Group")[bcol].mean().sort_values().index.tolist()
            rank = {g: i + 1 for i, g in enumerate(order)}
            s7 = f"Fig07_{short}_Entropy_state_raw"
            _emit(exp_dir, s7,
                  bar_frame({"Plot_order_x": sub["Group"].map(rank),
                              f"y_{short}_raw": sub[bcol],
                              "Entropy_reference_only": sub["Entropy"]}, order=order),
                  burden_pts, multi_week, written, bar_value_col=f"y_{short}_raw")
            s8 = f"Fig08_{short}_Entropy_state_normalized"
            _emit(exp_dir, s8,
                  bar_frame({"Plot_order_x": sub["Group"].map(rank),
                              f"y_{short}_linear_0_100": disp_b_lin,
                              "Entropy_reference_only": sub["Entropy"]}, order=order),
                  _collect_points(sample_level, exp, base, groups, weeks, multi_week,
                                  normalize="linear"),
                  multi_week, written, bar_value_col=f"y_{short}_linear_0_100")
            for s, ycol, note in ((s7, f"y_{short}_raw", "原始值"),
                                  (s8, f"y_{short}_linear_0_100", "线性锚定 0–100（无 gamma）")):
                _map(exp, "图07" if s is s7 else "图08", "x 轴（组的先后顺序）",
                     f"{s}_bars.csv", "Plot_order_x",
                     f"按各组原始 {short} 均值升序排列（healthy→diseased），x 不是熵")
                _map(exp, "图07" if s is s7 else "图08", "大圆点 y 值", f"{s}_bars.csv",
                     ycol, note)
                _map(exp, "图07" if s is s7 else "图08", "小灰点（每个受试对象）",
                     f"{s}_points.csv / _prism.csv", "Value", f"对应尺度的样本值（{note}）")
                _map(exp, "图07" if s is s7 else "图08", "连线", "（由大圆点生成）", "-",
                     "通过各组大圆点的 PCHIP 单调样条，无额外数据")
        else:
            # ============ Figure 08 : weekly state panels =================== #
            s8 = f"Fig08_{short}_Entropy_weekly_state_panels"
            _emit(exp_dir, s8,
                  bar_frame({"x_Entropy": sub["Entropy"], f"y_{short}": sub[bcol]}),
                  burden_pts, multi_week, written, bar_value_col=f"y_{short}")
            _map(exp, "图08", "每周每组的气泡（x=熵, y=负担）", f"{s8}_bars.csv",
                 f"x_Entropy / y_{short}", "该周该组的熵与负担均值")
            _map(exp, "图08", "小散点", f"{s8}_points.csv / _prism.csv", "Value",
                 "该周该组每个受试对象的负担值")
            ann8 = _annotation_frame(comparisons, weeks,
                                     [f"{base}_diff_treat_minus_comp", f"{base}_diff"],
                                     [f"{base}_p_mannwhitney", f"{base}_p_boot"])
            if _write(ann8, exp_dir / f"{s8}_annotations.csv"):
                written.append(f"{s8}_annotations.csv")
                _map(exp, "图08", "每周面板上方的显著性星号", f"{s8}_annotations.csv",
                     "Significance_label / p_value", "来自 03_治疗比较.xlsx 的每周 p 值")

        # ================= mechanism drivers ================================ #
        mech = exp_res.get("feature_mechanism", pd.DataFrame())
        if isinstance(mech, pd.DataFrame) and not mech.empty:
            m = mech.copy()
            if "Experiment" in m.columns:
                m = m[m["Experiment"].astype(str) == str(exp)]
            fig_no = "图07" if multi_week else "图09"
            stem = ("Fig07_feature_mechanism_drivers" if multi_week
                    else "Fig09_feature_mechanism_drivers")
            if _write(m, exp_dir / f"{stem}_bars.csv"):
                written.append(f"{stem}_bars.csv")
                _map(exp, fig_no, "各指标的机制分解条", f"{stem}_bars.csv", "（见表内各列）",
                     "按指标给出各组相对参考的偏离/改善分解")

        files.extend(exp_dir / w for w in written)
        manifest.append({"experiment": str(exp), "burden": short,
                         "multi_week": multi_week, "files": written})

    if manifest:
        readme = out_root / "README_图形数据说明.md"
        readme.write_text(_readme_text(manifest), encoding="utf-8")
        files.append(readme)
    if mapping_rows:
        mdf = pd.DataFrame(mapping_rows)
        p_csv = out_root / "图表-数据对照表.csv"
        mdf.to_csv(p_csv, index=False, encoding="utf-8-sig")
        files.append(p_csv)
        p_md = out_root / "图表-数据对照表.md"
        p_md.write_text(_mapping_md(mdf), encoding="utf-8")
        files.append(p_md)
    return files


def _mapping_md(mdf: pd.DataFrame) -> str:
    lines: List[str] = []
    a = lines.append
    a("# 图表 → 数据 对照表")
    a("")
    a("本表把**图中出现的每一个视觉元素**（柱高、折线点、散点、误差棒、显著性标注、")
    a("坐标轴取值）对应到导出文件中的具体列，并说明该数值如何由原始数据算出。")
    a("")
    a("> 说明：所有导出数值均由绘图所用的同一套函数生成，因此与图形完全一致。")
    a("> 样本级原始指标见结果表 `09_样本级结果.xlsx`；组间统计见 `03_治疗比较.xlsx`。")
    a("")
    for exp in mdf["实验"].unique():
        sub = mdf[mdf["实验"] == exp]
        a(f"## 实验：{exp}")
        a("")
        for fig in sub["图"].unique():
            s = sub[sub["图"] == fig]
            a(f"### {fig}")
            a("")
            a("| 图中元素 | 文件 | 列名 | 计算方式 |")
            a("|---|---|---|---|")
            for _, r in s.iterrows():
                a(f"| {r['图中元素']} | `{r['文件']}` | `{r['列名']}` | {r['计算方式']} |")
            a("")
    return "\n".join(lines)


def _readme_text(manifest: List[Dict[str, Any]]) -> str:
    lines: List[str] = []
    a = lines.append
    a("# 图形数据（GraphPad Prism 可直接使用）")
    a("")
    a("本文件夹给出主图中**实际绘制的每一个数值**。逐元素的对照见 `图表-数据对照表.md`。")
    a("")
    a("## ⚠ 最重要的一点：柱高 ≠ 散点均值（除图 01 外）")
    a("")
    a("每个 `*_bars.csv` 都带有 `Mean_of_plotted_points` 与 `Bar_equals_mean_of_points` 两列：")
    a("")
    a("| 图 | 柱高 = 散点均值？ | 原因 |")
    a("|---|---|---|")
    a("| 图 01（原始 NMD/NRBS） | **YES** | 柱＝各组均值，点＝各样本值 |")
    a("| 图 02 / 04 / 05 | **NO** | 柱＝normalize(组均值)，点＝normalize(样本值)；gamma=3 非线性 |")
    a("| 图 03（熵） | **NO** | 柱＝组水平协方差熵；点＝样本参考态熵贡献，二者不同量 |")
    a("| 图 08（状态图） | **NO** | 同上，y 为线性锚定 0–100 |")
    a("")
    a("**因此：把散点粘进 Prism 后，不要让 Prism 自动算均值当柱高**（图 01 除外）；")
    a("请把 `*_bars.csv` 中的柱高作为**独立的一组数据**输入。")
    a("")
    a("## 文件类型")
    a("")
    a("| 后缀 | 内容 |")
    a("|---|---|")
    a("| `*_bars.csv` | 柱/折线/标记的高度，含 CI、SD、n |")
    a("| `*_points.csv` | 每个散点（长表：Group / Week / Replicate / Value） |")
    a("| `*_prism.csv` | 同样的散点，每组一列，可直接粘贴进 Prism |")
    a("| `*_annotations.csv` | 图上显示的 Δ、p 值与显著性标记（ns/*/**/***） |")
    a("")
    a("## 在 Prism 中重建")
    a("")
    a("- **图 01**：Column 表粘贴 `_prism.csv`，选 *Scatter dot plot with mean*；")
    a("  误差棒用 `_bars.csv` 的 CI 列（自助法 200 次、seed=42），不要用 Prism 自算的 SD/SEM。")
    a("- **图 02/03/04**：柱高单独输入（`_bars.csv`），散点作为第二数据集叠加。")
    a("- **图 05**：两组柱高 + 两套散点（`*Dots_prism.csv`）。")
    a("- **图 07/08（单时间点）**：XY 图，x=`Plot_order_x`，y=`y_*`，平滑连线；散点另加。")
    a("")
    a("## 本次导出的实验")
    a("")
    for m in manifest:
        a(f"- **{m['experiment']}**：负担指标 {m['burden']}，"
          f"{'多时间点' if m['multi_week'] else '单时间点'}，共 {len(m['files'])} 个文件")
    a("")
    return "\n".join(lines)
