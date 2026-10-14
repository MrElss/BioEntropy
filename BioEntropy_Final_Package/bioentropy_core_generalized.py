
from __future__ import annotations

import logging
import zipfile
from pathlib import Path
from typing import Any, Dict, List, Optional

import openpyxl
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from openpyxl.utils import get_column_letter

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
)


import bioentropy_display_state as _dstate
from bioentropy_plots import *  # noqa: F401,F403  (plotting + group-order helpers)


# Color palette and DISPLAY_CONTRAST_GAMMA now live in bioentropy_constants
# (imported above). Mutable display state remains module-local here.


def build_markdown_report(results: Dict[str, Any]) -> str:
    base = _legacy.build_markdown_report(results)
    addon = [
        "",
        "## NMD × Entropy 联合展示建议",
        "- 对临床纵向实验，多时间点图继续以 Week 为横坐标。",
        "- 对动物实验或其他单时间点实验，软件会自动改用 **给药组/处理组** 作为横坐标，不再把 0W 当作真实时间维度解释。",
        "- 多时间点实验通常输出 8 个主图/文件，单时间点实验通常输出 9 个；原始值、归一化展示、联合展示、状态图和机制图均保留。",
        "- 多时间点 Fig02/04 先对齐各组基线，再按统一锚点进行0–100缩放与gamma=3展示，起点相同且非零。仅有下降或不变轨迹时起点为100；缺少0W时标明最早有效周。展示指数不是原始值或百分比。",
        "- 疾病偏离度个体散点和组曲线使用同一份组汇总数据确定缩放锚点；截断及非线性变换后，组曲线不必等于散点的平均值。熵仅展示组估计，不叠加个体对数距离。",
        f"- 归一化展示同时覆盖 NMD/BRI 和 Entropy：结果表保留线性锚定值（Normal=0、model=100）与 gamma={DISPLAY_CONTRAST_GAMMA:g} 治疗窗口增强展示值；正式判断仍回到原始值、线性锚定值、置信区间和统计检验。",
        "- 桌面端会自动检测每个实验的组别，并在界面中生成可拖拽排序列表；拖拽后的顺序只影响图和展示表的排列，不改变计算结果。",
    ]
    return base + "\n".join(addon)


def build_input_guide_html(out_path: Path) -> Path:
    """Write the current input guide, overriding legacy manual-order wording."""
    _legacy.build_input_guide_html(out_path)
    path = Path(out_path)
    text = path.read_text(encoding="utf-8")
    old_norm = (
        "<li>归一化疾病偏离度仅用于展示：有 Normal/健康组和 model/模型组时，0–100 尺度采用生物学锚点归一化，"
        "即 normal control = 0、model control = 100；真实判断仍看原始数值、置信区间和统计检验。</li>"
    )
    new_norm = (
        f"<li>归一化展示同时覆盖疾病偏离度和 Entropy：有 Normal/健康组和 model/模型组时，先生成可追溯的线性锚定值"
        f"（normal control = 0、model control = 100），主展示图再使用 gamma={DISPLAY_CONTRAST_GAMMA:g} 的治疗窗口增强值拉开高疾病负担区间；"
        "真实判断仍看原始数值、线性锚定值、置信区间和统计检验。</li>"
    )
    old_order = (
        "<li>组别显示顺序可以自定义。桌面端“自定义组排序”可填写全局顺序，也可按实验号分别填写，例如："
        "<code>exp1: normal control, treatment, comparator, model control; exp2: normal control, low-dose, mid-dose, high-dose, comparator, model control</code>。</li>"
    )
    new_order = (
        "<li>组别显示顺序无需手动输入：桌面端会自动检测每个实验有哪些组，并在文件识别区域生成“图中组别顺序”列表；"
        "按住组名上下拖动即可调整主图、辅助图和归一化展示表的显示顺序。</li>"
    )
    old_figures = (
        "<li><code>figures</code> 文件夹只保留 7 个主图/主文件：原始马氏距离/疾病偏离度、归一化疾病偏离度展示图、"
        "Entropy、二维柱形图、可旋转 3D HTML、二维状态图、药物机制驱动指标图。</li>"
    )
    new_figures = (
        "<li><code>figures</code> 文件夹只保留 9 个主图/主文件，顺序为：马氏距离/疾病偏离度真实值、马氏距离/疾病偏离度归一化值、"
        "熵真实值、熵归一化值、二维柱形图、可旋转 3D HTML、真实值可点击走势点图、归一化可点击走势点图、药物机制驱动指标图。</li>"
    )
    old_aux = "<li><code>figures_auxiliary</code> 文件夹会额外输出辅助图：相图、治疗反应排序图、交互式 response landscape。</li>"
    old_aux2 = "<li><code>figures_auxiliary</code> 文件夹会额外输出更多辅助图：相图、治疗反应排序图、交互式 response landscape、静态状态图、双热图、归一化状态图等。</li>"
    # figures_auxiliary is no longer produced; drop any stale aux-folder note.
    new_aux = ""
    text = (
        text.replace(old_figures, new_figures)
        .replace(old_aux, new_aux)
        .replace(old_aux2, new_aux)
        .replace(old_norm, new_norm)
        .replace(old_order, new_order)
    )
    text = text.replace(
        "后面为各指标列。</p>",
        "后面为各指标列；每行一个样本，数值型 ID 在同一文件中须全局唯一，包括不同组之间。</p>",
    ).replace(
        "没有任何 Normal 参考时：主要输出 Entropy 和 BRI（相对锚定周变化）。",
        "没有真实 Normal 样本时：有可匹配的上传/内置参考范围则使用 NRBS；否则使用相对最早时间点的 BRI。",
    ).replace(
        "如果没有真实 Normal 参考，第一张主图会使用 BRI/HDI 作为基线参考疾病偏离度；只有存在 Normal 参考时才严格解释为 True-Normal NMD。",
        "第一张主图的负担指标优先级为 NMD > NRBS > BRI；只有存在真实 Normal 参考时才解释为 True-Normal NMD。",
    ).replace(
        new_figures,
        "<li><code>figures</code> 根据实验类型保留主图：多时间点通常为 8 个，单时间点通常为 9 个；包含原值、归一化展示、联合展示、3D HTML、状态图和机制驱动图。"
        "多时间点 Fig02/04 展示基线对齐后的归一化疾病偏离度/熵，起点相同且非零；仅有下降或不变轨迹时起点为100。缺失0W时注明最早有效周；统一锚点缩放并使用gamma=3，非原始值或下降百分比。</li>"
        "<li>疾病偏离度散点与组值共用组汇总数据的归一化锚点，展示值截断到 0–100。截断与 gamma 变换后，曲线不必等于散点均值，0/100 也不代表生物学极限。熵图仅显示组估计。</li>",
    )
    text = text.replace("</body></html>", (
        "<h2>8. 运行记录与质量提示</h2>"
        "<p>结果同级的 <code>*_audit</code> 目录保存运行记录、质量提示和指标字典；不修改已有结果，也不写入原 ZIP。</p>"
        "<p>reference-state entropy 使用 <code>mean(NMD²)</code>，与 NMD 存在结构性耦合；协方差熵保留在 "
        "<code>Entropy_robust_covariance</code>。请核查全缺失样本、参考单位和成人范围适用性。</p>"
        "</body></html>"
    ))
    path.write_text(text, encoding="utf-8")
    return path


def _write_dataframe_sheet(path: Path, sheet_name: str, df: pd.DataFrame) -> None:
    wb = openpyxl.load_workbook(path)
    if sheet_name in wb.sheetnames:
        wb.remove(wb[sheet_name])
    ws = wb.create_sheet(sheet_name)
    if df is None or df.empty:
        ws["A1"] = "No data"
        wb.save(path)
        return
    rows = [list(df.columns)] + df.astype(object).where(pd.notna(df), "").values.tolist()
    for row in rows:
        ws.append(row)
    for cell in ws[1]:
        cell.font = openpyxl.styles.Font(bold=True)
        cell.fill = openpyxl.styles.PatternFill("solid", fgColor="D9EAF7")
    for col_idx in range(1, ws.max_column + 1):
        max_len = 0
        for row_idx in range(1, ws.max_row + 1):
            v = ws.cell(row_idx, col_idx).value
            max_len = max(max_len, len(str(v)) if v is not None else 0)
        ws.column_dimensions[get_column_letter(col_idx)].width = min(max(12, max_len + 2), 34)
    wb.save(path)


def _refresh_zip(output_dir: Path, zip_path: Path) -> None:
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for path in output_dir.rglob("*"):
            if path.is_file():
                zf.write(path, arcname=str(path.relative_to(output_dir.parent)))


def _delete_unselected_figure_files(figures_dir: Path, keep_paths: List[Path]) -> None:
    keep = {p.resolve() for p in keep_paths if p.exists()}
    for path in figures_dir.glob("*"):
        if path.is_file() and path.resolve() not in keep:
            try:
                path.unlink()
            except Exception:
                pass


def _append_joint_visual_artifacts(results: Dict[str, Any], bundle: Dict[str, Any]) -> None:
    output_dir = Path(bundle["output_dir"])
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    figure_paths: Dict[str, str] = {}
    keep_main_files: List[Path] = []

    # Transparent display-normalization table, written as its own workbook in the
    # curated 数据表 folder (visual readability only; inference stays on raw metrics).
    norm_df = build_display_normalized_metric_table(results)
    if isinstance(norm_df, pd.DataFrame) and not norm_df.empty:
        tables_dir = output_dir / "数据表"
        idx = _legacy.NECESSARY_TABLE_KEYS.index("Display_Normalized_Index") + 1
        tables_dir.mkdir(parents=True, exist_ok=True)
        _legacy.save_table_bundle_excel(
            {"Display_Normalized_Index": norm_df},
            tables_dir / f"{idx:02d}_Display_Normalized_Index.xlsx",
        )

    for exp, exp_res in results.get("experiments", {}).items():
        previous_display_experiment = _dstate.CURRENT_DISPLAY_EXPERIMENT
        previous_sample_level = getattr(_dstate, "CURRENT_SAMPLE_LEVEL", None)
        previous_summary = _dstate.CURRENT_SUMMARY
        _dstate.CURRENT_DISPLAY_EXPERIMENT = str(exp)
        _dstate.CURRENT_SAMPLE_LEVEL = exp_res.get("sample_level")
        summary_df = exp_res.get("summary", pd.DataFrame())
        _dstate.CURRENT_SUMMARY = summary_df
        if not isinstance(summary_df, pd.DataFrame) or summary_df.empty:
            _dstate.CURRENT_DISPLAY_EXPERIMENT = previous_display_experiment
            _dstate.CURRENT_SAMPLE_LEVEL = previous_sample_level
            _dstate.CURRENT_SUMMARY = previous_summary
            continue

        try:
            treatment_group = exp_res.get("treatment_group")
            comparator_group = exp_res.get("comparator_group")
            comparison_df = exp_res.get("comparisons", exp_res.get("comparison", pd.DataFrame()))
            burden_kind = _preferred_burden_kind(exp_res)
            meta = _burden_metadata(burden_kind)
            # Multi-week experiments use aligned displays in figures 02/04; a single
            # uploaded time point does not imply a longitudinal baseline.
            try:
                multi_week = int(pd.to_numeric(summary_df.get("Week"), errors="coerce").nunique()) > 1
            except Exception:
                multi_week = False

            # Main Figure 01: raw Mahalanobis / disease-deviation values.
            pburden = figures_dir / f"{exp}_main_01_{meta['prefix']}_raw_trajectory.png"
            _copy_or_regenerate_main_metric_plot(
                summary_df, exp, pburden, burden_kind, comparison_df, treatment_group, comparator_group
            )
            if pburden.exists():
                figure_paths[f"{exp}_main_01_{meta['prefix']}_raw_trajectory"] = str(pburden)
                keep_main_files.append(pburden)

            # Main Figure 02: aligned longitudinal display or single-visit levels.
            suffix = "baseline_aligned_normalized" if multi_week else "normalized_index"
            pnorm = figures_dir / f"{exp}_main_02_{meta['prefix']}_{suffix}.png"
            if multi_week:
                plot_burden_baseline_aligned(summary_df, exp, pnorm, burden_kind=burden_kind,
                                             treatment_group=treatment_group, comparator_group=comparator_group)
            else:
                plot_normalized_burden_index(summary_df, exp, pnorm, burden_kind=burden_kind,
                                             treatment_group=treatment_group, comparator_group=comparator_group)
            if pnorm.exists():
                figure_paths[pnorm.stem] = str(pnorm)
                keep_main_files.append(pnorm)

            # Main Figure 03: raw Entropy.
            pent = figures_dir / f"{exp}_main_03_entropy_trajectory.png"
            plot_entropy_trajectory(summary_df, exp, pent, comparison_df, treatment_group, comparator_group)
            if pent.exists():
                figure_paths[f"{exp}_main_03_entropy_trajectory"] = str(pent)
                keep_main_files.append(pent)

            # Main Figure 04: aligned longitudinal entropy or single-visit levels.
            pent_norm = figures_dir / f"{exp}_main_04_entropy_{suffix}.png"
            if multi_week:
                plot_entropy_baseline_aligned(summary_df, exp, pent_norm,
                                              treatment_group=treatment_group, comparator_group=comparator_group)
            else:
                plot_normalized_entropy_index(summary_df, exp, pent_norm,
                                              treatment_group=treatment_group, comparator_group=comparator_group)
            if pent_norm.exists():
                figure_paths[pent_norm.stem] = str(pent_norm)
                keep_main_files.append(pent_norm)

            # Main Figure 05: 2D grouped bar for raw-labeled burden x Entropy display.
            pbar = figures_dir / f"{exp}_main_05_{meta['prefix']}_entropy_grouped_bar.png"
            plot_burden_entropy_grouped_bar(summary_df, exp, pbar, burden_kind=burden_kind,
                                            treatment_group=treatment_group, comparator_group=comparator_group)
            if pbar.exists():
                figure_paths[f"{exp}_main_05_{meta['prefix']}_entropy_grouped_bar"] = str(pbar)
                keep_main_files.append(pbar)

            # Main Figure 06: interactive 3D grouped bar in perspective HTML.
            p3d = figures_dir / f"{exp}_main_06_{meta['prefix']}_entropy_3d_interactive.html"
            plot_burden_entropy_3d_perspective_html(summary_df, exp, p3d, burden_kind=burden_kind,
                                                    treatment_group=treatment_group, comparator_group=comparator_group)
            if p3d.exists():
                figure_paths[f"{exp}_main_06_{meta['prefix']}_entropy_3d_interactive"] = str(p3d)
                keep_main_files.append(p3d)

            # Main Figures 07/08: burden-vs-entropy state plots — only for
            # single-timepoint (animal) experiments. Each group is a marker along
            # the entropy ordering (healthy → diseased) joined by a smooth curve;
            # 07 shows raw values, 08 the normalized display index.
            if not multi_week:
                # Main Figure 07: raw-value state plot.
                pstate_raw = figures_dir / f"{exp}_main_07_{meta['prefix']}_entropy_state.png"
                plot_single_timepoint_nmd_entropy_state(summary_df, exp, pstate_raw, burden_kind=burden_kind,
                                                        treatment_group=treatment_group, comparator_group=comparator_group,
                                                        normalized=False)
                if pstate_raw.exists():
                    figure_paths[f"{exp}_main_07_{meta['prefix']}_entropy_state"] = str(pstate_raw)
                    keep_main_files.append(pstate_raw)

                # Main Figure 08: normalized-display state plot.
                pstate_norm = figures_dir / f"{exp}_main_08_{meta['prefix']}_entropy_state_normalized.png"
                plot_single_timepoint_nmd_entropy_state(summary_df, exp, pstate_norm, burden_kind=burden_kind,
                                                        treatment_group=treatment_group, comparator_group=comparator_group,
                                                        normalized=True)
                if pstate_norm.exists():
                    figure_paths[f"{exp}_main_08_{meta['prefix']}_entropy_state_normalized"] = str(pstate_norm)
                    keep_main_files.append(pstate_norm)

            # Feature mechanism drivers. Numbered so each experiment type stays
            # contiguous: single-timepoint uses 07/08 for the clickable maps, so
            # mechanism is 09 there; multi-week has no 07/08, so mechanism is 07.
            mech_no = "07" if multi_week else "09"
            pmechanism = figures_dir / f"{exp}_main_{mech_no}_feature_mechanism_drivers.png"
            plot_feature_mechanism_drivers(exp_res, exp, pmechanism, burden_kind=burden_kind)
            if pmechanism.exists():
                figure_paths[f"{exp}_main_{mech_no}_feature_mechanism_drivers"] = str(pmechanism)
                keep_main_files.append(pmechanism)

            # Main Figure 08 (multi-week only): weekly state small-multiples. One
            # panel per week; groups placed by (entropy, burden) on shared axes so
            # the week-to-week march is directly comparable — diseased top-right →
            # healthy bottom-left.
            psmall = figures_dir / f"{exp}_main_08_{meta['prefix']}_entropy_weekly_state_panels.png"
            plot_state_time_small_multiples(summary_df, exp, psmall, burden_kind=burden_kind,
                                            treatment_group=treatment_group, comparator_group=comparator_group,
                                            comparison_df=comparison_df,
                                            normal_band=exp_res.get("normal_reference_band"))
            if psmall.exists():
                figure_paths[f"{exp}_main_08_{meta['prefix']}_entropy_weekly_state_panels"] = str(psmall)
                keep_main_files.append(psmall)

        finally:
            _dstate.CURRENT_DISPLAY_EXPERIMENT = previous_display_experiment
            _dstate.CURRENT_SAMPLE_LEVEL = previous_sample_level
            _dstate.CURRENT_SUMMARY = previous_summary

    _delete_unselected_figure_files(figures_dir, keep_main_files)
    bundle["figure_paths"] = figure_paths


def generate_output_bundle(results: Dict[str, Any], package_name: str = "BioEntropy_Results", explanation_options: Optional[Dict[str, Any]] = None, *, audit_context: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    old_group_order_map = _dstate.DISPLAY_GROUP_ORDER_MAP
    _dstate.DISPLAY_GROUP_ORDER_MAP = dict((explanation_options or {}).get("group_order_map") or results.get("group_order_map") or {})
    try:
        bundle = _legacy.generate_output_bundle(results, package_name=package_name)
        output_dir = Path(bundle["output_dir"])

        _append_joint_visual_artifacts(results, bundle)

        # Analysis_Report / explanation-report generation intentionally disabled.
        bundle["explanation_paths"] = {}

        # Per-figure numbers (bar heights and every scatter point) so the main
        # figures can be redrawn in GraphPad Prism from the exported values.
        export_issues = []
        try:
            from bioentropy_figure_data import export_figure_data_tables
            bundle["figure_data_files"] = [
                str(x) for x in export_figure_data_tables(results, output_dir)
            ]
        except Exception as exc:
            bundle["figure_data_files"] = []
            export_issues.append(f"figure_data: {type(exc).__name__}: {exc}")
            logging.getLogger(__name__).exception("Optional figure-data export failed")

        _refresh_zip(output_dir, Path(bundle["zip_path"]))
        # Sidecar metadata is deliberately outside output_dir and the legacy ZIP.
        # Any audit failure is reported without discarding the completed analysis.
        try:
            from bioentropy_audit import write_audit_bundle
            context = dict(audit_context or {})
            context["display_options"] = dict(explanation_options or {})
            context["display_options"]["group_order_map"] = dict(_dstate.DISPLAY_GROUP_ORDER_MAP)
            context["export_issues"] = list(context.get("export_issues", [])) + export_issues
            bundle.update(write_audit_bundle(results, output_dir, context))
        except Exception as exc:
            bundle["audit_error"] = f"{type(exc).__name__}: {exc}"
            logging.getLogger(__name__).exception("Analysis completed, but audit export failed")
        return bundle
    finally:
        _dstate.DISPLAY_GROUP_ORDER_MAP = old_group_order_map
