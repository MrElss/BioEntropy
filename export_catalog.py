from pathlib import Path
import json
import numpy as np
import pandas as pd
from metric_selection import metrics, score_field, LABELS


def write_traceability(output, results):
    from analysis_engine import csv, write_json
    from table_exports import workbook
    from analysis_outputs import annotation

    output = Path(output)
    index = json.loads((output / "图表索引.json").read_text(encoding="utf-8"))
    lookup = {r["name"]: r for r in results}
    elements, plotted_tests, counts, different = ([], [], [], [])
    for item in index:
        result = lookup[item["Experiment"]]
        clinical = result["mode"] == "clinical"
        number = item["Number"]
        source = "分周原值" if number in (10, 11) else "常规"
        normalized = number in (2, 4, 5, 8) or number == 6
        contrast = "change" if clinical and source == "常规" and normalized else "raw"
        item.update(
            Category="auxiliary" if number in (6, 9) else "main",
            Figure_policy=source,
            P_contrast=contrast if number in (1, 2, 3, 4, 5, 10, 11) else "none",
            Scale="baseline-aligned linear display"
            if clinical and normalized
            else "linear display"
            if normalized
            else "raw",
            P_data="图形数据_GraphPad/每图显著性.csv",
        )
        detail = output / result["folder"]
        pfile = detail / ("探索性置换检验_" + source + "方案.csv")
        tests = (
            result.get("comparisons_visit", [])
            if source == "分周原值"
            else result.get("comparisons", [])
        )
        groupfile = output / item["Data"]

        def add(element, file, columns, calculation, note=""):
            path = Path(file)
            elements.append(
                dict(
                    实验=result["name"],
                    图=f"图{number:02d}",
                    图中元素=element,
                    文件=path.relative_to(output).as_posix(),
                    列名=columns,
                    计算方式=calculation,
                    图表方案=source,
                    软件版本="V1.0.2",
                    说明=note,
                )
            )

        selected = metrics(result)[1:]
        figure_metrics = (
            ["NMD"] if number in (1, 2) else selected if number in (3, 4) else metrics(result)
        )
        if number in (1, 2, 3, 4):
            metric = figure_metrics[0]
            value = (
                "有效样本NMD的算术平均。"
                if metric == "NMD"
                else "读取所选方法的组估计，非样本分数的独立检验。"
            )
            scale = (
                (
                    "各组减自身0W后按本次临床组估计变化范围共用线性尺度。"
                    if clinical
                    else "各组估计按本实验该指标的共同最小值和最大值映射0–100。"
                )
                if normalized
                else "使用原始计算值。"
            )
            add("柱高或折线点", groupfile, f"Group; Week; {metric}; Raw; Display", value + scale)
            pointfile = groupfile.with_name(groupfile.name.replace("_bars.csv", "_points.csv"))
            if pointfile.is_file():
                field = score_field(metric)
                description = (
                    "每个有效样本的NMD。"
                    if metric == "NMD"
                    else "Ci=−∫Kh(z−zi)ln f̂(z)dz，是样本信息分数。"
                    if metric == "Entropy_Shannon"
                    else "本方法的样本贡献分数，均值对应组能量距离。"
                )
                add(
                    "散点",
                    pointfile,
                    f"{field}; Display; Group; Week; USUBJID或Sample_ID",
                    description
                    + (
                        "与组统计量共用相同线性变换，不截断。"
                        if normalized
                        else "原始分数，不作显示归一化。"
                    ),
                )
            add(
                "误差棒",
                detail / "组结果.csv",
                "无生物学误差棒",
                "这些主图不使用数值积分误差作为生物学标准误或95%置信区间。",
            )
        elif number == 5:
            for metric in figure_metrics:
                add(
                    "组点 " + LABELS[metric],
                    groupfile,
                    f"Metric={metric}; Raw; Display",
                    "与图02或04使用相同的组值和显示变换，两个指标各有自己的尺度。",
                )
                n = 2 if metric == "NMD" else 4
                original = next(
                    x for x in index if x["Experiment"] == result["name"] and x["Number"] == n
                )
                points = output / Path(original["Data"]).with_name(
                    Path(original["Data"]).name.replace("_bars.csv", "_points.csv")
                )
                if points.is_file():
                    add(
                        "散点 " + LABELS[metric],
                        points,
                        f"{score_field(metric)}; Display; Group; Week",
                        f"与图{n:02d}完全相同的样本字段和显示变换。",
                    )
        elif number in (7, 8):
            add(
                "组状态点",
                groupfile,
                "NMD_plot; Entropy_plot; Group; Week",
                "横轴为NMD，纵轴为所选方法组估计；"
                + ("使用原始值。" if number == 7 else "使用线性显示指数。"),
                "每点为一个组，没有个体散点或自动P值；不作联合显著性检验。",
            )
        elif number in (10, 11):
            add(
                "大点",
                groupfile,
                "X; Y; NMD; " + selected[0] + "; NMD_P_Holm; Metric_P_Holm",
                "X为所选方法原始组估计，Y为原始NMD组均值；0W不平移。",
                "共同样本用于两项组估计。"
                if result.get("shared_shannon_sample")
                else "组估计使用各指标的有效记录。",
            )
            points = groupfile.with_name(groupfile.name.replace("_groups.csv", "_points.csv"))
            add(
                "小散点",
                points,
                "USUBJID; Week; X; Y; Joint_plot_included",
                "同一患者同次访视的原始样本分数与NMD，同时有效才进入联合图。",
                "按所选方法展示组估计及可用的样本分数。",
            )
            limits = groupfile.with_name(groupfile.name.replace("_groups.csv", "_counts.csv"))
            add(
                "箭头及完整分布小图",
                limits,
                "N_joint; N_outside_main; Left; Right; Above; Below",
                "统计有效联合点超出局部窗口的方向与数量，小图保留完整有效分布。",
                "一个角落外点可被两个方向计数；不表示从计算中删掉样本。",
            )
        elif number == 6:
            add(
                "交互组估计",
                groupfile,
                "week; group; metric; value; raw",
                "value为既定线性显示值，raw为原始组估计；图中无个体熵或自动P值。",
            )
        elif number == 9:
            add(
                "热图颜色",
                groupfile,
                "Feature; Group; Week; Mean_standardized; Center; SD",
                "本实验各指标统一标准化后计算组/访视均值。",
                "仅指标概况，展示指标概况。",
            )
        if number in (1, 2, 3, 4, 5, 10, 11):
            for metric in figure_metrics:
                add(
                    "显著性 " + LABELS[metric],
                    pfile,
                    f"Metric={metric}; Week; Contrast={contrast}; Inference_source={source}; P_Holm; Annotation",
                    (
                        "两组相对0W变化之差"
                        if contrast == "change"
                        else "当周或当次两组原始组统计量之差"
                    )
                    + "；探索性置换及Holm校正。",
                    "仅比较所选治疗组和比较组。使用P_Holm；保留†并查看Caution。",
                )
                for week in sorted({row["Week"] for row in result["summary"]}):
                    matched = next(
                        (
                            r
                            for r in tests
                            if r["Metric"] == metric
                            and r["Week"] == week
                            and (r["Contrast"] == contrast)
                        ),
                        None,
                    )
                    label = annotation(result, metric, week, normalized, source=source)
                    record = (
                        dict(matched)
                        if matched
                        else dict(
                            Experiment=result["name"],
                            Metric=metric,
                            Week=week,
                            Contrast=contrast,
                            Inference_source=source,
                            P_Holm=np.nan,
                        )
                    )
                    record.update(
                        Figure=number,
                        Plot_annotation=label,
                        Figure_file=item["File"],
                        P_source=pfile.relative_to(output).as_posix(),
                    )
                    plotted_tests.append(record)
    for result in results:
        obs = pd.DataFrame(result["observations"])
        if "Week" not in obs:
            obs["Week"] = 0
        selected = metrics(result)[1:]
        field = score_field(selected[0]) if selected else None
        for (group, week), part in obs.groupby(["Group", "Week"], sort=False):
            nmd = np.isfinite(part.NMD)
            row = dict(
                Experiment=result["name"],
                Group=group,
                Week=int(week),
                N_present=len(part),
                N_NMD=int(nmd.sum()),
            )
            if selected:
                usable = (
                    part.Included_entropy.to_numpy(bool)
                    if "Included_entropy" in part
                    else np.isfinite(part[field]).to_numpy()
                    if field in part
                    else np.zeros(len(part), bool)
                )
                score = (
                    np.isfinite(part[field]).to_numpy()
                    if field and field in part
                    else np.zeros(len(part), bool)
                )
                row.update(
                    Selected_metric=selected[0],
                    N_selected=int(usable.sum()),
                    N_scores=int(score.sum()),
                    N_common_analysis=int((nmd.to_numpy() & usable).sum()),
                    N_joint_scatter=int((nmd.to_numpy() & score).sum()),
                )
                for i, (_, original) in enumerate(part.iterrows()):
                    if bool(nmd.iloc[i]) != bool(usable[i]):
                        different.append(
                            dict(
                                Experiment=result["name"],
                                Group=group,
                                Week=int(week),
                                Subject=original.get("USUBJID", original.get("Sample_ID", "")),
                                NMD_valid=bool(nmd.iloc[i]),
                                Selected_valid=bool(usable[i]),
                                Missing_selected_features=original.get(
                                    "Missing_entropy_features", ""
                                ),
                                NMD_status=original.get("NMD_status", "见NMD面板及原始记录"),
                                Source_file=original.get("Source_file", ""),
                                Excel_row=original.get("Excel_row", ""),
                            )
                        )
            counts.append(row)
    order = {r["name"]: i for i, r in enumerate(results)}
    index.sort(key=lambda z: (order[z["Experiment"]], z["Category"] != "main", z["Number"]))
    write_json(output / "图表索引.json", index)
    csv(index, output / "图形数据_GraphPad/图表-数据对照表.csv")
    csv(elements, output / "图形数据_GraphPad/图表元素说明.csv")
    csv(plotted_tests, output / "图形数据_GraphPad/每图显著性.csv")
    workbook(
        output / "数据表/16_图表元素与来源.xlsx",
        {"图元素说明": elements, "每图显著性": plotted_tests, "图文件索引": index},
        translated=False,
    )
    workbook(
        output / "数据表/17_逐组纳入人数.xlsx",
        {"逐组人数": counts, "纳入差异记录": different},
        translated=False,
    )
