import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from metric_selection import score_field, method_name, LABELS, unit


def coordinates(result, metric, displays):
    groups = pd.DataFrame(result["summary"]).copy()
    points = pd.DataFrame(result["observations"]).copy()
    if "Week" not in points:
        points["Week"] = 0
    for axis, key, field in [("X", metric, score_field(metric)), ("Y", "NMD", "NMD")]:
        groups[axis] = groups[key]
        if score_field(key) is not None:
            points[axis] = points[field] if field in points else np.nan
        else:
            points[axis] = np.nan
    groups["Plot_scale"] = "raw"
    points["Plot_scale"] = "raw"
    groups["Contrast"] = "raw"
    points["Joint_plot_included"] = np.isfinite(points[["X", "Y"]]).all(axis=1)
    return (groups, points)


def limits(values, centers, quantiles=(0.0, 1.0)):
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    c = np.asarray(centers, float)
    c = c[np.isfinite(c)]
    q = (
        np.quantile(v, quantiles)
        if len(v)
        else np.array([c.min(), c.max()])
        if len(c)
        else np.array([0.0, 1.0])
    )
    lo = min(q[0], c.min() if len(c) else 0.0)
    hi = max(q[1], c.max() if len(c) else 100.0)
    span = hi - lo if hi > lo else max(abs(lo), 1.0) * 0.1
    return (float(lo - 0.05 * span), float(hi + 0.05 * span))


def draw(result, metric, displays, output, focus=False):
    from plots import group_color, LABELS, experiment_name
    from analysis_outputs import annotation
    from analysis_engine import csv

    centers, allpoints = coordinates(result, metric, displays)
    points = allpoints[allpoints.Joint_plot_included].copy()
    groups = list(dict.fromkeys(centers.Group))
    weeks = sorted(centers.Week.unique())
    scatter = score_field(metric) is not None
    xlim = limits([], centers.X) if focus else limits(points.X, centers.X, (0.1, 0.9))
    ylim = limits([], centers.Y) if focus else limits(points.Y, centers.Y, (0.1, 0.9))
    fullx = limits(points.X, centers.X)
    fully = limits(points.Y, centers.Y)
    rows = 2 if scatter else 1
    fig, axes = plt.subplots(
        rows,
        len(weeks),
        figsize=(max(12, 2.8 * len(weeks)), 9.8 if scatter else 7.5),
        squeeze=False,
        gridspec_kw={"height_ratios": [3.2, 1.1]} if scatter else {},
    )
    fig.subplots_adjust(
        left=0.055,
        right=0.985,
        top=0.73,
        bottom=0.22 if scatter else 0.26,
        hspace=0.77,
        wspace=0.12,
    )
    counts = []
    for col, week in enumerate(weeks):
        ax = axes[0, col]
        cw = centers[centers.Week.eq(week)]
        if len(cw) > 1:
            ax.plot(cw.X, cw.Y, color="#b2bac4", lw=1.1, zorder=2)
        for gi, group in enumerate(groups):
            color = group_color(group, gi)
            marker = "o" if gi % 2 == 0 else "^"
            part = cw[cw.Group.eq(group)]
            if part.empty:
                continue
            center = part.iloc[0]
            b = points[points.Group.eq(group) & points.Week.eq(week)]
            inside = b.X.between(*xlim) & b.Y.between(*ylim)
            ax.scatter(
                b.loc[inside, "X"],
                b.loc[inside, "Y"],
                s=8,
                alpha=0.23,
                color=color,
                marker=marker,
                zorder=1,
            )
            if np.isfinite(center[["X", "Y"]].to_numpy(float)).all():
                ax.scatter(
                    center.X,
                    center.Y,
                    s=135 if gi == 0 else 65,
                    facecolor="white" if gi == 0 else color,
                    edgecolor=color,
                    linewidth=2,
                    zorder=4 + gi,
                )
                other = cw[~cw.Group.eq(group)]
                above = (
                    gi == 0
                    if other.empty or np.isclose(center.Y, other.Y.iloc[0])
                    else center.Y > other.Y.iloc[0]
                )
                frac = (center.X - xlim[0]) / (xlim[1] - xlim[0])
                dx, ha = (
                    (4, "left") if frac < 0.25 else (-4, "right") if frac > 0.75 else (0, "center")
                )
                ax.annotate(
                    f"({center.X:.4g}, {center.Y:.4g})",
                    (center.X, center.Y),
                    xytext=(dx, 12 if above else -18),
                    textcoords="offset points",
                    ha=ha,
                    fontsize=7.5,
                    color=color,
                    zorder=8,
                    bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=0.3),
                )
            cc = dict(
                Experiment=result["name"],
                Group=group,
                Week=int(week),
                N_NMD=int(center.N_NMD),
                N_entropy=int(center.N_entropy),
                N_joint=len(b),
                N_outside_main=int((~inside).sum()),
                Left=int((b.X < xlim[0]).sum()),
                Right=int((b.X > xlim[1]).sum()),
                Below=int((b.Y < ylim[0]).sum()),
                Above=int((b.Y > ylim[1]).sum()),
                X_min=xlim[0],
                X_max=xlim[1],
                Y_min=ylim[0],
                Y_max=ylim[1],
            )
            counts.append(cc)
            if scatter:
                ax.text(
                    0.5,
                    -0.145 - 0.075 * gi,
                    f"←{cc['Left']}  →{cc['Right']}  ↑{cc['Above']}  ↓{cc['Below']}",
                    transform=ax.transAxes,
                    ha="center",
                    va="top",
                    color=color,
                    fontsize=7.5,
                )
                full = axes[1, col]
                full.scatter(b.X, b.Y, s=3, alpha=0.32, color=color, marker=marker)
                full.scatter(
                    center.X,
                    center.Y,
                    s=21,
                    facecolor="white" if gi == 0 else color,
                    edgecolor=color,
                    zorder=4,
                )
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.grid(alpha=0.15)
        ax.set_axisbelow(True)
        ax.set_title(str(week) + "W", fontweight="bold", fontsize=15, pad=49)
        note = (
            "NMD: "
            + annotation(result, "NMD", week, False, source="分周原值")
            + "\n"
            + LABELS[metric]
            + ": "
            + annotation(result, metric, week, False, source="分周原值")
        )
        ax.text(
            0.5,
            1.02,
            note,
            transform=ax.transAxes,
            ha="center",
            va="bottom",
            fontsize=9,
            color="#29374c",
            linespacing=1.6,
        )
        ax.tick_params(labelsize=8)
        if col:
            ax.tick_params(labelleft=False)
        if scatter:
            full = axes[1, col]
            full.set_xlim(*fullx)
            full.set_ylim(*fully)
            full.add_patch(
                Rectangle(
                    (xlim[0], ylim[0]),
                    xlim[1] - xlim[0],
                    ylim[1] - ylim[0],
                    fill=False,
                    edgecolor="#707b89",
                    lw=0.8,
                )
            )
            full.tick_params(labelsize=7)
            full.grid(alpha=0.12)
            if col:
                full.tick_params(labelleft=False)
    axes[0, 0].set_ylabel("NMD（原始值）", fontsize=12)
    if scatter:
        axes[1, 0].set_ylabel("完整分布\nNMD原始值", fontsize=9)
    fig.suptitle(
        f"Experiment {experiment_name(result['name'])}：分周 NMD–{LABELS[metric]} 联合状态图",
        fontsize=22,
        y=0.975,
    )
    subtitle = "趋势放大版" if focus else "兼顾散点分布版" if scatter else "组统计量视图"
    fig.text(
        0.5,
        0.925,
        subtitle + " · 原始数值 · 各周共用坐标 · 0W不平移、不对齐",
        ha="center",
        fontsize=12,
        color="#536273",
    )
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markeredgecolor=group_color(g, i),
            markerfacecolor="white" if i == 0 else group_color(g, i),
            markersize=9,
            markeredgewidth=2,
            label="原安慰剂组" if g == "安慰剂" else g,
        )
        for i, g in enumerate(groups)
    ]
    fig.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.905),
        ncol=min(4, len(groups)),
        frameon=False,
        fontsize=12,
    )
    xlabel = (
        LABELS[metric]
        + (
            " / 样本贡献分数"
            if metric == "Energy_distance"
            else " / 样本信息分数"
            if scatter
            else ""
        )
        + unit(metric)
        + " · 原始值 →"
    )
    fig.text(0.52, 0.157, xlabel, ha="center", fontsize=12)
    point_note = (
        "小点：同一患者同次访视的NMD与能量距离贡献分数。"
        if metric == "Energy_distance"
        else "小点：同一患者同次访视的NMD与信息分数。"
        if scatter
        else "本图显示所选方法的组统计量。"
    )
    fig.text(
        0.055, 0.112, "大点：组NMD与所选方法组估计。" + point_note, fontsize=9, color="#536273"
    )
    fig.text(
        0.055,
        0.08,
        "箭头为各方向超出主图的数量；角落外的同一点可能被两个方向计数。下方小图完整保留散点。"
        if scatter
        else "每个组点表示当周原始组估计。",
        fontsize=9,
        color="#536273",
    )
    roles = result.get("roles", {})
    fig.text(
        0.055,
        0.048,
        f"P值：{roles.get('treatment', '未指定')} vs {roles.get('comparator', '未指定')}，当周两组原始数值之差（含0W）；探索性患者轨迹置换，Holm校正；†见说明。",
        fontsize=9,
        color="#536273",
    )
    for key in ("NMD", metric):
        prefix = "NMD" if key == "NMD" else "Metric"
        found = {
            r["Week"]: r
            for r in result.get("comparisons_visit", [])
            if r["Metric"] == key and r["Contrast"] == "raw"
        }
        centers[prefix + "_P_Holm"] = centers.Week.map({w: r["P_Holm"] for w, r in found.items()})
        centers[prefix + "_Annotation"] = centers.Week.map(
            {w: r["Annotation"] for w, r in found.items()}
        )
    method = method_name(metric)
    num = 11 if focus else 10
    stem = f"{result['name']}_main_{num:02d}_weekly_state"
    target = output / "figures" / method / (stem + ".png")
    fig.savefig(target, dpi=180, facecolor="white")
    plt.close(fig)
    directory = output / "图形数据_GraphPad" / result["name"]
    csv(centers, directory / (stem + "_" + method + "_groups.csv"))
    csv(allpoints, directory / (stem + "_" + method + "_points.csv"))
    csv(counts, directory / (stem + "_" + method + "_counts.csv"))
    return dict(
        Experiment=result["name"],
        Method=method,
        Number=num,
        Title="分周联合图 · " + subtitle,
        File=target.relative_to(output).as_posix(),
        Data=(directory / (stem + "_" + method + "_groups.csv")).relative_to(output).as_posix(),
    )
