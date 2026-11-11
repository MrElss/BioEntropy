from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from metric_selection import LABELS, metrics, method_name, score_field, unit

COLORS = ["#315e9b", "#cf4939", "#8b749d", "#45998c", "#ad873c", "#636f83", "#cb765e"]


def setup():
    plt.rcParams.update(
        {
            "font.family": ["Microsoft YaHei", "DejaVu Sans"],
            "axes.unicode_minus": False,
            "font.size": 11,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "savefig.facecolor": "white",
        }
    )


def aligned_table(results):
    rows = []
    for r in results:
        if r["mode"] != "clinical":
            continue
        df = pd.DataFrame(r["summary"])
        for metric in metrics(r):
            for group, part in df.groupby("Group", sort=False):
                baseline = part.loc[part.Week == 0, metric].iloc[0]
                for row in part.itertuples():
                    value = getattr(row, metric)
                    rows.append(
                        dict(
                            Experiment=r["name"],
                            Group=group,
                            Week=row.Week,
                            Metric=metric,
                            Raw=value,
                            Baseline_raw=baseline,
                            Change=value - baseline,
                        )
                    )
    df = pd.DataFrame(rows)
    bounds = {}
    if df.empty:
        return (df, bounds)
    for metric in df.Metric.unique():
        mask = df.Metric.eq(metric)
        finite = df.loc[mask, "Change"].dropna()
        if finite.empty:
            lo, hi = (0.0, 0.0)
        else:
            lo, hi = (float(finite.min()), float(finite.max()))
        span = hi - lo
        if span > 0:
            df.loc[mask, "Display"] = 100 * (df.loc[mask, "Change"] - lo) / span
            baseline_display = -100 * lo / span
        else:
            df.loc[mask, "Display"] = np.where(df.loc[mask, "Change"].notna(), 50.0, np.nan)
            baseline_display = 50.0
        bounds[metric] = dict(
            delta_min=lo,
            delta_max=hi,
            span=span,
            baseline_display=baseline_display,
            scope="本次所选临床实验、两组、全部访视共用线性范围",
        )
    return (df, bounds)


def experiment_name(name):
    return {"105and106": "105＋106"}.get(name, name)


def group_color(group, i):
    s = str(group).lower()
    if "normal" in s:
        return "#919dbb"
    if "model" in s or "modle" in s or "安慰剂" in s:
        return "#e4191c"
    if "htd" in s:
        return "#4e6597"
    if s == "pm":
        return "#f5a18a"
    return ["#8b715a", "#b3a38d", "#45a0a3", "#8a76a6"][i % 4]


def display_table(results, aligned):
    frames = [aligned] if not aligned.empty else []
    for r in results:
        if r["mode"] != "animal":
            continue
        d = pd.DataFrame(r["summary"])
        for metric in metrics(r):
            v = d[metric].to_numpy(float)
            finite = v[np.isfinite(v)]
            lo = float(finite.min()) if len(finite) else np.nan
            hi = float(finite.max()) if len(finite) else np.nan
            span = hi - lo
            shown = 100 * (v - lo) / span if span > 0 else np.where(np.isfinite(v), 50.0, np.nan)
            frames.append(
                pd.DataFrame(
                    dict(
                        Experiment=r["name"],
                        Group=d.Group,
                        Week=d.Week,
                        Metric=metric,
                        Raw=v,
                        Display=shown,
                        Display_min=lo,
                        Display_max=hi,
                    )
                )
            )
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def data_for_plot(r, metric, normalized, displays):
    d = pd.DataFrame(r["summary"]).copy()
    if normalized:
        a = displays[(displays.Experiment == r["name"]) & (displays.Metric == metric)]
        d = d.merge(a[["Group", "Week", "Display"]], on=["Group", "Week"], validate="one_to_one")
    else:
        d["Display"] = d[metric]
    d["Raw"] = d[metric]
    d["Metric"] = metric
    d["N"] = d.N_NMD if metric == "NMD" else d.N_entropy
    return d


def draw_metric(ax, r, metric, normalized, displays):
    d = data_for_plot(r, metric, normalized, displays)
    clinical = r["mode"] == "clinical"
    if clinical:
        groups = list(dict.fromkeys(d.Group))
        series = {g: d[d.Group.eq(g)].sort_values("Week") for g in groups}
        for gi, g in enumerate(groups):
            b = series[g]
            color = group_color(g, gi)
            ax.plot(
                b.Week,
                b.Display,
                color=color,
                lw=2.7,
                marker="o",
                ms=7,
                label=g,
                markeredgecolor="white",
                markeredgewidth=1,
                zorder=4,
            )
            for row in b.itertuples():
                if not np.isfinite(row.Display):
                    continue
                others = [
                    float(z.loc[z.Week.eq(row.Week), "Display"].iloc[0])
                    for gg, z in series.items()
                    if gg != g
                ]
                above = (
                    gi == 0
                    if not others or np.isclose(row.Display, others[0])
                    else row.Display > others[0]
                )
                label = (
                    f"{row.Display:.4f}"
                    if metric != "NMD" and (not normalized)
                    else f"{row.Display:.3g}"
                )
                ax.annotate(
                    label,
                    (row.Week, row.Display),
                    xytext=(0, 10 if above else -18),
                    textcoords="offset points",
                    ha="center",
                    fontsize=9,
                    color=color,
                )
        ax.set_xticks(sorted(d.Week.unique()), [str(w) + "W" for w in sorted(d.Week.unique())])
        ax.set_xlabel("Week")
        ax.axvline(24, ls="--", lw=0.8, color="#b4bcc8")
        ax.legend(
            loc="upper center", bbox_to_anchor=(0.5, 1.19), ncol=min(3, len(groups)), frameon=False
        )
    else:
        colors = [group_color(g, i) for i, g in enumerate(d.Group)]
        ax.bar(np.arange(len(d)), d.Display, color=colors, width=0.64, edgecolor="white", zorder=3)
        ax.set_xticks(np.arange(len(d)), d.Group, rotation=24, ha="right")
        ax.set_xlabel("Treatment group")
        for i, row in enumerate(d.itertuples()):
            if np.isfinite(row.Display):
                label = (
                    f"{row.Display:.4f}"
                    if metric != "NMD" and (not normalized)
                    else f"{row.Display:.3g}"
                )
                ax.annotate(
                    label,
                    (i, row.Display),
                    xytext=(0, 7),
                    textcoords="offset points",
                    ha="center",
                    fontsize=10,
                )
            else:
                ax.text(i, 0, "NA", ha="center")
    finite = d.Display.dropna().to_numpy()
    if normalized:
        lo, hi = (-18.0, 124.0)
    elif len(finite):
        span = max(float(np.ptp(finite)), max(abs(finite)) * 0.08, 0.1)
        lo = (
            min(0.0, float(finite.min()) - 0.12 * span)
            if not clinical
            else float(finite.min()) - 0.18 * span
        )
        hi = float(finite.max()) + 0.25 * span
    else:
        lo, hi = (0.0, 1.0)
    ax.set_ylim(lo, hi)
    ax.grid(axis="y", alpha=0.18)
    ax.set_axisbelow(True)
    short = "NMD" if metric == "NMD" else LABELS[metric]
    ax.set_ylabel(
        short + (" 显示指数（0–100）" if normalized else "" if metric == "NMD" else unit(metric))
    )
    return d


def nmd_points(ax, r, normalized, displays, groups, metric="NMD"):
    obs = pd.DataFrame(r["observations"]).copy()
    if "Week" not in obs:
        obs["Week"] = 0
    field = score_field(metric)
    if field not in obs:
        return (pd.DataFrame(), False)
    obs = obs[np.isfinite(obs[field])].copy()
    obs["Display"] = obs[field]
    obs["Point_metric"] = field
    if normalized:
        a = displays[(displays.Experiment == r["name"]) & (displays.Metric == metric)]
        if r["mode"] == "clinical":
            bases = dict(
                a[["Group", "Baseline_raw"]].drop_duplicates().itertuples(index=False, name=None)
            )
            f = a.dropna(subset=["Change", "Display"])
            unique = f.drop_duplicates("Change")
            if len(unique) >= 2:
                u, v = (unique.iloc[0], unique.iloc[-1])
                slope = (v.Display - u.Display) / (v.Change - u.Change)
                offset = u.Display - slope * u.Change
            else:
                slope, offset = (0.0, 50.0)
            obs["Display"] = (obs[field] - obs.Group.map(bases)) * slope + offset
        else:
            lo = float(a.Display_min.iloc[0])
            hi = float(a.Display_max.iloc[0])
            obs["Display"] = 100 * (obs[field] - lo) / (hi - lo) if hi > lo else 50.0
    rng = np.random.default_rng(42)
    ymin, ymax = ax.get_ylim()
    outside = False
    for i, g in enumerate(groups):
        b = obs[obs.Group.eq(g)].copy()
        positions = b.Week.to_numpy(float) if r["mode"] == "clinical" else np.full(len(b), i, float)
        jitter = 0.34 if r["mode"] == "clinical" else 0.18
        b["x"] = positions + rng.uniform(-jitter, jitter, len(b))
        visible = b.Display.between(ymin, ymax)
        color = group_color(g, i)
        ax.scatter(
            b.loc[visible, "x"],
            b.loc[visible, "Display"],
            s=14 if r["mode"] == "clinical" else 20,
            color=color if r["mode"] == "clinical" else "#26374b",
            alpha=0.25 if r["mode"] == "clinical" else 0.55,
            edgecolor="none" if r["mode"] == "clinical" else "white",
            linewidth=0.5,
            zorder=2 if r["mode"] == "clinical" else 3.5,
        )
        for week, bb in b.groupby("Week"):
            pos = float(week) if r["mode"] == "clinical" else i
            for up, flag in [(True, bb.Display > ymax), (False, bb.Display < ymin)]:
                n = int(flag.sum())
                if n:
                    outside = True
                    shift = (-6 if i == 0 else 6) if r["mode"] == "clinical" else 0
                    ha = ("right" if i == 0 else "left") if r["mode"] == "clinical" else "center"
                    ax.annotate(
                        ("↑ " if up else "↓ ") + str(n),
                        (pos, ymax if up else ymin),
                        xytext=(shift, -7 if up else 7),
                        textcoords="offset points",
                        ha=ha,
                        va="top" if up else "bottom",
                        fontsize=7.5,
                        color=color,
                        zorder=7,
                        bbox=dict(facecolor="white", edgecolor="none", alpha=0.9, pad=1),
                    )
    if outside:
        inset = ax.inset_axes([1.065, 0.44, 0.3, 0.46])
        inset.set_title("全部散点", fontsize=9)
        for i, g in enumerate(groups):
            b = obs[obs.Group.eq(g)]
            x = b.Week if r["mode"] == "clinical" else np.full(len(b), i)
            inset.scatter(x, b.Display, s=3, color=group_color(g, i), alpha=0.4)
        inset.tick_params(labelsize=6)
        inset.grid(alpha=0.15)
    return (obs, outside)


def metric_figure(r, metric, normalized, displays, path):
    from analysis_outputs import annotate_metric

    clinical = r["mode"] == "clinical"
    is_nmd = metric == "NMD"
    has_points = score_field(metric) is not None
    fig, ax = plt.subplots(figsize=(14 if has_points else 12, 8.3))
    fig.subplots_adjust(
        left=0.08 if has_points else 0.1,
        right=0.74 if has_points else 0.97,
        top=0.72,
        bottom=0.29 if clinical else 0.33,
    )
    d = draw_metric(ax, r, metric, normalized, displays)
    points = pd.DataFrame()
    outside = False
    if has_points:
        points, outside = nmd_points(
            ax, r, normalized, displays, list(dict.fromkeys(d.Group)), metric
        )
    annotate_metric(ax, r, metric, normalized)
    label = LABELS[metric]
    view = "基线对齐显示" if normalized and clinical else "归一化显示" if normalized else "原始值"
    fig.suptitle(f"Experiment {experiment_name(r['name'])}：{label}", fontsize=20, y=0.965)
    fig.text(0.5, 0.9, view, ha="center", fontsize=12, color="#536273")
    note = (
        "各组减自身0W后共用线性尺度；起点重叠仅用于显示。"
        if normalized and clinical
        else "按各组估计值的共同范围线性映射。"
        if normalized
        else "原始估计值。"
        if is_nmd
        else "原始估计值。"
    )
    if outside:
        note += " 超出主图的散点用箭头计数，小图保留全部分布。"
    fig.text(0.08, 0.18, note, fontsize=8.5, color="#596879", wrap=True)
    point_note = (
        "小点：每个样本的NMD。"
        if is_nmd
        else "小点：样本信息分数Cᵢ；大点/柱高：组Shannon熵。"
        if has_points
        else "每组/访视一个所选方法的组估计。"
    )
    if metric == "Energy_distance":
        point_note = "小点：样本能量距离贡献分数；均值为组估计。"
    fig.text(0.08, 0.145, point_note, fontsize=9, color="#596879")
    contrast = "两组相对0W变化之差" if normalized and clinical else "当次两组原值之差"
    roles = r.get("roles", {})
    fig.text(
        0.08,
        0.11,
        f"标记：{roles.get('treatment', '未指定')} vs {roles.get('comparator', '未指定')}；{contrast}。探索性置换，Holm校正；†需谨慎核查。",
        fontsize=8.5,
        color="#596879",
    )
    counts = "；".join(
        str(g) + ": " + "/".join(str(int(x)) for x in b.N)
        for g, b in d.groupby("Group", sort=False)
    )
    fig.text(0.08, 0.065, "n：" + counts, fontsize=8, color="#596879")
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return (d, points)


def combined_figure(r, metric, displays, path):
    from analysis_outputs import annotate_metric

    fig, axes = plt.subplots(1, 2, figsize=(19, 7.6))
    fig.subplots_adjust(left=0.06, right=0.85, top=0.7, bottom=0.31, wspace=0.73)
    for ax, key in zip(axes, ["NMD", metric]):
        data = draw_metric(ax, r, key, True, displays)
        if score_field(key) is not None:
            nmd_points(ax, r, True, displays, list(dict.fromkeys(data.Group)), key)
        annotate_metric(ax, r, key, True)
    fig.suptitle(
        f"Experiment {experiment_name(r['name'])}：NMD 与 {LABELS[metric]}", fontsize=19, y=0.97
    )
    fig.text(
        0.5,
        0.89,
        "同一实验的两种指标 · 各自线性显示尺度",
        ha="center",
        fontsize=11,
        color="#596879",
    )
    fig.text(
        0.06,
        0.115,
        "同指标组统计量与对应样本分数共用线性变换。完整分布见小图。",
        fontsize=9,
        color="#596879",
    )
    fig.text(
        0.06,
        0.07,
        "各指标采用各自的线性显示尺度。标记仅针对指定两组；临床为两组相对0W变化之差，动物为原值之差；探索性置换、Holm校正。",
        fontsize=9,
        color="#596879",
    )
    fig.savefig(path, dpi=150)
    plt.close(fig)


def state_figure(r, metric, normalized, displays, path):
    x = data_for_plot(r, "NMD", normalized, displays)
    y = data_for_plot(r, metric, normalized, displays)
    d = (
        x[["Experiment", "Group", "Week", "Display"]]
        .rename(columns={"Display": "NMD_plot"})
        .merge(y[["Group", "Week", "Display"]], on=["Group", "Week"], validate="one_to_one")
        .rename(columns={"Display": "Entropy_plot"})
    )
    fig, ax = plt.subplots(figsize=(10, 7))
    fig.subplots_adjust(left=0.13, right=0.95, bottom=0.2, top=0.79)
    for i, (g, b) in enumerate(d.groupby("Group", sort=False)):
        color = group_color(g, i)
        b = b.sort_values("Week")
        ax.plot(
            b.NMD_plot,
            b.Entropy_plot,
            "o-" if r["mode"] == "clinical" else "o",
            color=color,
            label=g,
            lw=1.7,
            ms=7,
        )
        if r["mode"] == "clinical":
            for row in b.itertuples():
                ax.annotate(
                    str(row.Week) + "W",
                    (row.NMD_plot, row.Entropy_plot),
                    xytext=(5, 5),
                    textcoords="offset points",
                    fontsize=8,
                    color=color,
                )
    ax.set_xlabel("NMD" + (" 显示指数" if normalized else ""))
    ax.set_ylabel(LABELS[metric] + (" 显示指数" if normalized else unit(metric)))
    ax.grid(alpha=0.15)
    ax.margins(0.18)
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, 1.18),
        ncol=min(3, len(d.Group.unique())),
        fontsize=9,
        frameon=False,
    )
    fig.suptitle(
        f"Experiment {experiment_name(r['name'])}：NMD–{LABELS[metric]} 状态图", fontsize=17, y=0.97
    )
    fig.text(
        0.12,
        0.08,
        "每点对应一个组/时间点；用于查看两项指标的联合变化。",
        fontsize=9,
        color="#596879",
    )
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return d


def create_report(output, results, display="Shannon"):
    from analysis_engine import csv, write_json
    from table_exports import export_tables, write_viewer

    setup()
    output = Path(output)
    aligned, bounds = aligned_table(results)
    all_display = display_table(results, aligned)
    csv(all_display, output / "显示值与原始值.csv")
    if not aligned.empty:
        csv(aligned, output / "临床_原始值与显示值.csv")
        write_json(output / "临床绘图尺度.json", bounds)
    figures = output / "figures"
    figures.mkdir()
    graph = output / "图形数据_GraphPad"
    graph.mkdir()
    index = []
    for r in results:
        directory = graph / r["name"]
        directory.mkdir()
        for metric in metrics(r):
            method = method_name(metric)
            dest = figures / method
            dest.mkdir(exist_ok=True)
            nums = (1, 2) if metric == "NMD" else (3, 4)
            for normalized, num in zip((False, True), nums):
                stem = (
                    f"{r['name']}_main_{num:02d}_{method}_{('normalized' if normalized else 'raw')}"
                )
                target = dest / (stem + ".png")
                data, points = metric_figure(r, metric, normalized, all_display, target)
                csv(data, directory / (stem + "_bars.csv"))
                wide = (
                    data.pivot(index="Week", columns="Group", values="Display")
                    if r["mode"] == "clinical"
                    else data.set_index("Group")[["Display"]]
                )
                csv(
                    wide.reindex(columns=list(dict.fromkeys(data.Group))).reset_index()
                    if r["mode"] == "clinical"
                    else wide.reset_index(),
                    directory / (stem + "_prism.csv"),
                )
                if len(points):
                    csv(points, directory / (stem + "_points.csv"))
                index.append(
                    dict(
                        Experiment=r["name"],
                        Method=method,
                        Number=num,
                        Title=(
                            "原始值"
                            if not normalized
                            else "基线对齐"
                            if r["mode"] == "clinical"
                            else "归一化显示"
                        )
                        + " · "
                        + LABELS[metric],
                        File=target.relative_to(output).as_posix(),
                        Data=(directory / (stem + "_bars.csv")).relative_to(output).as_posix(),
                    )
                )
            if metric == "NMD":
                continue
            target = dest / f"{r['name']}_main_05_NMD_{method}_comparison.png"
            combined_figure(r, metric, all_display, target)
            index.append(
                dict(
                    Experiment=r["name"],
                    Method=method,
                    Number=5,
                    Title="NMD与" + LABELS[metric] + "对照图",
                    File=target.relative_to(output).as_posix(),
                    Data="显示值与原始值.csv",
                )
            )
            for normalized, num in [(False, 7), (True, 8)] if r["mode"] == "animal" else []:
                target = dest / f"{r['name']}_main_{num:02d}_NMD_{method}_state.png"
                data = state_figure(r, metric, normalized, all_display, target)
                dp = directory / (target.stem + ".csv")
                csv(data, dp)
                index.append(
                    dict(
                        Experiment=r["name"],
                        Method=method,
                        Number=num,
                        Title="NMD–"
                        + LABELS[metric]
                        + "状态图"
                        + ("（显示值）" if normalized else "（原始值）"),
                        File=target.relative_to(output).as_posix(),
                        Data=dp.relative_to(output).as_posix(),
                    )
                )
    csv(index, graph / "图表-数据对照表.csv")
    write_json(output / "图表索引.json", index)
    from analysis_engine import VERSION

    write_json(
        output / "结果索引.json",
        dict(
            version=VERSION,
            layout="desktop",
            experiments=[{k: r[k] for k in ("name", "mode", "folder")} for r in results],
        ),
    )
