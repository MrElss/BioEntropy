from pathlib import Path
import json, html, warnings
import numpy as np
import pandas as pd
from scipy.stats import chi2, ttest_ind
from sklearn.ensemble import IsolationForest
from sklearn.neighbors import LocalOutlierFactor
from sklearn.covariance import MinCovDet
from clinical_inputs import InputSelectionError
from clinical_patient_ids import flat_read
from experiment_config import NONE
from metric_selection import metrics, method_name, score_field


def range_bounds(paths):
    bounds = {}
    for path in paths:
        h, rows, _ = flat_read(path)
        names = [str(x).strip() for x in h]
        lower = [x.lower() for x in names]
        records = []
        if "normal_min" in lower and "normal_max" in lower:
            feature = next(
                (
                    i
                    for i, x in enumerate(lower)
                    if x in ("feature", "biomarker", "indicator", "指标")
                ),
                None,
            )
            if feature is None:
                raise InputSelectionError(Path(path).name + "：正常范围表缺少Feature列。")
            records = [
                (
                    str(r[feature]).strip(),
                    r[lower.index("normal_min")],
                    r[lower.index("normal_max")],
                )
                for _, r in rows
            ]
        else:
            limits = {str(r[0]).strip().lower(): r for _, r in rows}
            if not {"normal_min", "normal_max"} <= set(limits):
                raise InputSelectionError(Path(path).name + "：需要Normal_min和Normal_max。")
            records = [
                (names[j], limits["normal_min"][j], limits["normal_max"][j])
                for j in range(1, len(names))
                if names[j].lower() not in ("sample id", "sample_id", "none", "")
            ]
        for f, lo, hi in records:
            try:
                lo, hi = (float(lo), float(hi))
            except (TypeError, ValueError):
                raise InputSelectionError(Path(path).name + "：" + f + "范围需填写数值。")
            if not np.isfinite([lo, hi]).all() or lo >= hi:
                raise InputSelectionError(Path(path).name + "：" + f + "下限必须小于上限。")
            if f in bounds and bounds[f] != (lo, hi):
                raise InputSelectionError(f + "：多个正常范围文件的界限冲突。")
            bounds[f] = (lo, hi)
    return bounds


def review_flags(data):
    x = np.asarray(data["values"], float)
    groups = np.asarray(data["groups"])
    weeks = np.asarray(data["weeks"])
    n, p = x.shape
    values = {
        k: np.full(n, np.nan) for k in ["MAD4", "IQR1_5", "IsolationForest", "LOF", "MCD_chi2_999"]
    }
    status = []
    for g in dict.fromkeys(groups):
        for w in sorted(set(weeks)):
            ii = np.flatnonzero((groups == g) & (weeks == w))
            a = x[ii]
            ok = np.isfinite(a).all(1)
            for key in ("MAD4", "IQR1_5"):
                if len(ii) < 4:
                    status.append(
                        dict(Group=g, Week=int(w), Method=key, Status="NA：组内记录不足4")
                    )
                    continue
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", RuntimeWarning)
                    if key == "MAD4":
                        med = np.nanmedian(a, 0)
                        scale = 1.4826 * np.nanmedian(np.abs(a - med), 0)
                        valid = np.isfinite(scale) & (scale > 0)
                        flags = np.abs(a[:, valid] - med[valid]) > 4 * scale[valid]
                    else:
                        q1, q3 = np.nanpercentile(a, [25, 75], axis=0)
                        span = q3 - q1
                        valid = np.isfinite(span) & (span > 0)
                        flags = (a[:, valid] < q1[valid] - 1.5 * span[valid]) | (
                            a[:, valid] > q3[valid] + 1.5 * span[valid]
                        )
                if valid.any():
                    values[key][ii] = np.where(
                        np.isfinite(a[:, valid]).any(1), flags.any(1).astype(float), np.nan
                    )
                status.append(
                    dict(
                        Group=g,
                        Week=int(w),
                        Method=key,
                        Status="仅标记；纳入有效尺度列"
                        if valid.any()
                        else "NA：无可用非零尺度指标",
                    )
                )
            aa = a[ok]
            idx = ii[ok]
            if len(aa) >= 5:
                scale = aa.std(0, ddof=1)
                valid = scale > 0
                z = (aa[:, valid] - aa[:, valid].mean(0)) / scale[valid]
                if z.shape[1]:
                    values["IsolationForest"][idx] = (
                        IsolationForest(
                            random_state=42, contamination="auto", n_jobs=1
                        ).fit_predict(z)
                        == -1
                    ).astype(float)
                    values["LOF"][idx] = (
                        LocalOutlierFactor(n_neighbors=min(20, len(z) - 1)).fit_predict(z) == -1
                    ).astype(float)
                    for key in ["IsolationForest", "LOF"]:
                        status.append(
                            dict(
                                Group=g,
                                Week=int(w),
                                Method=key,
                                Status="仅标记；完整样本、共同标准化",
                            )
                        )
                    if (
                        len(z) > 2 * z.shape[1] + 1
                        and np.linalg.matrix_rank(z - z.mean(0)) == z.shape[1]
                    ):
                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore")
                            fit = MinCovDet(random_state=42).fit(z)
                        values["MCD_chi2_999"][idx] = (
                            fit.mahalanobis(z) > chi2.ppf(0.999, z.shape[1])
                        ).astype(float)
                        status.append(
                            dict(
                                Group=g,
                                Week=int(w),
                                Method="MCD_chi2_999",
                                Status="仅标记；自由度=有效指标数",
                            )
                        )
                    else:
                        status.append(
                            dict(
                                Group=g,
                                Week=int(w),
                                Method="MCD_chi2_999",
                                Status="NA：样本/维数或满秩条件不满足",
                            )
                        )
                else:
                    status.append(
                        dict(Group=g, Week=int(w), Method="多变量检测", Status="NA：无变异指标")
                    )
            else:
                status.append(
                    dict(Group=g, Week=int(w), Method="多变量检测", Status="NA：完整样本不足5")
                )
    table = pd.DataFrame(values)
    table["Available_methods"] = table.notna().sum(1)
    table["Flag_count"] = table[list(values)].fillna(0).sum(1).astype(int)
    table["Review_only"] = np.where(table.Flag_count >= 2, "多方法提示，待核查", "")
    return (table, pd.DataFrame(status))


def extra_tables(output, r, e, progress):
    from analysis_engine import csv
    from table_exports import workbook

    data = r["_data"]
    obs = pd.DataFrame(r["observations"])
    source_obs = (
        sorted(
            r["observations"],
            key=lambda z: (e.weeks.index(z.get("Week", 0)), z.get("Excel_row", 0)),
        )
        if e.kind == "animal"
        else None
    )
    identity = pd.DataFrame(r.get("_source_observations", r["observations"]))
    progress(e.name + "：生成异常记录提示（仅标记，不删值、不收缩）")
    flags, status = review_flags(data)
    if len(identity) != len(flags):
        raise InputSelectionError("样本核查行数不一致。")
    flags = pd.concat([identity.reset_index(drop=True), flags], axis=1)
    detail = output / r["folder"]
    csv(flags, detail / "异常记录提示.csv")
    csv(status, detail / "检测方法适用情况.csv")
    bounds = range_bounds(e.range_files) if e.range_files else {}
    ranges = []
    for f, (lo, hi) in bounds.items():
        if f not in data["features"]:
            ranges.append(dict(Feature=f, Status="NA：当前指标面板没有同名同单位列，未自动换算"))
            continue
        j = data["features"].index(f)
        for i, v in enumerate(np.asarray(data["values"])[:, j]):
            if np.isfinite(v):
                ranges.append(
                    dict(
                        Group=data["groups"][i],
                        Week=data["weeks"][i],
                        Feature=f,
                        Value=v,
                        Normal_min=lo,
                        Normal_max=hi,
                        Source_file=identity.iloc[i].get("Source_file", Path(e.files[0]).name),
                        Excel_row=identity.iloc[i].get("Excel_row"),
                        Status="低于范围" if v < lo else "高于范围" if v > hi else "范围内",
                    )
                )
    if not bounds:
        ranges = [dict(Status="可载入正常范围文件进行指标核查；NMD采用Normal样本。")]
    workbook(
        detail / "补充核查.xlsx",
        {"异常提示": flags, "方法适用情况": status, "正常范围核查": ranges},
    )
    contrasts = []
    if e.treatment != NONE and e.comparator != NONE:
        for w in e.weeks:
            mask = obs.Week.eq(w) if "Week" in obs else np.ones(len(obs), bool)
            a = obs.loc[mask & obs.Group.eq(e.treatment), "NMD"].dropna().to_numpy(float)
            b = obs.loc[mask & obs.Group.eq(e.comparator), "NMD"].dropna().to_numpy(float)
            g = p = np.nan
            reason = "两组各至少2个有效NMD才作比较"
            if min(len(a), len(b)) >= 2:
                dof = len(a) + len(b) - 2
                sp = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / dof)
                if sp > 0:
                    g = (1 - 3 / (4 * dof - 1)) * (a.mean() - b.mean()) / sp
                    p = float(ttest_ind(a, b, equal_var=False).pvalue)
                    reason = "独立组；按访视分别计算；未作多重比较校正"
            contrasts.append(
                dict(
                    Experiment=e.name,
                    Week=w,
                    Treatment=e.treatment,
                    Comparator=e.comparator,
                    N_treatment=len(a),
                    N_comparator=len(b),
                    NMD_Hedges_g=g,
                    NMD_Welch_p=p,
                    Note=reason,
                )
            )
    workbook(detail / "治疗组与比较组_NMD.xlsx", {"NMD比较": contrasts})


def extra_figures(output, results):
    from plots import setup, group_color, LABELS
    from analysis_engine import csv, write_json
    import matplotlib.pyplot as plt

    setup()
    output = Path(output)
    displays = pd.read_csv(output / "显示值与原始值.csv")
    index = json.loads((output / "图表索引.json").read_text("utf-8"))
    for r in results:
        name = r["name"]
        data = r["_data"]
        raw = np.asarray(data["values"], float)
        groups = np.asarray(data["groups"])
        weeks = np.asarray(data["weeks"])
        order = list(dict.fromkeys(x["Group"] for x in r["summary"]))
        center = np.nanmean(raw, 0)
        sd = np.nanstd(raw, 0, ddof=1)
        z = (raw - center) / np.where(sd > 0, sd, np.nan)
        means = []
        labels = []
        records = []
        for g in order:
            for w in sorted(set(weeks)):
                a = z[(groups == g) & (weeks == w)]
                v = np.nanmean(a, 0) if len(a) else np.full(raw.shape[1], np.nan)
                means.append(v)
                labels.append(g + (f" · {w}W" if r["mode"] == "clinical" else ""))
                records.extend(
                    dict(
                        Group=g,
                        Week=int(w),
                        Feature=f,
                        Mean_standardized=v[j],
                        Center=center[j],
                        SD=sd[j],
                    )
                    for j, f in enumerate(data["features"])
                )
        fig, ax = plt.subplots(
            figsize=(max(10, len(data["features"]) * 0.75), max(5, len(labels) * 0.27))
        )
        finite = np.asarray(means)
        limit = max(1, float(np.nanmax(np.abs(finite))) if np.isfinite(finite).any() else 1)
        im = ax.imshow(finite, aspect="auto", cmap="RdBu_r", vmin=-limit, vmax=limit)
        ax.set_yticks(range(len(labels)), labels)
        ax.set_xticks(
            range(len(data["features"])), data["features"], rotation=55, ha="right", fontsize=8
        )
        ax.set_title(name + "：指标概况（各指标均值）")
        fig.colorbar(im, ax=ax, label="本次实验共同尺度的标准化均值")
        fig.tight_layout()
        target = output / "figures/NMD" / f"{name}_main_09_feature_profile.png"
        fig.savefig(target, dpi=150)
        plt.close(fig)
        dp = output / "图形数据_GraphPad" / name / "09_指标概况.csv"
        csv(records, dp)
        index.append(
            dict(
                Experiment=name,
                Method="NMD",
                Number=9,
                Title="指标概况（不作熵贡献分解）",
                File=target.relative_to(output).as_posix(),
                Data=dp.relative_to(output).as_posix(),
            )
        )
        for metric in metrics(r)[1:]:
            method = method_name(metric)
            d = displays[
                (displays.Experiment.astype(str) == name) & displays.Metric.isin(["NMD", metric])
            ]
            records = []
            for w in sorted({row["Week"] for row in r["summary"]}):
                for g in order:
                    for m in ["NMD", metric]:
                        part = d[(d.Group == g) & (d.Week == w) & (d.Metric == m)]
                        if len(part) and np.isfinite(part.Display.iloc[0]):
                            records.append(
                                dict(
                                    week=int(w),
                                    group=g,
                                    metric=LABELS[m],
                                    value=float(part.Display.iloc[0]),
                                    raw=float(part.Raw.iloc[0]),
                                )
                            )
            target = output / "figures" / method / f"{name}_main_06_interactive.html"
            interactive(target, name, records)
            dp = output / "图形数据_GraphPad" / name / f"06_{method}_三维显示.csv"
            csv(records, dp)
            index.append(
                dict(
                    Experiment=name,
                    Method=method,
                    Number=6,
                    Title="三维交互对照（拖动旋转）",
                    File=target.relative_to(output).as_posix(),
                    Data=dp.relative_to(output).as_posix(),
                )
            )
            if r["mode"] == "clinical":
                from weekly_state import draw

                index.append(draw(r, metric, displays, output))
                if score_field(metric) is not None:
                    index.append(draw(r, metric, displays, output, focus=True))
    write_json(output / "图表索引.json", index)
    csv(index, output / "图形数据_GraphPad/图表-数据对照表.csv")
    export_context_tables(output, results)


def export_context_tables(output, results):
    from table_exports import workbook
    from scipy.stats import spearmanr

    within = []
    normal = []
    configuration = []
    signatures = []
    for r in results:
        d = pd.DataFrame(r["summary"])
        cfg = json.loads((output / r["folder"] / "计算配置与来源.json").read_text("utf-8"))
        for metric in metrics(r)[1:]:
            b = d[["NMD", metric]].dropna()
            rho = (
                float(spearmanr(b.NMD, b[metric]).statistic)
                if len(b) >= 3 and b.NMD.nunique() > 1 and (b[metric].nunique() > 1)
                else np.nan
            )
            within.append(
                dict(
                    Experiment=r["name"],
                    Metric=metric,
                    N_group_visits=len(b),
                    Spearman_rho=rho,
                    Note="仅描述组/访视点的关联；无独立样本p值，不构成独立药效证据",
                )
            )
        n = cfg.get("normal", {})
        normal.append(
            dict(
                Experiment=r["name"],
                Reference=r["roles"]["reference"],
                N=n.get("n", n.get("Normal_n")),
                Estimator="LedoitWolf",
                Status=n.get("status", "可计算"),
                Features="; ".join(
                    cfg.get("nmd_features", n.get("features", cfg.get("features", [])))
                ),
            )
        )
        user = cfg.get("user_configuration", {})
        configuration.append(
            dict(
                Experiment=r["name"],
                Treatment=r["roles"]["treatment"],
                Comparator=r["roles"]["comparator"],
                Normal_reference=r["roles"]["reference"],
                Display_order=" → ".join(user.get("order", [])),
            )
        )
        signatures.append(
            (
                r,
                normal[-1]["Features"],
                json.dumps(cfg.get("energy_scale", cfg.get("entropy_scale", {})), sort_keys=True),
            )
        )
    pooled = []
    for i, (a, af, as_) in enumerate(signatures):
        for b, bf, bs in signatures[i + 1 :]:
            reasons = []
            if "105and106" in (a["name"], b["name"]) and (
                a["name"] in ("105", "106") or b["name"] in ("105", "106")
            ):
                reasons.append("合并数据与单试验重叠")
            if af != bf:
                reasons.append("NMD指标面板不同")
            if as_ != bs:
                reasons.append("熵坐标尺度不同")
            reasons.append("未声明跨实验可交换性，不自动合并")
            pooled.append(
                dict(
                    Experiment_A=a["name"],
                    Experiment_B=b["name"],
                    Pooled_correlation=np.nan,
                    Reason="；".join(reasons),
                )
            )
    workbook(output / "数据表/05_实验内相关性.xlsx", {"描述性相关": within})
    workbook(output / "数据表/06_跨实验合并核查.xlsx", {"合并条件核查": pooled})
    workbook(output / "数据表/07_正常参照概览.xlsx", {"正常参照": normal})
    workbook(output / "数据表/11_人工设置.xlsx", {"组别与顺序": configuration})


def interactive(path, name, records):
    data = json.dumps(records, ensure_ascii=False, allow_nan=False).replace("</", "<\\/")
    metric_label = next(
        (x["metric"] for x in records if x["metric"] != "NMD（正常参照马氏距离）"), "所选方法"
    )
    page = (
        '<!doctype html><meta charset="utf-8"><title>BioEntropy '
        "三维交互</title><style>body{font-family:Microsoft "
        "YaHei;background:#edf1f7;color:#18334e;margin:24px}canvas{backgro"
        "und:white;max-width:100%;touch-action:none}button,select{padding:"
        "8px;margin:8px}p{line-height:1.7}</style><h2>TITLE · "
        "NMD与熵</h2><p>两项指标各用自己的线性显示尺度；柱高采用线性显示尺度。拖动旋转，选择访视查看；原始值见下表。</p><l"
        'abel>访视 <select id="visit"></select></label><button '
        'id="reset">重置视角</button><br><canvas id="chart" width="1100" '
        'height="650"></canvas><table id="table" '
        'cellpadding="8"></table><script>\nconst '
        "data=DATA,c=document.getElementById('chart'),ctx=c.getContext('2d"
        "'),v=document.getElementById('visit');let "
        "angle=-.65,tilt=.35,drag=null;\n[...new "
        "Set(data.map(x=>x.week))].forEach(w=>{const "
        "o=document.createElement('option');o.value=w;o.textContent=w+'W';"
        "v.append(o)});\nfunction proj(x,y,z){const "
        "a=x*Math.cos(angle)-y*Math.sin(angle),b=x*Math.sin(angle)+y*Math."
        "cos(angle);return [550+a*68,420+b*30*tilt-z*2.5]}\nfunction "
        "poly(points,color){ctx.beginPath();points.forEach((q,i)=>{let "
        "p=proj(...q);i?ctx.lineTo(...p):ctx.moveTo(...p)});ctx.closePath("
        ");ctx.fillStyle=color;ctx.fill();ctx.strokeStyle='#fff';ctx.strok"
        "e()}\nfunction render(){ctx.clearRect(0,0,c.width,c.height);const "
        "d=data.filter(x=>x.week==v.value),gs=[...new "
        "Set(d.map(x=>x.group))],ms=[...new Set(d.map(x=>x.metric))];\n"
        "const "
        "bars=d.map(x=>({...x,x:gs.indexOf(x.group)-(gs.length-1)/2,y:ms.i"
        "ndexOf(x.metric)*2-1})).sort((a,b)=>(a.x*Math.sin(angle)+a.y*Math"
        ".cos(angle))-(b.x*Math.sin(angle)+b.y*Math.cos(angle)));\n"
        "bars.forEach(b=>{const "
        "x=b.x,y=b.y,h=b.value,w=.65;poly([[x-w/2,y-w/2,0],[x+w/2,y-w/2,0]"
        ",[x+w/2,y-w/2,h],[x-w/2,y-w/2,h]],b.y<0?'#4e6597':'#d86753');poly"
        "([[x+w/2,y-w/2,0],[x+w/2,y+w/2,0],[x+w/2,y+w/2,h],[x+w/2,y-w/2,h]"
        "],b.y<0?'#334975':'#a74230');poly([[x-w/2,y-w/2,h],[x+w/2,y-w/2,h"
        "],[x+w/2,y+w/2,h],[x-w/2,y+w/2,h]],b.y<0?'#8294bb':'#efa191');let "
        "p=proj(x,y,h);ctx.fillStyle='#18334e';ctx.font='12px Microsoft "
        "YaHei';ctx.textAlign='center';ctx.fillText(h.toFixed(1),p[0],p[1]"
        "-8)});\ngs.forEach((g,i)=>{let "
        "p=proj(i-(gs.length-1)/2,2.4,0);ctx.fillStyle='#18334e';ctx.fillT"
        "ext(g,p[0],p[1]+18)});ctx.textAlign='left';ctx.fillText('蓝色：NMD    "
        "红色：熵    显示指数 0–100',20,30);\nconst "
        "t=document.getElementById('table');t.replaceChildren();let "
        "header=t.insertRow();['组别','指标','原始值','显示指数'].forEach(s=>header.i"
        "nsertCell().textContent=s);d.forEach(b=>{let "
        "r=t.insertRow();[b.group,b.metric,b.raw.toFixed(5),b.value.toFixe"
        "d(3)].forEach(s=>r.insertCell().textContent=s)});}\n"
        "c.onpointerdown=e=>{drag=[e.clientX,e.clientY];c.setPointerCaptur"
        "e(e.pointerId)};c.onpointermove=e=>{if(drag){angle+=(e.clientX-dr"
        "ag[0])*.01;tilt=Math.max(.1,Math.min(1.5,tilt+(e.clientY-drag[1])"
        "*.005));drag=[e.clientX,e.clientY];render()}};c.onpointerup=()=>d"
        "rag=null;v.onchange=render;document.getElementById('reset').oncli"
        "ck=()=>{angle=-.65;tilt=.35;render()};render();</script>"
    )
    page = page.replace("NMD与熵", "NMD与" + html.escape(metric_label)).replace(
        "红色：熵", "红色：" + html.escape(metric_label)
    )
    path.write_text(page.replace("TITLE", html.escape(name)).replace("DATA", data), "utf-8")
