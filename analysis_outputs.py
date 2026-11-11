from pathlib import Path
import html
import numpy as np
import pandas as pd
from sample_information import add_scores
from comparisons import compare_policies, all_comparisons


def enhance_result(result, detail, progress=print, permutations=9999):
    from analysis_engine import csv, write_json
    from table_exports import workbook

    method = result.get("selected_method", "Shannon")
    result["score_checks"] = []
    if method == "Shannon":
        add_scores(result, progress)
    elif method == "Energy":
        obs = pd.DataFrame(result["observations"])
        if "Week" not in obs:
            obs["Week"] = 0
        for group in result["summary"]:
            part = obs[obs.Group.eq(group["Group"]) & obs.Week.eq(group["Week"])]
            values = part.Energy_contribution_score.dropna()
            result["score_checks"].append(
                dict(
                    Experiment=result["name"],
                    Group=group["Group"],
                    Week=group["Week"],
                    N_entropy=group["N_entropy"],
                    N_scores=len(values),
                    Original_Energy=group["Energy_distance"],
                    Mean_contribution_score=float(values.mean()) if len(values) else np.nan,
                    Numerical_difference=float(values.mean() - group["Energy_distance"])
                    if len(values)
                    else np.nan,
                )
            )
    compare_policies(result, permutations, progress)
    detail = Path(detail)
    csv(result["observations"], detail / "逐样本NMD及纳入.csv")
    if method in ("Shannon", "Energy"):
        csv(result["observations"], detail / (method + "样本分数.csv"))
    if method in ("Shannon", "Energy"):
        csv(result["score_checks"], detail / "散点与组估计数值核对.csv")
    csv(all_comparisons(result), detail / "探索性置换检验.csv")
    csv(result["comparisons"], detail / "探索性置换检验_常规方案.csv")
    if result.get("comparisons_visit"):
        csv(result["comparisons_visit"], detail / "探索性置换检验_分周原值方案.csv")
        csv(
            result["permutation_numeric_checks_visit"], detail / "置换积分精度核对_分周原值方案.csv"
        )
    csv(result.get("permutation_numeric_checks", []), detail / "置换积分精度核对.csv")
    write_json(detail / "散点与显著性设置.json", result["inference_settings"])
    workbook(
        detail / "散点与显著性核查.xlsx",
        {
            "样本分数": result["observations"],
            "分数均值核对": result["score_checks"],
            "探索性置换检验": all_comparisons(result),
            "置换积分核对": result.get("permutation_numeric_checks", []),
            "原值联合图积分核对": result.get("permutation_numeric_checks_visit", []),
        },
        translated=False,
    )
    from openpyxl import load_workbook

    path = detail / "结果表.xlsx"
    book = load_workbook(path)
    ws = book["逐样本NMD及纳入"]
    ws.delete_rows(1, ws.max_row)
    frame = pd.DataFrame(result["observations"])
    ws.append(list(frame.columns))
    for row_index, row in enumerate(frame.itertuples(index=False, name=None), 2):
        values = [None if pd.isna(v) else v.item() if isinstance(v, np.generic) else v for v in row]
        ws.append(values)
        for column_index, value in enumerate(values, 1):
            if isinstance(value, str):
                ws.cell(row_index, column_index).data_type = "s"
    ws.freeze_panes = "A2"
    ws.auto_filter.ref = ws.dimensions
    book.save(path)
    book.close()


def annotation(result, metric, week, normalized, source="常规"):
    if source == "常规" and normalized and (result["mode"] == "clinical") and (week == 0):
        return "0W对齐"
    if result.get("inference_settings", {}).get("permutations") == 0:
        return "未检验"
    contrast = (
        "change" if source == "常规" and normalized and (result["mode"] == "clinical") else "raw"
    )
    rows = (
        result.get("comparisons_visit", [])
        if source == "分周原值"
        else result.get("comparisons", [])
    )
    match = [
        r
        for r in rows
        if r["Metric"] == metric and r["Week"] == week and (r["Contrast"] == contrast)
    ]
    return match[0]["Annotation"] if match else "NA"


def annotate_metric(ax, result, metric, normalized):
    roles = result.get("roles", {})
    if not result.get("comparisons"):
        return
    if result["mode"] == "clinical":
        for w in sorted({r["Week"] for r in result["summary"]}):
            ax.text(
                w,
                1.025,
                annotation(result, metric, w, normalized),
                transform=ax.get_xaxis_transform(),
                ha="center",
                va="bottom",
                fontsize=8,
                color="#29374c",
            )
    else:
        ax.text(
            0.5,
            1.025,
            str(roles.get("treatment"))
            + " vs "
            + str(roles.get("comparator"))
            + ": "
            + annotation(result, metric, 0, False),
            transform=ax.transAxes,
            ha="center",
            va="bottom",
            fontsize=9,
            color="#29374c",
        )


def guide_for(method):
    from output_help import guide_for as explain

    return explain(method)


GUIDE = guide_for("Shannon")


def write_guide(output, results):
    from analysis_engine import csv, write_json
    from table_exports import workbook

    output = Path(output)
    scores = pd.concat([pd.DataFrame(r["observations"]) for r in results], ignore_index=True)
    checks = [x for r in results for x in r.get("score_checks", [])]
    tests = [x for r in results for x in all_comparisons(r)]
    config = [dict(Experiment=r["name"], **r.get("inference_settings", {})) for r in results]
    csv(tests, output / "探索性显著性结果.csv")
    workbook(
        output / "数据表/15_散点与显著性.xlsx",
        {"样本分数": scores, "分数均值核对": checks, "探索性置换检验": tests},
        translated=False,
    )
    write_json(output / "散点与显著性设置.json", config)
    guide = guide_for(results[0].get("selected_method", "Shannon"))
    extra = "\n\n## 本次运行\n\n" + "\n".join(
        "- "
        + r["name"]
        + "："
        + r.get("inference_settings", {}).get("status", "未运行")
        + "；置换次数 "
        + str(r.get("inference_settings", {}).get("permutations", 0))
        for r in results
    )
    (output / "散点与显著性说明.md").write_text(guide + extra, "utf-8")
    body = []
    for paragraph in (guide + extra).split("\n\n"):
        text = html.escape(paragraph).replace("\n", "<br>")
        tag = "h1" if paragraph.startswith("# ") else "h2" if paragraph.startswith("## ") else "p"
        if tag in ("h1", "h2"):
            text = text.lstrip("# ")
        body.append("<" + tag + ">" + text + "</" + tag + ">")
    page = (
        (
            '<html lang="zh-CN"><meta '
            'charset="utf-8"><title>散点与显著性说明</title><style>body{font:16px/1.85 '
            '"Microsoft YaHei",sans-serif;max-width:1000px;margin:35px '
            "auto;padding:0 "
            "25px;color:#223344}h2{margin-top:2em;border-bottom:1px solid "
            "#ccd}p{overflow-wrap:anywhere}</style>"
        )
        + "".join(body)
        + "</html>"
    )
    (output / "散点与显著性说明.html").write_text(page, "utf-8")
