import html, json
import numpy as np
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter
from metric_selection import metrics, LABELS

TITLES = {
    "Energy_distance": "能量距离",
    "N_entropy": "所选方法纳入数",
    "p_entropy": "所选方法指标数",
    "Entropy_Shannon": "Shannon熵",
    "Entropy_Renyi2": "Rényi二阶熵",
    "NMD": "NMD",
    "Experiment": "实验",
    "Group": "组别",
    "Week": "周",
    "N_entropy": "所选方法纳入数",
    "N_NMD": "NMD纳入数",
    "N_present": "记录数",
    "Shannon_numeric_SE": "Shannon数值积分误差",
    "N_paired_baseline": "与基线同时完整人数",
}


def workbook(path, tables, translated=True):
    book = Workbook()
    book.remove(book.active)
    from localization import localize_frame, sheet_title

    for name, df in tables.items():
        ws = book.create_sheet(sheet_title(name))
        df = localize_frame(pd.DataFrame(df))
        headers = (
            [TITLES.get(str(k), str(k)) for k in df.columns] if translated else list(df.columns)
        )
        ws.append(headers or ["记录"])
        for row_index, row in enumerate(df.itertuples(index=False, name=None), 2):
            values = [
                None if pd.isna(v) else v.item() if isinstance(v, np.generic) else v for v in row
            ]
            ws.append(values)
            for column_index, value in enumerate(values, 1):
                if isinstance(value, str):
                    ws.cell(row_index, column_index).data_type = "s"
        for c in ws[1]:
            c.font = Font(bold=True, color="FFFFFF")
            c.fill = PatternFill("solid", fgColor="3D5488")
            c.alignment = Alignment(wrap_text=True)
        ws.freeze_panes = "A2"
        ws.auto_filter.ref = ws.dimensions
        for i, col in enumerate(headers, 1):
            ws.column_dimensions[get_column_letter(i)].width = min(
                42, max(16, len(str(col)) * 1.25)
            )
    book.save(path)


def export_tables(output, results, displays):
    folder = output / "数据表"
    folder.mkdir()
    summary = pd.concat([pd.DataFrame(r["summary"]) for r in results], ignore_index=True)
    overview = []
    contrasts = []
    trends = []
    samples = []
    scales = []
    checks = []
    selected_metrics = metrics(results[0])
    for r in results:
        d = pd.DataFrame(r["summary"])
        detail = output / r["folder"]
        obs = pd.DataFrame(r["observations"])
        samples.append(obs)
        overview.append(
            dict(
                Experiment=r["name"],
                类型="临床" if r["mode"] == "clinical" else "动物",
                组数=d.Group.nunique(),
                时间点数=d.Week.nunique(),
                指标数_熵=int(d.p_entropy.iloc[0]) if "p_entropy" in d else 0,
                数据文件数=len(r["sources"]),
                详情目录=r["folder"],
            )
        )
        for week, b in d.groupby("Week", sort=False):
            groups = list(b.Group)
            reference = r.get("roles", {}).get("comparator") or next(
                (
                    g
                    for g in groups
                    if g == "安慰剂" or "model" in g.lower() or "modle" in g.lower()
                ),
                None,
            )
            if reference not in groups:
                reference = None
            if reference is not None:
                ref = b[b.Group.eq(reference)].iloc[0]
                for row in b.itertuples():
                    if row.Group == reference:
                        continue
                    contrasts.append(
                        dict(
                            Experiment=r["name"],
                            Week=week,
                            Group=row.Group,
                            参照组=reference,
                            指定治疗组=row.Group == r.get("roles", {}).get("treatment"),
                            **{
                                m + "_difference": getattr(row, m) - ref[m]
                                for m in selected_metrics
                            },
                            说明="组估计值减参照组；描述性差值，无配对推断或显著性检验",
                        )
                    )
        if r["mode"] == "clinical":
            for group, b in d.groupby("Group", sort=False):
                baseline = b[b.Week.eq(0)].iloc[0]
                for row in b.itertuples():
                    trends.append(
                        dict(
                            Experiment=r["name"],
                            Group=group,
                            Week=row.Week,
                            **{
                                m + "_change_from_0W": getattr(row, m) - baseline[m]
                                for m in selected_metrics
                            },
                            说明="组估计相对0W的变化；不同访视纳入者可能不同",
                        )
                    )
        for filename, destination in [("固定尺度.csv", scales), ("缺失与记录核查.csv", checks)]:
            p = detail / filename
            try:
                df = pd.read_csv(p)
            except pd.errors.EmptyDataError:
                continue
            if len(df):
                df.insert(0, "Experiment", r["name"])
                destination.append(df)
    workbook(folder / "01_实验总览.xlsx", {"实验总览": overview})
    workbook(folder / "02_组时间点汇总.xlsx", {"组时间点汇总": summary})
    workbook(folder / "03_治疗比较.xlsx", {"描述性比较": contrasts})
    workbook(folder / "04_时间趋势.xlsx", {"相对0W变化": trends})
    workbook(folder / "08_Display_Normalized_Index.xlsx", {"原值与显示值": displays})
    workbook(
        folder / "09_样本级结果.xlsx", {"样本NMD及纳入": pd.concat(samples, ignore_index=True)}
    )
    workbook(
        folder / "10_尺度与数据核查.xlsx",
        {
            "固定尺度": pd.concat(scales, ignore_index=True) if scales else [],
            "数据核查": pd.concat(checks, ignore_index=True) if checks else [],
        },
    )
    from analysis_outputs import guide_for

    readme = guide_for(results[0].get("selected_method", "Shannon"))
    (output / "结果说明.md").write_text(readme, "utf-8")


def write_viewer(output, results, index, display):
    esc = html.escape
    blocks = []
    for r in results:
        cards = []
        for item in index:
            if item["Experiment"] != r["name"]:
                continue
            media = (
                f'<a href="{esc(item["File"], quote=True)}" target="_blank">打开三维交互图</a>'
                if item["File"].endswith(".html")
                else f'<img loading="lazy" src="{esc(item["File"], quote=True)}" alt="{esc(item["Title"], quote=True)}">'
            )
            kind = "辅助图" if item.get("Category") == "auxiliary" else "图"
            comparison = {
                "raw": "当周/当次原值之差",
                "change": "两组相对0W变化之差",
                "none": "本图无自动P值",
            }.get(item.get("P_contrast"), "")
            explanation = (
                "<p>图表方案："
                f"{esc(item.get('Figure_policy', '常规'))}"
                "；"
                f"{esc(comparison)}"
                '。<a href="数据表/16_图表元素与来源.xlsx">查看元素、散点与P值来源</a></p>'
            )
            cards.append(
                (
                    '<article data-method="'
                    f"{item['Method']}"
                    '"><h3>'
                    f"{kind}"
                    f"{item['Number']:02d}"
                    " · "
                    f"{esc(item['Title'])}"
                    "</h3>"
                    f"{explanation}"
                    '<a href="'
                    f"{esc(item['Data'], quote=True)}"
                    '">图形数据</a>'
                    f"{media}"
                    "</article>"
                )
            )
        blocks.append(
            f'<section data-experiment="{esc(r["name"], quote=True)}"><h2>{esc(r["name"])}</h2>'
            + "".join(cards)
            + "</section>"
        )
    options = "".join("<option>" + esc(r["name"]) + "</option>" for r in results)
    page = (
        '<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta '
        'name="viewport" '
        'content="width=device-width,initial-scale=1"><title>BioEntropy '
        "V1.0.2</title>\n"
        "<style>body{margin:0;background:#edf1f7;color:#18334e;font-family"
        ":'Microsoft "
        "YaHei',sans-serif}header{padding:28px;background:#18334e;color:wh"
        "ite;border-bottom:4px solid #00a087}header "
        "h1{margin:0}main{max-width:1150px;margin:auto;padding:22px}.contr"
        "ols{position:sticky;top:0;background:white;padding:15px;display:f"
        "lex;gap:20px;flex-wrap:wrap;z-index:2;box-shadow:0 2px 6px "
        "#0002}select{padding:7px;font-size:15px}article{background:white;"
        "border:1px solid #d9e0e8;padding:16px;margin:16px "
        "0}img{width:100%;height:auto}p{line-height:1.8}a{color:#3d5488}h3"
        "{margin:0 0 8px}[hidden]{display:none!important}</style>\n"
        "<header><h1>BioEntropy</h1><p>V1.0.2 · "
        "本次仅输出NMD与所选方法</p></header><main>\n<div "
        'class="controls"><label>本次方法 <select '
        'id="method">METHOD_OPTION</select></label><label>实验 <select '
        'id="experiment"><option '
        'value="all">全部实验</option>OPTIONS</select></label></div>\n<p><a '
        'href="数据表/02_组时间点汇总.xlsx">组时间点汇总</a> · <a '
        'href="数据表/08_Display_Normalized_Index.xlsx">原值与显示值</a> · <a '
        'href="图形数据_GraphPad/图表-数据对照表.csv">GraphPad数据索引</a> · <a '
        'href="结果说明.md">结果说明</a> · <a href="散点与显著性说明.html">散点怎么算、怎么看</a> · '
        '<a href="数据表/15_散点与显著性.xlsx">散点及显著性结果表</a></p>\nCONTENT\n'
        "<p>原始数值、显示尺度和纳入记录分别保存。仅生成本次所选方法；如需其他方法，请在软件中重新选择并计算到新目录。</p></mai"
        "n>\n<script>const "
        "m=document.getElementById('method'),e=document.getElementById('ex"
        "periment');function "
        "refresh(){document.querySelectorAll('article').forEach(a=>a.hidde"
        "n=a.dataset.method!=='NMD'&&m.value!=='both'&&a.dataset.method!=="
        "m.value);document.querySelectorAll('section').forEach(a=>a.hidden"
        "=e.value!=='all'&&a.dataset.experiment!==e.value);}const p=new "
        "URLSearchParams(location.hash.slice(1));m.value='DEFAULT';m.oncha"
        "nge=e.onchange=refresh;refresh();</script></html>"
    )
    (output / "查看结果.html").write_text(
        page.replace("OPTIONS", options)
        .replace("CONTENT", "\n".join(blocks))
        .replace("'DEFAULT'", json.dumps(display))
        .replace(
            "METHOD_OPTION",
            '<option value="'
            + html.escape(display)
            + '">'
            + html.escape(
                LABELS.get(
                    {
                        "Shannon": "Entropy_Shannon",
                        "Renyi2": "Entropy_Renyi2",
                        "Energy": "Energy_distance",
                    }.get(display),
                    "仅NMD",
                )
            )
            + "</option>",
        ),
        "utf-8",
    )
