from pathlib import Path
from tempfile import TemporaryDirectory, gettempdir
import datetime
import html
import json
import re
import shutil
import traceback
import numpy as np
import pandas as pd
from analysis_engine import VERSION, sha, write_json
from clinical_inputs import InputSelectionError
from comparisons import all_comparisons
from experiment_config import save_project
from metric_selection import metrics, score_field, LABELS
from output_paths import desktop_directory
from table_exports import workbook
from localization import translate, presentation, current_language, localize_html, install_rendering


def joint_figures(output, results):
    from weekly_state import draw

    output = Path(output)
    displays = pd.read_csv(output / "显示值与原始值.csv")
    index = json.loads((output / "图表索引.json").read_text("utf-8"))
    for result in results:
        if result["mode"] != "clinical":
            continue
        for metric in metrics(result)[1:]:
            index.append(draw(result, metric, displays, output))
            if score_field(metric) is not None:
                index.append(draw(result, metric, displays, output, focus=True))
    write_json(output / "图表索引.json", index)


def destination_path(path, experiments=()):
    target = Path(path).resolve()
    if target.exists() and (not target.is_dir() or any(target.iterdir())):
        raise InputSelectionError("请选择空的结果目录。")
    for experiment in experiments:
        if not experiment.enabled:
            continue
        for filename in [*experiment.files, *experiment.normal_files, *experiment.range_files]:
            parent = Path(filename).resolve().parent
            if target == parent or parent in target.parents:
                if not (parent == desktop_directory().resolve() and target.parent == parent):
                    raise InputSelectionError("结果目录应在原始数据目录之外，避免混入输入。")
    return target


def safe_name(value, used):
    stem = re.sub(r'[<>:"/\\|?*\x00-\x1f]', "_", str(value)).rstrip(" .") or "实验"
    if stem.split(".")[0].upper() in {
        "CON",
        "PRN",
        "AUX",
        "NUL",
        *(f"COM{i}" for i in range(1, 10)),
        *(f"LPT{i}" for i in range(1, 10)),
    }:
        stem = "_" + stem
    candidate = stem
    number = 2
    while candidate.casefold() in used:
        candidate = f"{stem}-{number}"
        number += 1
    used.add(candidate.casefold())
    return candidate


def read_csv(path):
    try:
        return pd.read_csv(
            path, dtype={"USUBJID": str, "Sample_ID": str, "Subject": str, "Trial": str}
        )
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def settings_rows(experiment, section, settings):
    rows = []
    for key, value in settings.items():
        if key == "entropy_uses_NMD":
            continue
        if isinstance(value, (list, dict)):
            value = json.dumps(value, ensure_ascii=False, allow_nan=True)
        text = str(value)
        if len(text) > 30000:
            rows.extend(
                dict(
                    Experiment=experiment,
                    Section=section,
                    Setting=f"{key}[{i // 30000 + 1}]",
                    Value=text[i : i + 30000],
                )
                for i in range(0, len(text), 30000)
            )
        else:
            rows.append(dict(Experiment=experiment, Section=section, Setting=key, Value=value))
    return rows


def write_viewer(target, index, method):
    esc = html.escape
    names = list(dict.fromkeys(item["Experiment"] for item in index))
    options = "".join(f"<option>{esc(name)}</option>" for name in names)
    sections = []
    for name in names:
        cards = []
        for item in index:
            if item["Experiment"] != name:
                continue
            title = (f"图{item['Number']:02d} · " if item.get("Number") is not None else "") + item[
                "Title"
            ]
            contrast = {"raw": "原始数值比较", "change": "相对0W变化比较", "none": ""}.get(
                item.get("P_contrast"), ""
            )
            data = item.get("Data_file", "结果表.xlsx")
            cards.append(
                (
                    "<article><h3>"
                    f"{esc(title)}"
                    "</h3><p>"
                    f"{esc(contrast)}"
                    ' <a href="'
                    f"{esc(data, quote=True)}"
                    '">绘图数据</a> · <a href="'
                    f"{esc(item['File'], quote=True)}"
                    '" target="_blank">打开图片</a></p><a href="'
                    f"{esc(item['File'], quote=True)}"
                    '" target="_blank"><img loading="lazy" src="'
                    f"{esc(item['File'], quote=True)}"
                    '" alt="'
                    f"{esc(title, quote=True)}"
                    '"></a></article>'
                )
            )
        sections.append(
            f'<section data-experiment="{esc(name, quote=True)}"><h2>{esc(name)}</h2>{"".join(cards)}</section>'
        )
    page = (
        '<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta '
        'name="viewport" '
        'content="width=device-width,initial-scale=1"><title>BioEntropy '
        "V1.0.2</title><style>body{margin:0;background:#f2f4f7;color:#2136"
        '4c;font:15px/1.65 "Microsoft '
        'YaHei",sans-serif}header{background:#18334e;color:white;padding:2'
        "2px 5vw}h1{margin:0;font-size:28px}header p{margin:4px "
        "0}main{max-width:1250px;margin:auto;padding:20px}.nav{position:st"
        "icky;top:0;background:white;padding:14px;border:1px solid "
        "#d9e0e8;z-index:1}select{padding:5px;font-size:15px}a{color:#3d54"
        "88}article{background:white;border:1px solid "
        "#d9e0e8;padding:16px;margin:15px "
        "0}h2{margin-top:25px}h3{margin:0}article p{margin:6px 0 "
        "12px}img{width:100%;height:auto}[hidden]{display:none!important}<"
        "/style><header><h1>BioEntropy</h1><p>V1.0.2 · "
        'METHOD</p></header><main><div class="nav"><a '
        'href="结果表.xlsx">打开结果表</a>\u3000实验：<select id="experiment"><option '
        'value="all">全部实验</option>OPTIONS</select></div>CONTENT</main><scr'
        "ipt>const "
        "select=document.getElementById('experiment');select.onchange=()=>"
        "document.querySelectorAll('section').forEach(s=>s.hidden=select.v"
        "alue!=='all'&&s.dataset.experiment!==select.value);</script></htm"
        "l>"
    )
    (target / "查看结果.html").write_text(
        localize_html(
            page.replace("METHOD", esc(method))
            .replace("OPTIONS", options)
            .replace("CONTENT", "".join(sections))
        ),
        "utf-8",
    )


def publish(stage, target, results, method, permutations):
    stage, target = Path(stage), Path(target)
    index = json.loads((stage / "图表索引.json").read_text("utf-8"))
    displays = read_csv(stage / "显示值与原始值.csv")
    chart_tests = read_csv(stage / "图形数据_GraphPad" / "每图显著性.csv")
    input_checks = json.loads((stage / "输入文件核验.json").read_text("utf-8"))
    used = set()
    names = {result["name"]: safe_name(result["name"], used) for result in results}
    settings, counts = [], []
    settings.extend(
        settings_rows(
            "",
            "运行",
            dict(
                Version=VERSION,
                Language=current_language(),
                Method=method,
                Permutations=permutations,
                Status="completed",
                Original_files_unchanged=True,
            ),
        )
    )
    for result in results:
        name = result["name"]
        folder = target / "图表" / names[name]
        folder.mkdir(parents=True, exist_ok=True)
        tables = {}
        for item in index:
            if item["Experiment"] != name:
                continue
            original = stage / item["File"]
            picture = folder / original.name
            shutil.copy2(original, picture)
            item["File"] = picture.relative_to(target).as_posix()
            number = item["Number"]
            data = read_csv(stage / item["Data"])
            if "Experiment" in data:
                data = data[data.Experiment.astype(str).eq(name)]
            tables[f"图{number:02d}_组值"] = data
            datafile = stage / item["Data"]
            if number in (1, 2, 3, 4):
                for suffix, label in [("_points.csv", "散点"), ("_prism.csv", "Prism")]:
                    path = datafile.with_name(datafile.name.replace("_bars.csv", suffix))
                    if path.is_file():
                        tables[f"图{number:02d}_{label}"] = read_csv(path)
            elif number in (10, 11):
                for suffix, label in [("_points.csv", "散点"), ("_counts.csv", "窗口计数")]:
                    path = datafile.with_name(datafile.name.replace("_groups.csv", suffix))
                    if path.is_file():
                        tables[f"图{number:02d}_{label}"] = read_csv(path)
            tests = (
                chart_tests[
                    chart_tests.Experiment.astype(str).eq(name) & chart_tests.Figure.eq(number)
                ].copy()
                if not chart_tests.empty
                else pd.DataFrame()
            )
            if not tests.empty:
                tests["Figure_file"] = item["File"]
                tests["P_source"] = "结果表.xlsx#显著性"
                tables[f"图{number:02d}_P值"] = tests
            item["Data_file"] = f"绘图数据/{names[name]}.xlsx"
            item["Data"] = item["Data_file"] + f"#图{number:02d}_组值"
            item["P_data"] = "结果表.xlsx#显著性"
        target.joinpath("绘图数据").mkdir(exist_ok=True)
        workbook(target / "绘图数据" / (names[name] + ".xlsx"), tables, translated=False)
        config = json.loads((stage / result["folder"] / "计算配置与来源.json").read_text("utf-8"))
        settings.extend(settings_rows(name, "计算设置", config))
        settings.extend(settings_rows(name, "统计设置", result.get("inference_settings", {})))
        observations = pd.DataFrame(result["observations"])
        if "Week" not in observations:
            observations["Week"] = 0
        selected = metrics(result)[1:]
        field = score_field(selected[0]) if selected else None
        for (group, week), part in observations.groupby(["Group", "Week"], sort=False):
            nmd = np.isfinite(part.NMD.to_numpy(float))
            row = dict(
                Experiment=name,
                Group=group,
                Week=int(week),
                N_present=len(part),
                N_NMD=int(nmd.sum()),
            )
            if selected:
                usable = (
                    part.Included_entropy.to_numpy(bool)
                    if "Included_entropy" in part
                    else np.isfinite(part[field].to_numpy(float))
                    if field in part
                    else np.zeros(len(part), bool)
                )
                score = (
                    np.isfinite(part[field].to_numpy(float))
                    if field and field in part
                    else np.zeros(len(part), bool)
                )
                row.update(
                    Selected_metric=selected[0],
                    N_selected=int(usable.sum()),
                    N_scores=int(score.sum()),
                    N_common_analysis=int((nmd & usable).sum()),
                    N_joint_scatter=int((nmd & score).sum()),
                )
            counts.append(row)
    if not chart_tests.empty:
        file_lookup = {(item["Experiment"], item["Number"]): item["File"] for item in index}
        chart_tests["Figure_file"] = [
            file_lookup[(str(row.Experiment), row.Figure)] for row in chart_tests.itertuples()
        ]
        chart_tests["P_source"] = "结果表.xlsx#显著性"
    common_audit = []
    reasons = {
        "Included": "纳入",
        "Incomplete entropy panel": "所选熵面板存在缺失指标",
        "NMD unavailable": "NMD不可计算",
        "Fewer than two common records in group/visit": "本组本访视共同有效记录不足2例",
    }
    for result in results:
        for row in result["observations"]:
            if "Included_common" in row:
                common_audit.append(
                    {
                        "实验": result["name"],
                        "组别": row["Group"],
                        "周": row.get("Week", 0),
                        "Trial": row.get("Trial", ""),
                        "个体编号": row.get("USUBJID", row.get("Sample_ID", "")),
                        "共同纳入": bool(row["Included_common"]),
                        "原因": reasons[row["Selection_reason"]],
                        "缺失指标": row.get("Missing_entropy_features", ""),
                        "来源文件": row.get("Source_file", ""),
                        "Excel行号": row.get("Excel_row"),
                    }
                )
    workbook(
        target / "结果表.xlsx",
        {
            "组结果": pd.concat(
                [pd.DataFrame(result["summary"]) for result in results], ignore_index=True
            ),
            "样本结果": pd.concat(
                [pd.DataFrame(result["observations"]) for result in results], ignore_index=True
            ),
            "显著性": [row for result in results for row in all_comparisons(result)],
            "原值与显示值": displays,
            "纳入人数": counts,
            "各图P值": chart_tests,
            "图表索引": index,
            "运行设置": settings,
            "输入文件": input_checks,
            **({"共同样本清单": common_audit} if common_audit else {}),
        },
        translated=False,
    )
    title = {
        "Shannon": "NMD + Shannon熵",
        "Renyi2": "NMD + Rényi二阶熵",
        "Energy": "NMD + 能量距离",
        "none": "NMD",
    }[method]
    write_viewer(target, index, title)


def logger(target, progress):
    def write(message):
        with (target / "运行日志.txt").open("a", encoding="utf-8") as handle:
            handle.write(translate(str(message)) + "\n")
        progress(message)

    return write


def _run_compact(calculate, experiments, output, display, progress, permutations):
    target = destination_path(output, experiments)
    target.mkdir(parents=True, exist_ok=True)
    log = logger(target, progress)
    save_project(target / "本次分析设置.json", experiments, display, str(target), permutations)
    project = json.loads((target / "本次分析设置.json").read_text("utf-8"))
    project["language"] = current_language()
    write_json(target / "本次分析设置.json", project)
    log(VERSION + " · " + datetime.datetime.now().isoformat(timespec="seconds"))
    try:
        with TemporaryDirectory(prefix="BioEntropy_V1_0_0_") as temporary:
            staging = Path(temporary).resolve()
            assert staging.parent == Path(gettempdir()).resolve() and staging.name.startswith(
                "BioEntropy_V1_0_0_"
            )
            stage = staging / "results"
            results = calculate(experiments, stage, display, log, permutations)
            log("整理结果图、绘图数据及结果表")
            publish(stage, target, results, display, permutations)
        log("全部完成：" + str(target / "查看结果.html"))
        return target
    except Exception:
        detail = traceback.format_exc()
        (target / "错误信息.txt").write_text(detail, "utf-8")
        log("运行未完成；详细原因见错误信息.txt。")
        raise


def _run_gvalues(source, destination, progress=print):
    import gvalue_recompute as engine

    root = Path(engine.find_gvalue_dir(str(source))).resolve()
    target = Path(destination).resolve()
    if target == root or root in target.parents:
        raise InputSelectionError("G值结果目录应位于原始G值数据目录之外。")
    target = destination_path(target)
    target.mkdir(parents=True, exist_ok=True)
    log = logger(target, progress)
    sources = [
        dict(path=str(p), sha256=sha(p))
        for p in root.rglob("*.xlsx")
        if not p.name.startswith("~$")
    ]
    log(VERSION + " · Hedges’ g · " + datetime.datetime.now().isoformat(timespec="seconds"))
    try:
        with TemporaryDirectory(prefix="BioEntropy_V1_0_0_G_") as temporary:
            staging = Path(temporary).resolve()
            assert staging.parent == Path(gettempdir()).resolve() and staging.name.startswith(
                "BioEntropy_V1_0_0_G_"
            )
            log("计算G值并生成森林图")
            result = engine.analyze_gvalues(str(root), str(staging))
            for item in sources:
                if sha(item["path"]) != item["sha256"]:
                    raise InputSelectionError("G值输入文件在运行期间发生变化。")
            folder = Path(result["output_dir"])
            titles = {
                "gvalue_part1_individual_endpoint_G": "单指标G值",
                "gvalue_part2_individual_endpoints": "模型内指标",
                "gvalue_part2_within_model_pooled": "模型内合并",
                "gvalue_part3_per_model": "跨模型各模型",
                "gvalue_part3_pooled": "跨模型合并",
                "gvalue_part4_distances": "NMD汇总",
                "gvalue_part4_g": "NMD_G值",
                "gvalue_part4_individual_NMD": "样本NMD",
                "gvalue_part4_Normal_reference": "正常参照",
            }
            tables = {
                titles.get(p.stem, p.stem[:31]): read_csv(p) for p in sorted(folder.glob("*.csv"))
            }
            for path in sorted((folder / "QC").glob("*.csv")):
                tables[path.stem[:31]] = read_csv(path)
            tables["输入文件"] = pd.DataFrame(sources)
            tables["运行设置"] = pd.DataFrame(
                [
                    dict(
                        Version=VERSION,
                        Method="Hedges_g",
                        NMD_sample_selection="All animals with a valid NMD",
                        NMD_reference_uncertainty="Conditional on the fitted Normal reference",
                        Status="completed",
                        Original_files_unchanged=True,
                    )
                ]
            )
            workbook(target / "结果表.xlsx", tables, translated=False)
            figures = target / "图表"
            figures.mkdir(exist_ok=True)
            index = []
            for path in sorted((folder / "figures").glob("*.png")):
                shutil.copy2(path, figures / path.name)
                index.append(
                    dict(
                        Experiment="G值",
                        Title=path.stem,
                        File="图表/" + path.name,
                        Data_file="结果表.xlsx",
                    )
                )
            write_viewer(target, index, "Hedges’ g 与森林图")
        log("全部完成：" + str(target / "查看结果.html"))
        return target
    except Exception:
        (target / "错误信息.txt").write_text(traceback.format_exc(), "utf-8")
        log("运行未完成；详细原因见错误信息.txt。")
        raise


def run_compact(calculate, experiments, output, display, progress, permutations, language="zh"):
    names = [value for experiment in experiments for value in [experiment.name, *experiment.groups]]
    install_rendering()
    with presentation(language, names):
        return _run_compact(calculate, experiments, output, display, progress, permutations)


def run_gvalues(source, destination, progress=print, language="zh"):
    install_rendering()
    with presentation(language):
        return _run_gvalues(source, destination, progress)
