from pathlib import Path
import json, datetime, platform
import pandas as pd
from analysis_engine import clinical, animal, profiles, sha, csv, write_json, VERSION
from clinical_inputs import InputBundle, WEEKS, InputSelectionError
from clinical_patient_ids import load_patient_data
from animal_inputs import read_animal
from experiment_config import discover, save_project
from output_paths import desktop_directory
from metric_selection import validate_method


def known_clinical(e):
    return e.name in profiles() and set(e.weeks) <= set(WEEKS)


def bundle(e):
    return InputBundle(
        e.name,
        Path(e.files[0]).parent,
        Path(e.normal_files[0]) if e.normal_files else None,
        visit_files=tuple(Path(p) for p in e.files),
    )


def preflight(experiments):
    from generic_longitudinal import prepare
    from supplementary import range_bounds
    from reference_data import read_reference
    from experiment_config import EXTERNAL

    active = [e for e in experiments if e.enabled]
    if not active:
        raise InputSelectionError("请至少勾选一个实验。")
    if len({e.name for e in active}) != len(active):
        raise InputSelectionError("重复实验名称。")
    notes = []
    for e in active:
        e.validate()
        if e.kind == "animal":
            d = read_animal(e.files[0], e.reference)
            if e.reference == EXTERNAL:
                read_reference(e.normal_files, d["features"])
        elif known_clinical(e):
            load_patient_data(e.name, profiles()[e.name], bundle(e), e.record())
        else:
            prepare(e)
        if e.range_files:
            range_bounds(e.range_files)
        notes.append(
            e.name
            + "："
            + ("单时间点分组数据" if e.kind == "animal" else "多时间点/患者数据")
            + "，"
            + str(len(e.groups))
            + "组，"
            + str(len(e.weeks))
            + "个时间点"
        )
    return notes


def _calculate_project(experiments, output, display="Shannon", progress=print, permutations=9999):
    from plots import create_report
    from compact_outputs import joint_figures
    from table_exports import write_viewer
    from generic_longitudinal import analyze

    validate_method(display)
    active = [e for e in experiments if e.enabled]
    for note in preflight(active):
        progress("输入核查：" + note)
    output = Path(output).resolve()
    if output.exists() and any(output.iterdir()):
        raise InputSelectionError("请选择空的结果目录。")
    paths = list(
        dict.fromkeys(p for e in active for p in [*e.files, *e.normal_files, *e.range_files])
    )
    for p in paths:
        parent = Path(p).resolve().parent
        if output == parent or parent in output.parents:
            if not (parent == desktop_directory().resolve() and output.parent == parent):
                raise InputSelectionError("结果目录应在原始数据目录之外，避免混入输入。")
    sources = {p: sha(p) for p in paths}
    output.mkdir(parents=True, exist_ok=True)
    detail = output / "计算明细"
    detail.mkdir()
    write_json(output / "运行状态.json", dict(status="running", version=VERSION))
    save_project(output / "本次分析设置.json", experiments, display, str(output), permutations)
    results = []
    try:
        for e in active:
            progress(e.name + "：按所选方法计算，并按设定顺序排列图表")
            options = dict(e.record(), _method=display)
            if display == "Energy":
                from energy_analysis import animal_energy, clinical_energy, generic_energy

                if e.kind == "animal":
                    r = animal_energy(e.files[0], detail / e.name, progress, options)
                elif known_clinical(e):
                    r = clinical_energy(bundle(e), detail / e.name, progress, options)
                else:
                    r = generic_energy(e, detail / e.name, progress)
            elif e.kind == "animal":
                r = animal(e.files[0], detail / e.name, progress, options)
            elif known_clinical(e):
                r = clinical(bundle(e), detail / e.name, progress, options)
            else:
                r = analyze(e, detail / e.name, progress, display)
            r["selected_method"] = display
            if display == "none":
                for row in r["summary"]:
                    for key in ("N_entropy", "p_entropy", "Bandwidth", "N_paired_baseline"):
                        row.pop(key, None)
            r["roles"] = dict(treatment=e.treatment, comparator=e.comparator, reference=e.reference)
            from analysis_outputs import enhance_result

            enhance_result(r, detail / e.name, progress, permutations)
            r["_source_observations"] = [dict(x) for x in r["observations"]]
            rank = {g: i for i, g in enumerate(e.order)}
            r["summary"].sort(key=lambda row: (rank[row["Group"]], row["Week"]))
            r["observations"].sort(key=lambda row: (rank[row["Group"]], row.get("Week", 0)))
            df = pd.DataFrame(r["summary"])
            csv(df, detail / e.name / "组结果.csv")
            from openpyxl import load_workbook

            book = load_workbook(detail / e.name / "结果表.xlsx")
            ws = book["组结果"]
            ws.delete_rows(2, ws.max_row)
            for row in df.itertuples(index=False, name=None):
                ws.append([None if pd.isna(x) else x for x in row])
                for c in ws[ws.max_row]:
                    if isinstance(c.value, str):
                        c.data_type = "s"
            book.save(detail / e.name / "结果表.xlsx")
            book.close()
            r["folder"] = "计算明细/" + r["folder"]
            cfgpath = output / r["folder"] / "计算配置与来源.json"
            cfg = json.loads(cfgpath.read_text("utf-8"))
            cfg["user_configuration"] = e.record()
            cfg["selected_method"] = display
            write_json(cfgpath, cfg)
            results.append(r)
        csv(
            pd.concat([pd.DataFrame(r["summary"]) for r in results], ignore_index=True),
            output / "全部实验_原始结果.csv",
        )
        progress("生成主结果图及绘图数据")
        create_report(output, results, display)
        joint_figures(output, results)
        from export_catalog import write_traceability

        (output / "数据表").mkdir(exist_ok=True)
        write_traceability(output, results)
        for path, digest in sources.items():
            if sha(path) != digest:
                raise InputSelectionError("计算期间原始文件发生变化：" + path)
        write_json(
            output / "输入文件核验.json",
            [dict(path=p, sha256=v, unchanged=True) for p, v in sources.items()],
        )
        write_json(
            output / "运行状态.json",
            dict(
                status="completed",
                version=VERSION,
                mode="auto",
                experiments=[e.name for e in active],
                original_files_unchanged=True,
                python=platform.python_version(),
                completed=datetime.datetime.now().isoformat(),
            ),
        )
        return results
    except Exception as exc:
        write_json(output / "运行状态.json", dict(status="failed", version=VERSION, error=str(exc)))
        raise


def run_auto(inputs, output, display="Shannon", progress=print, permutations=9999):
    return run_project(discover(inputs), output, display, progress, permutations)


def run_project(
    experiments, output, display="Shannon", progress=print, permutations=9999, language="zh"
):
    from compact_outputs import run_compact

    return run_compact(
        _calculate_project, experiments, output, display, progress, permutations, language
    )
