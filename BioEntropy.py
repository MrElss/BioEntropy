import os

for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[name] = "1"
import sys, argparse, traceback
from pathlib import Path
from analysis_engine import VERSION, write_json, sha
from clinical_inputs import InputSelectionError

MODES = {
    "临床实验（带患者编号的分周表）": "clinical",
    "动物实验（分组工作簿）": "animal",
    "G值计算与森林图": "gvalues",
}
SCOPES = {
    "自动识别所选目录": "auto",
    "105、106及合并数据": "all",
    "仅105": "105",
    "仅106": "106",
    "仅105＋106合并": "105and106",
}


def app_root():
    return (
        Path(sys.executable).parent
        if getattr(sys, "frozen", False)
        else Path(__file__).resolve().parent
    )


def g_run(source, destination, progress, language="zh"):
    from compact_outputs import run_gvalues

    return run_gvalues(source, destination, progress, language)


def execute(
    mode,
    inputs,
    destination,
    scope="auto",
    display="Shannon",
    progress=print,
    permutations=9999,
    language="zh",
):
    if not inputs:
        raise InputSelectionError("尚未选择输入文件或目录。")
    if mode == "gvalues":
        return g_run(inputs[0], destination, progress, language)
    from analysis_workflow import run_project
    from experiment_config import discover

    experiments = discover(inputs)
    if mode in ("animal", "clinical"):
        experiments = [e for e in experiments if e.kind == mode]
    if mode == "clinical" and scope not in ("auto", "all"):
        experiments = [e for e in experiments if e.name == scope]
    return run_project(experiments, destination, display, progress, permutations, language=language)


def gui(smoke=None, language=None):
    from desktop_app import launch

    launch(execute, smoke, language)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--run", action="store_true")
    p.add_argument("--mode", choices=["auto", "clinical", "animal", "gvalues"], default="auto")
    p.add_argument("--project")
    p.add_argument("--input", action="append")
    p.add_argument("--output")
    p.add_argument("--scope", choices=list(SCOPES.values()), default="auto")
    p.add_argument("--display", choices=["Shannon", "Renyi2", "Energy", "none"], default=None)
    p.add_argument("--language", choices=["zh", "en"], default=None)
    p.add_argument("--ui-smoke-test")
    p.add_argument("--permutations", type=int, choices=[0, 999, 9999], default=None)
    a = p.parse_args()
    if a.ui_smoke_test:
        gui(a.ui_smoke_test, a.language)
        return
    if not a.run:
        gui(language=a.language)
        return
    if not a.input and (not a.project) or not a.output:
        p.error("--run requires --input or --project, and --output")
    target = Path(a.output)

    def log(s):
        if sys.stdout:
            try:
                print(s, flush=True)
            except (OSError, UnicodeError):
                pass

    try:
        if a.project:
            from experiment_config import load_project
            from analysis_workflow import run_project

            experiments, config = load_project(a.project)
            run_project(
                experiments,
                target,
                a.display or config.get("entropy_display", "Shannon"),
                log,
                a.permutations if a.permutations is not None else config.get("permutations", 9999),
                language=a.language or config.get("language", "zh"),
            )
        else:
            execute(
                a.mode,
                a.input,
                target,
                a.scope,
                a.display or "Shannon",
                log,
                a.permutations if a.permutations is not None else 9999,
                language=a.language or "zh",
            )
    except Exception:
        detail = traceback.format_exc()
        log(detail)
        target.mkdir(parents=True, exist_ok=True)
        (target / "错误信息.txt").write_text(detail, "utf-8")
        write_json(target / "运行状态.json", dict(status="failed", version=VERSION))
        sys.exit(1)


if __name__ == "__main__":
    main()
