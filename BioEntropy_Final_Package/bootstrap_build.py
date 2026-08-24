
from __future__ import annotations
import shutil, subprocess, sys, venv
from pathlib import Path

ROOT = Path(__file__).resolve().parent
VENV = ROOT / ".venv_build"
REQ = ROOT / "requirements_desktop.txt"
APP = ROOT / "BioEntropy_Desktop_App.py"

def vpy() -> Path:
    if sys.platform.startswith("win"):
        return VENV / "Scripts" / "python.exe"
    return VENV / "bin" / "python"

def ensure_venv():
    if not vpy().exists():
        print("Creating virtual environment:", VENV)
        venv.EnvBuilder(with_pip=True).create(VENV)

def run(cmd):
    print(">", " ".join(str(x) for x in cmd))
    subprocess.check_call([str(c) for c in cmd], cwd=str(ROOT))

def main():
    ensure_venv()
    py = vpy()
    run([py, "-m", "pip", "install", "--upgrade", "pip", "setuptools", "wheel"])
    run([py, "-m", "pip", "install", "-r", REQ, "pyinstaller>=6.6"])
    for folder in ["build", "dist"]:
        shutil.rmtree(ROOT / folder, ignore_errors=True)
    sep = ";" if sys.platform.startswith("win") else ":"
    cmd = [
        py, "-m", "PyInstaller",
        "--noconfirm",
        "--clean",
        "--onefile",
        "--windowed",
        "--name", "BioEntropy_Desktop",
        "--add-data", f"bioentropy_core_generalized.py{sep}.",
        "--add-data", f"bioentropy_core_legacy.py{sep}.",
        "--add-data", f"gvalue_recompute.py{sep}.",
        "--add-data", f"bioentropy_constants.py{sep}.",
        "--add-data", f"bioentropy_entropy.py{sep}.",
        "--add-data", f"bioentropy_display_state.py{sep}.",
        "--add-data", f"bioentropy_plots.py{sep}.",
        "--collect-all", "matplotlib",
        "--collect-all", "sklearn",
        "--collect-all", "scipy",
        "--collect-all", "openpyxl",
        "--collect-all", "pandas",
        "--collect-all", "plotly",
        str(APP),
    ]
    run(cmd)
    print("Build complete. EXE should be in:", ROOT / "dist" / "BioEntropy_Desktop.exe")

if __name__ == "__main__":
    main()
