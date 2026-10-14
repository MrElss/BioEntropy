
from __future__ import annotations
import subprocess, sys, venv
from pathlib import Path

ROOT = Path(__file__).resolve().parent
VENV = ROOT / ".venv_run"
REQ = ROOT / "requirements_desktop.txt"
LOCK = ROOT / "requirements_lock.txt"
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
    # Preserve an existing working runtime; ordinary launches never upgrade it.
    probe = subprocess.run(
        [str(py), "-c", "import numpy, pandas, scipy, sklearn, matplotlib, openpyxl, xlsxwriter, plotly, tkinter"],
        cwd=str(ROOT), capture_output=True, text=True,
    )
    if probe.returncode != 0:
        # The captured lock is Python 3.12+; Python 3.11 retains the supported
        # flexible requirements instead of installing incompatible wheels.
        supports_lock = subprocess.run(
            [str(py), "-c", "import sys; sys.exit(0 if sys.version_info >= (3, 12) else 1)"],
            capture_output=True,
        ).returncode == 0
        requirement = LOCK if supports_lock and LOCK.exists() else REQ
        print("Preparing runtime from:", requirement.name)
        run([py, "-m", "pip", "install", "-r", requirement])
    run([py, APP])

if __name__ == "__main__":
    main()
