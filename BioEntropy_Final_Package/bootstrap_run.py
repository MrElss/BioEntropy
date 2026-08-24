
from __future__ import annotations
import subprocess, sys, venv
from pathlib import Path

ROOT = Path(__file__).resolve().parent
VENV = ROOT / ".venv_run"
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
    run([py, "-m", "pip", "install", "-r", REQ])
    run([py, APP])

if __name__ == "__main__":
    main()
