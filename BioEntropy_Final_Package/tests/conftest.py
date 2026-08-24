"""Shared pytest config: make the package importable and force a headless
matplotlib backend so plotting code never needs a display."""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

PKG = Path(__file__).resolve().parent.parent
if str(PKG) not in sys.path:
    sys.path.insert(0, str(PKG))
