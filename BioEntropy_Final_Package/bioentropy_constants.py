"""Shared constants: version string, numerical parameters, figure palette."""
from __future__ import annotations

# Keep in sync with APP_VERSION in the desktop app.
SOFTWARE_VERSION = "V1.6"

# Numerical / algorithm configuration (shared by the core and the entropy module).
RNG_SEED = 42
WINSOR_Q_DEFAULT = 0.01
ENTROPY_STRONG_OUTLIER_METHODS = 4
ENTROPY_OUTLIER_SHRINK = 0.5
NMD_Z_CAP = 8.0

# Nature Publishing Group (NPG) qualitative palette used across all figures.
_NPG_COLORS = [
    "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
    "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
]
_FALLBACK_COLORS = _NPG_COLORS

# Semantic colors for the core indices / group roles.
_NMD_COLOR = "#3C5488"
_ENTROPY_COLOR = "#E64B35"
_BRI_COLOR = "#00A087"
_MODEL_COLOR = "#DC0000"
_NORMAL_COLOR = "#8491B4"
_NORMALIZED_COLOR = "#4DBBD5"
_GRID_COLOR = "#D8DEE8"
_TEXT_COLOR = "#202A35"
_LIGHT_BG = "#FAFBFC"

# Gamma applied to the 0-100 anchored display values to spread high-burden bins.
DISPLAY_CONTRAST_GAMMA = 3.0
