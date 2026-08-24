"""Display state shared by the orchestrator and the plotting code.

Holds the group display order and the experiment being rendered. Set by
bioentropy_core_generalized around a render pass; read by bioentropy_plots.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

DISPLAY_GROUP_ORDER_MAP: Dict[str, List[str]] = {}
CURRENT_DISPLAY_EXPERIMENT: Optional[str] = None
# Per-sample table for the experiment currently being rendered. Set by the
# orchestrator around a render pass so single-timepoint bar charts can overlay
# each subject's raw value on top of the group bar without threading the frame
# through the whole plotting call stack. None ⇒ no overlay.
CURRENT_SAMPLE_LEVEL: Optional[Any] = None
