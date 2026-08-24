"""Recalculation of Hedges' g (G-values) for the forest-plot datasets.

This standalone analysis utility recomputes the G-values from the refreshed
raw data under ``data/G值森林图数据``.  It covers the first two deliverables:

    1) 所有单个指标G值  – per-metric Hedges' g for every treatment arm
                          (efficacy of each arm vs the disease Model, expressed
                          in "improvement" direction = movement toward Control).
    2) Pooled G值计算    – per-metric directional Hedges' g for HTD1801 vs PM,
                          plus inverse-variance pooled (fixed + random effects)
                          G across the metrics of each animal model.

The Hedges' g formula follows the standard small-sample-corrected definition:

    Sp        = sqrt(((n1-1)·SD1² + (n2-1)·SD2²) / (n1+n2-2))
    Cohen's d = (mean1 - mean2) / Sp
    J         = 1 - 3 / (4·(n1+n2-2) - 1)
    Hedges' g = J · d
    Var(g)    = J² · ( (n1+n2)/(n1·n2) + d²/(2·(n1+n2)) )
    SE(g)     = sqrt(Var(g)),  95% CI = g ± 1.96·SE

Requires openpyxl only. Run from the repository root:

    python BioEntropy_Final_Package/gvalue_recompute.py
"""

from __future__ import annotations

import csv
import math
import os
import re
from dataclasses import dataclass, field
from datetime import datetime
from typing import Dict, List, Optional, Sequence, Set, Tuple

import openpyxl

# --------------------------------------------------------------------------- #
# Configuration – the contested decisions are exposed here so they can be
# changed without touching the statistics below.
# --------------------------------------------------------------------------- #

# How to treat values stored as text with a trailing "*" (e.g. "515*").
#   "include" – strip the "*" and keep the value (DEFAULT: never auto-exclude;
#               the meaning of the asterisk must be confirmed by the user).
#   "exclude" – drop the value from its group (only after explicit confirmation).
# Every starred value is listed in the flagged-values QC table regardless.
ASTERISK_MODE = "include"

# Columns whose canonical name starts with this prefix are treated as
# wild-type / drug-on-healthy safety controls and excluded from efficacy G.
SAFETY_PREFIX = "WT"

# z value for the 95% confidence interval (normal approximation).
Z95 = 1.959963984540054

# Context endpoints CONFIRMED by the user as higher-in-Model (disease) and
# reduced by treatment, i.e. lower_better. They are now treated as objective
# higher_worse and enter the main analysis (still validity-checked: included
# only where the Model is actually worse than Control).
CONFIRMED_HIGHER_WORSE_TOKENS = (
    "insulin", "胰岛", "体重", "bodyweight", "bodywt", "weight",
)
# Endpoints whose disease direction still depends on context (NOT confirmed):
# food / water intake and body-weight CHANGE RATE. These stay
# "needs_manual_review" and are kept OUT of the headline/main results.
CONTEXT_ENDPOINT_TOKENS = (
    "摄食", "feed", "foodintake", "饮水", "waterintake", "water",
    "变化率", "changerate",
)
# Objective endpoints whose disease direction is fixed; main analysis keeps a
# (model, endpoint) only if the observed Model-vs-Control change matches this.
LOWER_WORSE_TOKENS = ("gfr",)            # lower = worse (must fall in disease)
HIGHER_WORSE_TOKENS = (
    "glu", "血糖", "hba1c", "hbalc", "ogtt", "ptt", "itt", "auc",
    "tg", "tc", "cho", "ldl", "alt", "ast", "livertg", "uacr", "kim",
    "ngal", "cre", "bun", "crp", "尿糖", "urineglu",
)

# The 5 DB metrics with an explicitly specified direction.
#   "lower_better"  : improvement = J·(mean_PM - mean_HTD1801)/Sp
#   "higher_better" : improvement = J·(mean_HTD1801 - mean_PM)/Sp
DB_DIRECTION = {
    "HbA1c": "lower_better",
    "ITT-AUC": "lower_better",
    "UACR": "lower_better",
    "KIM-1": "lower_better",
    "GFR": "higher_better",
}

# Which metric file (folder 2) supplies each DB pooled metric, and which sheet.
# DB-KIM-1.xlsx bundles all five metrics as sheets; the others are single-sheet.
DB_POOLED_METRICS = ["HbA1c", "ITT-AUC", "UACR", "KIM-1", "GFR"]


# --------------------------------------------------------------------------- #
# Group-name canonicalisation
# --------------------------------------------------------------------------- #

CONTROL_ALIASES = {"control", "ct", "normal", "normalcontrol"}
MODEL_ALIASES = {"model", "m", "t2dm", "modelcontrol"}
PM_ALIASES = {"pm", "physicalmixture"}


def canonical_group(name: object) -> str:
    """Map a raw column header to a canonical group label."""
    raw = str(name or "").replace("\xa0", " ").strip()
    key = raw.lower().replace("_", "-").replace(" ", "")
    if key in CONTROL_ALIASES:
        return "Control"
    if key in MODEL_ALIASES:
        return "Model"
    if key in {"1801", "htd1801"}:
        return "HTD1801"
    if key in {"htd1801-l", "1801-l"}:
        return "HTD1801-L"
    if key in {"htd1801-m", "1801-m"}:
        return "HTD1801-M"
    if key in {"htd1801-h", "1801-h"}:
        return "HTD1801-H"
    if key in PM_ALIASES:
        return "PM"
    if key in {"met", "metformin"}:
        return "MET"
    if key in {"udca"}:
        return "UDCA"
    if key in {"bbr"}:
        return "BBR"
    # Wild-type / safety controls: WT1801, WTBBR, WT UDCA, ...
    if key.startswith("wt"):
        return SAFETY_PREFIX + raw[2:].strip().lstrip(" \xa0")
    return raw


def _norm_token(name: object) -> str:
    return str(name or "").strip().lower().replace(" ", "").replace("_", "").replace("-", "")


def canonical_endpoint(name: object) -> str:
    """Standardise an endpoint label so the same marker matches across models.

    e.g. 血糖 / GLU / Blood glucose → "Blood glucose"; 血清胰岛素 / insulin →
    "Serum insulin"; CHO / TC → "TC"; HBA1C / 糖化血红蛋白 → "HbA1c";
    Urine Glucos / 尿糖 → "Urine glucose"; liver TG → "Liver TG".
    """
    raw = str(name or "").strip()
    k = _norm_token(name)
    if "变化率" in raw or "changerate" in k:
        return "Body weight change rate"
    if "尿糖" in raw or "urineglu" in k:
        return "Urine glucose"
    if "livertg" in k or ("liver" in k and "tg" in k):
        return "Liver TG"
    if "糖化" in raw or "hba1c" in k or "hbalc" in k:
        return "HbA1c"
    if "血糖" in raw or k == "glu" or "bloodglucose" in k or k == "bg":
        return "Blood glucose"
    if "胰岛素" in raw or "insulin" in k:
        return "Serum insulin"
    if "体重" in raw or k in ("bodyweight", "bw", "weight"):
        return "Body weight"
    if "摄食" in raw or "feed" in k or "foodintake" in k:
        return "Food intake"
    if "饮水" in raw or "waterintake" in k or k == "water":
        return "Water intake"
    if k in ("cho", "tc", "totalcholesterol"):
        return "TC"
    if k in ("ldl", "ldlc"):
        return "LDL-c"
    if "ogtt" in k:
        return "OGTT-AUC"
    if "ptt" in k:
        return "PTT-AUC"
    if "itt" in k:
        return "ITT-AUC"
    if "kim" in k:
        return "KIM-1"
    for tok, disp in (("uacr", "UACR"), ("ngal", "NGAL"), ("alt", "ALT"),
                      ("ast", "AST"), ("bun", "BUN"), ("cre", "CRE"),
                      ("crp", "CRP"), ("gfr", "GFR"), ("tg", "TG")):
        if tok in k:
            return disp
    return raw


def exclusion_reason(status: str) -> str:
    """Plain-language reason an endpoint is/ isn't in the within-model pool."""
    return {
        "valid": "",
        "model_not_valid_for_this_endpoint": "Model not worse than Control",
        "needs_manual_review": "context endpoint — direction needs confirmation",
        "direction_uncertain": "Control≈Model / direction undefined",
    }.get(status, status)


def classify_endpoint(name: object) -> Tuple[str, Optional[str]]:
    """Classify an endpoint by disease-direction handling.

    Returns (category, expected_change) where:
      category in {"objective", "context", "unknown"}
      expected_change in {"up", "down", None}  – the Model-vs-Control change
          that signals disease for an objective endpoint ("up" = higher_worse,
          "down" = lower_worse).
    Insulin and body weight are user-confirmed higher_worse objective endpoints
    (see CONFIRMED_HIGHER_WORSE_TOKENS) and enter the main analysis when the
    Model is actually worse than Control. The remaining "context" endpoints
    (food/water intake, body-weight change-rate) still need manual review and are
    kept out of headline results by default.
    """
    k = _norm_token(name)
    # body-weight CHANGE RATE stays context (not confirmed) — check before 体重.
    if any(tok in k for tok in ("变化率", "changerate")):
        return "context", None
    if "gfr" in k:
        return "objective", "down"
    # user-confirmed higher_worse endpoints (insulin, body weight).
    if any(tok in k for tok in CONFIRMED_HIGHER_WORSE_TOKENS):
        return "objective", "up"
    # remaining context (food / water intake) — needs manual review.
    if any(tok in k for tok in CONTEXT_ENDPOINT_TOKENS):
        return "context", None
    if any(tok in k for tok in HIGHER_WORSE_TOKENS):
        return "objective", "up"
    return "unknown", None


def coerce_value(value: object) -> Tuple[Optional[float], Optional[str], Optional[str]]:
    """Coerce one cell to float.

    Returns (numeric_or_None, flag, original_text):
      flag is None for clean numbers, "asterisk" for starred values, or
      "nonnumeric" for unparseable text. For starred values the numeric is
      returned when ASTERISK_MODE == "include" (default) and None otherwise;
      either way the caller logs it in the flagged-values table.
    """
    if value is None:
        return None, None, None
    if isinstance(value, (int, float)):
        return float(value), None, None
    text = str(value).strip()
    if not text:
        return None, None, None
    if text.endswith("*"):
        stripped = text[:-1].strip()
        try:
            num = float(stripped)
        except ValueError:
            return None, "nonnumeric", text
        return (num if ASTERISK_MODE == "include" else None), "asterisk", text
    try:
        return float(text), None, None
    except ValueError:
        return None, "nonnumeric", text


@dataclass
class FlaggedValue:
    model: str
    endpoint: str
    group: str
    animal: str
    original: str
    kind: str       # "asterisk" | "nonnumeric"
    action: str     # "retained" | "excluded"


@dataclass
class ValidationRecord:
    model: str
    family: str
    file: str
    status: str         # "ok" | "warning" | "error"
    groups_found: str
    messages: List[str] = field(default_factory=list)


# Groups the core HTD1801-vs-PM analysis needs.
REQUIRED_GROUPS = ("Control", "Model", "PM")


def validate_group_file(model: str, family: str, fname: str,
                        groups: Dict[str, List[float]]) -> ValidationRecord:
    """Plain-language validation of a per-metric (column = group) file."""
    found = sorted(groups.keys())
    msgs: List[str] = []
    status = "ok"
    has_htd = any(g.split("#")[0].startswith("HTD1801") for g in groups)
    missing = [g for g in REQUIRED_GROUPS if g not in groups]
    if not has_htd:
        missing = missing + ["HTD1801"]
    if missing:
        status = "warning"
        msgs.append("缺少分析所需的分组：" + "、".join(missing) +
                    "（请检查表头是否拼写一致，例如 Control / Model / HTD1801 / PM）。")
    unknown = [g for g in groups if g.split("#")[0] not in (
        "Control", "Model", "PM", "HTD1801", "HTD1801-L", "HTD1801-M",
        "HTD1801-H", "MET", "UDCA", "BBR") and not g.startswith(SAFETY_PREFIX)]
    if unknown:
        msgs.append("未识别的分组列（按原名保留，未纳入标准比较）：" + "、".join(unknown) + "。")
    return ValidationRecord(model=model, family=family, file=fname,
                            status=status, groups_found="、".join(found), messages=msgs)


# --------------------------------------------------------------------------- #
# Workbook reading
# --------------------------------------------------------------------------- #

# Positional group order assumed when a file has NO header row (the first row is
# already numeric data). Matches the standard 4-arm animal layout used throughout
# the dataset; only applied to 4-column headerless sheets, and always flagged.
_HEADERLESS_GROUP_ORDER = ("Control", "Model", "HTD1801", "PM")


def read_groups(path: str, sheet: Optional[str] = None, model: str = "",
                endpoint: Optional[str] = None) -> Tuple[Dict[str, List[float]], List[FlaggedValue]]:
    """Read a workbook sheet into {canonical_group: [values]}.

    Row 1 is normally the header row (group names). If row 1 is already numeric
    data — i.e. the file has no header — the columns are mapped positionally to
    the standard ``Control / Model / HTD1801 / PM`` layout (4-column sheets only)
    so a header-less export is not silently dropped; the assumption is recorded
    as a FlaggedValue for the QC table.

    Duplicate canonical names keep the first column and append a suffix to the
    rest so nothing is silently merged. Returns the groups plus a list of
    FlaggedValue records (starred / non-numeric cells) for the QC table.
    """
    if endpoint is None:
        endpoint = _metric_name(os.path.basename(path))
    wb = openpyxl.load_workbook(path, data_only=True, read_only=True)
    ws = wb[sheet] if sheet else wb[wb.sheetnames[0]]
    headers = [ws.cell(row=1, column=c).value for c in range(1, ws.max_column + 1)]
    groups: Dict[str, List[float]] = {}
    flagged: List[FlaggedValue] = []

    # Detect a header-less sheet: every non-empty cell in row 1 is numeric, so
    # there are no group labels to read. In that case row 1 is data, not a header.
    nonempty = [h for h in headers if not (h is None or str(h).strip() == "")]
    headerless = bool(nonempty) and all(coerce_value(h)[0] is not None for h in nonempty)

    if headerless and len(headers) == len(_HEADERLESS_GROUP_ORDER):
        labels = list(_HEADERLESS_GROUP_ORDER)
        data_start = 1
        flagged.append(FlaggedValue(
            model=model, endpoint=endpoint, group="/".join(labels),
            animal="header row", original="(no header row)", kind="headerless",
            action="assumed Control/Model/HTD1801/PM by column order",
        ))
    else:
        labels = [canonical_group(h) if not (h is None or str(h).strip() == "") else None
                  for h in headers]
        data_start = 2

    seen: Dict[str, int] = {}
    for c, (raw_header, label) in enumerate(zip(headers, labels), start=1):
        if label is None:
            continue
        if label in seen:
            seen[label] += 1
            label = f"{label}#{seen[label]}"
        else:
            seen[label] = 1
        group_name = label if headerless else str(raw_header).strip()
        col_values: List[float] = []
        for r in range(data_start, ws.max_row + 1):
            num, flag, original = coerce_value(ws.cell(row=r, column=c).value)
            if flag is not None:
                flagged.append(FlaggedValue(
                    model=model, endpoint=endpoint, group=group_name,
                    animal=f"row {r}", original=original or "", kind=flag,
                    action="retained" if num is not None else "excluded",
                ))
            if num is not None:
                col_values.append(num)
        if col_values:
            groups[label] = col_values
    wb.close()
    return groups, flagged


# --------------------------------------------------------------------------- #
# Statistics
# --------------------------------------------------------------------------- #

@dataclass
class GroupStat:
    name: str
    n: int
    mean: float
    sd: float


def describe(name: str, values: Sequence[float]) -> Optional[GroupStat]:
    n = len(values)
    if n < 2:
        return None
    mean = sum(values) / n
    var = sum((x - mean) ** 2 for x in values) / (n - 1)
    return GroupStat(name=name, n=n, mean=mean, sd=math.sqrt(var))


@dataclass
class GResult:
    """Hedges' g for sample1 - sample2 (raw signed)."""
    g: float
    se: float
    ci_low: float
    ci_high: float
    d: float
    pooled_sd: float
    s1: GroupStat
    s2: GroupStat


def _finite(values: Sequence[float]) -> List[float]:
    """Drop None / NaN / inf so a stray non-finite cell cannot poison the g."""
    out: List[float] = []
    for v in values:
        if v is None:
            continue
        try:
            f = float(v)
        except (TypeError, ValueError):
            continue
        if math.isfinite(f):
            out.append(f)
    return out


def hedges_g(values1: Sequence[float], values2: Sequence[float],
             name1: str = "g1", name2: str = "g2") -> Optional[GResult]:
    """Hedges' g for (mean1 - mean2). Formula identical to the core package.

    Non-finite inputs are dropped first so a single NaN/inf cannot propagate
    into the effect size.
    """
    s1 = describe(name1, _finite(values1))
    s2 = describe(name2, _finite(values2))
    if s1 is None or s2 is None:
        return None
    df = s1.n + s2.n - 2
    if df <= 0:
        return None
    pooled_var = ((s1.n - 1) * s1.sd ** 2 + (s2.n - 1) * s2.sd ** 2) / df
    pooled = math.sqrt(pooled_var)
    if not math.isfinite(pooled) or pooled <= 0:
        return None
    d = (s1.mean - s2.mean) / pooled
    j = 1.0 - 3.0 / (4.0 * df - 1.0) if df > 1 else 1.0
    g = j * d
    variance = j * j * ((s1.n + s2.n) / (s1.n * s2.n) + d * d / (2.0 * (s1.n + s2.n)))
    se = math.sqrt(max(variance, 0.0))
    return GResult(
        g=g, se=se, ci_low=g - Z95 * se, ci_high=g + Z95 * se,
        d=d, pooled_sd=pooled, s1=s1, s2=s2,
    )


def directed(res: GResult, multiplier: float) -> GResult:
    """Flip a raw g (and its CI) by +1/-1 to express it in improvement units."""
    lo, hi = sorted((multiplier * res.ci_low, multiplier * res.ci_high))
    return GResult(
        g=multiplier * res.g, se=res.se, ci_low=lo, ci_high=hi,
        d=multiplier * res.d, pooled_sd=res.pooled_sd, s1=res.s1, s2=res.s2,
    )


# two-sided t critical values at 0.975 by degrees of freedom (small-sample,
# used for the Hartung-Knapp adjustment; falls back to 1.96 for large df).
_T975 = {1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447,
         7: 2.365, 8: 2.306, 9: 2.262, 10: 2.228, 11: 2.201, 12: 2.179,
         13: 2.160, 14: 2.145, 15: 2.131, 16: 2.120, 17: 2.110, 18: 2.101,
         19: 2.093, 20: 2.086, 25: 2.060, 30: 2.042}


def _t975(df: int) -> float:
    if df <= 0:
        return float("nan")
    if df in _T975:
        return _T975[df]
    if df > 30:
        return Z95
    keys = sorted(_T975)
    return _T975[min(keys, key=lambda k: abs(k - df))]


@dataclass
class Pooled:
    label: str
    k: int
    estimate: float
    se: float
    ci_low: float
    ci_high: float
    q: float = 0.0
    tau2: float = 0.0
    i2: float = 0.0
    # Hartung-Knapp adjusted random-effects CI (recommended for few studies).
    hk_se: Optional[float] = None
    hk_ci_low: Optional[float] = None
    hk_ci_high: Optional[float] = None

    @property
    def hk_crosses_zero(self) -> Optional[bool]:
        if self.hk_ci_low is None or self.hk_ci_high is None:
            return None
        return self.hk_ci_low <= 0.0 <= self.hk_ci_high


def inverse_variance_pool(effects: Sequence[float], ses: Sequence[float]) -> Tuple[Pooled, Pooled]:
    """Return (fixed_effect, random_effect) pooled estimates.

    The random-effects estimate is DerSimonian-Laird; its Pooled also carries a
    Hartung-Knapp adjusted CL (t-based, df=k-1) which is the recommended honest
    interval when only a few studies/models are pooled.
    """
    eff = [e for e, s in zip(effects, ses) if s and math.isfinite(s) and s > 0
           and math.isfinite(e)]
    se = [s for e, s in zip(effects, ses) if s and math.isfinite(s) and s > 0
          and math.isfinite(e)]
    k = len(eff)
    w = [1.0 / (s * s) for s in se]
    sw = sum(w)
    fixed = sum(wi * ei for wi, ei in zip(w, eff)) / sw
    fixed_se = math.sqrt(1.0 / sw)
    q = sum(wi * (ei - fixed) ** 2 for wi, ei in zip(w, eff)) if k > 1 else 0.0
    dfq = max(k - 1, 0)
    c = sw - sum(wi * wi for wi in w) / sw if k > 1 else 0.0
    tau2 = max(0.0, (q - dfq) / c) if c > 0 else 0.0
    wr = [1.0 / (s * s + tau2) for s in se]
    swr = sum(wr)
    rand = sum(wi * ei for wi, ei in zip(wr, eff)) / swr
    rand_se = math.sqrt(1.0 / swr)
    i2 = max(0.0, (q - dfq) / q * 100.0) if q > 0 and dfq > 0 else 0.0
    # Hartung-Knapp variance: (1/(k-1)) * Σ wr_i (y_i-μ)² / Σ wr_i, t(k-1) CI.
    hk_se = hk_lo = hk_hi = None
    if k >= 2:
        qhk = sum(wi * (ei - rand) ** 2 for wi, ei in zip(wr, eff)) / (k - 1)
        hk_se = math.sqrt(qhk / swr)
        # Sidik-Jonkman safeguard: do not let HK be narrower than classic RE.
        hk_se = max(hk_se, rand_se)
        tcrit = _t975(k - 1)
        hk_lo, hk_hi = rand - tcrit * hk_se, rand + tcrit * hk_se
    fixed_p = Pooled("fixed", k, fixed, fixed_se,
                     fixed - Z95 * fixed_se, fixed + Z95 * fixed_se, q, 0.0, i2)
    rand_p = Pooled("random", k, rand, rand_se,
                    rand - Z95 * rand_se, rand + Z95 * rand_se, q, tau2, i2,
                    hk_se=hk_se, hk_ci_low=hk_lo, hk_ci_high=hk_hi)
    return fixed_p, rand_p


def two_sided_p(g: float, se: float) -> float:
    """Two-sided p-value for an effect size via the normal approximation."""
    if se is None or not math.isfinite(se) or se <= 0:
        return float("nan")
    z = abs(g) / se
    # erfc-based normal tail; p = 2*(1 - Phi(z)) = erfc(z/sqrt2)
    return math.erfc(z / math.sqrt(2.0))


def fmt_p(p: float) -> str:
    """Compact p-value string for tables/plots (e.g. '<0.001', '0.021', 'n/a')."""
    if p is None or not math.isfinite(p):
        return "n/a"
    if p < 0.001:
        return "<0.001"
    if p < 0.01:
        return f"{p:.3f}"
    return f"{p:.2f}"


def benjamini_hochberg(pvals: Sequence[float]) -> List[float]:
    """Benjamini-Hochberg FDR-adjusted q-values (NaNs preserved)."""
    idx = [i for i, p in enumerate(pvals) if p is not None and math.isfinite(p)]
    m = len(idx)
    q = [float("nan")] * len(pvals)
    if m == 0:
        return q
    order = sorted(idx, key=lambda i: pvals[i])
    prev = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        r = m - rank + 1
        val = min(prev, pvals[i] * m / r)
        q[i] = val
        prev = val
    return q


# --------------------------------------------------------------------------- #
# Direction / validity checks
# --------------------------------------------------------------------------- #

@dataclass
class DirectionCheck:
    analysis: str          # which family raised it (Part 2 / Mahalanobis / ...)
    model: str
    endpoint: str
    n_control: int
    n_model: int
    control_mean: float
    model_mean: float
    observed: str          # "up" / "down" / "flat"
    category: str          # objective / context / unknown
    expected: Optional[str]    # up / down / None
    status: str            # valid / model_not_valid_for_this_endpoint /
                           # needs_manual_review / direction_uncertain
    included_in_main: bool


def assess_endpoint(model: str, endpoint: str, control_vals: Sequence[float],
                    model_vals: Sequence[float], analysis: str = "") -> DirectionCheck:
    """Infer disease direction from Control vs Model and decide main inclusion.

    Rules:
      * objective endpoints (fixed disease direction) are INCLUDED in the main
        analysis only if the Model deteriorates in the expected direction;
        otherwise -> model_not_valid_for_this_endpoint.
      * insulin and body weight are user-confirmed objective higher_worse
        endpoints and are INCLUDED when the Model deteriorates as expected.
      * remaining context endpoints (food/water intake, body-weight change rate)
        -> needs_manual_review and excluded from the main analysis by default.
      * if Control/Model are missing or equal -> direction_uncertain.
    """
    cs = describe("Control", _finite(control_vals))
    ms = describe("Model", _finite(model_vals))
    category, expected = classify_endpoint(endpoint)
    cmean = cs.mean if cs else float("nan")
    mmean = ms.mean if ms else float("nan")
    nC = cs.n if cs else 0
    nM = ms.n if ms else 0
    if cs is None or ms is None or cmean == mmean:
        observed = "flat"
        status = "direction_uncertain"
        return DirectionCheck(analysis, model, endpoint, nC, nM, cmean, mmean,
                              observed, category, expected, status, False)
    observed = "up" if mmean > cmean else "down"
    if category == "context":
        status = "needs_manual_review"
        included = False
    elif category == "objective":
        if observed == expected:
            status = "valid"
            included = True
        else:
            status = "model_not_valid_for_this_endpoint"
            included = False
    else:  # unknown endpoint name
        status = "direction_uncertain"
        included = False
    return DirectionCheck(analysis, model, endpoint, nC, nM, cmean, mmean,
                          observed, category, expected, status, included)


def control_anchor_multiplier(groups: Dict[str, List[float]]) -> Optional[float]:
    """+1 if higher=better (Control>Model), -1 if lower=better (Control<Model)."""
    control = describe("Control", groups.get("Control", []))
    model = describe("Model", groups.get("Model", []))
    if control is None or model is None:
        return None
    if control.mean == model.mean:
        return None
    return 1.0 if control.mean > model.mean else -1.0


def direction_multiplier(metric: str, groups: Dict[str, List[float]]) -> Tuple[Optional[float], str]:
    """Return (+1/-1 multiplier for HTD1801-vs-PM improvement, basis-string)."""
    if metric in DB_DIRECTION:
        if DB_DIRECTION[metric] == "lower_better":
            return -1.0, "lower_better (explicit)"
        return 1.0, "higher_better (explicit)"
    mult = control_anchor_multiplier(groups)
    if mult is None:
        return None, "Control-anchor (unavailable: missing Control/Model)"
    basis = "higher_better" if mult > 0 else "lower_better"
    return mult, f"{basis} (Control-anchor)"


# --------------------------------------------------------------------------- #
# Path discovery (avoids hard-coding non-ASCII literals)
# --------------------------------------------------------------------------- #

def find_gvalue_dir(repo_root: str) -> str:
    # Accept the G值森林图数据 directory itself (it holds the 1）..5） subfolders).
    if os.path.isdir(repo_root) and any(
            k[:1] in "12345" and os.path.isdir(os.path.join(repo_root, k))
            for k in os.listdir(repo_root)):
        return repo_root
    # Otherwise look one level down under data/.
    for base in (os.path.join(repo_root, "data"), repo_root):
        if not os.path.isdir(base):
            continue
        for d in sorted(os.listdir(base)):
            full = os.path.join(base, d)
            if os.path.isdir(full) and any(k[:1] in "12345" for k in os.listdir(full)):
                return full
    raise FileNotFoundError("Could not locate the G值森林图数据 directory under data/")


def _optional_subdir(parent: str, prefix: str) -> Optional[str]:
    """Return the path of the first subdir whose name starts with ``prefix``,
    or None if absent. The five G-value families are independent, so any
    subset (or a single one) may be present."""
    for d in sorted(os.listdir(parent)):
        if d.startswith(prefix) and os.path.isdir(os.path.join(parent, d)):
            return os.path.join(parent, d)
    return None


def _subdir_startswith(parent: str, prefix: str) -> str:
    found = _optional_subdir(parent, prefix)
    if found is None:
        raise FileNotFoundError(f"No subdir starting with {prefix!r} in {parent}")
    return found


def model_subdirs(folder: str) -> Dict[str, str]:
    """Map a normalised model key -> path for the model subfolders."""
    out: Dict[str, str] = {}
    for d in sorted(os.listdir(folder)):
        full = os.path.join(folder, d)
        if not os.path.isdir(full):
            continue
        low = d.lower()
        if low.startswith("db") or "db" in low:
            out["DB"] = full
        elif low.startswith("kk") or "kk" in low:
            out["KK-ay"] = full
        else:
            out["地鼠"] = full
    return out


# --------------------------------------------------------------------------- #
# Part 1 – per-metric, per-arm efficacy G
# --------------------------------------------------------------------------- #

TREATMENT_ARMS = ("HTD1801", "HTD1801-L", "HTD1801-M", "HTD1801-H",
                  "PM", "MET", "UDCA", "BBR")


@dataclass
class Part1Record:
    """One individual model × endpoint HTD1801-vs-PM efficacy G (no pooling).

    This is the Part-1 unit: a single endpoint's effect size within one model,
    g>0 ⇒ HTD1801 better than PM. Validity (does the Model deteriorate vs
    Control?) decides whether the endpoint feeds the Part-2 within-model pool.
    """
    model: str
    endpoint: str
    dose: str
    category: str            # objective / context
    status: str              # valid / model_not_valid_for_this_endpoint / ...
    included: bool           # contributes to the Part-2 within-model pool
    g: float
    se: float
    ci_low: float
    ci_high: float
    htd_mean: float
    htd_n: int
    pm_mean: float
    pm_n: int
    model_vs_control_g: Optional[float]


# --------------------------------------------------------------------------- #
# Part 2 – pooled HTD1801 vs PM
# --------------------------------------------------------------------------- #

@dataclass
class Part2Metric:
    model: str
    metric: str
    direction_basis: str
    stats: Dict[str, GroupStat] = field(default_factory=dict)
    htd_vs_pm: Optional[GResult] = None          # directional (improvement)
    model_vs_control: Optional[GResult] = None   # raw signed
    htd_vs_model: Optional[GResult] = None        # improvement
    pm_vs_model: Optional[GResult] = None         # improvement
    htd_label: str = "HTD1801"


def _improve_vs_model(groups, arm, mult) -> Optional[GResult]:
    if arm not in groups or "Model" not in groups:
        return None
    res = hedges_g(groups[arm], groups["Model"], arm, "Model")
    return directed(res, mult) if res else None


def compute_part2_metric(model: str, metric: str, groups: Dict[str, List[float]],
                         htd_group: str = "HTD1801") -> Tuple[Optional[Part2Metric], List[str]]:
    warns: List[str] = []
    mult, basis = direction_multiplier(metric, groups)
    if mult is None:
        warns.append(f"{model}/{metric}: direction undefined ({basis}); skipped.")
        return None, warns
    if htd_group not in groups or "PM" not in groups:
        warns.append(f"{model}/{metric}: missing {htd_group} or PM; skipped.")
        return None, warns

    out = Part2Metric(model=model, metric=metric, direction_basis=basis, htd_label=htd_group)
    for key in ("Control", "Model", htd_group, "PM"):
        st = describe(key, groups.get(key, []))
        if st is not None:
            out.stats[key] = st

    raw = hedges_g(groups[htd_group], groups["PM"], htd_group, "PM")
    if raw is None:
        warns.append(f"{model}/{metric}: {htd_group} vs PM not computable.")
        return None, warns
    out.htd_vs_pm = directed(raw, mult)

    mvc = hedges_g(groups.get("Model", []), groups.get("Control", []), "Model", "Control")
    out.model_vs_control = mvc
    cmult = control_anchor_multiplier(groups)
    if cmult is not None:
        out.htd_vs_model = _improve_vs_model(groups, htd_group, cmult)
        out.pm_vs_model = _improve_vs_model(groups, "PM", cmult)
    return out, warns


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #

def _list_xlsx(folder: str) -> List[str]:
    return sorted(f for f in os.listdir(folder) if f.endswith(".xlsx"))


def _metric_name(filename: str) -> str:
    stem = os.path.splitext(filename)[0]
    for pref in ("DB-", "DB鼠-", "kk-", "KK-", "地鼠-", "地鼠- "):
        if stem.startswith(pref):
            stem = stem[len(pref):]
    return stem.strip()


def _kk_repr_dose(groups_or_labels) -> str:
    """The representative HTD1801 dose label present (medium for KK-ay)."""
    if PART3_PRIMARY_KK_DOSE in groups_or_labels:
        return PART3_PRIMARY_KK_DOSE
    for d in ("HTD1801-M", "HTD1801-H", "HTD1801-L", "HTD1801"):
        if d in groups_or_labels:
            return d
    return "HTD1801"


def run(repo_root: str) -> Dict[str, object]:
    gdir = find_gvalue_dir(repo_root)
    folder1 = _optional_subdir(gdir, "1")
    folder2 = _optional_subdir(gdir, "2")

    all_flagged: List[FlaggedValue] = []
    all_warns: List[str] = []
    direction_checks: List[DirectionCheck] = []
    validations: List[ValidationRecord] = []

    # ---- Part 1 = individual model × endpoint G (HTD1801 vs PM, no pooling) --
    # Source: folder 1 (所有单个指标) — the comprehensive per-endpoint data.
    part1_records: List[Part1Record] = []
    for model, sub in (model_subdirs(folder1).items() if folder1 else []):
        for fname in _list_xlsx(sub):
            endpoint = canonical_endpoint(_metric_name(fname))
            groups, flagged = read_groups(os.path.join(sub, fname), model=model, endpoint=endpoint)
            all_flagged += flagged
            validations.append(validate_group_file(model, "1）单个指标", fname, groups))
            dc = assess_endpoint(model, endpoint, groups.get("Control", []),
                                 groups.get("Model", []), "1）单个指标")
            direction_checks.append(dc)
            dose = _kk_repr_dose(list(groups)) if model == "KK-ay" else (
                "HTD1801" if "HTD1801" in groups else _kk_repr_dose(list(groups)))
            res, warns = compute_part2_metric(model, endpoint, groups, htd_group=dose)
            all_warns += warns
            if res is None or res.htd_vs_pm is None:
                continue
            htd = res.stats.get(dose)
            pm = res.stats.get("PM")
            hp = res.htd_vs_pm
            part1_records.append(Part1Record(
                model=model, endpoint=endpoint, dose=dose, category=dc.category,
                status=dc.status, included=dc.included_in_main,
                g=hp.g, se=hp.se, ci_low=hp.ci_low, ci_high=hp.ci_high,
                htd_mean=htd.mean if htd else float("nan"), htd_n=htd.n if htd else 0,
                pm_mean=pm.mean if pm else float("nan"), pm_n=pm.n if pm else 0,
                model_vs_control_g=res.model_vs_control.g if res.model_vs_control else None))

    # ---- Part 2 = within-model pooled G (PRIMARY: folder 2 「pooled G值计算」) --
    # The dedicated within-model pooling set. One pooled G per model
    # (random-effects + Hartung-Knapp) over its valid endpoints.
    part2_within: Dict[str, Tuple[Pooled, Pooled]] = {}
    part2_within_members: Dict[str, List[str]] = {}
    part2_within_points: Dict[str, List[Dict[str, float]]] = {}
    for model, sub in (model_subdirs(folder2).items() if folder2 else []):
        items: List[Part2Metric] = []
        for fname in _list_xlsx(sub):
            metric_raw = _metric_name(fname)
            sheet = "KIM-1" if (model == "DB" and metric_raw.upper() == "KIM-1") else None
            metric = canonical_endpoint(metric_raw)
            groups, flagged = read_groups(os.path.join(sub, fname), sheet=sheet,
                                          model=model, endpoint=metric)
            all_flagged += flagged
            validations.append(validate_group_file(model, "2）Pooled", fname, groups))
            dc = assess_endpoint(model, metric, groups.get("Control", []),
                                 groups.get("Model", []), "2）Pooled")
            direction_checks.append(dc)
            dose = _kk_repr_dose(list(groups)) if model == "KK-ay" else (
                "HTD1801" if "HTD1801" in groups else _kk_repr_dose(list(groups)))
            res, warns = compute_part2_metric(model, metric, groups, htd_group=dose)
            all_warns += warns
            if res is not None and res.htd_vs_pm is not None and dc.included_in_main:
                items.append(res)
        if items:
            part2_within[model] = inverse_variance_pool(
                [r.htd_vs_pm.g for r in items], [r.htd_vs_pm.se for r in items])
            part2_within_members[model] = [r.metric for r in items]
            # keep each endpoint's individual g so the within-model figure can show
            # the per-indicator points above the pooled diamond (Part-3 style).
            part2_within_points[model] = [{
                "metric": r.metric, "g": r.htd_vs_pm.g, "se": r.htd_vs_pm.se,
                "lo": r.htd_vs_pm.ci_low, "hi": r.htd_vs_pm.ci_high,
            } for r in items]

    # ---- Part 2 (sensitivity) = comprehensive pool of ALL valid Part-1 endpoints
    part2_comprehensive: Dict[str, Tuple[Pooled, Pooled]] = {}
    part2_comprehensive_members: Dict[str, List[str]] = {}
    for model in dict.fromkeys(r.model for r in part1_records):
        members = [r for r in part1_records if r.model == model and r.included]
        if members:
            part2_comprehensive[model] = inverse_variance_pool(
                [r.g for r in members], [r.se for r in members])
            part2_comprehensive_members[model] = [r.endpoint for r in members]

    # ---- Part 3 -----------------------------------------------------------
    part3 = run_part3(repo_root, all_warns, all_flagged, direction_checks, validations)

    # ---- Part 4 -----------------------------------------------------------
    part4 = run_part4(repo_root, all_warns, all_flagged, direction_checks, validations)

    # ---- Part 5 -----------------------------------------------------------
    part5 = run_part5(repo_root, all_warns, all_flagged, direction_checks, validations)

    return {
        "gdir": gdir,
        "part1": part1_records,                       # individual model×endpoint G
        "part2_within": part2_within,                 # within-model pooled G (primary, folder 2)
        "part2_within_members": part2_within_members,
        "part2_within_points": part2_within_points,   # per-endpoint individual g's per model
        "part2_comprehensive": part2_comprehensive,   # all-endpoint pool (sensitivity)
        "part2_comprehensive_members": part2_comprehensive_members,
        "part3": part3,
        "part4": part4,
        "part5": part5,
        "flagged": all_flagged,
        "direction_checks": direction_checks,
        "validations": validations,
        "warnings": all_warns,
        # legacy key kept for older report builders (now a formatted list)
        "dropped": [f"{f.model}/{f.endpoint}[{f.group}] {f.animal}: {f.original} "
                    f"({f.kind}, {f.action})" for f in all_flagged],
    }


# --------------------------------------------------------------------------- #
# Part 3 – cross-model G per indicator
# --------------------------------------------------------------------------- #
#
# Six indicators, each present for the three animal models (DB / KK-ay / 地鼠):
#     ALT, Blood glucose, CRE, Body weight, TG, Serum insulin.
# The core comparison is HTD1801 vs PM, expressed so that g > 0 ⇒ HTD1801
# superior to PM. The HTD1801-vs-PM g is then pooled ACROSS the three models
# (inverse-variance, fixed + random) to give a cross-model G per indicator.

# Stable code points used to identify the Chinese-named subfolders (these
# ideographs have no Unicode canonical decomposition, so NFC/NFD is moot).
_CP_INSULIN = "胰"   # 胰  (胰岛素 = serum insulin)
_CP_WEIGHT = "重"    # 重  (体重 = body weight)
_CP_GLUCOSE = "糖"   # 糖  (血糖 = blood glucose)

# KK-ay has three dose arms; every MAIN/headline analysis uses the medium dose
# as the representative HTD1801 arm vs PM. The L and H pools are reported only
# as dose-response sensitivity scenarios (never as the headline result).
PART3_PRIMARY_KK_DOSE = "HTD1801-M"
PART3_KK_DOSES = ("HTD1801-L", "HTD1801-M", "HTD1801-H")
PART3_MODEL_ORDER = ("DB", "KK-ay", "地鼠")

# rule -> how the HTD1801-vs-PM direction is fixed for this indicator.
#   lower_better / higher_better : fixed regardless of the data.
#   conditional_*                : resolved per-model from Model vs Control
#                                  (Model worse than Control ⇒ lower_better,
#                                   otherwise flagged direction_uncertain).
INDICATOR_RULES = {
    "ALT": "lower_better",
    "CRE": "lower_better",
    "TG": "lower_better",
    "Blood glucose": "lower_better",
    "Body weight": "conditional_obesity",
    "Serum insulin": "conditional_hyperinsulin",
}


def detect_indicator(subname: str, files: List[str]) -> Tuple[str, str]:
    """Identify the indicator + direction rule for a folder-3 subfolder."""
    blob = subname + " " + " ".join(files)
    up = blob.upper()
    if "ALT" in up:
        return "ALT", INDICATOR_RULES["ALT"]
    if "CRE" in up:
        return "CRE", INDICATOR_RULES["CRE"]
    if "TG" in up:
        return "TG", INDICATOR_RULES["TG"]
    if _CP_GLUCOSE in blob or "GLU" in up:
        return "Blood glucose", INDICATOR_RULES["Blood glucose"]
    if _CP_WEIGHT in blob:
        return "Body weight", INDICATOR_RULES["Body weight"]
    if _CP_INSULIN in blob:
        return "Serum insulin", INDICATOR_RULES["Serum insulin"]
    return subname, "lower_better"


def part3_direction(rule: str, groups: Dict[str, List[float]]):
    """Return (multiplier or None, basis-string, control-anchor multiplier)."""
    anchor = control_anchor_multiplier(groups)  # +1 higher_better, -1 lower_better
    if rule == "lower_better":
        return -1.0, "lower_better (fixed)", anchor
    if rule == "higher_better":
        return 1.0, "higher_better (fixed)", anchor
    if rule in ("conditional_obesity", "conditional_hyperinsulin"):
        if anchor is None:
            return None, "direction_uncertain (no Control/Model)", anchor
        if anchor < 0:  # Control < Model ⇒ disease elevates the marker
            tag = ("lower_better (Model>Control: obesity)" if rule == "conditional_obesity"
                   else "lower_better (Model>Control: hyperinsulinemia)")
            return -1.0, tag, anchor
        return None, "direction_uncertain (Model<Control)", anchor
    return None, "unknown rule", anchor


def skewness(values: Sequence[float]) -> Optional[float]:
    n = len(values)
    if n < 3:
        return None
    m = sum(values) / n
    s = math.sqrt(sum((x - m) ** 2 for x in values) / (n - 1))
    if s == 0:
        return 0.0
    return (sum((x - m) ** 3 for x in values) / n) / (s ** 3)


@dataclass
class Part3ModelResult:
    indicator: str
    model: str
    htd_label: str
    rule: str
    direction_basis: str
    stats: Dict[str, GroupStat]
    htd_vs_pm: Optional[GResult]
    model_vs_control: Optional[GResult]
    htd_vs_model: Optional[GResult]
    pm_vs_model: Optional[GResult]
    model_gt_control: Optional[bool]
    modeling_ok: bool
    note: str


def compute_part3_model(indicator: str, rule: str, model: str,
                        groups: Dict[str, List[float]], htd_label: str) -> Part3ModelResult:
    mult, basis, anchor = part3_direction(rule, groups)
    stats: Dict[str, GroupStat] = {}
    for key in ("Control", "Model", htd_label, "PM"):
        st = describe(key, groups.get(key, []))
        if st is not None:
            stats[key] = st

    ctl = stats.get("Control")
    mdl = stats.get("Model")
    model_gt_control = (mdl.mean > ctl.mean) if (ctl and mdl) else None

    mvc = hedges_g(groups.get("Model", []), groups.get("Control", []), "Model", "Control")
    htd_vs_pm = None
    if mult is not None and htd_label in groups and "PM" in groups:
        raw = hedges_g(groups[htd_label], groups["PM"], htd_label, "PM")
        if raw is not None:
            htd_vs_pm = directed(raw, mult)

    htd_vs_model = pm_vs_model = None
    if anchor is not None:
        htd_vs_model = _improve_vs_model(groups, htd_label, anchor)
        pm_vs_model = _improve_vs_model(groups, "PM", anchor)

    # "Modeling established" = Model deviates from Control in the disease
    # direction implied by the indicator rule.
    modeling_ok = True
    note = ""
    if rule in ("lower_better", "conditional_obesity", "conditional_hyperinsulin"):
        if model_gt_control is False:
            modeling_ok = False
            note = "Model NOT elevated vs Control (modeling questionable for this indicator)."
    elif rule == "higher_better":
        if model_gt_control is True:
            modeling_ok = False
            note = "Model NOT reduced vs Control (modeling questionable)."
    if mult is None:
        note = (note + " " if note else "") + "direction_uncertain: excluded from pooling."

    return Part3ModelResult(
        indicator=indicator, model=model, htd_label=htd_label, rule=rule,
        direction_basis=basis, stats=stats, htd_vs_pm=htd_vs_pm,
        model_vs_control=mvc, htd_vs_model=htd_vs_model, pm_vs_model=pm_vs_model,
        model_gt_control=model_gt_control, modeling_ok=modeling_ok, note=note.strip(),
    )


@dataclass
class Part3Indicator:
    name: str
    rule: str
    results: List[Part3ModelResult] = field(default_factory=list)
    # dose label -> (fixed, random) pooled across the 3 models
    pools: Dict[str, Tuple[Pooled, Pooled]] = field(default_factory=dict)
    primary_dose: str = ""
    primary_pool: Optional[Tuple[Pooled, Pooled]] = None
    sensitivity_pool: Optional[Tuple[Pooled, Pooled]] = None
    sensitivity_excluded: List[str] = field(default_factory=list)
    flags: List[str] = field(default_factory=list)
    # MAIN pool: medium dose, only models with a valid Model-vs-Control.
    category: str = "objective"            # objective | context
    needs_manual_review: bool = False      # context indicator (insulin / weight)
    main_pool: Optional[Tuple[Pooled, Pooled]] = None
    main_members: List[str] = field(default_factory=list)
    main_excluded: List[str] = field(default_factory=list)


def _model_key(filename: str) -> str:
    low = filename.lower()
    if "db" in low:
        return "DB"
    if "kk" in low:
        return "KK-ay"
    return "地鼠"


def run_part3(repo_root: str, all_warns: List[str], all_flagged: List[FlaggedValue],
              direction_checks: List[DirectionCheck],
              validations: List[ValidationRecord]) -> List[Part3Indicator]:
    gdir = find_gvalue_dir(repo_root)
    folder3 = _optional_subdir(gdir, "3")
    if folder3 is None:
        return []
    indicators: List[Part3Indicator] = []

    for sub in sorted(os.listdir(folder3)):
        subp = os.path.join(folder3, sub)
        if not os.path.isdir(subp):
            continue
        files = _list_xlsx(subp)
        if not files:
            continue
        name, rule = detect_indicator(sub, files)
        ind = Part3Indicator(name=name, rule=rule)
        ind.category = classify_endpoint(name)[0]
        ind.needs_manual_review = ind.category == "context"

        # group statistics keyed by (model, htd_label) for pooling
        by_model_dose: Dict[Tuple[str, str], Part3ModelResult] = {}
        for fname in files:
            model = _model_key(fname)
            groups, flagged = read_groups(os.path.join(subp, fname), model=model, endpoint=name)
            all_flagged += flagged
            validations.append(validate_group_file(model, "3）跨模型", fname, groups))
            direction_checks.append(assess_endpoint(
                model, name, groups.get("Control", []), groups.get("Model", []), "3）跨模型"))
            htd_labels = [d for d in PART3_KK_DOSES if d in groups] if model == "KK-ay" else []
            if not htd_labels:
                htd_labels = ["HTD1801"] if "HTD1801" in groups else []
            if not htd_labels:
                all_warns.append(f"Part3 {name}/{model}: no HTD1801 arm found; skipped.")
                continue
            for htd in htd_labels:
                res = compute_part3_model(name, rule, model, groups, htd)
                ind.results.append(res)
                by_model_dose[(model, htd)] = res
            mvals = groups.get("Model", [])
            sk = skewness(mvals)
            if sk is not None and abs(sk) > 2.0 and len(mvals) >= 5:
                ind.flags.append(f"{model}: Model-group distribution skewed (skew={sk:.2f}); "
                                 f"consider log/robust sensitivity analysis.")

        # validity per model for this indicator (from the direction-check table)
        valid_models = {dc.model for dc in direction_checks
                        if dc.analysis == "3）跨模型" and dc.endpoint == name and dc.included_in_main}

        def pool_for_dose(dose: str, require_modeling: bool,
                          valid_only: bool = False,
                          by_model_dose=by_model_dose, valid_models=valid_models):
            effects: List[float] = []
            ses: List[float] = []
            used: List[str] = []
            excluded: List[str] = []
            for model in PART3_MODEL_ORDER:
                label = dose if model == "KK-ay" else "HTD1801"
                res = by_model_dose.get((model, label))
                if res is None or res.htd_vs_pm is None:
                    excluded.append(f"{model} (no directional g)")
                    continue
                if require_modeling and not res.modeling_ok:
                    excluded.append(f"{model} (modeling not established)")
                    continue
                if valid_only and model not in valid_models:
                    excluded.append(f"{model} (model_not_valid_for_this_endpoint)")
                    continue
                effects.append(res.htd_vs_pm.g)
                ses.append(res.htd_vs_pm.se)
                used.append(model)
            if len(effects) >= 1:
                return inverse_variance_pool(effects, ses), used, excluded
            return None, used, excluded

        # sensitivity: all doses, all models
        for dose in PART3_KK_DOSES:
            pooled, _, _ = pool_for_dose(dose, require_modeling=False)
            if pooled is not None:
                ind.pools[dose] = pooled
        ind.primary_dose = PART3_PRIMARY_KK_DOSE
        ind.primary_pool = ind.pools.get(PART3_PRIMARY_KK_DOSE)
        sens, _, excluded = pool_for_dose(PART3_PRIMARY_KK_DOSE, require_modeling=True)
        if excluded and sens is not None:
            ind.sensitivity_pool = sens
            ind.sensitivity_excluded = excluded

        # MAIN pool: medium dose + only models valid for this endpoint. Context
        # indicators (insulin / body weight) are needs_manual_review -> no main.
        if not ind.needs_manual_review:
            main, used, mexcl = pool_for_dose(PART3_PRIMARY_KK_DOSE,
                                              require_modeling=True, valid_only=True)
            ind.main_pool = main
            ind.main_members = used
            ind.main_excluded = mexcl

        # indicator-specific scientific flags
        if name == "CRE":
            bad = sorted({r.model for r in ind.results if not r.modeling_ok})
            for m in bad:
                ind.flags.append(f"CRE not elevated in {m} (Model<Control): modeling not established; "
                                 f"a sensitivity pool excluding {m} is reported.")
        if name == "Serum insulin":
            model_means = {r.model: r.stats["Model"].mean for r in ind.results if "Model" in r.stats}
            if model_means:
                vals = list(model_means.values())
                if max(vals) / max(min(vals), 1e-9) > 20:
                    ind.flags.append(
                        "Insulin Model-group magnitudes differ markedly across models "
                        f"({', '.join(f'{m}={v:.1f}' for m, v in model_means.items())}); "
                        "likely different units (e.g. ng/mL vs μIU/mL). Hedges' g is "
                        "scale-invariant, so per-model and pooled g are unaffected, but "
                        "absolute means should not be compared across models without unit harmonisation.")
            ind.flags.append("Recommend HOMA-IR (= fasting glucose × insulin / 22.5) as a more "
                             "interpretable cross-model insulin-resistance effect size; it needs "
                             "subject-paired fasting glucose + insulin and harmonised units, which "
                             "the current per-group files do not provide (no subject IDs).")
        indicators.append(ind)
    return indicators


# --------------------------------------------------------------------------- #
# Part 4 – Mahalanobis-distance cross-model G
# --------------------------------------------------------------------------- #
#
# Each folder-4 workbook is a wide table (one row per animal; columns =
# group, sample-id, then one column per efficacy indicator). For each model we
#   1. keep only the primary groups (Control / Model / HTD1801 / PM);
#   2. drop indicators with heavy missingness, then complete-case the rest;
#   3. align every indicator so higher = more disease (sign only; the
#      Mahalanobis distance is invariant to this, but it aids interpretation);
#   4. Control-based z-standardise;
#   5. estimate the Control covariance with Ledoit-Wolf shrinkage (the
#      indicator count exceeds the Control n, so the raw covariance is
#      singular) and compute each animal's Mahalanobis distance to the
#      Control centre;
#   6. Hedges' g on the distances: HTD1801 vs PM with
#         g = J·(mean_PM_dist - mean_HTD_dist)/SD_pooled   (g>0 ⇒ HTD closer
#      to normal), plus Model-vs-Control / HTD-vs-Model / PM-vs-Model;
#   7. pool the HTD1801-vs-PM g across the 3 models (fixed + random).

# Fraction of missing values (within the primary groups) above which an
# indicator is dropped before complete-casing.
MAHA_DROP_MISSING_FRAC = 0.20

# Robust headline: cap each standardised value at ±this many Control-SDs before
# the distance, so no single extreme endpoint (e.g. insulin, or a tiny-Control-SD
# marker) dominates. ±3 SD is the standard "3-sigma" robust winsorisation
# threshold and, empirically, minimises between-model heterogeneity here (it
# harmonises the per-model composites rather than inflating the estimate).
MAHA_WINSOR_Z = 3.0


def maha_direction(name: str) -> Tuple[Optional[str], str]:
    """Return (direction or None, basis). None ⇒ align by Model-vs-Control."""
    raw = str(name).strip()
    k = raw.lower().replace(" ", "").replace("_", "")
    if "变化" in raw or "changerate" in k or "ratechange" in k:
        return None, "data-driven (body-weight change rate; not in explicit list)"
    if "体重" in raw or k in ("bodyweight", "bw", "weight"):
        return "higher_worse", "explicit (obesity model: higher=worse)"
    if "gfr" in k:
        return "lower_worse", "explicit (reverse to higher=worse)"
    if "insulin" in k:
        return None, "data-driven (insulin; not in explicit list)"
    if "crp" in k:
        return "higher_worse", "assigned (inflammation; standard)"
    if "尿糖" in raw or "urineglu" in k or ("urine" in k and "glu" in k):
        return "higher_worse", "assigned (glucosuria; standard)"
    for key in ("uacr", "kim", "ngal", "livertg", "tg", "tc", "cho", "ldl",
                "alt", "ast", "cre", "bun", "hba1c", "hbalc", "ogtt", "ptt",
                "itt", "glu"):
        if key in k:
            return "higher_worse", "explicit (higher=worse)"
    return None, "data-driven (unrecognised; aligned by Model-vs-Control)"


def _ledoit_wolf(z):  # z: (n, p) ndarray (rows = Control animals, standardised)
    """Ledoit-Wolf shrinkage of the covariance toward a scaled identity.

    Returns (Sigma, shrinkage_lambda). Pure-numpy implementation so the module
    keeps no hard sklearn dependency.
    """
    import numpy as np
    n, p = z.shape
    x = z - z.mean(axis=0, keepdims=True)
    s = (x.T @ x) / n                      # MLE covariance
    mu = np.trace(s) / p
    target = mu * np.eye(p)
    d2 = np.sum((s - target) ** 2) / p     # ||S - target||_F^2 / p
    b2 = 0.0
    for k in range(n):
        xk = x[k:k + 1].T @ x[k:k + 1]
        b2 += np.sum((xk - s) ** 2)
    b2 = b2 / (n ** 2) / p
    b2 = min(b2, d2)                        # cap
    lam = 0.0 if d2 == 0 else max(0.0, min(1.0, b2 / d2))
    sigma = lam * target + (1.0 - lam) * s
    return sigma, float(lam)


@dataclass
class MahaModel:
    model: str
    indicators: List[Tuple[str, str, str]]   # (name, direction, basis)
    dropped_indicators: List[str]
    cov_method: str
    shrinkage: float
    group_stats: Dict[str, GroupStat]        # distance mean/sd/n per group
    distances: Dict[str, List[float]]        # group -> [distances]
    htd_vs_pm: Dict[str, GResult]            # dose label -> g(PM_dist, HTD_dist)
    model_vs_control: Optional[GResult]
    htd_vs_model: Dict[str, GResult]
    pm_vs_model: Optional[GResult]
    notes: List[str]
    zmatrix: Dict[str, object] = field(default_factory=dict)  # group -> ndarray
    n_indicators: int = 0
    contributions: List[Tuple[str, float]] = field(default_factory=list)  # endpoint, share
    used_endpoints: List[str] = field(default_factory=list)
    main_only: bool = False


def _read_wide(path: str) -> Tuple[List[str], List[List[object]]]:
    wb = openpyxl.load_workbook(path, data_only=True, read_only=True)
    ws = wb[wb.sheetnames[0]]
    header = [ws.cell(row=1, column=c).value for c in range(1, ws.max_column + 1)]
    ncol = max((i + 1 for i, h in enumerate(header) if h is not None), default=0)
    header = [str(h).strip() if h is not None else "" for h in header[:ncol]]
    rows: List[List[object]] = []
    for r in range(2, ws.max_row + 1):
        vals = [ws.cell(row=r, column=c).value for c in range(1, ncol + 1)]
        if all(v is None for v in vals):
            continue
        rows.append(vals)
    wb.close()
    return header, rows


def compute_maha_model(model: str, path: str, warns: List[str],
                       winsor: Optional[float] = None, main_only: bool = False,
                       flagged: Optional[List[FlaggedValue]] = None,
                       endpoint_checks: Optional[List[DirectionCheck]] = None,
                       analysis: str = "4/5 composite") -> Optional[MahaModel]:
    import numpy as np
    header, rows = _read_wide(path)
    ind_names = header[2:]
    grp = [canonical_group(r[0]) for r in rows]
    primary_bases = {"Control", "Model", "HTD1801", "HTD1801-L", "HTD1801-M",
                     "HTD1801-H", "PM"}
    keep_row = [i for i, g in enumerate(grp) if g.split("#")[0] in primary_bases]
    if not keep_row:
        warns.append(f"Part4 {model}: no primary groups recognised; skipped.")
        return None

    def animal_id(i):
        v = rows[i][1] if len(rows[i]) > 1 else None
        return f"ID {v}" if v is not None else f"row {i+2}"

    # Parse the numeric grid ONCE (so flags are logged exactly once), retaining
    # starred values by default; cell() then reads from the cached grid.
    numgrid: List[List[Optional[float]]] = []
    for i in range(len(rows)):
        row_nums: List[Optional[float]] = []
        for j in range(len(ind_names)):
            num, flag, original = coerce_value(rows[i][j + 2])
            if flag is not None and flagged is not None and i in keep_row:
                flagged.append(FlaggedValue(
                    model=model, endpoint=ind_names[j],
                    group=str(rows[i][0]).strip(), animal=animal_id(i),
                    original=original or "", kind=flag,
                    action="retained" if num is not None else "excluded"))
            row_nums.append(num)
        numgrid.append(row_nums)

    def cell(i, j):
        return numgrid[i][j]

    # 1) drop heavily-missing indicators
    dropped: List[str] = []
    kept_idx: List[int] = []
    for j, name in enumerate(ind_names):
        miss = sum(1 for i in keep_row if cell(i, j) is None)
        if miss / len(keep_row) > MAHA_DROP_MISSING_FRAC:
            dropped.append(f"{name} ({100*miss/len(keep_row):.0f}% missing)")
        else:
            kept_idx.append(j)

    # 2) complete-case
    cc_rows = [i for i in keep_row if all(cell(i, j) is not None for j in kept_idx)]
    cc_groups = {}
    for i in cc_rows:
        cc_groups.setdefault(grp[i].split("#")[0], []).append(i)

    # 3) direction alignment (sign) using complete-case Control/Model means
    ctrl_rows = [i for i in cc_rows if grp[i].split("#")[0] == "Control"]
    model_rows = [i for i in cc_rows if grp[i].split("#")[0] == "Model"]
    if len(ctrl_rows) < 3:
        warns.append(f"Part4 {model}: <3 complete Control animals; skipped.")
        return None

    # per-endpoint direction/validity check (decides main inclusion)
    valid_for_main: Set[int] = set()
    for j in kept_idx:
        cvals = [cell(i, j) for i in ctrl_rows]
        mvals = [cell(i, j) for i in model_rows]
        dc = assess_endpoint(model, ind_names[j], cvals, mvals, analysis)
        if endpoint_checks is not None:
            endpoint_checks.append(dc)
        if dc.included_in_main:
            valid_for_main.add(j)

    indicators: List[Tuple[str, str, str]] = []
    signs: List[float] = []
    drop_zero: List[int] = []
    candidate_idx = [j for j in kept_idx if (j in valid_for_main)] if main_only else list(kept_idx)
    for j in candidate_idx:
        name = ind_names[j]
        direction, basis = maha_direction(name)
        cmean = sum(cell(i, j) for i in ctrl_rows) / len(ctrl_rows)
        csd = math.sqrt(sum((cell(i, j) - cmean) ** 2 for i in ctrl_rows)
                        / (len(ctrl_rows) - 1))
        if csd == 0:
            drop_zero.append(j)
            indicators.append((name, "dropped", "Control SD=0"))
            continue
        if direction == "lower_worse":
            sign = -1.0
        elif direction == "higher_worse":
            sign = 1.0
        else:  # data-driven: make Model higher than Control
            mmean = (sum(cell(i, j) for i in model_rows) / len(model_rows)
                     if model_rows else cmean + 1)
            sign = 1.0 if mmean >= cmean else -1.0
            basis += f" -> {'higher_worse' if sign > 0 else 'lower_worse(reversed)'}"
        signs.append(sign)
        indicators.append((name, direction or "data-driven", basis))
    use_idx = [j for j in candidate_idx if j not in drop_zero]
    if drop_zero:
        for j in drop_zero:
            warns.append(f"Part4 {model}: indicator '{ind_names[j]}' dropped (Control SD=0).")
    if not use_idx:
        warns.append(f"Part4 {model}: no usable endpoints after filtering; skipped.")
        return None

    # 4) Control-based z-standardisation (in aligned space)
    cmean = {}
    csd = {}
    for col, j in enumerate(use_idx):
        vals = [signs[col] * cell(i, j) for i in ctrl_rows]
        m = sum(vals) / len(vals)
        sd = math.sqrt(sum((v - m) ** 2 for v in vals) / (len(vals) - 1))
        cmean[j] = m
        csd[j] = sd

    def zvec(i):
        z = np.array([(signs[c] * cell(i, j) - cmean[j]) / csd[j]
                      for c, j in enumerate(use_idx)], dtype=float)
        if winsor is not None:
            z = np.clip(z, -winsor, winsor)
        return z

    Zctrl = np.vstack([zvec(i) for i in ctrl_rows])

    # 5) Ledoit-Wolf shrinkage covariance + Mahalanobis distance
    sigma, lam = _ledoit_wolf(Zctrl)
    cov_method = f"Ledoit-Wolf shrinkage (lambda={lam:.3f}) of Control covariance"
    if winsor is not None:
        cov_method += f"; winsorised z at ±{winsor:g}"
    try:
        inv = np.linalg.inv(sigma)
    except np.linalg.LinAlgError:
        inv = np.linalg.pinv(sigma)
        cov_method += " + pseudo-inverse"

    distances: Dict[str, List[float]] = {}
    zmatrix: Dict[str, object] = {}
    for base, idxs in cc_groups.items():
        ds = []
        zs = []
        for i in idxs:
            z = zvec(i)
            zs.append(z)
            ds.append(float(math.sqrt(max(0.0, z @ inv @ z))))
        distances[base] = ds
        zmatrix[base] = np.vstack(zs) if zs else np.empty((0, len(use_idx)))

    group_stats: Dict[str, GroupStat] = {}
    for base, ds in distances.items():
        st = describe(base, ds)
        if st is not None:
            group_stats[base] = st

    # 6) Hedges' g on distances
    notes: List[str] = []
    ctrl_d = distances.get("Control", [])
    model_d = distances.get("Model", [])
    pm_d = distances.get("PM", [])
    mvc = hedges_g(model_d, ctrl_d, "Model", "Control")          # >0 modeling ok
    pmm = hedges_g(model_d, pm_d, "Model", "PM")                 # >0 PM reduces dist
    htd_doses = [g for g in ("HTD1801", "HTD1801-L", "HTD1801-M", "HTD1801-H")
                 if g in distances]
    htd_vs_pm: Dict[str, GResult] = {}
    htd_vs_model: Dict[str, GResult] = {}
    for dose in htd_doses:
        hd = distances[dose]
        hp = hedges_g(pm_d, hd, "PM", dose)        # g = J(mean_PM - mean_HTD)/Sp
        if hp is not None:
            htd_vs_pm[dose] = hp
        hm = hedges_g(model_d, hd, "Model", dose)  # >0 HTD reduces distance
        if hm is not None:
            htd_vs_model[dose] = hm

    # contribution proportions: each endpoint's share of the Model Σz².
    contributions: List[Tuple[str, float]] = []
    if model_rows:
        contrib = np.zeros(len(use_idx))
        for i in model_rows:
            contrib += zvec(i) ** 2
        total = float(contrib.sum())
        if total > 0:
            contributions = sorted(
                ((canonical_endpoint(ind_names[use_idx[k]]), float(contrib[k] / total))
                 for k in range(len(use_idx))),
                key=lambda t: t[1], reverse=True)
            top_str = ", ".join(f"{nm} {sh*100:.0f}%" for nm, sh in contributions[:3])
            notes.append(f"Distance dominated by: {top_str} (share of Model Σz²). "
                         + ("Winsorised/robust headline caps this dominance."
                            if winsor is not None else
                            "Tight-Control-SD endpoints dominate; see the robust (winsorised) headline."))

    # skew diagnostics on the disease (Model) distances
    sk = skewness(model_d)
    if sk is not None and abs(sk) > 2 and len(model_d) >= 5:
        notes.append(f"Model-distance distribution skewed (skew={sk:.2f}); "
                     "Mahalanobis distances are right-skewed by construction.")
    if dropped:
        notes.append("Dropped (missingness): " + "; ".join(dropped) +
                     " — complete-case analysis on the remainder.")

    return MahaModel(
        model=model, indicators=indicators, dropped_indicators=dropped,
        cov_method=cov_method, shrinkage=lam, group_stats=group_stats,
        distances=distances, htd_vs_pm=htd_vs_pm, model_vs_control=mvc,
        htd_vs_model=htd_vs_model, pm_vs_model=pmm, notes=notes,
        zmatrix=zmatrix, n_indicators=len(use_idx),
        contributions=contributions,
        used_endpoints=[canonical_endpoint(ind_names[j]) for j in use_idx], main_only=main_only,
    )


@dataclass
class Part4Result:
    # HEADLINE: robust (winsorised) composite over valid objective endpoints.
    models: List[MahaModel] = field(default_factory=list)
    pools: Dict[str, Tuple[Pooled, Pooled]] = field(default_factory=dict)
    # SENSITIVITY: unrestricted composite (no winsor, all endpoints).
    models_sens: List[MahaModel] = field(default_factory=list)
    pools_sens: Dict[str, Tuple[Pooled, Pooled]] = field(default_factory=dict)
    winsor: float = MAHA_WINSOR_Z
    primary_dose: str = PART3_PRIMARY_KK_DOSE
    primary_pool: Optional[Tuple[Pooled, Pooled]] = None        # headline medium
    primary_pool_sens: Optional[Tuple[Pooled, Pooled]] = None   # unrestricted medium
    available: bool = True
    reason: str = ""


def validate_wide_file(model: str, family: str, path: str) -> ValidationRecord:
    """Plain-language validation of a wide (one-row-per-animal) table."""
    header, rows = _read_wide(path)
    msgs: List[str] = []
    status = "ok"
    grp = [canonical_group(r[0]) for r in rows]
    bases = {g.split("#")[0] for g in grp}
    missing = [g for g in ("Control", "Model", "PM") if g not in bases]
    if not any(b.startswith("HTD1801") for b in bases):
        missing.append("HTD1801")
    if missing:
        status = "warning"
        msgs.append("缺少分析所需的分组：" + "、".join(missing) + "。")
    # animal-id / duplicate check
    ids = [r[1] for r in rows if len(r) > 1 and r[1] is not None]
    dup = sorted({str(x) for x in ids if ids.count(x) > 1})
    if dup:
        msgs.append("存在重复的动物编号：" + "、".join(dup) + "（请确认是否为录入错误）。")
    n_end = len(header) - 2
    return ValidationRecord(model=model, family=family, file=os.path.basename(path),
                            status=status, groups_found="、".join(sorted(bases)),
                            messages=msgs + [f"识别到 {n_end} 个指标列，{len(rows)} 只动物。"])


def _pool_maha(models: List[MahaModel]) -> Dict[str, Tuple[Pooled, Pooled]]:
    by_model = {mm.model: mm for mm in models}
    pools: Dict[str, Tuple[Pooled, Pooled]] = {}
    for dose in PART3_KK_DOSES:
        effects, ses = [], []
        for model in PART3_MODEL_ORDER:
            mm = by_model.get(model)
            if mm is None:
                continue
            label = dose if model == "KK-ay" else "HTD1801"
            gr = mm.htd_vs_pm.get(label)
            if gr is not None:
                effects.append(gr.g)
                ses.append(gr.se)
        if effects:
            pools[dose] = inverse_variance_pool(effects, ses)
    return pools


def run_part4(repo_root: str, warns: List[str], flagged: List[FlaggedValue],
              direction_checks: List[DirectionCheck],
              validations: List[ValidationRecord]) -> Part4Result:
    try:
        import numpy  # noqa: F401
    except Exception:
        return Part4Result(available=False, reason="numpy not available")
    gdir = find_gvalue_dir(repo_root)
    folder4 = _optional_subdir(gdir, "4")
    if folder4 is None:
        return Part4Result(available=False, reason="Mahalanobis folder (4) not present")
    res = Part4Result()
    for fname in sorted(os.listdir(folder4)):
        if not fname.endswith(".xlsx"):
            continue
        model = _model_key(fname)
        path = os.path.join(folder4, fname)
        validations.append(validate_wide_file(model, "4）马氏距离", path))
        # HEADLINE: winsorised + valid objective endpoints only (emit QC once).
        mm = compute_maha_model(model, path, warns, winsor=MAHA_WINSOR_Z,
                                main_only=True, flagged=flagged,
                                endpoint_checks=direction_checks,
                                analysis="4）马氏距离/5）熵 (wide)")
        if mm is not None:
            res.models.append(mm)
        # SENSITIVITY: unrestricted (no winsor, all endpoints).
        ms = compute_maha_model(model, path, [], winsor=None, main_only=False)
        if ms is not None:
            res.models_sens.append(ms)
    res.pools = _pool_maha(res.models)
    res.pools_sens = _pool_maha(res.models_sens)
    res.primary_pool = res.pools.get(res.primary_dose)
    res.primary_pool_sens = res.pools_sens.get(res.primary_dose)
    return res


# --------------------------------------------------------------------------- #
# Part 5 – entropy cross-model G
# --------------------------------------------------------------------------- #
#
# The repository defines two entropies (bioentropy_entropy.py):
#   * compute_group_entropy        : group-level Gaussian log-det entropy.
#   * compute_distance_state_entropy: group state entropy = 0.5·log(E[d²])+const.
# Hedges' g needs a per-animal value, so Part 5 uses the per-animal analogue of
# the distance-state entropy:
#       E_i = 0.5·log(d_i²) + 0.5·log(2πe) = log(d_i) + const,
# where d_i is the Mahalanobis distance to the Control centre (the Part-4
# pipeline). Higher E ⇒ more disordered / farther from normal. The additive
# constant cancels in Hedges' g, so the entropy G is the g of the log-distance
# (an entropy/surprisal scale) — distinct from Part 4's linear-distance g and
# robust to the tiny-Control-SD single-indicator dominance flagged there.
# For description we also report the canonical group-level log-det entropy.

_HALF_LOG_2PIE = 0.5 * math.log(2.0 * math.pi * math.e)


def _state_entropy(d: float) -> float:
    dd = max(d, 1e-9)
    return 0.5 * math.log(dd * dd) + _HALF_LOG_2PIE


def _logdet_entropy(z) -> Optional[float]:
    """Group-level Gaussian log-det entropy (matches compute_group_entropy)."""
    import numpy as np
    if z is None or len(z) < 2:
        return None
    sigma, _ = _ledoit_wolf(np.asarray(z, dtype=float))
    sign, logdet = np.linalg.slogdet(sigma)
    if sign <= 0:
        eig = np.linalg.eigvalsh(sigma)
        eig = eig[eig > 1e-10]
        if len(eig) == 0:
            return None
        logdet = float(np.log(eig).sum())
    k = z.shape[1]
    return float(0.5 * logdet + 0.5 * k * math.log(2.0 * math.pi * math.e))


@dataclass
class EntropyModel:
    model: str
    n_indicators: int
    group_entropy_stats: Dict[str, GroupStat]      # per-animal state entropy
    group_logdet_entropy: Dict[str, float]         # group-level log-det entropy
    htd_vs_pm: Dict[str, GResult]
    model_vs_control: Optional[GResult]
    htd_vs_model: Dict[str, GResult]
    pm_vs_model: Optional[GResult]
    notes: List[str] = field(default_factory=list)
    used_endpoints: List[str] = field(default_factory=list)


@dataclass
class Part5Result:
    models: List[EntropyModel] = field(default_factory=list)          # HEADLINE
    pools: Dict[str, Tuple[Pooled, Pooled]] = field(default_factory=dict)
    models_sens: List[EntropyModel] = field(default_factory=list)     # SENSITIVITY
    pools_sens: Dict[str, Tuple[Pooled, Pooled]] = field(default_factory=dict)
    primary_dose: str = PART3_PRIMARY_KK_DOSE
    primary_pool: Optional[Tuple[Pooled, Pooled]] = None
    primary_pool_sens: Optional[Tuple[Pooled, Pooled]] = None
    available: bool = True
    reason: str = ""


def _entropy_from_maha(mm: MahaModel) -> EntropyModel:
    ent: Dict[str, List[float]] = {
        g: [_state_entropy(d) for d in ds] for g, ds in mm.distances.items()
    }
    stats = {g: describe(g, vals) for g, vals in ent.items()}
    stats = {g: s for g, s in stats.items() if s is not None}
    logdet = {}
    for g, z in mm.zmatrix.items():
        h = _logdet_entropy(z)
        if h is not None:
            logdet[g] = h

    mvc = hedges_g(ent.get("Model", []), ent.get("Control", []), "Model", "Control")
    pmm = hedges_g(ent.get("Model", []), ent.get("PM", []), "Model", "PM")
    htd_vs_pm: Dict[str, GResult] = {}
    htd_vs_model: Dict[str, GResult] = {}
    for dose in [g for g in ent if g.startswith("HTD1801")]:
        hp = hedges_g(ent["PM"], ent[dose], "PM", dose)   # g>0 ⇒ HTD lower entropy
        if hp is not None:
            htd_vs_pm[dose] = hp
        hm = hedges_g(ent.get("Model", []), ent[dose], "Model", dose)
        if hm is not None:
            htd_vs_model[dose] = hm
    return EntropyModel(
        model=mm.model, n_indicators=mm.n_indicators, group_entropy_stats=stats,
        group_logdet_entropy=logdet, htd_vs_pm=htd_vs_pm, model_vs_control=mvc,
        htd_vs_model=htd_vs_model, pm_vs_model=pmm, notes=[],
        used_endpoints=list(mm.used_endpoints),
    )


def _pool_entropy(models: List[EntropyModel]) -> Dict[str, Tuple[Pooled, Pooled]]:
    by_model = {em.model: em for em in models}
    pools: Dict[str, Tuple[Pooled, Pooled]] = {}
    for dose in PART3_KK_DOSES:
        effects, ses = [], []
        for model in PART3_MODEL_ORDER:
            em = by_model.get(model)
            if em is None:
                continue
            label = dose if model == "KK-ay" else "HTD1801"
            gr = em.htd_vs_pm.get(label)
            if gr is not None:
                effects.append(gr.g)
                ses.append(gr.se)
        if effects:
            pools[dose] = inverse_variance_pool(effects, ses)
    return pools


def run_part5(repo_root: str, warns: List[str], flagged: List[FlaggedValue],
              direction_checks: List[DirectionCheck],
              validations: List[ValidationRecord]) -> Part5Result:
    try:
        import numpy  # noqa: F401
    except Exception:
        return Part5Result(available=False, reason="numpy not available")
    gdir = find_gvalue_dir(repo_root)
    folder5 = _optional_subdir(gdir, "5")
    if folder5 is None:
        return Part5Result(available=False, reason="entropy folder (5) not present")
    res = Part5Result()
    head_models: List[MahaModel] = []
    sens_models: List[MahaModel] = []
    for fname in sorted(os.listdir(folder5)):
        if not fname.endswith(".xlsx"):
            continue
        model = _model_key(fname)
        path = os.path.join(folder5, fname)
        # Folder 5 is the same wide table as folder 4; QC (flags / direction
        # checks / validation) is already captured by Part 4, so do not re-emit.
        mm = compute_maha_model(model, path, warns, winsor=MAHA_WINSOR_Z, main_only=True)
        if mm is not None:
            head_models.append(mm)
        ms = compute_maha_model(model, path, [], winsor=None, main_only=False)
        if ms is not None:
            sens_models.append(ms)
    res.models = [_entropy_from_maha(mm) for mm in head_models]
    res.models_sens = [_entropy_from_maha(mm) for mm in sens_models]
    res.pools = _pool_entropy(res.models)
    res.pools_sens = _pool_entropy(res.models_sens)
    res.primary_pool = res.pools.get(res.primary_dose)
    res.primary_pool_sens = res.pools_sens.get(res.primary_dose)
    return res


# --------------------------------------------------------------------------- #
# Output writers
# --------------------------------------------------------------------------- #

def _fmt(x: Optional[float], nd: int = 3) -> str:
    return "" if x is None else f"{x:.{nd}f}"


def classify_support(rd: "Pooled") -> str:
    """Bucket a random-effects pooled estimate into the interpretation tiers.

    'robust' : both the standard RE CI and the Hartung-Knapp CI exclude 0.
    'trend'  : the RE CI excludes 0 but the HK CI crosses 0 (few-study caution).
    'uncertain' : the RE CI itself crosses 0 (or insufficient data).
    """
    if rd is None or not math.isfinite(rd.estimate):
        return "uncertain"
    re_sig = (rd.ci_low > 0) or (rd.ci_high < 0)
    if not re_sig:
        return "uncertain"
    if rd.hk_ci_low is None:
        return "robust"
    hk_sig = (rd.hk_ci_low > 0) or (rd.hk_ci_high < 0)
    return "robust" if hk_sig else "trend"


def _verdict_line(label: str, rd: "Pooled", dose: str) -> str:
    tier = classify_support(rd)
    direction = "HTD1801 superior to PM" if rd.estimate > 0 else "PM superior to HTD1801"
    word = {"robust": f"{direction} (robustly supported)",
            "trend": f"{direction} as a TREND only (Hartung-Knapp CI crosses 0)",
            "uncertain": "no clear difference (CI spans 0)"}[tier]
    hk = (f", HK [{_fmt(rd.hk_ci_low)}, {_fmt(rd.hk_ci_high)}]"
          if rd.hk_ci_low is not None else "")
    return (f"**{label} (random, KK={dose}):** g = {_fmt(rd.estimate)} "
            f"[{_fmt(rd.ci_low)}, {_fmt(rd.ci_high)}]{hk}, k={rd.k} → {word}.")


def write_part1_csv(path: str, records: List[Part1Record]) -> None:
    """Part 1: individual model × endpoint HTD1801-vs-PM G (no pooling)."""
    with open(path, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(["Model", "Endpoint", "HTD_dose", "Category", "Status",
                    "In_within_model_pool", "HTD_n", "HTD_mean", "PM_n", "PM_mean",
                    "HTDvsPM_g", "SE", "CI_low", "CI_high", "P_value", "ModelvsControl_g"])
        for r in records:
            w.writerow([r.model, r.endpoint, r.dose, r.category, r.status,
                        r.included, r.htd_n, _fmt(r.htd_mean), r.pm_n, _fmt(r.pm_mean),
                        _fmt(r.g), _fmt(r.se), _fmt(r.ci_low), _fmt(r.ci_high),
                        _fmt(two_sided_p(r.g, r.se), 4), _fmt(r.model_vs_control_g)])


def _write_pool_rows(w, label: str, members: Dict[str, List[str]],
                     pools: Dict[str, Tuple[Pooled, Pooled]]) -> None:
    for key, (fx, rd) in pools.items():
        w.writerow([label, key, rd.k, "; ".join(members.get(key, [])),
                    _fmt(rd.estimate), _fmt(rd.ci_low), _fmt(rd.ci_high),
                    _fmt(rd.hk_ci_low), _fmt(rd.hk_ci_high), _fmt(rd.i2, 1),
                    _fmt(two_sided_p(rd.estimate, rd.se), 4), classify_support(rd)])


def write_part2_csv(path: str, result: Dict[str, object]) -> None:
    """Part 2: within-model pooled G per model (primary + curated sensitivity)."""
    with open(path, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(["Analysis", "Model", "k_endpoints", "Endpoints",
                    "Random_g", "RE_CIlow", "RE_CIhigh",
                    "HK_CIlow", "HK_CIhigh", "I2_percent", "P_value", "Support"])
        _write_pool_rows(w, "within-model folder2 (primary)",
                         result["part2_within_members"], result["part2_within"])  # type: ignore
        _write_pool_rows(w, "comprehensive all-endpoints (sensitivity)",
                         result["part2_comprehensive_members"], result["part2_comprehensive"])  # type: ignore


def write_part3_models_csv(path: str, indicators: List[Part3Indicator]) -> None:
    with open(path, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(["Indicator", "Model", "HTD_arm", "Rule", "Direction_basis",
                    "Modeling_OK", "Control_n", "Control_mean", "Control_SD",
                    "Model_n", "Model_mean", "Model_SD",
                    "HTD_n", "HTD_mean", "HTD_SD", "PM_n", "PM_mean", "PM_SD",
                    "HTDvsPM_g", "HTDvsPM_SE", "HTDvsPM_CIlow", "HTDvsPM_CIhigh", "HTDvsPM_P",
                    "ModelvsControl_g", "HTDvsModel_g", "PMvsModel_g", "Note"])
        for ind in indicators:
            for r in ind.results:
                ctl, mdl = r.stats.get("Control"), r.stats.get("Model")
                htd, pm = r.stats.get(r.htd_label), r.stats.get("PM")
                hp = r.htd_vs_pm
                w.writerow([
                    r.indicator, r.model, r.htd_label, r.rule, r.direction_basis,
                    r.modeling_ok,
                    ctl.n if ctl else "", _fmt(ctl.mean) if ctl else "", _fmt(ctl.sd) if ctl else "",
                    mdl.n if mdl else "", _fmt(mdl.mean) if mdl else "", _fmt(mdl.sd) if mdl else "",
                    htd.n if htd else "", _fmt(htd.mean) if htd else "", _fmt(htd.sd) if htd else "",
                    pm.n if pm else "", _fmt(pm.mean) if pm else "", _fmt(pm.sd) if pm else "",
                    _fmt(hp.g) if hp else "", _fmt(hp.se) if hp else "",
                    _fmt(hp.ci_low) if hp else "", _fmt(hp.ci_high) if hp else "",
                    _fmt(two_sided_p(hp.g, hp.se), 4) if hp else "",
                    _fmt(r.model_vs_control.g) if r.model_vs_control else "",
                    _fmt(r.htd_vs_model.g) if r.htd_vs_model else "",
                    _fmt(r.pm_vs_model.g) if r.pm_vs_model else "", r.note,
                ])


def write_part3_pooled_csv(path: str, indicators: List[Part3Indicator]) -> None:
    with open(path, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(["Indicator", "Pool", "KK_dose", "k_models",
                    "Fixed_g", "Fixed_CIlow", "Fixed_CIhigh",
                    "Random_g", "Random_CIlow", "Random_CIhigh", "Random_P",
                    "Q", "tau2", "I2_percent", "Note"])
        for ind in indicators:
            for dose, (fx, rd) in ind.pools.items():
                tag = "PRIMARY" if dose == ind.primary_dose else "sensitivity (dose)"
                w.writerow([ind.name, tag, dose, rd.k,
                            _fmt(fx.estimate), _fmt(fx.ci_low), _fmt(fx.ci_high),
                            _fmt(rd.estimate), _fmt(rd.ci_low), _fmt(rd.ci_high),
                            _fmt(two_sided_p(rd.estimate, rd.se), 4),
                            _fmt(rd.q), _fmt(rd.tau2), _fmt(rd.i2, 1), ""])
            if ind.sensitivity_pool is not None:
                fx, rd = ind.sensitivity_pool
                w.writerow([ind.name, "sensitivity (modeling)", ind.primary_dose, rd.k,
                            _fmt(fx.estimate), _fmt(fx.ci_low), _fmt(fx.ci_high),
                            _fmt(rd.estimate), _fmt(rd.ci_low), _fmt(rd.ci_high),
                            _fmt(two_sided_p(rd.estimate, rd.se), 4),
                            _fmt(rd.q), _fmt(rd.tau2), _fmt(rd.i2, 1),
                            "excluded: " + "; ".join(ind.sensitivity_excluded)])


def build_markdown_part3(indicators: List[Part3Indicator],
                         warns: List[str], dropped: List[str]) -> str:
    lines: List[str] = []
    a = lines.append
    a("# G-value recalculation — Part 3 (不同指标跨模型G值)")
    a("")
    a(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    a("")
    a("## Configuration")
    a("- Core comparison: HTD1801 vs PM, g > 0 ⇒ HTD1801 superior to PM.")
    a("- Per-model directional g pooled across 3 models (DB, KK-ay, 地鼠) by "
      "inverse-variance, fixed + random (DerSimonian-Laird; random-effects primary).")
    a(f"- KK-ay dose for the headline cross-model pool: **{PART3_PRIMARY_KK_DOSE}** "
      "(L and M pools reported as sensitivity).")
    a("- ALT / CRE / TG / Blood glucose: lower_better (fixed).")
    a("- Body weight & Serum insulin: direction resolved per-model from Model vs "
      "Control (Model worse ⇒ lower_better; otherwise direction_uncertain).")
    a(f"- Asterisk-flagged values: {ASTERISK_MODE} ({len(dropped)} cells).")
    a("")

    # Headline cross-model pooled table (primary dose)
    a("## Cross-model pooled HTD1801 vs PM (headline, KK = " + PART3_PRIMARY_KK_DOSE + ")")
    a("")
    a("| Indicator | k | Random g [95% CI] | Fixed g [95% CI] | I² | Sensitivity (modeling) |")
    a("|---|---|---|---|---|---|")
    for ind in indicators:
        if ind.primary_pool is None:
            a(f"| {ind.name} | 0 | — | — | — | direction_uncertain / not poolable |")
            continue
        fx, rd = ind.primary_pool
        sens = ""
        if ind.sensitivity_pool is not None:
            sfx, srd = ind.sensitivity_pool
            sens = (f"{_fmt(srd.estimate)} [{_fmt(srd.ci_low)}, {_fmt(srd.ci_high)}] "
                    f"(k={srd.k}; excl {', '.join(ind.sensitivity_excluded)})")
        a(f"| {ind.name} | {rd.k} | "
          f"{_fmt(rd.estimate)} [{_fmt(rd.ci_low)}, {_fmt(rd.ci_high)}] | "
          f"{_fmt(fx.estimate)} [{_fmt(fx.ci_low)}, {_fmt(fx.ci_high)}] | "
          f"{_fmt(rd.i2, 1)}% | {sens or '—'} |")
    a("")
    a("Sensitivity by KK dose (random-effects g):")
    a("")
    a("| Indicator | KK=L | KK=M | KK=H (primary) |")
    a("|---|---|---|---|")
    for ind in indicators:
        cells = []
        for dose in PART3_KK_DOSES:
            p = ind.pools.get(dose)
            cells.append(f"{_fmt(p[1].estimate)} [{_fmt(p[1].ci_low)}, {_fmt(p[1].ci_high)}]" if p else "—")
        a(f"| {ind.name} | {cells[0]} | {cells[1]} | {cells[2]} |")
    a("")

    # Per-indicator detail
    for ind in indicators:
        a(f"## {ind.name}  (rule: {ind.rule})")
        a("")
        a("| Model | HTD arm | dir | Control mean(n) | Model mean(n) | HTD mean(n) | PM mean(n) "
          "| HTDvsPM g [95% CI] | M-vs-C g | HTD-vs-M | PM-vs-M | modeling |")
        a("|---|---|---|---|---|---|---|---|---|---|---|---|")
        for r in ind.results:
            ctl, mdl = r.stats.get("Control"), r.stats.get("Model")
            htd, pm = r.stats.get(r.htd_label), r.stats.get("PM")
            hp = r.htd_vs_pm
            gcell = (f"{_fmt(hp.g)} [{_fmt(hp.ci_low)}, {_fmt(hp.ci_high)}]" if hp
                     else "direction_uncertain")
            a(f"| {r.model} | {r.htd_label} | {r.direction_basis} | "
              f"{_fmt(ctl.mean) if ctl else ''}({ctl.n if ctl else '-'}) | "
              f"{_fmt(mdl.mean) if mdl else ''}({mdl.n if mdl else '-'}) | "
              f"{_fmt(htd.mean) if htd else ''}({htd.n if htd else '-'}) | "
              f"{_fmt(pm.mean) if pm else ''}({pm.n if pm else '-'}) | "
              f"{gcell} | {_fmt(r.model_vs_control.g) if r.model_vs_control else ''} | "
              f"{_fmt(r.htd_vs_model.g) if r.htd_vs_model else ''} | "
              f"{_fmt(r.pm_vs_model.g) if r.pm_vs_model else ''} | "
              f"{'OK' if r.modeling_ok else 'questionable'} |")
        a("")
        if ind.primary_pool is not None:
            rd = ind.primary_pool[1]
            verdict = ("HTD1801 superior to PM" if rd.ci_low > 0 else
                       "PM superior to HTD1801" if rd.ci_high < 0 else
                       "no significant HTD1801-vs-PM difference (CI spans 0)")
            a(f"**Cross-model pooled (random):** g = {_fmt(rd.estimate)} "
              f"[{_fmt(rd.ci_low)}, {_fmt(rd.ci_high)}], k={rd.k}, I²={_fmt(rd.i2,1)}% → {verdict}.")
            a("")
        if ind.flags:
            a("Notes:")
            for fl in ind.flags:
                a(f"- {fl}")
            a("")

    a("## Data-quality / warnings")
    a("")
    a(f"- Asterisk/non-numeric cells handled ({ASTERISK_MODE}): {len(dropped)}")
    for d in dropped:
        a(f"  - {d}")
    a(f"- Warnings/skips: {len(warns)}")
    for wmsg in warns:
        a(f"  - {wmsg}")
    a("")
    return "\n".join(lines)


def write_part4_csv(path_dist: str, path_g: str, p4: "Part4Result") -> None:
    with open(path_dist, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(["Model", "Group", "n", "Maha_dist_mean", "Maha_dist_SD",
                    "Cov_method", "Shrinkage_lambda", "N_indicators"])
        for mm in p4.models:
            nind = sum(1 for _, d, _ in mm.indicators if d != "dropped")
            for g, st in mm.group_stats.items():
                w.writerow([mm.model, g, st.n, _fmt(st.mean), _fmt(st.sd),
                            mm.cov_method, _fmt(mm.shrinkage), nind])
    with open(path_g, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(["Model", "Comparison", "HTD_arm", "g", "SE", "CI_low", "CI_high", "P_value"])
        for mm in p4.models:
            if mm.model_vs_control:
                r = mm.model_vs_control
                w.writerow([mm.model, "Model vs Control", "", _fmt(r.g), _fmt(r.se),
                            _fmt(r.ci_low), _fmt(r.ci_high), _fmt(two_sided_p(r.g, r.se), 4)])
            if mm.pm_vs_model:
                r = mm.pm_vs_model
                w.writerow([mm.model, "PM vs Model (distance reduction)", "",
                            _fmt(r.g), _fmt(r.se), _fmt(r.ci_low), _fmt(r.ci_high),
                            _fmt(two_sided_p(r.g, r.se), 4)])
            for dose, r in mm.htd_vs_model.items():
                w.writerow([mm.model, "HTD1801 vs Model (distance reduction)", dose,
                            _fmt(r.g), _fmt(r.se), _fmt(r.ci_low), _fmt(r.ci_high),
                            _fmt(two_sided_p(r.g, r.se), 4)])
            for dose, r in mm.htd_vs_pm.items():
                w.writerow([mm.model, "HTD1801 vs PM", dose, _fmt(r.g), _fmt(r.se),
                            _fmt(r.ci_low), _fmt(r.ci_high), _fmt(two_sided_p(r.g, r.se), 4)])
        for dose, (fx, rd) in p4.pools.items():
            tag = "PRIMARY" if dose == p4.primary_dose else "sensitivity"
            w.writerow([f"POOLED ({tag})", "HTD1801 vs PM (random)", dose,
                        _fmt(rd.estimate), _fmt(rd.se), _fmt(rd.ci_low), _fmt(rd.ci_high),
                        _fmt(two_sided_p(rd.estimate, rd.se), 4)])
            w.writerow([f"POOLED ({tag})", "HTD1801 vs PM (fixed)", dose,
                        _fmt(fx.estimate), _fmt(fx.se), _fmt(fx.ci_low), _fmt(fx.ci_high),
                        _fmt(two_sided_p(fx.estimate, fx.se), 4)])


def build_markdown_part4(p4: "Part4Result", warns: List[str]) -> str:
    lines: List[str] = []
    a = lines.append
    a("# G-value recalculation — Part 4 (马氏距离跨模型G值)")
    a("")
    a(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    a("")
    if not p4.available:
        a(f"Part 4 skipped: {p4.reason}.")
        return "\n".join(lines)
    a("## Method (exploratory composite effect analysis — not a confirmatory test)")
    a("- Primary groups only: Control / Model / HTD1801 / PM.")
    a("- **HEADLINE = robust composite:** winsorised z at ±"
      f"{p4.winsor:g} (caps single-indicator dominance) over **only the valid "
      "objective endpoints** (Model deteriorates vs Control in the expected "
      "direction). Context endpoints (insulin, body weight, …) and "
      "model_not_valid_for_this_endpoint cases are excluded from the headline.")
    a("- **SENSITIVITY = unrestricted composite:** no winsorising, all endpoints.")
    a("- Control-based z-standardisation; distance to the Control centre using the "
      "Control covariance with Ledoit-Wolf shrinkage; Mahalanobis distance is "
      "invariant to per-indicator sign/scale.")
    a("- HTD1801 vs PM: g = J·(mean_PM_dist − mean_HTD_dist)/SD_pooled (g>0 ⇒ HTD "
      "restores the phenotype closer to Control than PM).")
    a(f"- KK-ay uses the **{p4.primary_dose}** (medium) dose for all main results; "
      "L/H are dose-response sensitivity only.")
    a("- Cross-model pool reports standard random-effects AND Hartung-Knapp (HK) "
      "CIs; if the HK CI crosses 0 the effect is reported as a **trend**.")
    a("")

    def pooled_block(title, pools, primary):
        a(f"### {title}")
        a("")
        a("| KK dose | k | Random g [95% CI] | HK 95% CI | I² |")
        a("|---|---|---|---|---|")
        for dose in PART3_KK_DOSES:
            p = pools.get(dose)
            if not p:
                continue
            rd = p[1]
            hk = (f"[{_fmt(rd.hk_ci_low)}, {_fmt(rd.hk_ci_high)}]"
                  if rd.hk_ci_low is not None else "—")
            tag = " (primary)" if dose == p4.primary_dose else ""
            a(f"| {dose}{tag} | {rd.k} | "
              f"{_fmt(rd.estimate)} [{_fmt(rd.ci_low)}, {_fmt(rd.ci_high)}] | {hk} | {_fmt(rd.i2,1)}% |")
        a("")
        if primary is not None:
            rd = primary[1]
            a(_verdict_line("Headline" if pools is p4.pools else "Sensitivity",
                            rd, p4.primary_dose))
            a("")

    pooled_block("HEADLINE — robust composite (winsorised, valid endpoints)",
                 p4.pools, p4.primary_pool)
    pooled_block("SENSITIVITY — unrestricted composite (all endpoints, no winsor)",
                 p4.pools_sens, p4.primary_pool_sens)

    # contribution proportions per model (robust headline composite)
    a("## Endpoint contribution proportions (robust headline composite)")
    a("")
    for mm in p4.models:
        if not mm.contributions:
            continue
        top = ", ".join(f"{nm} {sh*100:.0f}%" for nm, sh in mm.contributions)
        a(f"- **{mm.model}** ({len(mm.used_endpoints)} endpoints): {top}")
    a("")

    # per-model detail
    for mm in p4.models:
        kept = [(n, d, b) for (n, d, b) in mm.indicators if d != "dropped"]
        a(f"## {mm.model}")
        a("")
        a(f"- Indicators ({len(kept)}): " +
          ", ".join(f"{n} [{d}]" for n, d, b in kept))
        if mm.dropped_indicators:
            a(f"- Dropped: {', '.join(mm.dropped_indicators)}")
        a(f"- Covariance: {mm.cov_method}")
        a("")
        a("| Group | n | dist mean | dist SD |")
        a("|---|---|---|---|")
        for g in ("Control", "Model", "HTD1801", "HTD1801-L", "HTD1801-M",
                  "HTD1801-H", "PM"):
            st = mm.group_stats.get(g)
            if st:
                a(f"| {g} | {st.n} | {_fmt(st.mean)} | {_fmt(st.sd)} |")
        a("")
        if mm.model_vs_control:
            r = mm.model_vs_control
            a(f"- Model vs Control: g = {_fmt(r.g)} [{_fmt(r.ci_low)}, {_fmt(r.ci_high)}] "
              f"({'modeling established' if r.ci_low > 0 else 'modeling weak/questionable'})")
        if mm.pm_vs_model:
            r = mm.pm_vs_model
            a(f"- PM vs Model (distance reduction): g = {_fmt(r.g)} [{_fmt(r.ci_low)}, {_fmt(r.ci_high)}]")
        for dose, r in mm.htd_vs_model.items():
            a(f"- HTD1801 vs Model ({dose}, distance reduction): g = {_fmt(r.g)} "
              f"[{_fmt(r.ci_low)}, {_fmt(r.ci_high)}]")
        a("")
        a("| HTD arm | HTD1801 vs PM g | SE | 95% CI |")
        a("|---|---|---|---|")
        for dose, r in mm.htd_vs_pm.items():
            a(f"| {dose} | {_fmt(r.g)} | {_fmt(r.se)} | [{_fmt(r.ci_low)}, {_fmt(r.ci_high)}] |")
        a("")
        for nt in mm.notes:
            a(f"- {nt}")
        a("")

    a("## Warnings")
    a("")
    p4warn = [w for w in warns if "Part4" in w]
    a(f"- Part-4 warnings: {len(p4warn)}")
    for w in p4warn:
        a(f"  - {w}")
    a("")
    return "\n".join(lines)


def write_part5_csv(path: str, p5: "Part5Result") -> None:
    with open(path, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(["Model", "Group", "n", "StateEntropy_mean", "StateEntropy_SD",
                    "GroupLogDetEntropy", "N_indicators"])
        for em in p5.models:
            for g in ("Control", "Model", "HTD1801", "HTD1801-L", "HTD1801-M",
                      "HTD1801-H", "PM"):
                st = em.group_entropy_stats.get(g)
                if st:
                    w.writerow([em.model, g, st.n, _fmt(st.mean), _fmt(st.sd),
                                _fmt(em.group_logdet_entropy.get(g)), em.n_indicators])
        w.writerow([])
        w.writerow(["Model", "Comparison", "HTD_arm", "g", "SE", "CI_low", "CI_high", "P_value"])
        for em in p5.models:
            if em.model_vs_control:
                r = em.model_vs_control
                w.writerow([em.model, "Model vs Control", "", _fmt(r.g), _fmt(r.se),
                            _fmt(r.ci_low), _fmt(r.ci_high), _fmt(two_sided_p(r.g, r.se), 4)])
            for dose, r in em.htd_vs_pm.items():
                w.writerow([em.model, "HTD1801 vs PM", dose, _fmt(r.g), _fmt(r.se),
                            _fmt(r.ci_low), _fmt(r.ci_high), _fmt(two_sided_p(r.g, r.se), 4)])
        for dose, (fx, rd) in p5.pools.items():
            tag = "PRIMARY" if dose == p5.primary_dose else "sensitivity"
            w.writerow([f"POOLED ({tag})", "HTD1801 vs PM (random)", dose,
                        _fmt(rd.estimate), _fmt(rd.se), _fmt(rd.ci_low), _fmt(rd.ci_high),
                        _fmt(two_sided_p(rd.estimate, rd.se), 4)])


def build_markdown_part5(p5: "Part5Result") -> str:
    lines: List[str] = []
    a = lines.append
    a("# G-value recalculation — Part 5 (熵跨模型G值)")
    a("")
    a(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    a("")
    if not p5.available:
        a(f"Part 5 skipped: {p5.reason}.")
        return "\n".join(lines)
    a("## Method")
    a("- Per-animal state entropy E = 0.5·log(d²) + 0.5·log(2πe) (the per-animal "
      "form of bioentropy_entropy.compute_distance_state_entropy), where d is the "
      "Mahalanobis distance to the Control centre (Part-4 pipeline).")
    a("- Higher entropy ⇒ more disordered / farther from normal. HTD1801 vs PM: "
      "g = J·(mean_PM_E − mean_HTD_E)/SD_pooled, g>0 ⇒ HTD1801 lowers the entropy "
      "(restores order) more than PM.")
    a("- Group-level log-det Gaussian entropy (compute_group_entropy) reported per "
      "group as a descriptive check.")
    a(f"- Cross-model inverse-variance pool (fixed + random; random primary); "
      f"KK-ay dose = **{p5.primary_dose}** primary, L/M sensitivity.")
    a("")

    a("## Cross-model pooled entropy G (HTD1801 vs PM)")
    a("")
    a("| KK dose | k | Random g [95% CI] | Fixed g [95% CI] | I² |")
    a("|---|---|---|---|---|")
    for dose in PART3_KK_DOSES:
        p = p5.pools.get(dose)
        if not p:
            continue
        fx, rd = p
        tag = " (primary)" if dose == p5.primary_dose else ""
        a(f"| {dose}{tag} | {rd.k} | "
          f"{_fmt(rd.estimate)} [{_fmt(rd.ci_low)}, {_fmt(rd.ci_high)}] | "
          f"{_fmt(fx.estimate)} [{_fmt(fx.ci_low)}, {_fmt(fx.ci_high)}] | {_fmt(rd.i2,1)}% |")
    a("")
    if p5.primary_pool is not None:
        rd = p5.primary_pool[1]
        verdict = ("HTD1801 superior to PM" if rd.ci_low > 0 else
                   "PM superior" if rd.ci_high < 0 else
                   "no significant difference (CI spans 0)")
        a(f"**Headline (random, KK={p5.primary_dose}):** g = {_fmt(rd.estimate)} "
          f"[{_fmt(rd.ci_low)}, {_fmt(rd.ci_high)}], k={rd.k}, I²={_fmt(rd.i2,1)}% → {verdict}.")
        a("")

    for em in p5.models:
        a(f"## {em.model}  ({em.n_indicators} indicators)")
        a("")
        a("| Group | n | state-entropy mean | state-entropy SD | group log-det entropy |")
        a("|---|---|---|---|---|")
        for g in ("Control", "Model", "HTD1801", "HTD1801-L", "HTD1801-M",
                  "HTD1801-H", "PM"):
            st = em.group_entropy_stats.get(g)
            if st:
                a(f"| {g} | {st.n} | {_fmt(st.mean)} | {_fmt(st.sd)} | "
                  f"{_fmt(em.group_logdet_entropy.get(g))} |")
        a("")
        if em.model_vs_control:
            r = em.model_vs_control
            a(f"- Model vs Control: g = {_fmt(r.g)} [{_fmt(r.ci_low)}, {_fmt(r.ci_high)}] "
              f"({'modeling established' if r.ci_low > 0 else 'weak/questionable'})")
        if em.pm_vs_model:
            r = em.pm_vs_model
            a(f"- PM vs Model (entropy reduction): g = {_fmt(r.g)} [{_fmt(r.ci_low)}, {_fmt(r.ci_high)}]")
        for dose, r in em.htd_vs_model.items():
            a(f"- HTD1801 vs Model ({dose}, entropy reduction): g = {_fmt(r.g)} "
              f"[{_fmt(r.ci_low)}, {_fmt(r.ci_high)}]")
        a("")
        a("| HTD arm | HTD1801 vs PM g | SE | 95% CI |")
        a("|---|---|---|---|")
        for dose, r in em.htd_vs_pm.items():
            a(f"| {dose} | {_fmt(r.g)} | {_fmt(r.se)} | [{_fmt(r.ci_low)}, {_fmt(r.ci_high)}] |")
        a("")
    return "\n".join(lines)


def build_markdown(result: Dict[str, object]) -> str:
    lines: List[str] = []
    a = lines.append
    a("# Part 1 (individual model × endpoint G) & Part 2 (within-model pooled G)")
    a("")
    a(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    a("")
    a("All effects are HTD1801 vs PM, g>0 ⇒ HTD1801 better than PM. KK-ay uses the "
      f"medium dose {PART3_PRIMARY_KK_DOSE}. Exploratory composite analysis, not a "
      "confirmatory test.")
    a("")

    # ---- Part 1: individual endpoint G, grouped by model ----
    records: List[Part1Record] = result["part1"]  # type: ignore
    a("## Part 1 — individual per-endpoint G (no pooling)")
    a("")
    for model in dict.fromkeys(r.model for r in records):
        a(f"### {model}")
        a("")
        a("| Endpoint | dose | HTD mean (n) | PM mean (n) | g [95% CI] | dir (M vs C) | "
          "in Part-2 pool | reason if excluded |")
        a("|---|---|---|---|---|---|---|---|")
        for r in [r for r in records if r.model == model]:
            mvc = r.model_vs_control_g
            dirn = ("Model higher" if (mvc or 0) > 0 else "Model lower") if mvc is not None else "—"
            a(f"| {r.endpoint} | {r.dose} | {_fmt(r.htd_mean)} ({r.htd_n}) | "
              f"{_fmt(r.pm_mean)} ({r.pm_n}) | {_fmt(r.g)} [{_fmt(r.ci_low)}, {_fmt(r.ci_high)}] "
              f"| {dirn} | {'yes' if r.included else 'no'} | {exclusion_reason(r.status)} |")
        a("")

    # ---- Part 2: within-model pooled G (3 values) ----
    a("## Part 2 — within-model pooled G (folder 2 「pooled G值计算」, primary)")
    a("")
    a("Random-effects (DerSimonian-Laird) pool of the curated within-model endpoint "
      "set (folder 2); 95% CI is the standard random-effects interval, with a "
      "Hartung-Knapp (HK) CI alongside.")
    a("")
    a("| Model | k endpoints | Random g [95% CI] | HK 95% CI | I² | support |")
    a("|---|---|---|---|---|---|")
    within = result["part2_within"]  # type: ignore
    members = result["part2_within_members"]  # type: ignore
    for model, (_fx, rd) in within.items():
        hk = (f"[{_fmt(rd.hk_ci_low)}, {_fmt(rd.hk_ci_high)}]" if rd.hk_ci_low is not None else "—")
        a(f"| {model} | {rd.k} | {_fmt(rd.estimate)} [{_fmt(rd.ci_low)}, {_fmt(rd.ci_high)}] "
          f"| {hk} | {_fmt(rd.i2, 1)}% | {classify_support(rd)} |")
    a("")
    for model in within:
        a(f"- **{model}** pooled endpoints: {', '.join(members.get(model, []))}")
    a("")
    comp = result["part2_comprehensive"]  # type: ignore
    if comp:
        a("Sensitivity — comprehensive pool of all valid Part-1 endpoints (folder 1):")
        a("")
        a("| Model | k | Random g [95% CI] | HK 95% CI |")
        a("|---|---|---|---|")
        for model, (_fx, rd) in comp.items():
            hk = (f"[{_fmt(rd.hk_ci_low)}, {_fmt(rd.hk_ci_high)}]" if rd.hk_ci_low is not None else "—")
            a(f"| {model} | {rd.k} | {_fmt(rd.estimate)} [{_fmt(rd.ci_low)}, {_fmt(rd.ci_high)}] | {hk} |")
        a("")
    return "\n".join(lines)


# --------------------------------------------------------------------------- #
# Forest plots + software entry point
# --------------------------------------------------------------------------- #

_NPG = ["#3C5488", "#E64B35", "#00A087", "#4DBBD5", "#F39B7F", "#7E6148", "#8491B4"]
_TEXT = "#202A35"
_GRID = "#D8DEE8"
_CJK_FONTS = ["Microsoft YaHei", "SimHei", "Noto Sans CJK SC", "Source Han Sans SC",
              "PingFang SC", "WenQuanYi Zen Hei", "Arial Unicode MS"]


def _setup_fonts(plt) -> None:
    try:
        from matplotlib import font_manager
        available = {f.name for f in font_manager.fontManager.ttflist}
        for name in _CJK_FONTS:
            if name in available:
                plt.rcParams["font.sans-serif"] = [name, "DejaVu Sans", "Arial"]
                break
        else:
            plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Arial"]
        plt.rcParams["font.family"] = "sans-serif"
        plt.rcParams["axes.unicode_minus"] = False
    except Exception:
        pass


# Nature-style palette: restrained, high-contrast, no gradients.
_FC_POINT = "#3A3A3A"     # individual per-model G — dark gray filled circle
_FC_PRIMARY = "#1F3B66"   # primary pooled G — dark navy diamond
_FC_HEADLINE = "#8C1D1D"  # composite headline pooled G — dark red diamond
_FC_SENS = "#FFFFFF"      # sensitivity pooled — open (white) diamond
_FC_CI = "#4A4A4A"        # CI line — thin gray-black
_FC_ZERO = "#BBBBBB"      # zero reference — light gray dashed
_FC_SEP = "#E6E6E6"       # data/number separator
_XLABEL_DIR = "Favors PM  ←  0  →  Favors HTD1801"
_XLIM = (-1.0, 4.0)       # baseline data x-range; the upper bound auto-expands to fit CIs


def _auto_hi(rows, margin: float = 0.3, cap: float = 8.0) -> float:
    """Smallest integer upper x-bound that fits every finite CI in ``rows``.

    Keeps the forest CIs on-axis (no truncation arrowheads) while staying at a
    tidy integer tick. Returns 0.0 when there is nothing finite to fit, so the
    caller's baseline bound wins.
    """
    his = [r.get("hi") for r in rows
           if r.get("hi") is not None and math.isfinite(r.get("hi"))]
    if not his:
        return 0.0
    return float(min(cap, math.ceil(max(his) + margin)))


def _ensure_dir_path(path: str) -> None:
    if path:
        os.makedirs(path, exist_ok=True)


def _safe(text: str) -> str:
    return re.sub(r"[^0-9A-Za-z_一-鿿]+", "_", str(text)).strip("_") or "item"


def _diamond(ax, x, y, face, edge, half_h=0.34, half_w=0.06, z=4):
    from matplotlib.patches import Polygon
    ax.add_patch(Polygon([[x - half_w, y], [x, y + half_h], [x + half_w, y],
                          [x, y - half_h]], closed=True, facecolor=face,
                         edgecolor=edge, lw=1.1, zorder=z))


def _forest_on_ax(ax, rows, *, xlim=_XLIM, title=None, method_note=None,
                  show_numbers=True, direction=True, label_fs=8.5,
                  num_fs=8.0, title_fs=11.0, ann_frac=1.05):
    """Draw a consistent Nature-style forest plot on the given axis.

    rows: list of dicts {label, g, lo, hi, kind, k?, hk?, tier?, note?, excluded?}
      with kind in {"point", "pool_primary", "pool_headline", "pool_sens"}.
    The pooled G value, 95% CI (and HK CI / tier when present) are annotated as
    a right-hand column.
    """
    import math as _m
    n = len(rows)
    ys = list(range(n))[::-1]
    data_l, data_r = xlim
    data_r = max(data_r, _auto_hi(rows))   # expand upper bound so CIs stay on-axis
    span = data_r - data_l
    ann_w = span * ann_frac if show_numbers else 0.0
    full_r = data_r + ann_w
    half_w = span * 0.016
    ax.axvline(0, color=_FC_ZERO, lw=1.0, ls=(0, (4, 3)), zorder=0)
    for y, row in zip(ys, rows):
        kind = row.get("kind", "point")
        g, lo, hi = row["g"], row["lo"], row["hi"]
        clo, chi = max(lo, data_l), min(hi, data_r)
        if chi >= clo:
            ax.plot([clo, chi], [y, y], color=_FC_CI, lw=1.1, zorder=2,
                    solid_capstyle="butt")
        # arrowheads when the CI runs off the unified axis
        if lo < data_l:
            ax.plot([data_l], [y], marker="<", ms=4, color=_FC_CI, zorder=2)
        if hi > data_r:
            ax.plot([data_r], [y], marker=">", ms=4, color=_FC_CI, zorder=2)
        gx = min(max(g, data_l), data_r)
        excluded = row.get("excluded", False)
        if kind == "point":
            if excluded:  # per-model point not in the pool — open gray circle
                ax.scatter([gx], [y], marker="o", s=30, facecolor="white",
                           edgecolor="#9AA0A6", linewidth=1.0, zorder=3)
            else:
                ax.scatter([gx], [y], marker="o", s=30, facecolor=_FC_POINT,
                           edgecolor="none", zorder=3)
        elif kind == "pool_sens":
            _diamond(ax, gx, y, _FC_SENS, _FC_PRIMARY, half_w=half_w)
        elif kind == "pool_headline":
            _diamond(ax, gx, y, _FC_HEADLINE, "#000000", half_w=half_w)
        else:
            _diamond(ax, gx, y, _FC_PRIMARY, "#000000", half_w=half_w)
        if show_numbers:
            txt = f"{g:.2f} [{lo:.2f}, {hi:.2f}]"
            if row.get("hk") is not None:
                hkl, hkh = row["hk"]
                txt += f"; HK [{hkl:.2f}, {hkh:.2f}]"
            if row.get("p") is not None:
                _ps = fmt_p(row["p"])
                txt += ("; p" + _ps) if _ps.startswith("<") else ("; p=" + _ps)
            if row.get("k") is not None:
                txt += f" k={row['k']}"
            if row.get("tier"):
                txt += f" · {row['tier']}"
            if row.get("note"):
                txt += f"  ({row['note']})"
            ax.text(full_r - span * 0.03, y, txt, ha="right", va="center",
                    fontsize=num_fs, color=_TEXT)
    ax.set_yticks(ys)
    ax.set_yticklabels([r["label"] for r in rows], fontsize=label_fs, color=_TEXT)
    for y, row in zip(ys, rows):
        if row.get("kind", "point") != "point":
            ax.get_yticklabels()[ys.index(y)].set_fontweight("bold")
    ax.set_xlim(data_l, full_r)
    ax.set_ylim(-0.7, n - 0.3)
    ax.set_xticks(list(range(int(_m.ceil(data_l)), int(data_r) + 1)))
    if show_numbers:
        ax.axvline(data_r + span * 0.05, color=_FC_SEP, lw=0.9, zorder=0)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color("#9AA0A6")
    ax.spines["bottom"].set_bounds(data_l, data_r)
    ax.tick_params(left=False, labelsize=8, colors=_TEXT)
    mid = (data_l + data_r) / 2.0
    if direction:
        ax.text(mid, -0.16, _XLABEL_DIR, transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=8, color="#5A6472")
    # Title and method note stacked ABOVE the axes (no overlap).
    if title:
        ax.text(data_l, 1.135, title, transform=ax.get_xaxis_transform(),
                ha="left", va="bottom", fontsize=title_fs, color=_TEXT, fontweight="bold")
    if method_note:
        ax.text(data_l, 1.02, method_note, transform=ax.get_xaxis_transform(),
                ha="left", va="bottom", fontsize=7.2, color="#8A929C", style="italic")


def _pool_kind(rd, headline=False):
    """Diamond style by support tier (headline composites use the red diamond)."""
    return "pool_headline" if headline else "pool_primary"


def _supp_word(rd):
    return {"robust": "robust", "trend": "trend", "uncertain": "trend"}[classify_support(rd)]


# ---- individual standalone forest figure ---------------------------------- #

def _save_forest(rows, out_path, title, method_note, figs, height=None, ann_frac=1.05):
    if not rows:
        return
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return
    _setup_fonts(plt)
    h = height or max(2.4, 0.46 * len(rows) + 1.5)
    fig, ax = plt.subplots(figsize=(11.0, h))
    _forest_on_ax(ax, rows, title=title, method_note=method_note, ann_frac=ann_frac)
    fig.tight_layout()
    _ensure_dir_path(os.path.dirname(out_path))
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    figs.append(out_path)


def _kk_dose_for(model):
    return PART3_PRIMARY_KK_DOSE if model == "KK-ay" else "HTD1801"

_TIER = {"robust": "robust", "trend": "trend", "uncertain": "uncertain"}


def _pool_annot(rd):
    """Common annotation fields for a pooled random-effects estimate."""
    return {"k": rd.k, "tier": _TIER[classify_support(rd)],
            "hk": (rd.hk_ci_low, rd.hk_ci_high) if rd.hk_ci_low is not None else None}


def _part3_rows(result, full=True):
    """Cross-model rows. full=True → per-model points + pooled diamond per
    endpoint; full=False → only the 6 pooled diamonds (clean composite panel)."""
    rows = []
    for ind in result["part3"]:
        if ind.needs_manual_review or ind.main_pool is None:
            continue
        if full:
            for r in ind.results:
                if r.htd_vs_pm and (r.model != "KK-ay" or r.htd_label == _kk_dose_for(r.model)):
                    rows.append({"label": f"  {ind.name} · {r.model}", "g": r.htd_vs_pm.g,
                                 "lo": r.htd_vs_pm.ci_low, "hi": r.htd_vs_pm.ci_high,
                                 "kind": "point", "excluded": r.model not in ind.main_members,
                                 "p": two_sided_p(r.htd_vs_pm.g, r.htd_vs_pm.se)})
        rd = ind.main_pool[1]
        excl = [e.split(" (")[0] for e in ind.main_excluded]
        row = {"label": f"{ind.name} — pooled", "g": rd.estimate, "lo": rd.ci_low,
               "hi": rd.ci_high, "kind": "pool_primary", "p": two_sided_p(rd.estimate, rd.se)}
        row.update(_pool_annot(rd))
        if excl:
            row["note"] = "excl " + ", ".join(excl)
        rows.append(row)
    return rows


def _within_rows(result):
    """Within-model rows: each model's per-endpoint g points followed by its
    pooled diamond (mirrors the Part-3 per-model-points-plus-pool layout)."""
    rows = []
    points = result.get("part2_within_points", {})
    for model, (_fx, rd) in result["part2_within"].items():
        for pt in points.get(model, []):
            rows.append({"label": f"  {model} · {pt['metric']}", "g": pt["g"],
                         "lo": pt["lo"], "hi": pt["hi"], "kind": "point",
                         "p": two_sided_p(pt["g"], pt["se"])})
        row = {"label": f"{model} — within-model pooled", "g": rd.estimate,
               "lo": rd.ci_low, "hi": rd.ci_high, "kind": "pool_primary",
               "p": two_sided_p(rd.estimate, rd.se)}
        row.update(_pool_annot(rd))
        rows.append(row)
    return rows


def _composite_rows(result, families=("Mahalanobis", "Entropy")):
    """Parts 4 & 5: per-model points + winsorised primary + unrestricted sens."""
    rows = []
    src = {"Mahalanobis": result["part4"], "Entropy": result["part5"]}
    for label in families:
        res = src[label]
        if not getattr(res, "available", False):
            continue
        for mm in res.models:
            r = mm.htd_vs_pm.get(_kk_dose_for(mm.model))
            if r:
                rows.append({"label": f"  {label} · {mm.model}", "g": r.g, "lo": r.ci_low,
                             "hi": r.ci_high, "kind": "point", "p": two_sided_p(r.g, r.se)})
        if res.primary_pool:
            rd = res.primary_pool[1]
            row = {"label": f"{label} · winsorised (primary)", "g": rd.estimate,
                   "lo": rd.ci_low, "hi": rd.ci_high, "kind": "pool_headline",
                   "p": two_sided_p(rd.estimate, rd.se)}
            row.update(_pool_annot(rd))
            rows.append(row)
        if res.primary_pool_sens:
            rd = res.primary_pool_sens[1]
            row = {"label": f"{label} · unrestricted (sensitivity)", "g": rd.estimate,
                   "lo": rd.ci_low, "hi": rd.ci_high, "kind": "pool_sens",
                   "p": two_sided_p(rd.estimate, rd.se)}
            row.update(_pool_annot(rd))
            rows.append(row)
    return rows


def _part1_rows(records, model):
    """Part-1 individual endpoint points for one model (included=filled,
    excluded=open), annotated with g, CI, n and the exclusion reason."""
    rows = []
    for r in [r for r in records if r.model == model]:
        reason = exclusion_reason(r.status)
        note = f"n={r.htd_n}/{r.pm_n}" + (f"; {reason}" if reason else "")
        rows.append({"label": r.endpoint, "g": r.g, "lo": r.ci_low, "hi": r.ci_high,
                     "kind": "point", "excluded": not r.included, "note": note,
                     "p": two_sided_p(r.g, r.se)})
    return rows


def _panelA_pipeline(ax):
    from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 1)
    ax.axis("off")
    steps = ["Raw data\n(per animal)", "Part 1\nper-model\nendpoint G",
             "Part 2\nwithin-model\npooled G", "Part 3\ncross-model\nendpoint pooled G",
             "Parts 4–5\nMahalanobis /\nentropy pooled G"]
    n = len(steps)
    bw, gap = 1.62, (10 - 1.62 * n) / (n - 1)
    for i, s in enumerate(steps):
        x = i * (bw + gap)
        ax.add_patch(FancyBboxPatch((x, 0.28), bw, 0.46, boxstyle="round,pad=0.02,rounding_size=0.04",
                                    facecolor="#F4F6F8", edgecolor="#9AA0A6", lw=0.9))
        ax.text(x + bw / 2, 0.51, s, ha="center", va="center", fontsize=8.0, color=_TEXT)
        if i < n - 1:
            ax.add_patch(FancyArrowPatch((x + bw + 0.04, 0.51), (x + bw + gap - 0.04, 0.51),
                                         arrowstyle="-|>", mutation_scale=11, color="#7A828C", lw=1.0))


def _panelB_heatmap(ax, result):
    import numpy as np
    records = result["part1"]
    models = list(dict.fromkeys(r.model for r in records))
    endpoints = list(dict.fromkeys(r.endpoint for r in records))
    M = np.full((len(endpoints), len(models)), np.nan)
    lookup = {(r.model, r.endpoint): r.g for r in records}
    for i, ep in enumerate(endpoints):
        for j, m in enumerate(models):
            if (m, ep) in lookup:
                M[i, j] = lookup[(m, ep)]
    import matplotlib
    try:
        cmap = matplotlib.colormaps["RdBu_r"].copy()
    except (AttributeError, KeyError):
        import matplotlib.pyplot as _plt
        cmap = _plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#ECEEF0")
    M = np.ma.masked_invalid(M)
    im = ax.imshow(M, cmap=cmap, vmin=-3, vmax=3, aspect="auto")
    ax.set_xticks(range(len(models)))
    ax.set_xticklabels(models, fontsize=8.5)
    ax.set_yticks(range(len(endpoints)))
    ax.set_yticklabels(endpoints, fontsize=7.5)
    ax.tick_params(length=0)
    for i, ep in enumerate(endpoints):
        for j, m in enumerate(models):
            v = lookup.get((m, ep))
            if v is not None:
                ax.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=6.4,
                        color="white" if abs(v) > 1.6 else "#202A35")
            else:  # explicitly label missing cells (not measured in this model)
                ax.text(j, i, "ND", ha="center", va="center", fontsize=6.0,
                        color="#AEB4BC")
    for s in ax.spines.values():
        s.set_visible(False)
    ax.text(0.0, 1.012, "ND = no data / not measured in this model",
            transform=ax.transAxes, fontsize=6.8, color="#8A929C", va="bottom")
    cb = ax.figure.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
    cb.set_label("Hedges' g (HTD1801 vs PM)", fontsize=7.5)
    cb.ax.tick_params(labelsize=7)


def _make_panel_figure(result, out_path, figs):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec
        from matplotlib.lines import Line2D
        from matplotlib.patches import Patch
    except Exception:
        return
    _setup_fonts(plt)
    fig = plt.figure(figsize=(14.4, 14.8), facecolor="white")
    gs = GridSpec(4, 2, height_ratios=[0.55, 1.0, 1.0, 1.0], width_ratios=[0.92, 1.7],
                  hspace=0.82, wspace=0.5, left=0.135, right=0.985, top=0.93, bottom=0.08)
    axA = fig.add_subplot(gs[0, :])
    axB = fig.add_subplot(gs[1:, 0])
    axC = fig.add_subplot(gs[1, 1])
    axD = fig.add_subplot(gs[2, 1])
    axE = fig.add_subplot(gs[3, 1])

    def tag(ax, letter):
        ax.text(-0.02, 1.06, letter, transform=ax.transAxes, fontsize=15,
                fontweight="bold", va="bottom", ha="right", color=_TEXT)

    _panelA_pipeline(axA); axA.set_title("Analysis pipeline", fontsize=11, loc="left",
                                         color=_TEXT, pad=4); tag(axA, "a")
    _panelB_heatmap(axB, result); axB.set_title("Part 1 · individual G by model × endpoint",
                                                fontsize=11, loc="left", color=_TEXT, pad=12)
    tag(axB, "b")
    # Keep panels C/D/E on one shared x-range (auto-expanded to fit every CI).
    rows_c = _within_rows(result)
    rows_d = _part3_rows(result, full=False)
    rows_e = _composite_rows(result)
    shared_r = max(_XLIM[1], _auto_hi(rows_c), _auto_hi(rows_d), _auto_hi(rows_e))
    shared_xlim = (_XLIM[0], shared_r)
    _forest_on_ax(axC, rows_c, xlim=shared_xlim, title="Part 2 · within-model pooled G",
                  method_note="Random-effects (DerSimonian–Laird); RE & HK 95% CI",
                  direction=False, ann_frac=1.5, num_fs=7.2)
    tag(axC, "c")
    # Panel D: only the 6 cross-model pooled diamonds (clean); per-model points
    # are in the standalone Part-3 figure.
    _forest_on_ax(axD, rows_d, xlim=shared_xlim,
                  title="Part 3 · cross-model pooled G (6 shared endpoints)",
                  method_note="Random-effects across 3 models; RE & HK 95% CI · tier",
                  direction=False, ann_frac=1.5, num_fs=7.2)
    tag(axD, "d")
    _forest_on_ax(axE, rows_e, xlim=shared_xlim,
                  title="Parts 4–5 · cross-model pooled G (Mahalanobis & entropy)",
                  method_note="Per-model g → random-effects meta; winsorised = primary",
                  ann_frac=1.5, num_fs=7.2)
    tag(axE, "e")

    legend = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=_FC_POINT, markersize=7,
               label="Per-model G"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="white",
               markeredgecolor="#9AA0A6", markersize=7, label="Per-model G (excluded from pool)"),
        Patch(facecolor=_FC_PRIMARY, edgecolor="black", label="Pooled G (primary)"),
        Patch(facecolor=_FC_HEADLINE, edgecolor="black", label="Winsorised composite (primary)"),
        Patch(facecolor=_FC_SENS, edgecolor=_FC_PRIMARY, label="Pooled G (sensitivity)"),
    ]
    fig.legend(handles=legend, loc="lower center", ncol=5, frameon=False, fontsize=8.0,
               bbox_to_anchor=(0.5, 0.028))
    fig.text(0.5, 0.004, "Pooling: random-effects (DerSimonian–Laird). Annotations give the "
             "standard RE 95% CI and the Hartung–Knapp (HK) 95% CI with a tier: robust = HK "
             "excludes 0; trend = RE excludes 0 but HK crosses 0; uncertain = RE crosses 0. "
             "KK-ay uses the medium dose.", ha="center", va="bottom", fontsize=7.2,
             color="#76808C")
    fig.suptitle("HTD1801 vs PM — hierarchical G-value summary    "
                 "(positive G favors HTD1801: lower disease burden / closer to control)",
                 fontsize=12.5, color=_TEXT, x=0.5, y=0.975)
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    figs.append(out_path)


def _save_part1_figure(result, out_path, figs):
    """Part 1 standalone: per-model forest of every endpoint's individual G."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec
    except Exception:
        return
    _setup_fonts(plt)
    records = result["part1"]
    models = list(dict.fromkeys(r.model for r in records))
    per = {m: _part1_rows(records, m) for m in models}
    heights = [max(1, len(per[m])) for m in models]
    fig = plt.figure(figsize=(11.5, max(6.0, 0.34 * sum(heights) + 2.4)), facecolor="white")
    # Extra top headroom so the first panel's title clears the figure suptitle.
    gs = GridSpec(len(models), 1, height_ratios=heights, hspace=0.32,
                  left=0.2, right=0.97, top=0.90, bottom=0.055)
    # Shared x-range across the per-model panels (auto-expanded to fit every CI).
    shared_xlim = (_XLIM[0], max(_XLIM[1], *(_auto_hi(per[m]) for m in models)))
    for i, m in enumerate(models):
        ax = fig.add_subplot(gs[i])
        _forest_on_ax(ax, per[m], xlim=shared_xlim,
                      title=f"Part 1 · {m} — individual endpoint G (HTD1801 vs PM)",
                      method_note=("filled = included in Part-2 pool; open = excluded "
                                   "(reason at right)") if i == 0 else None,
                      direction=(i == len(models) - 1), ann_frac=1.35, num_fs=7.0, label_fs=8.0)
    fig.suptitle("Part 1 · individual model × endpoint Hedges' g", fontsize=12.5,
                 color=_TEXT, x=0.5, y=0.995)
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    figs.append(out_path)


def _contribution_plot(models, out_path, figs):
    """Supplementary QC: each endpoint's contribution share to the PART-4/5
    Mahalanobis/entropy composite (a different endpoint set than Part 2)."""
    rows = [(mm.model, mm.contributions) for mm in models if mm.contributions]
    if not rows:
        return
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return
    _setup_fonts(plt)
    fig, ax = plt.subplots(figsize=(11.5, max(2.4, 0.9 * len(rows) + 1.6)), facecolor="white")
    y = list(range(len(rows)))[::-1]
    greys = ["#3A3A3A", "#5C5C5C", "#7E7E7E", "#9E9E9E", "#BDBDBD", "#D6D6D6"]
    ep_order = []
    for _, c in rows:
        for nm, _s in c:
            if nm not in ep_order:
                ep_order.append(nm)
    cmap = {ep: greys[i % len(greys)] for i, ep in enumerate(ep_order)}
    for yi, (_m, contribs) in zip(y, rows):
        left = 0.0
        for nm, sh in contribs:
            ax.barh(yi, sh, left=left, color=cmap[nm], edgecolor="white", linewidth=0.5)
            if sh > 0.07:
                ax.text(left + sh / 2, yi, f"{nm}\n{sh*100:.0f}%", ha="center", va="center",
                        fontsize=6.6, color="white")
            left += sh
    ax.set_yticks(y)
    ax.set_yticklabels([f"{m} (k={len(c)})" for m, c in rows])
    ax.set_xlim(0, 1)
    ax.set_xlabel("Contribution share of the Mahalanobis/entropy composite (Model Σz²)", fontsize=9)
    ax.set_title("Supplementary · endpoint contributions to the Part-4/5 composite",
                 fontsize=11, color=_TEXT, pad=8, loc="left")
    ax.text(0.0, 1.10, "NB: this is the multi-indicator Mahalanobis/entropy composite "
            "(folder-4 wide table, all valid objective endpoints) — a different and larger "
            "endpoint set than the Part-2 within-model pool (folder 2, k shown in Part 2).",
            transform=ax.transAxes, fontsize=7.2, color="#8A929C", va="bottom")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    figs.append(out_path)


def _make_figures(result: Dict[str, object], fig_dir: str) -> List[str]:
    figs: List[str] = []
    # Nature-style multi-panel composite (primary deliverable).
    _make_panel_figure(result, os.path.join(fig_dir, "Figure_composite_main.png"), figs)
    # One standalone results figure per Part (1–5).
    _save_part1_figure(result, os.path.join(fig_dir, "Figure_Part1_individual.png"), figs)
    _save_forest(_within_rows(result), os.path.join(fig_dir, "Figure_Part2_within_model.png"),
                 "Part 2 · within-model pooled G (HTD1801 vs PM)",
                 "Random-effects pool of folder-2 endpoints; RE & HK 95% CI · tier",
                 figs, ann_frac=1.6)
    _save_forest(_part3_rows(result, full=True),
                 os.path.join(fig_dir, "Figure_Part3_cross_model_endpoints.png"),
                 "Part 3 · cross-model pooled G · 6 shared endpoints",
                 "Per-model g (gray; open = excluded from pool) → random-effects pool; RE & HK CI",
                 figs, ann_frac=1.6)
    _save_forest(_composite_rows(result, ("Mahalanobis",)),
                 os.path.join(fig_dir, "Figure_Part4_mahalanobis.png"),
                 "Part 4 · cross-model pooled G · Mahalanobis distance",
                 "Per-model g → RE meta; winsorised (±3 SD, valid endpoints) = primary, "
                 "unrestricted = sensitivity", figs, ann_frac=1.6)
    _save_forest(_composite_rows(result, ("Entropy",)),
                 os.path.join(fig_dir, "Figure_Part5_entropy.png"),
                 "Part 5 · cross-model pooled G · entropy",
                 "Per-model g → RE meta; winsorised (±3 SD, valid endpoints) = primary, "
                 "unrestricted = sensitivity", figs, ann_frac=1.6)
    p4 = result["part4"]
    if getattr(p4, "available", False):
        _contribution_plot(p4.models, os.path.join(fig_dir, "Figure_Supp_contributions.png"), figs)
    return figs


# --------------------------------------------------------------------------- #
# Quality-control tables + main/sensitivity interpretation
# --------------------------------------------------------------------------- #

def write_qc_tables(result: Dict[str, object], out_dir: str) -> None:
    qc = os.path.join(out_dir, "QC")
    _ensure_dir_path(qc)
    # 1) input validation
    with open(os.path.join(qc, "QC_input_validation.csv"), "w", newline="",
              encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(["Model", "Family", "File", "Status", "Groups_found", "Messages"])
        for v in result["validations"]:  # type: ignore
            w.writerow([v.model, v.family, v.file, v.status, v.groups_found,
                        " ".join(v.messages)])
    # 2) flagged / starred values
    with open(os.path.join(qc, "QC_flagged_values.csv"), "w", newline="",
              encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(["Model", "Endpoint", "Group", "Animal", "Original", "Kind", "Action"])
        for f in result["flagged"]:  # type: ignore
            w.writerow([f.model, f.endpoint, f.group, f.animal, f.original, f.kind, f.action])
    # 3) direction / validity checks
    with open(os.path.join(qc, "QC_direction_check.csv"), "w", newline="",
              encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(["Analysis", "Model", "Endpoint", "n_Control", "n_Model",
                    "Control_mean", "Model_mean", "Observed_change", "Category",
                    "Expected_disease_change", "Status", "Included_in_main"])
        for d in result["direction_checks"]:  # type: ignore
            w.writerow([d.analysis, d.model, d.endpoint, d.n_control, d.n_model,
                        _fmt(d.control_mean), _fmt(d.model_mean), d.observed,
                        d.category, d.expected or "", d.status, d.included_in_main])
    # 4) contribution proportions (robust composite)
    p4 = result["part4"]
    with open(os.path.join(qc, "QC_indicator_contributions.csv"), "w", newline="",
              encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(["Model", "Endpoint", "Contribution_share"])
        if getattr(p4, "available", False):
            for mm in p4.models:
                for nm, sh in mm.contributions:
                    w.writerow([mm.model, nm, _fmt(sh, 4)])


def _fdr_family(result: Dict[str, object]) -> List[Tuple[str, "Pooled"]]:
    """The family of headline cross-model comparisons for FDR adjustment."""
    fam: List[Tuple[str, Pooled]] = []
    for ind in result["part3"]:  # type: ignore
        if not ind.needs_manual_review and ind.main_pool is not None:
            fam.append((f"Part3:{ind.name}", ind.main_pool[1]))
    p4, p5 = result["part4"], result["part5"]
    if getattr(p4, "available", False) and p4.primary_pool:
        fam.append(("Part4:Mahalanobis", p4.primary_pool[1]))
    if getattr(p5, "available", False) and p5.primary_pool:
        fam.append(("Part5:Entropy", p5.primary_pool[1]))
    return fam


def write_fdr_table(result: Dict[str, object], out_dir: str) -> List[Tuple[str, Pooled, float, float]]:
    fam = _fdr_family(result)
    pvals = [two_sided_p(rd.estimate, rd.se) for _, rd in fam]
    qvals = benjamini_hochberg(pvals)
    rows = [(name, rd, p, q) for (name, rd), p, q in zip(fam, pvals, qvals)]
    with open(os.path.join(out_dir, "QC", "QC_FDR_main_comparisons.csv"), "w",
              newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(["Comparison", "Random_g", "SE", "p_raw", "q_FDR_BH",
                    "Significant_q<0.05"])
        for name, rd, p, q in rows:
            w.writerow([name, _fmt(rd.estimate), _fmt(rd.se), _fmt(p, 4), _fmt(q, 4),
                        (q < 0.05) if math.isfinite(q) else ""])
    return rows


def asterisk_sensitivity(repo_root: str) -> Dict[str, float]:
    """Re-run headline composites with starred values EXCLUDED, for comparison."""
    global ASTERISK_MODE
    saved = ASTERISK_MODE
    out: Dict[str, float] = {}
    try:
        ASTERISK_MODE = "exclude"
        r = run(repo_root)
        p4, p5 = r["part4"], r["part5"]
        if getattr(p4, "available", False) and p4.primary_pool:
            out["mahalanobis"] = p4.primary_pool[1].estimate
        if getattr(p5, "available", False) and p5.primary_pool:
            out["entropy"] = p5.primary_pool[1].estimate
        for m, pool in r["part2_within"].items():  # type: ignore
            out[f"part2_{m}"] = pool[1].estimate
    except Exception:
        pass
    finally:
        ASTERISK_MODE = saved
    return out


def build_main_summary(result: Dict[str, object], fdr_rows, ast_sens: Dict[str, float]) -> str:
    L: List[str] = []
    a = L.append
    a("# Hedges' g 值计算 — 主分析与稳健性总报告")
    a("")
    a(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    a("")
    a("> 本报告为**探索性复合效应分析（exploratory composite effect analysis）**，"
      "用于趋势与可视化，**不是确证性统计检验**；未对全部比较做严格多重校正（见 FDR 敏感性）。")
    a("")
    a("## 方法约定（主分析）")
    a(f"- KK-Ay 模型一律使用 **中剂量 {PART3_PRIMARY_KK_DOSE}** 代表 HTD1801（高/低剂量仅作剂量-反应敏感性）。")
    a(f"- 星号数据默认**保留**（ASTERISK_MODE={ASTERISK_MODE}）；是否剔除需经确认，见 QC 表与排除/保留敏感性。")
    a("- 复合评分（马氏距离/熵）**主结果=稳健版**（winsor 截尾 + 仅纳入“有效客观指标”），"
      "不受限版作为敏感性。")
    a("- 仅当 Model vs Control 的变化方向与疾病逻辑一致时，该指标才纳入主分析；"
      "胰岛素/体重/摄食/饮水等需人工确认（needs_manual_review），暂不进入主结果。")
    a("- 跨模型合并同时给出标准随机效应与 **Hartung-Knapp（HK）** 95% CI；HK CI 跨 0 时按“趋势”解读。")
    a("")

    # ---- MAIN results table ----
    a("## 一、主分析结果（HTD1801 中剂量 vs PM）")
    a("")
    a("| 分析 | g (随机效应) | 95% CI | HK 95% CI | 证据等级 |")
    a("|---|---|---|---|---|")

    def row(name, pool):
        if pool is None:
            a(f"| {name} | — | — | — | 不可计算 |")
            return
        rd = pool[1]
        tier = {"robust": "稳健支持", "trend": "趋势", "uncertain": "不确定"}[classify_support(rd)]
        hk = (f"[{_fmt(rd.hk_ci_low)}, {_fmt(rd.hk_ci_high)}]" if rd.hk_ci_low is not None else "—")
        a(f"| {name} | {_fmt(rd.estimate)} | [{_fmt(rd.ci_low)}, {_fmt(rd.ci_high)}] | {hk} | {tier} |")

    for m, pool in result["part2_within"].items():  # type: ignore
        members = result["part2_within_members"].get(m, [])  # type: ignore
        row(f"Part2 模型内合并 · {m}（{len(members)} 指标）", pool)
    for ind in result["part3"]:  # type: ignore
        if ind.needs_manual_review:
            continue
        excl = f"（排除 {', '.join(ind.main_excluded)}）" if ind.main_excluded else ""
        row(f"Part3 跨模型 · {ind.name}{excl}", ind.main_pool)
    p4, p5 = result["part4"], result["part5"]
    row("Part4 马氏距离复合（稳健）", p4.primary_pool if getattr(p4, "available", False) else None)
    row("Part5 熵复合（稳健）", p5.primary_pool if getattr(p5, "available", False) else None)
    a("")

    # ---- needs_manual_review ----
    nmr = [ind for ind in result["part3"] if ind.needs_manual_review]  # type: ignore
    if nmr:
        a("## 二、需人工确认的指标（needs_manual_review，未计入主结果）")
        a("")
        a("以下指标的疾病方向取决于具体疾病背景，请确认后再决定是否纳入主分析：")
        for ind in nmr:
            pool = ind.primary_pool
            note = (f"（参考：剂量 {PART3_PRIMARY_KK_DOSE} 跨模型 g="
                    f"{_fmt(pool[1].estimate)}）" if pool else "")
            a(f"- **{ind.name}** {note}")
        a("")

    # ---- direction-check summary ----
    a("## 三、方向/有效性检查（Model vs Control）")
    a("")
    a("| 模型 | 指标 | Control 均值 | Model 均值 | 观察变化 | 类别 | 判定 | 计入主分析 |")
    a("|---|---|---|---|---|---|---|---|")
    seen = set()
    for d in result["direction_checks"]:  # type: ignore
        key = (d.model, d.endpoint, d.analysis)
        if key in seen:
            continue
        seen.add(key)
        a(f"| {d.model} | {d.endpoint} | {_fmt(d.control_mean)} | {_fmt(d.model_mean)} "
          f"| {d.observed} | {d.category} | {d.status} | {'是' if d.included_in_main else '否'} |")
    a("")

    # ---- sensitivity ----
    a("## 四、敏感性分析")
    a("")
    a("**(a) 复合评分：稳健（主）vs 不受限。**")
    a("")
    a("| 复合 | 稳健 g (HK CI) | 不受限 g |")
    a("|---|---|---|")
    if getattr(p4, "available", False) and p4.primary_pool:
        rd, rs = p4.primary_pool[1], (p4.primary_pool_sens[1] if p4.primary_pool_sens else None)
        a(f"| 马氏距离 | {_fmt(rd.estimate)} [{_fmt(rd.hk_ci_low)}, {_fmt(rd.hk_ci_high)}] "
          f"| {_fmt(rs.estimate) if rs else '—'} |")
    if getattr(p5, "available", False) and p5.primary_pool:
        rd, rs = p5.primary_pool[1], (p5.primary_pool_sens[1] if p5.primary_pool_sens else None)
        a(f"| 熵 | {_fmt(rd.estimate)} [{_fmt(rd.hk_ci_low)}, {_fmt(rd.hk_ci_high)}] "
          f"| {_fmt(rs.estimate) if rs else '—'} |")
    a("")
    a("**(b) KK-Ay 剂量-反应（跨模型合并 g，随机效应）。**")
    a("")
    a("| 指标 | KK 低剂量 | KK 中剂量（主） | KK 高剂量 |")
    a("|---|---|---|---|")
    for ind in result["part3"]:  # type: ignore
        cells = []
        for dose in PART3_KK_DOSES:
            p = ind.pools.get(dose)
            cells.append(_fmt(p[1].estimate) if p else "—")
        a(f"| {ind.name} | {cells[0]} | {cells[1]} | {cells[2]} |")
    a("")
    if ast_sens:
        a("**(c) 星号数据：保留（主）vs 剔除。**")
        a("")
        a("| 复合/合并 | 保留（主） | 剔除 |")
        a("|---|---|---|")
        p4p = p4.primary_pool[1].estimate if getattr(p4, "available", False) and p4.primary_pool else None
        p5p = p5.primary_pool[1].estimate if getattr(p5, "available", False) and p5.primary_pool else None
        a(f"| 马氏距离 | {_fmt(p4p)} | {_fmt(ast_sens.get('mahalanobis'))} |")
        a(f"| 熵 | {_fmt(p5p)} | {_fmt(ast_sens.get('entropy'))} |")
        for m, pool in result["part2_within"].items():  # type: ignore
            a(f"| Part2 {m} | {_fmt(pool[1].estimate)} | {_fmt(ast_sens.get(f'part2_{m}'))} |")
        a("")

    # ---- FDR ----
    a("**(d) 多重比较 FDR（Benjamini-Hochberg）敏感性（主跨模型比较族）。**")
    a("")
    a("| 比较 | g | p (原始) | q (FDR) | q<0.05 |")
    a("|---|---|---|---|---|")
    for name, rd, p, q in fdr_rows:
        a(f"| {name} | {_fmt(rd.estimate)} | {_fmt(p,4)} | {_fmt(q,4)} | "
          f"{'是' if (math.isfinite(q) and q < 0.05) else '否'} |")
    a("")

    # ---- interpretation framework ----
    a("## 五、解释框架（证据分级）")
    a("")
    fam = _fdr_family(result)
    robust = [n for n, rd in fam if classify_support(rd) == "robust"]
    trend = [n for n, rd in fam if classify_support(rd) == "trend"]
    uncertain = [n for n, rd in fam if classify_support(rd) == "uncertain"]
    nmr_names = [ind.name for ind in result["part3"] if ind.needs_manual_review]  # type: ignore
    a("- **稳健支持（robustly supported）**：随机效应与 HK CI 均不跨 0 — "
      + ("、".join(robust) if robust else "（无）"))
    a("- **趋势（trend-supported）**：随机效应显著但 HK CI 跨 0（仅 3 个模型，谨慎）— "
      + ("、".join(trend) if trend else "（无）"))
    a("- **不确定 / 未纳入主分析（uncertain）**：CI 跨 0，或 needs_manual_review / "
      "model_not_valid_for_this_endpoint — "
      + ("、".join(uncertain + nmr_names) if (uncertain or nmr_names) else "（无）"))
    a("")
    a("> 结论以“稳健支持”为主，“趋势”仅作提示，“不确定”不支撑主要结论。")
    a("")
    return "\n".join(L)


def _write_all_outputs(result: Dict[str, object], out_dir: str) -> Dict[str, str]:
    """Write every CSV + markdown report into out_dir. Returns key paths."""
    write_part1_csv(os.path.join(out_dir, "gvalue_part1_individual_endpoint_G.csv"), result["part1"])  # type: ignore
    write_part2_csv(os.path.join(out_dir, "gvalue_part2_within_model_pooled.csv"), result)
    write_part3_models_csv(os.path.join(out_dir, "gvalue_part3_per_model.csv"), result["part3"])  # type: ignore
    write_part3_pooled_csv(os.path.join(out_dir, "gvalue_part3_pooled.csv"), result["part3"])  # type: ignore
    write_part4_csv(os.path.join(out_dir, "gvalue_part4_distances.csv"),
                    os.path.join(out_dir, "gvalue_part4_g.csv"), result["part4"])  # type: ignore
    write_part5_csv(os.path.join(out_dir, "gvalue_part5_entropy.csv"), result["part5"])  # type: ignore
    reports = {
        "gvalue_part1_2_report.md": build_markdown(result),
        "gvalue_part3_report.md": build_markdown_part3(result["part3"], result["warnings"], result["dropped"]),  # type: ignore
        "gvalue_part4_report.md": build_markdown_part4(result["part4"], result["warnings"]),  # type: ignore
        "gvalue_part5_report.md": build_markdown_part5(result["part5"]),  # type: ignore
    }
    for name, text in reports.items():
        with open(os.path.join(out_dir, name), "w", encoding="utf-8") as fh:
            fh.write(text)
    return reports


def analyze_gvalues(input_dir: str, output_parent: str) -> Dict[str, object]:
    """Software entry point: compute all 5 G-value families from the
    G值森林图数据 directory (or a repo root containing data/), then export
    CSV tables, markdown reports and forest-plot PNGs into a timestamped
    bundle (+ zip). Mirrors analyze_hedges_forest's return contract.
    """
    gdir = find_gvalue_dir(str(input_dir))
    result = run(gdir)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = os.path.join(str(output_parent), f"BioEntropy_GValues_{stamp}")
    fig_dir = os.path.join(out_dir, "figures")
    _ensure_dir_path(fig_dir)
    _write_all_outputs(result, out_dir)
    figures = _make_figures(result, fig_dir)

    # QC tables + FDR + asterisk-exclude sensitivity + comprehensive summary
    write_qc_tables(result, out_dir)
    fdr_rows = write_fdr_table(result, out_dir)
    ast_sens = asterisk_sensitivity(gdir)
    summary_md = build_main_summary(result, fdr_rows, ast_sens)
    with open(os.path.join(out_dir, "00_MAIN_SUMMARY.md"), "w", encoding="utf-8") as fh:
        fh.write(summary_md)

    p4 = result["part4"]
    p5 = result["part5"]
    headline = {}
    if getattr(p4, "available", False) and p4.primary_pool:
        headline["mahalanobis_pooled_random_g"] = p4.primary_pool[1].estimate
        headline["mahalanobis_support"] = classify_support(p4.primary_pool[1])
    if getattr(p5, "available", False) and p5.primary_pool:
        headline["entropy_pooled_random_g"] = p5.primary_pool[1].estimate
        headline["entropy_support"] = classify_support(p5.primary_pool[1])

    # Which analyses actually ran (based on the data present).
    parts_done: List[str] = []
    if result["part1"]:
        parts_done.append("1）单个指标 G 值")
    if result["part2_within"]:
        parts_done.append("2）模型内合并 G 值")
    if result["part3"]:
        parts_done.append("3）不同指标跨模型 G 值")
    if getattr(p4, "available", False) and getattr(p4, "models", []):
        parts_done.append("4）马氏距离跨模型 G 值")
    if getattr(p5, "available", False) and getattr(p5, "models", []):
        parts_done.append("5）熵跨模型 G 值")

    # No zip archive is produced: results are delivered as the output folder
    # only (kept as an empty string for backward-compatible callers).
    zip_path = ""

    return {
        "output_dir": out_dir,
        "zip_path": zip_path,
        "figures": figures,
        "figure_count": len(figures),
        "warning_count": len(result["warnings"]),  # type: ignore
        "dropped_count": len(result["dropped"]),  # type: ignore
        "headline": headline,
        "parts_done": parts_done,
        "part_counts": {
            "part1_rows": len(result["part1"]),  # type: ignore
            "part3_indicators": len(result["part3"]),  # type: ignore
            "part4_models": len(getattr(p4, "models", [])),
            "part5_models": len(getattr(p5, "models", [])),
        },
    }


def main() -> None:
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    result = run(repo_root)
    out_dir = os.path.join(repo_root, "BioEntropy_Output")
    os.makedirs(out_dir, exist_ok=True)
    write_part1_csv(os.path.join(out_dir, "gvalue_part1_individual_endpoint_G.csv"), result["part1"])  # type: ignore
    write_part2_csv(os.path.join(out_dir, "gvalue_part2_within_model_pooled.csv"), result)
    md = build_markdown(result)
    with open(os.path.join(out_dir, "gvalue_recompute_report.md"), "w", encoding="utf-8") as fh:
        fh.write(md)

    part3 = result["part3"]  # type: ignore
    write_part3_models_csv(os.path.join(out_dir, "gvalue_part3_per_model.csv"), part3)
    write_part3_pooled_csv(os.path.join(out_dir, "gvalue_part3_pooled.csv"), part3)
    md3 = build_markdown_part3(part3, result["warnings"], result["dropped"])  # type: ignore
    with open(os.path.join(out_dir, "gvalue_part3_report.md"), "w", encoding="utf-8") as fh:
        fh.write(md3)

    part4 = result["part4"]  # type: ignore
    write_part4_csv(os.path.join(out_dir, "gvalue_part4_distances.csv"),
                    os.path.join(out_dir, "gvalue_part4_g.csv"), part4)
    md4 = build_markdown_part4(part4, result["warnings"])  # type: ignore
    with open(os.path.join(out_dir, "gvalue_part4_report.md"), "w", encoding="utf-8") as fh:
        fh.write(md4)

    part5 = result["part5"]  # type: ignore
    write_part5_csv(os.path.join(out_dir, "gvalue_part5_entropy.csv"), part5)
    md5 = build_markdown_part5(part5)
    with open(os.path.join(out_dir, "gvalue_part5_report.md"), "w", encoding="utf-8") as fh:
        fh.write(md5)

    for chunk in (md, md3, md4, md5):
        print("\n\n" + "=" * 78 + "\n\n")
        print(chunk)


if __name__ == "__main__":
    main()
