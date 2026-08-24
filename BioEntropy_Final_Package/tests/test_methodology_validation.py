"""Synthetic-data validation of the core claim machinery.

These tests construct data with a *known* truth and confirm the pipeline:
  - recovers a positive entropy-deviation association when group dispersion
    genuinely grows with disease deviation (positive control);
  - does NOT manufacture a strong association when dispersion is constant while
    deviation varies (null control);
  - loses the association when group labels are permuted (negative control).

They deliberately use plain group names (no "normal"/"model" keywords) so the
reported Entropy is the covariance entropy -- the measure that is *independent*
of deviation. (When a Normal reference exists the pipeline instead reports the
normal-reference state entropy, which tracks deviation almost tautologically;
see METHODS.md and the explicit test below.)
"""
import numpy as np
import pandas as pd

import bioentropy_core_generalized as G


def _synth_raw(coupled, seed, p=6, n=16, shifts=(0.0, 0.8, 1.6, 2.4, 3.2),
               group_names=None, permute=False):
    rng = np.random.default_rng(seed)
    names = group_names or [f"G{i}" for i in range(len(shifts))]
    rows = []
    sid = 0
    for gi, shift in enumerate(shifts):
        disp = (1.0 + 0.7 * gi) if coupled else 1.8
        for _ in range(n):
            sid += 1
            x = rng.normal(shift, disp, size=p)
            rows.append({"Group": names[gi], "Sample ID": sid,
                         **{f"f{j}": x[j] for j in range(p)}})
    df = pd.DataFrame(rows)
    if permute:
        df["Group"] = rng.permutation(df["Group"].to_numpy())
    return {("SYN", 0): df}


def _rho(raw):
    res = G.analyze_all(raw, role_map=None, winsor_q=G.WINSOR_Q_DEFAULT,
                        external_range_map={}, true_normal_map={})
    return float(res["tables"]["实验内相关性"].iloc[0].get(
        "Spearman_rho_BRI_vs_Entropy", float("nan")))


def _median_rho(coupled, seeds=(1, 2, 3, 4, 5), **kw):
    return float(np.median([_rho(_synth_raw(coupled, s, **kw)) for s in seeds]))


def test_positive_control_recovers_association():
    # Dispersion grows with deviation -> covariance entropy must track deviation.
    assert _median_rho(coupled=True) >= 0.8


def test_null_control_no_strong_association():
    # Constant dispersion -> entropy should not strongly track deviation.
    coupled = _median_rho(coupled=True)
    null = _median_rho(coupled=False)
    assert null < 0.8
    assert coupled - null >= 0.3   # the method clearly separates signal from null


def test_label_permutation_negative_control():
    # Permuting group labels destroys the dispersion-deviation coupling.
    truth = np.median([_rho(_synth_raw(True, s)) for s in (1, 2, 3)])
    permuted = np.median([_rho(_synth_raw(True, s, permute=True)) for s in (1, 2, 3)])
    assert truth - permuted >= 0.3


def test_reference_state_entropy_coupling_is_documented():
    # KNOWN behaviour (see METHODS.md): with a Normal reference present, the
    # reported Entropy is the reference-state entropy, which tracks deviation
    # even when covariance dispersion is constant. Pin it so it can't change
    # silently and be misread as an independent finding.
    raw = _synth_raw(coupled=False, seed=1,
                     group_names=["normal control", "model control", "T1", "T2", "T3"])
    res = G.analyze_all(raw, role_map=None, winsor_q=G.WINSOR_Q_DEFAULT,
                        external_range_map={}, true_normal_map={})
    summary = res["tables"]["组时间点汇总"]
    assert (summary["Entropy_method"] == "normal_reference_state_entropy").all()
    rho = res["tables"]["实验内相关性"].iloc[0]["Spearman_rho_BRI_vs_Entropy"]
    assert rho >= 0.8  # tautological coupling under the null -> interpret with care
