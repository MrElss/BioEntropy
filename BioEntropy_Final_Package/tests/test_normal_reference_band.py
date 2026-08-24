"""Contract tests for the Normal-reference NMD band.

The band is the distance-from-Normal that the healthy reference subjects
themselves carry (same robust-capped Mahalanobis transform). It is the range a
treatment arm must reach to be indistinguishable from Normal, and is drawn as a
shaded region on the NMD trajectory. These tests pin that:

  * the band is computed and exposed on the experiment result when a true Normal
    reference exists, with sane ordering of its summary statistics;
  * the healthy band sits well below the (still-diseased) patient NMD values, so
    the plot honestly shows the residual gap rather than implying normalization.
"""
import glob
from pathlib import Path

import pandas as pd
import pytest

import bioentropy_core_legacy as L

DATA = str(Path(__file__).resolve().parent.parent.parent / "data"
           / "临床实验105和106" / "临床实验105")


@pytest.fixture(scope="module")
def result_105():
    files = glob.glob(f"{DATA}/*.xlsx")
    if not files:
        pytest.skip("clinical 105 data not available")
    raw, _rr, nrf, *_ = L.load_input_bundle(files)
    sub = {k: v for k, v in raw.items() if k[0] == "105"}
    return L.analyze_experiment(sub, "105", treatment_group="HTD1801",
                                comparator_group="安慰剂", normal_sample_df=nrf.get("105"))


def test_band_present_and_ordered(result_105):
    band = result_105["normal_reference_band"]
    assert band is not None and band["metric"] == "NMD"
    assert band["n"] >= 3
    # percentiles and CI must be monotone / well-formed
    assert band["p2_5"] <= band["p25"] <= band["median"] <= band["p75"] <= band["p97_5"]
    assert band["mean_ci_low"] <= band["mean"] <= band["mean_ci_high"]


def test_band_sits_below_patient_nmd(result_105):
    band = result_105["normal_reference_band"]
    s = result_105["summary"]
    patient_min = pd.to_numeric(
        s[s["Group"].isin(["HTD1801", "安慰剂"])]["NMD_mean"], errors="coerce"
    ).min()
    # even the lowest treated group-mean NMD stays above the healthy 95% range top
    assert band["p97_5"] < patient_min


def test_band_none_without_true_normal():
    files = glob.glob(f"{DATA}/*.xlsx")
    if not files:
        pytest.skip("clinical 105 data not available")
    raw, *_ = L.load_input_bundle(files)
    sub = {k: v for k, v in raw.items() if k[0] == "105"}
    # no normal_sample_df → NRBS/BRI path, no band
    res = L.analyze_experiment(sub, "105", treatment_group="HTD1801", comparator_group="安慰剂")
    assert res["normal_reference_band"] is None
