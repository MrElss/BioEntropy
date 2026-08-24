"""Contract tests pinning the endpoint-inclusion policy to what PROJECT_NOTES.md and the
module docstrings promise.

These guard the two places where the code and the docs previously drifted:

  * Insulin and body weight are user-confirmed *objective higher_worse*
    endpoints (``CONFIRMED_HIGHER_WORSE_TOKENS``) and ENTER the main analysis
    when the Model is actually worse than Control.
  * Food/water intake and body-weight change-rate remain ``context`` endpoints
    (``needs_manual_review``) and stay OUT of the main results.

If the inclusion policy is intentionally changed again, update PROJECT_NOTES.md and the
``classify_endpoint`` / ``assess_endpoint`` docstrings together with this test so
the doc/code contract stays honest.
"""
import pytest

from gvalue_recompute import assess_endpoint, classify_endpoint


@pytest.mark.parametrize(
    "name",
    ["Serum insulin", "血清胰岛素", "insulin", "Body weight", "体重", "bodyweight"],
)
def test_insulin_and_body_weight_are_objective_higher_worse(name):
    category, expected = classify_endpoint(name)
    assert category == "objective"
    assert expected == "up"


@pytest.mark.parametrize(
    "name",
    ["摄食量", "food intake", "饮水", "water intake", "体重变化率", "weight change rate"],
)
def test_food_water_and_change_rate_stay_context(name):
    category, expected = classify_endpoint(name)
    assert category == "context"
    assert expected is None


def test_gfr_is_objective_lower_worse():
    assert classify_endpoint("GFR") == ("objective", "down")


def test_insulin_included_when_model_worse_than_control():
    # Model > Control matches the expected "up" disease direction -> included.
    dc = assess_endpoint("DB", "Serum insulin",
                         control_vals=[1.0, 1.1, 0.9, 1.05],
                         model_vals=[3.0, 3.2, 2.8, 3.1])
    assert dc.status == "valid"
    assert dc.included_in_main is True


def test_body_weight_included_when_model_worse_than_control():
    dc = assess_endpoint("DB", "Body weight",
                         control_vals=[20.0, 21.0, 19.5, 20.5],
                         model_vals=[35.0, 36.0, 34.0, 35.5])
    assert dc.status == "valid"
    assert dc.included_in_main is True


def test_insulin_excluded_when_model_not_worse():
    # Model < Control contradicts the fixed disease direction -> excluded.
    dc = assess_endpoint("DB", "Serum insulin",
                         control_vals=[3.0, 3.2, 2.8, 3.1],
                         model_vals=[1.0, 1.1, 0.9, 1.05])
    assert dc.status == "model_not_valid_for_this_endpoint"
    assert dc.included_in_main is False


def test_food_intake_always_needs_manual_review():
    # Even with Model > Control, a context endpoint stays out of the main pool.
    dc = assess_endpoint("DB", "food intake",
                         control_vals=[10.0, 11.0, 9.5, 10.5],
                         model_vals=[20.0, 21.0, 19.0, 20.5])
    assert dc.status == "needs_manual_review"
    assert dc.included_in_main is False
