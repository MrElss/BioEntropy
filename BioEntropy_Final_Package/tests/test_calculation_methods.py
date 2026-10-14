"""Independent scientific checks for the intended V1.0 method revision."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.spatial.distance import mahalanobis
from scipy.stats import chi2, mannwhitneyu

import bioentropy_core_legacy as L
import bioentropy_entropy as E
from bioentropy_constants import SOFTWARE_VERSION, NMD_PRIMARY_METHOD
import _public_regression as R


def _reference(covariance, m=20):
    covariance = np.asarray(covariance, dtype=float)
    p = len(covariance)
    features = [f'f{i}' for i in range(p)]
    return L.TrueNormalReference(
        experiment='SYN', features=features, mean=np.zeros(p),
        precision=np.linalg.inv(covariance), covariance=covariance,
        source_feature_map={f: f for f in features},
        unit_factors={f: 1.0 for f in features}, use_log_map={f: False for f in features},
        per_feature_mean={f: 0.0 for f in features},
        per_feature_sd=dict(zip(features, np.sqrt(np.diag(covariance)))),
        n_normal=m, source_files=[], clip_bounds={f: (-4.0, 4.0) for f in features},
    )


def test_full_distance_distinguishes_correlated_directions():
    # Both points have the same marginal standardized magnitude. A full
    # covariance distance must penalize the rare, anti-correlated direction.
    rho = 0.8
    ref = _reference([[1, rho], [rho, 1]])
    result = L.apply_true_normal_reference(pd.DataFrame({'f0': [1., 1.], 'f1': [1., -1.]}), ref)
    assert result['NMD_sq'].tolist() == pytest.approx([2/(1+rho), 2/(1-rho)])
    assert result['NMD_diagonal_capped'].tolist() == pytest.approx([np.sqrt(2)]*2)
    assert (result['NMD_method'] == NMD_PRIMARY_METHOD).all()


def test_primary_extreme_is_uncapped_but_sensitivity_keeps_cap():
    result = L.apply_true_normal_reference(pd.DataFrame({'f0': [20.], 'f1': [20.]}), _reference(np.eye(2)))
    assert result['NMD'].iloc[0] == pytest.approx(np.sqrt(800))
    assert result['NMD_diagonal_capped'].iloc[0] == pytest.approx(np.sqrt(128))
    assert np.isnan(result['NMD_z_cap'].iloc[0])
    assert result['NMD_diagonal_z_cap'].iloc[0] == 8
    assert result['NMD_p_chi2'].iloc[0] == pytest.approx(chi2.sf(800, 2))


def test_full_distance_agrees_with_scipy_and_consistent_sign_change():
    cov = np.array([[2., .7], [.7, 1.]])
    ref = _reference(cov)
    points = np.array([[1., 3.], [-2., .5]])
    result = L.apply_true_normal_reference(pd.DataFrame(points, columns=ref.features), ref)
    expected = [mahalanobis(point, ref.mean, ref.precision) for point in points]
    assert result['NMD'].tolist() == pytest.approx(expected)
    sign = np.diag([1., -1.])
    flipped = L.apply_true_normal_reference(pd.DataFrame(points @ sign, columns=ref.features), _reference(sign @ cov @ sign))
    assert flipped['NMD'].tolist() == pytest.approx(expected)


def test_missing_features_and_small_reference_are_explicit():
    result = L.apply_true_normal_reference(pd.DataFrame({'f0': [np.nan, 1.], 'f1': [np.nan, 2.]}), _reference(np.eye(2), m=3))
    assert result['NMD'].iloc[0] == 0
    assert result['NMD_feature_count'].tolist() == [0., 2.]
    assert result['NMD_p_smallsample_F'].isna().all()


def test_ledoitwolf_reference_supports_more_features_than_samples():
    rng = np.random.default_rng(181)
    features = [f'f{i}' for i in range(7)]
    study = pd.DataFrame(rng.normal(size=(12, 7)), columns=features)
    study['Group'] = ['A']*6+['B']*6
    study['Sample ID'] = np.arange(12)
    normal = pd.DataFrame(rng.normal(size=(4, 7)), columns=features)
    ref, _ = L.fit_true_normal_reference({('SYN', 0): study}, 'SYN', normal)
    actual = L.apply_true_normal_reference(study, ref)
    assert len(ref.features) > ref.n_normal
    expected = [mahalanobis(row, ref.mean, ref.precision) for row in study[ref.features].to_numpy()]
    assert actual['NMD'].tolist() == pytest.approx(expected)
    assert np.isfinite(actual['NMD']).all()
    assert actual['NMD_p_smallsample_F'].isna().all()


def test_entropy_keeps_extremes_and_does_not_use_outlier_votes(monkeypatch):
    points = np.array([[0., 0.], [1., 2.], [2., 1.], [100., -100.]])
    before = points.copy()
    def forbid_votes(*args):
        raise AssertionError('Outlier flags must not affect covariance entropy')
    monkeypatch.setattr(E, 'detect_outliers_methods', forbid_votes)
    value, used, info = E.compute_group_entropy_robust(points)
    np.testing.assert_array_equal(points, before)
    np.testing.assert_array_equal(used, points)
    assert value == pytest.approx(E.compute_group_entropy(points))
    assert info['n_entropy_outliers'] == 0
    assert 'no subject-level attenuation' in info['entropy_outlier_rule']


def test_qc_strong_flags_are_not_adjusted_subjects(monkeypatch, tmp_path):
    def flag_everyone(X):
        return {k: np.ones(len(X), dtype=bool) for k in ['robust_z', 'iqr', 'iso', 'lof', 'rmd']}
    monkeypatch.setattr(L, 'detect_outliers_methods', flag_everyone)
    results, *_ = R.run_case(tmp_path, 'bri_single')
    summary = results['tables']['组时间点汇总']
    np.testing.assert_allclose(summary['Entropy_raw'], summary['Entropy_covariance'])
    np.testing.assert_allclose(summary['Entropy'], summary['Entropy_covariance'])
    assert (summary['Entropy_n_adjusted_subjects'] == 0).all()
    assert (summary['Entropy_n_high_confidence_outliers'] == 0).all()
    assert (summary['Entropy_n_strong_qc_flags'] == summary['n']).all()
    assert (summary['Calculation_version'] == SOFTWARE_VERSION).all()


@pytest.fixture(scope='module', params=R.CASES)
def cross_version_case(request, tmp_path_factory):
    result, *_ = R.run_case(tmp_path_factory.mktemp('metric_reference'), request.param)
    return request.param, result


def _legacy_frame(case, table):
    data = json.loads(Path(__file__).with_name('reference_metric_baseline.json').read_text(encoding='utf-8'))
    snap = data['cases'][case][table]
    return pd.DataFrame(snap['rows'], columns=snap['columns'])


def test_public_results_route_both_distances_against_reference(cross_version_case):
    case, results = cross_version_case
    sample = results['tables']['样本级结果']
    old = _legacy_frame(case, '样本级结果')
    assert sample[['Experiment', 'Week', 'Group', 'Sample ID']].values.tolist() == old[['Experiment', 'Week', 'Group', 'Sample ID']].values.tolist()
    for current, previous in [('NMD', 'NMD_covariance_raw'), ('NMD_diagonal_capped', 'NMD'), ('BRI', 'BRI'), ('NRBS', 'NRBS')]:
        np.testing.assert_allclose(pd.to_numeric(sample[current]), pd.to_numeric(old[previous]), rtol=1e-10, atol=1e-12, equal_nan=True)
    summary = results['tables']['组时间点汇总']
    old_summary = _legacy_frame(case, '组时间点汇总')
    np.testing.assert_allclose(summary['Entropy_covariance'], pd.to_numeric(old_summary['Entropy_raw']), rtol=1e-10)
    np.testing.assert_allclose(summary['Entropy_covariance'], summary['Entropy_raw'])
    np.testing.assert_allclose(summary['Entropy_covariance'], summary['Entropy_robust_covariance'])
    for current, previous in [('NMD_diagonal_capped_mean_ci_low', 'NMD_mean_ci_low'),
                              ('NMD_diagonal_capped_mean_ci_high', 'NMD_mean_ci_high')]:
        np.testing.assert_allclose(summary[current], pd.to_numeric(old_summary[previous]), rtol=1e-10, atol=1e-12, equal_nan=True)


def test_nmd_downstream_and_sensitivity_use_correct_samples(cross_version_case):
    case, results = cross_version_case
    if case != 'nmd_longitudinal':
        assert 'NMD方法敏感性' not in results['tables']
        return
    sample = results['tables']['样本级结果']
    summary = results['tables']['组时间点汇总']
    for _, row in summary.iterrows():
        sub = sample[(sample['Week'] == row['Week']) & (sample['Group'] == row['Group'])]
        assert row['NMD_mean'] == pytest.approx(sub['NMD'].mean())
        assert row['NMD_diagonal_capped_mean'] == pytest.approx(sub['NMD_diagonal_capped'].mean())
        assert row['Entropy'] == pytest.approx(.5*np.log(np.mean(sub['NMD'].values**2)) + .5*np.log(2*np.pi*np.e))
    assert len(results['tables']['NMD方法敏感性']) == len(summary)
    for _, row in results['tables']['治疗比较'].iterrows():
        week = sample[sample['Week'] == row['Week']]
        t, c = week[week['Group'] == 'Treatment'], week[week['Group'] == 'Comparator']
        assert row['NMD_diff_treat_minus_comp'] == pytest.approx(t['NMD'].mean() - c['NMD'].mean())
        assert row['NMD_diagonal_capped_diff_treat_minus_comp'] == pytest.approx(t['NMD_diagonal_capped'].mean() - c['NMD_diagonal_capped'].mean())
        assert row['NMD_p_mannwhitney'] == pytest.approx(mannwhitneyu(t['NMD'], c['NMD']).pvalue)
        assert row['NMD_diagonal_capped_p_mannwhitney'] == pytest.approx(mannwhitneyu(t['NMD_diagonal_capped'], c['NMD_diagonal_capped']).pvalue)


def test_desktop_version_settings_and_generated_guide_are_distinct(tmp_path):
    import BioEntropy_Desktop_App as D
    import bioentropy_core_generalized as G
    app = D.BioEntropyDesktopApp.__new__(D.BioEntropyDesktopApp)
    assert D.APP_VERSION == SOFTWARE_VERSION == 'V1.0'
    assert app._settings_path().name == '.bioentropy_v1_0_desktop.json'
    guide = G.build_input_guide_html(tmp_path / 'guide.html').read_text(encoding='utf-8')
    assert 'Ledoit–Wolf 全协方差 NMD' in guide
    assert '取消额外 0.5 个体收缩' in guide
    assert '输出稳健封顶 NMD' not in guide
