"""Opt-in source/frozen runtime smoke check using public synthetic inputs only."""
from pathlib import Path
import json
import sys
import tempfile
import traceback
import zipfile

import numpy as np
import pandas as pd

from bioentropy_constants import SOFTWARE_VERSION, NMD_PRIMARY_METHOD


def run_self_test(report_path: str) -> int:
    """Exercise the shipped loader, analysis, XLSX/figure exports and audit."""
    destination = Path(report_path).resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    report = {'software_version': SOFTWARE_VERSION, 'frozen': bool(getattr(sys, 'frozen', False))}
    bundle = None
    try:
        import bioentropy_core_generalized as G
        from bioentropy_audit import capture_input_files, cleanup_generated_bundle
        with tempfile.TemporaryDirectory(prefix='bioentropy_selftest_') as folder:
            root = Path(folder)
            observations = np.array([[-1.2,-.3,.8],[-.8,.9,-.4],[-.4,-1.1,-.7],[-.1,.4,1.2],
                                     [.2,-.6,-1.],[.5,1.3,.1],[.8,-.8,.5],[1.,.2,-.3]])
            features = ['f0', 'f1', 'f2']
            paths = []
            for week in [0, 6]:
                rows = []
                for gi, group in enumerate(['Treatment', 'Comparator']):
                    shift = 1.5 if week == 0 else (.4 if gi == 0 else 1.2)
                    for i, point in enumerate(observations):
                        rows.append({'Group': group, 'Sample ID': gi*8+i+1,
                                     **dict(zip(features, point+shift+gi*.1))})
                path = root / f'SMOKE-{week}W.xlsx'
                pd.DataFrame(rows).to_excel(path, index=False)
                paths.append(path)
            normal = pd.DataFrame(np.tile(observations, (3,1)), columns=features)
            normal.insert(0, 'Sample ID', np.arange(100,124))
            path = root / 'SMOKE-Normal.xlsx'
            normal.to_excel(path, index=False)
            paths.append(path)
            uploads, provenance = capture_input_files(paths)
            raw, ranges, normals, *_ = G.load_input_bundle(uploads)
            roles = {'SMOKE': {'treatment': 'Treatment', 'comparator': 'Comparator'}}
            result = G.analyze_all(raw, role_map=roles, external_range_map=ranges, true_normal_map=normals)
            summary = result['tables']['组时间点汇总']
            sample = result['tables']['样本级结果']
            assert len(sample) == 32
            assert (sample['NMD_method'] == NMD_PRIMARY_METHOD).all()
            np.testing.assert_allclose(sample['NMD'], sample['NMD_covariance_raw'])
            np.testing.assert_allclose(summary['Entropy_covariance'], summary['Entropy_raw'])
            assert (summary['Entropy_n_adjusted_subjects'] == 0).all()
            original_summary = summary.copy(deep=True)
            bundle = G.generate_output_bundle(result, package_name='BioEntropySelfTest',
                audit_context={'inputs': provenance, 'role_map': roles, 'winsor_q': G.WINSOR_Q_DEFAULT,
                               'normal_reference_experiments': list(normals)})
            assert not bundle.get('audit_error')
            assert len(bundle['figure_paths']) > 0 and len(bundle['figure_data_files']) > 0
            assert len(bundle['figure_paths']) == 8
            assert not any('_main_09_' in p or '_main_10_' in p for p in bundle['figure_paths'])
            assert all(Path(p).is_file() for p in bundle['figure_paths'].values())
            pd.testing.assert_frame_equal(summary, original_summary, check_exact=True)
            change_files = [
                Path(bundle['output_dir']) / '图形数据_GraphPad' / 'SMOKE' / name
                for name in ['Fig02_NMD_baseline_aligned_normalized_bars.csv',
                             'Fig04_Entropy_baseline_aligned_normalized_bars.csv']
            ]
            for change_file in change_files:
                changes = pd.read_csv(change_file)
                baseline_rows = changes.loc[changes['Week'] == 0]
                assert len(baseline_rows) == 2
                assert (baseline_rows['Change_from_baseline'] == 0).all()
                assert (baseline_rows['Baseline_aligned_display_0_100_gamma3'] == 100).all()
                assert (baseline_rows['Baseline_week'] == 0).all()
            with zipfile.ZipFile(bundle['zip_path']) as z:
                assert z.testzip() is None
                assert any('NMD方法敏感性' in name for name in z.namelist())
            audit = json.loads((Path(bundle['audit_dir'])/'run_manifest.json').read_text(encoding='utf-8'))
            assert audit['software_version'] == SOFTWARE_VERSION
            assert audit['method_policy']['primary_nmd_z_cap'] is None
            report.update(status='passed', samples=len(sample), groups=len(summary),
                          figure_count=len(bundle['figure_paths']),
                          figure_data_files=len(bundle['figure_data_files']),
                          method_policy=audit['method_policy'],
                          environment=audit['environment'],
                          nmd_mean=summary['NMD_mean'].tolist(),
                          diagonal_mean=summary['NMD_diagonal_capped_mean'].tolist())
            report.update(entropy=summary['Entropy'].tolist(),
                          baseline_aligned_files=len(change_files),
                          baseline_starts_equal_nonzero=True,
                          output_preserves_summary=True,
                          display_revision=audit.get('code', {}).get('build_info', {}).get('display_revision'))
            cleanup_generated_bundle(bundle)
        status = 0
    except Exception:
        report.update(status='failed', error=traceback.format_exc())
        status = 1
    destination.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    return status
