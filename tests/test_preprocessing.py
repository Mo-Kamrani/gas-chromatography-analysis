"""Numerical checks for signed baselines, disjoint areas, and exported decisions."""

from copy import deepcopy
from pathlib import Path
from tempfile import TemporaryDirectory
import ast
import csv
import json
import unittest

import numpy as np
from scipy.integrate import trapezoid
from scipy.signal import find_peaks

from test_notebook import ROOT, NOTEBOOK, namespace as functions


def default_settings():
    for cell in NOTEBOOK['cells']:
        if cell['cell_type'] != 'code':
            continue
        for node in ast.parse(''.join(cell['source'])).body:
            if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'PREPROCESSING_SETTINGS' for t in node.targets):
                return ast.literal_eval(node.value)
    raise AssertionError('Missing preprocessing settings')


class PreprocessingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.settings = default_settings()
        cls.time, cls.raw = functions['extract_chromatogram'](functions['read_file'](ROOT / 'data/Gaschroma.txt'))
        cls.result = functions['preprocess_chromatogram'](cls.time, cls.raw, deepcopy(cls.settings))

    def test_supplied_export_and_signed_measurements(self):
        result, summary = self.result, self.result['summary']
        self.assertEqual(summary['input_points'], 25500)
        self.assertEqual(summary['peak_candidate_count'], 32)
        self.assertAlmostEqual(summary['candidate']['peak_time_min'], 1.70883)
        self.assertAlmostEqual(summary['first_selected_time_min'], 2.63817)
        self.assertEqual(summary['negative_values']['raw_full_negative_count'], 3070)
        self.assertEqual(summary['negative_values']['corrected_negative_count'], 7191)
        self.assertEqual(summary['baseline']['selected_point_count'], 3385)
        np.testing.assert_array_equal(result['raw_full'], self.raw)
        np.testing.assert_allclose(result['corrected'], result['raw_selected'] - result['baseline'])
        self.assertFalse(np.shares_memory(result['raw_full'], self.raw))
        for left, right in zip(result['integration_records'], result['integration_records'][1:]):
            self.assertLessEqual(left['integration_end_min'], right['integration_start_min'])

    def test_quadratic_baseline_and_transformed_coefficients(self):
        time = np.linspace(2, 12, 501)
        known = 4 + 0.4 * time - 0.03 * time**2
        baseline, info = functions['estimate_polynomial_baseline'](time, known, np.ones(len(time), dtype=bool), 2)
        np.testing.assert_allclose(baseline, known, atol=1e-12)
        scaled = (time - info['time_center_min']) / info['time_scale_min']
        np.testing.assert_allclose(np.polyval(info['coefficients_highest_power_first'], scaled), known, atol=1e-12)

    def test_isolated_gaussian_area_with_finite_boundaries(self):
        time = np.linspace(0, 8, 16001)
        signal = 100 * np.exp(-0.5 * ((time - 4) / 0.1)**2)
        regions, _ = functions['define_integration_regions'](time, signal, [8000], 0.001, self.settings)
        analytical = 100 * 0.1 * np.sqrt(2 * np.pi)
        for method in ('trapezoid', 'simpson'):
            with self.subTest(method=method):
                result = functions['integrate_peak_regions'](time, signal, regions, method)[0]
                self.assertLess(abs(result['area_signal_min'] - analytical) / analytical, 0.0025)
                self.assertAlmostEqual(result['area_signal_s'], result['area_signal_min'] * 60)

    def test_overlapping_peak_partitions_do_not_duplicate_area(self):
        time = np.linspace(0, 10, 20001)
        signal = 100 * (np.exp(-0.5*((time-4)/0.55)**2) + np.exp(-0.5*((time-5.3)/0.55)**2))
        peaks, _ = find_peaks(signal, prominence=1)
        self.assertEqual(len(peaks), 2)
        regions, valleys = functions['define_integration_regions'](time, signal, peaks, 0.001, self.settings)
        self.assertTrue(valleys[0]['overlap_review_required'])
        self.assertEqual(regions[0]['integration_end_min'], regions[1]['integration_start_min'])
        integrated = functions['integrate_peak_regions'](time, signal, regions, 'trapezoid')
        start, end = regions[0]['integration_start_min'], regions[-1]['integration_end_min']
        times = np.r_[start, time[(time > start) & (time < end)], end]
        union = trapezoid(np.interp(times, time, signal), x=times)
        self.assertAlmostEqual(sum(row['area_signal_min'] for row in integrated), union, places=10)

    def test_signed_negative_area_and_invalid_regions(self):
        time = np.linspace(0, 1, 101)
        signal = np.full(101, -2.0)
        region = {'integration_start_min': 0.25, 'integration_end_min': 0.75}
        for method in ('trapezoid', 'simpson'):
            self.assertAlmostEqual(functions['integrate_peak_regions'](time, signal, [region], method)[0]['area_signal_min'], -1.0)
        for regions in ([{'integration_start_min': -1, 'integration_end_min': 0.5}], [region, region]):
            with self.assertRaisesRegex(ValueError, 'inside the data range'):
                functions['integrate_peak_regions'](time, signal, regions, 'trapezoid')

    def test_manual_baseline_anchors_and_missing_selection(self):
        time = np.linspace(0, 6, 601)
        settings = {**self.settings, 'baseline_point_selection': 'manual', 'baseline_anchor_times_min': [0, 2, 4, 6]}
        mask, _ = functions['select_baseline_points'](time, time**2, settings)
        np.testing.assert_array_equal(time[mask], [0, 2, 4, 6])
        with self.assertRaisesRegex(ValueError, 'needs baseline intervals or anchor times'):
            functions['select_baseline_points'](time, time**2, {**settings, 'baseline_anchor_times_min': []})

    def test_sampling_and_detection_copy(self):
        time = np.linspace(0, 1, 101)
        raw = np.sin(time)
        time_copy, raw_copy, step = functions['validate_signal'](time, raw, 0.05)
        functions['smooth_for_detection'](raw_copy, step, 0.05, 2)
        np.testing.assert_array_equal(raw, np.sin(time))
        bad = time.copy()
        bad[50] += 0.005
        with self.assertRaisesRegex(ValueError, 'irregular'):
            functions['validate_signal'](bad, raw, 0.05)
        with self.assertRaises(ValueError):
            functions['validate_signal'](time, raw, float('nan'))

    def test_unconfirmed_recovery_requires_manual_start(self):
        time = np.linspace(0, 5, 501)
        candidate = {'found': True, 'recovered': False}
        with self.assertRaisesRegex(ValueError, 'no validated return'):
            functions['select_analysis_region'](time, candidate, self.settings)
        mask, selection = functions['select_analysis_region'](time, candidate, {**self.settings, 'analysis_start_min': 2.0})
        self.assertEqual(selection['start_source'], 'manual')
        self.assertGreaterEqual(time[mask][0], 2)

    def test_asls_alternative_and_constant_negative_input(self):
        time = np.linspace(0, 6, 601)
        raw = np.full(len(time), -5.0)
        settings = {**self.settings, 'detect_early_peak': False, 'baseline_method': 'asls'}
        result = functions['preprocess_chromatogram'](time, raw, settings)
        self.assertTrue(result['summary']['baseline']['converged'])
        np.testing.assert_allclose(result['baseline'], raw, atol=1e-10)
        self.assertEqual(result['summary']['peak_candidate_count'], 0)
        np.testing.assert_array_equal(result['raw_full'], raw)

    def test_summary_keeps_a_snapshot_of_settings(self):
        settings = deepcopy(self.settings)
        result = functions['preprocess_chromatogram'](self.time, self.raw, settings, run_peak_preview=False)
        settings['polynomial_degree'] = 99
        self.assertEqual(result['summary']['settings']['polynomial_degree'], 3)

    def test_signal_export_retains_original_points_and_blank_unselected_fields(self):
        with TemporaryDirectory() as directory:
            paths = functions['save_preprocessing_results'](self.result, directory, 'check')
            with paths[0].open(newline='') as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(len(rows), len(self.raw))
            np.testing.assert_array_equal([float(row['raw_intensity']) for row in rows], self.raw)
            self.assertEqual(rows[0]['selected_for_analysis'], '0')
            self.assertEqual(rows[0]['corrected_intensity_signed'], '')
            summary = json.loads(paths[-1].read_text())
            self.assertEqual(summary['integration']['integrated_peak_count'], 32)


if __name__ == '__main__':
    unittest.main()
