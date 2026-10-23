"""Regression checks for actual export parsing and misleading comparison cases."""

from pathlib import Path
from tempfile import TemporaryDirectory
import ast
import json
import unittest

import matplotlib
matplotlib.use('Agg')
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK = json.loads((ROOT / 'gas_chromatography.ipynb').read_text())
# Exercise the notebook's functions directly without rerunning plot cells on import.
namespace = {}
for cell in NOTEBOOK['cells']:
    if cell['cell_type'] != 'code':
        continue
    module = ast.parse(''.join(cell['source']))
    module.body = [node for node in module.body if isinstance(node, (ast.Import, ast.ImportFrom, ast.FunctionDef))]
    exec(compile(module, 'gas_chromatography.ipynb', 'exec'), namespace)


class AnalysisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.lines = namespace['read_file'](ROOT / 'data/Gaschroma.txt')
        cls.time, cls.signal = namespace['extract_chromatogram'](cls.lines)
        cls.compounds = namespace['extract_compound_results'](cls.lines)

    def test_complete_reference_export(self):
        self.assertEqual(len(self.time), 25500)
        np.testing.assert_allclose(self.time[[0, -1]], [0.00017, 16.99950])
        self.assertTrue(np.any(self.signal < 0))
        self.assertEqual(np.count_nonzero((self.time >= 4) & (self.time <= 16)), 18000)
        self.assertEqual(len(self.compounds), 35)
        self.assertEqual(sum(row['reported_percent'] > 0 for row in self.compounds), 31)
        self.assertAlmostEqual(sum(row['reported_percent'] for row in self.compounds), 100.0)
        by_name = {row['name']: row for row in self.compounds}
        self.assertEqual(by_name['C18:1 cis/trans']['reported_percent'], 5.568)
        self.assertEqual(by_name['C24:0']['reported_percent'], 7.39)

    def test_incomplete_signal_is_rejected(self):
        lines = self.lines.copy()
        header = lines.index('R.Time (min)\tIntensity')
        del lines[header + 1]
        with self.assertRaisesRegex(ValueError, 'Expected 25500 signal points'):
            namespace['extract_chromatogram'](lines)

    def test_signal_header_and_count_are_required(self):
        for marker in ('R.Time (min)\tIntensity', '# of Points\t25500'):
            with self.subTest(marker=marker), self.assertRaises(ValueError):
                namespace['extract_chromatogram']([line for line in self.lines if line != marker])

    def test_nonfinite_and_duplicate_time_are_rejected(self):
        for time, signal in (([0, 1, 1], [1, 2, 3]), ([0, np.nan, 2], [1, 2, 3]), ([0, 1, 2], [1, np.inf, 3])):
            with self.subTest(time=time, signal=signal), self.assertRaises(ValueError):
                namespace['validate_chromatogram'](time, signal)

    def test_normalization_uses_window_maximum(self):
        time, normalized = namespace['normalize_signal']([0, 1, 2, 3], [1000, 2, 4, -1], 1, 3)
        np.testing.assert_array_equal(time, [1, 2, 3])
        np.testing.assert_allclose(normalized, [0.5, 1, -0.25])

    def test_invalid_or_partial_normalization_windows(self):
        for start, end in ((3, 1), (1, np.nan), (-1, 2), (1, 4), (1.2, 1.3)):
            with self.subTest(start=start, end=end), self.assertRaises(ValueError):
                namespace['normalize_signal']([0, 1, 2, 3], [0, 1, 2, 3], start, end)
        with self.assertRaisesRegex(ValueError, 'positive signal maximum'):
            namespace['normalize_signal']([0, 1, 2], [-2, -1, 0], 0, 2)

    def test_bad_compound_rows_raise_clear_errors(self):
        header = self.lines.index(next(line for line in self.lines if line.startswith('ID#\tName')))
        for malformed in ('1\tC4:0', self.lines[header + 1].replace('0,000', 'NaN', 1)):
            lines = self.lines.copy()
            lines[header + 1] = malformed
            with self.subTest(malformed=malformed), self.assertRaises(ValueError):
                namespace['extract_compound_results'](lines)

    def test_duplicate_compounds_are_rejected(self):
        lines = self.lines.copy()
        header = lines.index(next(line for line in lines if line.startswith('ID#\tName')))
        lines[header + 2] = lines[header + 2].replace('C6:0', 'C4:0')
        with self.assertRaisesRegex(ValueError, 'unique'):
            namespace['extract_compound_results'](lines)

    def test_missing_values_are_distinct_from_explicit_zero(self):
        samples = {
            'a': {'compounds': [{'name': 'x', 'reported_percent': 0, 'area': 0}, {'name': 'y', 'reported_percent': 100, 'area': 10}]},
            'b': {'compounds': [{'name': 'x', 'reported_percent': 100, 'area': 20}]},
        }
        _, names, matrix, _ = namespace['build_compound_matrix'](samples)
        self.assertEqual(names, ['x', 'y'])
        self.assertEqual(matrix[0, 0], 0)
        self.assertTrue(np.isnan(matrix[1, 1]))
        _, _, area_matrix, _ = namespace['build_compound_matrix'](samples, 'area_percent')
        np.testing.assert_allclose(np.nansum(area_matrix, axis=1), [100, 100])

    def test_empty_profiles_and_invalid_metrics(self):
        with self.assertRaises(ValueError):
            namespace['build_compound_matrix']({})
        with self.assertRaises(ValueError):
            namespace['build_compound_matrix']({'a': {'compounds': []}})
        with self.assertRaises(ValueError):
            namespace['build_compound_matrix']({'a': {'compounds': self.compounds}}, 'unknown')

    def test_duplicate_file_copies_are_not_independent_samples(self):
        with TemporaryDirectory() as directory:
            first = Path(directory) / 'a.txt'
            second = Path(directory) / 'b.txt'
            raw = (ROOT / 'data/Gaschroma.txt').read_bytes()
            first.write_bytes(raw)
            second.write_bytes(raw)
            with self.assertRaisesRegex(ValueError, 'Byte-identical'):
                namespace['load_samples']({'a': first, 'b': second}, 'utf-8-sig', '[Chromatogram (Ch1)]', '[Compound Results(Ch1)]')

    def test_top_n_must_be_integer(self):
        for top_n in (0, -1, 1.5, True):
            with self.subTest(top_n=top_n), self.assertRaisesRegex(ValueError, 'positive integer'):
                namespace['plot_compound_bars'](['a'], ['x'], np.array([[1.0]]), 'Value', top_n)

    def test_plotting_preserves_raw_arrays(self):
        original_time, original_signal = self.time.copy(), self.signal.copy()
        figure = namespace['plot_chromatogram'](self.time, self.signal, 4, 16)
        namespace['plt'].close(figure)
        np.testing.assert_array_equal(self.time, original_time)
        np.testing.assert_array_equal(self.signal, original_signal)


if __name__ == '__main__':
    unittest.main()
