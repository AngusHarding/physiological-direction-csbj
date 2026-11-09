"""Tests for the prospective full-cohort transformation selector."""

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import unittest

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1] / "1_prepare_nhanes"  # public_code: the preparation stage
MODULE_PATH = ROOT / "7_marker_panel/02_select_transformations.py"  # public_code
module_spec = spec_from_file_location("transform_selector", MODULE_PATH)
selector = module_from_spec(module_spec)
module_spec.loader.exec_module(selector)


class TransformationFunctionTest(unittest.TestCase):
    def test_candidate_formulas(self):
        x = np.array([1.0, 8.0, 27.0])
        np.testing.assert_allclose(selector.transform_values(x, "identity"), x)
        np.testing.assert_allclose(selector.transform_values(x, "cbrt"), [1, 2, 3])
        np.testing.assert_allclose(selector.transform_values(x, "sqrt"), np.sqrt(x))
        np.testing.assert_allclose(selector.transform_values(x, "log"), np.log(x))
        np.testing.assert_allclose(selector.transform_values(x, "reciprocal"), 1 / x)

    def test_zero_is_not_artificially_shifted(self):
        x = np.array([0.0, 1.0, 4.0])
        self.assertIsNotNone(selector.transform_values(x, "sqrt"))
        self.assertIsNone(selector.transform_values(x, "log"))
        self.assertIsNone(selector.transform_values(x, "reciprocal"))

    def test_selector_is_deterministic(self):
        x = np.arange(1.0, 101.0) ** 2
        first = selector.evaluate_marker(x, "example", retain_identity=False)
        second = selector.evaluate_marker(x, "example", retain_identity=False)
        pd.testing.assert_frame_equal(first, second)
        self.assertEqual(int(first.selected.sum()), 1)

    def test_identity_is_locked_when_gate_passes(self):
        x = np.arange(1.0, 101.0) ** 2
        ledger = selector.evaluate_marker(x, "example", retain_identity=True)
        self.assertEqual(ledger.loc[ledger.selected, "candidate"].iloc[0], "identity")
        self.assertEqual(ledger.selection_basis.iloc[0], "raw_identity_passed_all_three")


class TransformationOutputTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.raw = pd.read_csv(ROOT / "outputs/05_panel/final_38_biomarker_panel_raw.csv")
        cls.spec = pd.read_csv(ROOT / "config/final_38_panel_spec.csv")
        cls.selection_path = ROOT / "outputs/06_transform/transformation_selection_v2.csv"
        cls.transformed_path = ROOT / "outputs/06_transform/final_38_panel_transformed_unstandardised_v2.csv"

    def test_outputs_exist_and_have_one_choice_per_marker(self):
        self.assertTrue(self.selection_path.exists())
        selection = pd.read_csv(self.selection_path)
        self.assertListEqual(selection.marker.tolist(), self.spec.marker.tolist())
        self.assertEqual(selection.marker.nunique(), 38)

    def test_transformed_panel_preserves_identity_and_shape(self):
        transformed = pd.read_csv(self.transformed_path)
        self.assertEqual(transformed.shape, self.raw.shape)
        self.assertListEqual(list(transformed.columns), list(self.raw.columns))
        pd.testing.assert_series_equal(transformed.SEQN, self.raw.SEQN)
        self.assertTrue(np.isfinite(transformed.iloc[:, 3:].to_numpy()).all())

    def test_saved_values_match_declared_formulas(self):
        transformed = pd.read_csv(self.transformed_path)
        selection = pd.read_csv(self.selection_path).set_index("marker")
        for marker in self.spec.marker:
            expected = selector.transform_values(
                self.raw[marker].to_numpy(), selection.loc[marker, "selected_transform"])
            np.testing.assert_allclose(transformed[marker], expected, rtol=1e-12, atol=1e-12)

    def test_two_step_counts_are_locked(self):
        selection = pd.read_csv(self.selection_path)
        self.assertEqual(int(selection.selection_basis.eq("raw_identity_passed_all_three").sum()), 13)
        self.assertEqual(int(selection.selection_basis.eq("least_aggressive_adequate").sum()), 13)
        self.assertEqual(int(selection.selection_basis.eq("best_available_not_adequate").sum()), 12)


if __name__ == "__main__":
    unittest.main()
