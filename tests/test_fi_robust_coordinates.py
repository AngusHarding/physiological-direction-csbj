"""Tests for FI raw-oriented robust coordinates."""

import importlib.util
from pathlib import Path
import unittest

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1] / "1_prepare_nhanes"  # public_code: the preparation stage
MODULE_PATH = ROOT / "8_adverse_distance/03_build_robust_coordinates.py"  # public_code
MODULE_SPEC = importlib.util.spec_from_file_location("build_fi_robust_coordinates", MODULE_PATH)
module = importlib.util.module_from_spec(MODULE_SPEC)
MODULE_SPEC.loader.exec_module(module)


class TestFIRobustCoordinates(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.transformed = pd.read_csv(module.TRANSFORMED)
        cls.membership = pd.read_csv(module.MEMBERSHIP)
        cls.parameters = pd.read_csv(module.PARAMETERS)
        cls.spec = pd.read_csv(module.SCORING_SPEC)
        cls.markers = cls.spec.marker.tolist()
        cls.coordinates = module.build_coordinates(
            cls.transformed, cls.membership, cls.parameters, cls.spec
        )

    def assert_direct_formula(self, marker, row_index):
        sex_code = int(self.transformed.loc[row_index, "RIAGENDR"])
        parameter = self.parameters.loc[
            self.parameters.marker.eq(marker) & self.parameters.sex_code.eq(sex_code)
        ].iloc[0]
        expected = (
            int(parameter.orientation)
            * (self.transformed.loc[row_index, marker] - parameter.reference_median_transformed)
            / parameter.reference_robust_scale
        )
        self.assertAlmostEqual(self.coordinates.loc[row_index, marker], expected, places=13)

    def test_direct_formula_preserving_transform(self):
        self.assert_direct_formula("LBXSCA", 100)

    def test_direct_formula_reciprocal_transform(self):
        self.assert_direct_formula("LBXSATSI", 100)

    def test_reference_location_and_scale(self):
        checks = module.build_checks(self.coordinates, self.markers)
        self.assertLess(checks.maximum_absolute_reference_median.max(), 1e-12)
        self.assertLess(checks.maximum_absolute_reference_iqr_minus_1_349.max(), 1e-12)

    def test_structure_and_finite_values(self):
        self.assertEqual(self.coordinates.shape, (29053, 42))
        self.assertTrue(np.isfinite(self.coordinates[self.markers].to_numpy(float)).all())
        self.assertTrue(self.coordinates.SEQN.equals(self.transformed.SEQN))


if __name__ == "__main__":
    unittest.main()
