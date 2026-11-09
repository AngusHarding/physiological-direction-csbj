"""Tests for sex-specific FI healthy-reference parameter estimation."""

import importlib.util
from pathlib import Path
import unittest

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1] / "1_prepare_nhanes"  # public_code: the preparation stage
MODULE_PATH = ROOT / "8_adverse_distance/02_estimate_reference_parameters.py"  # public_code
MODULE_SPEC = importlib.util.spec_from_file_location("estimate_fi_reference_parameters", MODULE_PATH)
module = importlib.util.module_from_spec(MODULE_SPEC)
MODULE_SPEC.loader.exec_module(module)


class TestFIReferenceParameters(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.transformed = pd.read_csv(module.TRANSFORMED)
        cls.membership = pd.read_csv(module.MEMBERSHIP)
        cls.spec = pd.read_csv(module.SCORING_SPEC)
        cls.parameters = module.estimate_parameters(cls.transformed, cls.membership, cls.spec)

    def test_parameter_structure_and_scales(self):
        self.assertEqual(self.parameters.shape, (76, 15))
        self.assertFalse(self.parameters.duplicated(["marker", "sex_code"]).any())
        self.assertTrue(np.isfinite(self.parameters.reference_robust_scale).all())
        self.assertTrue(self.parameters.reference_robust_scale.gt(0).all())

    def test_direct_formula_reproduction(self):
        marker, sex_code = "LBXSATSI", 1
        reference = (
            self.membership.is_fi_reference.astype(bool)
            & self.membership.RIAGENDR.eq(sex_code)
        )
        values = self.transformed.loc[reference, marker].to_numpy(float)
        q1, median, q3 = np.quantile(values, [0.25, 0.50, 0.75], method="linear")
        observed = self.parameters.loc[
            self.parameters.marker.eq(marker) & self.parameters.sex_code.eq(sex_code)
        ].iloc[0]
        self.assertAlmostEqual(observed.reference_q1_transformed, q1, places=14)
        self.assertAlmostEqual(observed.reference_median_transformed, median, places=14)
        self.assertAlmostEqual(observed.reference_q3_transformed, q3, places=14)
        self.assertAlmostEqual(observed.reference_robust_scale, (q3 - q1) / 1.349, places=14)

    def test_metadata_matches_scoring_spec(self):
        first_sex = self.parameters[self.parameters.sex_code.eq(1)].reset_index(drop=True)
        self.assertEqual(first_sex.marker.tolist(), self.spec.marker.tolist())
        self.assertEqual(first_sex.orientation.tolist(), self.spec.orientation.tolist())
        self.assertEqual(first_sex.adverse_direction.tolist(), self.spec.adverse_direction.tolist())


if __name__ == "__main__":
    unittest.main()
