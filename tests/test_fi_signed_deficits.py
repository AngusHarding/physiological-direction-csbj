"""Tests for FI signed graded deficits."""

import importlib.util
from pathlib import Path
import unittest

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1] / "1_prepare_nhanes"  # public_code: the preparation stage
MODULE_PATH = ROOT / "8_adverse_distance/06_build_signed_deficits.py"  # public_code
MODULE_SPEC = importlib.util.spec_from_file_location("build_fi_signed_deficits", MODULE_PATH)
module = importlib.util.module_from_spec(MODULE_SPEC)
MODULE_SPEC.loader.exec_module(module)


class TestFISignedDeficits(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.coordinates = pd.read_csv(module.COORDINATES)
        cls.severities = pd.read_csv(module.SEVERITIES)
        cls.spec = pd.read_csv(module.SCORING_SPEC)
        cls.markers = cls.spec.marker.tolist()
        cls.signed = module.build_signed_deficits(
            cls.coordinates, cls.severities, cls.spec
        )

    def test_pure_sign_rules(self):
        z = np.array([-2.0, 0.0, 2.0])
        d = np.array([0.5, 0.0, 0.5])
        np.testing.assert_array_equal(module.signed_deficit(z, d, "high_only"), d)
        np.testing.assert_array_equal(module.signed_deficit(z, d, "low_only"), -d)
        np.testing.assert_array_equal(
            module.signed_deficit(z, d, "two_sided"), [-0.5, 0.0, 0.5]
        )

    def test_marker_examples(self):
        np.testing.assert_array_equal(self.signed.systolic, self.severities.systolic)
        np.testing.assert_array_equal(self.signed.LBXSAL, -self.severities.LBXSAL)
        expected_calcium = np.sign(self.coordinates.LBXSCA) * self.severities.LBXSCA
        np.testing.assert_array_equal(self.signed.LBXSCA, expected_calcium)

    def test_magnitude_bounds_and_identity(self):
        """Version 2: unbounded magnitude; only |signed| == severity is invariant.

        The former [-1, 1] bound came from the version-1 deficit cap, removed as
        a specification defect (FI_METHODS.md Amendment 1).
        """
        g = self.signed[self.markers].to_numpy(float)
        d = self.severities[self.markers].to_numpy(float)
        self.assertTrue(np.isfinite(g).all())
        np.testing.assert_array_equal(np.abs(g), d)
        # Positive check that the cap is gone rather than merely unasserted.
        self.assertGreater(np.abs(g).max(), 1.0)
        self.assertTrue(self.signed.SEQN.equals(self.coordinates.SEQN))


if __name__ == "__main__":
    unittest.main()
