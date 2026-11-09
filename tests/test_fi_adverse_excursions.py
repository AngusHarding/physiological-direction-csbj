"""Tests for direction-specific FI adverse excursions."""

import importlib.util
from pathlib import Path
import unittest

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1] / "1_prepare_nhanes"  # public_code: the preparation stage
MODULE_PATH = ROOT / "8_adverse_distance/04_build_adverse_excursions.py"  # public_code
MODULE_SPEC = importlib.util.spec_from_file_location("build_fi_adverse_excursions", MODULE_PATH)
module = importlib.util.module_from_spec(MODULE_SPEC)
MODULE_SPEC.loader.exec_module(module)


class TestFIAdverseExcursions(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.coordinates = pd.read_csv(module.COORDINATES)
        cls.spec = pd.read_csv(module.SCORING_SPEC)
        cls.markers = cls.spec.marker.tolist()
        cls.excursions = module.build_excursions(cls.coordinates, cls.spec)

    def test_pure_direction_rules(self):
        values = np.array([-2.0, 0.0, 3.0])
        np.testing.assert_array_equal(
            module.adverse_excursion(values, "high_only"), [0.0, 0.0, 3.0]
        )
        np.testing.assert_array_equal(
            module.adverse_excursion(values, "low_only"), [2.0, 0.0, 0.0]
        )
        np.testing.assert_array_equal(
            module.adverse_excursion(values, "two_sided"), [2.0, 0.0, 3.0]
        )

    def test_marker_examples(self):
        expected = {
            "systolic": np.maximum(0.0, self.coordinates.systolic),
            "LBXSAL": np.maximum(0.0, -self.coordinates.LBXSAL),
            "LBXSCA": np.abs(self.coordinates.LBXSCA),
        }
        for marker, values in expected.items():
            np.testing.assert_array_equal(self.excursions[marker], values)

    def test_structure_and_range(self):
        self.assertEqual(self.excursions.shape, (29053, 42))
        self.assertTrue(self.excursions.SEQN.equals(self.coordinates.SEQN))
        values = self.excursions[self.markers].to_numpy(float)
        self.assertTrue(np.isfinite(values).all())
        self.assertGreaterEqual(values.min(), 0.0)
        self.assertGreater(values.max(), 3.0)  # confirms this step has not capped excursions


if __name__ == "__main__":
    unittest.main()
