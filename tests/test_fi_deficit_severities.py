"""Tests for bounded FI deficit-severity grading."""

import importlib.util
from pathlib import Path
import unittest

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1] / "1_prepare_nhanes"  # public_code: the preparation stage
MODULE_PATH = ROOT / "8_adverse_distance/05_build_deficit_severities.py"  # public_code
MODULE_SPEC = importlib.util.spec_from_file_location("build_fi_deficit_severities", MODULE_PATH)
module = importlib.util.module_from_spec(MODULE_SPEC)
MODULE_SPEC.loader.exec_module(module)


class TestFIDeficitSeverities(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.excursions = pd.read_csv(module.EXCURSIONS)
        cls.spec = pd.read_csv(module.SCORING_SPEC)
        cls.markers = cls.spec.marker.tolist()
        cls.severities = module.build_severities(cls.excursions, cls.spec)

    def test_severity_is_the_excursion(self):
        """Version 2: no dead zone, no cap, no free parameters."""
        excursions = np.array([0.0, 0.5, 1.0, 1.5, 3.0, 8.0, 17.7])
        np.testing.assert_allclose(module.grade_excursion(excursions), excursions)

    def test_no_dead_zone(self):
        """Every non-zero excursion must score strictly positive.

        Regression test for FI_METHODS.md Amendment 1. Under the version-1 rule
        every one of these excursions scored exactly zero.
        """
        excursions = np.array([1e-9, 0.01, 0.25, 0.5, 0.75, 0.999, 1.0])
        self.assertTrue((module.grade_excursion(excursions) > 0).all(), "dead zone reintroduced")

    def test_no_cap(self):
        """Distinct excursions above the former cap must stay distinct.

        Under the version-1 rule all of these collapsed to exactly 1.0, which
        flattened 68.7% of waist-to-height values.
        """
        excursions = np.array([3.0, 4.0, 8.0, 15.4, 17.7])
        graded = module.grade_excursion(excursions)
        self.assertEqual(len(np.unique(graded)), len(excursions), "cap reintroduced")
        self.assertTrue((np.diff(graded) > 0).all(), "grading is not strictly increasing")

    def test_proportionality(self):
        """Doubling the excursion doubles the deficit, at every magnitude."""
        for base in (0.3, 0.7, 2.0, 6.0):
            self.assertAlmostEqual(
                float(module.grade_excursion(np.array([2 * base]))[0]),
                2 * float(module.grade_excursion(np.array([base]))[0]),
            )

    def test_version_one_rule_still_reproducible_for_sensitivity(self):
        """The superseded bounded rule must remain available as a sensitivity."""
        excursions = np.array([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0])
        expected = np.array([0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0])
        np.testing.assert_allclose(
            module.grade_excursion(excursions, dead_zone=1.0, cap=3.0), expected
        )

    def test_direct_matrix_formula(self):
        expected = self.excursions[self.markers].to_numpy(float)
        np.testing.assert_allclose(self.severities[self.markers].to_numpy(float), expected)

    def test_structure_bounds_and_identity(self):
        self.assertEqual(self.severities.shape, (29053, 42))
        self.assertTrue(self.severities.SEQN.equals(self.excursions.SEQN))
        values = self.severities[self.markers].to_numpy(float)
        self.assertTrue(np.isfinite(values).all())
        self.assertGreaterEqual(values.min(), 0.0)
        # Version 2: non-negative but unbounded above. Assert the maximum is NOT
        # a ceiling, which is what distinguishes this from the version-1 rule.
        self.assertGreater(values.max(), 1.0)
        self.assertAlmostEqual(
            values.max(), self.excursions[self.markers].to_numpy(float).max()
        )


if __name__ == "__main__":
    unittest.main()
