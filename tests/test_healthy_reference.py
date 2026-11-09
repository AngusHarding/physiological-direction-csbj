"""Locks for the independently rebuilt availability-aware reference."""

from pathlib import Path
import unittest

import pandas as pd

ROOT = Path(__file__).resolve().parents[1] / "1_prepare_nhanes"  # public_code: the preparation stage


class HealthyReferenceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.flags = pd.read_csv(ROOT / "outputs/08_reference/healthy_reference_flags.csv")

    def test_primary_reference_counts(self):
        ref = self.flags[self.flags.is_reference_primary]
        self.assertEqual(len(ref), 1238)
        self.assertEqual(int(ref.RIAGENDR.eq(1).sum()), 547)
        self.assertEqual(int(ref.RIAGENDR.eq(2).sum()), 691)

    def test_historical_logic_is_independently_reproduced(self):
        self.assertEqual(int(self.flags.is_reference_historical_compatible.sum()), 1241)

    def test_primary_excludes_only_three_genuine_unknowns_from_historical(self):
        historical_only = (self.flags.is_reference_historical_compatible
                           & ~self.flags.is_reference_primary)
        self.assertEqual(int(historical_only.sum()), 3)
        self.assertTrue(self.flags.loc[historical_only, "primary_first_exclusion"].eq(
            "disease_unknown").all())

    def test_primary_is_subset_of_historical_compatible(self):
        self.assertFalse((self.flags.is_reference_primary
                         & ~self.flags.is_reference_historical_compatible).any())

    def test_reference_has_sufficient_sex_specific_ratio(self):
        ref = self.flags[self.flags.is_reference_primary]
        self.assertGreater(int(ref.RIAGENDR.eq(1).sum()) / 38, 10)
        self.assertGreater(int(ref.RIAGENDR.eq(2).sum()) / 38, 10)


if __name__ == "__main__":
    unittest.main()
