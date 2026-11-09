"""Locks for the raw fixed-width linked-mortality reconstruction."""

from pathlib import Path
import unittest

import pandas as pd

ROOT = Path(__file__).resolve().parents[1] / "1_prepare_nhanes"  # public_code: the preparation stage


class MortalityTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mortality = pd.read_csv(ROOT / "outputs/01_imported/mortality.csv")

    def test_locked_rows_deaths_and_identity(self):
        self.assertEqual(len(self.mortality), 59064)
        self.assertEqual(int(self.mortality.mortstat.eq(1).sum()), 9249)
        self.assertEqual(int(self.mortality.mortstat.eq(0).sum()), 49815)
        self.assertFalse(self.mortality.SEQN.duplicated().any())

    def test_lossless_fields_are_retained(self):
        expected = ["SEQN", "mortstat", "ucod_leading", "ucod_label",
                    "diabetes_mcod", "hyperten_mcod", "permth_int",
                    "permth_exm", "cycle"]
        self.assertListEqual(list(self.mortality.columns), expected)

    def test_accidental_deaths_are_not_deleted(self):
        accidents = self.mortality.mortstat.eq(1) & self.mortality.ucod_leading.eq(4)
        self.assertGreater(int(accidents.sum()), 0)

    def test_alive_cause_is_missing_not_synthetic_code(self):
        alive = self.mortality.mortstat.eq(0)
        self.assertTrue(self.mortality.loc[alive, "ucod_leading"].isna().all())
        self.assertTrue(self.mortality.loc[alive, "ucod_label"].eq("Alive").all())

    def test_both_followup_clocks_are_preserved(self):
        self.assertIn("permth_int", self.mortality.columns)
        self.assertIn("permth_exm", self.mortality.columns)
        self.assertEqual(int(self.mortality.permth_int.isna().sum()), 0)
        # Public files contain unavailable exam clocks; these remain missing and
        # are handled by the downstream survival-analysis eligibility rule.
        self.assertEqual(int(self.mortality.permth_exm.isna().sum()), 2811)


if __name__ == "__main__":
    unittest.main()
