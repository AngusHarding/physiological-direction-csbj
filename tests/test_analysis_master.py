"""Structural locks for the inclusion-preserving analysis master."""

from pathlib import Path
import unittest

import pandas as pd

ROOT = Path(__file__).resolve().parents[1] / "1_prepare_nhanes"  # public_code: the preparation stage
MASTER = ROOT / "outputs/07_master/analysis_master_raw.csv"
RAW = ROOT / "outputs/05_panel/final_38_biomarker_panel_raw.csv"


class AnalysisMasterTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.master = pd.read_csv(MASTER)
        cls.raw = pd.read_csv(RAW)

    def test_master_preserves_locked_participants(self):
        self.assertEqual(len(self.master), 29053)
        self.assertFalse(self.master.SEQN.duplicated().any())
        pd.testing.assert_series_equal(self.master.SEQN, self.raw.SEQN)

    def test_raw_panel_is_preserved_exactly(self):
        pd.testing.assert_frame_equal(self.master[self.raw.columns], self.raw)

    def test_every_participant_passed_shared_eligibility(self):
        required = ["eligible_demographic_complete", "eligible_recent_illness_matched",
                    "eligible_recent_illness_complete", "eligible_age_18_79",
                    "eligible_not_acute", "eligible_not_pregnant", "eligible_adult"]
        self.assertFalse(self.master[required].isna().any().any())
        self.assertTrue(self.master[required].all().all())

    def test_optional_covariates_do_not_define_master_membership(self):
        # These missing values must remain visible rather than causing silent
        # complete-case deletion during the master merge.
        self.assertGreater(int(self.master.Smoking_status.isna().sum()), 0)
        self.assertGreater(int(self.master.drinks_per_week.isna().sum()), 0)
        self.assertGreater(int(self.master.mortstat.isna().sum()), 0)

    def test_reference_is_not_silently_inherited(self):
        self.assertNotIn("is_reference", self.master.columns)


if __name__ == "__main__":
    unittest.main()
