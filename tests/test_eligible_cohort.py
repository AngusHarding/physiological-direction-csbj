from pathlib import Path
import unittest

import pandas as pd


ROOT = Path(__file__).resolve().parents[1] / "1_prepare_nhanes"  # public_code: the preparation stage


class EligibleCohortTests(unittest.TestCase):
    """Protect the certified cohort boundary against accidental code changes."""

    @classmethod
    def setUpClass(cls):
        # Load large CSVs once for the whole class rather than before every test.
        cls.cohort = pd.read_csv(ROOT / "outputs" / "04_cohorts" / "eligible_adult_cohort.csv")
        cls.flags = pd.read_csv(ROOT / "outputs" / "04_cohorts" / "eligibility_flags.csv")

    def test_certified_counts(self):
        self.assertEqual(len(self.cohort), 34785)
        self.assertEqual(int((self.cohort.RIAGENDR == 1).sum()), 17924)
        self.assertEqual(int((self.cohort.RIAGENDR == 2).sum()), 16861)

    def test_unique_complete_key(self):
        self.assertFalse(self.cohort.SEQN.isna().any())
        self.assertFalse(self.cohort.SEQN.duplicated().any())

    def test_all_output_participants_are_flagged_eligible(self):
        eligible_ids = set(self.flags.loc[self.flags.eligible_adult, "SEQN"])
        self.assertEqual(set(self.cohort.SEQN), eligible_ids)


if __name__ == "__main__":
    unittest.main()
