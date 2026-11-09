from pathlib import Path
import unittest
import pandas as pd

ROOT = Path(__file__).resolve().parents[1] / "1_prepare_nhanes"  # public_code: the preparation stage


class HeiTests(unittest.TestCase):
    def test_total_output_is_unique_and_bounded(self):
        total = pd.read_csv(ROOT / "outputs/02_harmonised/hei_2020_total.csv")
        self.assertEqual(len(total), 49_621)
        self.assertTrue(total.SEQN.is_unique)
        self.assertTrue(total.HEI_2020_total.between(0, 100).all())

    def test_component_schema_and_sum(self):
        frame = pd.read_csv(ROOT / "outputs/02_harmonised/hei_2020_component_scores.csv")
        scores = [c for c in frame if c.startswith("score_")]
        self.assertEqual(len(scores), 13)
        complete = frame.missing_component_count.eq(0)
        difference = frame.loc[complete, scores].sum(axis=1) - frame.loc[complete, "HEI_2020_total"]
        self.assertLess(difference.abs().max(), 1e-10)
        self.assertEqual(int(frame.missing_component_count.gt(0).sum()), 8)
        self.assertTrue(frame.loc[frame.missing_component_count.gt(0), "AGE"].le(3).all())


if __name__ == "__main__":
    unittest.main()

