"""Structural locks for the authoritative untransformed 38-marker panel."""

from pathlib import Path
import unittest
import pandas as pd

ROOT = Path(__file__).resolve().parents[1] / "1_prepare_nhanes"  # public_code: the preparation stage
PANEL = ROOT / "outputs/05_panel/final_38_biomarker_panel_raw.csv"


class Final38PanelTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.panel = pd.read_csv(PANEL)
        cls.spec = pd.read_csv(ROOT / "config/final_38_panel_spec.csv")

    def test_locked_dimensions_and_sex_counts(self):
        self.assertEqual(self.panel.shape, (29053, 41))
        self.assertEqual(int(self.panel.RIAGENDR.eq(1).sum()), 14916)
        self.assertEqual(int(self.panel.RIAGENDR.eq(2).sum()), 14137)

    def test_column_order_follows_specification(self):
        expected = ["SEQN", "RIAGENDR", "RIDAGEYR"] + self.spec.marker.tolist()
        self.assertListEqual(list(self.panel.columns), expected)

    def test_panel_is_unique_complete_case(self):
        self.assertFalse(self.panel.SEQN.duplicated().any())
        self.assertFalse(self.panel.isna().any().any())

    def test_values_are_untransformed(self):
        # These anchors are clinically interpretable raw units. A transformed
        # or z-standardised table would fail them immediately.
        self.assertGreater(self.panel.LBXTC.median(), 100)       # mg/dL
        self.assertGreater(self.panel.systolic.median(), 80)    # mmHg
        self.assertGreater(self.panel.eGFR.median(), 50)        # mL/min/1.73m2

    def test_auxiliary_fields_are_not_panel_columns(self):
        excluded = {"BMI", "LBXSCR", "LBXSOSSI", "LBXHCT", "LBXMCHSI", "LBXRBCSI"}
        self.assertTrue(excluded.isdisjoint(self.panel.columns))


if __name__ == "__main__":
    unittest.main()
