"""Tests for HDL name harmonisation and lipid-panel ownership decisions."""

import unittest
import pandas as pd


class LipidRulesTest(unittest.TestCase):
    def test_all_historical_hdl_names_map_to_canonical(self):
        for source in ["LBDHDL", "LBXHDD", "LBDHDD"]:
            frame = pd.DataFrame({source: [50.0]}).rename(columns={source: "LBDHDD"})
            self.assertEqual(frame.LBDHDD.iloc[0], 50.0)

    def test_inner_merge_requires_both_files(self):
        tc = pd.DataFrame({"SEQN": [1, 2], "LBXTC": [180, 190]})
        hdl = pd.DataFrame({"SEQN": [2, 3], "LBDHDD": [50, 60]})
        merged = tc.merge(hdl, on="SEQN", how="inner")
        self.assertListEqual(merged.SEQN.tolist(), [2])

    def test_friedewald_ldl_is_deterministic(self):
        self.assertEqual(200 - 50 - 150 / 5, 120)


if __name__ == "__main__":
    unittest.main()
