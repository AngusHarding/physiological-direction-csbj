"""Regression tests for the FI reference-membership gate."""

import importlib.util
from pathlib import Path
import unittest

import pandas as pd

ROOT = Path(__file__).resolve().parents[1] / "1_prepare_nhanes"  # public_code: the preparation stage
MODULE_PATH = ROOT / "8_adverse_distance/01_lock_reference.py"  # public_code
MODULE_SPEC = importlib.util.spec_from_file_location("lock_fi_reference", MODULE_PATH)
module = importlib.util.module_from_spec(MODULE_SPEC)
MODULE_SPEC.loader.exec_module(module)


def inputs():
    transformed = pd.read_csv(module.TRANSFORMED)
    reference = pd.read_csv(module.REFERENCE)
    markers = pd.read_csv(module.SCORING_SPEC).marker.tolist()
    return transformed, reference, markers


class TestFIReference(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.transformed, cls.reference, cls.markers = inputs()

    def test_locked_reference_passes(self):
        membership, summary = module.build_membership(
            self.transformed, self.reference, self.markers
        )
        self.assertEqual(int(membership.is_fi_reference.sum()), 1238)
        self.assertEqual(summary.reference_n.tolist(), [547, 691])

    def test_participant_order_is_preserved(self):
        membership, _ = module.build_membership(
            self.transformed, self.reference, self.markers
        )
        self.assertTrue(membership.SEQN.equals(self.transformed.SEQN))

    def test_changed_reference_membership_fails(self):
        changed = self.reference.copy()
        changed.loc[changed.index[0], "is_reference"] = True
        with self.assertRaisesRegex(AssertionError, "counts changed|age 30 or older"):
            module.build_membership(self.transformed, changed, self.markers)


if __name__ == "__main__":
    unittest.main()
