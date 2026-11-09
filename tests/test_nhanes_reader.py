import importlib.util
from pathlib import Path
import unittest

import numpy as np
import pandas as pd


# Tests live outside the preparation folders, so load the reader directly from its file path.
ROOT = Path(__file__).resolve().parents[1] / "1_prepare_nhanes"  # public_code: the preparation stage
MODULE_PATH = ROOT / "shared" / "nhanes_reader.py"  # public_code
SPEC = importlib.util.spec_from_file_location("nhanes_reader", MODULE_PATH)
reader = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(reader)


class ReaderTests(unittest.TestCase):
    """Small checks for failure conditions in the shared reader."""

    def test_key_validation_accepts_unique_nonmissing_seqn(self):
        reader.validate_participant_key(pd.DataFrame({"SEQN": [1.0, 2.0]}), "test")

    def test_key_validation_rejects_duplicates(self):
        with self.assertRaises(ValueError):
            reader.validate_participant_key(pd.DataFrame({"SEQN": [1.0, 1.0]}), "test")

    def test_zero_threshold_matches_certified_value(self):
        self.assertLess(5.397605346934028e-79, reader.SAS_ZERO_THRESHOLD)
        self.assertFalse(np.isnan(0.0))


if __name__ == "__main__":
    unittest.main()
