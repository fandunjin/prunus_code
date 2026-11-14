from pathlib import Path
import sys
import unittest
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from enrichment_utils import bh_group, partial_conjunction


class EnrichmentTests(unittest.TestCase):
    def test_partial_conjunction_includes_missing_species(self):
        self.assertAlmostEqual(partial_conjunction([0.01, 0.02, 0.8]), 0.04)
        self.assertAlmostEqual(partial_conjunction([0.01, 0.02]), 0.04)
        self.assertEqual(partial_conjunction([0.01]), 1)

    def test_bh_preserves_full_universe_and_row_order(self):
        values = pd.Series([0.04, 0.01, 1.0, 0.03], index=["a", "b", "c", "d"])
        result = bh_group(values)
        self.assertEqual(list(result.index), list(values.index))
        for got, want in zip(
            result, [0.05333333333333334, 0.04, 1.0, 0.05333333333333334]
        ):
            self.assertAlmostEqual(got, want)


if __name__ == "__main__":
    unittest.main()
