"""Regression tests for MATLAB singleton JSON and independent summary checking."""
import copy
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_shadow_selection as check


class WireFormatTests(unittest.TestCase):
    def test_numeric_array_shapes(self):
        self.assertEqual(check.matlab_array([]), [])
        self.assertEqual(check.matlab_array(20033), [20033])
        self.assertEqual(check.matlab_array([20033, 20034]), [20033, 20034])

    def test_schema_normalization_preserves_provenance_and_values(self):
        raw = {"seeds": 20033, "perSeed": {"seed": 20033, "arms": {"name": "standard"}},
               "pairedDifferences": {"p": {"seeds": 20033, "gainCorrected": -0.5, "raw": 1.25}},
               "recordBootstrapCI95": {"p": {"replicates": [0.1, 0.2]}},
               "provenance": {"commit": "original-commit", "workflowRunId": "original-run"}}
        original = copy.deepcopy(raw)
        normalized = check.normalize_run_arrays(raw)
        self.assertEqual(raw, original)
        self.assertEqual(normalized["provenance"], original["provenance"])
        self.assertEqual(normalized["pairedDifferences"]["p"]["seeds"], [20033])
        self.assertEqual(normalized["pairedDifferences"]["p"]["gainCorrected"], [-0.5])
        self.assertEqual(normalized["pairedDifferences"]["p"]["raw"], [1.25])
        self.assertEqual(normalized["perSeed"][0]["arms"], [{"name": "standard"}])

    def test_missing_pick_encodings_preserve_numeric_ids(self):
        self.assertIsNone(check.normalized_pick(None))
        self.assertIsNone(check.normalized_pick([]))
        self.assertEqual(check.normalized_pick(1), 1)
        self.assertEqual(check.normalized_pick(0), 0)

    def test_empty_arrays_remain_empty(self):
        raw = {"seeds": [], "perSeed": [], "pairedDifferences": {"p": {"seeds": [], "gainCorrected": [], "raw": []}}}
        self.assertEqual(check.normalize_run_arrays(raw), raw)

    def test_matlab_fixtures_and_tamper_detection(self):
        root = Path(os.environ["SHADOW_FIXTURES"])
        checker = Path(check.__file__).resolve()
        for case, count in (("zero", 1), ("single", 1), ("multiple", 3)):
            with self.subTest(case=case):
                folder = root / case
                run = next(p for p in folder.glob("shadow_selection_*.json"))
                csv_file = run.with_suffix(".csv")
                summary = run.with_name(run.stem + "_summary.csv")
                command = [sys.executable, str(checker), str(folder / "fragments"),
                           str(run), str(csv_file), str(summary), "fixture1234567890", str(count)]
                result = subprocess.run(command, text=True, capture_output=True)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                with tempfile.TemporaryDirectory() as dest:
                    bad = Path(dest) / run.name
                    data = json.loads(run.read_text())
                    data["pairedMeanDifference"]["guardedMinusStandard"] = 12345
                    bad.write_text(json.dumps(data))
                    result = subprocess.run(command[:3] + [str(bad)] + command[4:], text=True, capture_output=True)
                    self.assertNotEqual(result.returncode, 0)
                    self.assertIn("mean difference mismatch", result.stdout)


if __name__ == "__main__":
    unittest.main()
