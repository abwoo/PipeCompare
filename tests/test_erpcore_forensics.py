"""Unit tests for the PR51 frozen-run read-only forensic parser."""
import importlib.util
from pathlib import Path
import unittest

MODULE = Path(__file__).resolve().parents[1] / "benchmarks" / "diagnose_erpcore_logs.py"
spec = importlib.util.spec_from_file_location("pr51_forensics", MODULE)
f = importlib.util.module_from_spec(spec)
spec.loader.exec_module(f)


def sample():
    return {
        "subject": 18, "completed": True,
        "algorithmCommit": "original", "executionCommit": "frozen",
        "workflowRun": "42809789847", "protocolSha256": "protocol",
        "candidates": 24, "feasible": 0,
        "noRecommendationReasons": [
            "icremove: Unable to perform assignment because the left and right sides have a different number of elements."
        ] * 18 + ["filters alone: component amplitude changed by 16%"] * 6,
        "erp": {"available": False},
        "minimal": {"available": False},
        "common": {"available": False},
        "errors": [],
    }


class FrozenLogTests(unittest.TestCase):
    def test_parse_timestamped_matlab_marker(self):
        import json
        r = sample()
        body = "2026-11-10T12:00:00Z " + f.MARKER + json.dumps(r) + "\n"
        self.assertEqual(f.result_from_log(body), r)
        with self.assertRaises(ValueError):
            f.result_from_log(body + body)
        with self.assertRaises(ValueError):
            f.result_from_log("no result marker")

    def test_runtime_error_is_not_relabelled_as_data_quality(self):
        r = sample()
        summary = f.summarize(r, "log-hash")
        self.assertEqual(summary["primary_diagnostic"], "icremove-execution-error-and-other-gates")
        self.assertEqual(summary["feasible_count"], 0)
        self.assertEqual(sorted(summary["reasons"].values()), [6, 18])

    def test_shape_gate_is_distinct_from_exception(self):
        r = sample()
        r["subject"] = 17
        r["noRecommendationReasons"] = ["filters alone: peak latency shifted by 31 ms"] * 24
        self.assertEqual(f.classify(r), "signal-shape-gate")

    def test_effective_support_differs_from_retained_counts(self):
        r = sample()
        r["subject"] = 8
        r["feasible"] = 19
        r["erp"] = {"available": True, "kept": [16, 71]}
        r["selectedSupport"] = [
            {"condition": "11", "referenceIds": list(range(22)),
             "flagged": [True] * 17 + [False] * 5,
             "kept": 16, "effectiveUnflaggedN": 5.0},
            {"condition": "31", "referenceIds": list(range(95)),
             "flagged": [True] * 75 + [False] * 20,
             "kept": 71, "effectiveUnflaggedN": 20.0},
        ]
        r["minimal"] = {"available": True, "support": [{
            "condition": "target", "referenceIds": list(range(22)),
            "flagged": [True] * 17 + [False] * 5,
            "kept": 4, "effectiveUnflaggedN": 4.0
        }]}
        r["common"] = {"available": True, "support": [{
            "condition": "target", "referenceIds": list(range(22)),
            "flagged": [True] * 17 + [False] * 5,
            "kept": 4, "effectiveUnflaggedN": 4.0
        }]}
        r["icaSeedCheck"] = {"recommended": [25, 25, 25], "ids": [25],
                             "objective": [[1.2, 1.2, 1.2]], "same": True}
        self.assertEqual(f.classify(r), "effective-unflagged-support-below-frozen-floor")
        self.assertEqual(f.support(r, "erp")[0]["effective_unflagged"], 5.0)
        self.assertEqual(f.support(r, "minimal")[0]["effective_unflagged"], 4.0)

    def test_seed_identity_and_objective_are_separate(self):
        r = sample()
        r["icaSeedCheck"] = {"recommended": [12, 25, 25],
                             "ids": [12, 25],
                             "objective": [[0.7, 0.8, 0.8], [0.9, 0.72, 0.73]],
                             "same": False}
        values = f.seed_summary(r)
        self.assertEqual(values["recommended_candidate_ids"], [12, 25, 25])
        self.assertEqual(values["recommended_objective"], [0.7, 0.72, 0.73])
        self.assertFalse(values["same"])


if __name__ == "__main__":
    unittest.main()
