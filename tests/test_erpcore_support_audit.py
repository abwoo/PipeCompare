"""Scientific guardrails for the exploratory PR51 ERP CORE support auditor."""
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))
from audit_erpcore_support import audit, effective


def case():
    return {
        "completed": True, "notRun": 0, "subject": 8,
        "selectedSupport": [
            {"condition": "target", "referenceN": 3, "kept": 2,
             "referenceIds": [100, 101, 102], "scores": [1., 2., None],
             "weights": [1., 1., 0.], "flagged": [False, True, False],
             "positiveN": 2, "flaggedN": 1, "unflaggedN": 1,
             "effectiveUnflaggedN": 1., "eyeBiasAssessed": True},
            {"condition": "standard", "referenceN": 3, "kept": 3,
             "referenceIds": [200, 201, 202], "scores": [2., 2., 2.],
             "weights": [1., 1., 1.], "flagged": [False, False, False],
             "positiveN": 3, "flaggedN": 0, "unflaggedN": 3,
             "effectiveUnflaggedN": 3., "eyeBiasAssessed": True}
        ],
        "minimal": {"support": None},
        "common": {"support": None}
    }


class SupportAuditTests(unittest.TestCase):
    def test_effective_n_respects_weights(self):
        self.assertAlmostEqual(effective([1., 2.]), 1.8)
        self.assertEqual(effective([0., 0.]), 0.)

    def test_flags_are_not_retained_trials(self):
        summaries, per_trial = audit(case())
        self.assertEqual(len(summaries), 6)
        self.assertEqual(len(per_trial), 6)
        self.assertEqual(summaries[0]["kept"], 2)
        self.assertEqual(summaries[0]["effectiveUnflaggedN"], 1.)
        self.assertFalse(summaries[0]["meets_original_floor"])
        self.assertTrue(per_trial[1]["blink_proxy"])

    def test_absent_comparator_is_not_imputed(self):
        summaries, _ = audit(case())
        self.assertEqual(sum(x["available"] for x in summaries), 2)
        self.assertTrue(all(not x["available"] for x in summaries[2:]))

    def test_does_not_allow_incomplete_result(self):
        x = case()
        x["completed"] = False
        with self.assertRaises(ValueError):
            audit(x)

    def test_rejects_misaligned_trial_vectors(self):
        x = case()
        x["selectedSupport"][0]["weights"] = [1., 1.]
        with self.assertRaises(ValueError):
            audit(x)

    def test_detects_bad_effective_support_claim(self):
        x = case()
        x["selectedSupport"][0]["effectiveUnflaggedN"] = 2.5
        with self.assertRaises(ValueError):
            audit(x)


if __name__ == "__main__":
    unittest.main()
