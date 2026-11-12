"""Condition-labelled independent eye-proxy support checks, no gate mutation."""
import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))
from validate_veog_support import canonical, crosscheck, eye_max


def sample():
    # fs=10, sample windows [-2,+8], 1-based event positions 11 and 30.
    veog = np.zeros(60)
    veog[10] = 80.0
    veog[29] = 10.0
    mat = {"fs": 10.0, "eventLatency": np.array([11, 30]),
           "eventTypes": np.array(["11", "31"]),
           "veog": veog, "heog": np.zeros(60)}
    def support(name, ref_id, flagged):
        return {"condition": name, "referenceN": 1, "kept": 1,
                "referenceIds": [ref_id], "scores": [5.0], "weights": [1.0],
                "flagged": [flagged],
                "effectiveUnflaggedN": 0.0 if flagged else 1.0}
    result = {"subject": 4, "completed": True, "notRun": 0,
              "selectedSupport": [support("11", 1, True),
                                  support("31", 2, False)],
              "minimal": {"support": [support("target", 1, False),
                                      support("standard", 2, False)]}}
    return result, mat


class EyeProxyTests(unittest.TestCase):
    def test_stimulus_codes_are_not_implicitly_standard(self):
        self.assertEqual(canonical("11"), "target")
        self.assertEqual(canonical("31"), "standard")
        self.assertEqual(canonical("target"), "target")
        with self.assertRaises(ValueError):
            canonical("21")

    def test_alignment_and_condition_specific_spikes(self):
        result, mat = sample()
        rows, trials = crosscheck(result, mat)
        self.assertEqual(len(rows), 6)
        self.assertEqual(len(trials), 4)
        std, target = rows[1], rows[0]
        self.assertEqual(target["hmm_flagged"], 1)
        self.assertEqual(target["both_positive"], 1)
        self.assertEqual(target["veog_proxy_positive"], 1)
        self.assertEqual(std["veog_proxy_positive"], 0)
        self.assertAlmostEqual(trials[0]["veog_max_uv"], 80.0)
        self.assertAlmostEqual(trials[1]["veog_max_uv"], 10.0)
        self.assertFalse(rows[4]["available"])

    def test_rejects_swapped_source_conditions(self):
        result, mat = sample()
        result["selectedSupport"][0]["referenceIds"] = [2]
        with self.assertRaisesRegex(ValueError, "source event/condition mismatch"):
            crosscheck(result, mat)

    def test_rejects_misaligned_trial_arrays(self):
        result, mat = sample()
        result["selectedSupport"][1]["flagged"] = []
        with self.assertRaisesRegex(ValueError, "misaligned support"):
            crosscheck(result, mat)

    def test_never_imputes_incomplete_run(self):
        result, mat = sample()
        result["completed"] = False
        with self.assertRaisesRegex(ValueError, "incomplete"):
            crosscheck(result, mat)

    def test_dc_baseline_correction(self):
        signal = np.full(30, 100.0)
        signal[10] = 170.0
        self.assertAlmostEqual(eye_max(signal, 11, 10), 70.0)


if __name__ == "__main__":
    unittest.main()
