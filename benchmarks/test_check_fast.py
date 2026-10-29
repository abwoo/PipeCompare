#!/usr/bin/env python3
"""Tests for Round 10 fast-loop Python guards (stdlib only).

Run: python3 -m unittest test_check_fast -v
These tests never run MATLAB, never touch the repo, and never
use network. They test classification logic only: a fast fragment
must be refused as held-out evidence, a full fragment must pass
the guard, and a calibration report must never carry a ratio
without coverage.
"""
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import check_shadow_selection as cs
import fast_replay as fr


def full_frag(seed=20031):
    return {"seed": seed, "arms": [],
            "labels": {"fullyHeldOut": True, "preprocessingRefitted": False,
                       "adaptiveLeavesIncluded": False}}


def micro_frag(seed=40101):
    f = full_frag(seed)
    f["labels"] = {"fullyHeldOut": False, "preprocessingRefitted": False,
                    "adaptiveLeavesIncluded": False}
    f["fastMode"] = {"mode": "slice-diagnostic", "partialSearch": True,
                     "slicedValidation": False, "designSeedDiagnostic": True}
    return f


class TestFastGuards(unittest.TestCase):
    def test_full_fragment_is_held_out(self):
        self.assertFalse(cs.is_fast_fragment(full_frag()))
        self.assertEqual(fr.evidence_class(full_frag()), "held-out")
        fr.require_paired_eligible([full_frag()])  # must not raise

    def test_micro_fragment_is_refused_as_held_out(self):
        self.assertTrue(cs.is_fast_fragment(micro_frag()))
        self.assertEqual(fr.evidence_class(micro_frag()), "partial")
        with self.assertRaises(ValueError):
            fr.require_paired_eligible([full_frag(), micro_frag()])

    def test_micro_fragment_valid_on_its_own_terms(self):
        self.assertEqual(cs.validate_micro_fragment(micro_frag()), [])

    def test_micro_fragment_claiming_held_out_fails_validation(self):
        f = micro_frag()
        f["labels"]["fullyHeldOut"] = True
        self.assertIn("must not claim fullyHeldOut=true",
                      " ".join(cs.validate_micro_fragment(f)))

    def test_micro_fragment_on_held_out_seed_fails_validation(self):
        f = micro_frag(seed=20031)
        self.assertTrue(cs.validate_micro_fragment(f))

    def test_sliced_validation_is_never_paired(self):
        f = micro_frag()
        f["fastMode"] = {"mode": "slice-diagnostic", "partialSearch": False,
                         "slicedValidation": True, "designSeedDiagnostic": True}
        self.assertEqual(fr.evidence_class(f), "sliced-validation")
        with self.assertRaises(ValueError):
            fr.require_paired_eligible([f])

    def test_calibration_report_always_carries_coverage(self):
        recs = [{"assessable": True, "actual_sq": 4.0, "reported_sq": 1.0},
                {"assessable": False, "actual_sq": 9.0, "reported_sq": 1.0}]
        rep = fr.calibration_report(recs)
        self.assertEqual(rep["coverage"], 0.5)
        self.assertFalse(rep["complete"])
        self.assertIsNotNone(rep["ratio"])
        full = fr.calibration_report(
            [{"assessable": True, "actual_sq": 4.0, "reported_sq": 1.0}])
        self.assertTrue(full["complete"])
        self.assertEqual(full["coverage"], 1.0)


if __name__ == "__main__":
    unittest.main()
