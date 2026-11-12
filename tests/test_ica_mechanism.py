"""Proof that no-op ICA and unavailable metadata are distinguished."""
import math
import sys
import unittest
from types import SimpleNamespace
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"benchmarks"))

from diagnose_ica_mechanism import removal_snapshot, condition_summary


def cand(steps, n):
    return SimpleNamespace(steps=[SimpleNamespace(type=t,icsRemoved=count)
        for t,count in steps], icsRemoved=n)


class MechanismAudit(unittest.TestCase):
    def test_fitted_ica_without_removing_anything_is_not_ica_stability(self):
        c=cand([("ica",float("nan")),("icremove",0)],0)
        z=removal_snapshot(c)
        self.assertEqual(z["seed0RemovedComponents"],0)
        self.assertEqual(z["seed0Mechanism"],"ica-fitted-no-components-removed-at-seed0")

    def test_fitted_ica_with_removed_components(self):
        c=cand([("ica",float("nan")),("icremove",3)],3)
        z=removal_snapshot(c)
        self.assertEqual(z["seed0RemovedComponents"],3)
        self.assertEqual(z["seed0Mechanism"],"ica-fitted-components-removed-at-seed0")

    def test_missing_component_counts_are_unavailable_not_zero(self):
        c=cand([("ica",float("nan")),("icremove",float("nan"))],float("nan"))
        z=removal_snapshot(c)
        self.assertIsNone(z["seed0RemovedComponents"])
        self.assertEqual(z["seed0Mechanism"],"ica-removal-unreported")

    def test_no_ica_step_is_not_ica_stability(self):
        c=cand([("highpass",float("nan"))],0)
        z=removal_snapshot(c)
        self.assertFalse(z["seed0IcaPresent"])
        self.assertEqual(z["seed0Mechanism"],"no-ica-step")

    def test_effective_trial_count_uses_frozen_floor_without_rounding(self):
        x={"hasMeasurement":True,"correctedDifferenceEstimate":2.0,
           "effectiveUnflagged":[4.999,50]}
        self.assertFalse(condition_summary(x)["bothConditionsAtLeast5"])
        x["effectiveUnflagged"]=[5,50]
        self.assertTrue(condition_summary(x)["bothConditionsAtLeast5"])

    def test_missing_trial_counts_fail_closed(self):
        for eff in ([],None,[float("nan"),50],[5]):
            z=condition_summary({"hasMeasurement":True,
              "correctedDifferenceEstimate":2.0,"effectiveUnflagged":eff})
            self.assertFalse(z["bothConditionsAtLeast5"])

if __name__=="__main__":
    unittest.main()
