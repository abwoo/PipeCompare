"""Regression tests for stage-specific ICA P300 contribution conservation."""
import sys
import unittest
from copy import deepcopy
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"benchmarks"))
from ica_stage_impact import (matched_stage_delta, fixed_trial_comparison,
                              stage_contributions)


def fixture():
    d={"nComponents":2,"nChannels":4,
       "channelLabels":["Fz","Cz","Pz","Oz"],
       "mixing":[[1.,0.],[0.,1.],[-1.,0.],[0.,-1.]],
       "removed":[1],"classNames":["Brain","Eye"],
       "classProbabilities":[[.8,.2],[.1,.9]]}
    s0={"stageComponentPzDifferenceUv":[1.,2.],
        "stageSensorPzDifferenceUv":4.,"stageResidualPzDifferenceUv":1.,
        "joinedFinalTrialIdentities":True,"equalsFinalP300Estimand":False}
    s1={"stageComponentPzDifferenceUv":[1.3,2.4],
        "stageSensorPzDifferenceUv":5.2,"stageResidualPzDifferenceUv":1.5,
        "joinedFinalTrialIdentities":True,"equalsFinalP300Estimand":False}
    a={"science":{"status":"complete","decompositions":[d]}}
    b=deepcopy(a)
    return s0,s1,a,b


class StageProjection(unittest.TestCase):
    def test_reconstruction_changes_and_residual_are_conserved(self):
        s0,s1,a,b=fixture()
        out=matched_stage_delta(s0,s1,a,b)
        self.assertEqual(out["weakMatchCount"],0)
        self.assertAlmostEqual(out["observedPreRemovalSensorPzDeltaUv"],1.2)
        self.assertAlmostEqual(out["signedStrongMatchedContributionDeltaUv"],.7)
        self.assertAlmostEqual(out["sensorReconstructionResidualDeltaUv"],.5)
        self.assertTrue(out["stageSensorResidualConserved"])
        self.assertFalse(out["componentSpecificFinalERPCausalEffectKnown"])

    def test_component_sign_permutation_preserves_pre_projected_amplitudes(self):
        s0,s1,a,b=fixture()
        dec=b["science"]["decompositions"][0]
        dec["mixing"]=[[-r[1],r[0]] for r in dec["mixing"]]
        dec["removed"]=[2]
        dec["classProbabilities"]=list(reversed(dec["classProbabilities"]))
        s1["stageComponentPzDifferenceUv"]=[2.4,1.3]
        out=matched_stage_delta(s0,s1,a,b)
        self.assertAlmostEqual(out["signedStrongMatchedContributionDeltaUv"],.7)
        self.assertEqual(out["weakMatchCount"],0)

    def test_stage_sensor_unequal_source_plus_residual_aborts(self):
        s0,s1,a,b=fixture()
        s1["stageSensorPzDifferenceUv"]=3.
        with self.assertRaisesRegex(ValueError,"mismatch"):
            matched_stage_delta(s0,s1,a,b)

    def test_missing_component_does_not_get_imputed_to_zero(self):
        s0,s1,a,b=fixture()
        s1["stageComponentPzDifferenceUv"]=[]
        with self.assertRaisesRegex(ValueError,"Missing"):
            matched_stage_delta(s0,s1,a,b)

    def test_unjustified_final_erp_claim_rejected(self):
        s0,s1,a,b=fixture()
        s0["equalsFinalP300Estimand"]=True
        with self.assertRaisesRegex(ValueError,"incorrectly"):
            matched_stage_delta(s0,s1,a,b)

    def test_weak_retained_component_is_separately_reported(self):
        s0,s1,a,b=fixture()
        a["science"]["decompositions"][0]["removed"]=[]
        b["science"]["decompositions"][0]["removed"]=[]
        # Change one topography enough to weaken the match but keep all
        # components in both decompositions; Pz contributions still conserve.
        b["science"]["decompositions"][0]["mixing"]=[[1.,0.],[0.,.7],
            [-1.,.7],[0.,-1.4]]
        got=matched_stage_delta(s0,s1,a,b)
        self.assertGreaterEqual(got["weakMatchCount"],1)
        self.assertGreaterEqual(got["weakRetainedBothCount"],1)
        self.assertGreater(got["absoluteRetainedWeakDeltaBudgetUv"],0.)

    def test_identical_retained_trials_and_weights_are_comparable(self):
        a={"conditions":[{"referenceIds":[1,2],"retainedIds":[1,2],"trialWeights":[.5,1.]},
                         {"referenceIds":[3,4],"retainedIds":[3,4],"trialWeights":[1.,1.]}]}
        b=deepcopy(a)
        z=fixed_trial_comparison(a,b)
        self.assertTrue(z["sameOriginalTrialWeights"])
        self.assertFalse(z["startDifferenceConfoundedByTrialSelectionOrWeights"])
        b["conditions"][0]["trialWeights"][0]=.25
        self.assertTrue(fixed_trial_comparison(a,b)["startDifferenceConfoundedByTrialSelectionOrWeights"])
        b=deepcopy(a)
        b["conditions"][0]["retainedIds"]=[1]
        self.assertFalse(fixed_trial_comparison(a,b)["sameRetainedTrialIdentity"])

    def test_failed_rank_component_cannot_prove_neural_truth(self):
        s0,s1,a,b=fixture()
        s0["stageComponentPzDifferenceUv"]=[float("nan"),2]
        with self.assertRaisesRegex(ValueError,"Missing"):
            stage_contributions(s0)


if __name__=="__main__":
    unittest.main()
