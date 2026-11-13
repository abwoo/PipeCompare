"""Mechanism-level ERP stage contrasts must conserve original event and IC signal."""
import sys
import unittest
from copy import deepcopy
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"benchmarks"))
from ica_postremoval_mechanisms import assess


def samples():
    conditions=[
        {"referenceIds":[1,2],"retainedIds":[1,2],"trialWeights":[1.,1.]},
        {"referenceIds":[3,4],"retainedIds":[3,4],"trialWeights":[1.,1.]}]
    pre0={"status":"complete","nEvents":4,"srate":100.,
          "eventIds":[1,2,3,4],"conditionIndex":[1,1,2,2],
          "observedPzUv":[2.,4.,1.,3.],
          "componentUv":[[0.,0.,0.,0.],[2.,4.,1.,3.]]}
    post0={"status":"complete","nEvents":4,"srate":100.,
           "eventIds":[1,2,3,4],"conditionIndex":[1,1,2,2],
           "observedPzUv":[2.,4.,1.,3.]}
    pre1=deepcopy(pre0)
    pre1["componentUv"]=[[0.,2.,0.,0.],[2.,2.,1.,3.]]
    post1=deepcopy(post0)
    post1["observedPzUv"]=[2.,2.,1.,3.]
    old={"science":{"status":"complete","conditions":deepcopy(conditions),
       "preRemovalContributions":[pre0],"postRemovalSignals":[post0],
       "decompositions":[{"nComponents":2,"removed":[1]}]}}
    new={"science":{"status":"complete","conditions":deepcopy(conditions),
       "preRemovalContributions":[pre1],"postRemovalSignals":[post1],
       "decompositions":[{"nComponents":2,"removed":[1]}]}}
    oldstage={"preRemovalFinalSelectedTrialPzDifferenceUv":1.,
              "postRemovalFinalSelectedTrialPzDifferenceUv":1.}
    newstage={"preRemovalFinalSelectedTrialPzDifferenceUv":1.,
              "postRemovalFinalSelectedTrialPzDifferenceUv":0.}
    return old,new,oldstage,newstage


class ActualPostMechanism(unittest.TestCase):
    def test_fixed_event_removal_effect_and_conservation(self):
        x=assess(*samples(),1.,1.4,.2)
        self.assertAlmostEqual(x["preRemovalFixedEventsSignalDriftUv"],0.)
        self.assertAlmostEqual(x["postRemovalFixedEventsSignalDriftUv"],-1.)
        self.assertAlmostEqual(x["actualICDeletionEffectDriftUv"],-1.)
        self.assertAlmostEqual(x["predictedICDeletionEffectDriftUv"],-1.)
        self.assertAlmostEqual(x["postRemovalTrialSelectionOrWeightTermUv"],0.)
        self.assertAlmostEqual(x["originalFinalGainCorrectedP300DriftUv"],.4)
        self.assertFalse(x["fullERPCausalEffectOrSuperiorRecipeEstablished"])

    def test_changing_final_weights_is_separate_from_physical_removal(self):
        a,b,c,d=samples()
        b["science"]["conditions"][1]["trialWeights"]=[1.,3.]
        d["postRemovalFinalSelectedTrialPzDifferenceUv"]=-.5
        x=assess(a,b,c,d,1.,1.,0.)
        self.assertAlmostEqual(x["actualICDeletionEffectDriftUv"],-1.)
        self.assertNotEqual(x["postRemovalTrialSelectionOrWeightTermUv"],0.)
        self.assertTrue(x["startDifferenceConfoundedByTrialSelectionOrWeights"])

    def test_stage_predictor_fails_when_source_does_not_match_actual_removed_ica(self):
        a,b,c,d=samples()
        b["science"]["preRemovalContributions"][0]["componentUv"][0][1]=1.
        with self.assertRaisesRegex(ValueError,"not predicted"):
            assess(a,b,c,d,1.,1.,0.)

    def test_undeclared_event_fails_closed(self):
        a,b,c,d=samples()
        b["science"]["postRemovalSignals"][0]["eventIds"][0]=99
        with self.assertRaisesRegex(ValueError,"physical source events"):
            assess(a,b,c,d,1.,1.,0.)

    def test_missing_source_measured_stage_fails_closed(self):
        a,b,c,d=samples()
        b["science"]["postRemovalSignals"]=[]
        with self.assertRaisesRegex(ValueError,"Absent/duplicate"):
            assess(a,b,c,d,1.,1.,0.)

    def test_wrong_old_validated_source_stage_fails_closed(self):
        a,b,c,d=samples()
        c["preRemovalFinalSelectedTrialPzDifferenceUv"]=99.
        with self.assertRaisesRegex(ValueError,"not reproduced"):
            assess(a,b,c,d,1.,1.,0.)

    def test_original_final_trial_missing_from_both_stages_fails(self):
        a,b,c,d=samples()
        b["science"]["conditions"][0]["retainedIds"]=[1,99]
        with self.assertRaisesRegex(ValueError,"Invalid original retained"):
            assess(a,b,c,d,1.,1.,0.)

if __name__=="__main__":
    unittest.main()
