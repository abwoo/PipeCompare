"""Fixed retained event/weight decomposition: exact arithmetic and missingness tests."""
import sys
import unittest
from copy import deepcopy
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"benchmarks"))
from ica_stage_selection_decomp import stage_map, final_trial_reference, compare, difference


def fixture():
    z={
        "nEvents":4,"status":"complete","eventIds":[1,2,3,4],
        "conditionIndex":[1,1,2,2],
        "observedPzUv":[1.,3.,0.,2.]
    }
    conditions=[
        {"referenceIds":[1,2],"retainedIds":[1,2],"trialWeights":[1.,1.]},
        {"referenceIds":[3,4],"retainedIds":[3,4],"trialWeights":[1.,1.]},
    ]
    old={"science":{"preRemovalContributions":[z],"conditions":conditions}}
    new=deepcopy(old)
    new["science"]["preRemovalContributions"][0]["observedPzUv"]=[2.,4.,1.,3.]
    return old,new


class EventDecomposition(unittest.TestCase):
    def test_unchanged_event_weights_have_zero_selection_effect(self):
        old,new=fixture()
        x=compare(old,new,1.,1.,5.,5.2)
        self.assertAlmostEqual(x["observedTotalStagePzDeltaUv"],0.)
        self.assertAlmostEqual(x["fixedSeed0TrialsSignalDeltaUv"],0.)
        self.assertAlmostEqual(x["changedTrialSelectionOrWeightsDeltaUv"],0.)
        self.assertTrue(x["arithmeticallyConserved"])
        self.assertFalse(x["startDifferenceConfoundedByTrialSelectionOrWeights"])
        self.assertAlmostEqual(x["originalFinalGainCorrectedP300DeltaUv"],.2)

    def test_retained_trial_change_explains_stage_delta(self):
        old,new=fixture()
        new["science"]["conditions"][0]["retainedIds"]=[2]
        new["science"]["conditions"][1]["retainedIds"]=[3]
        x=compare(old,new,1.,3.,5.,5.)
        self.assertAlmostEqual(x["observedTotalStagePzDeltaUv"],2.)
        self.assertAlmostEqual(x["fixedSeed0TrialsSignalDeltaUv"],0.)
        self.assertAlmostEqual(x["changedTrialSelectionOrWeightsDeltaUv"],2.)
        self.assertTrue(x["startDifferenceConfoundedByTrialSelectionOrWeights"])

    def test_weight_change_only_is_distinct_from_signal_change(self):
        old,new=fixture()
        new["science"]["conditions"][0]["trialWeights"]=[1.,3.]
        x=compare(old,new,1.,1.5,3.,3.2)
        self.assertAlmostEqual(x["observedTotalStagePzDeltaUv"],.5)
        self.assertAlmostEqual(x["fixedSeed0TrialsSignalDeltaUv"],0.)
        self.assertAlmostEqual(x["changedTrialSelectionOrWeightsDeltaUv"],.5)

    def test_stage_signal_change_with_fixed_trials(self):
        old,new=fixture()
        new["science"]["preRemovalContributions"][0]["observedPzUv"][0]+=2.
        x=compare(old,new,1.,2.,5.,5.)
        self.assertAlmostEqual(x["fixedSeed0TrialsSignalDeltaUv"],1.)
        self.assertAlmostEqual(x["changedTrialSelectionOrWeightsDeltaUv"],0.)
        self.assertAlmostEqual(x["observedTotalStagePzDeltaUv"],1.)

    def test_original_retained_event_missing_fails(self):
        old,new=fixture()
        new["science"]["preRemovalContributions"][0]["eventIds"][0]=99
        with self.assertRaisesRegex(ValueError,"missing"):
            compare(old,new,1.,1.,5.,5.)

    def test_original_event_wrong_condition_fails(self):
        old,new=fixture()
        new["science"]["preRemovalContributions"][0]["conditionIndex"][0]=2
        with self.assertRaisesRegex(ValueError,"wrong condition"):
            compare(old,new,1.,1.,5.,5.)

    def test_original_event_duplicate_fails(self):
        old,new=fixture()
        new["science"]["preRemovalContributions"][0]["eventIds"][0]=2
        with self.assertRaisesRegex(ValueError,"misaligned"):
            stage_map(new)

    def test_original_sensor_earlier_verified_stage_mismatch_fails(self):
        old,new=fixture()
        with self.assertRaisesRegex(ValueError,"differs"):
            compare(old,new,1.,2.,5.,5.)

    def test_missing_original_weight_fails(self):
        old,new=fixture()
        new["science"]["conditions"][0]["trialWeights"][0]=0.
        with self.assertRaisesRegex(ValueError,"invalid"):
            final_trial_reference(new)

    def test_no_zero_fill_for_unobserved_starts(self):
        old,new=fixture()
        new["science"]["preRemovalContributions"]=[]
        with self.assertRaisesRegex(ValueError,"Missing"):
            compare(old,new,1.,1.,5.,5.)

    def test_bad_condition_count_raises(self):
        with self.assertRaisesRegex(ValueError,"Two-condition"):
            difference([1.])

if __name__=="__main__":
    unittest.main()
