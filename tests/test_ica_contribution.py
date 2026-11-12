"""Read-only pre-removal ICA contribution validator: event and residual parity."""
import sys
import unittest
from copy import deepcopy
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"benchmarks"))
from validate_ica_contribution import contribution


def fixture():
    raw=[[1.,2.,3.,4.],[.5,.5,.5,.5]]
    stage={
       "status":"complete","eventIds":[1,2,3,4],
       "conditionIndex":[1,1,2,2],"nComponents":2,
       "nEvents":4,"nICAChannels":2,"rank":2,
       "componentUv":raw,"observedPzUv":[1.5,2.5,3.5,4.5],
       "residualPzUv":[0.,0.,0.,0.],"maxAbsResidualUv":0.,
    }
    final={"conditions":[
        {"referenceIds":[1,2],"retainedIds":[1,2],"trialWeights":[1,1]},
        {"referenceIds":[3,4],"retainedIds":[3,4],"trialWeights":[1,1]},
    ]}
    return stage,final


class AttributionIntegrity(unittest.TestCase):
    def test_full_rank_original_sensor_reconstruction(self):
        stage,final=fixture()
        report=contribution(stage,final)
        self.assertTrue(report["joinedFinalTrialIdentities"])
        self.assertFalse(report["equalsFinalP300Estimand"])
        self.assertEqual(report["stageComponentPzDifferenceUv"],[-2.,0.])
        self.assertEqual(report["stageSensorPzDifferenceUv"],-2.)
    def test_rank_deficiency_residual_is_kept_not_hidden(self):
        stage,final=fixture()
        stage["rank"]=1
        stage["residualPzUv"]=[.25]*4
        stage["observedPzUv"]=[x+.25 for x in stage["observedPzUv"]]
        stage["maxAbsResidualUv"]=.25
        result=contribution(stage,final)
        self.assertEqual(result["maxAbsStageReconstructionResidualUv"],.25)
    def test_full_rank_missing_sensor_signal_fails_closed(self):
        stage,final=fixture()
        stage["residualPzUv"]=[.25]*4
        stage["observedPzUv"]=[x+.25 for x in stage["observedPzUv"]]
        stage["maxAbsResidualUv"]=.25
        with self.assertRaisesRegex(ValueError,"Full-rank"):
            contribution(stage,final)
    def test_original_candidate_trial_missing_from_stage_fails(self):
        stage,final=fixture()
        final["conditions"][0]["retainedIds"]=[1,99]
        with self.assertRaisesRegex(ValueError,"absent"):
            contribution(stage,final)
    def test_missing_or_misaligned_component_matrix_fails(self):
        stage,final=fixture()
        stage["componentUv"]=[[1.,2.],[3.,4.]]
        with self.assertRaisesRegex(ValueError,"Misaligned"):
            contribution(stage,final)
    def test_sensor_identity_mismatch_fails(self):
        stage,final=fixture()
        stage["observedPzUv"][0]+=1
        with self.assertRaisesRegex(ValueError,"reconstruction"):
            contribution(stage,final)
    def test_not_assume_pre_ica_measure_is_final_erp(self):
        stage,final=fixture()
        # Artificially choose only one trial in the first condition:
        # stage metric follows the ORIGINAL fixed selection, not truth.
        final["conditions"][0]["retainedIds"]=[2]
        report=contribution(stage,final)
        self.assertEqual(report["stageComponentPzDifferenceUv"],[-1.5,0.])
        self.assertFalse(report["equalsFinalP300Estimand"])

if __name__=="__main__":
    unittest.main()
