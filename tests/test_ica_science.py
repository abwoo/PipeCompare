"""Scientific recorder comparison and alignment regression tests."""
import sys
import unittest
from copy import deepcopy
from pathlib import Path
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"benchmarks"))
from validate_ica_science import validate_science, match_components, validate_identity


def synth():
    times=np.arange(-200,801,10).astype(float)
    wave0=2+np.exp(-((times-400)/70)**2)
    wave1=np.zeros_like(times)
    p3=(times>=300)&(times<=600)
    score0=float(wave0[p3].mean())
    base={"nComponents":2,"nChannels":4,
          "channelLabels":["Fz","Cz","Pz","Oz"],
          "mixing":[[1.,2.],[2.,-1.],[3.,0.],[4.,1.]],
          "removed":[1],"classNames":["Brain","Eye"],
          "classProbabilities":[[.2,.8],[.8,.2]]}
    cond=[
      {"name":"target","referenceIds":[12,13,14],"retainedIds":[12,13,14],
       "referenceFlagged":[False,False,False],
       "trialScores":[score0-.2,score0,score0+.2],"trialWeights":[1,1,1],
       "usedCount":3,"effectiveN":3,"meanWaveform":wave0.tolist(),
       "windowMean":score0,"originalEstimate":score0},
      {"name":"standard","referenceIds":[21,22,23],"retainedIds":[21,22,23],
       "referenceFlagged":[False,False,False],
       "trialScores":[-.1,0,.1],"trialWeights":[1,1,1],
       "usedCount":3,"effectiveN":3,"meanWaveform":wave1.tolist(),
       "windowMean":0,"originalEstimate":0}]
    measure={"candidateId":4,"seed":0,"gain":1.,
             "sourceUnflagged":[3,3],"effectiveUnflagged":[3,3],
             "correctedDifferenceEstimate":score0,"science":{
                "status":"complete","decompositions":[base],
                "timeMs":times.tolist(),"channel":"Pz","conditions":cond}}
    return measure


class ScienceRecorderTest(unittest.TestCase):
    def test_finite_science_snapshot_and_bootstrap(self):
        x=validate_science(synth(),11,4,0)
        self.assertTrue(x["complete"],x["errors"])
        self.assertGreater(x["bootstrapSE"],0.)
    def test_modified_p3_waveform_rejected(self):
        x=synth();x["science"]["conditions"][0]["meanWaveform"][60]+=999.
        got=validate_science(x,11,4,0)
        self.assertFalse(got["complete"])
    def test_missing_trial_identity_rejected(self):
        x=synth();x["science"]["conditions"][0]["retainedIds"]=[12,14]
        got=validate_science(x,11,4,0)
        self.assertFalse(got["complete"])
    def test_scientific_source_support_disagreement_is_rejected(self):
        x=synth()
        x["sourceUnflagged"]=[4,3]
        result=validate_science(x,11,4,0)
        self.assertFalse(result["complete"])
        self.assertIn("physical-hmm-support-disagrees-with-source:0",result["errors"])

    def test_weighted_unflagged_support_disagreement_is_rejected(self):
        x=synth()
        x["effectiveUnflagged"]=[2,3]
        result=validate_science(x,11,4,0)
        self.assertFalse(result["complete"])
        self.assertIn("weighted-hmm-support-disagrees-with-source:0",result["errors"])

    def test_iclabel_shape_incoherent(self):
        x=synth();x["science"]["decompositions"][0]["classProbabilities"]=[[1.,0.]]
        got=validate_science(x,11,4,0)
        self.assertFalse(got["complete"])
        self.assertIn("iclabel-probability-misaligned:0",got["errors"])
    def test_mixing_permutation_and_sign_are_matched(self):
        x=validate_science(synth(),11,4,0)
        y=deepcopy(x)
        topo=y["topology"][0]
        topo["mixing"]=np.column_stack([-topo["mixing"][:,1],topo["mixing"][:,0]])
        topo["removed"]=[2]
        match=match_components(x,y)
        self.assertTrue(match["available"])
        self.assertAlmostEqual(match["medianAbsoluteTopographicCorrelation"],1.,places=10)
        self.assertAlmostEqual(match["removedComponentMatchedJaccard"],1.,places=10)
    def test_removed_component_mismatch_remains_visible(self):
        x=validate_science(synth(),11,4,0);y=deepcopy(x)
        y["topology"][0]["removed"]=[2]
        match=match_components(x,y)
        self.assertEqual(match["removedComponentMatchedJaccard"],0.)
    def test_cross_montage_not_compared(self):
        x=validate_science(synth(),11,4,0);y=deepcopy(x)
        y["topology"][0]["labels"][0]="other"
        self.assertFalse(match_components(x,y)["available"])
    def test_original_record_identity_not_assumed(self):
        prior={"subject":11,"algorithmCommit":"a","protocolSha256":"p",
               "measureSha256":"m","completed":True,"notRun":0,
               "icaSeedCheck":{"ids":[4],"recommended":[4]*6,
                    "objective":[[1.]*6],"rankStatus":[["feasible"]*6],
                    "strata":[["s"]*6],"candidateMeasurements":[[{
                      "candidateId":4,"seed":s,"hasMeasurement":True,
                      "gain":1.,"correctedDifferenceEstimate":2.,
                      "effectiveUnflagged":[10,10],
                      "sourceUnflagged":[10,10]} for s in range(6)]]},
               "erp":{"differenceAmplitude":2.}}
        new=deepcopy(prior)
        self.assertEqual(validate_identity(new,prior),[])
        new["icaSeedCheck"]["recommended"][5]=9
        self.assertIn("recommended-changed",validate_identity(new,prior))


if __name__=="__main__":
    unittest.main()
