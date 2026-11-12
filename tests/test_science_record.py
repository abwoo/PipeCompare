"""Read-only science recorder independent validation regression."""
import sys
import unittest
from copy import deepcopy
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"benchmarks"))
from validate_science_record import inspect_one


def example():
    scores0=[1.,2.,3.,4.,5.]
    scores1=[0.,0.,1.,0.,1.]
    c=[]
    for name,scores in (("target",scores0),("standard",scores1)):
        c.append({"condition":name,"sourceIds":[11,12,13,14,15],
                  "keptIds":[11,12,13,14,15],
                  "weights":[1.,1.,1.,1.,1.],"scores":scores,
                  "waveform":[0.,1.,2.,3.,4.,5.,6.,7.,8.,9.],
                  "meanAmplitude":sum(scores)/5,"analyticSME":0.1})
    return {"seed":3,"candidateId":7,"gain":1.,
            "hasMeasurement":True,"correctedDifferenceEstimate":2.6,
            "science":{"available":True,"seed":3,"candidateId":7,
             "ica":[{"weights":[[1.,0.],[0.,1.]],"sphere":[[1.,0.],[0.,1.]],
                      "mixing":[[1.,0.],[0.,1.]],"componentCount":2,
                      "removedIndices":[1],"channelIndices":[1,2]}],
             "waveform":{"timesMs":list(range(10)),"difference":[0.]*10,
                         "condition":c}}}


class RecordIntegrity(unittest.TestCase):
    def test_complete_observational_record(self):
        z=inspect_one(example(),7,3)
        self.assertTrue(z["valid"],z["issues"])
        self.assertGreater(z["bootstrap"]["se"],0)
        self.assertEqual(z["bootstrap"]["draws"],1000)
    def test_component_identity_failure(self):
        x=example()
        x["science"]["ica"][0]["removedIndices"]=[3]
        z=inspect_one(x,7,3)
        self.assertIn("removed-component-id-invalid",z["issues"])
    def test_trial_id_alignment(self):
        x=example()
        x["science"]["waveform"]["condition"][0]["sourceIds"].pop()
        self.assertIn("trial-array-misaligned-0",inspect_one(x,7,3)["issues"])
    def test_gain_and_erp_difference_consistency(self):
        x=example()
        x["correctedDifferenceEstimate"]=99.
        self.assertIn("bootstrap-point-mismatch",inspect_one(x,7,3)["issues"])
    def test_no_record_not_mistaken_for_validity(self):
        x=example()
        x["science"]={"available":False,"reason":"absent"}
        self.assertFalse(inspect_one(x,7,3)["valid"])
    def test_deterministic_bootstrap(self):
        x=example()
        a=inspect_one(x,7,3)
        b=inspect_one(deepcopy(x),7,3)
        self.assertEqual(a["bootstrap"],b["bootstrap"])


if __name__=="__main__":
    unittest.main()
