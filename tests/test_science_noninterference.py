"""Scientific recorder must not silently change saved ranking and P300 results."""
import sys
import unittest
from copy import deepcopy
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"benchmarks"))
from verify_science_noninterference import one_pair


def result():
    return {"subject":11,"completed":True,"feasible":20,"candidates":28,"notRun":0,
            "errors":[],"source":{"source_set_sha256":"a","source_fdt_sha256":"b",
                                  "source_events_sha256":"c"},
            "erp":{"differenceAmplitude":1.5,"gain":0.8},
            "icaSeedCheck":{"ids":[12,25],"recommended":[12,25,25,25,25,25],"same":False,
                            "rankStatus":[["feasible"]*6,["feasible"]*6],
                            "strata":[["average"]*6,["average"]*6],
                            "objective":[[1.,1.,1.,1.,1.,1.],[1.5]*6],
                            "measurements":[{"seed":i,"candidateId":12 if i==0 else 25,
                                             "gain":0.8,"hasMeasurement":True,
                                             "correctedDifferenceEstimate":1.5} for i in range(6)]}}


class NonIntervention(unittest.TestCase):
    def test_recorder_addition_cannot_change_scientific_results(self):
        old=result();new=deepcopy(old)
        new["icaSeedCheck"]["measurements"][0]["science"]={"available":True}
        self.assertTrue(one_pair(old,new)["noninterventionPassed"])
    def test_changed_recommendation_detected(self):
        old=result();new=deepcopy(old)
        new["icaSeedCheck"]["recommended"][3]=12
        self.assertIn("ica-recommended-changed",one_pair(old,new)["reasons"])
    def test_changed_p300_detected(self):
        old=result();new=deepcopy(old)
        new["icaSeedCheck"]["measurements"][3]["correctedDifferenceEstimate"]=1.7
        self.assertIn("seed-3-correctedDifferenceEstimate-changed",
                      one_pair(old,new)["reasons"])
    def test_changed_source_identity_detected(self):
        old=result();new=deepcopy(old)
        new["source"]["source_set_sha256"]="different"
        self.assertIn("source-source_set_sha256-changed",one_pair(old,new)["reasons"])


if __name__=="__main__":
    unittest.main()
