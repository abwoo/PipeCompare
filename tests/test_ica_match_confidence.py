"""P1/P2 matching-confidence forensic tests: deterministic, no EEG fit."""
import sys
import unittest
from copy import deepcopy
from pathlib import Path

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"benchmarks"))
from ica_match_confidence import assignment, _waveform, _cycle

def test_maps():
    a={"nComponents":3,"nChannels":5,
       "channelLabels":["Fp1","Fz","Cz","Pz","Oz"],
       "mixing":[[1.,0.,1.],[0.,1.,1.],[-1.,0.,-1.],[0.,-1.,-1.],[0.,0.,0.]],
       "removed":[1],"classNames":["Brain","Eye","Muscle"],
       "classProbabilities":[[.9,.05,.05],[.05,.9,.05],[.1,.1,.8]]}
    return a

class MatchConfidence(unittest.TestCase):
    def test_sign_permutation_normalization(self):
        a=test_maps()
        b=deepcopy(a)
        b["mixing"]=[[-row[1],row[0],row[2]] for row in b["mixing"]]
        b["removed"]=[2]
        b["classProbabilities"]=[a["classProbabilities"][1],a["classProbabilities"][0],
                                  a["classProbabilities"][2]]
        z=assignment(a,b)
        self.assertEqual(z["lowCorrelationComponents"],0)
        self.assertEqual(z["assignmentConflicts"],0)
        self.assertEqual(z["weakWithRemovalDiscordance"],0)
        self.assertAlmostEqual(min(x["absoluteCorrelation"] for x in z["allPairs"]),1.,places=9)

    def test_no_component_relabeling_claim_from_class_disagreement(self):
        a=test_maps();b=deepcopy(a)
        b["classProbabilities"][0]=[0.01,.98,.01]
        z=assignment(a,b)
        self.assertEqual(z["lowCorrelationComponents"],0)
        self.assertEqual(z["weakWithICLabelClassChange"],0)
        self.assertTrue(any(x["classDiscordant"] for x in z["allPairs"]))

    def test_assignment_conflict_is_not_called_intrinsic_weakness(self):
        a=test_maps();b=deepcopy(a)
        base=np.asarray(a["mixing"],dtype=float)
        v0=base[:,0];v1=base[:,1]
        # Two seed-0 sources both prefer the same other-seed component;
        # one-to-one assignment can force a weak match.
        a["mixing"]=np.column_stack([v0,.8*v0+.6*v1,base[:,2]]).tolist()
        b["mixing"]=np.column_stack([.9*v0+.43*v1,v1,base[:,2]]).tolist()
        z=assignment(a,b)
        self.assertGreaterEqual(z["assignmentConflicts"],1)
        self.assertGreaterEqual(z["lowCorrelationComponents"],1)
        self.assertTrue(any(row["weakType"]=="assignment-conflict" for row in z["weakPairs"]))

    def test_intrinsic_weakness_is_distinguished(self):
        a=test_maps();b=deepcopy(a)
        # Duplicate source maps -> underdetermined assignment is still audited.
        b["mixing"]=[[row[0],row[0],row[2]] for row in b["mixing"]]
        z=assignment(a,b)
        self.assertTrue(any(x["weakType"]=="intrinsic-weak" for x in z["weakPairs"]))

    def test_missing_pz_label_is_unavailable_not_effect(self):
        a=test_maps();b=deepcopy(a)
        a["channelLabels"][3]="CPz"
        b["channelLabels"][3]="CPz"
        z=assignment(a,b)
        self.assertTrue(all(x["normalizedPzMapLoadingOnly"] is None for x in z["allPairs"]))

    def test_bad_montage_fails_closed(self):
        a=test_maps();b=deepcopy(a)
        b["channelLabels"][0]="T7"
        with self.assertRaisesRegex(ValueError,"electrode"):
            assignment(a,b)

    def test_json_artifact_is_lossless_and_serializable(self):
        import json
        z=assignment(test_maps(),deepcopy(test_maps()))
        blob=json.dumps(z,allow_nan=False)
        self.assertEqual(json.loads(blob)["lowCorrelationComponents"],0)

    def test_missing_class_probabilities_are_unknown(self):
        a=test_maps();b=deepcopy(a)
        b["classProbabilities"]=[]
        z=assignment(a,b)
        self.assertTrue(all(x["classDiscordant"] is None for x in z["allPairs"]))

    def test_cycle_transitivity(self):
        good=_cycle({1:2,2:1,3:3},{2:3,1:1,3:2},{1:3,2:1,3:2})
        self.assertEqual((good["assessed"],good["consistent"]),(3,3))
        bad=_cycle({1:2,2:1},{2:3,1:1},{1:1,2:3})
        self.assertEqual(bad["consistent"],0)

    def test_waveform_uses_gain_corrected_full_pz(self):
        tm=np.linspace(-200.,800.,101)
        m={"gain":2.,"correctedDifferenceEstimate":.5,
           "science":{"channel":"Pz","timeMs":tm.tolist(),
              "conditions":[{"meanWaveform":[1.]*len(tm)},
                            {"meanWaveform":[0.]*len(tm)}]}}
        t,w=_waveform(m)
        self.assertTrue(np.array_equal(t,tm))
        self.assertAlmostEqual(float(w.mean()),.5)
        m["correctedDifferenceEstimate"]=1.0
        with self.assertRaisesRegex(ValueError,"parity"):
            _waveform(m)

if __name__=="__main__":
    unittest.main()
