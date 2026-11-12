"""ICA weak-component audit must align maps before comparing removal decisions."""
import sys
import unittest
from copy import deepcopy
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"benchmarks"))
from diagnose_ica_component_outliers import map_pairs


def science():
    a={"nComponents":3,"nChannels":5,
       "channelLabels":["Fp1","Fz","Cz","Pz","Oz"],
       "mixing":[[1,0,1],[0,1,1],[0,0,-1],[0,0,0],[-1,-1,-1]],
       "removed":[1],"classNames":["Brain","Eye","Muscle"],
       "classProbabilities":[[.9,.05,.05],[.05,.9,.05],[.1,.1,.8]]}
    return {"decompositions":[a]}


class WeakIC(unittest.TestCase):
    def test_sign_and_permutation_do_not_fake_instability(self):
        a=science(); b=deepcopy(a)
        d=b["decompositions"][0]
        d["mixing"]=[[r[1],-r[0],r[2]] for r in d["mixing"]]
        d["removed"]=[2]
        d["classProbabilities"]=[d["classProbabilities"][1],
                                 d["classProbabilities"][0],
                                 d["classProbabilities"][2]]
        z=map_pairs(a,b)
        self.assertAlmostEqual(z["minimumCorrelation"],1.0,places=10)
        self.assertEqual(z["belowPoint90"],0)
        self.assertEqual(z["matchedRemovedJaccard"],1.0)

    def test_removed_decision_disagreement_remains_visible(self):
        a=science();b=deepcopy(a)
        b["decompositions"][0]["removed"]=[2]
        z=map_pairs(a,b)
        self.assertEqual(z["matchedRemovedJaccard"],0.0)

    def test_weak_mixing_map_is_reported_not_hidden_by_median(self):
        a=science();b=deepcopy(a)
        # A third map changed but two other maps are constant.
        b["decompositions"][0]["mixing"]=[[row[0],row[1],replacement]
            for row,replacement in zip(b["decompositions"][0]["mixing"],[0,1,0,-1,0])]
        z=map_pairs(a,b)
        self.assertGreater(z["belowPoint90"],0)
        self.assertEqual(z["lowCorrelationRemovalDiscordant"],0)
        self.assertTrue(z["medianCorrelation"]>=z["minimumCorrelation"])

    def test_channel_mismatch_is_not_aligned_by_position(self):
        a=science();b=deepcopy(a)
        b["decompositions"][0]["channelLabels"][0]="T7"
        with self.assertRaisesRegex(ValueError,"montages"):
            map_pairs(a,b)

    def test_invalid_iclabel_dimensions_fail_closed(self):
        a=science();b=deepcopy(a)
        b["decompositions"][0]["classProbabilities"]=[[1,0,0]]
        with self.assertRaisesRegex(ValueError,"ICLabel"):
            map_pairs(a,b)


if __name__=="__main__":
    unittest.main()
