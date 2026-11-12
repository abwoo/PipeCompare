"""Scientifically correct one-to-one ICA map matching across fit starts."""
import sys
import unittest
from copy import deepcopy
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"benchmarks"))

from compare_ica_topography import alignment


def pairs():
    baseline={"available":True,"channels":[1,2,3],
      "maps":[[1,0,0],[0,1,0],[0,0,1]],"removedIndices":[1,3]}
    # Matched seed maps: new columns = [-old3, old1, -old2].
    shuffled={"available":True,"channels":[1,2,3],
      "maps":[[0,1,0],[0,0,-1],[-1,0,0]],"removedIndices":[1,2]}
    return baseline, shuffled


class TopographicMatching(unittest.TestCase):
    def test_permutation_and_polarity_do_not_create_instability(self):
        a,b=pairs()
        z=alignment(a,b)
        self.assertTrue(z["available"])
        self.assertAlmostEqual(z["medianMatchedAbsTopologyCorrelation"],1)
        self.assertAlmostEqual(z["minMatchedAbsTopologyCorrelation"],1)
        self.assertAlmostEqual(z["alignedRemovedJaccard"],1)
        self.assertEqual(z["matchedColumns"],[[1,2],[2,3],[3,1]])

    def test_changed_removed_components_are_not_hidden_by_map_matching(self):
        a,b=pairs()
        b["removedIndices"]=[1]
        z=alignment(a,b)
        self.assertEqual(z["alignedRemovedJaccard"],0.5)

    def test_matlab_singleton_removed_indices_are_not_missing(self):
        a,b=pairs()
        a["removedIndices"]=3
        b["removedIndices"]=1.0
        z=alignment(a,b)
        self.assertTrue(z["available"])
        self.assertEqual(z["removedCountBaseline"],1)
        self.assertEqual(z["alignedRemovedJaccard"],1)

    def test_nothing_removed_cannot_be_called_perfect(self):
        a,b=pairs()
        a["removedIndices"]=[];b["removedIndices"]=[]
        z=alignment(a,b)
        self.assertTrue(z["zeroRemovalBothSeeds"])
        self.assertIsNone(z["alignedRemovedJaccard"])

    def test_channel_order_reconciled_by_explicit_channel_indices(self):
        a,b=pairs()
        b["channels"]=[3,1,2]
        b["maps"]=[b["maps"][2],b["maps"][0],b["maps"][1]]
        z=alignment(a,b)
        self.assertTrue(z["available"])
        self.assertAlmostEqual(z["medianMatchedAbsTopologyCorrelation"],1)

    def test_channel_set_drift_abstains(self):
        a,b=pairs()
        b["channels"]=[1,2,4]
        z=alignment(a,b)
        self.assertFalse(z["available"])

    def test_missing_or_invalid_maps_abstain(self):
        a,b=pairs()
        b["available"]=False
        self.assertFalse(alignment(a,b)["available"])
        a,b=pairs()
        b["maps"]=[[1,0],[0,1],[1,1]]
        self.assertFalse(alignment(a,b)["available"])

    def test_invalid_component_labels_abstain(self):
        a,b=pairs()
        b["removedIndices"]=[0,3]
        self.assertFalse(alignment(a,b)["available"])
        a,b=pairs()
        b["removedIndices"]=[1,1]
        self.assertFalse(alignment(a,b)["available"])

if __name__=="__main__":
    unittest.main()
