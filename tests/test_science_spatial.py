"""ICA maps must be matched independent of sign and component ordering."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"benchmarks"))
from analyze_science_recorder import ic_match,waveform_compare


def snap(mix,removed,curve=None):
    if curve is None:curve=[0,1,2,3,4,5,6,7,8,9]
    return {"ica":{"mixing":mix,"removedIndices":removed,"channelIndices":[1,2,3,4]},
            "waveform":{"timesMs":list(range(10)),"gainApplied":False,"difference":curve}}


class SpatialComparison(unittest.TestCase):
    def test_permutation_and_sign_are_matched(self):
        # Two orthogonal-ish nonconstant independent components.
        a=snap([[1,0],[0,1],[0,0],[-1,-1]],[1])
        b=snap([[0,-1],[-1,0],[0,0],[1,1]],[2])
        z=ic_match(a,b)
        self.assertAlmostEqual(z["medianAbsTopoCorrelation"],1,places=8)
        self.assertAlmostEqual(z["removedMatchJaccard"],1,places=8)
    def test_true_removal_difference_is_preserved(self):
        a=snap([[1,0],[0,1],[0,0],[-1,-1]],[1])
        b=snap([[1,0],[0,1],[0,0],[-1,-1]],[2])
        self.assertEqual(ic_match(a,b)["removedMatchJaccard"],0)
    def test_waveform_correlation_and_difference(self):
        a=snap([[1],[0],[0],[-1]],[],list(range(10)))
        b=snap([[1],[0],[0],[-1]],[],list(range(10)))
        z=waveform_compare(a,b)
        self.assertAlmostEqual(z["fullEpochPearson"],1)
        self.assertEqual(z["fullEpochRmsDifferenceUv"],0)
    def test_waveform_time_mismatch_fails(self):
        a=snap([[1],[0],[0],[-1]],[])
        b=snap([[1],[0],[0],[-1]],[])
        b["waveform"]["timesMs"]=[1]*10
        with self.assertRaisesRegex(ValueError,"time axes"):
            waveform_compare(a,b)


if __name__=="__main__":
    unittest.main()
