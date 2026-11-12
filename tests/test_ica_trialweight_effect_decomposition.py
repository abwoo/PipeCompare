"""Development ERP CORE P300 candidate effect shift decomposition tests."""
import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"benchmarks"))
from ica_trialweight_effect_decomposition import bootstrap_decomposition,weighted_mean


def event_data(oldw=(1,1,1,1,1,1),neww=(1,2,2,1,1,1),effect=0):
    a=np.array([1.,2.,3.,4.,5.,6.])
    b=np.array([2.,4.,6.,8.,10.,12.])
    old=[(a,np.asarray(oldw,float),np.ones(6,dtype=bool)),
         (b,np.ones(6),np.ones(6,dtype=bool))]
    new=[(a+effect,np.asarray(neww,float),np.ones(6,dtype=bool)),
         (b,np.ones(6),np.ones(6,dtype=bool))]
    return old,new


class ScientificEffectReconciliation(unittest.TestCase):
    def test_identical_scores_weight_change_only(self):
        old,new=event_data()
        r=bootstrap_decomposition(old,new,1.,1.,draws=1000,seed=515)
        self.assertTrue(r["allSourceScoresIdentical"])
        self.assertAlmostEqual(r["signalTermRaw"],0.,places=12)
        self.assertAlmostEqual(r["gainTermCorrected"],0.,places=12)
        self.assertNotAlmostEqual(r["weightTermRaw"],0.)
        self.assertAlmostEqual(r["netCorrectedDifference"],r["weightTermRaw"])
        self.assertEqual(r["pairedBootstrap"]["drawsValid"],1000)

    def test_same_weights_different_signal(self):
        old,new=event_data(neww=(1,1,1,1,1,1),effect=2.)
        r=bootstrap_decomposition(old,new,1.,1.,draws=1000,seed=515)
        self.assertFalse(r["allSourceScoresIdentical"])
        self.assertAlmostEqual(r["weightTermRaw"],0.,places=12)
        self.assertAlmostEqual(r["signalTermRaw"],2.,places=12)

    def test_gain_difference_only(self):
        old,new=event_data(neww=(1,1,1,1,1,1))
        r=bootstrap_decomposition(old,new,1.5,.75,draws=1000,seed=515)
        self.assertAlmostEqual(r["signalTermRaw"],0.,places=12)
        self.assertAlmostEqual(r["weightTermRaw"],0.,places=12)
        self.assertAlmostEqual(r["netCorrectedDifference"],r["gainTermCorrected"],places=12)

    def test_paired_bootstrap_is_repeatable(self):
        old,new=event_data()
        p=bootstrap_decomposition(old,new,1.,1.,draws=250,seed=515)
        q=bootstrap_decomposition(old,new,1.,1.,draws=250,seed=515)
        self.assertEqual(p["pairedBootstrap"]["ci95"],q["pairedBootstrap"]["ci95"])

    def test_invalid_gain_fails(self):
        old,new=event_data()
        with self.assertRaisesRegex(ValueError,"gain"):
            bootstrap_decomposition(old,new,0.,1.,draws=100)

    def test_bad_source_alignment_fails(self):
        old,new=event_data()
        new[0]=(new[0][0][:-1],new[0][1][:-1],new[0][2][:-1])
        with self.assertRaisesRegex(ValueError,"axes"):
            bootstrap_decomposition(old,new,1.,1.,draws=100)

    def test_wrong_conditions_fail(self):
        old,new=event_data()
        with self.assertRaisesRegex(ValueError,"two ERP"):
            bootstrap_decomposition(old[:1],new[:1],1.,1.,draws=100)

    def test_no_positive_weights_fails_closed(self):
        with self.assertRaisesRegex(ValueError,"empty weighted"):
            weighted_mean(np.ones(3),np.zeros(3))


if __name__=="__main__":
    unittest.main()
