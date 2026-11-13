"""Eventwise post-IC removal Pz scientific bridge must fail closed."""
import sys
import unittest
from copy import deepcopy
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"benchmarks"))
from validate_ica_postremoval import verify_bridge, true_event_table


def data():
    pre={"status":"complete","sourceStage":"immediately-before-icremove",
         "eventIds":[1,2,3,4],"conditionIndex":[1,1,2,2],"nEvents":4,
         "nComponents":2,"srate":100.,
         "componentUv":[[1.,2.,1.,2.],[.5,.5,.5,.5]],
         "observedPzUv":[1.5,2.5,1.5,2.5]}
    post={"status":"complete","sourceStage":"immediately-after-icremove",
          "eventIds":[1,2,3,4],"conditionIndex":[1,1,2,2],
          "observedPzUv":[1.,2.,1.,2.],"nEvents":4,"srate":100.}
    decomposition={"nComponents":2,"removed":[2]}
    science={"conditions":[
       {"referenceIds":[1,2],"retainedIds":[1,2],"trialWeights":[1.,1.]},
       {"referenceIds":[3,4],"retainedIds":[3,4],"trialWeights":[1.,1.]}]}
    return pre,post,decomposition,science


class PostICAStage(unittest.TestCase):
    def test_actual_removed_component_conserves_physical_pz(self):
        d=verify_bridge(*data())
        self.assertEqual(d["removedComponentCount"],1)
        self.assertTrue(d["directICStepLinearBridgeValid"])
        self.assertAlmostEqual(d["linearPzRemovalMaxResidualUv"],0.)
        self.assertFalse(d["downstreamCausalityOrAccuracyEstablished"])

    def test_no_removed_components_holds_original_signal(self):
        a,b,c,s=data()
        c["removed"]=[]
        b["observedPzUv"]=a["observedPzUv"][:]
        d=verify_bridge(a,b,c,s)
        self.assertEqual(d["removedComponentCount"],0)

    def test_real_signal_change_without_removed_component_fails(self):
        a,b,c,s=data()
        c["removed"]=[]
        with self.assertRaisesRegex(ValueError,"does not match"):
            verify_bridge(a,b,c,s)

    def test_bad_removed_component_index_fails(self):
        a,b,c,s=data()
        c["removed"]=[3]
        with self.assertRaisesRegex(ValueError,"Removed IC"):
            verify_bridge(a,b,c,s)

    def test_nonzero_reconstruction_residual_does_not_corrupt_delta(self):
        a,b,c,s=data()
        a["observedPzUv"]=[x+.8 for x in a["observedPzUv"]]
        b["observedPzUv"]=[x+.8 for x in b["observedPzUv"]]
        d=verify_bridge(a,b,c,s)
        self.assertAlmostEqual(d["linearPzRemovalMaxResidualUv"],0.)

    def test_different_event_order_is_permutation_invariant(self):
        a,b,c,s=data()
        b["eventIds"]=[2,1,4,3]
        b["conditionIndex"]=[1,1,2,2]
        b["observedPzUv"]=[2.,1.,2.,1.]
        d=verify_bridge(a,b,c,s)
        self.assertTrue(d["sourceEventIDsPreserved"])

    def test_missing_event_fails_even_if_counts_match(self):
        a,b,c,s=data()
        b["eventIds"][0]=99
        with self.assertRaisesRegex(ValueError,"physical event"):
            verify_bridge(a,b,c,s)

    def test_wrong_event_condition_fails(self):
        a,b,c,s=data()
        b["conditionIndex"][0]=2
        with self.assertRaisesRegex(ValueError,"physical event"):
            verify_bridge(a,b,c,s)

    def test_stage_not_completed_fails(self):
        a,b,c,s=data()
        b={"status":"failed","reason":"Missing sampling"}
        with self.assertRaisesRegex(ValueError,"Missing sampling"):
            verify_bridge(a,b,c,s)

    def test_wrong_original_final_retained_id_fails(self):
        a,b,c,s=data()
        s["conditions"][0]["retainedIds"]=[1,999]
        with self.assertRaisesRegex(ValueError,"retained"):
            verify_bridge(a,b,c,s)

    def test_bad_sample_rate_fails(self):
        a,b,c,s=data()
        b["srate"]=250.
        with self.assertRaisesRegex(ValueError,"sampling rate"):
            verify_bridge(a,b,c,s)

    def test_numerical_threshold_is_frozen(self):
        a,b,c,s=data()
        b["observedPzUv"][0]+=2.0e-4
        with self.assertRaisesRegex(ValueError,"does not match"):
            verify_bridge(a,b,c,s)

    def test_nonfinite_stage_fails(self):
        a,b,c,s=data()
        b["observedPzUv"][0]=float("nan")
        with self.assertRaisesRegex(ValueError,"malformed"):
            true_event_table(b,"post")

if __name__=="__main__":
    unittest.main()
