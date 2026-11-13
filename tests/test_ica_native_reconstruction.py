"""Native EEG.icawinv, not pinv(W), must model actual pop_subcomp."""
import sys
import unittest
from copy import deepcopy
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"benchmarks"))
from validate_ica_native_reconstruction import native_reconstruction


def example():
    native_first=[2.,4.,3.,5.]
    native_second=[1.,1.,1.,1.]
    pre={"status":"complete","nEvents":4,"srate":100.,
         "eventIds":[1,2,3,4],"conditionIndex":[1,1,2,2],
         "nComponents":2,
         "observedPzUv":[3.5,5.5,4.5,6.5],
         "nativeICComponentUv":[native_first,native_second],
         "nativeResidualPzUv":[.5,.5,.5,.5],
         "nativeMaxAbsResidualUv":.5,
         "nativeVsPinvMixingMaxAbs":.2,
         "componentUv":[[1.8,3.8,2.8,4.8],[1.7,1.7,1.7,1.7]]}
    post={"status":"complete","nEvents":4,"srate":100.,
          "eventIds":[1,2,3,4],"conditionIndex":[1,1,2,2],
          "observedPzUv":[1.,1.,1.,1.]}
    d={"nComponents":2,"removed":[1]}
    science={"conditions":[
       {"referenceIds":[1,2],"retainedIds":[1,2],"trialWeights":[1.,1.]},
       {"referenceIds":[3,4],"retainedIds":[3,4],"trialWeights":[1.,1.]}]}
    return pre,post,d,science


class NativeEEGLABProof(unittest.TestCase):
    def test_real_keep_reconstruction_can_pass_when_old_pinv_rule_fails(self):
        z=native_reconstruction(*example())
        self.assertTrue(z["nativeRealEEGLABProjectionVerified"])
        self.assertAlmostEqual(z["actualNativeKeptProjectionMaxErrorUv"],0.)
        self.assertAlmostEqual(z["nativePzSubspaceResidualMaxAbsUv"],.5)
        self.assertGreater(z["originalPseudoinverseHypothesisMaxErrorUv"],.5)
        self.assertFalse(z["fullFinalP300CausalityOrSuperiorAlgorithmShown"])

    def test_noop_when_zero_comps_removed_preserves_physical_pz(self):
        pre,post,d,s=example()
        d["removed"]=[]
        post["observedPzUv"]=pre["observedPzUv"][:]
        x=native_reconstruction(pre,post,d,s)
        self.assertEqual(x["numberRemovedICs"],0)

    def test_rank_deficient_residual_not_erased_from_provenance(self):
        pre,post,d,s=example()
        pre["nativeResidualPzUv"]=[2.,2.,2.,2.]
        pre["nativeMaxAbsResidualUv"]=2.
        pre["observedPzUv"]=[5.,7.,6.,8.]
        x=native_reconstruction(pre,post,d,s)
        self.assertEqual(x["nativePzSubspaceResidualMaxAbsUv"],2.)

    def test_unexplained_new_post_signal_must_fail_closed(self):
        pre,post,d,s=example()
        post["observedPzUv"][0]+=0.0002
        with self.assertRaisesRegex(ValueError,"fails"):
            native_reconstruction(pre,post,d,s)

    def test_missing_original_native_mixing_fails(self):
        pre,post,d,s=example()
        del pre["nativeICComponentUv"]
        with self.assertRaisesRegex(ValueError,"unavailable"):
            native_reconstruction(pre,post,d,s)

    def test_wrong_native_subspace_residual_fails(self):
        pre,post,d,s=example()
        pre["nativeResidualPzUv"]=[0.]*4
        pre["nativeMaxAbsResidualUv"]=0.
        with self.assertRaisesRegex(ValueError,"does not reproduce"):
            native_reconstruction(pre,post,d,s)

    def test_wrong_event_identity_fails(self):
        pre,post,d,s=example()
        post["eventIds"][0]=99
        with self.assertRaisesRegex(ValueError,"not preserved"):
            native_reconstruction(pre,post,d,s)

    def test_wrong_removed_ic_identity_fails(self):
        pre,post,d,s=example()
        d["removed"]=[2]
        with self.assertRaisesRegex(ValueError,"fails"):
            native_reconstruction(pre,post,d,s)

    def test_wrong_rate_fails(self):
        pre,post,d,s=example()
        post["srate"]=250.
        with self.assertRaisesRegex(ValueError,"sampling rate"):
            native_reconstruction(pre,post,d,s)

    def test_changing_trial_weights_not_equal_to_source_change(self):
        pre,post,d,s=example()
        s["conditions"][0]["trialWeights"]=[1.,3.]
        z=native_reconstruction(pre,post,d,s)
        self.assertAlmostEqual(z["afterStageOriginalFinalTrialsUv"],0.)
        self.assertNotEqual(z["beforeStageOriginalFinalTrialsUv"],-1.)

if __name__=="__main__":
    unittest.main()
