"""Confirmed source HMM support does not imply downstream eligibility."""
import sys
import unittest
from copy import deepcopy
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"benchmarks"))
from validate_subject8_extended_pipeline import analyze

def support(n,k,condition):
    return {"condition":condition,"referenceN":n,
        "effectiveUnflaggedN":k,
        "scores":[float((-1)**i) for i in range(n)],
        "weights":[1.0]*n}
def goodarm():
    return {"available":True,"differenceAmplitude":0.0,
        "gain":1.0,"reportedDifferenceSE":1.0,
        "support":[support(30,6,"target"),support(126,24,"standard")]}
def fixture():
    e=goodarm()
    return {"subject":8,"developmentOnly":True,"extended460Seconds":True,
      "originalFrozenVerdict":"PARTIAL","promotionAuthorized":False,
      "source":{"subject":8,"duration_seconds":460,"original_frozen_duration_seconds":300,
          "exploratory_duration_override":True,
          "protocol_sha256":"p","source_set_sha256":"s",
          "source_fdt_sha256":"f","source_events_sha256":"v"},
      "completed":True,"errors":[],"notRun":0,
      "replayMaxError":0,"exportMaxError":0,
      "protocolSha256":"p","reference":{"available":True},
      "erp":e,"minimal":goodarm(),"common":goodarm(),
      "selectedSupport":e["support"]}
def original_hmm():
    def row(n,k): return {"referenceN":n,"unflaggedN":k}
    return {"subject":8,"originalSeconds":300,"fullSeconds":460,
      "fullWindowUnflaggedFloorMet":True,
      "baseline":[row(22,5),row(95,20)],
      "extended":[row(30,6),row(126,24)],
      "sourceSetSha256":"s","sourceFdtSha256":"f","sourceEventsSha256":"v"}

class FullPipelineGate(unittest.TestCase):
    def test_good_pair_is_still_only_one_dev_subject(self):
        r=analyze(fixture(),original_hmm())
        self.assertTrue(r["allThreeArmsScientificallyEligible"])
        self.assertFalse(r["productionPromotionAuthorized"])
        self.assertFalse(r["statisticalGeneralizationClaimSupported"])
        self.assertEqual(r["original20SubjectValidPairs"],15)

    def test_source_support_does_not_rescue_downstream_trial_rejection(self):
        r=fixture()
        r["minimal"]["support"][0]["effectiveUnflaggedN"]=4
        o=analyze(r,original_hmm())
        self.assertFalse(o["allThreeArmsScientificallyEligible"])
        self.assertIn("effective-support-below-five:target",o["arms"][1]["failureReasons"])

    def test_absent_standard_cannot_be_counted_valid(self):
        r=fixture()
        r["erp"]["available"]=False
        r["selectedSupport"]=None
        r["replayMaxError"]=float("nan")
        r["exportMaxError"]=float("nan")
        x=analyze(r,original_hmm())
        self.assertFalse(x["allThreeArmsScientificallyEligible"])
        self.assertIn("candidate-unavailable",x["arms"][0]["failureReasons"])

    def test_changed_hmm_original_counts_fail(self):
        h=original_hmm()
        h["baseline"][0]["unflaggedN"]=6
        with self.assertRaises(ValueError):
            analyze(fixture(),h)

    def test_wrong_source_sha_fails(self):
        h=original_hmm()
        h["sourceEventsSha256"]="wrong"
        with self.assertRaises(ValueError):
            analyze(fixture(),h)

    def test_wrong_cohort_rejected(self):
        r=fixture()
        r["subject"]=24
        with self.assertRaises(ValueError):
            analyze(r,original_hmm())

if __name__=="__main__":
    unittest.main()
