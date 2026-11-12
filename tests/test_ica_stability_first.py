"""Stability-first choice cannot depend on unseen starts."""
import sys
import unittest
from copy import deepcopy
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"benchmarks"))
from ica_stability_first import training_only,choose,compare,unseen
from test_ica_unseen_seeds import fixture


def prepared():
    x=fixture()
    for row in x["icaSeedCheck"]["candidateMeasurements"]:
        for m in row:
            m["sourceUnflagged"]=[12,22]
    return x


class StabilityFirst(unittest.TestCase):
    def test_selects_minimum_training_range(self):
        x=prepared()
        x["icaSeedCheck"]["candidateMeasurements"][0][2]["correctedDifferenceEstimate"]=2.50
        t=training_only(x)
        chosen,_=choose(t)
        self.assertEqual(chosen,25)
    def test_late_seed_data_cannot_reselect(self):
        a=prepared()
        a["icaSeedCheck"]["candidateMeasurements"][1][1]["correctedDifferenceEstimate"]=2.8
        b=deepcopy(a)
        for seed in (3,4,5):
            b["icaSeedCheck"]["candidateMeasurements"][0][seed]["correctedDifferenceEstimate"]=999.
            b["icaSeedCheck"]["objective"][0][seed]=999.
        self.assertEqual(choose(training_only(a)),choose(training_only(b)))
    def test_unseen_support_failure_is_explicit(self):
        x=prepared()
        x["icaSeedCheck"]["candidateMeasurements"][0][4]["effectiveUnflagged"]=[3,20]
        r=unseen(x,4,training_only(x))
        self.assertFalse(r["allUnseenFeasible"])
        self.assertIn("effective-versus-physical-support-failed",r["unseenSeeds"][1]["reasons"])
    def test_unseen_stratum_drift_is_explicit(self):
        x=prepared()
        x["icaSeedCheck"]["strata"][0][4]="different"
        self.assertFalse(unseen(x,4,training_only(x))["allUnseenFeasible"])
    def test_missing_support_must_fail_closed(self):
        x=prepared()
        x["icaSeedCheck"]["candidateMeasurements"][0][5].pop("sourceUnflagged")
        with self.assertRaisesRegex(ValueError,"bound integrity"):
            compare(x)
    def test_objective_and_stability_tie_abstention(self):
        t={"eligible":[{"id":4,"p3Range":0.,"maxObjective":1.},
                       {"id":25,"p3Range":0.,"maxObjective":1.}]}
        self.assertEqual(choose(t),(None,"unresolved-stability-and-objective-tie"))
    def test_training_only_stability_flag(self):
        x=prepared()
        x["icaSeedCheck"]["recommended"]=[4,4,4,25,25,25]
        x["icaSeedCheck"]["same"]=False
        self.assertTrue(training_only(x)["sourceCandidateStable"])
    def test_report_cannot_authorize_promotion(self):
        x=prepared()
        result=compare(x)
        self.assertFalse(result["scientificPromotionAuthorized"])
        self.assertFalse(result["newParticipantValidation"])


if __name__=="__main__":
    unittest.main()
