"""Real-ICA scientific gate rechecker: false success and missingness regressions."""
import sys
import unittest
from copy import deepcopy
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"benchmarks"))
from validate_ica_robust_real import validate


def sample():
    ms = [
        {"seed":j, "candidateId":4, "effectiveUnflagged":[8,9],
         "differenceEstimate":1.2, "gain":0.8, "correctedDifferenceEstimate":1.5}
        for j in range(3)
    ]
    data = {"subject":21,"completed":True,"errors":[],"notRun":0,
            "icaSeedCheck":{
                "recommended":[25,4,4],"same":False,
                "measurements":[{"seed":i} for i in range(3)],
                "ids":[4,25],
                "objective":[[1.12,1.127,1.133],[1.121,1.350,1.354]],
                "statusBySeed":[["feasible"]*3,["feasible"]*3],
                "stratumBySeed":[["P3"]*3,["P3"]*3],
                "supportPassBySeed":[[True]*3,[True]*3],
                "robust":{"status":"selected","candidateId":4,"worstObjective":1.133,
                    "originalRecommendationUnchanged":True,
                    "allSeedsScientificallyFeasible":True,
                    "allSeedsSupportFloorPassed":True,
                    "reason":"all-seed-Rank-run-feasible-and-finite",
                    "seedMeasurements":ms}
            }}
    c=data["icaSeedCheck"]
    c["candidateSeedAudit"]=[
        {"candidateId":idx,"objective":objective,
         "status":status,"supportPass":support,"stratum":group}
        for idx,objective,status,support,group in zip(
            c["ids"],c["objective"],c["statusBySeed"],
            c["supportPassBySeed"],c["stratumBySeed"])]
    return data


class GateAuditTests(unittest.TestCase):
    def test_valid_feasible_candidate(self):
        x=validate(sample())
        self.assertEqual(x["robustId"],4)
        self.assertFalse(x["scientificGateChanged"])
        self.assertFalse(x["confirmation"])

    def test_single_failed_shape_gate_invalidates_recommendation(self):
        x=sample()
        x["icaSeedCheck"]["candidateSeedAudit"][0]["status"][2]="rejected"
        with self.assertRaisesRegex(ValueError,"scientific gate"):
            validate(x)

    def test_single_insufficient_support_gate_invalidates_recommendation(self):
        x=sample()
        x["icaSeedCheck"]["candidateSeedAudit"][0]["supportPass"][2]=False
        with self.assertRaisesRegex(ValueError,"scientific gate"):
            validate(x)

    def test_support_below_five_not_accepted(self):
        x=sample()
        x["icaSeedCheck"]["robust"]["seedMeasurements"][0]["effectiveUnflagged"][0]=4
        with self.assertRaisesRegex(ValueError,"below fixed floor"):
            validate(x)

    def test_ambiguous_strata_cannot_be_ranked_together(self):
        x=sample()
        x["icaSeedCheck"]["candidateSeedAudit"][1]["stratum"]=["OTHER"]*3
        with self.assertRaisesRegex(ValueError,"incomparable"):
            validate(x)

    def test_all_failed_support_abstain(self):
        x=sample()
        x["icaSeedCheck"]["candidateSeedAudit"][0]["supportPass"]=[False]*3
        x["icaSeedCheck"]["candidateSeedAudit"][1]["supportPass"]=[False]*3
        x["icaSeedCheck"]["robust"].update({
            "status":"abstain","candidateId":None,"reason":"no-common-feasible-candidate",
            "worstObjective":None,"seedMeasurements":[],
            "allSeedsScientificallyFeasible":False,"allSeedsSupportFloorPassed":False})
        z=validate(x)
        self.assertEqual(z["robustState"],"abstain")

    def test_scientific_failed_but_selected_fails(self):
        x=sample()
        x["icaSeedCheck"]["candidateSeedAudit"][0]["status"]=["failed"]*3
        x["icaSeedCheck"]["candidateSeedAudit"][1]["status"]=["rejected"]*3
        with self.assertRaisesRegex(ValueError,"scientific gate"):
            validate(x)

    def test_candidate_with_bad_gain_fails(self):
        x=sample()
        x["icaSeedCheck"]["robust"]["seedMeasurements"][1]["gain"]=0
        with self.assertRaisesRegex(ValueError,"invalid gain"):
            validate(x)

    def test_incomplete_run_never_valid(self):
        x=sample();x["completed"]=False
        with self.assertRaisesRegex(ValueError,"incomplete"):
            validate(x)

    def test_flattened_matlab_status_without_explicit_records_fails(self):
        x=sample()
        del x["icaSeedCheck"]["candidateSeedAudit"]
        with self.assertRaisesRegex(ValueError,"explicit seed audit missing"):
            validate(x)

    def test_swapped_scientific_record_order_fails(self):
        x=sample()
        x["icaSeedCheck"]["candidateSeedAudit"].reverse()
        with self.assertRaisesRegex(ValueError,"not aligned"):
            validate(x)

    def test_different_objective_minimum_fails(self):
        x=sample()
        x["icaSeedCheck"]["robust"]["candidateId"]=25
        x["icaSeedCheck"]["robust"]["seedMeasurements"]=[]
        with self.assertRaisesRegex(ValueError,"not unique minimax"):
            validate(x)


if __name__=="__main__":
    unittest.main()
