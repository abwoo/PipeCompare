"""Six ICA starts must be two genuinely disjoint, science-gated triplets."""
import sys
import unittest
from copy import deepcopy
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"benchmarks"))
from validate_ica_independent_groups import validate


def sample():
    ids=[3,5]
    scores=[[1.11,1.12,1.12,1.20,1.21,1.22],
            [1.10,1.11,1.10,1.25,1.26,1.25]]
    by=[]
    for idx,values in zip(ids,scores):
        by.append({"candidateId":idx, "objective":values,
                   "status":["feasible"]*6,"stratum":["ERP_CORE_P3"]*6,
                   "supportPass":[True]*6,
                   "estimationFamily":["weighted"]*6})
    def m(seed,idx):
        return {"seed":seed,"candidateId":idx,"hasMeasurement":True,
                "estimationFamily":"weighted","effectiveUnflagged":[8,15],
                "gain":1.0,"differenceEstimate":1.0,
                "correctedDifferenceEstimate":1.0}
    groups=[]
    for number,choice in enumerate([5,3]):
        start=number*3
        groups.append({"independentGroup":number+1,
                       "seedNumbers":list(range(start,start+3)),
                       "baselineId":5,"baselineFamily":"weighted",
                       "status":"selected","candidateId":choice,
                       "seedMeasurements":[m(j,choice) for j in range(start,start+3)]})
    return {"subject":15,"completed":True,"notRun":0,"errors":[],
            "searchOptions":{"icaSeeds":5},
            "erp":{"differenceAmplitude":1.0,"available":True},
            "icaSeedCheck":{
                "seeds":[1,2,3,4,5],"ids":ids,"objective":scores,
                "recommended":[5]*6,"measurements":[m(j,5) for j in range(6)],
                "candidateSeedAudit":by,
                "independentGroups":{"sourceSeedLabels":list(range(6)),
                      "originalSearchRankingUnchanged":True,
                      "isIndependentHeldoutStudy":False,
                      "sameCandidate":False,
                      "status":"different-across-disjoint-triplets",
                      "groups":groups}}}


class IndependentGroupsTests(unittest.TestCase):
    def test_different_groups_true_independent_starts(self):
        x=validate(sample())
        self.assertEqual(x["groups"][0]["seeds"],[0,1,2])
        self.assertEqual(x["groups"][1]["seeds"],[3,4,5])
        self.assertFalse(x["sameAcrossDisjointIcaTriplets"])
        self.assertFalse(x["newFrozenStudy"])

    def test_replayed_same_triplet_is_not_independent(self):
        x=sample()
        x["icaSeedCheck"]["independentGroups"]["groups"][1]["seedNumbers"]=[0,1,2]
        with self.assertRaisesRegex(ValueError,"overlap"):validate(x)

    def test_non_minimax_must_not_be_reported_as_selected(self):
        x=sample()
        group=x["icaSeedCheck"]["independentGroups"]["groups"][0]
        group["candidateId"]=3
        group["seedMeasurements"]=[dict(z,candidateId=3) for z in group["seedMeasurements"]]
        with self.assertRaisesRegex(ValueError,"not unique"):validate(x)

    def test_failed_seed_support_cannot_pass(self):
        x=sample()
        x["icaSeedCheck"]["candidateSeedAudit"][0]["supportPass"][4]=False
        with self.assertRaisesRegex(ValueError,"not unique"):validate(x)

    def test_no_ica_seed5_rejected(self):
        x=sample()
        x["icaSeedCheck"]["seeds"]=[1,2,3,4]
        with self.assertRaisesRegex(ValueError,"planned"):validate(x)

    def test_five_seeds_not_six_rejected(self):
        x=sample()
        x["searchOptions"]["icaSeeds"]=2
        with self.assertRaisesRegex(ValueError,"six distinct"):validate(x)

    def test_change_weight_estimator_within_group_fails(self):
        x=sample()
        x["icaSeedCheck"]["independentGroups"]["groups"][1]["seedMeasurements"][0]["estimationFamily"]="unweighted"
        with self.assertRaisesRegex(ValueError,"seed or estimator family"):validate(x)

    def test_no_science_evaluable_group_abstains_honestly(self):
        x=sample()
        for r in x["icaSeedCheck"]["candidateSeedAudit"]:
            r["supportPass"][3:6]=[False]*3
        g=x["icaSeedCheck"]["independentGroups"]["groups"][1]
        g["status"]="abstain"
        g["candidateId"]=None
        g["seedMeasurements"]=[]
        p=x["icaSeedCheck"]["independentGroups"]
        p["status"]="at-least-one-group-abstained"
        y=validate(x)
        self.assertEqual(y["status"],"at-least-one-group-abstained")


if __name__=="__main__":
    unittest.main()
