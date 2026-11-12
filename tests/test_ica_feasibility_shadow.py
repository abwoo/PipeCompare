"""Frozen scientific gates are necessary for all-seed ICA robust ranking."""
import sys
import unittest
from copy import deepcopy
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"benchmarks"))
from ica_feasibility_shadow import shadow


def row(candidate, seed, amp=2.0, eff=(10,25)):
    return {"seed":seed,"candidateId":candidate,"hasMeasurement":True,
            "correctedDifferenceEstimate":amp,"gain":1.0,
            "effectiveUnflagged":list(eff)}


def fixture():
    return {"subject":21,"completed":True,"errors":[],
      "icaSeedCheck":{"ids":[4,25],"recommended":[25,4,4],"same":False,
      "objective":[[1.12,1.13,1.14],[1.10,1.35,1.36]],
      "rankStatus":[["feasible"]*3,["feasible"]*3],
      "strata":[["standard"]*3,["standard"]*3],
      "candidateMeasurements":[
          [row(4,j,2+.01*j) for j in range(3)],
          [row(25,j,2+.02*j) for j in range(3)]]}}


class ShadowTest(unittest.TestCase):
    def test_common_feasible_worst_case(self):
        z=shadow(fixture())
        self.assertEqual(z["shadowCandidate"],4)
        self.assertFalse(z["scientificPromotionAuthorized"])
        self.assertAlmostEqual(z["shadowP3Range"],.02)

    def test_matlab_column_major_cell_array_serialization(self):
        x=fixture()
        check=x["icaSeedCheck"]
        # MATLAB jsonencode flattens Nx3 cell matrices in column-major order:
        # seed0/candidate0, seed0/candidate1, seed1/candidate0, ...
        for key in ("rankStatus","strata","candidateMeasurements"):
            rows=check[key]
            check[key]=[rows[i][seed] for seed in range(3) for i in range(2)]
        z=shadow(x)
        self.assertEqual(z["shadowCandidate"],4)
        self.assertEqual(z["eligibleCount"],2)
        self.assertEqual(z["sourceRecommended"],[25,4,4])

    def test_matlab_flattened_rejection_remains_per_correct_candidate(self):
        x=fixture()
        check=x["icaSeedCheck"]
        check["rankStatus"][0][2]="rejected"
        for key in ("rankStatus","strata","candidateMeasurements"):
            rows=check[key]
            check[key]=[rows[i][seed] for seed in range(3) for i in range(2)]
        z=shadow(x)
        self.assertEqual(z["shadowCandidate"],25)

    def test_excludes_best_score_when_rejected_in_one_seed(self):
        x=fixture()
        x["icaSeedCheck"]["rankStatus"][0][2]="rejected"
        z=shadow(x)
        self.assertEqual(z["shadowCandidate"],25)
        self.assertEqual(z["eligibleCount"],1)

    def test_abstains_when_original_trial_floor_missing(self):
        x=fixture()
        for measurements in x["icaSeedCheck"]["candidateMeasurements"]:
            measurements[1]["effectiveUnflagged"]=[4,10]
        z=shadow(x)
        self.assertEqual(z["eligibleCount"],0)
        self.assertIsNone(z["shadowCandidate"])
        self.assertEqual(z["reason"],"no-common-feasible-candidate")

    def test_no_cross_strata_comparison(self):
        x=fixture()
        x["icaSeedCheck"]["strata"][0][1]="other"
        z=shadow(x)
        self.assertEqual(z["shadowCandidate"],25)

    def test_tie_abstains_instead_of_candidate_order(self):
        x=fixture()
        x["icaSeedCheck"]["objective"]=[[1,1,1],[1,1,1]]
        z=shadow(x)
        self.assertIsNone(z["shadowCandidate"])
        self.assertEqual(z["reason"],"objective-tie-requires-predefined-tiebreak")

    def test_missing_measurement_never_promotes(self):
        x=fixture()
        x["icaSeedCheck"]["candidateMeasurements"][0][1]["hasMeasurement"]=False
        z=shadow(x)
        self.assertEqual(z["shadowCandidate"],25)

    def test_misaligned_candidate_is_error_or_ineligible(self):
        x=fixture()
        x["icaSeedCheck"]["candidateMeasurements"][0][0]["candidateId"]=9
        z=shadow(x)
        self.assertFalse(z["audit"][0]["eligible"])

    def test_incomplete_run_rejected(self):
        x=fixture()
        x["completed"]=False
        with self.assertRaises(ValueError):
            shadow(x)


if __name__=="__main__":
    unittest.main()
