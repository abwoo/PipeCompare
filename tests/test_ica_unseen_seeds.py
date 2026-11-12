"""Selection on starts 0-2 must be immune to observations on starts 3-5."""
import sys
import unittest
from copy import deepcopy
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"benchmarks"))
from validate_ica_unseen_seeds import inspect


def fixture():
    ids=[4,25]
    objective=[[1.12,1.13,1.14,1.12,1.13,1.15],
               [1.10,1.35,1.36,1.11,1.11,1.11]]
    status=[["feasible"]*6 for _ in ids]
    strata=[["reref"]*6 for _ in ids]
    measurements=[]
    for ident in ids:
        measurements.append([{"seed":j,"candidateId":ident,"hasMeasurement":True,
          "gain":1.0,"correctedDifferenceEstimate":2+.01*j,
          "effectiveUnflagged":[10,20]} for j in range(6)])
    return {"subject":21,"completed":True,"errors":[],
            "icaSeedCheck":{"ids":ids,"objective":objective,"rankStatus":status,
              "strata":strata,"candidateMeasurements":measurements,
              "recommended":[25,4,4,4,25,25],"same":False}}


class UnseenStarts(unittest.TestCase):
    def test_selection_uses_only_first_three(self):
        a=fixture();b=deepcopy(a)
        # Late data would favor baseline under an invalid leaky max-all-six rule.
        b["icaSeedCheck"]["objective"][0][3:]=[999,999,999]
        z0=inspect(a);z1=inspect(b)
        self.assertEqual(z0["trainingShadow"]["shadowCandidate"],4)
        self.assertEqual(z1["trainingShadow"]["shadowCandidate"],4)
        self.assertTrue(z0["shadowUnseen"]["validAllThree"])

    def test_unseen_rejection_invalidates_selected_recipe(self):
        x=fixture()
        x["icaSeedCheck"]["rankStatus"][0][4]="rejected"
        z=inspect(x)
        self.assertEqual(z["trainingShadow"]["shadowCandidate"],4)
        self.assertFalse(z["shadowUnseen"]["validAllThree"])
        self.assertTrue(z["originalBaseline"]["validAllThree"])

    def test_unseen_support_failure_not_rescued(self):
        x=fixture()
        x["icaSeedCheck"]["candidateMeasurements"][0][5]["effectiveUnflagged"]=[4,20]
        z=inspect(x)
        self.assertFalse(z["shadowUnseen"]["validAllThree"])

    def test_matlab_cell_column_major_serialization(self):
        x=fixture()
        c=x["icaSeedCheck"]
        for name in ("rankStatus","strata","candidateMeasurements"):
            rows=c[name]
            c[name]=[rows[i][j] for j in range(6) for i in range(2)]
        z=inspect(x)
        self.assertEqual(z["trainingShadow"]["shadowCandidate"],4)
        self.assertTrue(z["shadowUnseen"]["validAllThree"])

    def test_no_unseen_data_is_needed_for_selection_but_required_for_audit(self):
        x=fixture()
        x["icaSeedCheck"]["candidateMeasurements"][0][4]=None
        z=inspect(x)
        self.assertEqual(z["trainingShadow"]["shadowCandidate"],4)
        self.assertFalse(z["shadowUnseen"]["validAllThree"])

    def test_unseen_nonfinite_objective_fails_closed(self):
        x=fixture()
        x["icaSeedCheck"]["objective"][0][4]=float("nan")
        z=inspect(x)
        self.assertEqual(z["trainingShadow"]["shadowCandidate"],4)
        self.assertFalse(z["shadowUnseen"]["validAllThree"])
        self.assertIn("nonfinite-objective",
                      z["shadowUnseen"]["seeds"][1]["failureReasons"])
        self.assertTrue(z["originalBaseline"]["validAllThree"])

    def test_unseen_cross_stratum_drift_fails_closed(self):
        x=fixture()
        x["icaSeedCheck"]["strata"][0][5]="different-reference"
        z=inspect(x)
        self.assertEqual(z["trainingShadow"]["shadowCandidate"],4)
        self.assertFalse(z["shadowUnseen"]["validAllThree"])
        self.assertIn("stratum-drift",
                      z["shadowUnseen"]["seeds"][2]["failureReasons"])
        self.assertTrue(z["originalBaseline"]["validAllThree"])

    def test_invalid_unseen_seed_cannot_reselect_training_recipe(self):
        x=fixture()
        x["icaSeedCheck"]["rankStatus"][0][4]="rejected"
        x["icaSeedCheck"]["objective"][0][4]=float("inf")
        x["icaSeedCheck"]["strata"][0][4]="other-stratum"
        result=inspect(x)
        self.assertEqual(result["trainingShadow"]["shadowCandidate"],4)
        self.assertFalse(result["shadowUnseen"]["validAllThree"])
        self.assertEqual(result["shadowUnseen"]["seeds"][1]["failureReasons"],
                         ["rank-not-feasible","stratum-drift","nonfinite-objective"])
        self.assertFalse(result["scientificPromotionAuthorized"])

    def test_bad_dimensions_are_rejected(self):
        x=fixture()
        x["icaSeedCheck"]["objective"][1].pop()
        with self.assertRaisesRegex(ValueError,"six"):
            inspect(x)


if __name__=="__main__":
    unittest.main()
