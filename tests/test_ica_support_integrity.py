"""Physical trial count is a hard upper bound on weighted effective N."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"benchmarks"))
from ica_support_integrity import check_measurement, audit


def measurement(candidate=4,seed=0,raw=(10,20),eff=(7,12)):
    return {"candidateId":candidate,"seed":seed,"sourceUnflagged":list(raw),
            "effectiveUnflagged":list(eff),"hasMeasurement":True}


def study():
    ids=[4,25]
    cells=[[measurement(candidate,seed) for seed in range(6)] for candidate in ids]
    return {"subject":11,"completed":True,"errors":[],
            "icaSeedCheck":{"ids":ids,"candidateMeasurements":cells}}


class SupportIntegrity(unittest.TestCase):
    def test_valid_nonuniform_weights(self):
        self.assertEqual(check_measurement(measurement(),4,0),[])
    def test_exact_physical_upper_bound_allowed(self):
        self.assertEqual(check_measurement(measurement(raw=(5,5),eff=(5,5)),4,0),[])
    def test_effective_more_than_physical_is_not_accepted(self):
        bad=check_measurement(measurement(raw=(4,10),eff=(5,6)),4,0)
        self.assertIn("target-effective-exceeds-physical",bad)
        self.assertIn("target-eligible-with-insufficient-physical-trials",bad)
    def test_negative_or_nan_effective_n_is_invalid(self):
        bad=check_measurement(measurement(eff=(-1,float("nan"))),4,0)
        self.assertIn("target-invalid-effective-n",bad)
        self.assertIn("standard-invalid-effective-n",bad)
    def test_no_missing_support_as_zero(self):
        m=measurement()
        m.pop("sourceUnflagged")
        self.assertEqual(check_measurement(m,4,0),["support-pair-missing"])
    def test_seed_candidate_alignment(self):
        self.assertEqual(check_measurement(measurement(),4,2),["candidate-seed-misaligned"])
    def test_matlab_flattened_cell_grid(self):
        x=study()
        cells=x["icaSeedCheck"]["candidateMeasurements"]
        x["icaSeedCheck"]["candidateMeasurements"]=[cells[i][s] for s in range(6) for i in range(2)]
        self.assertTrue(audit(x)["supportIntegrityPassed"])
        self.assertEqual(audit(x)["examinedCandidateStarts"],12)
    def test_invalid_one_seed_preserves_negative_evidence(self):
        x=study()
        x["icaSeedCheck"]["candidateMeasurements"][0][5]["effectiveUnflagged"]=[11,12]
        z=audit(x)
        self.assertFalse(z["supportIntegrityPassed"])
        self.assertEqual(z["findings"][0]["candidateId"],4)
        self.assertEqual(z["findings"][0]["seed"],5)
        self.assertFalse(z["productionPromotionAuthorized"])
    def test_absent_physical_count_not_inferred(self):
        x=study()
        x["icaSeedCheck"]["candidateMeasurements"][1][3]["sourceUnflagged"]=[]
        self.assertFalse(audit(x)["supportIntegrityPassed"])
    def test_incomplete_scientific_record_rejected(self):
        x=study()
        x["completed"]=False
        with self.assertRaisesRegex(ValueError,"Incomplete"):
            audit(x)


if __name__=="__main__":
    unittest.main()
