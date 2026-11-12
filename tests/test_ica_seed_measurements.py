"""Regression: ICA seed P300 calibration must match the published pilot estimand."""
import sys
import unittest
from copy import deepcopy
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))
from validate_ica_seed_measurements import validate


def fixture():
    gain = 1.080371669
    raw = 1.7599033039905114
    first = {"seed":0,"candidateId":25,"hasMeasurement":True,
             "conditionEstimates":[2.9932702825145108,1.2333669785239993],
             "differenceEstimate":raw,"gain":gain,
             "correctedDifferenceEstimate":raw/gain,
             "sme":[1.608,0.574],"sourceUnflagged":[11,68],
             "effectiveUnflagged":[9,53]}
    second=deepcopy(first);second.update({"seed":1,"candidateId":4})
    third=deepcopy(first);third.update({"seed":2,"candidateId":4})
    return {"subject":21,"completed":True,"errors":[],
            "erp":{"available":True,"differenceAmplitude":raw/gain},
            "icaSeedCheck":{"recommended":[25,4,4],"same":False,
                            "measurements":[first,second,third]}}


class ValidateSeedResults(unittest.TestCase):
    def test_exact_corrected_estimand(self):
        r=validate(fixture())
        self.assertEqual(r["subject"],21)
        self.assertEqual(r["candidateIds"],[25,4,4])
        self.assertFalse(r["stableCandidateId"])
        self.assertAlmostEqual(r["correctedDeltaFromFirst"][1],0)
        self.assertTrue(r["noScientificGatesChanged"])

    def test_raw_is_not_corrected(self):
        x=fixture()
        x["erp"]["differenceAmplitude"]=x["icaSeedCheck"]["measurements"][0]["differenceEstimate"]
        with self.assertRaisesRegex(ValueError,"gain-corrected"):
            validate(x)

    def test_calibration_inconsistent(self):
        x=fixture()
        x["icaSeedCheck"]["measurements"][1]["gain"]=0.25
        with self.assertRaisesRegex(ValueError,"inconsistent"):
            validate(x)

    def test_invalid_effective_source_trial_bound(self):
        x=fixture()
        x["icaSeedCheck"]["measurements"][2]["effectiveUnflagged"]=[12,53]
        with self.assertRaisesRegex(ValueError,"exceed"):
            validate(x)

    def test_wrong_seed_candidate_id(self):
        x=fixture()
        x["icaSeedCheck"]["measurements"][1]["candidateId"]=25
        with self.assertRaisesRegex(ValueError,"does not match ranking"):
            validate(x)

    def test_missing_different_seed_is_not_imputed(self):
        x=fixture()
        x["icaSeedCheck"]["measurements"][1]={"seed":1,"candidateId":4,
                                                "hasMeasurement":False}
        z=validate(x)
        self.assertFalse(z["bySeed"][1]["available"])
        self.assertIsNone(z["correctedDeltaFromFirst"][1])

    def test_technical_failure_not_scientific_result(self):
        x=fixture()
        x["errors"]=[{"stage":"standard","message":"failed"}]
        with self.assertRaisesRegex(ValueError,"uncompleted"):
            validate(x)


if __name__=="__main__":
    unittest.main()
