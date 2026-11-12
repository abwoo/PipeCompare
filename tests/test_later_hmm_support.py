"""No accidental eligibility promotion from longer raw recordings."""
import sys
import unittest
from copy import deepcopy
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"benchmarks"))

from validate_later_hmm_support import check

def fixture():
    def cond(name,n,k):
        return {"condition":name,"referenceN":n,"unflaggedN":k,"flaggedN":n-k,
                "unflaggedIds":list(range(1,k+1)),
                "flaggedIds":list(range(k+1,n+1))}
    return {"subject":4,"originalSeconds":300,"fullSeconds":544,
        "sameHMMScienceAndTrialRules":True,
        "original20SubjectVerdict":"PARTIAL","productionPromotionAuthorized":False,
        "baseline":[cond("target",18,4),cond("standard",95,3)],
        "extended":[cond("target",30,5),cond("standard",159,6)],
        "fullWindowUnflaggedFloorMet":True}

class HMMGate(unittest.TestCase):
    def test_complete_later_data_can_meet_reference_floor_only(self):
        z=check(fixture())
        self.assertTrue(z["fullHMMFloorMet"])
        self.assertFalse(z["pairedMeasurementNowValid"])
        self.assertEqual(z["originalUnflagged"],[4,3])

    def test_full_window_negative_scientific_result_is_preserved(self):
        z=fixture()
        z["extended"][0]["unflaggedN"]=4
        z["extended"][0]["flaggedN"]=26
        z["extended"][0]["unflaggedIds"]=list(range(1,5))
        z["extended"][0]["flaggedIds"]=list(range(5,31))
        z["fullWindowUnflaggedFloorMet"]=False
        self.assertFalse(check(z)["fullHMMFloorMet"])

    def test_original_source_support_cannot_be_rewritten(self):
        z=fixture()
        z["baseline"][0]["unflaggedN"]=5
        with self.assertRaises(ValueError):
            check(z)

    def test_invented_extra_support_rejected(self):
        z=fixture()
        z["extended"][1]["unflaggedN"]=7
        with self.assertRaises(ValueError):
            check(z)

    def test_modified_science_rules_rejected(self):
        z=fixture();z["sameHMMScienceAndTrialRules"]=False
        with self.assertRaises(ValueError):
            check(z)

if __name__=="__main__":
    unittest.main()
