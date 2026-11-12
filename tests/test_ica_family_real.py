"""Scientific gate and provenance regression for family-restricted ICA minimax."""
import sys
import unittest
from copy import deepcopy
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"tests"))
sys.path.insert(0,str(ROOT/"benchmarks"))
from test_ica_robust_real import sample
from validate_ica_family_real import validate


def fixture():
    x=sample()
    a=x["icaSeedCheck"]
    for row in a["candidateSeedAudit"]:
        row["estimationFamily"]=["weighted"]*3 if row["candidateId"]==25 else ["unweighted"]*3
    selected=deepcopy(a["robust"]["seedMeasurements"])
    for z in selected:
        z["candidateId"]=25
        z["estimationFamily"]="weighted"
    a["familyRobust"]={
        "status":"selected","candidateId":25,
        "baselineId":25,"baselineFamily":"weighted",
        "estimatorFamilyLocked":True,"originalRecommendationUnchanged":True,
        "scientific":{"worstObjective":1.354},
        "seedMeasurements":selected}
    return x


class FamilyResearchGateTests(unittest.TestCase):
    def test_only_original_weighting_family_may_be_used(self):
        r=validate(fixture())
        self.assertEqual(r["familyCandidate"],25)
        self.assertEqual(r["baselineFamily"],"weighted")
        self.assertEqual(r["oldCommonFeasibleCandidate"],4)
        self.assertEqual(r["allSeedScientificAndFamilyFeasibleCount"],1)
        self.assertFalse(r["isNewFrozenExperiment"])

    def test_wrong_original_family_fails(self):
        x=fixture()
        x["icaSeedCheck"]["familyRobust"]["baselineFamily"]="unweighted"
        with self.assertRaisesRegex(ValueError,"baseline estimator-family"):
            validate(x)

    def test_cross_family_shortcut_fails(self):
        x=fixture()
        x["icaSeedCheck"]["familyRobust"]["candidateId"]=4
        with self.assertRaisesRegex(ValueError,"family-constrained minimax"):
            validate(x)

    def test_weighting_change_under_seed_fails(self):
        x=fixture()
        x["icaSeedCheck"]["familyRobust"]["seedMeasurements"][1]["estimationFamily"]="unweighted"
        with self.assertRaisesRegex(ValueError,"changed scientific estimator"):
            validate(x)

    def test_scientific_gate_fails_closed(self):
        x=fixture()
        x["icaSeedCheck"]["candidateSeedAudit"][1]["status"][2]="rejected"
        with self.assertRaisesRegex(ValueError,"family-constrained minimax"):
            validate(x)

    def test_insufficient_HMM_support_fails_closed(self):
        x=fixture()
        x["icaSeedCheck"]["familyRobust"]["seedMeasurements"][0]["effectiveUnflagged"][0]=4
        with self.assertRaisesRegex(ValueError,"trial floor"):
            validate(x)

    def test_missing_family_alignment_fails(self):
        x=fixture()
        x["icaSeedCheck"]["candidateSeedAudit"][1]["estimationFamily"]=["weighted"]
        with self.assertRaisesRegex(ValueError,"not aligned"):
            validate(x)


if __name__=="__main__":
    unittest.main()
