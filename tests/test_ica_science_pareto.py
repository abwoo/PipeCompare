"""Training-only Pareto comparison must abstain on scientific tradeoffs."""
import sys
import unittest
from copy import deepcopy
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"benchmarks"))
sys.path.insert(0,str(Path(__file__).resolve().parent))
from test_ica_unseen_seeds import fixture as original_fixture
from test_ica_science import synth
from ica_science_pareto import dominates, train_selection, evaluate


def make_study():
    x=original_fixture()
    x["subject"]=21
    ids=x["icaSeedCheck"]["ids"]
    for i,ident in enumerate(ids):
        for seed in range(6):
            m=synth()
            m["candidateId"]=ident
            m["seed"]=seed
            m["gain"]=1.0
            m["hasMeasurement"]=True
            m["sourceUnflagged"]=[6,6]
            m["effectiveUnflagged"]=[6,6]
            amplitude=m["correctedDifferenceEstimate"]
            for j,cond in enumerate(m["science"]["conditions"]):
                center=amplitude if j==0 else 0.
                spread=0.1 if ident==4 else 1.0
                cond["referenceIds"]=[10*j+k for k in range(1,7)]
                cond["retainedIds"]=cond["referenceIds"][:]
                cond["referenceFlagged"]=[False]*6
                cond["trialScores"]=[center+spread*v for v in (-3,-2,-1,1,2,3)]
                cond["trialWeights"]=[1.]*6
                cond["usedCount"]=6
                cond["effectiveN"]=6
                cond["windowMean"]=center
                cond["originalEstimate"]=center
            x["icaSeedCheck"]["candidateMeasurements"][i][seed]=m
    return x


class ParetoScience(unittest.TestCase):
    def test_dominance_four_axes(self):
        low={"worstBootstrapSE":1.,"maxFullWaveRms":0.1,
             "minusMinTopologyCorrelation":-.95,"minusMinRemovalJaccard":-1.}
        high=dict(low,worstBootstrapSE=2.)
        self.assertTrue(dominates(low,high))
        self.assertFalse(dominates(high,low))

    def test_training_selection_prefers_lower_error_when_other_axes_tie(self):
        x=make_study()
        selected=train_selection(x)
        self.assertEqual(selected["candidateId"],4)
        self.assertEqual(selected["frontierCount"],1)

    def test_unseen_scores_cannot_change_training_choice(self):
        x=make_study()
        before=train_selection(x)["candidateId"]
        for row in x["icaSeedCheck"]["candidateMeasurements"]:
            for s in (3,4,5):
                row[s]["science"]["conditions"][0]["trialScores"]=[10,10,10,10,10,10]
                row[s]["correctedDifferenceEstimate"]=100.
        self.assertEqual(before,train_selection(x)["candidateId"])

    def test_real_tradeoff_returns_abstention(self):
        x=make_study()
        # The lower-error candidate changes removed IC identity at one
        # training start, while the higher-error candidate stays consistent.
        x["icaSeedCheck"]["candidateMeasurements"][0][1]["science"]["decompositions"][0]["removed"]=[2]
        selected=train_selection(x)
        self.assertIsNone(selected["candidateId"])
        self.assertEqual(selected["reason"],"multiobjective-tradeoff-abstain")
        self.assertEqual(selected["frontierCount"],2)

    def test_missing_science_cannot_be_inferred(self):
        x=make_study()
        x["icaSeedCheck"]["candidateMeasurements"][0][2]["science"]["status"]="failed"
        selected=train_selection(x)
        self.assertEqual(selected["candidateId"],25)
        self.assertEqual(selected["scienceCompleteCandidates"],1)

    def test_never_authorizes_production_or_new_holdout(self):
        x=make_study()
        z=evaluate(x)
        self.assertFalse(z["productionPromotionAuthorized"])
        self.assertFalse(z["open14SubjectConfirmationAuthorized"])
        self.assertEqual(z["training"]["selectionSeeds"],[0,1,2])


if __name__=="__main__":
    unittest.main()
