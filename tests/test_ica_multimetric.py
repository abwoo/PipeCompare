"""Multi-metric ICA screen: scientific gates, no leakage, explicit abstention."""
import sys
import unittest
from copy import deepcopy
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"benchmarks"))
from ica_multimetric import evaluate, source_gate, compare_to_original, COHORT, report


def fixture(subject=11):
    ids=COHORT[subject]
    n=len(ids)
    checks=[]
    reports=[]
    for candidate in ids:
        baseline=(candidate==ids[0])
        # Alternative has superior precision, less waveform variability,
        # and equal or better preserved IC-removal decisions.
        p3step=.015 if baseline else .006
        se=1.0 if baseline else .95
        wave=.020 if baseline else .010
        removal=.80 if baseline else 1.0
        measures=[{
            "candidateId":candidate,"seed":seed,"hasMeasurement":True,
            "gain":1.0,"correctedDifferenceEstimate":2.0+seed*p3step,
            "effectiveUnflagged":[10.,20.],"sourceUnflagged":[12.,24.]
        } for seed in range(6)]
        checks.append(measures)
        science_starts=[{
            "seed":seed,"scienceValid":True,"errors":[],
            "P3":2.0+seed*p3step,"bootstrapSE":se+seed*.0001,
            "positiveTrials":[12,24]
        } for seed in range(6)]
        science_matches=[{
            "seed":seed,"available":True,
            "minimumMatchedTopographicCorrelation":.96,
            "removedComponentMatchedJaccard":removal,
            "fullWaveformRmsDifference":wave+seed*.0001
        } for seed in range(1,6)]
        reports.append({"candidateId":candidate,"startReports":science_starts,
                        "componentMatchesVsSeed0":science_matches})
    raw={"subject":subject,"completed":True,"errors":[],
         "workflowRun":"43169562528",
         "executionCommit":"5c9388562974c8305da95cb59e7249982788d589",
         "icaSeedCheck":{"ids":ids,"recommended":[ids[0]]*6,
             "objective":[[2.0+.01*seed+.01*i for seed in range(6)] for i in range(n)],
             "rankStatus":[["feasible"]*6 for _ in ids],
             "strata":[["reref:average"]*6 for _ in ids],
             "candidateMeasurements":checks}}
    science={"subject":subject,"identityErrors":[],"scienceFailures":[],
             "unalteredFrozenRankAndEstimates":True,
             "allPredeclaredFeasibleScienceComplete":True,
             "predeclaredCandidateScience":reports}
    return raw,science


class ScientificMultiMetric(unittest.TestCase):
    def test_strict_dominating_candidate_is_shadow_only(self):
        raw,science=fixture()
        result=evaluate(raw,science)
        self.assertEqual(result["decision"]["kind"],"shadow-switch")
        self.assertEqual(result["decision"]["candidateId"],25)
        self.assertTrue(result["assessment"]["passesProtectedNonregression"])
        self.assertFalse(result["productionPromotionAuthorized"])

    def test_late_seed_changes_never_influence_training_decision(self):
        raw,science=fixture()
        before=evaluate(raw,science)["decision"]
        altered=deepcopy(science)
        for item in altered["predeclaredCandidateScience"]:
            for seed in (3,4,5):
                item["startReports"][seed]["P3"]=1000.
                item["startReports"][seed]["bootstrapSE"]=1000.
        after=evaluate(raw,altered)
        self.assertEqual(before,after["decision"])
        self.assertFalse(after["assessment"]["passesProtectedNonregression"])

    def test_less_stable_p300_is_rejected(self):
        raw,science=fixture()
        alt=science["predeclaredCandidateScience"][1]
        alt["startReports"][1]["P3"]=2.1
        result=evaluate(raw,science)
        self.assertEqual(result["decision"]["kind"],"retain-original")
        self.assertIn("p300-amplitude-drift-worse",result["alternativeTraining"][0]["failures"])

    def test_removed_ic_worsening_is_rejected(self):
        raw,science=fixture()
        science["predeclaredCandidateScience"][1]["componentMatchesVsSeed0"][1]["removedComponentMatchedJaccard"]=.5
        result=evaluate(raw,science)
        self.assertIn("removal-consistency-worse",result["alternativeTraining"][0]["failures"])

    def test_waveform_worsening_is_rejected(self):
        raw,science=fixture()
        science["predeclaredCandidateScience"][1]["componentMatchesVsSeed0"][1]["fullWaveformRmsDifference"]=.4
        result=evaluate(raw,science)
        self.assertIn("full-waveform-drift-worse",result["alternativeTraining"][0]["failures"])

    def test_under_2_percent_precision_improvement_is_rejected(self):
        raw,science=fixture()
        for z in science["predeclaredCandidateScience"][1]["startReports"]:
            z["bootstrapSE"]=.99
        result=evaluate(raw,science)
        self.assertIn("mean-bootstrap-se-gain-less-than-2pct",result["alternativeTraining"][0]["failures"])

    def test_physical_and_effective_trial_floor_never_relaxed(self):
        raw,science=fixture()
        raw["icaSeedCheck"]["candidateMeasurements"][1][2]["effectiveUnflagged"][0]=4.99
        result=evaluate(raw,science)
        self.assertIn("candidate-ineligible",result["alternativeTraining"][0]["failures"])

    def test_missing_science_fails_closed_not_zero(self):
        raw,science=fixture()
        science["predeclaredCandidateScience"][1]["componentMatchesVsSeed0"][0]["available"]=False
        self.assertIn("candidate-ineligible",evaluate(raw,science)["alternativeTraining"][0]["failures"])

    def test_bad_original_baseline_forces_abstention(self):
        raw,science=fixture()
        raw["icaSeedCheck"]["rankStatus"][0][0]="rejected"
        self.assertEqual(evaluate(raw,science)["decision"]["kind"],"abstain")

    def test_baseline_roster_and_provenance_must_match(self):
        raw,science=fixture()
        raw["icaSeedCheck"]["recommended"][0]=25
        with self.assertRaisesRegex(ValueError,"baseline ID"):
            evaluate(raw,science)
        raw,science=fixture()
        raw["workflowRun"]="mismatched"
        with self.assertRaisesRegex(ValueError,"identity"):
            evaluate(raw,science)

    def test_matlab_flat_cell_arrays_keep_candidate_seed_alignment(self):
        raw,science=fixture()
        n=len(raw["icaSeedCheck"]["ids"])
        for field in ("rankStatus","strata","candidateMeasurements"):
            rows=raw["icaSeedCheck"][field]
            raw["icaSeedCheck"][field]=[rows[i][s] for s in range(6) for i in range(n)]
        self.assertEqual(evaluate(raw,science)["decision"]["kind"],"shadow-switch")

    def test_no_tie_guessed(self):
        raw,science=fixture(15)
        # two alternative scientific options have exactly equal worst objectives.
        raw["icaSeedCheck"]["objective"][1]=raw["icaSeedCheck"]["objective"][2].copy()
        self.assertEqual(evaluate(raw,science)["decision"]["kind"],"abstain")

    def test_no_hypothesis_or_holdout_claim(self):
        raw,science=fixture()
        result=evaluate(raw,science)
        self.assertFalse(result["prospectiveValidationAuthorized"])
        self.assertEqual(result["originalPR51Verdict"],"PARTIAL")


if __name__=="__main__":
    unittest.main()
