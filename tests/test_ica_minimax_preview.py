"""Worst-across-three-seed candidate score is NOT scientific pass evidence."""
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))
from preview_ica_minimax import preview


def sample():
    return {"subject":21,"ids":[4,25],
            "objective_by_candidate_and_seed":[
                [1.12,1.127,1.133],
                [1.121,1.350,1.354]
            ],"seed_pipelines":[25,4,4]}


class MinimaxPreviewTests(unittest.TestCase):
    def test_worst_case_proxy_prefers_stable_objective_not_seed0_id(self):
        x=preview(sample())
        self.assertEqual(x["proposal"],4)
        self.assertEqual(x["original_first_seed_candidate"],25)
        self.assertGreater(x["proxy_worst_case_reduction"],0.2)
        self.assertFalse(x["scientific_feasibility_checked"])
        self.assertFalse(x["eligible_for_real_algorithm_promotion"])

    def test_absent_seed_objective_excludes_candidate(self):
        x=sample()
        x["objective_by_candidate_and_seed"][0][2]=None
        out=preview(x)
        self.assertEqual(out["proposal"],25)
        self.assertEqual(out["usable_finite_rows"],1)
        self.assertFalse(out["eligible_for_real_algorithm_promotion"])

    def test_no_finite_rows_abstains(self):
        x=sample()
        x["objective_by_candidate_and_seed"]=[[None,1,1],[1,None,1]]
        out=preview(x)
        self.assertIsNone(out["proposal"])
        self.assertEqual(out["reason"],"no-all-seed-finite-objective-candidate")

    def test_exact_ties_not_arbitrarily_resolved(self):
        x=sample()
        x["objective_by_candidate_and_seed"]=[[1,2,3],[1,2,3]]
        out=preview(x)
        self.assertIsNone(out["proposal"])
        self.assertEqual(out["tie_ids"],[4,25])

    def test_misaligned_candidate_ids_fail(self):
        x=sample()
        x["ids"]=[4,4]
        with self.assertRaises(ValueError):
            preview(x)


if __name__=="__main__":
    unittest.main()
