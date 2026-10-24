import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"benchmarks"))
from diagnose_selection import candidates, aggregate

def row(record, condition, selected, error, label="same label"):
    return dict(recording=str(record),condition=str(condition),selected=str(selected),
        steps=label,gain="1",actualError=str(error),rawError=str(error),
        reportedError="1",sme="1",eyeBias="0",shapeSpread="0",correctedEstimate="2")

class AttributionTests(unittest.TestCase):
    def test_duplicate_labels_do_not_merge_distinct_candidates(self):
        data=[row(1,1,1,3),row(1,2,1,3),row(1,1,0,1),row(1,2,0,1)]
        cs=candidates(data)["1"]
        self.assertEqual(len(cs),2)
        self.assertEqual([c["mse"] for c in cs],[9,1])
        self.assertEqual([c["selected"] for c in cs],[True,False])

    def test_incomplete_condition_blocks_are_rejected(self):
        with self.assertRaises(AssertionError):
            candidates([row(1,1,1,3),row(1,1,0,1),row(1,2,0,1)])

    def test_truth_lower_bound_and_record_balancing(self):
        def record(selected, oracle):
            return dict(selected=selected,oracle=oracle,score_minimum=selected,
                association_minimum=selected,selected_raw=selected,raw_oracle=oracle,
                reported=1,association=4)
        result=aggregate([record(9,1),record(1,1)])
        self.assertEqual(result["records"],2)
        self.assertAlmostEqual(result["selected_rmse"],5**0.5)
        self.assertEqual(result["oracle_rmse"],1)
        self.assertEqual(result["mean_selected_minus_oracle_mse"],4)
        self.assertAlmostEqual(result["association_ratio_at_unchanged_selection"],(5/4)**0.5)

if __name__ == "__main__":
    unittest.main()
