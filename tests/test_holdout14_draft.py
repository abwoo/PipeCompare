"""The prospective ERP CORE holdout MUST remain entirely unseen and unexecuted."""
import json
import sys
import unittest
from copy import deepcopy
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"benchmarks"))
from validate_holdout14_draft import validate


def fixture():
    p=json.loads((ROOT/"benchmarks"/"erpcore_holdout14_preanalysis_20261112.json").read_text())
    f=json.loads((ROOT/"benchmarks"/"erpcore_protocol.json").read_text())
    return p,f


class HoldoutRules(unittest.TestCase):
    def test_nonoverlap_and_original_exclusions(self):
        p,f=fixture()
        r=validate(p,f)
        self.assertEqual(r["development_count"],20)
        self.assertEqual(len(r["holdout_subjects"]),14)
        self.assertFalse(r["outcomes_reviewed"])

    def test_does_not_treat_preanalysis_as_frozen_eeg(self):
        p,f=fixture();p["status"]="FROZEN_AND_VALIDATED"
        with self.assertRaises(ValueError): validate(p,f)

    def test_overlap_invalidates_cohort(self):
        p,f=fixture();p["previously_unseen_14"][0]=11
        with self.assertRaisesRegex(ValueError,"contamination"):validate(p,f)

    def test_original_public_exclusions_cannot_be_changed(self):
        p,f=fixture();p["excluded_subjects"]=[6,9]
        with self.assertRaisesRegex(ValueError,"exclusions"):validate(p,f)

    def test_cannot_drop_poor_participants(self):
        p,f=fixture();p["frozen_after_algorithm_selection_gates"]["matched_pair_subjects"]=13
        with self.assertRaisesRegex(ValueError,"All 14"):validate(p,f)

    def test_cannot_weaken_five_unflagged_trials(self):
        p,f=fixture();p["frozen_after_algorithm_selection_gates"]["minimum_effective_unflagged_each_arm_each_condition"]=4
        with self.assertRaisesRegex(ValueError,"support"):validate(p,f)


if __name__=="__main__":
    unittest.main()
