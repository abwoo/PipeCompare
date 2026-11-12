"""Tests for stage-2 read-only P300 support audit; no EEG processing."""
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SNAP=ROOT/"benchmarks"/"results"/"erpcore_rootcause_20261112"
SPEC=importlib.util.spec_from_file_location("audit_p3_effective_support",ROOT/"benchmarks"/"audit_p3_effective_support.py")
AUDIT=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)

class SupportAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.subjects=AUDIT.load_snapshots(SNAP)
        cls.output=AUDIT.make_audit(cls.subjects)

    def test_exhaustive_identity_and_arm_coverage(self):
        self.assertEqual(len(self.subjects),20)
        self.assertEqual(tuple(sorted(self.subjects)),AUDIT.EXPECTED)
        self.assertEqual(len(self.output["arm_condition_rows"]),60)
        self.assertEqual(self.output["protocol_sha256"],AUDIT.PROTOCOL)

    def test_original_frozen_support_floor_for_priority_participants(self):
        def two(sid,arm):
            return [v["effective_unflagged"] for v in
                    AUDIT.audit_arm(sid,arm,self.subjects[sid])["conditions"]]
        self.assertEqual(two(4,"standard"),[4,3])
        self.assertEqual(two(4,"minimal"),[1,1])
        self.assertEqual(two(8,"minimal"),[4,18])
        self.assertEqual(two(8,"common"),[4,18])
        self.assertEqual(two(17,"minimal"),[4,27])
        self.assertEqual(two(18,"minimal"),[0,9])
        self.assertEqual(two(23,"common"),[0,11])
        self.assertFalse(AUDIT.audit_arm(4,"standard",self.subjects[4])["all_condition_support_meets_floor"])

    def test_no_recommendation_does_not_mean_measured_zero_standard_support(self):
        row=AUDIT.audit_arm(18,"standard",self.subjects[18])
        self.assertFalse(row["readout_available"])
        self.assertFalse(row["support_available"])
        self.assertEqual(row["conditions"],[])
        self.assertEqual(row["unassessable_without_trial_data"][0],
                         "which exact urevents were blink-flagged")

    def test_keep_fractional_weighted_effective_support(self):
        s=self.subjects[15]
        values=[v["effective_unflagged"] for v in AUDIT.audit_arm(15,"standard",s)["conditions"]]
        self.assertTrue(any(x!=int(x) for x in values))
        self.assertAlmostEqual(values[0],21.858364864849502)

    def test_provenance_mismatch_and_duplicate_subject_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            blobs=[json.loads(p.read_text()) for p in sorted(SNAP.glob("source_*.json"))]
            blobs[0]["subjects"][0]["protocol_sha256"]="wrong"
            for i,blob in enumerate(blobs):
                (root/f"source_{i:02}.json").write_text(json.dumps(blob))
            with self.assertRaisesRegex(ValueError,"mixed algorithm/protocol"):
                AUDIT.load_snapshots(root)

    def test_three_output_files_keep_completeness(self):
        with tempfile.TemporaryDirectory() as d:
            AUDIT.write_artifacts(self.output,Path(d))
            self.assertEqual(sorted(x.name for x in Path(d).iterdir()),
                ["effective_support.csv","effective_support.json","effective_support.md"])
            obj=json.loads((Path(d)/"effective_support.json").read_text())
            self.assertEqual(obj["unique_subjects"],20)
            self.assertEqual(len(obj["arm_condition_rows"]),60)
            self.assertIn("zero retained target trials",(Path(d)/"effective_support.md").read_text())

if __name__=="__main__":
    unittest.main()
