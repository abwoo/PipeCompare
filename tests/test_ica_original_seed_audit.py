"""Frozen P3 ICA seed objective audit: provenance, uncertainty and integrity."""
import copy
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location("p3ica",ROOT/"benchmarks/audit_p3_ica_original_seeds.py")
A=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(A)
DATA=ROOT/"benchmarks/results/pr51_ica_original_20261112"

class FrozenIcaSeedTests(unittest.TestCase):
    def load(self,subject):
        return json.loads((DATA/f"subject_{subject:03}.json").read_text())

    def test_known_frozen_recommendations_and_shapes(self):
        for subject,(recommended,count) in A.EXPECTED.items():
            d=A.audit(self.load(subject))
            self.assertEqual(d["recommended"],recommended)
            self.assertEqual(d["n_candidates_compared"],count)
            self.assertFalse(d["same_recommendation"])
            self.assertTrue(d["no_alternate_seed_p3_estimate_in_original_source"])
            self.assertTrue(d["no_inference_on_P3_amplitude_equivalence"])

    def test_original_subject21_candidate25_has_large_seed_sensitivity(self):
        d=A.audit(self.load(21))
        row=next(x for x in d["highest_candidate_relative_variation"] if x["id"]==25)
        self.assertGreater(row["max_relative_change_pct"],15)
        self.assertEqual(d["recommended"],[25,4,4])

    def test_sub11_selected_objectives_and_fixed_recommended_id(self):
        d=A.audit(self.load(11))
        v=[x["objective"] for x in d["selected_objective_by_seed"]]
        self.assertAlmostEqual(v[0],0.7499797775422635)
        self.assertAlmostEqual(v[1],0.7307544288840147)
        self.assertEqual(v[1],v[2])

    def test_reject_mixed_original_algorithm_or_protocol(self):
        x=self.load(11)
        x["algorithm_commit"]="altered"
        with self.assertRaisesRegex(ValueError,"identity"):
            A.audit(x)
        x=self.load(11);x["protocol_sha256"]="new_protocol"
        with self.assertRaisesRegex(ValueError,"identity"):
            A.audit(x)

    def test_reject_candidate_matrix_or_invalid_seed_selected(self):
        d=self.load(15)
        d["icaSeedCheck"]["objective"][0].pop()
        with self.assertRaisesRegex(ValueError,"matrix"):
            A.audit(d)
        d=self.load(15)
        d["icaSeedCheck"]["recommended"][0]=999
        with self.assertRaisesRegex(ValueError,"matrix"):
            A.audit(d)

    def test_emit_three_subject_report_without_inferred_p3_amplitude(self):
        with tempfile.TemporaryDirectory() as dirname:
            out=Path(dirname)/"report"
            data=A.run(DATA,out)
            self.assertEqual(len(data),3)
            self.assertEqual(json.loads((out/"seed_objective_audit.json").read_text())[0]["subject"],11)
            report=(out/"seed_objective_audit.md").read_text()
            self.assertIn("not",report.lower())
            self.assertIn("11",report)
            self.assertIn("21",report)

if __name__=="__main__":
    unittest.main()
