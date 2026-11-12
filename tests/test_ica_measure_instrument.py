"""Ensure the stage-3 ICA instrument edits only a temporary checked-out Executor."""
import importlib.util
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location("ica_instrument",ROOT/".github/scripts/pr51-instrument-ica-measure.py")
MOD=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MOD)

class IcaInstrumentation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.original=(ROOT/"+pipecompare/+run/Executor.m").read_text()

    def test_patch_frozen_source_unique_and_nonmutating(self):
        patched=MOD.patch(self.original)
        self.assertEqual(self.original.count("out.seedMeasures=seedMeasures;"),0)
        self.assertEqual(patched.count("out.seedMeasures=seedMeasures;"),1)
        self.assertEqual(patched.count("seedMeasures{seed+1}"),1)
        self.assertIn("correctedDifferenceAmplitude",patched)
        self.assertIn("differenceSEAnalytic",patched)
        self.assertIn("nICRemoved",patched)
        self.assertIn("diagnosticError",patched)
        self.assertEqual((ROOT/"+pipecompare/+run/Executor.m").read_text(),self.original)

    def test_patch_only_expected_seed_observation_anchors(self):
        patched=MOD.patch(self.original)
        self.assertEqual(patched.count(MOD.START),1)
        self.assertEqual(patched.count(MOD.PERSEED),1)
        self.assertEqual(patched.count(MOD.ENDING),1)
        self.assertEqual(patched.count("function c = rerunWithIca"),1)

    def test_double_instrumentation_is_rejected(self):
        patched=MOD.patch(self.original)
        with self.assertRaisesRegex(ValueError,"unique"):
            MOD.patch(patched)

    def test_expected_anchor_absent_is_rejected(self):
        with self.assertRaisesRegex(ValueError,"unique"):
            MOD.patch(self.original.replace(MOD.PERSEED,"changed-score-function"))

if __name__ == "__main__":
    unittest.main()
