"""Pinned P3 VEOG proxy tests; diagnostic, not ground-truth blink accuracy."""
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
import numpy as np
from scipy.io import savemat

ROOT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location("p3veog",ROOT/"benchmarks/probe_p3_veog.py")
M=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)
SAVED=ROOT/"benchmarks/results/erpcore_support_trial_20261112"

class VeogProbeTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.folder=Path(self.temp.name)
        self.source={"subject":4,"protocol_sha256":M.PROTOCOL,"sampling_rate":256,
                     "n_valid_targets":1,"n_valid_standards":1,
                     "source_set_sha256":"original-set","source_fdt_sha256":"original-fdt",
                     "source_events_sha256":"original-events"}
        self.veog=np.zeros(2300)
        # Event indices 1 and 2 have latencies 500 and 1500 (MATLAB 1-based).
        self.veog[500-1+20]=75
        self.save_input()

    def save_input(self):
        self.input=self.folder/"input.mat"
        savemat(self.input,{"fs":256.,"veog":self.veog,
            "eventTypes":np.array(["11","31"],dtype=object),
            "eventLatency":np.array([500,1500])})
        self.source_path=self.folder/"source.json"
        self.source_path.write_text(json.dumps(self.source))

    def test_baseline_epoch_boundary_and_proxy_class(self):
        d=M.analyze(self.input,self.source_path,4)
        self.assertEqual(d["n_events"],2)
        self.assertTrue(d["proxy_is_not_ground_truth"])
        self.assertEqual(d["conditions"][0]["veog_proxy_above50"],1)
        self.assertEqual(d["conditions"][1]["veog_proxy_above50"],0)
        self.assertAlmostEqual(d["events"][0]["veog_peak_abs_baseline_corrected_uv"],75)
        self.assertEqual(d["events"][0]["event_index_1based"],1)
        self.assertIsNone(d["events"][0]["selected_standard_hmm_flag"])

    def test_frozen_selected_support_join_on_event_and_condition(self):
        packet={"subject":4,"source_run_id":M.FROZEN_RUN,"protocol_sha256":M.PROTOCOL,
             "source_job_id":1234,
             "source_files_sha256":{"set":"original-set","fdt":"original-fdt","events":"original-events"},
             "selectedSupport":[{"condition":"11","referenceN":1,"referenceIds":[1],
                                 "flagged":[True],"weights":[0.],"scores":[2.5]},
                                 {"condition":"31","referenceN":1,"referenceIds":[2],
                                  "flagged":[False],"weights":[1.],"scores":[1.5]}]}
        p=self.folder/"selected.json";p.write_text(json.dumps(packet))
        d=M.analyze(self.input,self.source_path,4,p)
        self.assertEqual(d["original_selected_support_job"],1234)
        self.assertEqual(d["conditions"][0]["hmm_flagged_and_veog_above50"],1)
        self.assertEqual(d["conditions"][1]["hmm_unflagged_and_veog_below_equal50"],1)
        self.assertEqual(d["events"][0]["selected_standard_weight"],0)
        self.assertEqual(d["events"][1]["selected_standard_weight"],1)

    def test_source_and_event_count_mismatch_rejected(self):
        self.source["n_valid_targets"]=0
        self.save_input()
        with self.assertRaisesRegex(ValueError,"source event counts"):
            M.analyze(self.input,self.source_path,4)

    def test_invalid_selected_support_source_sha_rejected(self):
        packet=json.loads((SAVED/"selected_004.json").read_text())
        source={"subject":4,"source_set_sha256":"invalid",
            "source_fdt_sha256":packet["source_files_sha256"]["fdt"],
            "source_events_sha256":packet["source_files_sha256"]["events"]}
        with self.assertRaisesRegex(ValueError,"source hash mismatch"):
            M.read_selected(SAVED/"selected_004.json",source,4)

    def test_real_original_4_8_event_vectors_are_self_consistent(self):
        for subject in (4,8):
            packet=json.loads((SAVED/f"selected_{subject:03}.json").read_text())
            source={f"source_{a}_sha256":packet["source_files_sha256"][a]
                    for a in ("set","fdt","events")}
            lookup,job=M.read_selected(SAVED/f"selected_{subject:03}.json",source,subject)
            self.assertGreater(job,0)
            self.assertEqual(len(lookup),113 if subject==4 else 117)

    def test_output_csv_and_json_keep_event_provenance(self):
        d=M.analyze(self.input,self.source_path,4)
        M.save(d,self.folder/"evidence")
        result=json.loads((self.folder/"evidence/subject_004.json").read_text())
        self.assertEqual(len(result["events"]),2)
        self.assertEqual(result["conditions"][0]["condition"],"11")
        self.assertTrue((self.folder/"evidence/subject_004.csv").exists())

if __name__=="__main__":
    unittest.main()
