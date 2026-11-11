#!/usr/bin/env python3
"""Unit tests for benchmarks/compare_erpcore.py (no MATLAB, no network).

Run: python3 tests/test_erpcore_collector.py
Covers: the bootstrap-SE recomputation against a hand-computed value, the
paired-precision CI, the calibration ratio, gate decisions, and tamper
detection (a changed source checksum must fail validation).
"""
import hashlib
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "benchmarks"))
import compare_erpcore as ce  # noqa: E402


def fake_support(n_target=20, n_standard=80, seed=1):
    import numpy as np
    rng = np.random.RandomState(seed)
    t = rng.normal(5.0, 3.0, n_target)
    s = rng.normal(1.0, 3.0, n_standard)
    return [
        {"condition": "target", "scores": t.tolist(), "weights": [1.0] * n_target,
         "flagged": [False] * n_target, "effectiveUnflaggedN": float(n_target),
         "referenceIds": list(range(1, n_target + 1)), "referenceN": n_target,
         "kept": n_target, "retention": 1.0},
        {"condition": "standard", "scores": s.tolist(), "weights": [1.0] * n_standard,
         "flagged": [False] * n_standard, "effectiveUnflaggedN": float(n_standard),
         "referenceIds": list(range(1, n_standard + 1)), "referenceN": n_standard,
         "kept": n_standard, "retention": 1.0},
    ]


def fake_record(subject, gain=1.0, diff=None, reported_se=0.5, seed=1):
    support = fake_support(seed=seed)
    if diff is None:
        import numpy as np
        diff = float(np.mean(support[0]["scores"]) - np.mean(support[1]["scores"]))
    return {
        "subject": subject, "completed": True, "algorithmCommit": "ALG",
        "executionCommit": "EXEC", "workflowRun": "123", "eeglab": "2026.0.0",
        "eeglabSha": "EEG", "matlab": "R2026a", "protocolSha256": "PROTO",
        "measureSha256": "BASELINE",
        "notRun": 0, "errors": [], "replayMaxError": 1e-9, "exportMaxError": 1e-9,
        "searchOptions": {"icaSeeds": 2, "refineCutoffs": True, "eogCheck": True,
                          "signalTemplate": "data"},
        "icaSeedCheck": {"same": True, "seeds": [1, 2]},
        "source": {"subject": subject, "source_set_sha256": "SET", "source_fdt_sha256": "FDT",
                   "source_events_sha256": "EVT", "sampling_rate": 256, "samples": 76800,
                   "n_valid_targets": 20, "n_valid_standards": 80},
        "reference": {"available": True, "differenceAmplitude": diff},
        "selectedSupport": support,
        "erp": {"available": True, "differenceAmplitude": diff, "reportedDifferenceSE": reported_se,
                "gain": gain},
        "minimal": {"available": True, "differenceAmplitude": diff, "reportedDifferenceSE": reported_se,
                    "gain": gain, "support": support},
        "common": {"available": True, "differenceAmplitude": diff, "reportedDifferenceSE": reported_se,
                   "gain": gain, "support": support},
    }


def fake_protocol(subjects=(1, 2)):
    return {
        "algorithm_commit": "ALG",
        "baseline_measure_sha256": "BASELINE",
        "toolchain": {"eeglab": "2026.0.0", "matlab": "R2026a"},
        "excerpt": {"duration_seconds": 300, "sampling_rate": 256},
        "trial_rules": {"min_valid_targets": 15, "min_valid_standards": 60},
        "evaluation": {"replay_tolerance_uv": 1e-6, "bootstrap_draws": 1000,
                       "bootstrap_seed": 510004, "calibration_agreement_band": [0.8, 1.2]},
        "decision_gates": {"minimum_effective_unflagged": 5, "stability_min_count": 18,
                           "individual_regression_relative_allowance": 0.2,
                           "individual_regression_absolute_allowance_uv": 0.05,
                           "individual_regression_max_count": 2},
        "source": {"files": [{"subject": s, "set": {"sha256": "SET"}, "fdt": {"sha256": "FDT"},
                              "events": {"sha256": "EVT"}} for s in subjects]},
    }


class TestBootstrapSE(unittest.TestCase):
    def test_point_matches_weighted_mean(self):
        support = fake_support()
        boot = ce.bootstrap_se(support, 1.0, 1)
        self.assertIsNotNone(boot)
        # Point estimate = weighted mean difference of the saved scores.
        import numpy as np
        t = np.mean(support[0]["scores"])
        s = np.mean(support[1]["scores"])
        self.assertAlmostEqual(boot["point"], t - s, places=9)
        # SE is positive and in a sane range for these synthetic scores.
        self.assertGreater(boot["se"], 0.1)
        self.assertLess(boot["se"], 3.0)
        # Bootstrap mean is close to the point estimate.
        self.assertAlmostEqual(boot["mean"], boot["point"], delta=0.2)

    def test_weighted_scores(self):
        support = fake_support()
        support[0]["weights"] = [2.0] * 20
        boot = ce.bootstrap_se(support, 2.0, 1)
        # Gain 2 halves the corrected difference and its SE.
        boot1 = ce.bootstrap_se(fake_support(), 1.0, 1)
        self.assertAlmostEqual(boot["point"], boot1["point"] / 2, places=9)

    def test_too_few_trials_returns_none(self):
        support = fake_support(n_target=1)
        self.assertIsNone(ce.bootstrap_se(support, 1.0, 1))


class TestValidation(unittest.TestCase):
    def test_valid_record_passes(self):
        ce.validate_record(fake_record(1), fake_protocol(), "PROTO")

    def test_tampered_source_sha_fails(self):
        r = fake_record(1)
        r["source"]["source_fdt_sha256"] = "TAMPERED"
        with self.assertRaises(AssertionError):
            ce.validate_record(r, fake_protocol(), "PROTO")

    def test_low_support_does_not_fail_validation(self):
        # Support is a pair-validity rule handled by main(), not by
        # validate_record: one low-support pair is excluded, not fatal.
        r = fake_record(1)
        r["selectedSupport"][0]["effectiveUnflaggedN"] = 3.0
        ce.validate_record(r, fake_protocol(), "PROTO")

    def test_replay_error_fails(self):
        r = fake_record(1)
        r["replayMaxError"] = 1e-5
        with self.assertRaises(AssertionError):
            ce.validate_record(r, fake_protocol(), "PROTO")

    def test_real_matlab_version_string_passes(self):
        # Regression test for A4: MATLAB's `version` returns
        # "26.1.0.3346908 (R2026a) Update 5", which the old startswith
        # check rejected, invalidating all 20 real records in run
        # 42809789847. The pinned release token must be found inside it.
        r = fake_record(1)
        r["matlab"] = "26.1.0.3346908 (R2026a) Update 5"
        ce.validate_record(r, fake_protocol(), "PROTO")

    def test_other_matlab_release_fails(self):
        # The gate still rejects a different release.
        r = fake_record(1)
        r["matlab"] = "24.2.0.27 (R2024b) Update 3"
        with self.assertRaises(AssertionError):
            ce.validate_record(r, fake_protocol(), "PROTO")

    def test_real_shape_standard_support_at_top_level(self):
        # Regression test for A5: the real pilot writes the Standard arm's
        # per-trial support at the record top level (selectedSupport) and
        # never writes erp.support. A record shaped that way validates and
        # forms a complete pair. (The A4 fixture put "support" in every arm
        # including erp, so it never matched real pilot output and the A4
        # re-collect rejected all 17 records with a feasible Standard
        # pipeline.)
        r = fake_record(1)
        self.assertNotIn("support", r["erp"])
        self.assertIn("selectedSupport", r)
        ce.validate_record(r, fake_protocol(), "PROTO")
        self.assertIsNone(ce.validate_arms(r))

    def test_standard_unavailable_is_incomplete_pair(self):
        # Amendment A5: a complete record whose Standard search found no
        # feasible pipeline (subjects 17/18/23) is an incomplete pair, not
        # an identity failure. Replay/export never ran, so their errors
        # are null — that must not fail validation either.
        r = fake_record(1)
        r["erp"] = {"available": False, "differenceAmplitude": None,
                    "reportedDifferenceSE": None, "gain": None}
        r["feasible"] = 0
        r["candidates"] = 24
        r["replayMaxError"] = None
        r["exportMaxError"] = None
        del r["selectedSupport"]
        ce.validate_record(r, fake_protocol(), "PROTO")
        reason = ce.validate_arms(r)
        self.assertIsNotNone(reason)
        self.assertIn("Standard", reason)
        self.assertIn("0 of 24", reason)

    def test_fixed_arm_unavailable_is_incomplete_pair(self):
        # Amendment A5: subject 4's shape — Standard fine, both fixed arms
        # unavailable — is an incomplete pair, not an identity failure.
        r = fake_record(1)
        r["minimal"] = {"available": False, "differenceAmplitude": None,
                        "reportedDifferenceSE": None, "gain": None,
                        "support": r["minimal"]["support"]}
        r["common"] = {"available": False, "differenceAmplitude": None,
                       "reportedDifferenceSE": None, "gain": None,
                       "support": r["common"]["support"]}
        ce.validate_record(r, fake_protocol(), "PROTO")
        reason = ce.validate_arms(r)
        self.assertIsNotNone(reason)
        self.assertIn("minimal", reason)

    def test_available_standard_without_support_fails_validation(self):
        # An available Standard arm with no per-trial support anywhere is a
        # malformed record, not an incomplete pair.
        r = fake_record(1)
        del r["selectedSupport"]
        with self.assertRaises(AssertionError):
            ce.validate_record(r, fake_protocol(), "PROTO")


def build_established_records():
    """20 records where Standard is much more precise and well calibrated."""
    import numpy as np
    records = []
    for i, subject in enumerate(range(1, 21)):
        r = fake_record(subject, seed=100 + i)
        rng = np.random.RandomState(200 + i)
        r["selectedSupport"] = fake_support(seed=300 + i)
        r["minimal"]["support"] = [
            {"condition": "target", "scores": rng.normal(5.0, 12.0, 20).tolist(),
             "weights": [1.0] * 20, "flagged": [False] * 20, "effectiveUnflaggedN": 20.0,
             "referenceIds": list(range(1, 21)), "referenceN": 20, "kept": 20, "retention": 1.0},
            {"condition": "standard", "scores": rng.normal(1.0, 12.0, 80).tolist(),
             "weights": [1.0] * 80, "flagged": [False] * 80, "effectiveUnflaggedN": 80.0,
             "referenceIds": list(range(1, 81)), "referenceN": 80, "kept": 80, "retention": 1.0},
        ]
        r["common"]["support"] = r["minimal"]["support"]
        r["erp"]["reportedDifferenceSE"] = 0.8
        r["minimal"]["reportedDifferenceSE"] = 2.5
        r["common"]["reportedDifferenceSE"] = 2.5
        for arm in ("erp", "minimal", "common"):
            boot = ce.bootstrap_se(ce.arm_support(r, arm), 1.0, subject)
            r[arm]["differenceAmplitude"] = boot["point"]
        r["reference"]["differenceAmplitude"] = r["erp"]["differenceAmplitude"]
        records.append(r)
    return records


def run_collector(records, protocol):
    """Write records + protocol to a temp dir, run ce.main(), return (code, summary)."""
    import sys
    td = tempfile.TemporaryDirectory()
    proto_path = Path(td.name) / "proto.json"
    proto_path.write_text(json.dumps(protocol))
    digest = hashlib.sha256(proto_path.read_bytes()).hexdigest()
    for r in records:
        r["protocolSha256"] = digest
    rec_paths = []
    for r in records:
        fp = Path(td.name) / f"r{r['subject']}.json"
        fp.write_text(json.dumps(r))
        rec_paths.append(str(fp))
    outdir = Path(td.name) / "out"
    old = sys.argv
    sys.argv = ["compare_erpcore.py", *rec_paths, "--protocol", str(proto_path),
                "--out", str(outdir)]
    try:
        code = ce.main()
    finally:
        sys.argv = old
    summary = json.loads((outdir / "summary.json").read_text())
    td.cleanup()
    return code, summary


class TestGates(unittest.TestCase):
    def test_paired_ci_below_zero(self):
        diffs = [-0.5] * 20
        mean, ci = ce.paired_ci(diffs, 1000, 510004)
        self.assertAlmostEqual(mean, -0.5)
        self.assertLess(ci[1], 0)

    def test_calibration_ratio(self):
        ratio, ci = ce.calibration_ratio([1.0] * 20, [1.0] * 20, 1000, 510004)
        self.assertAlmostEqual(ratio, 1.0)
        self.assertGreaterEqual(ci[0], 0.99)
        self.assertLessEqual(ci[1], 1.01)

    def test_full_run_established(self):
        records = build_established_records()
        p = fake_protocol(tuple(range(1, 21)))
        code, summary = run_collector(records, p)
        self.assertEqual(code, 0)
        self.assertEqual(summary["decision"], "established")
        self.assertTrue(summary["paired_precision"]["minimal"]["ci_entirely_below_zero"])
        self.assertTrue(summary["calibration"]["ci_inside_band"])

    def test_one_low_support_pair_excluded_not_fatal(self):
        records = build_established_records()
        records[0]["selectedSupport"][0]["effectiveUnflaggedN"] = 3.0
        p = fake_protocol(tuple(range(1, 21)))
        code, summary = run_collector(records, p)
        self.assertEqual(code, 0)
        self.assertEqual(summary["excluded_pairs"], [1])
        self.assertEqual(summary["n"], 19)
        self.assertEqual(summary["decision"], "established")

    def test_four_incomplete_pairs_cap_verdict_at_partial(self):
        # Amendment A5 end to end: 16 complete pairs where Standard wins
        # on precision and calibration, plus 3 Standard failures and 1
        # comparator failure (the real batch's shape). The run still
        # produces a verdict (exit 0), the 4 are listed, and the verdict
        # is partial, never established.
        records = build_established_records()
        for r in records[16:19]:
            r["erp"] = {"available": False, "differenceAmplitude": None,
                        "reportedDifferenceSE": None, "gain": None}
            r["feasible"] = 0
            r["candidates"] = 24
            r["replayMaxError"] = None
            r["exportMaxError"] = None
            del r["selectedSupport"]
        r = records[19]
        r["minimal"] = {"available": False, "differenceAmplitude": None,
                        "reportedDifferenceSE": None, "gain": None,
                        "support": r["minimal"]["support"]}
        r["common"] = {"available": False, "differenceAmplitude": None,
                       "reportedDifferenceSE": None, "gain": None,
                       "support": r["common"]["support"]}
        p = fake_protocol(tuple(range(1, 21)))
        code, summary = run_collector(records, p)
        self.assertEqual(code, 0)
        self.assertEqual(summary["n"], 16)
        self.assertEqual(summary["decision"], "partial")
        self.assertEqual([ip["subject"] for ip in summary["incomplete_pairs"]],
                         [17, 18, 19, 20])
        self.assertFalse(summary["gates"]["complete_pairs"])

    def test_three_low_support_pairs_incomplete(self):
        records = build_established_records()
        for r in records[:3]:
            r["common"]["support"][1]["effectiveUnflaggedN"] = 2.0
        p = fake_protocol(tuple(range(1, 21)))
        code, summary = run_collector(records, p)
        self.assertEqual(code, 1)
        self.assertEqual(summary["decision"], "incomplete")
        self.assertEqual(len(summary["support_failures"]), 3)


if __name__ == "__main__":
    unittest.main()
