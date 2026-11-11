#!/usr/bin/env python3
"""Collector for the ERP CORE real-P3 validation.

Reads one result.json per participant (plus faults.csv / candidates.csv when
present), verifies execution identity, independently recomputes the primary
endpoint -- the per-participant bootstrap SE (1000 trial resamples) of the
P3 difference amplitude -- from the saved per-trial scores and weights, and
applies the frozen decision gates. No known truth exists in real data, so
there is no MSE here: precision is compared through bootstrap SEs, and
calibration through reported SE vs bootstrap SE.

Exit code 0 when all 20 planned records are present and valid, whatever the
scientific verdict; non-zero on missing/invalid records or identity failure.
The verdict is recorded in summary.json (decision tier) and summary.md.
"""
import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

ARMS = ["erp", "minimal", "common"]
BOOTSTRAP_DRAWS_TRIAL = 1000
BOOTSTRAP_SEED_TRIAL = 520004


def finite(x):
    return x is not None and isinstance(x, (int, float)) and math.isfinite(x)


def interval(draws, q):
    """Order statistic ceil((B+1)*q), clipped to [1,B]; same rule as PR51."""
    b = len(draws)
    k = int(math.ceil((b + 1) * q))
    k = max(1, min(b, k))
    return float(np.sort(draws)[k - 1])


def paired_ci(diffs, draws, seed):
    """Mean of paired diffs with a participant-bootstrap 95% CI."""
    diffs = np.asarray(diffs, dtype=float)
    rng = np.random.RandomState(seed)
    idx = rng.randint(0, len(diffs), (draws, len(diffs)))
    means = diffs[idx].mean(axis=1)
    return float(diffs.mean()), [interval(means, 0.025), interval(means, 0.975)]


def bootstrap_se(support_rows, gain, subject):
    """Bootstrap SE of the gain-corrected P3 difference amplitude.

    Resample trials within each condition with replacement (1000 draws),
    recompute the weighted-mean difference each time, take the SD. Uses only
    the saved per-trial scores and weights, so it is independent of the
    MATLAB SME and of the tool's reported error. Also returns the bootstrap
    mean (must agree with the MATLAB point estimate) and the per-condition
    trial counts used.
    """
    rng = np.random.RandomState(BOOTSTRAP_SEED_TRIAL + subject)
    good_data = []
    for row in support_rows:
        scores = np.asarray(row["scores"], dtype=float)
        weights = np.asarray(row["weights"], dtype=float)
        good = np.isfinite(scores) & np.isfinite(weights) & (weights > 0)
        if int(good.sum()) < 2:
            return None
        good_data.append((scores[good], weights[good]))
    diffs = []
    for s, w in good_data:
        idx = rng.randint(0, len(s), (BOOTSTRAP_DRAWS_TRIAL, len(s)))
        diffs.append((w[idx] * s[idx]).sum(axis=1) / w[idx].sum(axis=1))
    boot = (diffs[0] - diffs[1]) / gain
    point = (np.average(good_data[0][0], weights=good_data[0][1])
             - np.average(good_data[1][0], weights=good_data[1][1])) / gain
    return {"se": float(boot.std(ddof=1)), "mean": float(boot.mean()),
            "point": float(point),
            "counts": {row["condition"]: int(len(good_data[i][0]))
                       for i, row in enumerate(support_rows)}}


def validate_record(r, p, digest):
    """Identity and completeness checks. Raises AssertionError on failure."""
    assert r["algorithmCommit"] == p["algorithm_commit"], "Algorithm revision mismatch"
    assert r["protocolSha256"] == digest, "Protocol mismatch"
    assert r.get("measureSha256") == p["baseline_measure_sha256"], "Measure.m identity mismatch"
    assert r["eeglab"] == p["toolchain"]["eeglab"], "EEGLAB version mismatch"
    # MATLAB's `version` returns e.g. "26.1.0.3346908 (R2026a) Update 5", so the
    # pinned release token must be looked for inside the string, not at its
    # start. (A3's collect rejected all 20 real records because the old
    # startswith check assumed the string begins with the release name.)
    assert p["toolchain"]["matlab"] in r["matlab"], "MATLAB version mismatch"
    entry = next(x for x in p["source"]["files"] if x["subject"] == r["subject"])
    src = r["source"]
    assert src["subject"] == r["subject"], "Source subject mismatch"
    assert src["source_set_sha256"] == entry["set"]["sha256"], "Source .set identity mismatch"
    assert src["source_fdt_sha256"] == entry["fdt"]["sha256"], "Source .fdt identity mismatch"
    assert src["source_events_sha256"] == entry["events"]["sha256"], "Source events identity mismatch"
    assert src["sampling_rate"] == p["excerpt"]["sampling_rate"], "Excerpt sampling rate mismatch"
    assert src["samples"] == p["excerpt"]["duration_seconds"] * p["excerpt"]["sampling_rate"], "Excerpt length mismatch"
    assert src["n_valid_targets"] >= p["trial_rules"]["min_valid_targets"], "Valid-target floor missed"
    assert src["n_valid_standards"] >= p["trial_rules"]["min_valid_standards"], "Valid-standard floor missed"
    if not r["completed"]:
        return
    assert not r.get("errors"), f"Completed record contains errors: {r.get('errors')}"
    assert r.get("notRun") == 0, "Incomplete candidate family"
    opts = r["searchOptions"]
    assert opts["icaSeeds"] == 2 and opts["refineCutoffs"] and opts["eogCheck"] \
        and opts["signalTemplate"] == "data", "Standard options changed"
    assert finite(r.get("replayMaxError")) and r["replayMaxError"] < p["evaluation"]["replay_tolerance_uv"], "Replay error too large"
    assert finite(r.get("exportMaxError")) and r["exportMaxError"] < p["evaluation"]["replay_tolerance_uv"], "Export error too large"
    for arm in ARMS:
        a = r[arm]
        assert a.get("available"), f"{arm} readout unavailable"
        assert finite(a.get("differenceAmplitude")) and finite(a.get("reportedDifferenceSE")), f"{arm} difference/SE missing"
        assert a.get("support"), f"{arm} per-trial support missing"
        assert len(a["support"]) == 2, f"{arm} support must have two conditions"
    # Support is NOT asserted here: the frozen protocol allows up to 2 pairs
    # with effective unflagged support < 5 (listed separately, excluded from
    # the paired analysis); more than 2 makes the batch incomplete. See main().
    # Independent recomputation: the collector's bootstrap mean must agree
    # with the MATLAB difference amplitude; the bootstrap SE must be finite.
    boot = bootstrap_se(r["erp"]["support"], r["erp"]["gain"], r["subject"])
    assert boot is not None and finite(boot["se"]) and boot["se"] > 0, "Bootstrap SE unavailable"
    assert abs(boot["point"] - r["erp"]["differenceAmplitude"]) <= 1e-6, \
        f"Point estimate mismatch: collector {boot['point']:.6f} vs MATLAB {r['erp']['differenceAmplitude']:.6f}"


def calibration_ratio(actual, reported, draws, seed):
    """RMS(actual)/RMS(reported) with a participant-bootstrap 95% CI."""
    a = np.asarray(actual, dtype=float)
    b = np.asarray(reported, dtype=float)
    ratio = float(np.sqrt((a ** 2).mean() / (b ** 2).mean())) if (b ** 2).mean() > 0 else float("inf")
    rng = np.random.RandomState(seed)
    idx = rng.randint(0, len(a), (draws, len(a)))
    aa = (a[idx] ** 2).mean(axis=1)
    bb = (b[idx] ** 2).mean(axis=1)
    rr = np.full(len(aa), np.nan)
    pos = bb > 0
    rr[pos] = np.sqrt(aa[pos] / bb[pos])
    rr[(bb == 0) & (aa > 0)] = np.inf
    lo = rr.copy()
    hi = rr.copy()
    lo[np.isnan(lo)] = 0
    hi[np.isnan(hi)] = np.inf
    return ratio, [interval(lo, 0.025), interval(hi, 0.975)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("records", nargs="+", help="result.json files, one per participant")
    ap.add_argument("--protocol", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    p = json.loads(Path(args.protocol).read_text())
    # The protocol JSON does not pin its own hash (a file cannot contain its
    # own SHA256); the digest below is computed over the file as frozen, and
    # every record carries the digest from fetch time.
    digest = hashlib.sha256(Path(args.protocol).read_bytes()).hexdigest()

    records = []
    for path in args.records:
        records.append(json.loads(Path(path).read_text()))
    # Deduplicate by subject: a retried job leaves one artifact per attempt;
    # keep the highest run attempt for each subject.
    by_subject = {}
    for r in records:
        s = r["subject"]
        att = r.get("runAttempt")
        try:
            att = int(att)
        except (TypeError, ValueError):
            att = -1
        if s not in by_subject or att > by_subject[s][0]:
            by_subject[s] = (att, r)
    records = [v[1] for v in by_subject.values()]
    records.sort(key=lambda r: r["subject"])
    planned = [f["subject"] for f in p["source"]["files"]]
    present = [r["subject"] for r in records]
    missing = [s for s in planned if s not in present]
    extra = [s for s in present if s not in planned]

    # One execution identity across the whole run: no mixed commits/runs/EEGLAB.
    exec_commits = {r.get("executionCommit") for r in records}
    runs = {r.get("workflowRun") for r in records}
    eeglab_shas = {r.get("eeglabSha") for r in records}

    valid = []
    problems = []
    for r in records:
        try:
            validate_record(r, p, digest)
            valid.append(r)
        except AssertionError as e:
            problems.append({"subject": r["subject"], "problem": str(e)})
    incomplete = [r["subject"] for r in records if not r.get("completed")]
    complete = [r for r in valid if r.get("completed")]

    out = {"protocol_sha256": digest, "planned": planned, "present": present,
           "missing": missing, "extra": extra, "incomplete": incomplete,
           "problems": problems,
           "execution_commits": sorted(c for c in exec_commits if c),
           "workflow_runs": sorted(c for c in runs if c),
           "eeglab_shas": sorted(c for c in eeglab_shas if c),
           "mixed_execution": len(exec_commits) > 1 or len(runs) > 1 or len(eeglab_shas) > 1}

    if missing or extra or problems or incomplete or out["mixed_execution"]:
        out["decision"] = "incomplete"
        write_out(args.out, out)
        print(json.dumps(out, indent=2))
        return 1

    ev = p["evaluation"]
    gates = p["decision_gates"]
    # Support gate (frozen protocol): a pair with effective unflagged support
    # < 5 in any arm/condition is not a valid pair for the main conclusion.
    # Up to 2 such pairs are listed separately and excluded from the paired
    # analysis; more than 2 makes the batch incomplete.
    min_sup = gates["minimum_effective_unflagged"]
    support_fail = []
    for r in complete:
        bad = []
        for arm in ARMS:
            for row in r[arm]["support"]:
                if row["effectiveUnflaggedN"] < min_sup:
                    bad.append(f"{arm}/{row['condition']}={row['effectiveUnflaggedN']:.2f}")
        if bad:
            support_fail.append({"subject": r["subject"], "details": bad})
    out["support_failures"] = support_fail
    if len(support_fail) > 2:
        out["decision"] = "incomplete"
        out["incomplete_reason"] = f"{len(support_fail)} pairs below the support floor (max 2 allowed)"
        write_out(args.out, out)
        print(json.dumps(out, indent=2))
        return 1
    valid = [r for r in complete
             if r["subject"] not in {f["subject"] for f in support_fail}]
    out["valid_pairs"] = [r["subject"] for r in valid]
    out["excluded_pairs"] = [f["subject"] for f in support_fail]

    # Primary endpoint per participant and arm (valid pairs only).
    se = {arm: [] for arm in ARMS}
    diff_amp = {arm: [] for arm in ARMS}
    reported = {arm: [] for arm in ARMS}
    for r in valid:
        for arm in ARMS:
            boot = bootstrap_se(r[arm]["support"], r[arm]["gain"], r["subject"])
            se[arm].append(boot["se"])
            diff_amp[arm].append(r[arm]["differenceAmplitude"])
            reported[arm].append(r[arm]["reportedDifferenceSE"])
    se = {k: np.array(v) for k, v in se.items()}
    out["n"] = len(valid)
    out["bootstrap_se_uv"] = {arm: [float(x) for x in se[arm]] for arm in ARMS}
    out["difference_amplitude_uv"] = {arm: [float(x) for x in diff_amp[arm]] for arm in ARMS}
    out["mean_bootstrap_se_uv"] = {arm: float(se[arm].mean()) for arm in ARMS}
    out["mean_difference_amplitude_uv"] = {arm: float(np.mean(diff_amp[arm])) for arm in ARMS}

    # Paired precision: Standard vs each fixed control.
    precision = {}
    for ctrl in ("minimal", "common"):
        d = se["erp"] - se[ctrl]
        mean, ci = paired_ci(d, ev["bootstrap_draws"], ev["bootstrap_seed"])
        precision[ctrl] = {"mean_diff_uv": mean, "ci95": ci,
                           "ci_entirely_below_zero": ci[1] < 0,
                           "wins": int((d < 0).sum()), "losses": int((d > 0).sum())}
        # Individual regressions: Standard worse by >20% relative AND beyond
        # the 0.05 uV absolute allowance.
        regs = [bool(se["erp"][i] > se[ctrl][i] * (1 + gates["individual_regression_relative_allowance"])
                      and se["erp"][i] - se[ctrl][i] > gates["individual_regression_absolute_allowance_uv"])
                for i in range(len(valid))]
        precision[ctrl]["individual_regressions"] = int(sum(regs))
        precision[ctrl]["regression_subjects"] = [valid[i]["subject"] for i, x in enumerate(regs) if x]
    out["paired_precision"] = precision

    # Calibration (Standard arm): reported SE vs bootstrap SE.
    ratio, ci = calibration_ratio(se["erp"], reported["erp"], ev["bootstrap_draws"], ev["bootstrap_seed"])
    lo, hi = ev["calibration_agreement_band"]
    out["calibration"] = {"arm": "erp", "ratio": ratio, "ci95": ci,
                          "band": [lo, hi], "ci_inside_band": ci[0] >= lo and ci[1] <= hi,
                          "ci_outside_band": ci[1] < lo or ci[0] > hi}

    # Stability: ICA-seed agreement across participants (same metric as PR51).
    stable = [bool(r.get("icaSeedCheck", {}).get("same")) for r in complete]
    out["stability"] = {"ica_seed_same_count": int(sum(stable)), "n": len(complete),
                        "gate_pass": int(sum(stable)) >= gates["stability_min_count"]}

    # Secondary: agreement of the Standard difference amplitude with the
    # published-pipeline reference amplitude.
    ref = [r["reference"]["differenceAmplitude"] for r in valid
           if r.get("reference", {}).get("available")]
    if len(ref) == len(valid):
        d = np.array(diff_amp["erp"]) - np.array(ref)
        out["reference_agreement"] = {"n": len(ref), "mean_abs_diff_uv": float(np.abs(d).mean()),
                                      "mean_diff_uv": float(d.mean()),
                                      "correlation": float(np.corrcoef(diff_amp["erp"], ref)[0, 1])}

    # Decision tiers (frozen draft wording):
    # established = precision vs both controls + calibration + stability +
    #   individual regressions all pass.
    # partial = individual regressions pass AND (won minimal but not common,
    #   OR won both but calibration or stability fails). Only what actually
    #   passed is claimed; never upgraded.
    # not-established = complete and valid but neither of the above
    #   (in particular both precision CIs failing to sit below 0).
    prec_min = precision["minimal"]["ci_entirely_below_zero"]
    prec_com = precision["common"]["ci_entirely_below_zero"]
    cal_ok = out["calibration"]["ci_inside_band"]
    stab_ok = out["stability"]["gate_pass"]
    reg_ok = all(precision[c]["individual_regressions"] <= gates["individual_regression_max_count"]
                 for c in ("minimal", "common"))
    if prec_min and prec_com and cal_ok and stab_ok and reg_ok:
        out["decision"] = "established"
    elif reg_ok and ((prec_min and not prec_com)
                     or (prec_min and prec_com and (not cal_ok or not stab_ok))):
        out["decision"] = "partial"
    else:
        out["decision"] = "not-established"
    out["gates"] = {"precision_vs_minimal": prec_min, "precision_vs_common": prec_com,
                    "calibration": cal_ok, "stability": stab_ok, "individual_regressions": reg_ok}
    write_out(args.out, out)
    print(json.dumps(out, indent=2))
    return 0


def write_out(out_dir, out):
    d = Path(out_dir)
    d.mkdir(parents=True, exist_ok=True)
    (d / "summary.json").write_text(json.dumps(out, indent=2))
    lines = ["# ERP CORE real-P3 validation — collector summary",
             "", f"Decision: **{out['decision']}**", "",
             f"Records: {out.get('n', 0)} valid pairs of {len(out['planned'])} planned.",
             f"Missing: {out['missing'] or 'none'}; problems: {len(out['problems'])}.",
             f"Excluded pairs (support < 5): {out.get('excluded_pairs') or 'none'}.",
             f"Mixed execution identity: {out['mixed_execution']}.", ""]
    if "paired_precision" in out:
        for ctrl, v in out["paired_precision"].items():
            lines.append(f"Precision Standard vs {ctrl}: mean diff {v['mean_diff_uv']:.3f} uV, "
                         f"95% CI [{v['ci95'][0]:.3f}, {v['ci95'][1]:.3f}], "
                         f"wins {v['wins']}, losses {v['losses']}, "
                         f"individual regressions {v['individual_regressions']}.")
        c = out["calibration"]
        lines.append(f"Calibration (Standard, reported vs bootstrap SE): ratio {c['ratio']:.3f}, "
                     f"95% CI [{c['ci95'][0]:.3f}, {c['ci95'][1]:.3f}], band {c['band']}.")
        s = out["stability"]
        lines.append(f"ICA-seed stability: {s['ica_seed_same_count']}/{s['n']} same recommendation.")
        if "reference_agreement" in out:
            ra = out["reference_agreement"]
            lines.append(f"Reference agreement (Standard vs published pipeline): mean |diff| "
                         f"{ra['mean_abs_diff_uv']:.3f} uV, correlation {ra['correlation']:.3f}.")
        lines.append("")
        lines.append(f"Gates: {json.dumps(out['gates'])}")
    (d / "summary.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    sys.exit(main())
