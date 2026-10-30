"""Independent checker for shadow-selection run output (variant A).

Recomputes, from the per-seed fragments (not from the run JSON's own
summary fields): paired differences and means on arm intersections,
the percentile bootstrap interval from the stored replicates, coverage
counts, verdicts, and the CSV contents. Compares every recomputed value
with the run JSON / CSVs and fails on any mismatch.

This implements the unavailable mode of benchmarks/
diagnose_paired_regression.py for shadow JSON: an arm with no unique
assessable recommendation is unavailable, a pair uses only seeds
assessable in BOTH arms, and single-arm-only seeds are listed, never
folded into a mean.

Percentile rule: linear interpolation, the same rule as MATLAB
prctile's default (0-based position (n-1)*p/100), implemented here so
the check needs no third-party package. Standard library only.

Usage:
  check_shadow_selection.py <fragments_dir> <final_json> <per_seed_csv>
      <summary_csv> [expected_head_sha]
Exit 0 on PASS, 1 on any FAIL (all failures are printed).
"""

import copy
import csv
import json
import math
import sys
from pathlib import Path

FAILURES = []


def fail(msg):
    FAILURES.append(msg)


def fast_mode_of(frag):
    """Round 10 fast-mode field of a fragment ({} when absent)."""
    fm = frag.get("fastMode") or {}
    return fm if isinstance(fm, dict) else {}


def is_fast_fragment(frag):
    """True for micro / slice-diagnostic / sliced-validation / design fragments."""
    fm = fast_mode_of(frag)
    mode = fm.get("mode") or "full"
    return mode != "full" or bool(fm.get("designSeedDiagnostic")) or bool(fm.get("partialSearch")) or bool(fm.get("slicedValidation"))


def validate_micro_fragment(frag):
    """Validate a fast fragment ON ITS OWN terms (never held-out).

    Returns a list of failures (empty = valid fast fragment). This
    does NOT certify the fragment as evidence: a valid fast fragment
    is still diagnostic only, and is_fast_fragment() must keep it
    out of every held-out paired computation.
    """
    failures = []
    fm = fast_mode_of(frag)
    mode = fm.get("mode")
    if mode not in ("micro", "slice-diagnostic", "resume-fast"):
        failures.append(f"fast fragment mode {mode!r} is not a known fast mode")
    if frag.get("seed") not in (40101, 40102, 40103, 40104) and mode != "resume-fast":
        failures.append(f"fast fragment seed {frag.get('seed')} is not a design seed 40101-40104")
    if not fm.get("designSeedDiagnostic") and mode != "resume-fast":
        failures.append("fast fragment must carry designSeedDiagnostic=true")
    labels = frag.get("labels") or {}
    if labels.get("fullyHeldOut") is True:
        failures.append("fast fragment must not claim fullyHeldOut=true")
    return failures


def is_nan(x):
    return x is None or (isinstance(x, float) and math.isnan(x))


def num(x):
    """JSON number-or-null to float (null means NaN here)."""
    if x is None:
        return float("nan")
    return float(x)


def close(a, b, tol=1e-9):
    """NaN-aware closeness: both NaN counts as equal."""
    if is_nan(a) and is_nan(b):
        return True
    if is_nan(a) or is_nan(b):
        return False
    return abs(a - b) <= tol + tol * abs(b)


def percentile_sorted(xs, p):
    """Linear percentile on an already-sorted list (MATLAB prctile default)."""
    n = len(xs)
    if n == 0:
        return float("nan")
    if n == 1:
        return xs[0]
    pos = (n - 1) * p / 100.0
    lo = math.floor(pos)
    hi = math.ceil(pos)
    if lo == hi:
        return xs[lo]
    return xs[lo] + (pos - lo) * (xs[hi] - xs[lo])


def percentile(xs, p):
    return percentile_sorted(sorted(xs), p)


def verdict_for(mean_diff, ci_lower, ci_upper, tol):
    if is_nan(mean_diff):
        return "failed", "no paired seeds"
    if ci_lower > tol:
        return "failed", "interval lower bound above +5% tolerance"
    if mean_diff <= 0 and ci_upper <= tol:
        return "established", ""
    return "not-established", "point or interval does not meet the frozen rule"


def parse_csv_num(text):
    if text is None or text.strip() in ("", "NaN", "nan", "null"):
        return float("nan")
    return float(text)



def matlab_array(value):
    """Accept MATLAB jsonencode's empty, singleton and vector representations."""
    if value is None:
        return []
    return value if isinstance(value, list) else [value]


def normalize_run_arrays(run):
    """Normalize only schema-defined array fields; preserve values/provenance."""
    run = copy.deepcopy(run)
    for field in ("seeds", "perSeed", "unavailable"):
        if field in run:
            run[field] = matlab_array(run[field])
    for entry in run.get("perSeed", []):
        if isinstance(entry, dict) and "arms" in entry:
            entry["arms"] = matlab_array(entry["arms"])
    for pair in run.get("pairedDifferences", {}).values():
        for field in ("seeds", "gainCorrected", "raw"):
            if field in pair:
                pair[field] = matlab_array(pair[field])
    for ci in run.get("recordBootstrapCI95", {}).values():
        if "replicates" in ci:
            ci["replicates"] = matlab_array(ci["replicates"])
    return run

def main():
    FAILURES.clear()
    if len(sys.argv) < 5:
        print("usage: check_shadow_selection.py <fragments_dir> <final_json> "
              "<per_seed_csv> <summary_csv> [expected_head_sha] [planned]")
        return 2
    frag_dir, final_json, per_seed_csv, summary_csv = map(Path, sys.argv[1:5])
    expected_head = sys.argv[5] if len(sys.argv) > 5 else None
    # Round 8 B1: optional 6th argument is the PLANNED seed count.
    # Without it, coverage.planned cannot be verified (the old check
    # compared it to the present count, which is tautological).
    expected_planned = int(sys.argv[6]) if len(sys.argv) > 6 else None

    run = normalize_run_arrays(json.loads(final_json.read_text()))
    # Round 6 arm fragments: shadow_selection_seed_<seed>_standard.json
    # and shadow_selection_seed_<seed>_balanced.json combine into one
    # per-seed fragment. Complete per-seed fragments
    # (shadow_selection_seed_<seed>.json) are used as-is. When both
    # arm fragments are present, their contract fingerprints must
    # agree, otherwise the run's Standard arm must be unavailable
    # (the merge step enforces this; the checker verifies it here).
    frags = {}
    arm_parts = {}
    for path in sorted(frag_dir.rglob("shadow_selection_seed_*.json")):
        frag = json.loads(path.read_text())
        name = path.name
        # Round 8 B1: filename seed must agree with content seed.
        import re as _re
        _m = _re.match(r"shadow_selection_seed_(\d+)_(standard|balanced|)\.json$", name)
        if _m and _m.group(1) and int(_m.group(1)) != frag.get("seed"):
            fail(f"fragment {name} filename seed disagrees with content seed {frag.get('seed')}")
        if name.endswith("_standard.json") or name.endswith("_balanced.json"):
            seed = frag["seed"]
            entry = arm_parts.setdefault(seed, {"arms": {}, "fingerprints": {}})
            # A single arm encodes as a JSON object, not an array.
            frag_arms = frag.get("arms", [])
            if isinstance(frag_arms, dict):
                frag_arms = [frag_arms]
            for arm in frag_arms:
                # Round 8 B1: duplicate seed-arm fragments must fail.
                if arm["name"] in entry["arms"]:
                    fail(f"duplicate fragment arm {arm['name']} for seed {seed}")
                entry["arms"][arm["name"]] = arm
            if "guardChanged" in frag:
                entry["guardChanged"] = frag["guardChanged"]
            if "guardDecision" in frag:
                entry["guardDecision"] = frag["guardDecision"]
            kind = "standard" if name.endswith("_standard.json") else "balanced"
            entry["fingerprints"][kind] = frag.get("contract")
            entry["labels"] = frag.get("labels")
            entry["provenance"] = frag.get("provenance")
            # Round 10: carry the fastMode field so the held-out
            # guard below can see it; without this, a micro arm
            # fragment would lose its label when combined.
            if frag.get("fastMode") is not None:
                entry["fastMode"] = frag.get("fastMode")
        else:
            frags[frag["seed"]] = frag
    for seed, entry in arm_parts.items():
        # Fill a missing arm with an unavailable placeholder, matching
        # the merge step's 'arm fragment missing' arm (checker
        # compares assessable/MSE/pick, not the reason string).
        arms = []
        for k in ("standard", "balancedGuarded", "balancedUnguarded"):
            if k in entry["arms"]:
                arms.append(entry["arms"][k])
            else:
                arms.append({"name": k, "pick": None, "selectedKey": "",
                             "mseGainCorrected": None, "mseRaw": None,
                             "assessable": False, "reason": "arm fragment missing"})
        combined = {"seed": seed,
                    "arms": arms,
                    "labels": entry.get("labels"),
                    "provenance": entry.get("provenance")}
        if entry.get("fastMode") is not None:
            combined["fastMode"] = entry.get("fastMode")
        if "guardChanged" in entry:
            combined["guardChanged"] = entry["guardChanged"]
        if "guardDecision" in entry:
            combined["guardDecision"] = entry["guardDecision"]
        # Fingerprint check: both present and differing means the
        # run's Standard arm must be unavailable.
        fps = entry["fingerprints"]
        if "standard" in fps and "balanced" in fps:
            if fps["standard"] != fps["balanced"]:
                std = combined["arms"] and next(
                    (a for a in combined["arms"] if a["name"] == "standard"), None)
                if std is not None and std.get("assessable"):
                    fail(f"seed {seed} contract fingerprints differ but "
                         "Standard arm is assessable in the run")
        frags[seed] = combined

    # ---- Round 10 fast guard: a micro / sliced / design fragment
    # is diagnostic only and must never enter a held-out paired
    # computation. Fail loudly here, in the checker, even if a
    # caller accidentally passes fast fragments to this script.
    for _seed, _frag in frags.items():
        if is_fast_fragment(_frag):
            fail(f"fragment seed {_seed} is a fast-mode fragment "
                 f"(fastMode={fast_mode_of(_frag)}); it is diagnostic "
                 "only and must never be checked as held-out evidence")
    if FAILURES:
        print("FAIL:")
        for f in FAILURES:
            print("  " + f)
        return 1

    # ---- required top-level fields ----
    required = ["commit", "seeds", "recipes", "perSeed", "pairedDifferences",
                "pairedMeanDifference", "recordBootstrapCI95", "coverage",
                "provenance", "labels", "frozenCriteria"]
    for field in required:
        if field not in run:
            fail(f"missing top-level field: {field}")
    if FAILURES:
        print("FAIL:")
        for f in FAILURES:
            print("  " + f)
        return 1

    # ---- labels (variant A constants) ----
    labels = run["labels"]
    for key, want in [("fullyHeldOut", True), ("preprocessingRefitted", False),
                      ("adaptiveLeavesIncluded", False)]:
        if labels.get(key) is not want:
            fail(f"label {key} is {labels.get(key)!r}, expected {want!r}")
    if not run.get("labelDefinitions"):
        fail("labelDefinitions missing or empty")

    # ---- recipes ----
    want_recipes = {
        "standard": ("standard", "legacy-total", True),
        "balancedGuarded": ("balanced", "pooled", True),
        "balancedUnguarded": ("balanced", "pooled", False),
    }
    for arm, (recipe, variance, protect) in want_recipes.items():
        got = run["recipes"].get(arm, {})
        if (got.get("recipe"), got.get("artifactVariance"),
                got.get("protectReferenceCleanup")) != (recipe, variance, protect):
            fail(f"recipe {arm} mismatch: {got}")

    # ---- provenance ----
    prov = run["provenance"]
    if prov.get("inputArtifactRunId") not in (None, []):
        fail("provenance.inputArtifactRunId must be null or an empty array for fresh seeds")
    if expected_head is not None:
        if prov.get("commit") != expected_head:
            fail(f"provenance.commit {prov.get('commit')} != expected {expected_head}")
        if expected_head[:7] not in final_json.name:
            fail(f"final JSON name {final_json.name} lacks head short SHA {expected_head[:7]}")
    if not prov.get("script"):
        fail("provenance.script missing")
    if not run["frozenCriteria"].get("text") or "+5%" not in run["frozenCriteria"]["text"]:
        fail("frozenCriteria text missing or lacks the +5% rule")

    # ---- seeds vs fragments ----
    if sorted(run["seeds"]) != sorted(frags):
        fail(f"run seeds {sorted(run['seeds'])} != fragment seeds {sorted(frags)}")
    if run["seeds"] != sorted(run["seeds"]):
        fail("run seeds are not sorted")

    # ---- per-seed arms: fragments must match the run JSON ----
    per_seed = {}
    for entry in run["perSeed"]:
        arms = {a["name"]: a for a in entry["arms"]}
        per_seed[entry["seed"]] = arms
        frag = frags.get(entry["seed"])
        if frag is None:
            continue
        frag_arms = {a["name"]: a for a in matlab_array(frag["arms"])}
        for arm in want_recipes:
            fa, ra = frag_arms.get(arm), arms.get(arm)
            if fa is None or ra is None:
                fail(f"seed {entry['seed']} arm {arm} missing in fragment or run JSON")
                continue
            if bool(fa["assessable"]) != bool(ra["assessable"]):
                fail(f"seed {entry['seed']} arm {arm} assessable mismatch")
            for key in ("mseGainCorrected", "mseRaw"):
                if not close(num(fa.get(key)), num(ra.get(key))):
                    fail(f"seed {entry['seed']} arm {arm} {key} fragment/run mismatch")
            if fa.get("pick") != ra.get("pick"):
                fail(f"seed {entry['seed']} arm {arm} pick fragment/run mismatch")

    # ---- paired recomputation ----
    pairs = {"guardedMinusStandard": ("balancedGuarded", "standard"),
             "guardedMinusUnguarded": ("balancedGuarded", "balancedUnguarded")}
    for pair, (first, second) in pairs.items():
        paired = [s for s in per_seed
                  if per_seed[s][first]["assessable"] and per_seed[s][second]["assessable"]]
        diffs_gc = [num(per_seed[s][first]["mseGainCorrected"])
                    - num(per_seed[s][second]["mseGainCorrected"]) for s in paired]
        diffs_raw = [num(per_seed[s][first]["mseRaw"])
                     - num(per_seed[s][second]["mseRaw"]) for s in paired]
        got_diff = run["pairedDifferences"][pair]
        if sorted(got_diff["seeds"]) != sorted(paired):
            fail(f"pair {pair} paired seeds mismatch: {got_diff['seeds']} vs {paired}")
        for got, want, name in [(got_diff["gainCorrected"], diffs_gc, "gainCorrected"),
                                (got_diff["raw"], diffs_raw, "raw")]:
            # Order by seed for a stable comparison.
            order = sorted(range(len(paired)), key=lambda i: paired[i])
            want_ordered = [want[i] for i in order]
            if len(got) != len(want_ordered) or any(
                    not close(g, w) for g, w in zip(got, want_ordered)):
                fail(f"pair {pair} {name} differences mismatch")
        mean_diff = sum(diffs_gc) / len(diffs_gc) if diffs_gc else float("nan")
        if not close(num(run["pairedMeanDifference"][pair]), mean_diff):
            fail(f"pair {pair} mean difference mismatch")
        ci = run["recordBootstrapCI95"][pair]
        if ci.get("draws") != 10000 or ci.get("seed") != 1:
            fail(f"pair {pair} bootstrap draws/seed mismatch: {ci.get('draws')}/{ci.get('seed')}")
        reps = ci.get("replicates", [])
        if diffs_gc:
            if len(reps) != 10000:
                fail(f"pair {pair} replicate count {len(reps)} != 10000")
            else:
                lo, hi = percentile(reps, 2.5), percentile(reps, 97.5)
                if not close(num(ci["lower"]), lo) or not close(num(ci["upper"]), hi):
                    fail(f"pair {pair} CI mismatch: JSON [{ci['lower']}, {ci['upper']}] "
                         f"vs recomputed [{lo}, {hi}]")
        ref_vals = [num(per_seed[s][second]["mseGainCorrected"]) for s in paired]
        ref_mean = sum(ref_vals) / len(ref_vals) if ref_vals else float("nan")
        tol = 0.05 * ref_mean if not is_nan(ref_mean) else float("nan")
        verdict, _ = verdict_for(mean_diff, num(ci["lower"]), num(ci["upper"]), tol)
        # The run JSON stores verdicts under 'verdicts' when present.
        stored = run.get("verdicts", {}).get(pair, {})
        if stored.get("verdict") != verdict:
            fail(f"pair {pair} verdict mismatch: {stored.get('verdict')} vs {verdict}")
        if stored and not close(num(stored.get("referenceMean")), ref_mean):
            fail(f"pair {pair} referenceMean mismatch")
        if stored and not close(num(stored.get("tolerance")), tol):
            fail(f"pair {pair} tolerance mismatch")
        if run["coverage"].get("paired_" + pair) != len(paired):
            fail(f"pair {pair} coverage mismatch")
        # Negative control: identical picks must give a zero difference.
        for s in paired:
            gp = per_seed[s]["balancedGuarded"].get("pick")
            up = per_seed[s]["balancedUnguarded"].get("pick")
            if pair == "guardedMinusUnguarded" and gp is not None and gp == up:
                d = num(per_seed[s]["balancedGuarded"]["mseGainCorrected"]) \
                    - num(per_seed[s]["balancedUnguarded"]["mseGainCorrected"])
                if abs(d) > 1e-12:
                    fail(f"seed {s} negative control violated: identical picks, diff {d}")

    # ---- coverage counts ----
    for arm in want_recipes:
        want_n = sum(1 for s in per_seed if per_seed[s][arm]["assessable"])
        if run["coverage"].get(arm) != want_n:
            fail(f"coverage {arm} mismatch: {run['coverage'].get(arm)} vs {want_n}")
    if expected_planned is not None:
        if run.get("nSeedsPlanned") != expected_planned:
            fail(f"nSeedsPlanned {run.get('nSeedsPlanned')} != planned {expected_planned}")
        if run["coverage"].get("planned") != expected_planned:
            fail(f"coverage.planned {run['coverage'].get('planned')} != planned {expected_planned}")
    # Without a planned argument, coverage.planned is NOT verified:
    # comparing it to len(per_seed) would be tautological (Round 8 B1).

    # ---- per-seed CSV ----
    with per_seed_csv.open(newline="") as f:
        rows = list(csv.DictReader(f))
    if len(rows) != len(per_seed) * 3:
        fail(f"per-seed CSV row count {len(rows)} != {len(per_seed) * 3}")
    for row in rows:
        seed, arm = int(row["seed"]), row["arm"]
        if seed not in per_seed or arm not in per_seed[seed]:
            fail(f"per-seed CSV unknown row: {seed} {arm}")
            continue
        got = per_seed[seed][arm]
        if not close(parse_csv_num(row["selectedMSE_gainCorrected"]),
                     num(got["mseGainCorrected"])):
            fail(f"per-seed CSV gainCorrected mismatch: {seed} {arm}")
        if not close(parse_csv_num(row["selectedMSE_raw"]), num(got["mseRaw"])):
            fail(f"per-seed CSV raw mismatch: {seed} {arm}")
        if int(row["assessable"]) != int(bool(got["assessable"])):
            fail(f"per-seed CSV assessable mismatch: {seed} {arm}")
        # Round 8 B1: the per-seed CSV must carry pick/selectedKey.
        if "pick" not in row or "selectedKey" not in row:
            fail("per-seed CSV missing pick/selectedKey columns")
        else:
            got_pick = got.get("pick")
            csv_pick = row["pick"]
            if (got_pick is None and csv_pick not in ("", "None")) or \
               (got_pick is not None and csv_pick != str(got_pick)):
                fail(f"per-seed CSV pick mismatch: {seed} {arm}")
            if row["selectedKey"] != (got.get("selectedKey") or ""):
                fail(f"per-seed CSV selectedKey mismatch: {seed} {arm}")

    # ---- summary CSV ----
    with summary_csv.open(newline="") as f:
        srows = {r["pair"]: r for r in csv.DictReader(f)}
    for pair in pairs:
        row = srows.get(pair)
        if row is None:
            fail(f"summary CSV missing pair {pair}")
            continue
        if int(row["nPaired"]) != run["coverage"]["paired_" + pair]:
            fail(f"summary CSV nPaired mismatch: {pair}")
        if not close(parse_csv_num(row["meanDifference"]),
                     num(run["pairedMeanDifference"][pair])):
            fail(f"summary CSV meanDifference mismatch: {pair}")
        ci = run["recordBootstrapCI95"][pair]
        if not close(parse_csv_num(row["ciLower"]), num(ci["lower"])) or \
                not close(parse_csv_num(row["ciUpper"]), num(ci["upper"])):
            fail(f"summary CSV CI mismatch: {pair}")
        stored = run.get("verdicts", {}).get(pair, {})
        if row["verdict"] != stored.get("verdict"):
            fail(f"summary CSV verdict mismatch: {pair}")

    if FAILURES:
        print("FAIL:")
        for failure in FAILURES:
            print("  " + failure)
        return 1
    print(f"PASS: {len(per_seed)} seeds, "
          f"{sum(1 for s in per_seed for a in per_seed[s] if per_seed[s][a]['assessable'])} "
          f"assessable arm records, all paired values, CIs, verdicts and CSVs match")
    return 0


if __name__ == "__main__":
    sys.exit(main())
