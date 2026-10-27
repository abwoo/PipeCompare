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

import csv
import json
import math
import sys
from pathlib import Path

FAILURES = []


def fail(msg):
    FAILURES.append(msg)


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


def main():
    if len(sys.argv) < 5:
        print("usage: check_shadow_selection.py <fragments_dir> <final_json> "
              "<per_seed_csv> <summary_csv> [expected_head_sha]")
        return 2
    frag_dir, final_json, per_seed_csv, summary_csv = map(Path, sys.argv[1:5])
    expected_head = sys.argv[5] if len(sys.argv) > 5 else None

    run = json.loads(final_json.read_text())
    frags = {}
    for path in sorted(frag_dir.rglob("shadow_selection_seed_*.json")):
        frag = json.loads(path.read_text())
        frags[frag["seed"]] = frag

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
    if prov.get("inputArtifactRunId") is not None:
        fail("provenance.inputArtifactRunId must be null for fresh seeds")
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
        frag_arms = {a["name"]: a for a in frag["arms"]}
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
    if run["coverage"].get("planned") != len(per_seed):
        fail("coverage.planned mismatch")

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
