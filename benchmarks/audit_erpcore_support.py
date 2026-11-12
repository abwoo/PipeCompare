#!/usr/bin/env python3
"""Read-only, per-trial ERP CORE support reconciliation for PR51.

Inputs are exploratory erpcore_pilot result.json files; no thresholds or
pipeline decisions are changed here. Missing arms remain explicit.
"""
import argparse
import csv
import json
import math
from pathlib import Path

ARMS = ("erp", "minimal", "common")
SUPPORT_FLOOR = 5  # Original frozen protocol; never tuned from results.

def effective(weights):
    s = sum(weights)
    q = sum(w * w for w in weights)
    return s * s / q if q else 0.0

def audit(result):
    if not result.get("completed"):
        raise ValueError("incomplete result must not enter support audit")
    if result.get("notRun", 0) != 0:
        raise ValueError("incomplete candidate search")
    subject = int(result["subject"])
    rows, trials = [], []
    for arm in ARMS:
        readout = result.get(arm) or {}
        support = (result.get("selectedSupport") if arm == "erp"
                   else readout.get("support"))
        if not support:
            for condition in ("target", "standard"):
                rows.append({"subject": subject, "arm": arm, "condition": condition,
                             "available": False, "measurementAvailable": bool(readout.get("available", False)),
                             "reason": "support-unavailable"})
            continue
        if isinstance(support, dict):
            support = [support]
        for record in support:
            name = str(record["condition"])
            ids = record["referenceIds"]
            scores = record["scores"]
            weights = record["weights"]
            blink = record["flagged"]
            if not all(len(a) == len(ids) for a in (scores, weights, blink)):
                raise ValueError(f"unaligned reference vectors: {subject}/{arm}/{name}")
            if len(set(map(str, ids))) != len(ids):
                raise ValueError(f"duplicate reference trial ids: {subject}/{arm}/{name}")
            clean, flagged = [], []
            count = 0
            for ident, score, weight, is_blink in zip(ids, scores, weights, blink):
                w = float(weight) if weight is not None else float("nan")
                s = float(score) if score is not None else float("nan")
                positive = math.isfinite(w) and w > 0 and math.isfinite(s)
                actual_w = w if positive else 0.0
                if positive:
                    count += 1
                if bool(is_blink):
                    flagged.append(actual_w)
                else:
                    clean.append(actual_w)
                trials.append({"subject": subject, "arm": arm, "condition": name,
                               "reference_id": ident, "score": score,
                               "weight": actual_w, "blink_proxy": bool(is_blink),
                               "positive_weight": positive,
                               "effective_contribution": actual_w})
            eff_clean = effective(clean)
            if count != int(record["positiveN"]):
                raise ValueError(f"positiveN mismatch: {subject}/{arm}/{name}")
            if abs(eff_clean - float(record["effectiveUnflaggedN"])) > 1e-7:
                raise ValueError(f"effective support mismatch: {subject}/{arm}/{name}")
            rows.append({"subject": subject, "arm": arm, "condition": name,
                         "available": True, "measurementAvailable": bool(readout.get("available", False)),
                         "referenceN": int(record["referenceN"]),
                         "referenceUnflaggedN": sum(not bool(x) for x in blink),
                         "kept": int(record["kept"]), "positiveN": count,
                         "flaggedN": int(record["flaggedN"]),
                         "unflaggedN": int(record["unflaggedN"]),
                         "effectiveN": effective(clean + flagged),
                         "effectiveUnflaggedN": eff_clean,
                         "eyeBiasAssessed": bool(record["eyeBiasAssessed"]),
                         "meets_original_floor": eff_clean >= SUPPORT_FLOOR,
                         "reason": ("below-original-effective-floor"
                                    if eff_clean < SUPPORT_FLOOR else "support-sufficient")})
    return rows, trials

def structural_triage(rows):
    """An immutable HMM mask bounds effective support by clean source-trial count.

    This is an upper bound, not an achieved accuracy estimate; no threshold or
    scientific gate is modified. A subject with fewer than five unflagged
    reference trials cannot attain five effective unflagged trials through
    weighting, rejection or a different candidate selection.
    """
    out = []
    for subject in sorted(set(x["subject"] for x in rows)):
        standard = [x for x in rows if x["subject"] == subject and x["arm"] == "erp"]
        if len(standard) != 2 or not all(x["available"] for x in standard):
            out.append({"subject": subject, "structurally_possible": None,
                        "status": "source-support-unavailable"})
            continue
        def label(value):
            return {"11": "target", "31": "standard"}.get(value, value)
        clean = {label(x["condition"]): x["referenceUnflaggedN"] for x in standard}
        if set(clean) != {"target", "standard"}:
            raise ValueError(f"unknown reference condition identities for {subject}: {clean}")
        possible = all(n >= SUPPORT_FLOOR for n in clean.values())
        out.append({"subject": subject, "source_unflagged": clean,
                    "structurally_possible": possible,
                    "status": ("source-upper-bound-below-frozen-floor"
                               if not possible else "not-structurally-excluded"),
                    "note": "Necessary condition only; not evidence of artifact-free ground truth"})
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--inputs", nargs="+", required=True, type=Path)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    rows, trials = [], []
    for path in args.inputs:
        r, t = audit(json.loads(path.read_text()))
        rows += r
        trials += t
    args.out.mkdir(parents=True, exist_ok=True)
    for filename, data in (("support_conditions.csv", rows), ("support_trials.csv", trials)):
        if not data:
            continue
        keys = list(dict.fromkeys(k for record in data for k in record.keys()))
        with (args.out / filename).open("w", newline="") as file:
            writer = csv.DictWriter(file, fieldnames=keys)
            writer.writeheader()
            writer.writerows(data)
    (args.out / "summary.json").write_text(json.dumps({
        "source": "exploratory-real-P3-results",
        "original_effective_floor": SUPPORT_FLOOR,
        "no_thresholds_changed": True,
        "do_not_pool_with_original_preregistered_run": True,
        "conditions": rows,
        "structural_triage": structural_triage(rows),
        "trial_records": len(trials)
    }, indent=2))
    print("PR51_SUPPORT_AUDIT", json.dumps({
        "subjects": sorted(set(x["subject"] for x in rows)),
        "condition_rows": len(rows), "trial_rows": len(trials),
        "available_rows": sum(bool(x["available"]) for x in rows)}))

if __name__ == "__main__":
    main()
