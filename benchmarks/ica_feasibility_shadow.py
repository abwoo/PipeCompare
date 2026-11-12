#!/usr/bin/env python3
"""Scientific-gate-aware ICA minimax SHADOW decision. Never mutates Rank.

Requires real per-seed candidate feasibility, candidate stratum, and
per-seed measurement/support. A finite objective without a feasible status
is NOT eligible. No scientific gate or original PR51 verdict is relaxed.
"""
import argparse
import json
import math
from pathlib import Path

FLOOR = 5.0
SEEDS = 3


def isfinite(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def seed_grid(values, candidates, label):
    """Accept MATLAB jsonencode's column-major flattened cell arrays.

    Numeric objective tables are true Nx3 arrays, but MATLAB encodes Nx3
    cell arrays as a flat list (seed0 candidates, seed1, seed2). Do not
    reshape row-major or cross candidate/seed scientific verdicts.
    """
    if not isinstance(values, list):
        raise ValueError(f"{label} must be a JSON array")
    if (len(values)==candidates
            and all(isinstance(row,list) and len(row)==SEEDS for row in values)):
        return values
    if len(values)==candidates*SEEDS:
        return [[values[seed*candidates+i] for seed in range(SEEDS)]
                for i in range(candidates)]
    raise ValueError(f"{label} candidate/seed alignment mismatch")


def shadow(study):
    if not study.get("completed") or study.get("errors"):
        raise ValueError("incomplete MATLAB source run")
    check = study.get("icaSeedCheck")
    if not isinstance(check, dict):
        raise ValueError("seed-check data unavailable")
    ids, scores = check["ids"], check["objective"]
    n = len(ids)
    if n == 0 or len(set(ids)) != n or len(scores) != n:
        raise ValueError("misaligned candidate IDs")
    if not all(isinstance(row,list) and len(row)==SEEDS for row in scores):
        raise ValueError("objective must have three aligned seed columns")
    statuses = seed_grid(check["rankStatus"], n, "rankStatus")
    strata = seed_grid(check["strata"], n, "strata")
    measure = seed_grid(check["candidateMeasurements"], n, "candidateMeasurements")
    initial = check["recommended"][0]
    if initial is None:
        raise ValueError("baseline recommendation absent")
    if initial not in ids:
        raise ValueError("baseline recommendation is outside re-evaluated ICA candidates")
    initial_row = ids.index(initial)
    initial_stratum = strata[initial_row][0]
    eligible = []
    audit = []
    for i,ident in enumerate(ids):
        errors = []
        if any(value != initial_stratum for value in strata[i]):
            errors.append("stratum-differs-from-baseline-recommended")
        if any(value != "feasible" for value in statuses[i]):
            errors.append("not-feasible-in-every-seed")
        if not all(isfinite(value) for value in scores[i]):
            errors.append("nonfinite-objective")
        values = []
        for seed,m in enumerate(measure[i]):
            if not isinstance(m, dict) or m.get("candidateId") != ident or m.get("seed") != seed:
                errors.append("mismatched-candidate-measurement")
                continue
            amp=m.get("correctedDifferenceEstimate")
            support=m.get("effectiveUnflagged")
            gain=m.get("gain")
            if not m.get("hasMeasurement") or not isfinite(amp) or not isfinite(gain) or gain<=0:
                errors.append("missing-finite-p3-measurement")
            if not isinstance(support,list) or len(support)!=2 or not all(
                    isfinite(v) and v+1e-9>=FLOOR for v in support):
                errors.append("frozen-effective-trial-floor-not-met")
            if isfinite(amp):
                values.append(amp)
        ok=not errors
        audit.append({"id":ident,"eligible":ok,"reasons":sorted(set(errors)),
                      "statuses":statuses[i],"objectives":scores[i]})
        if ok:
            eligible.append({"id":ident,"maxObjective":max(scores[i]),
                             "medianObjective":sorted(scores[i])[1],
                             "calibratedP3BySeed":values,
                             "p3Range":max(values)-min(values)})
    minWorst=min((a["maxObjective"] for a in eligible),default=None)
    tied=[x for x in eligible if x["maxObjective"]==minWorst]
    # Lacking a frozen cross-seed deterministic tie-breaker, explicitly abstain.
    recommended=tied[0]["id"] if len(tied)==1 else None
    out={"subject":study["subject"],"sourceRecommended":check["recommended"],
         "sourceCandidateStable":check["same"],"baselineStratum":initial_stratum,
         "eligibleCount":len(eligible),"eligible":eligible,
         "shadowCandidate":recommended,
         "shadowCandidateUnique":len(tied)==1,
         "reason":("unique-feasible-worst-seed-objective"
                   if recommended is not None else
                   "no-common-feasible-candidate" if not eligible else
                   "objective-tie-requires-predefined-tiebreak"),
         "audit":audit, "method":"exploratory-shadow-only",
         "originalRecommendationUnchanged":True,
         "scientificThresholdsUnchanged":True,
         "scientificPromotionAuthorized":False}
    if recommended is not None:
        choice=tied[0]
        out["shadowP3Range"]=choice["p3Range"]
        out["shadowP3BySeed"]=choice["calibratedP3BySeed"]
    return out


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--inputs",type=Path,nargs="+",required=True)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    results=[shadow(json.loads(path.read_text())) for path in args.inputs]
    if len({a["subject"] for a in results})!=len(results):
        raise ValueError("duplicated participant")
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps({
        "status":"exploratory-not-scientific-confirmation",
        "frozen_floor":FLOOR,"subjects":results},indent=2))
    for result in results:
        print("PR51_ICA_FEASIBILITY_SHADOW",json.dumps({
            k:result[k] for k in ("subject","eligibleCount","shadowCandidate",
             "shadowCandidateUnique","reason","sourceRecommended",
             "sourceCandidateStable")}))

if __name__=="__main__":
    main()
