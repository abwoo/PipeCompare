#!/usr/bin/env python3
"""Independent science-gated estimator-family ICA robustness cross-check.

Development-only P300. Original Rank-run gates and effective >=5 support are
not modified, no candidate is promoted by this audit.
"""
import argparse
import json
import math
from pathlib import Path
from validate_ica_robust_real import validate as validate_common


def finite(x):
    return isinstance(x,(float,int)) and not isinstance(x,bool) and math.isfinite(x)


def validate(record):
    original=validate_common(record)
    s=record["icaSeedCheck"]
    family=s["familyRobust"]
    rows=s["candidateSeedAudit"]
    initial=s["recommended"][0]
    by_id={r["candidateId"]:r for r in rows}
    if initial not in by_id:
        raise ValueError("initial ICA recommendation absent from evaluated candidates")
    if len(by_id)!=len(rows):
        raise ValueError("candidate ids are duplicated")
    if family.get("estimatorFamilyLocked") is not True or family.get("originalRecommendationUnchanged") is not True:
        raise ValueError("method provenance or original rank immutability broken")
    baseline_family=by_id[initial]["estimationFamily"][0]
    if baseline_family not in ("weighted","unweighted"):
        raise ValueError("unknown original P300 estimator family")
    if family["baselineId"]!=initial or family["baselineFamily"]!=baseline_family:
        raise ValueError("baseline estimator-family mismatch")
    candidates=[]
    for r in rows:
        (statuses,strata,support,objective,families)=(r["status"],r["stratum"],r["supportPass"],r["objective"],r["estimationFamily"])
        if not all(isinstance(v,list) and len(v)==3 for v in (statuses,strata,support,objective,families)):
            raise ValueError("per-seed candidate status and estimator provenance not aligned")
        if not all(type(v) is bool for v in support):
            raise ValueError("independent HMM-support mask must be explicit booleans")
        if not all(isinstance(f,str) and f in {"weighted","unweighted","unavailable"} for f in families):
            raise ValueError("unknown per-seed estimator family")
        if (all(f==baseline_family for f in families) and
            all(z=="feasible" for z in statuses) and all(support) and
            len(set(strata))==1 and
            all(finite(v) and v>=0 for v in objective)):
            candidates.append((r["candidateId"],max(objective),strata[0]))
    feasible_count=len(candidates)
    groups={c[2] for c in candidates}
    result={"subject":record["subject"],
            "originalId":initial,"baselineFamily":baseline_family,
            "allSeedScientificAndFamilyFeasibleCount":feasible_count,
            "oldCommonFeasibleCandidate":original["robustId"],
            "familyCandidate":family["candidateId"],
            "familyStatus":family["status"],
            "familyMeasurements":family.get("seedMeasurements",[]),
            "noRetrospectiveGateChange":True,
            "isNewFrozenExperiment":False}
    if family["status"]=="selected":
        if not candidates:
            raise ValueError("family-constrained minimax selected candidate fails scientific gates")
        if len(groups)!=1:
            raise ValueError("cannot select across incomparable scientific strata")
        best=min(z[1] for z in candidates)
        winners=[z[0] for z in candidates if abs(z[1]-best)<=1e-12*max(1,abs(best))]
        if len(winners)!=1 or family["candidateId"]!=winners[0]:
            raise ValueError("family-constrained minimax candidate is incorrect or not feasible")
        if not finite(family["scientific"]["worstObjective"]) or not math.isclose(
            family["scientific"]["worstObjective"],best,abs_tol=1e-9,rel_tol=0):
            raise ValueError("family minimax objective inconsistent with all-seed science matrix")
        measures=family.get("seedMeasurements",[])
        if len(measures)!=3:
            raise ValueError("missing family-locked candidate repeat measurements")
        for i,m in enumerate(measures):
            if m["seed"]!=i or m["candidateId"]!=family["candidateId"] or m["estimationFamily"]!=baseline_family:
                raise ValueError("selected recipe changed scientific estimator family across seeds")
            effective=m["effectiveUnflagged"]
            if not (len(effective)==2 and all(finite(n) and n>=5 for n in effective)):
                raise ValueError("scientific HMM-unflagged effective trial floor changed")
            if not (finite(m["gain"]) and m["gain"]>0 and
                    finite(m["differenceEstimate"]) and finite(m["correctedDifferenceEstimate"]) and
                    math.isclose(m["differenceEstimate"]/m["gain"],m["correctedDifferenceEstimate"],
                                 abs_tol=1e-8,rel_tol=0)):
                raise ValueError("invalid original-versus-family calibrated P300")
    elif family["status"]=="abstain":
        if family["candidateId"] is not None or family.get("seedMeasurements"):
            raise ValueError("abstention cannot silently carry a recommendation")
        if candidates and len(groups)==1:
            best=min(z[1] for z in candidates)
            winners=[z for z in candidates if abs(z[1]-best)<=1e-12*max(1,abs(best))]
            if len(winners)==1:
                raise ValueError("unique science-feasible family candidate was incorrectly abstained")
    else:
        raise ValueError("unrecognized family robustness state")
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--input",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    result=validate(json.loads(args.input.read_text()))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2))
    print("PR51_ICA_ESTIMATOR_FAMILY_RESULT",json.dumps({
        "subject":result["subject"],"baselineFamily":result["baselineFamily"],
        "oldCommonFeasible":result["oldCommonFeasibleCandidate"],
        "familyCandidate":result["familyCandidate"],"eligible":result["allSeedScientificAndFamilyFeasibleCount"],
        "P3":[m["correctedDifferenceEstimate"] for m in result["familyMeasurements"]]}))

if __name__=="__main__":
    main()
