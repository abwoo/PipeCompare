#!/usr/bin/env python3
"""Read-only audit of physical support versus weighted effective support.

No scientific gate is relaxed, no trial is restored, no candidate is selected.
Used on the immutable six-start real-EEG development artifacts.
"""
import argparse
import json
import math
from pathlib import Path

SEEDS = 6
FLOOR = 5.
TOL = 1e-8


def finite_number(x):
    return isinstance(x,(int,float)) and not isinstance(x,bool) and math.isfinite(x)


def matrix(x,n,name):
    if not isinstance(x,list):
        raise ValueError(f"{name}: missing array")
    if len(x)==n and all(isinstance(row,list) and len(row)==SEEDS for row in x):
        return x
    if len(x)==n*SEEDS:
        return [[x[s*n+i] for s in range(SEEDS)] for i in range(n)]
    raise ValueError(f"{name}: misaligned MATLAB cell grid")


def check_measurement(m,candidate,seed,has_measurement=None):
    """Check unavoidable mathematical bounds, independent of candidate rank."""
    issues=[]
    if not isinstance(m,dict) or m.get("candidateId")!=candidate or m.get("seed")!=seed:
        return ["candidate-seed-misaligned"]
    raw=m.get("sourceUnflagged")
    eff=m.get("effectiveUnflagged")
    if not isinstance(raw,list) or not isinstance(eff,list) or len(raw)!=2 or len(eff)!=2:
        return ["support-pair-missing"]
    for i,(n,w) in enumerate(zip(raw,eff)):
        label="target" if i==0 else "standard"
        if not finite_number(n) or n<0 or n!=int(n):
            issues.append(f"{label}-invalid-physical-count")
            continue
        if not finite_number(w) or w < -TOL:
            issues.append(f"{label}-invalid-effective-n")
        elif w > n + TOL:
            issues.append(f"{label}-effective-exceeds-physical")
        # If a candidate is declared eligible, physical support must be >= floor.
        if finite_number(w) and w >= FLOOR-TOL and n<FLOOR:
            issues.append(f"{label}-eligible-with-insufficient-physical-trials")
    return issues


def audit(result):
    if result.get("completed") is not True or result.get("errors"):
        raise ValueError("Incomplete MATLAB study record")
    c=result["icaSeedCheck"]
    ids=c["ids"]
    if not isinstance(ids,list) or not ids or len(ids)!=len(set(ids)):
        raise ValueError("Candidate IDs invalid")
    ms=matrix(c["candidateMeasurements"],len(ids),"candidateMeasurements")
    findings=[]
    for i,candidate in enumerate(ids):
        for seed in range(SEEDS):
            problems=check_measurement(ms[i][seed],candidate,seed)
            if problems:
                findings.append({"candidateId":candidate,"seed":seed,"issues":problems})
    return {"subject":result["subject"],"ids":len(ids),
            "examinedCandidateStarts":len(ids)*SEEDS,
            "findings":findings,"supportIntegrityPassed":not findings,
            "originalRecommendationUnchanged":True,
            "productionPromotionAuthorized":False}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--evidence",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    paths=sorted(args.evidence.rglob("result.json"))
    if len(paths)!=3:
        raise SystemExit(f"Expected exactly 3 original real-EEG results, saw {len(paths)}")
    rows=[audit(json.loads(path.read_text())) for path in paths]
    if sorted(r["subject"] for r in rows)!=[11,15,21]:
        raise SystemExit("Unexpected development cohort")
    rows.sort(key=lambda x:x["subject"])
    output={"study":"read-only support bounds across all original ICA candidate/seed fits",
            "rawRunId":43169562528,"newEEGFits":False,"newSubjects":False,
            "original20SubjectVerdict":"PARTIAL","rows":rows,
            "allPassed":all(r["supportIntegrityPassed"] for r in rows),
            "decision":"audit-only; do not promote robust selector or merge PR51"}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(output,indent=2))
    for row in rows:
        print("PR51_ICA_SUPPORT_INTEGRITY",json.dumps({
            "subject":row["subject"],"candidateStarts":row["examinedCandidateStarts"],
            "passed":row["supportIntegrityPassed"],"findings":len(row["findings"])}))
    print("PR51_ICA_SUPPORT_INTEGRITY_SUMMARY",json.dumps({
        "allPassed":output["allPassed"],"totalFindings":sum(len(x["findings"]) for x in rows)}))
    # Preserve evidence artifact even with a scientific failure. The workflow
    # exit code must not be mistaken for evidence of improvement.
    if not output["allPassed"]:
        raise SystemExit("Scientific support integrity failure; do not promote")


if __name__=="__main__":
    main()
