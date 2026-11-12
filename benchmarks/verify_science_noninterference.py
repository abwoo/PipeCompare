#!/usr/bin/env python3
"""Compare recorder-enabled real EEG against original non-recorder run.

No source EEG refit here; paired checks refer to separately completed runs.
This detects intervention effects of the observational science recorder.
"""
import argparse
import json
import math
from pathlib import Path

SUBJECTS={11,15,21}


def close(a,b,tol=1e-7):
    if a is None or b is None:return a is None and b is None
    if isinstance(a,(int,float)) and isinstance(b,(int,float)):
        if math.isnan(a) or math.isnan(b):
            return math.isnan(a) and math.isnan(b)
        return math.isclose(a,b,rel_tol=0,abs_tol=tol)
    return a==b


def one_pair(old,new):
    s=old["subject"]
    reasons=[]
    if s!=new["subject"]:reasons.append("subject-mismatch")
    for key in ("completed","feasible","candidates","notRun"):
        if old.get(key)!=new.get(key):reasons.append(f"{key}-changed")
    if old.get("errors")!=new.get("errors"):reasons.append("execution-errors-changed")
    for key in ("source_set_sha256","source_fdt_sha256","source_events_sha256"):
        if old["source"].get(key)!=new["source"].get(key):
            reasons.append("source-"+key+"-changed")
    for key in ("differenceAmplitude","gain"):
        if not close(old["erp"].get(key),new["erp"].get(key)):
            reasons.append("erp-"+key+"-changed")
    a,b=old["icaSeedCheck"],new["icaSeedCheck"]
    for key in ("ids","recommended","same"):
        if a.get(key)!=b.get(key):reasons.append("ica-"+key+"-changed")
    for key in ("rankStatus","strata"):
        if a.get(key)!=b.get(key):reasons.append("ica-"+key+"-changed")
    ao,bo=a.get("objective"),b.get("objective")
    if len(ao)!=len(bo) or any(len(x)!=len(y) or any(not close(u,v) for u,v in zip(x,y))
                                  for x,y in zip(ao,bo)):
        reasons.append("ica-objective-changed")
    # Scalar measured effect is independent from observational snapshot.
    for j in range(len(a["measurements"])):
        x,y=a["measurements"][j],b["measurements"][j]
        for key in ("candidateId","hasMeasurement","seed"):
            if not close(x.get(key),y.get(key)):
                reasons.append(f"seed-{j}-{key}-changed")
        for key in ("correctedDifferenceEstimate","gain"):
            if not close(x.get(key),y.get(key)):
                reasons.append(f"seed-{j}-{key}-changed")
    numerical=[]
    for seed,(x,y) in enumerate(zip(a["measurements"],b["measurements"])):
        left=x.get("correctedDifferenceEstimate")
        right=y.get("correctedDifferenceEstimate")
        delta=(right-left) if isinstance(left,(int,float)) and isinstance(right,(int,float)) and math.isfinite(left) and math.isfinite(right) else None
        numerical.append({"seed":seed,"originalP300":left,"recorderP300":right,"p300DeltaUv":delta,
                          "originalCandidate":x.get("candidateId"),
                          "recorderCandidate":y.get("candidateId")})
    return {"subject":s,"noninterventionPassed":not reasons,"reasons":sorted(set(reasons)),
            "originalIds":a["recommended"],"recorderIds":b["recommended"],
            "seedP300Comparison":numerical}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original",type=Path,required=True)
    p.add_argument("--recorder",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    def load(folder):
        files=sorted(folder.rglob("result.json"))
        rows={}
        for f in files:
            x=json.loads(f.read_text())
            if x["subject"] in rows:raise ValueError("duplicate subject")
            rows[x["subject"]]=x
        if set(rows)!=SUBJECTS:raise ValueError(f"Expected unchanged three-subject cohort, got {set(rows)}")
        return rows
    a,b=load(args.original),load(args.recorder)
    rows=[one_pair(a[k],b[k]) for k in sorted(SUBJECTS)]
    out={"originalRun":43169562528,"recorderRun":43189639057,
         "newSubjects":False,"originalPR51Verdict":"PARTIAL",
         "allNonInterventional":all(x["noninterventionPassed"] for x in rows),
         "rows":rows,"promotionAuthorized":False}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,indent=2))
    for r in rows:print("PR51_SCIENCE_NONINTERFERENCE",json.dumps(r))
    if not out["allNonInterventional"]:raise SystemExit("Recorder changed one or more scientific outputs")


if __name__=="__main__":
    main()
