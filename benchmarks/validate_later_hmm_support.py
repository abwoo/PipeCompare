#!/usr/bin/env python3
"""Scientific completeness gate for longer-window EEG artifact support.

A negative HMM support result is valid evidence. Missing data, a changed
300-second baseline, invented support, and different source SHA fail closed.
"""
import argparse
import json
from pathlib import Path

EXPECTED={4:{"seconds":544,"before":[18,95],"unflagged":[4,3]},
          8:{"seconds":460,"before":[22,95],"unflagged":[5,20]}}


def check(x):
    sid=x.get("subject")
    if sid not in EXPECTED:
        raise ValueError("unplanned subject")
    want=EXPECTED[sid]
    if x.get("originalSeconds")!=300 or x.get("fullSeconds")!=want["seconds"]:
        raise ValueError("wrong source duration")
    if x.get("original20SubjectVerdict")!="PARTIAL" or x.get("productionPromotionAuthorized") is not False:
        raise ValueError("retrospective scientific promotion")
    if x.get("sameHMMScienceAndTrialRules") is not True:
        raise ValueError("scientific test changed its rules")
    a=x.get("baseline");b=x.get("extended")
    if not (isinstance(a,list) and isinstance(b,list) and len(a)==len(b)==2):
        raise ValueError("missing conditions")
    if [r.get("condition") for r in a]!=["target","standard"]:
        raise ValueError("original condition order mismatch")
    if [r.get("condition") for r in b]!=["target","standard"]:
        raise ValueError("extension condition order mismatch")
    if [r.get("referenceN") for r in a]!=want["before"]:
        raise ValueError("frozen 300-second source count not reproduced")
    if [r.get("unflaggedN") for r in a]!=want["unflagged"]:
        raise ValueError("original 300-second HMM result not reproduced")
    for group in (a,b):
        for item in group:
            n=item.get("referenceN"); k=item.get("unflaggedN"); flagged=item.get("flaggedN")
            ids=item.get("unflaggedIds")
            fi=item.get("flaggedIds")
            if (not all(isinstance(v,int) and not isinstance(v,bool) for v in (n,k,flagged))
                    or min(n,k,flagged)<0 or k+flagged!=n):
                raise ValueError("HMM flag cardinality inconsistent")
            def array(v):
                if isinstance(v,list): return v
                if isinstance(v,(int,float)) and v==int(v): return [int(v)]
                return None
            kept=array(ids); marked=array(fi)
            if kept is None or marked is None or len(kept)!=k or len(marked)!=flagged:
                raise ValueError("missing per-trial HMM identities")
            if len(set(kept+marked))!=n:
                raise ValueError("overlapping/duplicate reference trial identities")
    full_ok=all(r["unflaggedN"]>=5 for r in b)
    if full_ok!=x.get("fullWindowUnflaggedFloorMet"):
        raise ValueError("incorrect full-support verdict")
    return {"subject":sid,"originalUnflagged":[r["unflaggedN"] for r in a],
            "fullUnflagged":[r["unflaggedN"] for r in b],
            "fullReferenceN":[r["referenceN"] for r in b],
            "fullHMMFloorMet":full_ok,
            "thisIsOnlySourceHMMFeasibility":True,
            "pairedMeasurementNowValid":False,
            "independentParticipantConfirmation":False}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--evidence",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    sources=sorted(a.evidence.rglob("later_hmm_support.json"))
    if len(sources)!=2:
        raise SystemExit(f"Expected 2 preserved full-run HMM results, found {len(sources)}")
    rows=[check(json.loads(f.read_text())) for f in sources]
    if {r["subject"] for r in rows}!={4,8}:
        raise SystemExit("Duplicated/missing original development participants")
    rows.sort(key=lambda r:r["subject"])
    result={"status":"development-only HMM reference support",
            "originalFrozenVerdict":"PARTIAL",
            "original20SubjectValidPairs":15,
            "noThresholdRelaxation":True,"sourceSubjects":[4,8],
            "newHeldoutSubjects":[ ],"productionPromotionAuthorized":False,
            "rows":rows}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2))
    for row in rows:
        print("PR51_LONGER_HMM_COLLECT",json.dumps(row))
if __name__=="__main__":
    main()
