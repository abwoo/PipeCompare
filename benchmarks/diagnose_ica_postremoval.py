#!/usr/bin/env python3
"""Forensic readout of observed-vs-predicted post-ICA sensor change.

We do NOT change a gate or use an alternate prediction for a pass.
All outputs refer to original archival development subjects and seeds.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from validate_ica_unseen_seeds import grid
from validate_ica_postremoval import true_event_table, serial

FROZEN={11:[12,25],15:[5,2,6],21:[25,4]}
SOURCE_RUN=43405859824
GATE=1e-4

def stats(a):
    x=np.asarray(a,dtype=float)
    return {"maxAbs":float(np.max(np.abs(x))),"rmse":float(np.sqrt(np.mean(x*x))),
            "mean":float(np.mean(x)),"std":float(np.std(x))}

def diagnose(pre,post,decomp):
    before=true_event_table(pre,"pre")
    after=true_event_table(post,"post")
    if set(before)!=set(after) or any(before[i][0]!=after[i][0] for i in before):
        raise ValueError("Actual pre/post event mismatch")
    ids=serial(pre["eventIds"])
    count=int(pre["nComponents"])
    amplitude=np.asarray(pre["componentUv"],dtype=float)
    if amplitude.shape!=(count,len(ids)):
        raise ValueError("Source ICA component dimensions inconsistent")
    removed=[int(z) for z in serial(decomp.get("removed"))]
    if any(z<1 or z>count for z in removed):
        raise ValueError("Original removed IC index is invalid")
    observed=np.asarray([before[i][1]-after[i][1] for i in ids])
    projected=amplitude[np.asarray([i-1 for i in removed],dtype=int),:].sum(axis=0) if removed else np.zeros(len(ids))
    kept=amplitude[np.asarray([i for i in range(count) if i+1 not in removed],dtype=int),:].sum(axis=0)
    diff=observed-projected
    wrong_sign=observed+projected
    wrong_kept=observed-kept
    numerator=float(np.dot(observed-np.mean(observed),projected-np.mean(projected)))
    denom=float(np.sum((projected-np.mean(projected))**2))
    slope=numerator/denom if denom>1e-12 else None
    intercept=float(np.mean(observed)-(slope*np.mean(projected) if slope is not None else 0))
    ordered_deletes={idx:stats(observed-amplitude[idx-1]) for idx in removed}
    zero_kept=[int(i)+1 for i in range(count) if int(i)+1 not in removed and float(np.std(amplitude[i]))<1e-10]
    return {
        "actualRemovedIndices":removed,
        "icaComponents":count,
        "physicalSourceEvents":len(ids),
        "beforePz":stats([before[i][1] for i in ids]),
        "afterPz":stats([after[i][1] for i in ids]),
        "actualDeletionChange":stats(observed),
        "predictedDeletedComponentChange":stats(projected),
        "actualMinusPredicted":stats(diff),
        "differenceExceedsOriginalFrozenGate":bool(np.max(np.abs(diff))>GATE),
        "oppositeSignAlternativeError":stats(wrong_sign),
        "keptComponentAlternativeError":stats(wrong_kept),
        "bestAffineSlopeNotCausal":slope,
        "bestAffineInterceptNotCausal":intercept,
        "singleRemovedComponentAlternativeErrors":ordered_deletes,
        "keptNearConstantComponents":zero_kept,
        "preComponentSumVsPzMax":float(np.max(np.abs(amplitude.sum(axis=0)+np.asarray([v for v in pre["residualPzUv"]])-np.asarray([before[i][1] for i in ids])))),
        "postStageStatus":post.get("status"),
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--evidence",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    fs=sorted(a.evidence.rglob("result.json"))
    if not fs:raise ValueError("Expected actual archived EEG results")
    rows=[]
    for file in fs:
        source=json.loads(file.read_text())
        sid=source["subject"]
        if sid not in FROZEN:raise ValueError("Out-of-scope subject")
        roster=FROZEN[sid]
        c=source["icaSeedCheck"]
        grid0=grid(c["candidateMeasurements"],len(c["ids"]),"candidateMeasurements")
        candidates=[]
        for candidate in roster:
            i=c["ids"].index(candidate)
            per=[]
            for seed in range(6):
                m=grid0[i][seed]
                if m.get("candidateId")!=candidate or m.get("seed")!=seed:
                    raise ValueError("Original seed and candidate misaligned")
                z=m["science"]
                pre=serial(z.get("preRemovalContributions"))
                post=serial(z.get("postRemovalSignals"))
                dec=serial(z.get("decompositions"))
                if len(pre)!=len(post)!=1 or len(pre)!=1 or len(dec)!=1:
                    raise ValueError("Missing actual Pz observer")
                result=diagnose(pre[0],post[0],dec[0])
                result["seed"]=seed
                per.append(result)
                print("PR51_POST_REMOVAL_FORENSIC",json.dumps({
                    "subject":sid,"candidate":candidate,"seed":seed,
                    "removed":result["actualRemovedIndices"],
                    "frozenGateFailed":result["differenceExceedsOriginalFrozenGate"],
                    "observedMaxAbsUv":result["actualDeletionChange"]["maxAbs"],
                    "predictedMaxAbsUv":result["predictedDeletedComponentChange"]["maxAbs"],
                    "residualMaxAbsUv":result["actualMinusPredicted"]["maxAbs"],
                    "slope":result["bestAffineSlopeNotCausal"],
                    "intercept":result["bestAffineInterceptNotCausal"],
                    "oppositeSignResidualMaxUv":result["oppositeSignAlternativeError"]["maxAbs"],
                    "keptAlternativeResidualMaxUv":result["keptComponentAlternativeError"]["maxAbs"]}))
            candidates.append({"candidate":candidate,"seeds":per})
        rows.append({"subject":sid,"candidates":candidates})
    if len(set(row["subject"] for row in rows))!=len(rows):
        raise ValueError("Duplicate participant from archive")
    output={"study":"unmodified real ICA deletion prediction mismatch forensics",
            "sourceRun":SOURCE_RUN,"original20SubjectVerdict":"PARTIAL",
            "originalGateMaxAbsoluteAllowedUv":GATE,
            "gateModified":False,"productionCodeChanged":False,
            "noNewEEGFits":True,"noNewParticipants":True,"diagnosticNoPromotion":True,
            "results":rows}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(output,indent=2,allow_nan=False))
    print("PR51_POST_REMOVAL_FORENSIC_STOP",json.dumps({
        "originalScientificGateRelaxed":False,"passAuthorized":False,
        "PR51MergeAuthorized":False,"frozenHoldoutTouched":False}))

if __name__=="__main__":
    main()
