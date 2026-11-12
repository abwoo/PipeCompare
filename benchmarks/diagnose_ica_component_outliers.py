#!/usr/bin/env python3
"""Read-only decomposition of worst matched ICA scalp-map correlations.

Explains failed frozen multimetric-v1 topology screen; cannot re-select or
promote a candidate, adjust thresholds, or open new participant data.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment
from validate_ica_unseen_seeds import grid
from validate_ica_science import FROZEN_DEVELOPMENT_IDS

LOW_CORR=.90
SEEDS=6


def serial(v):
    if v is None: return []
    if isinstance(v,list): return v
    return [v]


def centered_maps(dec):
    labels=serial(dec.get("channelLabels"))
    A=np.asarray(dec.get("mixing",[]),dtype=float)
    if (A.ndim!=2 or A.shape!=(len(labels),int(dec.get("nComponents",0)))
            or len(labels)!=len(set(labels)) or len(labels)<4
            or not np.isfinite(A).all()):
        raise ValueError("Missing/misaligned source mixing maps")
    A=A-A.mean(axis=0)
    norm=np.linalg.norm(A,axis=0)
    if np.any(norm<1e-12):
        raise ValueError("Degenerate source scalp map")
    return labels,A/norm


def labels_for(dec):
    names=serial(dec.get("classNames",[]))
    probs=np.asarray(dec.get("classProbabilities",[]),dtype=float)
    n=int(dec["nComponents"])
    if probs.size==0:
        return [None]*n
    if probs.shape!=(n,len(names)) or not np.isfinite(probs).all():
        raise ValueError("Incorrect ICLabel probability dimensions")
    return [names[int(i)] for i in np.argmax(probs,axis=1)]


def map_pairs(first,later):
    """One matched decomposition; sign-invariant and permutation-invariant."""
    a=serial(first.get("decompositions"))
    b=serial(later.get("decompositions"))
    if len(a)!=1 or len(b)!=1:
        raise ValueError("Expected one ICA fit in the prespecified ERP CORE recipes")
    a,b=a[0],b[0]
    labs,A=centered_maps(a)
    labs2,B=centered_maps(b)
    if set(labs)!=set(labs2):
        raise ValueError("Different ICA electrode montages")
    B=B[[labs2.index(label) for label in labs],:]
    corr=np.abs(A.T @ B)
    ii,jj=linear_sum_assignment(-corr)
    ra=set(int(x) for x in serial(a.get("removed",[])))
    rb=set(int(x) for x in serial(b.get("removed",[])))
    if any(x<1 or x>A.shape[1] for x in ra) or any(x<1 or x>B.shape[1] for x in rb):
        raise ValueError("Removed component outside decomposition")
    label_a=labels_for(a)
    label_b=labels_for(b)
    rows=[]
    for i,j in zip(ii,jj):
        ia,ib=int(i)+1,int(j)+1
        removed_a=ia in ra
        removed_b=ib in rb
        rows.append({
            "baselineComponent":ia,"otherComponent":ib,
            "absoluteMapCorrelation":float(corr[i,j]),
            "baselineRemoved":removed_a,"otherRemoved":removed_b,
            "removalConcordant":removed_a==removed_b,
            "baselineTopClass":label_a[i],
            "otherTopClass":label_b[j],
            "classConcordant":label_a[i]==label_b[j] if label_a[i] is not None and label_b[j] is not None else None,
        })
    rms=[r["absoluteMapCorrelation"] for r in rows]
    rm_a={k for k,r in enumerate(rows) if r["baselineRemoved"]}
    rm_b={k for k,r in enumerate(rows) if r["otherRemoved"]}
    union=rm_a|rm_b
    jacc=1. if not union else len(rm_a & rm_b)/len(union)
    bad=[r for r in rows if r["absoluteMapCorrelation"] < LOW_CORR]
    return {
        "matchedComponents":len(rows),"minimumCorrelation":min(rms),
        "medianCorrelation":float(np.median(rms)),
        "matchedRemovedJaccard":jacc,
        "belowPoint90":len(bad),
        "lowCorrelationBothKept":sum(not r["baselineRemoved"] and not r["otherRemoved"] for r in bad),
        "lowCorrelationBothRemoved":sum(r["baselineRemoved"] and r["otherRemoved"] for r in bad),
        "lowCorrelationRemovalDiscordant":sum(not r["removalConcordant"] for r in bad),
        "lowCorrelationClassDiscordant":sum(r["classConcordant"] is False for r in bad),
        "lowCorrelationPairs":bad,
    }


def examine(original,validation):
    if (original.get("subject")!=validation.get("subject") or
        original.get("subject") not in FROZEN_DEVELOPMENT_IDS or
        not original.get("completed") or original.get("errors") or
        validation.get("identityErrors") or validation.get("scienceFailures") or
        not validation.get("unalteredFrozenRankAndEstimates")):
        raise ValueError("Invalid/inconsistent archived source evidence")
    c=original["icaSeedCheck"]
    ids=c["ids"]
    m=grid(c["candidateMeasurements"],len(ids),"candidateMeasurements")
    results=[]
    for candidate in FROZEN_DEVELOPMENT_IDS[original["subject"]]:
        if candidate not in ids:
            raise ValueError("Predeclared candidate missing")
        r=ids.index(candidate)
        snaps=[m[r][s] for s in range(SEEDS)]
        if any(not isinstance(x,dict) or x.get("candidateId")!=candidate or
               x.get("seed")!=i or not isinstance(x.get("science"),dict) or
               x["science"].get("status")!="complete" for i,x in enumerate(snaps)):
            raise ValueError("Incomplete/misaligned scientific decomposition")
        expected=next((x for x in validation["predeclaredCandidateScience"]
                       if x["candidateId"]==candidate),None)
        if expected is None:raise ValueError("Candidate science summary absent")
        summaries={x["seed"]:x for x in expected["componentMatchesVsSeed0"]}
        per_seed=[]
        for seed in range(1,SEEDS):
            z=map_pairs(snaps[0]["science"],snaps[seed]["science"])
            summary=summaries[seed]
            if (not summary.get("available") or
                not math.isclose(z["minimumCorrelation"],summary["minimumMatchedTopographicCorrelation"],abs_tol=1e-8) or
                not math.isclose(z["medianCorrelation"],summary["medianAbsoluteTopographicCorrelation"],abs_tol=1e-8) or
                not math.isclose(z["matchedRemovedJaccard"],summary["removedComponentMatchedJaccard"],abs_tol=1e-8)):
                raise ValueError("Component matching deviates from verified recorder summary")
            per_seed.append({"seed":seed,**z})
        results.append({"candidateId":candidate,"starts":per_seed})
    return {"subject":original["subject"],"candidates":results,
            "newSubjects":False,"candidateSelectionChanged":False,
            "scientificThresholdChanged":False,"promotionAuthorized":False}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--evidence",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    raws=sorted(a.evidence.rglob("result.json"))
    val=sorted(a.evidence.rglob("ica_science_validated.json"))
    if len(raws)!=3 or len(val)!=3:
        raise SystemExit(f"Expected three real science and validated records: {len(raws)}/{len(val)}")
    original={x["subject"]:x for p in raws for x in [json.loads(p.read_text())]}
    summaries={x["subject"]:x for p in val for x in [json.loads(p.read_text())]}
    if set(original)!={11,15,21} or set(summaries)!={11,15,21}:
        raise SystemExit("Wrong development cohort")
    rows=[examine(original[k],summaries[k]) for k in sorted(original)]
    out={"study":"ICA weak-component topography audit, read-only, pre-examined development cohort",
         "scienceRun":43184910782,"v1MultimetricRun":43193549225,
         "newEEGFits":False,"newParticipants":False,"original20SubjectVerdict":"PARTIAL",
         "scientificPromotionAuthorized":False,"results":rows}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,indent=2))
    for z in rows:
        for cand in z["candidates"]:
            train=cand["starts"][:2]
            bad=sum(s["belowPoint90"] for s in train)
            affected=sum(s["lowCorrelationRemovalDiscordant"] for s in train)
            print("PR51_COMPONENT_OUTLIER",json.dumps({
                "subject":z["subject"],"candidateId":cand["candidateId"],
                "weakMatchedComponentsTraining":bad,
                "weakComponentsWithRemovalDiscordanceTraining":affected,
                "weakBothKeptTraining":sum(s["lowCorrelationBothKept"] for s in train),
                "weakBothRemovedTraining":sum(s["lowCorrelationBothRemoved"] for s in train),
                "weakICLabelClassDiscordantTraining":sum(s["lowCorrelationClassDiscordant"] for s in train),
                "lowestCorrelationTraining":min(s["minimumCorrelation"] for s in train),
                "allSeedWeakComponents":sum(s["belowPoint90"] for s in cand["starts"])}))
    print("PR51_COMPONENT_OUTLIER_DECISION",json.dumps({
        "multimetricV1Verdict":"ABSTAIN","candidatePromotionAuthorized":False,
        "originalPR51Verdict":"PARTIAL"}))


if __name__=="__main__":
    main()
