#!/usr/bin/env python3
"""Compare exactly observed ERP CORE ICA maps across six starts.

This is development-only, not a new scientific estimator. Hungarian
one-to-one component matching is sign invariant; component numbers alone
cannot be compared across independent ICA fits. Jaccard of removed sets
is *undefined* when both sets are empty, not artificially perfect.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment
from validate_ica_unseen_seeds import grid, inspect

SUBJECTS = {11,15,21}


def alignment(reference, target):
    if not (isinstance(reference,dict) and isinstance(target,dict)
            and reference.get("available") and target.get("available")):
        return {"available":False,"reason":"missing-ICA-fingerprint"}
    ch0=reference.get("channels")
    ch1=target.get("channels")
    if not (isinstance(ch0,list) and isinstance(ch1,list)
            and len(ch0)>1 and len(ch0)==len(ch1)
            and len(set(ch0))==len(ch0)
            and set(ch0)==set(ch1)):
        return {"available":False,"reason":"different-or-invalid-channel-subsets"}
    aa=np.asarray(reference.get("maps"),dtype=float)
    bb=np.asarray(target.get("maps"),dtype=float)
    if (aa.ndim!=2 or bb.ndim!=2 or aa.shape[0]!=len(ch0)
            or bb.shape[0]!=len(ch1) or aa.shape[1]!=bb.shape[1]
            or aa.shape[1]<1 or not np.isfinite(aa).all()
            or not np.isfinite(bb).all()):
        return {"available":False,"reason":"nonfinite-or-incompatible-ICA-dimensions"}
    bb=bb[[ch1.index(c) for c in ch0],:]
    # Pearson correlation of spatial topographies, not time courses.
    aa=aa-aa.mean(axis=0,keepdims=True)
    bb=bb-bb.mean(axis=0,keepdims=True)
    na=np.linalg.norm(aa,axis=0)
    nb=np.linalg.norm(bb,axis=0)
    if np.any(na<=1e-12) or np.any(nb<=1e-12):
        return {"available":False,"reason":"constant-or-degenerate-topography"}
    similarity=np.abs((aa/na).T @ (bb/nb))
    rows,cols=linear_sum_assignment(-similarity)
    matched=similarity[rows,cols]
    r0=reference.get("removedIndices")
    r1=target.get("removedIndices")
    if not (isinstance(r0,list) and isinstance(r1,list)):
        return {"available":False,"reason":"missing-component-removal-identities"}
    if (any(not isinstance(v,int) or isinstance(v,bool)
            or v<1 or v>aa.shape[1] for v in r0+r1)
            or len(set(r0))!=len(r0) or len(set(r1))!=len(r1)):
        return {"available":False,"reason":"invalid-component-removal-identities"}
    a_removed={v-1 for v in r0}
    aligned_removed={int(row) for row,col in zip(rows,cols) if int(col)+1 in r1}
    union=a_removed|aligned_removed
    removed_jaccard=len(a_removed & aligned_removed)/len(union) if union else None
    return {"available":True,"reason":"one-to-one-absolute-topographic-correlation",
            "nComponents":int(aa.shape[1]),"nChannels":int(aa.shape[0]),
            "medianMatchedAbsTopologyCorrelation":float(np.median(matched)),
            "minMatchedAbsTopologyCorrelation":float(np.min(matched)),
            "removedCountBaseline":len(a_removed),
            "removedCountOtherStart":len(r1),
            "alignedRemovedJaccard":removed_jaccard,
            "zeroRemovalBothSeeds":not union,
            "matchedColumns":[[int(a)+1,int(b)+1] for a,b in zip(rows,cols)]}


def analyze(subject_result):
    if not subject_result.get("completed") or subject_result.get("errors"):
        raise ValueError("Incomplete real-data participant")
    s=subject_result["subject"]
    if s not in SUBJECTS:
        raise ValueError("Unexpected development participant")
    ic=subject_result["icaSeedCheck"]
    ids=ic["ids"]
    if len(ic.get("recommended",[]))!=6 or len(ic.get("objective",[]))!=len(ids):
        raise ValueError("Not all six planned ICA starts recorded")
    grids=grid(ic["icaComponentFingerprints"],len(ids),"icaComponentFingerprints")
    if not all(len(x)==6 for x in grids):
        raise ValueError("ICA maps do not align to six starts")
    vetted=inspect(subject_result)
    chosen=[("baseline",ic["recommended"][0]),
            ("minimax",vetted["shadowUnseen"]["candidateId"])]
    results=[]
    for role,ident in chosen:
        if ident is None:
            results.append({"role":role,"id":None,"abstain":True})
            continue
        if ident not in ids:
            raise ValueError(f"Candidate {ident} not in recorded ICA set")
        at=ids.index(ident)
        ref=grids[at][0]
        checks=[]
        for seed in range(1,6):
            checks.append({"seed":seed,**alignment(ref,grids[at][seed])})
        if any(not t["available"] for t in checks):
            # An unavailable fit cannot be counted as a stable one.
            unavailable=[x["seed"] for x in checks if not x["available"]]
        else:
            unavailable=[]
        results.append({"role":role,"id":ident,
            "seed0FingerprintAvailable":bool(ref.get("available")),
            "seed0RemovedCount":ref.get("removedCount"),
            "matchedStartCount":5-len(unavailable),
            "unmatchedStartIndices":unavailable,"startComparisons":checks,
            "topographicComparisonIsNotSourceIdentity":True,
            "fixedIDByConstructionIsNotValidation":role=="minimax"})
    cross=None
    if len(results)==2 and not any(x.get("abstain") for x in results):
        x=grids[ids.index(results[0]["id"])][0]
        y=grids[ids.index(results[1]["id"])][0]
        cross=alignment(x,y)
    return {"subject":s,"sourceRun":subject_result.get("workflowRun"),
        "seedIndices":[0,1,2,3,4,5],"referenceSeed":0,
        "independentAuditSeeds":[3,4,5],"candidates":results,
        "sameSeedCrossRecipeAlignment":cross,
        "newParticipantsEvaluated":False,"scientificPromotionAuthorized":False}


def run(folder):
    files=sorted(folder.rglob("result.json"))
    if len(files)!=3:
        raise ValueError(f"Expected 3 complete real-EEG records; observed {len(files)}")
    rows=[analyze(json.loads(p.read_text())) for p in files]
    if {r["subject"] for r in rows}!=SUBJECTS:
        raise ValueError("Duplicate or unexpected participant")
    rows.sort(key=lambda r:r["subject"])
    return {"study":"PR51 six-start ICA component topography mechanistic development",
        "source":"ERP CORE P3, original 300 s per subject",
        "originalFrozenVerdict":"PARTIAL","validOriginalPairs":15,
        "method":"one-to-one maximization of absolute Pearson scalp-map correlations",
        "notTested":["cross-recording source identity","full ERP waveform shape",
                     "new participant generalization","improvement of P300 precision"],
        "rows":rows,"promotionAuthorized":False,"freeze14NewSubjectsAuthorized":False}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--evidence",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    result=run(a.evidence)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,allow_nan=False))
    for subject in result["rows"]:
        for record in subject["candidates"]:
            checks=record.get("startComparisons",[])
            print("PR51_ICA_TOPOLOGY",json.dumps({
                "subject":subject["subject"],"role":record["role"],
                "candidate":record["id"],"removedSeed0":record.get("seed0RemovedCount"),
                "matchedStarts":record.get("matchedStartCount"),
                "unseenTopographicMed":[c.get("medianMatchedAbsTopologyCorrelation")
                    for c in checks if c["seed"] in (3,4,5)],
                "unseenRemovalJaccard":[c.get("alignedRemovedJaccard")
                    for c in checks if c["seed"] in (3,4,5)]}))
if __name__=="__main__":
    main()
