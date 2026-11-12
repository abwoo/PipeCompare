#!/usr/bin/env python3
"""Post-capture observational ICA topology and ERP stability diagnostics.

No pipeline tuning, no new ranking, no scientific promotion. Match ICs by
absolute scalp-map correlation with Hungarian assignment, robust to
component permutation and sign; retain mismatches as negative evidence.
"""
import argparse
import json
import math
from pathlib import Path
import numpy as np
from scipy.optimize import linear_sum_assignment
from validate_ica_unseen_seeds import grid

RECIPES={11:{"original":12,"minimax":25,"stabilityFirst":25},
         15:{"original":5,"minimax":2,"stabilityFirst":6},
         21:{"original":25,"minimax":4,"stabilityFirst":4}}


def one_ica(s):
    ica=s.get("ica",[])
    if isinstance(ica,dict):ica=[ica]
    if not isinstance(ica,list) or len(ica)!=1:
        raise ValueError("Expected one ICA decomposition per Standard recipe")
    return ica[0]


def scalar_ids(x):
    if isinstance(x,list):return x
    if isinstance(x,(int,float)) and math.isfinite(x):return [int(x)]
    return []


def ic_match(a,b):
    x=one_ica(a);y=one_ica(b)
    ci=x["channelIndices"];cj=y["channelIndices"]
    if ci!=cj:raise ValueError("Different ICA channel sets; not comparable")
    A=np.asarray(x["mixing"],dtype=float)
    B=np.asarray(y["mixing"],dtype=float)
    if A.ndim!=2 or B.ndim!=2 or A.shape[0]!=B.shape[0] or A.shape[0]<3:
        raise ValueError("Insufficient common spatial channels")
    # Normalize maps across channels; sign and column scale are arbitrary.
    A=A-A.mean(axis=0,keepdims=True)
    B=B-B.mean(axis=0,keepdims=True)
    na=np.linalg.norm(A,axis=0);nb=np.linalg.norm(B,axis=0)
    if np.any(na<=1e-12) or np.any(nb<=1e-12):
        raise ValueError("Zero-norm IC map")
    C=np.abs((A/na).T@(B/nb))
    r,c=linear_sum_assignment(-C)
    matched=C[r,c]
    original_removed=set(scalar_ids(x.get("removedIndices",[])))
    alternate_removed=set(scalar_ids(y.get("removedIndices",[])))
    # MATLAB component indexing is one-based.
    removed_map={int(i+1):int(j+1) for i,j in zip(r,c)}
    corresponding={removed_map[i] for i in original_removed if i in removed_map}
    both=corresponding|alternate_removed
    jaccard=(len(corresponding&alternate_removed)/len(both)) if both else 1.
    return {"nComponentsA":A.shape[1],"nComponentsB":B.shape[1],
            "matchedCount":len(matched),"medianAbsTopoCorrelation":float(np.median(matched)),
            "minAbsTopoCorrelation":float(np.min(matched)),
            "removedA":sorted(original_removed),"removedB":sorted(alternate_removed),
            "removedMatchJaccard":float(jaccard),
            "topoMatchingIsSignAndPermutationInvariant":True}


def waveform_compare(a,b):
    X=a["waveform"];Y=b["waveform"]
    if X["timesMs"]!=Y["timesMs"]:
        raise ValueError("Different ERP time axes")
    if X["gainApplied"] or Y["gainApplied"]:
        raise ValueError("Unexpected already-gain-corrected waveform")
    x=np.asarray(X["difference"],dtype=float);y=np.asarray(Y["difference"],dtype=float)
    if x.shape!=y.shape or x.size<10 or not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError("ERP waveforms absent/nonfinite")
    delta=y-x
    centeredX=x-x.mean();centeredY=y-y.mean()
    norm=np.linalg.norm(centeredX)*np.linalg.norm(centeredY)
    corr=float(centeredX@centeredY/norm) if norm>1e-12 else None
    return {"samples":x.size,"fullEpochRmsDifferenceUv":float(np.sqrt(np.mean(delta**2))),
            "fullEpochPearson":corr,"gainCorrected":False}


def examine(result,science_report):
    subj=int(result["subject"])
    c=result["icaSeedCheck"];ids=c["ids"]
    ms=grid(c["candidateMeasurements"],len(ids),"candidateMeasurements")
    bootstrap_by_pair={(x["candidateId"],x["seed"]):x.get("bootstrap")
                       for x in science_report["rows"] if x["valid"]}
    rows=[]
    for name,ident in RECIPES[subj].items():
        if ident not in ids:raise ValueError("Named recipe not present in ICA candidate roster")
        values=ms[ids.index(ident)]
        first=values[0].get("science")
        if not isinstance(first,dict) or not first.get("available"):
            raise ValueError("Baseline seed ICA science unavailable")
        seed_zero_boot=bootstrap_by_pair.get((ident,0))
        if not isinstance(seed_zero_boot,dict):
            raise ValueError("Seed-zero bootstrap SME missing")
        seeds=[]
        for seed in (3,4,5):
            z=values[seed].get("science")
            if not isinstance(z,dict) or not z.get("available"):
                raise ValueError(f"Missing independent seed {seed} scientific snapshot")
            try:
                ic=ic_match(first,z)
                wp=waveform_compare(first,z)
                problems=[]
            except (ValueError,KeyError,TypeError) as exc:
                ic=None;wp=None;problems=[str(exc)]
            bootstrap=bootstrap_by_pair.get((ident,seed))
            if not isinstance(bootstrap,dict) or not isinstance(bootstrap.get("se"),(int,float)):
                problems.append("independent-seed-bootstrap-unavailable")
            seeds.append({"seed":seed,"topology":ic,"waveform":wp,
                          "bootstrapSMEuv":bootstrap.get("se") if bootstrap else None,
                          "bootstrapEffectiveCounts":bootstrap.get("sampleCounts") if bootstrap else None,
                          "unavailable":problems})
        rows.append({"recipe":name,"candidateId":ident,
                     "baseSeed":0,"seedZeroBootstrapSMEuv":seed_zero_boot["se"],
                     "independentSeeds":seeds})
    return {"subject":subj,"candidateStudies":rows,
            "selectionAlgorithmChanged":False,
            "newParticipantsUsed":False}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--evidence",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    files=sorted(a.evidence.rglob("result.json"))
    if len(files)!=3:raise SystemExit(f"Expected three real MATLAB artifacts, got {len(files)}")
    rows=[]
    for f in files:
        result=json.loads(f.read_text())
        report=json.loads((f.parent/"science_record_audit.json").read_text())
        if report["subject"]!=result["subject"] or report["problems"]:
            raise ValueError("Scientific recorder not validated")
        rows.append(examine(result,report))
    if sorted(x["subject"] for x in rows)!=sorted(RECIPES):raise ValueError("Unexpected development cohort")
    rows.sort(key=lambda r:r["subject"])
    out={"analysis":"descriptive post-capture component and waveform alignment only",
         "originalPR51Verdict":"PARTIAL","recorderRun":43189639057,
         "trainingOrTuningOnNewParticipants":False,
         "algorithmPromotionAuthorized":False,"rows":rows}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,indent=2))
    for row in rows:
        for rec in row["candidateStudies"]:
            for s in rec["independentSeeds"]:
                print("PR51_COMPONENT_WAVEFORM",json.dumps({
                    "subject":row["subject"],"recipe":rec["recipe"],
                    "seed":s["seed"],
                    "medianAbsTopoCorrelation":s["topology"]["medianAbsTopoCorrelation"] if s["topology"] else None,
                    "removedJaccard":s["topology"]["removedMatchJaccard"] if s["topology"] else None,
                    "fullEpochRmsDifference":s["waveform"]["fullEpochRmsDifferenceUv"] if s["waveform"] else None,
                    "seedZeroBootstrapSMEuv":rec["seedZeroBootstrapSMEuv"],
                    "unseenBootstrapSMEuv":s["bootstrapSMEuv"],
                    "unavailable":s["unavailable"]}))


if __name__=="__main__":
    main()
