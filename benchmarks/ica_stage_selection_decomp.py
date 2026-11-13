#!/usr/bin/env python3
"""Fixed-reference arithmetic decomposition of original pre-removal Pz events.

Development-only, does not refit or alter ICA/ERP trial decisions. Source:
docs/PR51_ICA_STAGE_SELECTION_PROTOCOL_20261113.md
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np

from ica_stage_impact import fixed_trial_comparison, stage_contributions
from validate_ica_science import FROZEN_DEVELOPMENT_IDS
from validate_ica_unseen_seeds import grid

SOURCE_RUN=43397978997
EPS=1e-6


def serial(x):
    if x is None:return []
    return x if isinstance(x,list) else [x]


def finite(x):
    return isinstance(x,(int,float)) and not isinstance(x,bool) and math.isfinite(x)


def stage_map(measure):
    science=measure.get("science",{})
    captures=serial(science.get("preRemovalContributions"))
    if len(captures)!=1 or not isinstance(captures[0],dict):
        raise ValueError("Missing/duplicate source pre-removal stage capture")
    z=captures[0]
    if z.get("status")!="complete":
        raise ValueError("Actual source pre-removal stage failed")
    ids=serial(z.get("eventIds"))
    cond=serial(z.get("conditionIndex"))
    volt=serial(z.get("observedPzUv"))
    n=z.get("nEvents")
    if (not isinstance(n,int) or n<2 or
            not (len(ids)==len(cond)==len(volt)==n) or
            len(set(ids))!=n or any(type(i) is not int for i in ids) or
            any(c not in (1,2) for c in cond) or
            any(not finite(y) for y in volt)):
        raise ValueError("Source stage sensor and event indices misaligned")
    if set(cond)!={1,2}:
        raise ValueError("Source stage lacks both ERP event conditions")
    return {i:(c,float(v)) for i,c,v in zip(ids,cond,volt)}


def final_trial_reference(measure):
    science=measure.get("science",{})
    conditions=serial(science.get("conditions"))
    if len(conditions)!=2:
        raise ValueError("Missing final two-condition trial/weight science")
    choices=[]
    for ci,condition in enumerate(conditions,1):
        ids=serial(condition.get("referenceIds"))
        weights=serial(condition.get("trialWeights"))
        kept=serial(condition.get("retainedIds"))
        if (len(ids)!=len(weights) or len(set(ids))!=len(ids) or
                len(set(kept))!=len(kept) or not kept):
            raise ValueError("Original reference trial/weight cardinality invalid")
        raw=dict(zip(ids,weights))
        selected={}
        for ident in kept:
            if ident not in raw or type(ident) is not int:
                raise ValueError("Final original trial missing from reference")
            w=raw[ident]
            if not finite(w) or w<=0:
                raise ValueError("Final source trial weight invalid")
            selected[ident]=float(w)
        choices.append(selected)
    return choices


def weighted(source,reference):
    results=[]
    for condition,selected in enumerate(reference,1):
        ids=list(selected)
        if not ids:
            raise ValueError("Empty final trial set")
        missing=[i for i in ids if i not in source or source[i][0]!=condition]
        if missing:
            raise ValueError("Original event missing from stage sensor or wrong condition: "+str(missing[:3]))
        x=np.asarray([source[i][1] for i in ids],dtype=float)
        w=np.asarray([selected[i] for i in ids],dtype=float)
        if not np.isfinite(x).all() or not np.isfinite(w).all() or np.any(w<=0):
            raise ValueError("Nonfinite or nonpositive stage data/weights")
        results.append(float(np.dot(w,x)/np.sum(w)))
    return results


def difference(v):
    if len(v)!=2:
        raise ValueError("Two-condition ERP only")
    return v[0]-v[1]


def compare(start0,start1,old_validated0,old_validated1,original_erp0,original_erp1):
    src0=stage_map(start0)
    src1=stage_map(start1)
    w0=final_trial_reference(start0)
    w1=final_trial_reference(start1)
    m0=weighted(src0,w0)
    ms_same0=weighted(src1,w0)
    ms_own=weighted(src1,w1)
    # For an explicitly symmetric provenance check, prove that source0 also
    # contains every later-retained event even though seed0 is the reference.
    weighted(src0,w1)
    t=[b-a for a,b in zip(m0,ms_own)]
    s=[b-a for a,b in zip(m0,ms_same0)]
    q=[b-a for a,b in zip(ms_same0,ms_own)]
    if not np.allclose(t,np.asarray(s)+q,rtol=1e-7,atol=EPS):
        raise ValueError("Observed stage selection/signal decomposition nonconservative")
    for measured,verified in [(difference(m0),old_validated0),(difference(ms_own),old_validated1)]:
        if not math.isclose(measured,verified,abs_tol=EPS,rel_tol=1e-7):
            raise ValueError("Event-level original stage source differs from verified stage Pz difference")
    status=fixed_trial_comparison(start0["science"],start1["science"])
    final=float(original_erp1-original_erp0)
    return {
        "observedTotalStagePzDeltaUv":difference(t),
        "fixedSeed0TrialsSignalDeltaUv":difference(s),
        "changedTrialSelectionOrWeightsDeltaUv":difference(q),
        "targetConditionTotalUv":t[0],"standardConditionTotalUv":t[1],
        "targetFixedSignalUv":s[0],"standardFixedSignalUv":s[1],
        "targetSelectionWeightUv":q[0],"standardSelectionWeightUv":q[1],
        "seed0StagePzUv":difference(m0),
        "otherStagePzUv":difference(ms_own),
        "originalFinalGainCorrectedP300DeltaUv":final,
        "originalFinalEndpointIsNotStageMeasurement":True,
        "subjectLevelCausalInferenceAuthorized":False,
        "referenceWeightsAreSeed0":True,
        "allPhysicalEventsPresentInBothStages":True,
        "arithmeticallyConserved":True,
        **status,
    }


def evaluate(real,science,record):
    sid=real.get("subject")
    if (sid not in FROZEN_DEVELOPMENT_IDS or sid!=science.get("subject") or
            sid!=record.get("subject") or not real.get("completed") or real.get("errors") or
            science.get("identityErrors") or science.get("scienceFailures") or
            not science.get("unalteredFrozenRankAndEstimates") or
            not science.get("allPredeclaredFeasibleScienceComplete")):
        raise ValueError("Original ERP science parity invalid")
    source=real["icaSeedCheck"]
    ids=source["ids"]
    if source["recommended"][0]!=FROZEN_DEVELOPMENT_IDS[sid][0]:
        raise ValueError("Original Standard selected candidate changed")
    measure=grid(source["candidateMeasurements"],len(ids),"candidateMeasurements")
    verified={z["candidateId"]:z for z in record["candidates"]}
    if set(verified)!=set(FROZEN_DEVELOPMENT_IDS[sid]):
        raise ValueError("Original pre-removal candidate roster incomplete")
    observations=[]
    for candidate in FROZEN_DEVELOPMENT_IDS[sid]:
        n=ids.index(candidate)
        row=[measure[n][s] for s in range(6)]
        external=sorted(verified[candidate]["startReports"],key=lambda z:z["seed"])
        if len(external)!=6 or [z["seed"] for z in external]!=list(range(6)):
            raise ValueError("Verified 6 independent ICA seed records incomplete")
        for seed in range(6):
            if row[seed].get("candidateId")!=candidate or row[seed].get("seed")!=seed:
                raise ValueError("Candidate seed identities changed")
            stage_contributions(external[seed])
        starts=[]
        for seed in range(1,6):
            x=compare(row[0],row[seed],
                      external[0]["stageSensorPzDifferenceUv"],
                      external[seed]["stageSensorPzDifferenceUv"],
                      row[0]["correctedDifferenceEstimate"],
                      row[seed]["correctedDifferenceEstimate"])
            x["seed"]=seed
            starts.append(x)
        observations.append({"candidateId":candidate,"baseline":candidate==FROZEN_DEVELOPMENT_IDS[sid][0],
                             "seed0FixedReference":True,"starts":starts})
    return {"subject":sid,"candidates":observations,
            "newParticipants":False,"icaAlgorithmModified":False,
            "causalEffectClaimed":False}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--evidence",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    def source(filename):
        files=sorted(args.evidence.rglob(filename))
        if len(files)!=3:
            raise ValueError("Three archived subject records required for "+filename)
        records=[json.loads(file.read_text()) for file in files]
        by_sid={x["subject"]:x for x in records}
        if set(by_sid)!={11,15,21}:
            raise ValueError("Missing, duplicated or unexpected science subjects")
        return by_sid
    original=source("result.json")
    science=source("ica_science_validated.json")
    stage=source("ica_contribution_validation.json")
    rows=[evaluate(original[s],science[s],stage[s]) for s in sorted(original)]
    result={"analysis":"development-only fixed seed0 trial/weight arithmetic separation of pre-removal Pz",
            "realEEGSourceRun":SOURCE_RUN,"original20SubjectVerdict":"PARTIAL",
            "original20SubjectValidPairs":15,"frozenHoldoutUsed":False,
            "ICARecipeChanged":False,"pr51MergeAuthorized":False,
            "candidatePromotionAuthorized":False,"results":rows}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,allow_nan=False,indent=2))
    for subj in rows:
        for c in subj["candidates"]:
            for z in c["starts"]:
                print("PR51_SELECTION_DECOMPOSITION",json.dumps({
                    "subject":subj["subject"],"candidate":c["candidateId"],"seed":z["seed"],
                    "totalStageUv":round(z["observedTotalStagePzDeltaUv"],9),
                    "fixedTrialStageSignalUv":round(z["fixedSeed0TrialsSignalDeltaUv"],9),
                    "trialSelectionWeightUv":round(z["changedTrialSelectionOrWeightsDeltaUv"],9),
                    "finalP300DriftUv":round(z["originalFinalGainCorrectedP300DeltaUv"],9),
                    "selectionWeightProvenanceChanged":z["startDifferenceConfoundedByTrialSelectionOrWeights"],
                    "conserved":z["arithmeticallyConserved"]
                }))
    print("PR51_SELECTION_DECOMPOSITION_DECISION",json.dumps({
        "pr51MergeAuthorized":False,"newParticipants":False,"causalityClaimed":False,
        "originalPR51Verdict":"PARTIAL"}))


if __name__=="__main__":
    main()
