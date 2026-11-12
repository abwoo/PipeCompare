#!/usr/bin/env python3
"""Check stage-specific ICA-component Pz signatures, without changing Rank.

The pre-IC-removal projected component amplitudes are NOT final ERP
counterfactual effects. Original candidate scores/ranks must pass the existing
separate ICA science identity/parity validator first.
"""
import argparse
import json
import math
from pathlib import Path
import numpy as np

from validate_ica_science import FROZEN_DEVELOPMENT_IDS
from validate_ica_unseen_seeds import grid

def number(x):
    return isinstance(x,(int,float)) and not isinstance(x,bool) and math.isfinite(x)

def sequence(v):
    if isinstance(v,list):return v
    return [v] if v is not None else []

def contribution(stage,final):
    errors=[]
    if not isinstance(stage,dict) or stage.get("status")!="complete":
        raise ValueError("Missing/failed pre-removal signal observer")
    ids=sequence(stage.get("eventIds"))
    cond=sequence(stage.get("conditionIndex"))
    nc=stage.get("nComponents"); n=stage.get("nEvents")
    try:
        mat=np.asarray(stage.get("componentUv",[]),dtype=float)
        observed=np.asarray(sequence(stage.get("observedPzUv")),dtype=float)
        residual=np.asarray(sequence(stage.get("residualPzUv")),dtype=float)
    except (ValueError,TypeError) as exc:
        raise ValueError("Non-numerical component signal") from exc
    if (not isinstance(n,int) or not isinstance(nc,int) or n<2 or nc<1
            or mat.shape!=(nc,n) or observed.shape!=(n,) or residual.shape!=(n,)
            or len(ids)!=n or len(cond)!=n or len(set(ids))!=n
            or any(c not in (1,2) for c in cond)
            or not np.isfinite(mat).all() or not np.isfinite(observed).all()
            or not np.isfinite(residual).all()):
        raise ValueError("Misaligned per-event component Pz signature")
    if not np.allclose(mat.sum(axis=0)+residual,observed,atol=1e-6,rtol=1e-7):
        raise ValueError("Sensor reconstruction and residual disagree")
    if not np.isclose(np.max(np.abs(residual)),stage.get("maxAbsResidualUv"),atol=1e-7):
        raise ValueError("Wrong stored sensor residual")
    rank=stage.get("rank"); nchan=stage.get("nICAChannels")
    if not isinstance(rank,int) or not isinstance(nchan,int) or rank>nchan or rank>nc:
        raise ValueError("Invalid ICA rank and dimension")
    if rank==nchan and np.max(np.abs(residual))>1e-4:
        raise ValueError("Full-rank decomposition unexpectedly misses Pz signal")
    conditions=sequence(final.get("conditions"))
    if len(conditions)!=2: raise ValueError("Missing final ERP condition evidence")
    by_id={int(event):k for k,event in enumerate(ids)}
    estimates=[]
    for condition_index,cs in enumerate(conditions,1):
        ref=sequence(cs.get("referenceIds"))
        final_ids=sequence(cs.get("retainedIds"))
        weights=sequence(cs.get("trialWeights"))
        if len(ref)!=len(weights) or len(final_ids)==0 or len(set(final_ids))!=len(final_ids):
            raise ValueError("Final original trial identities/weights invalid")
        weight_of={int(event):weight for event,weight in zip(ref,weights)}
        if len(weight_of)!=len(ref): raise ValueError("Duplicate original trial")
        if any(event not in by_id or cond[by_id[event]]!=condition_index
               or event not in weight_of or not number(weight_of[event])
               or weight_of[event]<=0 for event in final_ids):
            raise ValueError("Final retained event absent from pre-removal signature")
        idx=[by_id[event] for event in final_ids]
        weight=np.asarray([weight_of[event] for event in final_ids],dtype=float)
        wsum=weight.sum()
        projected=(mat[:,idx]@weight)/wsum
        sensor=float(observed[idx]@weight/wsum)
        residue=float(residual[idx]@weight/wsum)
        if not np.isclose(projected.sum()+residue,sensor,atol=1e-6,rtol=1e-7):
            raise ValueError("Weighted original trial reconstruction failed")
        estimates.append({"conditionIndex":condition_index,
                          "retainedFinalEvents":len(idx),
                          "preRemovalComponentPzMeanUv":projected.tolist(),
                          "preRemovalSensorPzMeanUv":sensor,
                          "preRemovalResidualPzMeanUv":residue})
    return {
        "componentCount":nc,"stageEvents":n,"icaRank":rank,"icaChannelCount":nchan,
        "maxAbsStageReconstructionResidualUv":float(np.max(np.abs(residual))),
        "stageComponentPzDifferenceUv":(
            np.asarray(estimates[0]["preRemovalComponentPzMeanUv"])-
            np.asarray(estimates[1]["preRemovalComponentPzMeanUv"])).tolist(),
        "stageSensorPzDifferenceUv":estimates[0]["preRemovalSensorPzMeanUv"]-
                                    estimates[1]["preRemovalSensorPzMeanUv"],
        "stageResidualPzDifferenceUv":estimates[0]["preRemovalResidualPzMeanUv"]-
                                      estimates[1]["preRemovalResidualPzMeanUv"],
        "joinedFinalTrialIdentities":True,
        "metricStage":"before original IC removal; fixed observed final trial identity and weights",
        "equalsFinalP300Estimand":False,
    }

def evaluate(source,old):
    sid=source.get("subject")
    if sid not in FROZEN_DEVELOPMENT_IDS or old.get("subject")!=sid:
        raise ValueError("Unexpected/mismatched research participant")
    if not source.get("completed") or source.get("errors") or old.get("identityErrors") or old.get("scienceFailures"):
        raise ValueError("Original source parity/science capture failure")
    c=source["icaSeedCheck"]
    ids=c["ids"]
    if c["recommended"][0]!=FROZEN_DEVELOPMENT_IDS[sid][0]:
        raise ValueError("Unexpected old Standard candidate")
    m=grid(c["candidateMeasurements"],len(ids),"candidateMeasurements")
    results=[]
    for candidate in FROZEN_DEVELOPMENT_IDS[sid]:
        if candidate not in ids:raise ValueError("Expected archived candidate absent")
        i=ids.index(candidate);starts=[]
        for seed in range(6):
            z=m[i][seed]
            if not isinstance(z,dict) or z.get("candidateId")!=candidate or z.get("seed")!=seed:
                raise ValueError("Candidate/start provenance mismatch")
            science=z.get("science")
            if not isinstance(science,dict) or science.get("status")!="complete":
                raise ValueError("Missing completed original science")
            captures=sequence(science.get("preRemovalContributions"))
            if len(captures)!=1:
                raise ValueError(f"Expected exactly one pre-removal observer: {sid}/{candidate}/{seed}, found {len(captures)}")
            measure=contribution(captures[0],science)
            starts.append({"seed":seed,**measure})
        results.append({"candidateId":candidate,"startReports":starts})
    return {"subject":sid,"candidates":results,"newParticipants":False,
            "candidateSelectionModified":False,"promotionAuthorized":False}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--new",type=Path,required=True)
    p.add_argument("--old",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    new=json.loads(a.new.read_text());old=json.loads(a.old.read_text())
    output=evaluate(new,old)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(output,indent=2,allow_nan=False))
    for c in output["candidates"]:
        print("PR51_ICA_PRE_REMOVAL",json.dumps({
            "subject":output["subject"],"candidateId":c["candidateId"],
            "seedCount":len(c["startReports"]),
            "preStageResidualMaxUv":max(x["maxAbsStageReconstructionResidualUv"]
                                       for x in c["startReports"]),
            "stageP3RangeUv":max(x["stageSensorPzDifferenceUv"] for x in c["startReports"])-
                             min(x["stageSensorPzDifferenceUv"] for x in c["startReports"]),
            "productionPromotionAuthorized":False}))
    print("PR51_ICA_PRE_REMOVAL_SCOPE",json.dumps({
        "originalPR51Verdict":"PARTIAL","newParticipants":False,
        "productionPromotionAuthorized":False,
        "fullERPComponentCausalityProven":False}))

if __name__=="__main__":
    main()
