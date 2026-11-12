#!/usr/bin/env python3
"""Fail-closed integrity and independent 1000-draw trial bootstrap of science records.

This reads already-computed experimental MATLAB results. It never refits ICA,
relabels source events, modifies a pipeline, or selects a new recommendation.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from validate_ica_unseen_seeds import grid, valid

N_SEEDS=6
DRAWS=1000
SEED=610927


def inspect_one(measure,ident,seed):
    failures=[]
    if not isinstance(measure,dict) or measure.get("seed")!=seed or measure.get("candidateId")!=ident:
        return {"candidateId":ident,"seed":seed,"valid":False,"issues":["measure-alignment"]}
    science=measure.get("science")
    if not isinstance(science,dict) or not science.get("available"):
        return {"candidateId":ident,"seed":seed,"valid":False,
                "issues":["science-snapshot-missing:"+str(science.get("reason") if isinstance(science,dict) else "absent")]}
    if science.get("seed")!=seed or science.get("candidateId")!=ident:
        failures.append("science-snapshot-alignment")
    wave=science.get("waveform",{})
    conditions=wave.get("condition",[])
    if not isinstance(conditions,list) or len(conditions)!=2:
        failures.append("waveform-conditions-missing")
        return {"candidateId":ident,"seed":seed,"valid":False,"issues":failures}
    times=np.asarray(wave.get("timesMs",[]),dtype=float)
    diff=np.asarray(wave.get("difference",[]),dtype=float)
    if len(times)<10 or len(diff)!=len(times) or not np.isfinite(times).all() or not np.isfinite(diff).all():
        failures.append("full-waveform-shape-invalid")
    ica=science.get("ica",[])
    # MATLAB jsonencode emits scalar 1x1 structs as objects, not arrays.
    if isinstance(ica,dict): ica=[ica]
    if not isinstance(ica,list) or not ica:
        failures.append("ica-decomposition-missing")
    else:
        for comp in ica:
            W=np.asarray(comp.get("weights",[]),dtype=float)
            sphere=np.asarray(comp.get("sphere",[]),dtype=float)
            mix=np.asarray(comp.get("mixing",[]),dtype=float)
            if W.ndim!=2 or sphere.ndim!=2 or mix.ndim!=2 or W.shape[0]!=mix.shape[1] or not np.isfinite(mix).all():
                failures.append("ica-topography-shape-invalid")
            ids=comp.get("removedIndices",[])
            if isinstance(ids,(int,float)) and math.isfinite(ids):ids=[ids]
            if not isinstance(ids,list) or any(not isinstance(k,int) or k<1 or k>W.shape[0] for k in ids):
                failures.append("removed-component-id-invalid")
    per_condition=[]
    for idx,cond in enumerate(conditions):
        ids=cond.get("sourceIds",[])
        kept=cond.get("keptIds",[])
        wt=cond.get("weights",[])
        scores=cond.get("scores",[])
        waveform=cond.get("waveform",[])
        if not all(isinstance(x,list) for x in (ids,kept,wt,scores,waveform)):
            failures.append(f"trial-array-missing-{idx}")
            continue
        if not (len(ids)==len(wt)==len(scores) and len(ids)>0 and len(waveform)==len(times)):
            failures.append(f"trial-array-misaligned-{idx}")
            continue
        if not set(kept).issubset(set(ids)) or len(set(ids))!=len(ids):
            failures.append(f"epoch-identity-invalid-{idx}")
        # MATLAB JSON encodes NaN scores and absent epoch weights as null.
        # Explicit float coercion maps null -> NaN; NumPy's implicit
        # object-array coercion cannot safely support isfinite().
        scalar_scores=np.asarray(scores,dtype=float)
        weights=np.asarray(wt,dtype=float)
        waveform_values=np.asarray(waveform,dtype=float)
        good=np.isfinite(scalar_scores)&np.isfinite(weights)&(weights>0)
        if int(np.count_nonzero(good))<2:
            failures.append(f"insufficient-bootstrap-scores-{idx}")
            continue
        s=scalar_scores[good];w=weights[good]
        per_condition.append((s,w))
        if not np.isfinite(waveform_values).all():failures.append(f"waveform-nonfinite-{idx}")
    boot=None
    gain=measure.get("gain")
    if len(per_condition)==2 and valid(gain) and gain>0:
        rng=np.random.RandomState(SEED+seed*10000+ident)
        estimates=[]
        for scores,weight in per_condition:
            ix=rng.randint(0,len(scores),size=(DRAWS,len(scores)))
            estimates.append((scores[ix]*weight[ix]).sum(axis=1)/weight[ix].sum(axis=1))
        draws=(estimates[0]-estimates[1])/gain
        point=(np.average(per_condition[0][0],weights=per_condition[0][1])-
               np.average(per_condition[1][0],weights=per_condition[1][1]))/gain
        recorded=measure.get("correctedDifferenceEstimate")
        if not valid(recorded) or abs(point-recorded)>1e-5:
            failures.append("bootstrap-point-mismatch")
        boot={"se":float(np.std(draws,ddof=1)),"point":float(point),
              "sampleCounts":[len(x[0]) for x in per_condition],"draws":DRAWS}
        if not valid(boot["se"]):failures.append("bootstrap-se-nonfinite")
    else:failures.append("bootstrap-input-missing")
    return {"candidateId":ident,"seed":seed,"valid":not failures,
            "issues":failures,"bootstrap":boot,
            "removedComponents":[x.get("removedIndices",[]) for x in ica] if isinstance(ica,list) else []}


def validate(result):
    if not result.get("completed") or result.get("errors"):raise ValueError("incomplete real EEG record")
    c=result["icaSeedCheck"]
    ids=c["ids"]
    measures=grid(c["candidateMeasurements"],len(ids),"candidateMeasurements")
    original=grid(c["rankStatus"],len(ids),"rankStatus")
    output=[]
    for i,ident in enumerate(ids):
        for seed in range(N_SEEDS):
            # All measured candidates should be recorded, including rejected
            # scientific candidates: Rank feasibility is not recorder success.
            snapshot=measures[i][seed]
            if isinstance(snapshot,dict) and snapshot.get("hasMeasurement"):
                z=inspect_one(snapshot,ident,seed)
                z["rankStatus"]=original[i][seed]
                output.append(z)
    return {"subject":result["subject"],"checked":len(output),
            "available":sum(z["valid"] for z in output),
            "problems":[x for x in output if not x["valid"]],
            "rows":output,"productionRecommendationUnchanged":True}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--input",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    raw=a.input.read_bytes()
    report=validate(json.loads(raw))
    report["rawResultSha256"]=hashlib.sha256(raw).hexdigest()
    report["newAlgorithmSelected"]=False
    report["originalPR51Verdict"]="PARTIAL"
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(report,indent=2))
    print("PR51_SCIENCE_RECORD_AUDIT",json.dumps({
        "subject":report["subject"],"checked":report["checked"],
        "available":report["available"],"failures":len(report["problems"])}))
    if not report["checked"] or report["problems"]:
        raise SystemExit("Recorder integrity failed; original scientific recommendations unchanged")


if __name__=="__main__":
    main()
