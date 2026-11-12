#!/usr/bin/env python3
"""Reconcile ICA robust versus original ERP CORE P3 effects at source trial level.

Development data only. Descriptive paired participant-internal bootstrap,
not an eligibility gate or new confirmatory result.
"""
import argparse
import io
import json
import math
from pathlib import Path
import zipfile

import numpy as np
from scipy.io import loadmat

BOOTSTRAP_DRAWS = 10000
BOOTSTRAP_SEED = 511515


def as_vector(x):
    return np.asarray(x, dtype=float).reshape(-1)


def objective_trials(candidate):
    o = candidate.m.objectives
    X = o.X
    W = o.w
    result = []
    for k in range(2):
        x = as_vector(X[k])
        if np.size(W) == 0:
            w = np.ones(len(x), dtype=float)
        else:
            w = as_vector(W[k])
        if len(w) != len(x):
            raise ValueError("candidate trial scores/weight vectors are misaligned")
        if np.any(np.isfinite(w) & (w < 0)):
            raise ValueError("negative ERP trial weights")
        ok = np.isfinite(x) & np.isfinite(w) & (w > 0)
        maskw = np.where(ok,w,0.)
        maskx = np.where(ok,x,0.)
        if not maskw.sum() > 0:
            raise ValueError("candidate has no finite positive weighted ERP trials")
        estimate = float(np.dot(maskx,maskw) / maskw.sum())
        claimed = as_vector(o.estimate)[k]
        if not np.isclose(estimate, claimed, rtol=0, atol=1e-8):
            raise ValueError(f"trial-level estimator does not reproduce MATLAB candidate: {estimate} vs {claimed}")
        result.append((maskx,maskw,ok))
    return result


def weighted_mean(x,w):
    if x.shape != w.shape or np.any(w < 0):
        raise ValueError("invalid weights")
    denom=float(w.sum())
    if denom <= 0:raise ValueError("empty weighted mean")
    return float(np.sum(x*w)/denom)


def bootstrap_decomposition(old, new, old_gain, new_gain, draws=BOOTSTRAP_DRAWS, seed=BOOTSTRAP_SEED):
    if not (math.isfinite(old_gain) and math.isfinite(new_gain) and old_gain>0 and new_gain>0):
        raise ValueError("invalid candidate gain")
    if len(old) != 2 or len(new) != 2:
        raise ValueError("two ERP conditions required")
    rng=np.random.default_rng(seed)
    diff_old=np.zeros(draws,dtype=float)
    diff_new=np.zeros(draws,dtype=float)
    enough=np.ones(draws,dtype=bool)
    old_means=[];new_means=[];cross_means=[];per_condition=[]
    for k,(oc,nc) in enumerate(zip(old,new)):
        old_x,old_w,_=oc;new_x,new_w,_=nc
        n=len(old_x)
        if n!=len(new_x):raise ValueError("different source trial axes")
        mo=weighted_mean(old_x,old_w)
        mn=weighted_mean(new_x,new_w)
        cross=weighted_mean(new_x,old_w)
        old_means.append(mo);new_means.append(mn);cross_means.append(cross)
        valid=(old_w>0) & (new_w>0)
        shared_delta=float(np.max(np.abs(old_x[valid]-new_x[valid]))) if any(valid) else None
        per_condition.append({
            "conditionIndex":k,
            "referenceTrials":n,"oldPositiveWeights":int(np.sum(old_w>0)),
            "newPositiveWeights":int(np.sum(new_w>0)),
            "allScoresIdentical": bool(np.allclose(old_x,new_x,rtol=0,atol=1e-10)),
            "allWeightsIdentical":bool(np.allclose(old_w,new_w,rtol=0,atol=1e-10)),
            "maxSharedScoreAbsDifference":shared_delta,
            "oldMean":mo,"newMean":mn,
            "newScoresWithOldWeightsMean":cross,
            "weightRangeOld":[float(old_w.min()),float(old_w.max())],
            "weightRangeNew":[float(new_w.min()),float(new_w.max())]
        })
        idx=rng.integers(0,n,size=(draws,n))
        old_num=np.sum(old_x[idx]*old_w[idx],axis=1)
        old_den=np.sum(old_w[idx],axis=1)
        new_num=np.sum(new_x[idx]*new_w[idx],axis=1)
        new_den=np.sum(new_w[idx],axis=1)
        good=(old_den>0)&(new_den>0)
        enough &= good
        old_est=np.divide(old_num,old_den,out=np.zeros(draws),where=old_den>0)
        new_est=np.divide(new_num,new_den,out=np.zeros(draws),where=new_den>0)
        sign=1 if k==0 else -1
        diff_old+=sign*old_est/old_gain
        diff_new+=sign*new_est/new_gain
    old_raw=old_means[0]-old_means[1]
    new_raw=new_means[0]-new_means[1]
    cross_raw=cross_means[0]-cross_means[1]
    signal_change=cross_raw-old_raw
    weight_change=new_raw-cross_raw
    total=new_raw/new_gain-old_raw/old_gain
    gain_change=new_raw/new_gain-new_raw/old_gain
    if not np.isclose(total,(signal_change+weight_change)/old_gain+gain_change,atol=1e-9,rtol=0):
        raise ValueError("signal/weight/gain decomposition is not additive")
    valid_delta=(diff_new-diff_old)[enough]
    if len(valid_delta)<draws*.9:
        raise ValueError("insufficient valid paired bootstrap resamples")
    ci=np.quantile(valid_delta,[.025,.975])
    return {
        "perCondition":per_condition,
        "allSourceScoresIdentical":all(x["allScoresIdentical"] for x in per_condition),
        "oldRawDifference":old_raw, "newRawDifference":new_raw,
        "oldCorrectedDifference":old_raw/old_gain,
        "newCorrectedDifference":new_raw/new_gain,
        "gainOriginal":old_gain,"gainRobust":new_gain,
        "signalTermRaw":signal_change,
        "weightTermRaw":weight_change,
        "gainTermCorrected":gain_change,
        "netCorrectedDifference":total,
        "pairedBootstrap": {
            "drawsPlanned":draws, "drawsValid":int(len(valid_delta)),
            "seed":seed, "meanDelta":float(valid_delta.mean()),
            "ci95":[float(ci[0]),float(ci[1])],
            "estimator":"participant-internal resampling of same source events within each P3 condition; original vs robust pipeline calibrated effect difference; exploratory only"
        }
    }


def analyze_matlab(mat_content, record):
    if not record["completed"] or record["errors"] or record.get("notRun")!=0:
        raise ValueError("failed original ERP CORE participant")
    s=record["icaSeedCheck"]
    chosen=s["robust"]["candidateId"]
    old=s["recommended"][0]
    if s["robust"]["status"]!="selected" or not s["robust"]["allSeedsScientificallyFeasible"]:
        return {"subject":record["subject"],"status":"not-eligible","robustId":None}
    a=loadmat(io.BytesIO(mat_content),squeeze_me=True,struct_as_record=False)["audit"]
    cands={int(c.id):c for c in np.atleast_1d(a.cands)}
    if old not in cands or chosen not in cands:
        raise ValueError(f"selected candidate absent from MATLAB score audit: {old} {chosen}")
    old_gain=float(record["erp"]["gain"])
    new_gain=float(s["robust"]["seedMeasurements"][0]["gain"])
    result=bootstrap_decomposition(objective_trials(cands[old]),
                                    objective_trials(cands[chosen]),
                                    old_gain,new_gain)
    for field,claimed in (("oldCorrectedDifference",record["erp"]["differenceAmplitude"]),
                          ("newCorrectedDifference",s["robust"]["seedMeasurements"][0]["correctedDifferenceEstimate"])):
        if not math.isclose(result[field],claimed,abs_tol=1e-8,rel_tol=0):
            raise ValueError(f"reconstruction differs from real pilot {field}")
    return {"subject":record["subject"],"status":"exploratory-effect-decomposition",
            "originalCandidateId":old,"robustCandidateId":chosen,
            "results":result,"scientificGatesChanged":False,
            "establishesSuperiority":False,"requiresNewFrozenHoldout":True}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--evidence-root",type=Path,required=True)
    p.add_argument("--subjects",nargs="+",type=int,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    rows=[]
    for subject in args.subjects:
        folder=args.evidence_root/str(subject)
        path=folder/"scoring_audit.mat"
        report=folder/"result.json"
        if not path.is_file() or not report.is_file():
            raise FileNotFoundError(f"missing complete original GH action artifact {subject}")
        rows.append(analyze_matlab(path.read_bytes(),json.loads(report.read_text())))
    if len({x["subject"] for x in rows})!=len(args.subjects):
        raise ValueError("duplicated subject evidence")
    result={"scope":"three original REAL-P3 development participants, NO heldout EEG",
            "notConfirmatory":True,"sourceRun":43162539734,
            "rows":rows}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2))
    for x in rows:
        if x["status"]!="exploratory-effect-decomposition":
            print("PR51_ICA_TRIAL_WEIGHTS",json.dumps(x))
            continue
        v=x["results"]
        print("PR51_ICA_TRIAL_WEIGHTS",json.dumps({
            "subject":x["subject"],"original":x["originalCandidateId"],"robust":x["robustCandidateId"],
            "same_scores":v["allSourceScoresIdentical"],
            "weight_term_raw":v["weightTermRaw"],
            "signal_term_raw":v["signalTermRaw"],
            "calibrated_effect_delta":v["netCorrectedDifference"],
            "pairedCI95":v["pairedBootstrap"]["ci95"]}))


if __name__=="__main__":
    main()
