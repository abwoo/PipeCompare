#!/usr/bin/env python3
"""P2c real EEGLAB native kept-component post-ICA Pz reconstruction.

Original P2 pinv(W)-based before-minus-removed prediction FAILED real EEG
and is retained in the GitHub evidence. This separate, preregistered
hypothesis tests actual pop_subcomp kept component reconstruction.
"""
import argparse
import json
import math
from pathlib import Path
import numpy as np
from validate_ica_postremoval import (
    ABS_TOL_UV,REL_TOL,finite,final_event_weights,serial,
    stage_weighted_diff,true_event_table)
from validate_ica_science import FROZEN_DEVELOPMENT_IDS
from validate_ica_unseen_seeds import grid

SEEDS=6

def native_reconstruction(pre,post,dec,science):
    before=true_event_table(pre,"before real IC removal")
    after=true_event_table(post,"after real IC removal")
    ids=serial(pre["eventIds"])
    if set(before)!=set(after) or any(before[i][0]!=after[i][0] for i in ids):
        raise ValueError("Original physical source time-locked event identities not preserved")
    if not math.isclose(pre["srate"],post["srate"],rel_tol=0,abs_tol=1e-9):
        raise ValueError("ICA removal altered sampling rate unexpectedly")
    n=pre.get("nComponents")
    if not isinstance(n,int) or n<1 or dec.get("nComponents")!=n:
        raise ValueError("Actual ICA component count inconsistent")
    try:
        source=np.asarray(pre["nativeICComponentUv"],dtype=float)
        residual=np.asarray(serial(pre["nativeResidualPzUv"]),dtype=float)
        old=np.asarray(pre["componentUv"],dtype=float)
    except (ValueError,TypeError,KeyError) as e:
        raise ValueError("Original actual EEG.icawinv source unavailable") from e
    if (source.shape!=(n,len(ids)) or old.shape!=(n,len(ids)) or
        residual.shape!=(len(ids),) or
        not np.isfinite(source).all() or not np.isfinite(residual).all() or
        not np.isfinite(old).all()):
        raise ValueError("Native component or residual dimensions/values invalid")
    pre_y=np.asarray([before[i][1] for i in ids])
    post_y=np.asarray([after[i][1] for i in ids])
    if not np.allclose(source.sum(axis=0)+residual,pre_y,atol=1e-6,rtol=1e-7):
        raise ValueError("Original before-removal native source residual does not reproduce actual Pz")
    if not math.isclose(float(np.max(np.abs(residual))),
                        float(pre["nativeMaxAbsResidualUv"]),abs_tol=1e-7):
        raise ValueError("Native subspace residual identity changed")
    r=serial(dec.get("removed"))
    if len(set(r))!=len(r) or any(not isinstance(i,int) or i<1 or i>n for i in r):
        raise ValueError("Removed IC indices not original/original class decision")
    keep=np.asarray([i for i in range(n) if i+1 not in r],dtype=int)
    dropped=np.asarray([i-1 for i in r],dtype=int)
    retained=source[keep,:].sum(axis=0) if len(keep) else np.zeros(len(ids))
    erased=source[dropped,:].sum(axis=0) if len(dropped) else np.zeros(len(ids))
    if len(r)>0:
        predicted_post=retained
        predicted_before_minus_after=erased+residual
    else:
        predicted_post=pre_y
        predicted_before_minus_after=np.zeros(len(ids))
    difference=post_y-predicted_post
    difference2=(pre_y-post_y)-predicted_before_minus_after
    max_abs=float(np.max(np.abs(difference)))
    if (not np.allclose(post_y,predicted_post,atol=ABS_TOL_UV,rtol=REL_TOL) or
        not np.allclose(pre_y-post_y,predicted_before_minus_after,
                        atol=ABS_TOL_UV,rtol=REL_TOL)):
        raise ValueError(f"Native kept-component Pz prediction fails real pop_subcomp: max={max_abs:.9g} uV")
    original_wrong=(pre_y-post_y)-(
        old[dropped,:].sum(axis=0) if len(dropped) else np.zeros(len(ids)))
    events=final_event_weights(science)
    before_ave=stage_weighted_diff(before,events)
    after_ave=stage_weighted_diff(after,events)
    return {
        "nComponents":n,"nativeICMixingSourceVerified":True,
        "nativeMixingVsPinvMaxAbs":float(pre["nativeVsPinvMixingMaxAbs"]),
        "nativePzSubspaceResidualMaxAbsUv":float(np.max(np.abs(residual))),
        "nativePzSubspaceResidualRMSUv":float(np.sqrt(np.mean(residual**2))),
        "originalPseudoinverseHypothesisMaxErrorUv":float(np.max(np.abs(original_wrong))),
        "actualNativeKeptProjectionMaxErrorUv":max_abs,
        "actualNativeDeletionPlusResidualMaxErrorUv":float(np.max(np.abs(difference2))),
        "actualICsRemoved":r,"numberRemovedICs":len(r),
        "physicalUreventCount":len(ids),
        "beforeStageOriginalFinalTrialsUv":before_ave,
        "afterStageOriginalFinalTrialsUv":after_ave,
        "directRemovalStageDeltaUv":after_ave-before_ave,
        "nativeRealEEGLABProjectionVerified":True,
        "fullFinalP300CausalityOrSuperiorAlgorithmShown":False,
    }

def examine(real,old,prior):
    sid=real.get("subject")
    if (sid not in FROZEN_DEVELOPMENT_IDS or old.get("subject")!=sid or
        prior.get("subject")!=sid or not real.get("completed") or real.get("errors") or
        old.get("identityErrors") or old.get("scienceFailures") or
        not old.get("unalteredFrozenRankAndEstimates") or
        not old.get("allPredeclaredFeasibleScienceComplete")):
        raise ValueError("Real original ICA, Rank and original P300 parity failed")
    c=real["icaSeedCheck"]
    if c["recommended"][0]!=FROZEN_DEVELOPMENT_IDS[sid][0]:
        raise ValueError("Original Standard candidate changed")
    candidates=grid(c["candidateMeasurements"],len(c["ids"]),"candidateMeasurements")
    verified={r["candidateId"]:r for r in prior["candidates"]}
    if set(verified)!=set(FROZEN_DEVELOPMENT_IDS[sid]):
        raise ValueError("Prior tested candidate roster incomplete")
    out=[]
    for candidate in FROZEN_DEVELOPMENT_IDS[sid]:
        idx=c["ids"].index(candidate)
        stages=sorted(verified[candidate]["startReports"],key=lambda x:x["seed"])
        if len(stages)!=SEEDS or [r["seed"] for r in stages]!=list(range(SEEDS)):
            raise ValueError("Old pre-deletion validation incomplete")
        starts=[]
        for seed in range(SEEDS):
            m=candidates[idx][seed]
            if m.get("candidateId")!=candidate or m.get("seed")!=seed:
                raise ValueError("Source candidate and ICA seed identity changed")
            science=m.get("science",{})
            if science.get("status")!="complete":
                raise ValueError("Original ICA science evidence not complete")
            pre=serial(science.get("preRemovalContributions"))
            post=serial(science.get("postRemovalSignals"))
            dec=serial(science.get("decompositions"))
            if len(pre)!=1 or len(post)!=1 or len(dec)!=1:
                raise ValueError(f"Missing true native pre/post ICA stage, {sid}/{candidate}/{seed}")
            stage=native_reconstruction(pre[0],post[0],dec[0],science)
            if not math.isclose(stage["beforeStageOriginalFinalTrialsUv"],
                                stages[seed]["stageSensorPzDifferenceUv"],
                                abs_tol=1e-6,rel_tol=1e-7):
                raise ValueError("Previous original pre-deletion stage observations changed")
            starts.append({"seed":seed,**stage,
                           "originalFinalGainCorrectedP300Uv":m["correctedDifferenceEstimate"]})
        out.append({"candidateId":candidate,"starts":starts})
    return {"subject":sid,"candidates":out,
            "newParticipants":False,"sourceAlgorithmOrRankChanged":False}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--new",type=Path,required=True)
    p.add_argument("--old",type=Path,required=True)
    p.add_argument("--pre",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    row=examine(json.loads(a.new.read_text()),json.loads(a.old.read_text()),
                json.loads(a.pre.read_text()))
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(row,indent=2,allow_nan=False))
    for x in row["candidates"]:
        print("PR51_NATIVE_POP_SUBCOMP_SOURCE",json.dumps({
            "subject":row["subject"],"candidateId":x["candidateId"],
            "starts":len(x["starts"]),
            "maxActualMixingVsPinv":max(z["nativeMixingVsPinvMaxAbs"] for z in x["starts"]),
            "maxNativeSubspaceResidualUv":max(z["nativePzSubspaceResidualMaxAbsUv"] for z in x["starts"]),
            "maxOldWrongHypothesisErrorUv":max(z["originalPseudoinverseHypothesisMaxErrorUv"] for z in x["starts"]),
            "maxActualNativeKeptProjectionErrorUv":max(z["actualNativeKeptProjectionMaxErrorUv"] for z in x["starts"]),
            "nativeLinearSensorReconstructionVerified":True,
            "promotionAuthorized":False}))
    print("PR51_NATIVE_RECONSTRUCTION_DECISION",json.dumps({
        "originalP2FailedHypothesisRetained":True,"newEEGLABEquationValid":True,
        "sourceAndRankParityPassed":True,"newMethodOutperformsStandard":False,
        "PR51MergeAuthorized":False,"untouchedHoldoutUsed":False}))

if __name__=="__main__":
    main()
