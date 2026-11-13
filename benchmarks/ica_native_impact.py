#!/usr/bin/env python3
"""Original real IC deletion vs downstream selection/weight P300 stage analysis.

A development-only observational before/after stage contrast, **not**
scientific evidence to promote an ICA recipe or a final-P300 causal claim.
Protocol: docs/PR51_ICA_POSTREMOVAL_MECHANISM_PROTOCOL_20261113.md.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from ica_stage_impact import fixed_trial_comparison
from validate_ica_postremoval import (
    ABS_TOL_UV, finite, final_event_weights, serial,
    stage_weighted_diff, true_event_table,
)
from validate_ica_science import FROZEN_DEVELOPMENT_IDS
from validate_ica_unseen_seeds import grid

SOURCE_RUN=43408389159
TOL=1e-6


def stage(snapshot,name):
    science=snapshot.get("science",{})
    if science.get("status")!="complete":
        raise ValueError("Original science record incomplete")
    key={"pre":"preRemovalContributions","post":"postRemovalSignals"}[name]
    captures=serial(science.get(key))
    if len(captures)!=1 or not isinstance(captures[0],dict):
        raise ValueError("Absent/duplicate original measured "+name+" removal stage")
    table=true_event_table(captures[0],name)
    return captures[0],table


def deleted_source_events(original,pre):
    ds=serial(original.get("science",{}).get("decompositions"))
    if len(ds)!=1:raise ValueError("Ambiguous ICA source decomposition")
    d=ds[0]
    ids=serial(pre["eventIds"])
    comps=np.asarray(pre.get("nativeICComponentUv"),dtype=float)
    count=d.get("nComponents")
    if not isinstance(count,int) or comps.shape!=(count,len(ids)) or not np.isfinite(comps).all():
        raise ValueError("Original pre-removal ICA source matrix does not match events")
    removed=serial(d.get("removed"))
    if len(set(removed))!=len(removed) or any(
        not isinstance(k,int) or k<1 or k>count for k in removed):
        raise ValueError("Actual removed IC indices missing/invalid")
    selected=np.asarray([k-1 for k in removed],dtype=int)
    residual=np.asarray(serial(pre.get("nativeResidualPzUv")),dtype=float)
    if residual.shape!=(len(ids),) or not np.isfinite(residual).all():
        raise ValueError("Actual native EEG.icawinv subspace residual missing")
    # pop_subcomp reconstructs from retained sources when >=1 component
    # is rejected: the lost signal includes the native subspace residual.
    # No rejected components => actual Steps.icRemove skips pop_subcomp.
    removed_signal=(comps[selected].sum(axis=0)+residual) if len(selected) else np.zeros(len(ids))
    condition=serial(pre["conditionIndex"])
    table={int(i):(int(c),float(v)) for i,c,v in zip(ids,condition,removed_signal)}
    if len(table)!=len(ids):raise ValueError("Duplicated ICA original event IDs")
    return table,removed


def difference(new,old):
    return float(new-old)


def assess(old_m,new_m,old_final_stage,new_final_stage,old_erp,new_erp,old_full_wave_rms):
    pre0,p0=stage(old_m,"pre")
    pre1,p1=stage(new_m,"pre")
    post0,q0=stage(old_m,"post")
    post1,q1=stage(new_m,"post")
    original_sets=[set(p0),set(p1),set(q0),set(q1)]
    if len(set(map(frozenset,original_sets)))!=1 or any(
        not (p0[k][0]==p1[k][0]==q0[k][0]==q1[k][0]) for k in p0):
        raise ValueError("Original physical source events are not identical across ICA fits")
    F0=final_event_weights(old_m["science"])
    F1=final_event_weights(new_m["science"])
    # Require every observed physical event needed by either frozen trial set.
    for table in (p0,p1,q0,q1):
        stage_weighted_diff(table,F0)
        stage_weighted_diff(table,F1)
    P0=stage_weighted_diff(p0,F0)
    P1=stage_weighted_diff(p1,F0)
    R0=stage_weighted_diff(q0,F0)
    R1=stage_weighted_diff(q1,F0)
    R1chosen=stage_weighted_diff(q1,F1)
    pre_drift=difference(P1,P0)
    post_frozen_drift=difference(R1,R0)
    deletion_drift=difference(post_frozen_drift,pre_drift)
    weight_selection=difference(R1chosen,R1)
    total_post_change=difference(R1chosen,R0)
    if not math.isclose(total_post_change,post_frozen_drift+weight_selection,
                        rel_tol=1e-7,abs_tol=TOL):
        raise ValueError("Post-removal fixed-event/selection budget not conserved")
    old_removed,old_ids=deleted_source_events(old_m,pre0)
    new_removed,new_ids=deleted_source_events(new_m,pre1)
    removal_predicted=-(stage_weighted_diff(new_removed,F0)-
                        stage_weighted_diff(old_removed,F0))
    if not math.isclose(removal_predicted,deletion_drift,rel_tol=1e-6,abs_tol=ABS_TOL_UV):
        raise ValueError("Actual Pz cross-ICA deletion drift not predicted from observed removed IC sources")
    if (not math.isclose(stage_weighted_diff(p0,F0),
                        old_final_stage["beforeStageOriginalFinalTrialsUv"],
                        abs_tol=TOL,rel_tol=1e-7) or
        not math.isclose(stage_weighted_diff(q0,F0),
                         old_final_stage["afterStageOriginalFinalTrialsUv"],
                         abs_tol=TOL,rel_tol=1e-7) or
        not math.isclose(R1chosen,
                         new_final_stage["afterStageOriginalFinalTrialsUv"],
                         abs_tol=TOL,rel_tol=1e-7)):
        raise ValueError("Original stage observation not reproduced by same events and weights")
    fixed=fixed_trial_comparison(old_m["science"],new_m["science"])
    return {
        "preRemovalFixedEventsSignalDriftUv":pre_drift,
        "postRemovalFixedEventsSignalDriftUv":post_frozen_drift,
        "actualNativeICRemovalEffectDriftUv":deletion_drift,
        "predictedNativeRemovalEffectDriftUv":removal_predicted,
        "postRemovalTrialSelectionOrWeightTermUv":weight_selection,
        "totalPostRemovalStageDriftUv":total_post_change,
        "baselineOriginalPostRemovalStageUv":R0,
        "otherOriginalPostRemovalStageUv":R1chosen,
        "originalICRemovedSeed0":old_ids,"originalICRemovedOther":new_ids,
        "oldRemovedICCount":len(old_ids),"otherRemovedICCount":len(new_ids),
        "originalFinalGainCorrectedP300DriftUv":difference(new_erp,old_erp),
        "originalFullPzERPWaveformRMSUv":old_full_wave_rms,
        "unchangedERPScienceSourceIdentities":True,
        "fixedICStageLinearBridgeVerified":True,
        "fullERPCausalEffectOrSuperiorRecipeEstablished":False,
        **fixed,
    }


def evaluate(original,science,postcheck):
    sid=original.get("subject")
    if (sid not in FROZEN_DEVELOPMENT_IDS or
        science.get("subject")!=sid or postcheck.get("subject")!=sid or
        not original.get("completed") or original.get("errors") or
        science.get("identityErrors") or science.get("scienceFailures") or
        not science.get("unalteredFrozenRankAndEstimates") or
        not science.get("allPredeclaredFeasibleScienceComplete")):
        raise ValueError("Frozen source Rank/P300 and source-trial parity failed")
    c=original["icaSeedCheck"];ids=c["ids"]
    if c["recommended"][0]!=FROZEN_DEVELOPMENT_IDS[sid][0]:
        raise ValueError("Original Standard candidate recommendation changed")
    grid_m=grid(c["candidateMeasurements"],len(ids),"candidateMeasurements")
    old={r["candidateId"]:r for r in science["predeclaredCandidateScience"]}
    post={r["candidateId"]:r for r in postcheck["candidates"]}
    roster=set(FROZEN_DEVELOPMENT_IDS[sid])
    if set(old)!=roster or set(post)!=roster:
        raise ValueError("Original candidate/source roster incomplete")
    subjects=[]
    for candidate in FROZEN_DEVELOPMENT_IDS[sid]:
        idx=ids.index(candidate)
        rows=[grid_m[idx][i] for i in range(6)]
        observed=sorted(post[candidate]["starts"],key=lambda x:x["seed"])
        if len(observed)!=6 or [x["seed"] for x in observed]!=list(range(6)):
            raise ValueError("Required six actual post-removal Pz records unavailable")
        matching={x["seed"]:x for x in old[candidate]["componentMatchesVsSeed0"]}
        if set(matching)!={1,2,3,4,5}:
            raise ValueError("Original full ERP scientific comparisons incomplete")
        comparisons=[]
        for seed in range(1,6):
            if (rows[0].get("candidateId")!=candidate or
                rows[seed].get("candidateId")!=candidate or
                rows[0].get("seed")!=0 or rows[seed].get("seed")!=seed or
                not observed[seed]["nativeRealEEGLABProjectionVerified"]):
                raise ValueError("One candidate/seed identity or actual stage check invalid")
            x=assess(rows[0],rows[seed],observed[0],observed[seed],
                     rows[0]["correctedDifferenceEstimate"],
                     rows[seed]["correctedDifferenceEstimate"],
                     matching[seed]["fullWaveformRmsDifference"])
            x["seed"]=seed
            x["minimumOriginalScalpMapCorrelation"]=matching[seed]["minimumMatchedTopographicCorrelation"]
            comparisons.append(x)
        subjects.append({"candidateId":candidate,"originalStandard":candidate==FROZEN_DEVELOPMENT_IDS[sid][0],
                         "comparisons":comparisons})
    return {"subject":sid,"candidates":subjects,"scientificImprovementValidated":False}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--evidence",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    def collect(name):
        files=sorted(a.evidence.rglob(name))
        if len(files)!=3:raise ValueError("Missing original source evidence: "+name)
        objs=[json.loads(f.read_text()) for f in files]
        by={x["subject"]:x for x in objs}
        if set(by)!={11,15,21}:raise ValueError("Missing or duplicated original participants")
        return by
    r=collect("result.json")
    sci=collect("ica_science_validated.json")
    p=collect("ica_native_validated.json")
    rows=[evaluate(r[s],sci[s],p[s]) for s in sorted(r)]
    output={"study":"P2d actual native EEGLAB IC removal plus residual and original ERP stage; frozen original events",
            "sourceRun":SOURCE_RUN,"originalFrozen20SubjectVerdict":"PARTIAL",
            "validOriginalPairs":15,"newParticipants":False,"frozenHoldoutOpened":False,
            "scientificPromotionAuthorized":False,"PR51MergeAuthorized":False,
            "results":rows}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(output,indent=2,allow_nan=False))
    for subj in rows:
        for cand in subj["candidates"]:
            for row in cand["comparisons"]:
                print("PR51_NATIVE_IC_REMOVAL_MECHANISM",json.dumps({
                    "subject":subj["subject"],"candidate":cand["candidateId"],
                    "seed":row["seed"],
                    "preFixedUv":round(row["preRemovalFixedEventsSignalDriftUv"],8),
                    "postFixedUv":round(row["postRemovalFixedEventsSignalDriftUv"],8),
                    "deletionEffectUv":round(row["actualNativeICRemovalEffectDriftUv"],8),
                    "postTrialWeightUv":round(row["postRemovalTrialSelectionOrWeightTermUv"],8),
                    "postStageDriftUv":round(row["totalPostRemovalStageDriftUv"],8),
                    "finalP300DriftUv":round(row["originalFinalGainCorrectedP300DriftUv"],8),
                    "removedICs":[row["oldRemovedICCount"],row["otherRemovedICCount"]],
                    "weightsChanged":row["startDifferenceConfoundedByTrialSelectionOrWeights"]}))
    print("PR51_NATIVE_IC_REMOVAL_DECISION",json.dumps({
        "scientificProgressIsStageLocalizationOnly":True,"newSelectorPromoted":False,
        "PR51MergeAuthorized":False,"original20SubjectVerdict":"PARTIAL"}))

if __name__=="__main__":
    main()
