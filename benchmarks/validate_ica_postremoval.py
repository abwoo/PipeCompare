#!/usr/bin/env python3
"""P2 science-only eventwise before/after IC removal bridge.

Runs only on formerly inspected subjects 11/15/21. Requires the unchanged
original Rank and science validation before this script in Actions.
Never modifies score, selector, source EEG, or post hoc thresholds.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np

from validate_ica_science import FROZEN_DEVELOPMENT_IDS
from validate_ica_unseen_seeds import grid

ABS_TOL_UV=1e-4
REL_TOL=1e-6
SEEDS=6

def serial(a):
    return a if isinstance(a,list) else [a] if a is not None else []

def finite(x):
    return isinstance(x,(int,float)) and not isinstance(x,bool) and math.isfinite(x)

def true_event_table(stage,stage_name):
    if not isinstance(stage,dict) or stage.get("status")!="complete":
        raise ValueError(stage_name+" missing/failed science observer: "+
                         str(stage.get("reason","") if isinstance(stage,dict) else stage))
    ids=serial(stage.get("eventIds"))
    cond=serial(stage.get("conditionIndex"))
    volts=serial(stage.get("observedPzUv"))
    n=stage.get("nEvents")
    if (not isinstance(n,int) or n<2 or
        len(ids)!=n or len(cond)!=n or len(volts)!=n or
        any(not isinstance(i,int) for i in ids) or len(set(ids))!=n or
        set(cond)!={1,2} or any(not finite(v) for v in volts) or
        not finite(stage.get("srate")) or stage["srate"]<=0):
        raise ValueError(stage_name+" event/sensor provenance is malformed")
    return {i:(c,float(y)) for i,c,y in zip(ids,cond,volts)}

def final_event_weights(science):
    z=serial(science.get("conditions"))
    if len(z)!=2:
        raise ValueError("Expected exactly two original validated ERP conditions")
    selections=[]
    for ci,condition in enumerate(z,1):
        ids=serial(condition.get("referenceIds"))
        kept=serial(condition.get("retainedIds"))
        weights=serial(condition.get("trialWeights"))
        if (not kept or len(ids)!=len(weights) or
            len(set(ids))!=len(ids) or len(set(kept))!=len(kept)):
            raise ValueError("Original final trial/weight alignment malformed")
        by_id=dict(zip(ids,weights))
        selected={}
        for i in kept:
            if i not in by_id or not finite(by_id[i]) or by_id[i]<=0:
                raise ValueError("Invalid original retained positive-weight trial")
            selected[int(i)]=float(by_id[i])
        selections.append(selected)
    return selections

def stage_weighted_diff(events,selected):
    ave=[]
    for ci,sel in enumerate(selected,1):
        if any(i not in events or events[i][0]!=ci for i in sel):
            raise ValueError("Original retained final trial missing from stage")
        vals=np.asarray([events[i][1] for i in sel],dtype=float)
        w=np.asarray(list(sel.values()),dtype=float)
        ave.append(float(np.dot(vals,w)/w.sum()))
    return ave[0]-ave[1]

def verify_bridge(pre,post,decomposition,science):
    p=true_event_table(pre,"pre")
    q=true_event_table(post,"post")
    if set(p)!=set(q) or any(p[i][0]!=q[i][0] for i in p):
        raise ValueError("Pre/post IC-removal physical event identities changed")
    if not math.isclose(pre["srate"],post["srate"],abs_tol=1e-9):
        raise ValueError("ICA removal unexpectedly changed sampling rate")
    nc=decomposition.get("nComponents")
    if not isinstance(nc,int) or nc<=0 or pre.get("nComponents")!=nc:
        raise ValueError("ICA component count mismatch with original decomposition")
    try:
        amplitude=np.asarray(pre.get("componentUv"),dtype=float)
    except (ValueError,TypeError) as exc:
        raise ValueError("Malformed pre-removal component signal") from exc
    ids=serial(pre["eventIds"])
    if amplitude.shape!=(nc,len(ids)) or not np.isfinite(amplitude).all():
        raise ValueError("Original ICA components do not match event-level measurements")
    raw_removed=serial(decomposition.get("removed"))
    if (len(set(raw_removed))!=len(raw_removed) or
        any(not isinstance(i,int) or i<1 or i>nc for i in raw_removed)):
        raise ValueError("Removed IC identities invalid")
    removed=np.asarray([int(x)-1 for x in raw_removed],dtype=int)
    actual=np.asarray([p[i][1]-q[i][1] for i in ids],dtype=float)
    predicted=amplitude[removed,:].sum(axis=0) if removed.size else np.zeros(len(ids))
    if not np.allclose(actual,predicted,atol=ABS_TOL_UV,rtol=REL_TOL):
        max_error=float(np.max(np.abs(actual-predicted)))
        raise ValueError(f"Actual post-ICA Pz deletion does not match removed IC activation: max={max_error:.9g} uV")
    joined=final_event_weights(science)
    pre_mean=stage_weighted_diff(p,joined)
    post_mean=stage_weighted_diff(q,joined)
    if not np.isfinite([pre_mean,post_mean]).all():
        raise ValueError("Invalid original retained-trial stage difference")
    return {
        "physicalUreventCount":len(ids),
        "sourceEventIDsPreserved":True,
        "sourceConditionIdentitiesPreserved":True,
        "removedComponentCount":int(removed.size),
        "removedComponentIDs":[int(v) for v in raw_removed],
        "stepPzSignalChangeMaxAbsUv":float(np.max(np.abs(actual))),
        "predictedRemovedICChangeMaxAbsUv":float(np.max(np.abs(predicted))),
        "linearPzRemovalMaxResidualUv":float(np.max(np.abs(actual-predicted))),
        "linearPzRemovalRMSEUv":float(np.sqrt(np.mean((actual-predicted)**2))),
        "preRemovalFinalSelectedTrialPzDifferenceUv":pre_mean,
        "postRemovalFinalSelectedTrialPzDifferenceUv":post_mean,
        "actualICRemovalFinalSelectedTrialPzDeltaUv":post_mean-pre_mean,
        "finalGainCorrectedP300MeasurementIdentical":False,
        "directICStepLinearBridgeValid":True,
        "downstreamCausalityOrAccuracyEstablished":False,
    }

def examine(new,scientific,prior_stage):
    sid=new.get("subject")
    if (sid not in FROZEN_DEVELOPMENT_IDS or scientific.get("subject")!=sid or
        prior_stage.get("subject")!=sid or not new.get("completed") or new.get("errors") or
        scientific.get("identityErrors") or scientific.get("scienceFailures") or
        not scientific.get("unalteredFrozenRankAndEstimates") or
        not scientific.get("allPredeclaredFeasibleScienceComplete")):
        raise ValueError("Original unmodified Rank/P300 and event evidence not verified")
    c=new["icaSeedCheck"]
    if c["recommended"][0]!=FROZEN_DEVELOPMENT_IDS[sid][0]:
        raise ValueError("Frozen original Standard recommendation drifted")
    ids=c["ids"]
    measurements=grid(c["candidateMeasurements"],len(ids),"candidateMeasurements")
    validated={x["candidateId"]:x for x in prior_stage["candidates"]}
    if set(validated)!=set(FROZEN_DEVELOPMENT_IDS[sid]):
        raise ValueError("Prior pre-removal candidate roster incomplete")
    results=[]
    for candidate in FROZEN_DEVELOPMENT_IDS[sid]:
        i=ids.index(candidate)
        archived=sorted(validated[candidate]["startReports"],key=lambda x:x["seed"])
        if len(archived)!=SEEDS or [v["seed"] for v in archived]!=list(range(SEEDS)):
            raise ValueError("Prior pre-removal stage not complete for all starts")
        starts=[]
        for seed in range(SEEDS):
            m=measurements[i][seed]
            if not isinstance(m,dict) or m.get("candidateId")!=candidate or m.get("seed")!=seed:
                raise ValueError("Actual candidate or seed archive identity invalid")
            science=m.get("science",{})
            pre=serial(science.get("preRemovalContributions"))
            post=serial(science.get("postRemovalSignals"))
            dec=serial(science.get("decompositions"))
            if len(pre)!=1 or len(post)!=1 or len(dec)!=1:
                raise ValueError(f"Missing pre/post real observer or one ICA fit {sid}/{candidate}/{seed}: {len(pre)}/{len(post)}/{len(dec)}")
            stage=verify_bridge(pre[0],post[0],dec[0],science)
            old=archived[seed]
            if (not math.isclose(stage["preRemovalFinalSelectedTrialPzDifferenceUv"],
                                old["stageSensorPzDifferenceUv"],abs_tol=1e-6,rel_tol=1e-7)):
                raise ValueError("Previously validated pre-removal Pz result no longer reproduces")
            starts.append({"seed":seed,**stage,
                           "actualFinalGainCorrectedP300Uv":m["correctedDifferenceEstimate"],
                           "originalEffectiveUnflagged":m["effectiveUnflagged"],
                           "originalPhysicalUnflagged":m["sourceUnflagged"]})
        results.append({"candidateId":candidate,"starts":starts})
    return {"subject":sid,"candidateReadouts":results,
            "standardUnchanged":True,"scientificGateAdjusted":False,
            "newSubjects":False,"candidateSelectionChanged":False}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--evidence",type=Path)
    p.add_argument("--new",type=Path)
    p.add_argument("--old",type=Path)
    p.add_argument("--pre",type=Path)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    if a.new is not None:
        if a.evidence is not None or a.old is None or a.pre is None:
            raise ValueError("Single subject requires --new, --old, --pre and no --evidence")
        real=json.loads(a.new.read_text())
        science=json.loads(a.old.read_text())
        former=json.loads(a.pre.read_text())
        report=[examine(real,science,former)]
    else:
        if a.evidence is None or a.old is not None or a.pre is not None:
            raise ValueError("Aggregate run requires exactly --evidence")
        def by_subject(name):
            files=sorted(a.evidence.rglob(name))
            if len(files)!=3:
                raise ValueError("Expected three source-artifacts for "+name)
            objs=[json.loads(path.read_text()) for path in files]
            tab={x["subject"]:x for x in objs}
            if set(tab)!={11,15,21}:raise ValueError("Missing/duplicated/out-of-scope participant")
            return tab
        result=by_subject("result.json")
        validated=by_subject("ica_science_validated.json")
        stage=by_subject("ica_contribution_validation.json")
        report=[examine(result[s],validated[s],stage[s]) for s in sorted(result)]
    output={"study":"read-only original actual post-ICA deletion Pz linear event bridge",
            "original20SubjectVerdict":"PARTIAL","validOriginalPairs":15,
            "realDevelopmentOnly":True,"frozen14HoldoutOpened":False,
            "pr51MergeAuthorized":False,"algorithmPromotionAuthorized":False,
            "results":report}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(output,indent=2,allow_nan=False))
    for row in report:
        for candidate in row["candidateReadouts"]:
            print("PR51_ICA_REMOVAL_BRIDGE",json.dumps({
                "subject":row["subject"],"candidateId":candidate["candidateId"],
                "starts":len(candidate["starts"]),
                "maxRemovalResidualUv":max(x["linearPzRemovalMaxResidualUv"] for x in candidate["starts"]),
                "removedCounts":[x["removedComponentCount"] for x in candidate["starts"]],
                "postRemovalStageRangeUv":max(x["postRemovalFinalSelectedTrialPzDifferenceUv"] for x in candidate["starts"])-
                                          min(x["postRemovalFinalSelectedTrialPzDifferenceUv"] for x in candidate["starts"]),
                "finalCalibratedP3RangeUv":max(x["actualFinalGainCorrectedP300Uv"] for x in candidate["starts"])-
                                          min(x["actualFinalGainCorrectedP300Uv"] for x in candidate["starts"]),
                "performanceImprovementProven":False}))
    print("PR51_ICA_REMOVAL_BRIDGE_VERDICT",json.dumps({
        "realStageRemovalsValidated":True,"preRemovalRankParityEnforced":True,
        "pr51MergeAuthorized":False,"holdoutUsed":False,"newICAMethodValid":False}))

if __name__=="__main__":
    main()
