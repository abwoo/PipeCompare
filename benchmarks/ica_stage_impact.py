#!/usr/bin/env python3
"""Retrospective ICA component stage signature vs final ERP association.

Source is the existing full six-start GitHub EEG experiment. This script
NEVER selects a recipe or estimates a causal effect on the final P300.
Read docs/PR51_ICA_STAGE_IMPACT_PROTOCOL_20261113.md first.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np

from ica_match_confidence import assignment
from validate_ica_science import FROZEN_DEVELOPMENT_IDS
from validate_ica_unseen_seeds import grid

SOURCE_RUN=43397978997
WEAK_CORRELATION=0.90
TOL=1e-6


def scalar(v,name):
    if not isinstance(v,(int,float)) or isinstance(v,bool) or not math.isfinite(v):
        raise ValueError("Missing/nonfinite "+name)
    return float(v)


def seq(v):
    return v if isinstance(v,list) else [v] if v is not None else []


def snapshot_decomposition(v):
    science=v.get("science")
    if not isinstance(science,dict) or science.get("status")!="complete":
        raise ValueError("Original science capture missing")
    dec=seq(science.get("decompositions"))
    if len(dec)!=1 or not isinstance(dec[0],dict):
        raise ValueError("Expected one real ICA decomposition")
    return dec[0]


def stage_contributions(v):
    comps=np.asarray(v.get("stageComponentPzDifferenceUv",[]),dtype=float)
    if comps.ndim!=1 or comps.size<1 or not np.isfinite(comps).all():
        raise ValueError("Missing original stage component amplitude")
    sensor=scalar(v.get("stageSensorPzDifferenceUv"),"stage sensor Pz difference")
    residual=scalar(v.get("stageResidualPzDifferenceUv"),"stage residual")
    if not math.isclose(float(comps.sum()+residual),sensor,rel_tol=1e-7,abs_tol=TOL):
        raise ValueError("Pre-removal Pz sensor/component/residual mismatch")
    if not v.get("joinedFinalTrialIdentities"):
        raise ValueError("Original final trial ID join not validated")
    if v.get("equalsFinalP300Estimand") is not False:
        raise ValueError("Stage contribution incorrectly represented as final ERP")
    return comps,sensor,residual


def used_trial_weights(science):
    conditions=seq(science.get("conditions"))
    if len(conditions)!=2:raise ValueError("Missing original two-condition ERP")
    rows=[]
    for c in conditions:
        refs=seq(c.get("referenceIds"))
        kept=seq(c.get("retainedIds"))
        weights=seq(c.get("trialWeights"))
        if not kept or len(set(kept))!=len(kept) or len(refs)!=len(weights) or len(refs)!=len(set(refs)):
            raise ValueError("Malformed original trial/weight provenance")
        by_id=dict(zip(refs,weights))
        row={}
        for ident in kept:
            if ident not in by_id:raise ValueError("Retained event missing from original reference")
            w=scalar(by_id[ident],"retained trial weight")
            if w<=0:raise ValueError("Nonpositive original retained trial weight")
            row[int(ident)]=w
        rows.append(row)
    return rows


def fixed_trial_comparison(a,b):
    aa=used_trial_weights(a)
    bb=used_trial_weights(b)
    event_ids=all(set(x)==set(y) for x,y in zip(aa,bb))
    same_weights=event_ids and all(
        math.isclose(x[k],y[k],rel_tol=1e-10,abs_tol=1e-12)
        for x,y in zip(aa,bb) for k in x)
    return {"sameRetainedTrialIdentity":event_ids,
            "sameOriginalTrialWeights":bool(same_weights),
            "startDifferenceConfoundedByTrialSelectionOrWeights":not same_weights,
            "retainedNSeed0":[len(x) for x in aa],
            "retainedNOther":[len(x) for x in bb]}


def matched_stage_delta(seed0,seedS,a,b):
    xa,sensor0,res0=stage_contributions(seed0)
    xb,sensor1,res1=stage_contributions(seedS)
    decA=snapshot_decomposition(a)
    decB=snapshot_decomposition(b)
    match=assignment(decA,decB)
    if len(xa)!=int(decA["nComponents"]) or len(xb)!=int(decB["nComponents"]):
        raise ValueError("Stage source component count differs from original ICA decomposition")
    rows=[]
    for q in match["allPairs"]:
        i=int(q["component0"])-1
        j=int(q["componentLater"])-1
        delta=float(xb[j]-xa[i])
        rows.append({
            "sourceICSeed0":i+1,"sourceICOther":j+1,
            "absoluteMapCorrelation":q["absoluteCorrelation"],
            "weakMap":bool(q["isWeak"]),
            "retainedInBoth":not q["removedSeed0"] and not q["removedLater"],
            "removedInBoth":q["removedSeed0"] and q["removedLater"],
            "removalDiscordant":q["removalDiscordant"],
            "icLabelClassDiscordant":q["classDiscordant"],
            "stagePzContributionSeed0Uv":float(xa[i]),
            "stagePzContributionOtherUv":float(xb[j]),
            "stageContributionChangeUv":delta,
        })
    weak=[z for z in rows if z["weakMap"]]
    retained_weak=[z for z in weak if z["retainedInBoth"]]
    strong=[z for z in rows if not z["weakMap"]]
    signed_weak=sum(z["stageContributionChangeUv"] for z in weak)
    signed_strong=sum(z["stageContributionChangeUv"] for z in strong)
    unmatched=float(xb.sum()-xa.sum()-(signed_weak+signed_strong))
    residual_delta=float(res1-res0)
    stage_delta=float(sensor1-sensor0)
    if not math.isclose(signed_weak+signed_strong+unmatched+residual_delta,
                        stage_delta,rel_tol=1e-7,abs_tol=TOL):
        raise ValueError("Stage decomposition change does not conserve the observed sensor Pz signal")
    return {
        "matchCount":len(rows),"weakMatchCount":len(weak),
        "weakRetainedBothCount":len(retained_weak),
        "weakWithRemovalDisagreement":sum(z["removalDiscordant"] for z in weak),
        "weakWithClassDisagreement":sum(z["icLabelClassDiscordant"] is True for z in weak),
        "signedWeakMatchedContributionDeltaUv":signed_weak,
        "signedStrongMatchedContributionDeltaUv":signed_strong,
        "unmatchedComponentDeltaUv":unmatched,
        "sensorReconstructionResidualDeltaUv":residual_delta,
        "observedPreRemovalSensorPzDeltaUv":stage_delta,
        "absoluteRetainedWeakDeltaBudgetUv":sum(abs(z["stageContributionChangeUv"]) for z in retained_weak),
        "absoluteOtherMatchedDeltaBudgetUv":sum(abs(z["stageContributionChangeUv"]) for z in rows if not (z["weakMap"] and z["retainedInBoth"])),
        "weakComponentDetails":weak,
        "matchedICCountDiff":int(xb.size-xa.size),
        "stageSensorResidualConserved":True,
        "componentSpecificFinalERPCausalEffectKnown":False,
    }


def evaluate(original,verified,stage):
    sid=int(original.get("subject"))
    if sid not in FROZEN_DEVELOPMENT_IDS or verified.get("subject")!=sid or stage.get("subject")!=sid:
        raise ValueError("Cohort or source mismatches")
    if (not original.get("completed") or original.get("errors") or
        verified.get("identityErrors") or verified.get("scienceFailures") or
        not verified.get("unalteredFrozenRankAndEstimates") or
        not verified.get("allPredeclaredFeasibleScienceComplete")):
        raise ValueError("Original ERP and scientific capture parity failed")
    cm=original["icaSeedCheck"]
    ids=cm["ids"]
    if cm["recommended"][0]!=FROZEN_DEVELOPMENT_IDS[sid][0]:
        raise ValueError("Original Standard choice drifted")
    measurement=grid(cm["candidateMeasurements"],len(ids),"candidateMeasurements")
    status=grid(cm["rankStatus"],len(ids),"rankStatus")
    old={x["candidateId"]:x for x in verified["predeclaredCandidateScience"]}
    new={x["candidateId"]:x for x in stage["candidates"]}
    if set(old)!=set(new) or set(old)!=set(FROZEN_DEVELOPMENT_IDS[sid]):
        raise ValueError("One or more source candidate IDs changed")
    results=[]
    for candidate in FROZEN_DEVELOPMENT_IDS[sid]:
        idx=ids.index(candidate)
        m=[measurement[idx][seed] for seed in range(6)]
        for seed,row in enumerate(m):
            if not isinstance(row,dict) or row.get("candidateId")!=candidate or row.get("seed")!=seed:
                raise ValueError("Candidate / seed alignment broken")
        latest=sorted(new[candidate]["startReports"],key=lambda z:z["seed"])
        if len(latest)!=6 or [x["seed"] for x in latest]!=list(range(6)):
            raise ValueError("New pre-removal stage capture incomplete")
        expected=old[candidate]["componentMatchesVsSeed0"]
        if len(expected)!=5:
            raise ValueError("Previously archived ERP matching observations incomplete")
        by_seed={x["seed"]:x for x in expected}
        comparisons=[]
        for seed in range(1,6):
            cp=matched_stage_delta(latest[0],latest[seed],m[0],m[seed])
            prev=by_seed[seed]
            if prev.get("available") is not True:
                raise ValueError("Prior full ERP waveform comparison absent")
            cp.update(fixed_trial_comparison(m[0]["science"],m[seed]["science"]))
            cp.update({
                "seed":seed,
                "oldRankStatus":status[idx][0],
                "otherRankStatus":status[idx][seed],
                "oldEffectiveUnflagged":m[0]["effectiveUnflagged"],
                "otherEffectiveUnflagged":m[seed]["effectiveUnflagged"],
                "oldPhysicalUnflagged":m[0]["sourceUnflagged"],
                "otherPhysicalUnflagged":m[seed]["sourceUnflagged"],
                "originalFinalGainCorrectedP300DriftUv":scalar(
                    m[seed]["correctedDifferenceEstimate"],"other original P300")-
                    scalar(m[0]["correctedDifferenceEstimate"],"old original P300"),
                "finalFullPzWaveformRMSUv":scalar(prev["fullWaveformRmsDifference"],"verified waveform RMS"),
            })
            comparisons.append(cp)
        results.append({
            "candidateId":candidate,"originalStandard":candidate==FROZEN_DEVELOPMENT_IDS[sid][0],
            "comparisons":comparisons,
            "lowCorrelationComponentsAppearInStage":any(c["weakMatchCount"]>0 for c in comparisons),
            "allTrialsAndWeightsIdenticalAcrossStarts":all(
                c["sameOriginalTrialWeights"] for c in comparisons),
        })
    return {"subject":sid,"candidateReports":results,
            "clinicalAccuracyImprovementDemonstrated":False,
            "newIndependentSubjectValidation":False,
            "recommendationChanged":False}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--evidence",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    arg=parser.parse_args()
    def records(file_name):
        fs=sorted(arg.evidence.rglob(file_name))
        if len(fs)!=3:raise ValueError("Expected exactly three archived "+file_name)
        obj=[json.loads(f.read_text()) for f in fs]
        result={int(x["subject"]):x for x in obj}
        if set(result)!={11,15,21}:raise ValueError("Missing/duplicated/out-of-scope evidence "+file_name)
        return result
    original=records("result.json")
    verified=records("ica_science_validated.json")
    staged=records("ica_contribution_validation.json")
    rows=[evaluate(original[s],verified[s],staged[s]) for s in sorted(original)]
    output={"analysis":"retrospective P1b ICA matched component Pz stage signature; not causal final ERP",
            "upstreamDevelopmentRun":SOURCE_RUN,
            "originalPR51Verdict":"PARTIAL","valid20SubjectMatchedPairs":15,
            "newSubjectsUsed":False,"holdout14Used":False,
            "PR51MergeAuthorized":False,"promotionAuthorized":False,
            "results":rows}
    arg.out.parent.mkdir(parents=True,exist_ok=True)
    arg.out.write_text(json.dumps(output,indent=2,allow_nan=False))
    for subj in rows:
        for cand in subj["candidateReports"]:
            train=cand["comparisons"][:2]
            print("PR51_ICA_STAGE_ASSOCIATION",json.dumps({
                "subject":subj["subject"],"candidateId":cand["candidateId"],
                "weakMatchedPairsSelectionStarts":sum(x["weakMatchCount"] for x in train),
                "retainedWeakPairsSelectionStarts":sum(x["weakRetainedBothCount"] for x in train),
                "weakStageDeltaAbsoluteBudgetUv":[round(x["absoluteRetainedWeakDeltaBudgetUv"],8) for x in train],
                "stageSensorPzDriftUv":[round(x["observedPreRemovalSensorPzDeltaUv"],8) for x in train],
                "finalP300DriftUv":[round(x["originalFinalGainCorrectedP300DriftUv"],8) for x in train],
                "trialSelectionOrWeightConfounding":[x["startDifferenceConfoundedByTrialSelectionOrWeights"] for x in train],
                "allSourceResidualConserved":all(x["stageSensorResidualConserved"] for x in train)}))
    print("PR51_ICA_STAGE_PROMOTION",json.dumps({
        "authorized":False,"sourceParitiesChecked":True,"originalVerdict":"PARTIAL",
        "causalComponentImpactClaimed":False}))


if __name__=="__main__":
    main()
