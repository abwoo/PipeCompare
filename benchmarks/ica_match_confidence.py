#!/usr/bin/env python3
"""Archived-only ICA topography assignment confidence and ERP co-variation.

This is P1/P2 diagnostics, NOT component contribution estimation, recommender
modification or participant-level confirmation. Protocol frozen in docs.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

from diagnose_ica_component_outliers import centered_maps, labels_for, map_pairs, serial
from validate_ica_science import FROZEN_DEVELOPMENT_IDS, validate_science
from validate_ica_unseen_seeds import grid

SOURCE_RUN=43184910782
WEAK=.90
SEEDS=6
TOL=1e-7


def one_dec(snapshot):
    science=snapshot.get("science")
    if not isinstance(science,dict) or science.get("status")!="complete":
        raise ValueError("Missing real ICA science evidence")
    d=serial(science.get("decompositions"))
    if len(d)!=1:
        raise ValueError("Expected exactly one actual ICA decomposition")
    return d[0]


def aligned(a,b):
    la,A=centered_maps(a)
    lb,B=centered_maps(b)
    if len(set(la))!=len(la) or set(la)!=set(lb):
        raise ValueError("ICA electrode labels do not align")
    B=B[[lb.index(label) for label in la]]
    corr=np.clip(np.abs(A.T@B),0.,1.)
    if not np.isfinite(corr).all():
        raise ValueError("Nonfinite topology correlation")
    return la,A,B,corr


def assignment(a,b):
    labels,A,B,C=aligned(a,b)
    ii,jj=linear_sum_assignment(-C)
    ref_classes=labels_for(a)
    next_classes=labels_for(b)
    ra=set(int(i) for i in serial(a.get("removed")))
    rb=set(int(i) for i in serial(b.get("removed")))
    if (len(ii)!=min(A.shape[1],B.shape[1]) or
        any(i<1 or i>A.shape[1] for i in ra) or
        any(i<1 or i>B.shape[1] for i in rb)):
        raise ValueError("Component count/removal mismatch")
    total=float(sum(C[i,j] for i,j in zip(ii,jj)))
    pz=labels.index("Pz") if "Pz" in labels else None
    diagnostics=[]
    for i,j in zip(ii,jj):
        i=int(i);j=int(j)
        corr=float(C[i,j])
        best_row=float(C[i,:].max())
        best_col=float(C[:,j].max())
        row_nearest=np.flatnonzero(np.isclose(C[i,:],best_row,atol=1e-12,rtol=0.))
        col_nearest=np.flatnonzero(np.isclose(C[:,j],best_col,atol=1e-12,rtol=0.))
        mutual=(len(row_nearest)==1 and len(col_nearest)==1
                and row_nearest[0]==j and col_nearest[0]==i)
        alternative=C.copy()
        alternative[i,j]=-1e9
        i2,j2=linear_sum_assignment(-alternative)
        alt_total=float(sum(alternative[k,l] for k,l in zip(i2,j2)))
        # If no alternate complete matching exists, report unavailable.
        margin=total-alt_total if alt_total>0 else None
        original_removed=i+1 in ra
        alternate_removed=j+1 in rb
        weak=corr<WEAK
        diagnostics.append({
            "component0":i+1,"componentLater":j+1,
            "absoluteCorrelation":corr,
            "nearestCorrelationForSeed0Component":best_row,
            "nearestCorrelationForOtherComponent":best_col,
            "mutualUniqueNearest":mutual,
            "alternativeCompleteAssignmentScoreMargin":margin,
            "isWeak":weak,
            "weakType":("intrinsic-weak" if best_row<WEAK else
                        "assignment-conflict" if weak else "not-weak"),
            "removedSeed0":original_removed,"removedLater":alternate_removed,
            "removalDiscordant":original_removed!=alternate_removed,
            "classSeed0":ref_classes[i],"classLater":next_classes[j],
            "classDiscordant":(ref_classes[i]!=next_classes[j]
                                if ref_classes[i] is not None and next_classes[j] is not None else None),
            "normalizedPzMapLoadingOnly":(float(abs(A[pz,i])) if pz is not None else None),
        })
    weak=[r for r in diagnostics if r["isWeak"]]
    return {
        "mapping":{int(i)+1:int(j)+1 for i,j in zip(ii,jj)},
        "matchedComponents":len(diagnostics),
        "lowCorrelationComponents":len(weak),
        "intrinsicallyWeak":sum(r["weakType"]=="intrinsic-weak" for r in weak),
        "assignmentConflicts":sum(r["weakType"]=="assignment-conflict" for r in weak),
        "weakUnambiguousMutualNearest":sum(r["mutualUniqueNearest"] for r in weak),
        "weakWithRemovalDiscordance":sum(r["removalDiscordant"] for r in weak),
        "weakWithICLabelClassChange":sum(r["classDiscordant"] is True for r in weak),
        "weakPairs":weak,
        "allPairs":diagnostics,
    }


def _waveform(measure):
    z=measure["science"]
    if z.get("channel")!="Pz":raise ValueError("Not Pz ERP")
    conditions=serial(z.get("conditions"))
    if len(conditions)!=2:raise ValueError("Not two-condition ERP")
    times=np.asarray(z.get("timeMs"),dtype=float)
    if times.ndim!=1 or len(times)<30 or not np.isfinite(times).all() or not np.all(np.diff(times)>0):
        raise ValueError("Waveform time axis unavailable")
    waves=[np.asarray(x.get("meanWaveform"),dtype=float) for x in conditions]
    if any(x.shape!=times.shape or not np.isfinite(x).all() for x in waves):
        raise ValueError("Missing complete ERP waveform")
    g=measure.get("gain")
    if not isinstance(g,(int,float)) or not np.isfinite(g) or g<=0:
        raise ValueError("Gain missing")
    diff=(waves[0]-waves[1])/g
    win=(times>=300-1e-7)&(times<=600+1e-7)
    if not win.any() or not math.isclose(float(diff[win].mean()),
           float(measure["correctedDifferenceEstimate"]),rel_tol=1e-7,abs_tol=1e-6):
        raise ValueError("Waveform/source P300 parity failed")
    return times,diff


def _near(a,b,label):
    if not np.isfinite(a) or not np.isfinite(b) or not np.isclose(a,b,atol=TOL,rtol=TOL):
        raise ValueError("Previously independently validated result changed: "+label)


def _cycle(m01,m12,m02):
    consistent=0; assessed=0
    for x,y in m01.items():
        if y in m12 and x in m02:
            assessed+=1
            consistent+=int(m12[y]==m02[x])
    return {"assessed":assessed,"consistent":consistent,
            "fraction":consistent/assessed if assessed else None}


def examine(source,previous):
    sid=source.get("subject")
    if (sid not in FROZEN_DEVELOPMENT_IDS or previous.get("subject")!=sid or
        not source.get("completed") or source.get("errors") or
        previous.get("identityErrors") or previous.get("scienceFailures") or
        not previous.get("unalteredFrozenRankAndEstimates") or
        not previous.get("allPredeclaredFeasibleScienceComplete")):
        raise ValueError("Invalid original/provenance science evidence")
    x=source["icaSeedCheck"];ids=x["ids"]
    if x["recommended"][0]!=FROZEN_DEVELOPMENT_IDS[sid][0]:
        raise ValueError("Original baseline ID changed")
    snapshot=grid(x["candidateMeasurements"],len(ids),"candidateMeasurements")
    prev={q["candidateId"]:q for q in previous["predeclaredCandidateScience"]}
    if set(prev)!=set(FROZEN_DEVELOPMENT_IDS[sid]):
        raise ValueError("Predeclared candidate roster changed")
    outputs=[]
    for candidate in FROZEN_DEVELOPMENT_IDS[sid]:
        i=ids.index(candidate)
        rows=[snapshot[i][seed] for seed in range(SEEDS)]
        checks=[]
        for seed,z in enumerate(rows):
            if not isinstance(z,dict) or z.get("candidateId")!=candidate or z.get("seed")!=seed:
                raise ValueError("Seed/candidate identity mismatch")
            check=validate_science(z,sid,candidate,seed)
            if not check.get("complete"):raise ValueError("Scientific evidence incomplete: "+str(check.get("errors")))
            checks.append(check)
        sr=prev[candidate]
        starts=sr["startReports"]
        if sorted(q["seed"] for q in starts)!=list(range(SEEDS)):
            raise ValueError("Old six-seed science records incomplete")
        starts=sorted(starts,key=lambda z:z["seed"])
        mats=[one_dec(z) for z in rows]
        all_erp=[]
        for seed in range(SEEDS):
            _near(checks[seed]["bootstrapSE"],starts[seed]["bootstrapSE"],"bootstrapSE")
            _near(rows[seed]["correctedDifferenceEstimate"],starts[seed]["P3"],"P300")
            tm,wave=_waveform(rows[seed])
            if seed>0:
                t0,_=_waveform(rows[0])
                if not np.array_equal(tm,t0):raise ValueError("ERP timebase changed")
            all_erp.append(wave)
        previous_matches={z["seed"]:z for z in sr["componentMatchesVsSeed0"]}
        per_seed=[];maps={}
        for seed in range(1,SEEDS):
            diag=assignment(mats[0],mats[seed])
            maps[seed]=diag["mapping"]
            old=previous_matches[seed]
            if not old.get("available"):raise ValueError("Prior component match unavailable")
            # This is a new diagnostic, not a replacement for verified original matching.
            mm=map_pairs({"decompositions":[mats[0]]},{"decompositions":[mats[seed]]})
            _near(mm["minimumCorrelation"],old["minimumMatchedTopographicCorrelation"],"minimumMapCorrelation")
            _near(mm["matchedRemovedJaccard"],old["removedComponentMatchedJaccard"],"removalJaccard")
            delta=all_erp[seed]-all_erp[0]
            rms=float(np.sqrt(np.mean(delta**2)))
            _near(rms,old["fullWaveformRmsDifference"],"fullPzWaveformRMS")
            p3delta=float(rows[seed]["correctedDifferenceEstimate"]-
                          rows[0]["correctedDifferenceEstimate"])
            per_seed.append({
                "seed":seed,
                "minimumAbsoluteScalpMapCorrelation":old["minimumMatchedTopographicCorrelation"],
                "medianAbsoluteScalpMapCorrelation":old["medianAbsoluteTopographicCorrelation"],
                "fullPzWaveformRMSChangeUv":rms,
                "P300DifferenceAmplitudeChangeUv":p3delta,
                "bootstrapSEChangeUv":starts[seed]["bootstrapSE"]-starts[0]["bootstrapSE"],
                "diagnostic":{key:value for key,value in diag.items()
                              if key not in ("mapping","allPairs")},
            })
        composed=assignment(mats[1],mats[2])["mapping"]
        cycle=_cycle(maps[1],composed,maps[2])
        outputs.append({"candidateId":candidate,"anchorSeed":0,"comparisons":per_seed,
                        "cycleConsistency0To1To2":cycle})
    return {"subject":sid,"candidateReadouts":outputs,
            "EEGRefitted":False,"rankUnchanged":True,"holdoutUsed":False,
            "causalComponentAttributionEstablished":False}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--evidence",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    sources=sorted(args.evidence.rglob("result.json"))
    validations=sorted(args.evidence.rglob("ica_science_validated.json"))
    if len(sources)!=3 or len(validations)!=3:
        raise SystemExit(f"Expected three source/validation files; got {len(sources)}/{len(validations)}")
    original={r["subject"]:r for f in sources for r in [json.loads(f.read_text())]}
    verified={r["subject"]:r for f in validations for r in [json.loads(f.read_text())]}
    if set(original)!={11,15,21} or set(verified)!={11,15,21}:
        raise SystemExit("Unexpected, duplicated or missing participants")
    outputs=[examine(original[s],verified[s]) for s in sorted(original)]
    result={"scope":"P1-P2 retrospective matching-confidence association audit",
            "scienceSourceRun":SOURCE_RUN,"allNewEEGFits":False,"newParticipants":False,
            "original20SubjectVerdict":"PARTIAL","pr51MergeAuthorized":False,
            "productionRecommendationChanged":False,"results":outputs}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False))
    for row in outputs:
        for candidate in row["candidateReadouts"]:
            training=candidate["comparisons"][:2]
            print("PR51_ICA_MATCH_CONFIDENCE",json.dumps({
                "subject":row["subject"],"candidateId":candidate["candidateId"],
                "weakPairsStarts1And2":sum(v["diagnostic"]["lowCorrelationComponents"] for v in training),
                "intrinsicWeakStarts1And2":sum(v["diagnostic"]["intrinsicallyWeak"] for v in training),
                "forcedAssignmentWeakStarts1And2":sum(v["diagnostic"]["assignmentConflicts"] for v in training),
                "mutualNearestWeakStarts1And2":sum(v["diagnostic"]["weakUnambiguousMutualNearest"] for v in training),
                "cycleConsistency":candidate["cycleConsistency0To1To2"],
                "waveformRMSStart1Uv":candidate["comparisons"][0]["fullPzWaveformRMSChangeUv"],
                "waveformRMSStart2Uv":candidate["comparisons"][1]["fullPzWaveformRMSChangeUv"]}))
    print("PR51_ICA_MATCH_CONFIDENCE_STOP",json.dumps({
        "cannotAttributeP300ToComponentsWithoutActivations":True,
        "minimumCorrelationGateRelaxed":False,"promote":False,"mergePR51":False}))


if __name__=="__main__":
    main()
