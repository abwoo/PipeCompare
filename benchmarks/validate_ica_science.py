#!/usr/bin/env python3
"""Audit newly recorded ICA and ERP science against the immutable six-start run.

This is a READOUT, not a new recommender. The prospective subject-level
holdout [24..39 excluding 30,35] is never opened. A failed recorder never
changes the candidate ranking; CI fails separately while preserving evidence.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

from validate_ica_unseen_seeds import grid, valid

SEEDS=6
BOOT_DRAWS=1000
BOOT_SEED=520004
# Declared from earlier development comparisons, before new science fields.
FROZEN_DEVELOPMENT_IDS={11:[12,25],15:[5,2,6],21:[25,4]}
REL_TOL=1e-7
ABS_TOL=1e-6


def array(x):
    return np.asarray(x,dtype=float)


def approx(a,b):
    if a is None or b is None:
        return a is None and b is None
    if isinstance(a,(list,tuple)) or isinstance(b,(list,tuple)):
        aa=array(a);bb=array(b)
        return aa.shape==bb.shape and np.allclose(aa,bb,atol=ABS_TOL,rtol=REL_TOL,equal_nan=True)
    return bool(np.isclose(float(a),float(b),atol=ABS_TOL,rtol=REL_TOL,equal_nan=True))


def validate_identity(new,prior):
    problems=[]
    for key in ("subject","algorithmCommit","protocolSha256","measureSha256"):
        if new.get(key)!=prior.get(key):
            problems.append(f"{key}-changed")
    for key in ("completed","notRun"):
        if new.get(key)!=prior.get(key):
            problems.append(f"{key}-changed")
    if new.get("errors"):
        problems.append("matlab-execution-errors")
    a,b=new["icaSeedCheck"],prior["icaSeedCheck"]
    for key in ("ids","recommended"):
        if a.get(key)!=b.get(key):
            problems.append(f"{key}-changed")
    for key in ("objective","rankStatus","strata"):
        if key=="objective":
            if not approx(a[key],b[key]):problems.append("objective-changed")
        else:
            n=len(a["ids"])
            if grid(a[key],n,key)!=grid(b[key],n,key):
                problems.append(f"{key}-changed")
    if not approx(new["erp"]["differenceAmplitude"],prior["erp"]["differenceAmplitude"]):
        problems.append("original-p3-estimate-changed")
    n=len(a["ids"])
    aa=grid(a["candidateMeasurements"],n,"new candidateMeasurements")
    bb=grid(b["candidateMeasurements"],n,"prior candidateMeasurements")
    for i,ident in enumerate(a["ids"]):
        for s in range(SEEDS):
            x,y=aa[i][s],bb[i][s]
            for field in ("candidateId","seed","hasMeasurement","gain",
                          "correctedDifferenceEstimate","effectiveUnflagged","sourceUnflagged"):
                if isinstance(x,dict) and isinstance(y,dict):
                    if not approx(x.get(field),y.get(field)):
                        problems.append(f"source-measurement-changed:{ident}:{s}:{field}")
                else:
                    problems.append(f"source-measurement-unavailable:{ident}:{s}")
    return problems


def serial(v):
    if v is None:return []
    if isinstance(v,list):return v
    return [v]


def validate_science(m,subject,candidate,seed):
    errors=[]
    z=m.get("science",{})
    if not isinstance(z,dict) or z.get("status")!="complete":
        return {"complete":False,"errors":["science-not-complete:"+str(z.get("reason")) if isinstance(z,dict) else "science-absent"]}
    dec=serial(z.get("decompositions"))
    if not dec:
        errors.append("no-fitted-decomposition")
    topology=[]
    for i,comp in enumerate(dec):
        chans=comp.get("channelLabels",[])
        mix=array(comp.get("mixing",[]))
        nComp=comp.get("nComponents")
        if mix.ndim!=2 or mix.shape!=(comp.get("nChannels"),nComp) or len(chans)!=len(set(chans)) or len(chans)!=mix.shape[0] or not np.isfinite(mix).all():
            errors.append(f"bad-mixing-dimensions:{i}")
            continue
        removed=[int(v) for v in serial(comp.get("removed"))]
        if len(set(removed))!=len(removed) or any(v<1 or v>nComp for v in removed):
            errors.append(f"removed-components-misaligned:{i}")
        probs=array(comp.get("classProbabilities",[]))
        classes=serial(comp.get("classNames",[]))
        if probs.size>0 and (probs.ndim!=2 or probs.shape!=(nComp,len(classes))
            or not np.isfinite(probs).all()):
            errors.append(f"iclabel-probability-misaligned:{i}")
        topology.append({"labels":chans,"mixing":mix,"removed":removed,
                         "nComp":nComp,"classNames":classes,"probabilities":probs.tolist()})
    times=array(z.get("timeMs",[]))
    conditions=serial(z.get("conditions"))
    if times.ndim!=1 or len(times)<30 or not np.all(np.diff(times)>0) or len(conditions)!=2:
        errors.append("invalid-full-waveform-timebase")
        return {"complete":False,"errors":errors,"topology":topology}
    if z.get("channel")!="Pz":
        errors.append("roi-not-pz")
    if len(conditions)==2:
        waveforms=[]
        for j,cond in enumerate(conditions):
            name=cond.get("name")
            ids=serial(cond.get("referenceIds"))
            retained=serial(cond.get("retainedIds"))
            wave=array(cond.get("meanWaveform",[]))
            scores=array(cond.get("trialScores",[]))
            weights=array(cond.get("trialWeights",[]))
            flags=serial(cond.get("referenceFlagged"))
            if wave.shape!=times.shape or not np.isfinite(wave).all():
                errors.append(f"waveform-missing:{j}")
            if len(ids)!=len(set(ids)) or scores.shape!=(len(ids),) or weights.shape!=(len(ids),) or len(flags)!=len(ids):
                errors.append(f"trial-identity-alignment:{j}")
                continue
            good=np.isfinite(scores)&np.isfinite(weights)&(weights>0)
            expected=[v for v,g in zip(ids,good) if g]
            if expected!=retained or np.count_nonzero(good)!=cond.get("usedCount"):
                errors.append(f"trial-retention-mismatch:{j}")
            if not good.any():
                errors.append(f"missing-positive-weight-trials:{j}")
                continue
            neff=weights[good].sum()**2/np.sum(weights[good]**2)
            if not approx(neff,cond.get("effectiveN")):
                errors.append(f"effective-support-mismatch:{j}")
            p3=(times>=300-1e-7)&(times<=600+1e-7)
            if not p3.any() or not approx(wave[p3].mean(),cond.get("originalEstimate")):
                errors.append(f"p3-window-waveform-does-not-reproduce-score:{j}")
            point=np.sum(scores[good]*weights[good])/np.sum(weights[good])
            if not approx(point,cond.get("windowMean")):
                errors.append(f"trial-score-waveform-inconsistent:{j}")
            waveforms.append({"condition":name,"waveform":wave.tolist(),
                              "scores":scores[good],"weights":weights[good],
                              "positiveCount":int(good.sum())})
        if len(waveforms)!=2:
            return {"complete":False,"errors":errors,"topology":topology}
        gain=m.get("gain")
        if not valid(gain) or gain<=0:
            errors.append("missing-finite-gain")
            return {"complete":False,"errors":errors,"topology":topology}
        diff=(array(waveforms[0]["waveform"])-array(waveforms[1]["waveform"]))/gain
        p3=(times>=300-1e-7)&(times<=600+1e-7)
        if not approx(float(diff[p3].mean()),m.get("correctedDifferenceEstimate")):
            errors.append("gain-corrected-p3-waveform-mismatch")
        rng=np.random.RandomState(BOOT_SEED+subject+1000*seed+candidate)
        draws=[]
        for row in waveforms:
            x,w=row["scores"],row["weights"]
            if len(x)<2:
                errors.append("bootstrap-insufficient-positive-trials")
                break
            idx=rng.randint(0,len(x),(BOOT_DRAWS,len(x)))
            draws.append(np.sum(w[idx]*x[idx],axis=1)/np.sum(w[idx],axis=1))
        boot_se=None
        if len(draws)==2:
            boot=(draws[0]-draws[1])/gain
            boot_se=float(boot.std(ddof=1))
            if not valid(boot_se) or boot_se<=0:
                errors.append("bootstrap-sme-invalid")
        return {"complete":not errors,"errors":errors,
                "bootstrapSE":boot_se, "topology":topology,
                "timeMs":times.tolist(),"p3DifferenceWave":diff.tolist(),
                "positiveTrials":[r["positiveCount"] for r in waveforms]}
    return {"complete":False,"errors":errors,"topology":topology}


def match_components(base,test):
    # Align channel labels, not original component IDs; IC sign is arbitrary.
    a,b=base["topology"],test["topology"]
    if len(a)!=1 or len(b)!=1:
        return {"available":False,"reason":"multiple-or-missing-decompositions"}
    a,b=a[0],b[0]
    if set(a["labels"])!=set(b["labels"]) or len(a["labels"])<4:
        return {"available":False,"reason":"incompatible-electrode-montage"}
    ix=[b["labels"].index(v) for v in a["labels"]]
    A=a["mixing"]; B=b["mixing"][ix,:]
    A=A-A.mean(axis=0,keepdims=True)
    B=B-B.mean(axis=0,keepdims=True)
    da=np.linalg.norm(A,axis=0);db=np.linalg.norm(B,axis=0)
    if np.any(da<1e-12) or np.any(db<1e-12):
        return {"available":False,"reason":"degenerate-mixing-topography"}
    C=np.abs((A/da).T@(B/db))
    ii,jj=linear_sum_assignment(-C)
    matched=C[ii,jj]
    ra=set(a["removed"]); rb=set(b["removed"])
    rma={int(i) for i,k in enumerate(ii) if int(k)+1 in ra}
    rmb={int(i) for i,k in enumerate(jj) if int(k)+1 in rb}
    union=len(rma|rmb);jaccard=1. if union==0 else len(rma&rmb)/union
    return {"available":True,"matchedComponents":len(ii),
            "unmatchedComponents":abs(A.shape[1]-B.shape[1]),
            "medianAbsoluteTopographicCorrelation":float(np.median(matched)),
            "minimumMatchedTopographicCorrelation":float(np.min(matched)),
            "removedComponentMatchedJaccard":float(jaccard),
            "originalRemovedCount":len(ra),"otherSeedRemovedCount":len(rb)}


def evaluate(new,prior):
    subject=new["subject"]
    if subject not in FROZEN_DEVELOPMENT_IDS:
        raise ValueError("Out-of-scope development participant")
    baseline_errors=validate_identity(new,prior)
    c=new["icaSeedCheck"]
    ids=c["ids"];n=len(ids)
    measures=grid(c["candidateMeasurements"],n,"candidateMeasurements")
    ranks=grid(c["rankStatus"],n,"rankStatus")
    findings=[]
    chosen=[]
    for candidate in FROZEN_DEVELOPMENT_IDS[subject]:
        if candidate not in ids:
            baseline_errors.append(f"predeclared-candidate-missing:{candidate}")
            continue
        idx=ids.index(candidate)
        per_seed=[]
        for seed in range(SEEDS):
            z=measures[idx][seed]
            if not isinstance(z,dict):
                per_seed.append({"seed":seed,"scienceValid":False,"errors":["measurement-absent"]})
                continue
            audit=validate_science(z,subject,candidate,seed)
            if ranks[idx][seed]=="feasible" and not audit.get("complete"):
                findings.append(f"feasible-science-capture-invalid:{candidate}:{seed}:{audit['errors']}")
            per_seed.append({"seed":seed,"scienceValid":audit.get("complete",False),
                             "bootstrapSE":audit.get("bootstrapSE"),
                             "P3":z.get("correctedDifferenceEstimate"),
                             "positiveTrials":audit.get("positiveTrials"),
                             "errors":audit.get("errors",[])})
            # Store in-memory only: no raw matrices in compact analysis output.
            per_seed[-1]["_audit"]=audit
        comparisons=[]
        for seed in range(1,SEEDS):
            x,y=per_seed[0]["_audit"],per_seed[seed]["_audit"]
            if x.get("complete") and y.get("complete"):
                comparisons.append({"seed":seed,**match_components(x,y),
                    "fullWaveformRmsDifference":float(np.sqrt(np.mean(
                    (array(x["p3DifferenceWave"])-array(y["p3DifferenceWave"]))**2)))})
            else:
                comparisons.append({"seed":seed,"available":False,
                                    "reason":"scientific-snapshot-incomplete"})
        for row in per_seed: row.pop("_audit",None)
        chosen.append({"candidateId":candidate,"startReports":per_seed,
                       "componentMatchesVsSeed0":comparisons})
    return {"subject":subject,"identityErrors":baseline_errors,"scienceFailures":findings,
        "predeclaredCandidateScience":chosen,
        "unalteredFrozenRankAndEstimates":not baseline_errors,
        "allPredeclaredFeasibleScienceComplete":not findings,
        "original20SubjectVerdict":"PARTIAL",
        "candidateSelectionModified":False,"newParticipants":False,
        "promotionAuthorized":False}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--new",type=Path,required=True)
    p.add_argument("--previous",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    new=json.loads(args.new.read_text())
    previous=json.loads(args.previous.read_text())
    result=evaluate(new,previous)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,allow_nan=False))
    print("PR51_ICA_SCIENCE_RECORD",json.dumps({
        "subject":result["subject"],
        "sourceIdentityMatched":result["unalteredFrozenRankAndEstimates"],
        "scientificCaptureComplete":result["allPredeclaredFeasibleScienceComplete"],
        "identityErrors":result["identityErrors"],
        "scienceFailures":result["scienceFailures"],
        "candidates":[v["candidateId"] for v in result["predeclaredCandidateScience"]]}))
    if result["identityErrors"] or result["scienceFailures"]:
        raise SystemExit("Scientific recorder validation failed closed")


if __name__=="__main__":
    main()
