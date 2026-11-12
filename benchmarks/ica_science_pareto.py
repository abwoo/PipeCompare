#!/usr/bin/env python3
"""Research-only multimetric ICA Pareto abstention study.

No weights or post-hoc numerical thresholds. Train seeds 0..2 only;
evaluate held-start 3..5 only AFTER a unique Pareto survivor is chosen.
Previously investigated participants 11/15/21 are development data, not
new confirmatory cases. This routine never touches production Rank.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np

from ica_stability_first import training_only, unseen
from validate_ica_science import validate_science, match_components
from validate_ica_unseen_seeds import inspect, grid

TRAIN=(0,1,2)
UNSEEN=(3,4,5)
EXPECTED={11,15,21}
EPS=1e-12


def finite(x):
    return isinstance(x,(int,float)) and not isinstance(x,bool) and math.isfinite(x)


def dominates(a,b):
    """Four objectives; maximize topology/removed stability, minimize error/drift."""
    keys=("worstBootstrapSE","maxFullWaveRms","minusMinTopologyCorrelation",
          "minusMinRemovalJaccard")
    av=[a[k] for k in keys]
    bv=[b[k] for k in keys]
    if not all(finite(x) for x in av+bv):
        raise ValueError("Nonfinite Pareto evidence")
    return all(x<=y+EPS for x,y in zip(av,bv)) and any(x<y-EPS for x,y in zip(av,bv))


def science_profile(ms, subject, ident, seeds):
    audit=[]
    for s in seeds:
        z=validate_science(ms[s],subject,ident,s)
        if not z["complete"]:
            return None, "missing-scientific-readout:"+str(z["errors"])
        audit.append(z)
    z0=audit[0]
    topo=[]
    waves=[]
    for k in range(1,len(audit)):
        t=match_components(z0,audit[k])
        if not t.get("available"):
            return None,"topographic-matching-unavailable:"+t.get("reason","unknown")
        topo.append(t)
        d=np.asarray(z0["p3DifferenceWave"])-np.asarray(audit[k]["p3DifferenceWave"])
        if d.shape!=np.asarray(z0["p3DifferenceWave"]).shape or not np.all(np.isfinite(d)):
            return None,"full-waveform-timebase-or-value-mismatch"
        waves.append(float(np.sqrt(np.mean(d*d))))
    se=[x["bootstrapSE"] for x in audit]
    if not all(finite(v) and v>0 for v in se):
        return None,"bootstrap-se-unavailable"
    vals=[a["medianAbsoluteTopographicCorrelation"] for a in topo]
    rem=[a["removedComponentMatchedJaccard"] for a in topo]
    return {"worstBootstrapSE":max(se),"meanBootstrapSE":sum(se)/len(se),
            "maxFullWaveRms":max(waves),"minusMinTopologyCorrelation":-min(vals),
            "minusMinRemovalJaccard":-min(rem),
            "topographyComparisons":len(topo),
            "waveformComparisons":len(waves)},None


def train_selection(study):
    subject=study["subject"]
    c=study["icaSeedCheck"];n=len(c["ids"])
    ids=c["ids"]
    data=grid(c["candidateMeasurements"],n,"candidateMeasurements")
    training=training_only(study)
    available={x["id"] for x in training["eligible"]}
    audit=[]
    candidates=[]
    for i,ident in enumerate(ids):
        if ident not in available:
            continue
        p,why=science_profile(data[i],subject,ident,TRAIN)
        if why:
            audit.append({"candidateId":ident,"scienceEligible":False,"reason":why})
            continue
        x={"candidateId":ident,**p}
        candidates.append(x)
        audit.append({"candidateId":ident,"scienceEligible":True})
    frontier=[x for x in candidates if not any(
        y["candidateId"]!=x["candidateId"] and dominates(y,x) for y in candidates)]
    if len(frontier)==1:
        choice=frontier[0]["candidateId"]
        reason="unique-scientific-four-axis-pareto-survivor"
    elif not frontier:
        choice=None
        reason="no-complete-feasible-scientific-candidate"
    else:
        choice=None
        reason="multiobjective-tradeoff-abstain"
    return {"subject":subject,"selectionSeeds":list(TRAIN),
            "originalCommonFeasibleCandidates":len(available),
            "scienceCompleteCandidates":len(candidates),
            "scienceExcluded":audit,
            "frontierCount":len(frontier),
            "frontier":frontier,
            "candidateId":choice,"reason":reason,
            "allCandidateTrainingMetrics":candidates}


def heldout(study,choice,training):
    if choice is None:
        return {"candidateId":None,"validAllUnseenSeeds":False,
                "reason":"training-abstention","unseenP3Range":None,
                "unseenBootstrapSE":None}
    base=unseen(study,choice,training)
    c=study["icaSeedCheck"];data=grid(c["candidateMeasurements"],len(c["ids"]),"candidateMeasurements")
    i=c["ids"].index(choice)
    new=[]
    for seed in UNSEEN:
        z=validate_science(data[i][seed],study["subject"],choice,seed)
        if not z["complete"]:
            return {"candidateId":choice,"validAllUnseenSeeds":False,
                    "reason":"unseen-science-invalid:"+str(z["errors"]),
                    "unseenP3Range":base["p3Range"]}
        new.append(z)
    # Unseen-only P3 range is never a training objective.
    return {"candidateId":choice,"validAllUnseenSeeds":base["allUnseenFeasible"],
            "reason":"all-strict-scientific-gates-pass" if base["allUnseenFeasible"] else "original-rank-or-support-failed",
            "unseenP3Range":base["p3Range"],
            "meanUnseenBootstrapSE":sum(q["bootstrapSE"] for q in new)/len(new),
            "unseenBootstrapSE":[q["bootstrapSE"] for q in new]}


def evaluate(study):
    if study.get("completed") is not True or study.get("errors"):
        raise ValueError("Incomplete archived real EEG evidence")
    if study["subject"] not in EXPECTED:
        raise ValueError("New subject forbidden in development analysis")
    original=inspect(study)
    train=train_selection(study)
    unseen_result=heldout(study,train["candidateId"],train)
    return {"subject":study["subject"],
        "training":train,
        "unseen":unseen_result,
        "baselineId":original["originalBaseline"]["candidateId"],
        "minimaxId":original["shadowUnseen"]["candidateId"],
        "originalBaselineUnseenP3Range":original["originalBaseline"]["unseenP3Range"],
        "minimaxUnseenP3Range":original["shadowUnseen"]["unseenP3Range"],
        "newParticipantsEvaluated":False,
        "original20SubjectVerdict":"PARTIAL",
        "productionPromotionAuthorized":False,
        "open14SubjectConfirmationAuthorized":False}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--evidence",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    files=sorted(a.evidence.rglob("result.json"))
    if len(files)!=3:raise SystemExit("Need exactly three archived science recorder subjects")
    rows=[evaluate(json.loads(f.read_text())) for f in files]
    if {v["subject"] for v in rows}!=EXPECTED:
        raise SystemExit("Duplicate or missing development subject")
    rows.sort(key=lambda x:x["subject"])
    out={"analysis":"four-axis Pareto scientific ICA development experiment",
         "sourceRun":43184910782,
         "noNewEEGFits":True,
         "noNewParticipants":True,
         "original20SubjectVerdict":"PARTIAL",
         "productionPromotionAuthorized":False,
         "freeze14SubjectValidationAuthorized":False,
         "rows":rows}
    a.out.parent.mkdir(exist_ok=True,parents=True)
    a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
    for r in rows:
        print("PR51_ICA_SCIENCE_PARETO",json.dumps({
            "subject":r["subject"],
            "originalId":r["baselineId"],
            "minimaxId":r["minimaxId"],
            "trainCommonFeasible":r["training"]["originalCommonFeasibleCandidates"],
            "scienceComplete":r["training"]["scienceCompleteCandidates"],
            "paretoFrontier":r["training"]["frontierCount"],
            "researchCandidate":r["training"]["candidateId"],
            "reason":r["training"]["reason"],
            "unseenFeasible":r["unseen"]["validAllUnseenSeeds"],
            "unseenP3Range":r["unseen"].get("unseenP3Range"),
            "unseenMeanBootstrapSE":r["unseen"].get("meanUnseenBootstrapSE")}))
    print("PR51_ICA_SCIENCE_PARETO_DECISION",json.dumps({
        "newParticipants":False,"promotion":False,
        "originalVerdict":"PARTIAL","holdoutNotOpened":True}))


if __name__=="__main__":
    main()
