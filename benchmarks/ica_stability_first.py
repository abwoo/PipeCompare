#!/usr/bin/env python3
"""Exploratory stability-FIRST versus minimax ICA recipe comparison.

Training seeds 0/1/2 exclusively choose a pipeline. Unseen seeds 3/4/5
are read strictly after choice. Already-observed subjects 11/15/21 are
development evidence, NOT a confirmatory participant-level validation.
"""
import argparse
import json
import math
from pathlib import Path
from ica_feasibility_shadow import shadow
from ica_support_integrity import audit as support_audit
from validate_ica_unseen_seeds import inspect, grid, valid

TRAIN = (0,1,2)
UNSEEN = (3,4,5)

def training_only(study):
    c=study["icaSeedCheck"]
    # Fresh three-seed projection, not the six-seed "same" flag.
    reduced=dict(study)
    d=dict(c)
    n=len(c["ids"])
    for k in ("rankStatus","strata","candidateMeasurements"):
        d[k]=[row[:3] for row in grid(c[k],n,k)]
    d["objective"]=[row[:3] for row in c["objective"]]
    d["recommended"]=c["recommended"][:3]
    d["same"]=len(set(d["recommended"]))==1
    reduced["icaSeedCheck"]=d
    return shadow(reduced)

def choose(training):
    eligible=training["eligible"]
    if not eligible:
        return None, "no-common-feasible-candidate"
    # Stability is the primary objective here; original numerical error
    # acts ONLY as a tie breaker. Do not tune a P3 range threshold.
    best=min(a["p3Range"] for a in eligible)
    tied=[a for a in eligible if a["p3Range"]==best]
    min_cost=min(a["maxObjective"] for a in tied)
    tied=[a for a in tied if a["maxObjective"]==min_cost]
    if len(tied)!=1:
        return None,"unresolved-stability-and-objective-tie"
    return tied[0]["id"],"unique-minimum-training-p3-range"

def unseen(study,ident,training):
    if ident is None:
        return {"candidateId":None,"allUnseenFeasible":False,
                "p3Range":None,"failureReasons":["training-abstention"]}
    c=study["icaSeedCheck"];ids=c["ids"];n=len(ids)
    if ident not in ids:
        raise ValueError("Chosen candidate absent from source")
    strata=grid(c["strata"],n,"strata")
    rank=grid(c["rankStatus"],n,"rankStatus")
    ms=grid(c["candidateMeasurements"],n,"candidateMeasurements")
    original=c["recommended"][0]
    if original not in ids: raise ValueError("Baseline ID absent")
    baseline_stratum=strata[ids.index(original)][0]
    i=ids.index(ident)
    observations=[]
    for seed in UNSEEN:
        m=ms[i][seed]
        p3=m.get("correctedDifferenceEstimate") if isinstance(m,dict) else None
        gain=m.get("gain") if isinstance(m,dict) else None
        eff=m.get("effectiveUnflagged") if isinstance(m,dict) else None
        original_n=m.get("sourceUnflagged") if isinstance(m,dict) else None
        problems=[]
        if rank[i][seed]!="feasible":problems.append("rank-not-feasible")
        if strata[i][seed]!=baseline_stratum:problems.append("stratum-drift")
        if not valid(c["objective"][i][seed]):problems.append("nonfinite-objective")
        if not (isinstance(m,dict) and m.get("candidateId")==ident
                and m.get("seed")==seed and m.get("hasMeasurement")
                and valid(p3) and valid(gain) and gain>0):
            problems.append("calibrated-measurement-invalid")
        if not (isinstance(eff,list) and len(eff)==2
                and isinstance(original_n,list) and len(original_n)==2
                and all(valid(e) and valid(s) and s>=5 and
                        e>=5-1e-9 and e<=s+1e-8 for e,s in zip(eff,original_n))):
            problems.append("effective-versus-physical-support-failed")
        observations.append({"seed":seed,"p3":p3,"passed":not problems,
                             "reasons":problems})
    feasible=all(o["passed"] for o in observations)
    values=[o["p3"] for o in observations]
    span=max(values)-min(values) if all(valid(v) for v in values) else None
    return {"candidateId":ident,"allUnseenFeasible":feasible,
            "p3Range":span,"unseenSeeds":observations}

def compare(study):
    if not support_audit(study)["supportIntegrityPassed"]:
        raise ValueError("Physical/effective trial bound integrity failed")
    historical=inspect(study)
    training=training_only(study)
    stability_id,reason=choose(training)
    baseline=historical["originalBaseline"]
    minimax=historical["shadowUnseen"]
    new=unseen(study,stability_id,training)
    return {"subject":study["subject"],
            "trainingSeeds":list(TRAIN),"unseenSeeds":list(UNSEEN),
            "originalId":baseline["candidateId"],
            "minimaxId":minimax["candidateId"],
            "stabilityFirstId":stability_id,
            "selectionReason":reason,
            "eligibleTrainingCandidates":training["eligibleCount"],
            "baselineUnseenRange":baseline.get("unseenP3Range"),
            "minimaxUnseenRange":minimax.get("unseenP3Range"),
            "stabilityFirstUnseenRange":new["p3Range"],
            "baselineUnseenFeasible":baseline["validAllThree"],
            "minimaxUnseenFeasible":minimax["validAllThree"],
            "stabilityFirstUnseenFeasible":new["allUnseenFeasible"],
            "stabilityFirstObservations":new["unseenSeeds"],
            "newParticipantValidation":False,
            "scientificPromotionAuthorized":False}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--evidence",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    files=sorted(args.evidence.rglob("result.json"))
    if len(files)!=3:raise SystemExit(f"Expected 3 archived records, got {len(files)}")
    rows=[compare(json.loads(file.read_text())) for file in files]
    if sorted(r["subject"] for r in rows)!=[11,15,21]:
        raise SystemExit("Unexpected participant IDs")
    rows.sort(key=lambda r:r["subject"])
    out={"study":"stability-first exploratory ICA fit-selection comparison",
         "rawRun":43169562528,"fitNewEEG":False,"newSubjects":False,
         "original20SubjectVerdict":"PARTIAL",
         "results":rows,"productionPromotionAuthorized":False,
         "freeze14SubjectValidationAuthorized":False}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,indent=2,allow_nan=False))
    for x in rows:
        print("PR51_ICA_STABILITY_FIRST",json.dumps({
            k:x[k] for k in ("subject","originalId","minimaxId","stabilityFirstId",
                "selectionReason","baselineUnseenRange","minimaxUnseenRange",
                "stabilityFirstUnseenRange","stabilityFirstUnseenFeasible")}))
    print("PR51_ICA_STABILITY_FIRST_SCOPE",json.dumps({
        "newSubjects":False,"productionPromotionAuthorized":False,
        "freeze14SubjectValidationAuthorized":False}))

if __name__=="__main__":
    main()
