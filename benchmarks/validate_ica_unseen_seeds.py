#!/usr/bin/env python3
"""Unseen ICA-start stress test: train seeds 0..2, audit starts 3..5.

This is a development-only research comparison on previously observed ERP
CORE subjects 11/15/21. It does NOT test new participants and cannot
confirm the protocol or justify relaxing original scientific eligibility.
"""
import argparse
import json
import math
from copy import deepcopy
from pathlib import Path
from ica_feasibility_shadow import shadow

SELECTION = (0,1,2)
UNSEEN = (3,4,5)
SEEDS = 6
FLOOR = 5.0


def valid(x):
    return isinstance(x,(int,float)) and not isinstance(x,bool) and math.isfinite(x)


def grid(x,n,name):
    if not isinstance(x,list): raise ValueError(f"{name} not JSON array")
    if len(x)==n and all(isinstance(r,list) and len(r)==SEEDS for r in x):
        return x
    if len(x)==n*SEEDS:
        return [[x[s*n+i] for s in range(SEEDS)] for i in range(n)]
    raise ValueError(f"{name} must align to six seeds for {n} candidate IDs")


def inspect(study):
    if not study.get("completed") or study.get("errors"):
        raise ValueError("cannot confirm incomplete real P300 run")
    c=study["icaSeedCheck"]
    ids=c["ids"]
    if len(c["objective"])!=len(ids) or len(c["recommended"])!=SEEDS:
        raise ValueError("ICA score or recommendation seed length mismatch")
    if not all(len(row)==SEEDS for row in c["objective"]):
        raise ValueError("Need six selection and independent ICA-start columns")
    status=grid(c["rankStatus"],len(ids),"rankStatus")
    strata=grid(c["strata"],len(ids),"strata")
    measures=grid(c["candidateMeasurements"],len(ids),"candidateMeasurements")
    first=deepcopy(study)
    d=first["icaSeedCheck"]
    d["objective"]=[r[:3] for r in c["objective"]]
    d["rankStatus"]=[r[:3] for r in status]
    d["strata"]=[r[:3] for r in strata]
    d["candidateMeasurements"]=[r[:3] for r in measures]
    d["recommended"]=c["recommended"][:3]
    chosen=shadow(first)
    target=chosen["shadowCandidate"]
    original=c["recommended"][0]
    if not valid(original) or original not in ids:
        raise ValueError("Original candidate missing")
    # Freeze cross-seed comparison to the selection seed-0 baseline stratum.
    # A feasible Rank.status alone does not establish same-estimand
    # comparability or a finite objective in a separately fitted seed.
    baseline_stratum=strata[ids.index(original)][0]
    def heldout(ident):
        if ident is None:
            return {"candidateId":None,"validAllThree":False,
                    "reason":"training-no-common-feasible-candidate"}
        row=ids.index(ident)
        observations=[]
        for seed in UNSEEN:
            z=measures[row][seed]
            eff=z.get("effectiveUnflagged",[]) if isinstance(z,dict) else []
            p3=z.get("correctedDifferenceEstimate") if isinstance(z,dict) else None
            gain=z.get("gain") if isinstance(z,dict) else None
            # No post-hoc relaxation: test independent fits against the
            # same frozen science gates and baseline-selected stratum.
            reasons=[]
            if status[row][seed]!="feasible":
                reasons.append("rank-not-feasible")
            if strata[row][seed]!=baseline_stratum:
                reasons.append("stratum-drift")
            if not valid(c["objective"][row][seed]):
                reasons.append("nonfinite-objective")
            if not (isinstance(z,dict) and z.get("candidateId")==ident
                    and z.get("seed")==seed):
                reasons.append("misaligned-measurement")
            if not (isinstance(z,dict) and z.get("hasMeasurement")
                    and valid(p3) and valid(gain) and gain>0):
                reasons.append("finite-calibrated-p3-missing")
            if not (isinstance(eff,list) and len(eff)==2 and
                    all(valid(v) and v+1e-9>=FLOOR for v in eff)):
                reasons.append("frozen-effective-trial-floor-not-met")
            observations.append({"seed":seed,"rankStatus":status[row][seed],
                  "stratum":strata[row][seed],
                  "objective":c["objective"][row][seed],
                  "p3":p3,"effectiveUnflagged":eff,
                  "passed":not reasons,"failureReasons":reasons})
        accepted=all(t["passed"] for t in observations)
        values=[t["p3"] for t in observations if valid(t["p3"])]
        return {"candidateId":ident,"validAllThree":accepted,
                "unseenP3Range":max(values)-min(values) if len(values)==3 else None,
                "seeds":observations}
    return {"subject":study["subject"],"selectionSeeds":list(SELECTION),
      "unseenFitSeeds":list(UNSEEN),
      "originalBaseline":heldout(original),
      "trainingShadow":chosen,
      "shadowUnseen":heldout(target),
      "newParticipantsEvaluated":False,
      "original20SubjectVerdictUnchanged":True,
      "scientificPromotionAuthorized":False}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--input",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    a=parser.parse_args()
    result=inspect(json.loads(a.input.read_text()))
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2))
    print("PR51_UNSEEN_SEED_STRESS",json.dumps({
      "subject":result["subject"],
      "original":result["originalBaseline"]["candidateId"],
      "shadow":result["shadowUnseen"]["candidateId"],
      "originalUnseenValid":result["originalBaseline"]["validAllThree"],
      "shadowUnseenValid":result["shadowUnseen"]["validAllThree"],
      "originalUnseenRange":result["originalBaseline"].get("unseenP3Range"),
      "shadowUnseenRange":result["shadowUnseen"].get("unseenP3Range")}))


if __name__=="__main__":
    main()
