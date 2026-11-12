#!/usr/bin/env python3
"""Independent P300 ICA triplet-by-triplet reproducibility validator.

Groups use distinct random ICA starts {0,1,2} and {3,4,5}; same EEG
source is shared, so this is NOT independent-participant replication.
Science gates are the original feasible status, source HMM >=5 effective
unflagged epochs, finite same-stratum objectives, same weighted/unweighted
estimand family within each group's 3 starts.
"""
import argparse
import json
import math
from pathlib import Path


def finite(value):
    return isinstance(value,(float,int)) and not isinstance(value,bool) and math.isfinite(value)


def check_group(rows, group, rec, slots, group_number):
    expected_seeds=[3*group_number+i for i in range(3)]
    if group["seedNumbers"] != expected_seeds or group["independentGroup"] != group_number+1:
        raise ValueError("group seed indices overlap or were reordered")
    baseline=rec[slots[0]]
    by_id={z["candidateId"]:z for z in rows}
    family=(by_id[baseline]["estimationFamily"][slots[0]]
            if baseline in by_id else "unavailable")
    if group["baselineId"] != baseline or group["baselineFamily"] != family:
        raise ValueError("baseline-recipe P300 estimator identity differs from seed group")
    eligible=[]
    if family in ("weighted","unweighted"):
        for row in rows:
            f=[row["estimationFamily"][j] for j in slots]
            statuses=[row["status"][j] for j in slots]
            strata=[row["stratum"][j] for j in slots]
            support=[row["supportPass"][j] for j in slots]
            objectives=[row["objective"][j] for j in slots]
            if (all(v==family for v in f) and all(v=="feasible" for v in statuses)
                    and all(v is True for v in support) and len(set(strata))==1
                    and all(finite(v) and v >= 0 for v in objectives)):
                eligible.append((row["candidateId"],max(objectives),strata[0]))
    groups={z[2] for z in eligible}
    if group["status"]=="selected":
        if not eligible or len(groups)!=1:
            raise ValueError("selected group has no comparable scientifically valid candidate")
        best=min(z[1] for z in eligible)
        winners=[z[0] for z in eligible if abs(z[1]-best) <= 1e-12*max(1,abs(best))]
        if len(winners)!=1 or group["candidateId"] != winners[0]:
            raise ValueError("selected group is not unique all-seed feasible minimax")
        if len(group.get("seedMeasurements",[]))!=3:
            raise ValueError("missing independently re-evaluated triplet")
        for j,m in zip(expected_seeds,group["seedMeasurements"]):
            if (m["candidateId"]!=group["candidateId"] or m["seed"]!=j
                    or m["estimationFamily"]!=family):
                raise ValueError("selected candidate, seed or estimator family changed")
            if not (m["hasMeasurement"] and
                    isinstance(m["effectiveUnflagged"],list) and
                    len(m["effectiveUnflagged"])==2 and
                    all(finite(v) and v>=5 for v in m["effectiveUnflagged"])):
                raise ValueError("P300 selected seed has insufficient HMM-unflagged support")
            if not all(finite(m[v]) for v in ("gain","differenceEstimate","correctedDifferenceEstimate")):
                raise ValueError("P300 measurement/gain missing")
            if m["gain"]<=0 or not math.isclose(
                    m["correctedDifferenceEstimate"],m["differenceEstimate"]/m["gain"],
                    abs_tol=1e-8,rel_tol=0):
                raise ValueError("corrected seed P300 is not calibrated")
    elif group["status"]=="abstain":
        if group.get("candidateId") is not None or group.get("seedMeasurements"):
            raise ValueError("scientific abstention cannot contain a selected pipeline")
        if eligible and len(groups)==1:
            best=min(z[1] for z in eligible)
            ties=[z for z in eligible if abs(z[1]-best)<=1e-12*max(1,abs(best))]
            if len(ties)==1:
                raise ValueError("abstained despite one eligible minimum")
    else:
        raise ValueError("unknown seed group outcome")
    return {
        "seeds":expected_seeds,
        "baseline":baseline,
        "estimationFamily":family,
        "scientificallyEligibleCandidates":len(eligible),
        "status":group["status"],
        "candidate":group.get("candidateId"),
        "correctedP3":[m["correctedDifferenceEstimate"] for m in group.get("seedMeasurements",[])],
        "originalScientificGatesUnchanged":True,
    }


def validate(result):
    if not result.get("completed") or result.get("errors") or result.get("notRun",0)!=0:
        raise ValueError("incomplete real-P3 source run")
    if result.get("searchOptions",{}).get("icaSeeds")!=5:
        raise ValueError("six distinct ICA starts were not used")
    data=result.get("icaSeedCheck")
    if not isinstance(data,dict):
        raise ValueError("missing ICA run evidence")
    rows=data.get("candidateSeedAudit")
    if not isinstance(rows,list) or not rows:
        raise ValueError("no candidate/seed source-aligned audit")
    ids=data.get("ids");recommended=data.get("recommended");ms=data.get("measurements")
    if not (isinstance(ids,list) and isinstance(recommended,list) and isinstance(ms,list)
            and len(ids)==len(rows)==len(data.get("objective",[]))
            and len(recommended)==len(ms)==6 and len(set(ids))==len(ids)):
        raise ValueError("source candidate/seed shape mismatch")
    if data.get("seeds")!=[1,2,3,4,5]:
        raise ValueError("ICA seed identities are not the planned five alternate starts")
    for i,row in enumerate(rows):
        if row["candidateId"]!=ids[i] or row["objective"]!=data["objective"][i]:
            raise ValueError("candidate rows re-ordered or objective altered")
        for field in ("status","stratum","supportPass","objective","estimationFamily"):
            if not isinstance(row[field],list) or len(row[field])!=6:
                raise ValueError("flattened MATLAB output or insufficient ICA starts")
        if not all(type(x) is bool for x in row["supportPass"]):
            raise ValueError("seed scientific support is not a boolean")
    first=ms[0]
    if not(first["hasMeasurement"] and finite(first["correctedDifferenceEstimate"]) and
            math.isclose(first["correctedDifferenceEstimate"],result["erp"]["differenceAmplitude"],
                         abs_tol=1e-8,rel_tol=0)):
        raise ValueError("first-seed calibrated P300 does not match unchanged original scientific report")
    pair=data.get("independentGroups")
    if not isinstance(pair,dict) or pair.get("sourceSeedLabels")!=list(range(6)):
        raise ValueError("six ICA starts were not partitioned with explicit provenance")
    if pair.get("originalSearchRankingUnchanged") is not True or pair.get("isIndependentHeldoutStudy") is not False:
        raise ValueError("misidentified original ranking or heldout context")
    group_records=pair.get("groups")
    if not isinstance(group_records,list) or len(group_records)!=2:
        raise ValueError("expected two disjoint independently selected triplets")
    group_out=[check_group(rows,g,recommended,list(range(3*i,3*(i+1))),i)
               for i,g in enumerate(group_records)]
    yes=all(g["status"]=="selected" for g in group_records)
    same=yes and group_records[0]["candidateId"]==group_records[1]["candidateId"]
    expected=("same-across-disjoint-triplets" if same else
              "different-across-disjoint-triplets" if yes else
              "at-least-one-group-abstained")
    if pair["sameCandidate"] is not same or pair["status"] != expected:
        raise ValueError("disjoint-group recommendation consistency incorrectly claimed")
    return {"subject":result["subject"],"originalSeeds":recommended,
            "groups":group_out,"sameAcrossDisjointIcaTriplets":same,
            "status":expected,"allOriginalScientificGatesPreserved":True,
            "heldoutParticipantEEGOpened":False,"newFrozenStudy":False}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--input",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args()
    result=validate(json.loads(a.input.read_text()))
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2))
    print("PR51_DISJOINT_ICA_SEED_RESULT",json.dumps(result))


if __name__=="__main__":
    main()
