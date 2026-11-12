#!/usr/bin/env python3
"""Fail-closed multi-metric ICA research screen on archived development data.

No ICA fitting, no Rank mutation, no hidden use of starts 3-5 in selection.
The scientific/precision criteria and thresholds are frozen in
docs/PR51_ICA_MULTIMETRIC_PROTOCOL_20261112.md.
"""
import argparse
import json
import math
from pathlib import Path

from validate_ica_unseen_seeds import grid, valid

COHORT = {11: [12, 25], 15: [5, 2, 6], 21: [25, 4]}
RAW_RUN = "43169562528"
RAW_HEAD = "5c9388562974c8305da95cb59e7249982788d589"
TRAIN = (0, 1, 2)
AUDIT = (3, 4, 5)
FLOOR = 5.
CORR_FLOOR = .90
PRECISION_RATIO = .98
EPSILON = 1e-6


def data_error(label):
    raise ValueError(label)


def metrics(starts, matches, seed_set):
    selected = [starts[s] for s in seed_set]
    m = {x["seed"]: x for x in matches}
    if any(s != 0 and s not in m for s in seed_set):
        data_error("Missing component-matching science")
    corr = [m[s]["minimumMatchedTopographicCorrelation"] for s in seed_set if s != 0]
    jacc = [m[s]["removedComponentMatchedJaccard"] for s in seed_set if s != 0]
    wave = [m[s]["fullWaveformRmsDifference"] for s in seed_set if s != 0]
    p3 = [v["P3"] for v in selected]
    se = [v["bootstrapSE"] for v in selected]
    for value in [*corr, *jacc, *wave, *p3, *se]:
        if not valid(value):
            data_error("Nonfinite topographic/ERP/precision readout")
    if any(v <= 0 for v in se):
        data_error("Nonpositive bootstrap SE")
    if any(not m[s].get("available", False) for s in seed_set if s != 0):
        data_error("Component match unavailable")
    return {
        "p300Range": max(p3) - min(p3),
        "meanBootstrapSE": sum(se)/len(se),
        "maxBootstrapSE": max(se),
        "minimumRemovedJaccard": min(jacc) if jacc else None,
        "maxFullWaveformRMS": max(wave) if wave else None,
        "minimumMatchedTopographyCorrelation": min(corr) if corr else None,
    }


def source_gate(source, candidate, seeds, base_stratum):
    c = source["icaSeedCheck"]
    ids = c["ids"]
    if candidate not in ids:
        return ["candidate-not-in-original-search"]
    n = len(ids)
    i = ids.index(candidate)
    status = grid(c["rankStatus"], n, "rankStatus")
    strata = grid(c["strata"], n, "strata")
    measure = grid(c["candidateMeasurements"], n, "candidateMeasurements")
    if len(c["objective"]) != n or any(len(x)!=6 for x in c["objective"]):
        data_error("Bad original objective shape")
    failures=[]
    for s in seeds:
        prefix=f"start-{s}:"
        if status[i][s] != "feasible": failures.append(prefix+"original-rank-not-feasible")
        if strata[i][s] != base_stratum: failures.append(prefix+"baseline-stratum-mismatch")
        if not valid(c["objective"][i][s]): failures.append(prefix+"nonfinite-original-objective")
        z = measure[i][s]
        if not isinstance(z,dict) or z.get("candidateId") != candidate or z.get("seed") != s:
            failures.append(prefix+"candidate-seed-misalignment")
            continue
        if not z.get("hasMeasurement") or not valid(z.get("gain")) or z["gain"] <= 0 or not valid(z.get("correctedDifferenceEstimate")):
            failures.append(prefix+"invalid-calibrated-p3")
        eff, phys = z.get("effectiveUnflagged"), z.get("sourceUnflagged")
        if not isinstance(eff,list) or not isinstance(phys,list) or len(eff)!=2 or len(phys)!=2:
            failures.append(prefix+"missing-two-condition-support")
        else:
            for condition,(e,p) in enumerate(zip(eff,phys)):
                if (not valid(e) or not valid(p) or p < FLOOR or e+1e-9 < FLOOR
                    or e < 0 or e > p+1e-8 or p != int(p)):
                    failures.append(prefix+f"invalid-condition-{condition}-support")
    return failures


def describe(choice, source, baseline_stratum, seeds):
    starts = choice["startReports"]
    matches = choice["componentMatchesVsSeed0"]
    if (len(starts)!=6 or sorted(s["seed"] for s in starts) != list(range(6))
            or len(matches)!=5 or sorted(m["seed"] for m in matches)!=[1,2,3,4,5]):
        data_error("Incomplete or duplicated six-start science snapshot")
    starts=sorted(starts,key=lambda x:x["seed"])
    failures=source_gate(source, choice["candidateId"], seeds, baseline_stratum)
    for s in seeds:
        if not starts[s].get("scienceValid") or starts[s].get("errors"):
            failures.append(f"start-{s}:science-capture-invalid")
    try:
        report=metrics(starts,matches,seeds)
    except (ValueError,KeyError,TypeError,IndexError) as exc:
        failures.append(f"science-metrics:{exc}")
        report=None
    if report and report["minimumMatchedTopographyCorrelation"] is not None:
        if report["minimumMatchedTopographyCorrelation"]+1e-12 < CORR_FLOOR:
            failures.append("matched-topography-below-0.90")
    return {"candidateId":choice["candidateId"],"eligible":not failures,
            "failures":failures,"metrics":report}


def compare_to_original(other, original):
    if not other["eligible"]:
        return False, ["candidate-ineligible"]+other["failures"]
    a=other["metrics"]; b=original["metrics"]
    checks={
        "removal-consistency-worse": a["minimumRemovedJaccard"]+EPSILON < b["minimumRemovedJaccard"],
        "full-waveform-drift-worse": a["maxFullWaveformRMS"] > b["maxFullWaveformRMS"]+EPSILON,
        "p300-amplitude-drift-worse": a["p300Range"] > b["p300Range"]+EPSILON,
        "mean-bootstrap-se-gain-less-than-2pct": a["meanBootstrapSE"] > b["meanBootstrapSE"]*PRECISION_RATIO+EPSILON,
        "worst-bootstrap-se-worse": a["maxBootstrapSE"] > b["maxBootstrapSE"]+EPSILON,
    }
    errors=[reason for reason,failed in checks.items() if failed]
    return not errors, errors


def evaluate(source, science):
    sid=source.get("subject")
    if sid not in COHORT or science.get("subject")!=sid:
        data_error("Unexpected or mismatched participant")
    if (not source.get("completed") or source.get("errors") or
        not science.get("unalteredFrozenRankAndEstimates") or
        science.get("identityErrors") or science.get("scienceFailures") or
        not science.get("allPredeclaredFeasibleScienceComplete")):
        data_error("Incomplete, changed or failed source evidence")
    # Identity may be encoded as integer/string in source JSON.
    if str(source.get("workflowRun")) != RAW_RUN or source.get("executionCommit") != RAW_HEAD:
        data_error("Raw EEG execution commit/run identity mismatch")
    c=source["icaSeedCheck"]
    baseline_id=COHORT[sid][0]
    if c["recommended"][0]!=baseline_id:
        data_error("Original Standard baseline ID does not match frozen roster")
    ids=c["ids"]
    base_stratum=grid(c["strata"],len(ids),"strata")[ids.index(baseline_id)][0]
    choices={x["candidateId"]:x for x in science["predeclaredCandidateScience"]}
    if set(choices)!=set(COHORT[sid]):
        data_error("Science candidate roster changed")
    candidates={i:describe(choices[i],source,base_stratum,TRAIN) for i in COHORT[sid]}
    orig=candidates[baseline_id]
    results=[]
    if not orig["eligible"]:
        # A failed baseline is an outcome, not a reason to hide independently
        # observed candidate eligibility. Never choose a candidate here.
        for ident in COHORT[sid][1:]:
            x=candidates[ident]
            results.append({"candidateId":ident,"passesTrainingConstraints":False,
                            "failures":["ineligible-original-comparator",*x["failures"]],
                            "trainingMetrics":x["metrics"]})
        decision={"kind":"abstain","candidateId":None,"reason":"ineligible-original-baseline"}
    else:
        eligible=[]
        for ident in COHORT[sid][1:]:
            x=candidates[ident]
            allowed,why=compare_to_original(x,orig)
            results.append({"candidateId":ident,"passesTrainingConstraints":allowed,
                            "failures":why,"trainingMetrics":x["metrics"]})
            if allowed:
                k=ids.index(ident)
                worst=max(c["objective"][k][s] for s in TRAIN)
                eligible.append((worst,ident))
        if not eligible:
            decision={"kind":"retain-original","candidateId":baseline_id,
                      "reason":"no-multimetric-noninferior-alternative"}
        else:
            eligible.sort()
            top=eligible[0][0]
            tied=[ident for cost,ident in eligible if cost==top]
            decision=({"kind":"shadow-switch","candidateId":tied[0],
                       "reason":"unique-multimetric-screen-and-worst-objective"}
                      if len(tied)==1 else
                      {"kind":"abstain","candidateId":None,
                       "reason":"unresolved-original-objective-tie"})
    # Independent-start numbers are *only* computed after the fixed decision.
    assessed={i:describe(choices[i],source,base_stratum,AUDIT) for i in COHORT[sid]}
    baseline_assess=assessed[baseline_id]
    picked_id=decision["candidateId"]
    outside={"baseline":baseline_assess,
             "candidate":assessed.get(picked_id) if picked_id is not None else None,
             "passesProtectedNonregression":None,"failureReasons":[]}
    if picked_id is not None and baseline_assess["eligible"] and assessed[picked_id]["eligible"]:
        ok,problems=compare_to_original(assessed[picked_id],baseline_assess)
        outside["passesProtectedNonregression"]=ok
        outside["failureReasons"]=problems
    else:
        outside["failureReasons"]=["missing-scientifically-eligible-comparison"]
    return {
        "subject":sid,"originalId":baseline_id,"selectionSeeds":list(TRAIN),
        "assessmentSeeds":list(AUDIT),"decision":decision,
        "originalTraining":orig,"alternativeTraining":results,
        "assessment":outside,
        "originalPR51Verdict":"PARTIAL","productionPromotionAuthorized":False,
        "prospectiveValidationAuthorized":False,
    }


def report(raw_dir,science_dir):
    raw=sorted(raw_dir.rglob("result.json"))
    sci=sorted(science_dir.rglob("ica_science_validated.json"))
    if len(raw)!=3 or len(sci)!=3:
        data_error(f"Exactly three raw and three validated science records required; {len(raw)}/{len(sci)}")
    raws={int((x:=json.loads(p.read_text()))["subject"]):x for p in raw}
    sciences={int((x:=json.loads(p.read_text()))["subject"]):x for p in sci}
    if set(raws)!=set(COHORT) or set(sciences)!=set(COHORT):
        data_error("Duplicated, missing or out-of-scope subjects")
    rows=[evaluate(raws[k],sciences[k]) for k in sorted(COHORT)]
    return {
        "scope":"development only; archived six ICA starts, no new EEG",
        "rawRun":43169562528,"scienceRun":43184910782,
        "trainingSeeds":list(TRAIN),"assessmentSeeds":list(AUDIT),
        "original20SubjectValidPairs":15,"original20SubjectVerdict":"PARTIAL",
        "newParticipantsEvaluated":False,"productionPromotionAuthorized":False,
        "newFrozen14ParticipantRunAuthorized":False,"results":rows,
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--raw",type=Path,required=True)
    p.add_argument("--science",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    result=report(a.raw,a.science)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,allow_nan=False))
    for row in result["results"]:
        print("PR51_MULTI_ICA_SUBJECT",json.dumps({
            "subject":row["subject"],"decision":row["decision"],
            "originalTrainingFailures":row["originalTraining"]["failures"],
            "originalTrainingMetrics":row["originalTraining"]["metrics"],
            "alternativeTraining":[{"candidateId":x["candidateId"],"passes":x["passesTrainingConstraints"],"failures":x["failures"]} for x in row["alternativeTraining"]],
            "assessmentPass":row["assessment"]["passesProtectedNonregression"],
            "assessmentFailures":row["assessment"]["failureReasons"]}))
    print("PR51_MULTI_ICA_PROMOTION",json.dumps({
        "authorized":False,"holdoutOpened":False,"originalPR51Verdict":"PARTIAL"}))


if __name__=="__main__":
    main()
