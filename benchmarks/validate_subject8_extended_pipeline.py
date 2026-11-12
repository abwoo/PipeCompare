#!/usr/bin/env python3
"""Prospectively scoped subject-8 complete full-record Standard+fixed pilot.

This is one previously evaluated development subject, NOT 14 held-out
people and NOT a new freeze. Scientific pair eligibility and bootstrap
are checked independently against archived trial scores; failed arms
remain explicit. The original 300-second results are never rewritten.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

from compare_erpcore import bootstrap_se, arm_support

SOURCE_RUN_HMM = "43179740278"
EXPECTED_SUBJECT = 8
ARMS = ("erp", "minimal", "common")
FLOOR = 5.0


def finite(x):
    return isinstance(x,(int,float)) and not isinstance(x,bool) and math.isfinite(x)


def analyze(r, support):
    if r.get("subject") != EXPECTED_SUBJECT:
        raise ValueError("Only known development subject 8 is permitted")
    if not (r.get("developmentOnly") is True and r.get("extended460Seconds") is True
            and r.get("originalFrozenVerdict") == "PARTIAL"
            and r.get("promotionAuthorized") is False):
        raise ValueError("Experiment is missing mandatory historical boundaries")
    source = r.get("source",{})
    if not (source.get("subject")==8 and source.get("duration_seconds")==460
            and source.get("original_frozen_duration_seconds")==300
            and source.get("exploratory_duration_override") is True):
        raise ValueError("Wrong extended source duration or protocol boundaries")
    if not (r.get("completed") is True and r.get("errors")==[] and r.get("notRun")==0):
        raise ValueError("Unfinished pipeline search cannot be scientifically interpreted")
    if not (finite(r.get("replayMaxError")) and r["replayMaxError"]<1e-5 and
            finite(r.get("exportMaxError")) and r["exportMaxError"]<1e-5):
        # Replay/export only exists if Standard actually recommends a recipe.
        if r.get("erp",{}).get("available"):
            raise ValueError("Selected Standard not accurately replayable")
    if r.get("protocolSha256")!=source.get("protocol_sha256"):
        raise ValueError("Protocol source identity mismatch")

    h = support
    if not (h.get("subject")==8 and h.get("originalSeconds")==300 and
            h.get("fullSeconds")==460 and h.get("fullWindowUnflaggedFloorMet") is True
            and [v["unflaggedN"] for v in h["baseline"]]==[5,20]
            and [v["unflaggedN"] for v in h["extended"]]==[6,24]):
        raise ValueError("Pinned HMM source support prerequisite not reproduced")
    for key in ("source_set_sha256","source_fdt_sha256","source_events_sha256"):
        alt={"source_set_sha256":"sourceSetSha256",
             "source_fdt_sha256":"sourceFdtSha256",
             "source_events_sha256":"sourceEventsSha256"}[key]
        if source.get(key)!=h.get(alt):
            raise ValueError(f"Mismatch with independently measured HMM source: {key}")
    if r.get("reference",{}).get("available") is not True:
        raise ValueError("Scientific P300 reference missing")
    original_floor = 5
    if original_floor != FLOOR:
        raise ValueError("Frozen floor altered")
    rows=[]
    for name in ARMS:
        arm=r.get(name,{})
        candidate=arm.get("available") is True
        sup=arm_support(r,name)
        reasons=[]
        if not candidate: reasons.append("candidate-unavailable")
        if not isinstance(sup,list) or len(sup)!=2:
            reasons.append("two-condition-support-unavailable")
            sup=[]
        elif [q.get("condition") for q in sup] != ["target","standard"]:
            reasons.append("condition-identity-invalid")
        effective=[]
        for row in sup:
            eff=row.get("effectiveUnflaggedN")
            effective.append(eff if finite(eff) else None)
            if not finite(eff) or eff+1e-9<FLOOR:
                reasons.append(f"effective-support-below-five:{row.get('condition')}")
        for metric in ("differenceAmplitude","gain","reportedDifferenceSE"):
            if not finite(arm.get(metric)) or (metric!="differenceAmplitude" and arm[metric]<=0):
                reasons.append(f"invalid-{metric}")
        se=None
        if not reasons:
            try:
                boot=bootstrap_se(sup,arm["gain"],r["subject"])
                if boot is None or not finite(boot["se"]) or boot["se"]<=0:
                    reasons.append("bootstrap-SE-unavailable")
                elif abs(boot["point"]-arm["differenceAmplitude"])>1e-6:
                    reasons.append("bootstrap-point-does-not-match-MATLAB")
                else:
                    se=boot["se"]
            except (ValueError,TypeError,KeyError,IndexError,ZeroDivisionError) as exc:
                reasons.append(f"bootstrap-error:{type(exc).__name__}")
        rows.append({"arm":name,"available":candidate,
            "effectiveUnflaggedN":effective,
            "referenceN":[q.get("referenceN") for q in sup],
            "selectedGainCorrectedP300":arm.get("differenceAmplitude"),
            "reportedDifferenceSE":arm.get("reportedDifferenceSE"),
            "bootstrapDifferenceSE":se,
            "individuallyScientificallyEligible":not reasons,
            "failureReasons":reasons})
    eligible=all(a["individuallyScientificallyEligible"] for a in rows)
    by={a["arm"]:a for a in rows}
    return {"scope":"single development participant 8; full 460-second P3 recording",
        "pipelineExecutionCommit":r.get("executionCommit"),
        "pipelineWorkflowRun":r.get("workflowRun"),
        "originalHMMRun":SOURCE_RUN_HMM,
        "rawEEGIdentityPinned":True,
        "HMMSourceFloorMet":True,
        "allThreeArmsScientificallyEligible":eligible,
        "pairedPrecisionComparisonJustifiedForThisOneSubject":eligible,
        "bootstrapSEStandardMinusMinimal": (
             by["erp"]["bootstrapDifferenceSE"]-by["minimal"]["bootstrapDifferenceSE"]
             if eligible else None),
        "bootstrapSEStandardMinusCommon": (
             by["erp"]["bootstrapDifferenceSE"]-by["common"]["bootstrapDifferenceSE"]
             if eligible else None),
        "arms":rows,"newHeldoutParticipants":0,
        "new14ParticipantFrozenStudyLaunched":False,
        "statisticalGeneralizationClaimSupported":False,
        "productionPromotionAuthorized":False,
        "original20SubjectVerdict":"PARTIAL",
        "original20SubjectValidPairs":15}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--record",type=Path,required=True)
    p.add_argument("--support",type=Path,required=True)
    p.add_argument("--protocol",type=Path,default=Path("benchmarks/erpcore_protocol.json"))
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    r=json.loads(a.record.read_text())
    original=json.loads(a.support.read_text())
    protocol=json.loads(a.protocol.read_text())
    source=next(e for e in protocol["source"]["files"] if e["subject"]==8)
    if any(r["source"][f"source_{k}_sha256"]!=source[k]["sha256"] for k in ("set","fdt","events")):
        raise ValueError("Archived source not pinned to original frozen manifest")
    if r.get("algorithmCommit")!=protocol["algorithm_commit"]:
        raise ValueError("Original algorithm identity claim mismatched")
    if r.get("measureSha256") != protocol["baseline_measure_sha256"]:
        raise ValueError("Measure source changed")
    digest=hashlib.sha256(a.protocol.read_bytes()).hexdigest()
    if r.get("protocolSha256")!=digest:
        raise ValueError("Original frozen protocol SHA changed")
    result=analyze(r,original)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2))
    for arm in result["arms"]:
        print("PR51_SUB8_EXTENDED_ARM",json.dumps(arm))
    print("PR51_SUB8_EXTENDED_SCIENCE",json.dumps({
        "pairedEligible":result["allThreeArmsScientificallyEligible"],
        "sourceFloorMet":result["HMMSourceFloorMet"],
        "newFrozenStudyLaunched":result["new14ParticipantFrozenStudyLaunched"],
        "productionPromotionAuthorized":result["productionPromotionAuthorized"]}))


if __name__=="__main__":
    main()
