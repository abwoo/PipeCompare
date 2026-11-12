#!/usr/bin/env python3
"""Audit cohort predeclaration without opening or downloading held-out EEG."""
import argparse
import json
from pathlib import Path


def validate(plan, frozen):
    if plan["status"]!="PRE-FREEZE-DESIGN-NOT_EXECUTED":
        raise ValueError("Cannot label this preanalysis design a completed freeze")
    if plan["source_manifest"]["status"]!="UNPINNED":
        raise ValueError("An independently checksum-pinned manifest needs a separate freeze review")
    if plan["algorithm"]["status"]!="UNSELECTED":
        raise ValueError("Algorithm must be approved and SHA-pinned outside the cohort draft")
    exc={int(x["subject"]) for x in frozen["source"]["excluded_from_p3"]}
    source={int(x["subject"]) for x in frozen["source"]["files"]}
    eligible=set(range(1,41))-exc
    unseen=sorted(eligible-source)
    if len(eligible)!=34 or len(source)!=20 or len(unseen)!=14:
        raise ValueError("Historical ERP CORE P3 cohorts do not match predeclared 34/20/14")
    if (sorted(plan["previously_used_20"])!=sorted(source) or
        sorted(plan["previously_unseen_14"])!=unseen or
        set(plan["previously_used_20"])&set(plan["previously_unseen_14"])):
        raise ValueError("Development/holdout participant contamination")
    if sorted(plan["excluded_subjects"])!=sorted(exc):
        raise ValueError("Original published P3 exclusions were changed")
    g=plan["frozen_after_algorithm_selection_gates"]
    if g["complete_subjects"]!=14 or g["matched_pair_subjects"]!=14:
        raise ValueError("All 14 independent participants must be reported")
    if g["minimum_effective_unflagged_each_arm_each_condition"]!=5:
        raise ValueError("Original source HMM effective support floor cannot be loosened")
    if g["max_individual_regression_per_control"]>1:
        raise ValueError("Individual SE regressions allowance was loosened")
    if g["bootstrap_draws"]!=10000 or g["bootstrap_seed"]!=510014:
        raise ValueError("Primary participant-level bootstrap is not pinned")
    if not plan["preflight_required"] or not plan["forbidden"]:
        raise ValueError("No fixed study guardrails")
    return {"status":"HOLDOUT_COHORT_DRAFT_VERIFIED_NOT_EXECUTED",
            "development_count":len(source),"holdout_subjects":unseen,
            "source_manifest_status":plan["source_manifest"]["status"],
            "algorithm_status":plan["algorithm"]["status"],
            "outcomes_reviewed":False,"algorithm_frozen":False}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--draft",type=Path,required=True)
    p.add_argument("--frozen-reference",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    result=validate(json.loads(args.draft.read_text()),json.loads(args.frozen_reference.read_text()))
    args.out.parent.mkdir(exist_ok=True,parents=True)
    args.out.write_text(json.dumps(result,indent=2))
    print("PR51_HOLDOUT_COHORT_AUDIT",json.dumps(result))


if __name__=="__main__":
    main()
