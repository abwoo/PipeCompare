#!/usr/bin/env python3
"""Read-only ERP CORE P3 effective-trial support audit for PR51.

Uses the original 20 frozen, connector-captured RESULT-marker snapshots.
Snapshots are derived summaries, not raw EEG, trial tables or unmodified
GitHub Actions log downloads. No trial recovery is claimed from this audit.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from pathlib import Path

EXPECTED = (1,2,3,4,5,7,8,11,12,13,14,15,16,17,18,19,20,21,22,23)
ALGORITHM = "7151821edf8a5a54c3167d7f18b66c530911b568"
PROTOCOL = "bb65c9a522e9f45210e5ee941a669630b75587a5255664aa0f1402b572e89c6b"
RUN = "42809789847"
ARMS = ("standard","minimal","common")
PRIMARY = (4,8,17,18,23)
FLOOR = 5.0

def load_snapshots(root: Path):
    files = sorted(root.glob("source_*.json"))
    if len(files) != 2:
        raise ValueError("expected exactly two frozen source snapshots")
    subjects = {}
    for path in files:
        doc = json.loads(path.read_text())
        if str(doc["source_run_id"]) != RUN:
            raise ValueError(f"wrong source run in {path.name}")
        for row in doc["subjects"]:
            sid = row["subject"]
            if sid in subjects:
                raise ValueError(f"duplicate subject {sid}")
            if row["algorithm_commit"] != ALGORITHM or row["protocol_sha256"] != PROTOCOL:
                raise ValueError(f"mixed algorithm/protocol subject {sid}")
            if row["workflow_run"] != RUN or not row["completed"] or row["errors"]:
                raise ValueError(f"incomplete original source subject {sid}")
            subjects[sid] = row
    if tuple(sorted(subjects)) != EXPECTED:
        raise ValueError("not the complete 20-person original frozen evaluation")
    return subjects

def audit_arm(sid, arm, record, floor=FLOOR):
    d = record[arm]
    support = d["support"]
    result = {"subject":sid,"arm":arm,"readout_available":bool(d["available"]),
              "conditions":[],"support_floor":floor,"support_available":bool(support)}
    for s in support:
        reference = s["reference_trials"]
        kept = s["kept"]
        flagged = s["flagged_trials"]
        effective = s["effective_unflagged"]
        numeric = [reference, kept, flagged, effective]
        if any(not isinstance(x,(int,float)) or isinstance(x,bool) or not math.isfinite(x) for x in numeric):
            raise ValueError(f"invalid support metric for subject {sid} arm {arm}")
        if reference < 0 or kept < 0 or kept > reference or not 0 <= flagged <= reference or not 0 <= effective <= kept + 1e-6:
            raise ValueError(f"out-of-range support metric for subject {sid} arm {arm}")
        result["conditions"].append({
            "condition":s["condition"],"reference":reference,"kept":kept,
            "flagged_in_reference":flagged,"effective_unflagged":effective,
            "kept_pct":round(100*kept/reference,2) if reference else None,
            "flagged_reference_pct":round(100*flagged/reference,2) if reference else None,
            "effective_deficit_to_five":round(max(floor-effective,0),6),
            "meets_original_support_floor":effective>=floor})
    result["all_condition_support_meets_floor"] = bool(support) and all(
        x["meets_original_support_floor"] for x in result["conditions"])
    result["readout_and_support_eligible"] = result["readout_available"] and result["all_condition_support_meets_floor"]
    result["unassessable_without_trial_data"] = ["which exact urevents were blink-flagged",
        "whether eye flags are true-positive or conservative",
        "whether the candidate biases the genuine P300 waveform",
        "per-trial non-brain artifact severity after ICA"]
    return result

def make_audit(subjects):
    rows = []
    for sid in EXPECTED:
        rec = subjects[sid]
        for arm in ARMS:
            rows.append(audit_arm(sid,arm,rec))
    focused = [r for r in rows if r["subject"] in PRIMARY]
    return {"source":"original frozen 20-subject connector-captured result-marker snapshots; not raw trial-level EEG",
            "source_run_id":RUN,"algorithm_commit":ALGORITHM,"protocol_sha256":PROTOCOL,
            "preexisting_effective_trial_floor":FLOOR,
            "unique_subjects":len(EXPECTED),"arm_condition_rows":rows,
            "primary_subjects":list(PRIMARY),
            "primary_rows":focused,
            "subject_outcomes":{str(sid):{"primary_diagnostic":subjects[sid]["primary_diagnostic"],
                "veog_proxy":subjects[sid]["veog_blink_diagnostic"]} for sid in EXPECTED},
            "methodological_nonclaim":"Low support is not proof of false eye detections. Do not lower the frozen floor or post-hoc tune blinks to restore valid pairs."}

def brief(data):
    return [
        {"subject":r["subject"],"arm":r["arm"],"available":r["readout_available"],
         "effective":[q["effective_unflagged"] for q in r["conditions"]],
         "reference":[q["reference"] for q in r["conditions"]],
         "flagged":[q["flagged_in_reference"] for q in r["conditions"]],
         "eligible_support":r["all_condition_support_meets_floor"]}
        for r in data["primary_rows"]
    ]

def write_artifacts(data,out):
    out.mkdir(parents=True,exist_ok=True)
    (out/"effective_support.json").write_text(json.dumps(data,ensure_ascii=False,indent=2)+"\n")
    fields=["subject","arm","readout_available","condition","reference","kept",
            "flagged_in_reference","effective_unflagged","kept_pct",
            "flagged_reference_pct","effective_deficit_to_five","meets_original_support_floor"]
    with (out/"effective_support.csv").open("w",newline="") as file:
        writer=csv.DictWriter(file,fieldnames=fields);writer.writeheader()
        for row in data["arm_condition_rows"]:
            for c in row["conditions"]:
                writer.writerow({"subject":row["subject"],"arm":row["arm"],
                                 "readout_available":row["readout_available"],**{k:c[k] for k in fields if k in c}})
    md=[
        "# ERP CORE original frozen P300: effective support audit",
        "",
        "Source: verified, derived result-marker snapshots from original Actions run "+RUN+". No raw trial-level flags are available here.",
        "",
        "This is descriptive forensics; it is not an algorithm-comparison rerun and does not revise frozen inclusion gates.",
        "",
        "| Participant | Arm | Readout | Effective target / standard | Original reference target / standard | Flagged in reference target / standard |",
        "|---|---|---|---|---|---|"]
    for x in brief(data):
        eff=x["effective"];ref=x["reference"];flag=x["flagged"]
        fmt=lambda a:" / ".join(str(v) for v in a) if a else "not recorded"
        md.append(f"| {x['subject']} | {x['arm']} | {x['available']} | {fmt(eff)} | {fmt(ref)} | {fmt(flag)} |")
    md+=["","## Mandatory interpretation",
         "",
         "- The value effective_unflagged is copied from the original recorder, including fractional values produced by weighting; do NOT use retained-minus-flagged as a substitute.",
         "- For subjects 18/23, Standard had no recommendation; absence of its support array is NOT a measured zero.",
         "- For subjects 18/23, each fixed comparator had **zero retained target trials**. These data cannot yield a complete paired P300 measurement under the frozen rules.",
         "- Subject 4 Standard had 4 target / 3 standard effective unflagged trials; subject 8 fixed comparators had 4 target trials.",
         "- Subject 17 has a separate signal-shape problem; its fixed comparator target support was 4.",
         "- VEOG threshold is an external *proxy*, not ground truth. A high VEOG reading among flagged trials does not prove every individual flag correct.",
         "- These snapshots cannot identify the underlying eye-flagged urevent IDs or reconstruct rejected EEG. That requires a new versioned, per-trial read-only recorder.",
         "- Keep preexisting >=5 effective-unflagged support floor and all signal-shape gates unchanged. No artificial rescue of valid-pair count.",
         ""]
    (out/"effective_support.md").write_text("\n".join(md))

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--snapshot-dir",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    subjects=load_snapshots(args.snapshot_dir)
    data=make_audit(subjects)
    write_artifacts(data,args.out)
    print("PR51_SUPPORT_AUDIT",json.dumps({"source_run_id":RUN,"subjects":data["unique_subjects"],
        "original_support_floor":FLOOR,"priority":brief(data)},separators=(",",":")))
    return 0

if __name__=="__main__":
    raise SystemExit(main())
