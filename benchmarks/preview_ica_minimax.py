#!/usr/bin/env python3
"""Descriptive *proxy* minimax ranking over three pre-existing ICA starts.

This is a scientific development thought experiment, not a deployable
recommendation. Matrix columns contain seed-specific ranking objectives but
NOT per-seed candidate feasibility or signal-preservation gate status;
therefore numerical minimax improvement does not establish scientific
eligibility, stability, or improved P300 measurements. Do not alter PR51.
"""
import argparse
import json
import math
from pathlib import Path


def finite(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v)


def preview(subject):
    ids = subject["ids"]
    rows = subject["objective_by_candidate_and_seed"]
    rec = subject["seed_pipelines"]
    if not (len(ids) == len(rows) and len(rec) == 3 and len(set(ids)) == len(ids)):
        raise ValueError("bad alignment of ICA candidate IDs / objectives / recommendations")
    if not all(isinstance(row, list) and len(row) == 3 for row in rows):
        raise ValueError("must have one frozen and two alternate ICA-start objectives")
    # The original ranking objective is an uncertainty proxy, not the final P3 amplitude.
    candidates = [{"id": int(idx), "worst": max(row), "median": sorted(row)[1],
                   "by_seed": row}
                  for idx, row in zip(ids, rows) if all(finite(v) for v in row)]
    if not candidates:
        return {"subject": subject["subject"], "proposal": None,
                "reason": "no-all-seed-finite-objective-candidate",
                "scientific_feasibility_checked": False}
    best = min(x["worst"] for x in candidates)
    # Exact objective ties: no arbitrary tie-breaking is claimed from only
    # these matrices (proper tree-complexity precedence is not recorded here).
    tied = [x for x in candidates if x["worst"] == best]
    prior = next((x for x in candidates if x["id"] == rec[0]), None)
    result = {"subject": subject["subject"],
              "original_three_seed_selected_ids": rec,
              "usable_finite_rows": len(candidates),
              "minimax_proxy_worst": best,
              "tie_ids": [x["id"] for x in tied],
              "proposal": tied[0]["id"] if len(tied) == 1 else None,
              "original_first_seed_candidate": rec[0],
              "original_first_seed_candidate_worst": prior["worst"] if prior else None,
              "scientific_feasibility_checked": False,
              "full_waveform_preservation_checked": False,
              "bootstrap_se_checked": False,
              "requires_per_seed_feasible_status_and_signal_gates": True,
              "eligible_for_real_algorithm_promotion": False}
    if prior is not None:
        result["proxy_worst_case_reduction"] = prior["worst"]-best
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--snapshot",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    data=json.loads(a.snapshot.read_text())
    if data["scientific_gate_changes"] is not False:
        raise ValueError("snapshot must not alter original study scientific rules")
    results=[preview(s) for s in data["subjects"]]
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps({
        "purpose":"DEVELOPMENT_ONLY_NUMERIC_PROXY_NO_SCIENTIFIC_FEASIBILITY",
        "source_workflow":data["original_github_run"],
        "subjects":results},indent=2))
    for r in results:
        print("PR51_ICA_MINIMAX_PROXY",json.dumps(r))


if __name__=="__main__":
    main()
