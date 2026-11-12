#!/usr/bin/env python3
"""Original real-ERP CORE P3 ICA selection identity / ranking-objective audit.

This is a READ-ONLY analysis of the saved original experiment markers.
Ranking objective is not the P300 amplitude nor the per-subject bootstrap SE.
No ICA reruns, new seeds, new acceptance criteria, or PR51 status changes.
"""
import argparse
import json
import math
from pathlib import Path

PLANNED = [1, 2, 3, 4, 5, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]
INCOMPLETE_PAIRS = {4, 17, 18, 23}  # original collector; support exclusion of 8 is subsequent

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    data = []
    for path in sorted(args.source_dir.glob("source_*.json")):
        data += json.loads(path.read_text())["subjects"]
    data.sort(key=lambda x:x["subject"])
    assert [x["subject"] for x in data] == PLANNED

    rows = []
    for s in data:
        ica = s["ica_seeds"]
        if not ica["available"]:
            rows.append({"subject":s["subject"],"fit_available":False,
                         "paired_control_available":False,
                         "identity_same":None,"picks":[],"ranking_objectives":[]})
            continue
        ids = ica["recommended_candidate_ids"]
        objs = ica["recommended_objective"]
        if len(ids) != 3 or len(objs) != 3:
            raise ValueError(f"Unexpected seeds for subject {s['subject']}: {ids}")
        if any(y is None or not math.isfinite(y) for y in objs):
            raise ValueError(f"Unassessable saved ranking objective for subject {s['subject']}")
        agrees = len(set(ids)) == 1
        if agrees != bool(ica["same"]):
            raise ValueError(f"Saved ICA agreement disagrees with ids for subject {s['subject']}")
        rows.append({
            "subject":s["subject"],"fit_available":True,
            "paired_control_available":s["subject"] not in INCOMPLETE_PAIRS,
            "identity_same":agrees,"picks":ids,
            "ranking_objectives":objs,
            "objective_range":max(objs)-min(objs),
            "objective_range_pct_of_original":100*(max(objs)-min(objs))/abs(objs[0]) if objs[0] else None,
            "note":"ranking objective only: cannot infer P3 amplitude/SE equivalence"
        })
    allfit=[r for r in rows if r["fit_available"]]
    paired=[r for r in allfit if r["paired_control_available"]]
    stable_all=sum(r["identity_same"] for r in allfit)
    stable_paired=sum(r["identity_same"] for r in paired)
    assert len(allfit)==17 and stable_all==14
    assert len(paired)==16 and stable_paired==13
    unstable=[r["subject"] for r in paired if not r["identity_same"]]
    assert unstable==[11,15,21]
    args.out.mkdir(parents=True,exist_ok=True)
    payload={
        "source":"original ERP CORE P3 frozen MATLAB result markers",
        "planned_n":20,"ICA_fits_available_n":17,"identity_same_among_fits_n":14,
        "paired_comparator_n":16,"identity_same_among_paired_n":13,
        "frozen_requirement":"at least 18 same recommendations out of planned 20",
        "unstable_subjects":unstable,
        "methodological_caveat":"recommendation identity and ranking objective are not P300 amplitude or bootstrap-SE stability",
        "records":rows
    }
    (args.out/"ica_seed_audit.json").write_text(json.dumps(payload,indent=2)+"\n")
    md=["# Frozen P3 ICA seed-choice diagnostics", "",
        "Original 20 planned, 17 with ICA selections, 14/17 same identities.",
        "Among 16 subjects with complete numerical comparator arms: **13/16** same (original collector).",
        "Prespecified full-cohort gate remains **18/20**; not met.",
        "**Different candidate IDs do not yet establish different P300 amplitude or SE.**",
        "**Same candidate ID does not establish numerical P300 stability either.**",
        "", "| Subject | Baseline / seed 1 / seed 2 IDs | Ranking objective range relative to baseline |",
        "|---:|---|---:|"]
    for r in allfit:
        md.append(f'| {r["subject"]} | {"/".join(str(i) for i in r["picks"])} | {r["objective_range_pct_of_original"]:.2f}% |')
    md += ["", "## Next confirmatory work", "",
           "Re-run prespecified ICA seeds with identical trial and channel definitions; save paired P3",
           "difference amplitude, bootstrap SE, identified IC topographies, and removed-component maps.",
           "Use independent component matching and a predetermined amplitude-equivalence margin;",
           "the saved ranking objective is insufficient for claiming P3 measurement reproducibility.",
           ""]
    (args.out/"ica_seed_audit.md").write_text("\n".join(md))
    print("PR51_FROZEN_ICA_CHOICES",json.dumps({k:payload[k] for k in ("ICA_fits_available_n","identity_same_among_fits_n","paired_comparator_n","identity_same_among_paired_n","unstable_subjects")}))
if __name__=="__main__":
    main()
