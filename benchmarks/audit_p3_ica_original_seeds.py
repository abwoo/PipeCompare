#!/usr/bin/env python3
"""Frozen original ERP CORE P3 ICA seed objective-matrix audit.

Read-only: the 0,1,2 columns are original search plus two ICA starts.
The original marker does NOT contain alternate-seed P300 estimates or SMEs.
Never interpret objective similarity as measurement-level robustness.
"""
import argparse
import json
import math
from pathlib import Path

EXPECTED={11:([12,25,25],20),15:([5,3,3],19),21:([25,4,4],19)}
RUN="42809789847"
ALG="7151821edf8a5a54c3167d7f18b66c530911b568"
PROTOCOL="bb65c9a522e9f45210e5ee941a669630b75587a5255664aa0f1402b572e89c6b"

def audit(original: dict) -> dict:
    subject=original["subject"]
    if subject not in EXPECTED or original["source_run_id"] != RUN or (
            original["algorithm_commit"] != ALG or original["protocol_sha256"] != PROTOCOL):
        raise ValueError("subject / frozen source identity mismatch")
    data=original["icaSeedCheck"]
    ids=list(data["ids"])
    objectives=data["objective"]
    recommended=data["recommended"]
    if (data["seeds"] != [1,2] or len(ids) != len(set(ids)) or
            len(ids) != EXPECTED[subject][1] or recommended != EXPECTED[subject][0] or
            len(objectives) != len(ids) or any(len(row) != 3 for row in objectives) or
            not all(v in ids for v in recommended)):
        raise ValueError("ICA seed objective matrix dimension or recommendation mismatch")
    chosen=[]
    for seed in range(3):
        objective=objectives[ids.index(recommended[seed])][seed]
        if not isinstance(objective,(int,float)) or not math.isfinite(objective):
            raise ValueError("non-finite selected objective")
        chosen.append({"seed":seed,"recommended_id":recommended[seed],"objective":objective})
    candidate_rows=[]
    for id,row in zip(ids,objectives):
        finite=[x for x in row if isinstance(x,(int,float)) and math.isfinite(x)]
        if not finite:continue
        first=row[0]
        variation=(max(finite)-min(finite))/abs(first)*100 if math.isfinite(first) and first!=0 else None
        candidate_rows.append({"id":id,"objective_by_seed":row,"max_relative_change_pct":variation})
    largest=sorted((r for r in candidate_rows if r["max_relative_change_pct"] is not None),
                   key=lambda r:-r["max_relative_change_pct"])[:3]
    vals=[x["objective"] for x in chosen]
    return {"subject":subject,"source_run_id":RUN,"algorithm_commit":ALG,
            "protocol_sha256":PROTOCOL,"n_candidates_compared":len(ids),
            "recommended":recommended,"same_recommendation":len(set(recommended))==1,
            "selected_objective_by_seed":chosen,
            "selected_objective_spread_pct":(max(vals)-min(vals))/min(vals)*100,
            "highest_candidate_relative_variation":largest,
            "no_alternate_seed_p3_estimate_in_original_source":True,
            "no_inference_on_P3_amplitude_equivalence":True}

def run(directory: Path,out: Path):
    originals=[json.loads((directory/f"subject_{subject:03}.json").read_text())
               for subject in sorted(EXPECTED)]
    reports=[audit(d) for d in originals]
    out.mkdir(parents=True,exist_ok=True)
    (out/"seed_objective_audit.json").write_text(json.dumps(reports,indent=2)+"\n")
    lines=["# Frozen original P3 ICA seed-objective audit","",
           "The original seed check stores selection objective values and candidate IDs, **not** alternate-seed P300 amplitudes, SEs, removed ICs or effective support.",
           "",
           "| Subject | Candidate IDs (original / seed 1 / seed 2) | Selected objectives | Spread (%) |",
           "|---|---|---|---:|"]
    for r in reports:
        ids="/".join(str(x) for x in r["recommended"])
        vals="/".join(f"{m['objective']:.5f}" for m in r["selected_objective_by_seed"])
        lines.append(f"| {r['subject']} | {ids} | {vals} | {r['selected_objective_spread_pct']:.2f} |")
    lines+=["","## Interpretation",
            "Different recommended IDs do **not** in themselves prove a significant P300 measurement difference.",
            "Conversely, similar scalar objectives do **not** prove P300 amplitude or SME equivalence; those quantities were not recorded per alternate seed.",
            "The next isolated, versioned real P3 reruns must measure seed-specific amplitude, uncertainty, IC removal and source-aligned HMM-clean support. Do not alter the original >=18/20 stability or >=5 effective support gates.",
            ""]
    (out/"seed_objective_audit.md").write_text("\n".join(lines))
    print("PR51_ORIGINAL_ICA_SEEDS "+json.dumps(
        [{"subject":r["subject"],"ids":r["recommended"],"relative_spread_pct":r["selected_objective_spread_pct"]}
         for r in reports],separators=(",",":")))
    return reports

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-dir",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    run(a.original_dir,a.out)

if __name__=="__main__":
    main()
