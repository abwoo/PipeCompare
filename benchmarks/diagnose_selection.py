"""Post-hoc attribution of selection and error-proxy failures from saved CSVs.
No pipeline is fitted or deployed; truth is used only for oracle diagnostics.
"""
import argparse
import csv
import json
import math
from pathlib import Path
from statistics import mean

def number(row, key):
    return float(row[key])

def candidates(rows):
    """Calibration.rows writes each candidate's conditions contiguously.
    Labels need not be unique; a condition reset starts a new candidate.
    """
    blocks = []
    for row in rows:
        condition = number(row, "condition")
        if not blocks or row["recording"] != blocks[-1][0]["recording"] or condition <= number(blocks[-1][-1], "condition"):
            blocks.append([])
        blocks[-1].append(row)
    result = {}
    for block in blocks:
        first = block[0]
        assert all(r["steps"] == first["steps"] for r in block), "Mixed candidate labels"
        conditions = [number(r, "condition") for r in block]
        assert conditions == list(range(1, len(block) + 1)), "Missing/reordered candidate conditions"
        assert len({r["selected"].lower() for r in block}) == 1, "Mixed selection flags"
        g = number(first, "gain")
        assert g > 0 and all(abs(number(r, "gain")-g) < 1e-9 for r in block), "Inconsistent gain"
        known = all(math.isfinite(number(r, k)) for r in block for k in
                    ("actualError", "rawError", "reportedError", "sme", "eyeBias", "shapeSpread", "correctedEstimate"))
        if not known:
            continue
        result.setdefault(first["recording"], []).append(dict(
            label=first["steps"], selected=first["selected"].lower() in ("1", "true"),
            mse=mean(number(r, "actualError")**2 for r in block),
            raw_mse=mean(number(r, "rawError")**2 for r in block),
            reported_mse=mean(number(r, "reportedError")**2 for r in block),
            # Holding estimates and gain fixed; no causal interpretation.
            association_mse=mean((number(r, "sme")/g)**2 +
                (abs(number(r, "eyeBias"))/g +
                 abs(number(r, "correctedEstimate"))*number(r, "shapeSpread"))**2 for r in block),
            conditions=len(block)))
    for record, cs in result.items():
        assert len({c["conditions"] for c in cs}) == 1, "Unequal condition coverage: "+record
    return result

def aggregate(rows):
    if not rows:
        return dict(records=0)
    result = dict(records=len(rows))
    for field in ("selected", "score_minimum", "association_minimum", "oracle", "selected_raw", "raw_oracle"):
        result[field+"_rmse"] = math.sqrt(mean(r[field] for r in rows))
    result["mean_selected_minus_oracle_mse"] = mean(r["selected"]-r["oracle"] for r in rows)
    result["mean_score_minimum_minus_oracle_mse"] = mean(r["score_minimum"]-r["oracle"] for r in rows)
    actual = mean(r["selected"] for r in rows)
    for name in ("reported", "association"):
        reported = mean(r[name] for r in rows)
        result[name+"_ratio_at_unchanged_selection"] = math.sqrt(actual/reported) if reported > 0 else None
    result["scope"] = "Same assessable recordings; equal recording weight; oracle is an unattainable truth-informed lower bound within the saved candidate family."
    return result

def diagnose(folder):
    with (folder/"calibration.csv").open(newline="") as f:
        grouped = candidates(list(csv.DictReader(f)))
    meta = json.loads((folder/"summary.json").read_text())
    prevalence = {}
    if (folder/"recordings.csv").exists():
        with (folder/"recordings.csv").open(newline="") as f:
            prevalence = {r["subject"]:float(r["blinkPrevalence"]) for r in csv.DictReader(f)}
    rows = []
    unavailable = []
    for record, cs in grouped.items():
        selected = [c for c in cs if c["selected"]]
        if len(selected) != 1:
            unavailable.append(dict(recording=record, reason="No unique assessable recommendation"))
            continue
        pick = selected[0]
        oracle = min(cs, key=lambda c:c["mse"])
        raw_oracle = min(cs, key=lambda c:c["raw_mse"])
        score = min(cs, key=lambda c:c["reported_mse"])
        assoc = min(cs, key=lambda c:c["association_mse"])
        rows.append(dict(recording=record, prevalence=prevalence.get(record),
            candidates=len(cs), selected_label=pick["label"], oracle_label=oracle["label"],
            selected=pick["mse"], selected_raw=pick["raw_mse"],
            reported=pick["reported_mse"], association=pick["association_mse"],
            oracle=oracle["mse"], raw_oracle=raw_oracle["raw_mse"],
            score_minimum=score["mse"], association_minimum=assoc["mse"]))
    planned = meta.get("plannedRecords", meta.get("planned"))
    return dict(source_commit=meta["commit"], planned=planned,
        assessable_selected=len(rows), unavailable_planned=planned-len(rows) if planned else None,
        aggregate=aggregate(rows), prevalence=[dict(prevalence=p, **aggregate([r for r in rows if r["prevalence"]==p]))
            for p in sorted({r["prevalence"] for r in rows if r["prevalence"] is not None})],
        records=rows, excluded_saved_records=unavailable)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    folders = sorted(args.root.glob("*validation-analysis"))
    assert len(folders) == 2, "Expected exactly the public and synthetic saved analyses"
    result = dict(scope="Post-hoc development diagnosis, not new independent confirmation or a deployment recommendation.",
        source_run=39974230047, analyses={p.name:diagnose(p) for p in folders})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, ensure_ascii=False, allow_nan=False, indent=2)+"\n")
    print("PIPECOMPARE_SELECTION_DIAGNOSIS "+json.dumps(result, allow_nan=False))

if __name__ == "__main__":
    main()
