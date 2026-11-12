#!/usr/bin/env python3
"""Descriptive frozen-cohort ERP CORE P3 trial-support audit.

Operates only on preserved original-job result markers. No EEG fitting,
trial relabeling, threshold relaxation, subject omission or new scientific
decisions. The 5-trial effective nonblink floor is copied from the original
pre-registered ERP CORE protocol.
"""
import argparse
from collections import Counter
import json
import math
from pathlib import Path

EXPECTED = [1, 2, 3, 4, 5, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]
FLOOR = 5
RUN = 42809789847
ALGORITHM = "7151821edf8a5a54c3167d7f18b66c530911b568"
PROTOCOL = "bb65c9a522e9f45210e5ee941a669630b75587a5255664aa0f1402b572e89c6b"
ARMS = ("standard", "minimal", "common")


def read_original_snapshots(folder):
    paths = sorted(Path(folder).glob("source_*.json"))
    if len(paths) != 2:
        raise ValueError(f"expected two original evidence snapshots; found {len(paths)}")
    records = {}
    identities = set()
    for path in paths:
        root = json.loads(path.read_text())
        if root["source_run_id"] != RUN:
            raise ValueError(f"not original run in {path}")
        for s in root["subjects"]:
            idx = s["subject"]
            if idx in records:
                raise ValueError(f"duplicate original subject {idx}")
            identity = (s["workflow_run"], s["algorithm_commit"], s["execution_commit"], s["protocol_sha256"])
            identities.add(identity)
            if s["workflow_run"] != str(RUN) or s["algorithm_commit"] != ALGORITHM or s["protocol_sha256"] != PROTOCOL:
                raise ValueError(f"wrong frozen identity at subject {idx}")
            if not s["completed"] or s["errors"]:
                raise ValueError(f"original incomplete source {idx}")
            records[idx] = s
    if sorted(records) != EXPECTED or len(identities) != 1:
        raise ValueError("cohort coverage or original execution consistency failed")
    return [records[k] for k in EXPECTED]


def analyze_arm(data):
    rows = data.get("support") or []
    if not rows:
        return {"available": bool(data.get("available")), "support_recorded": False, "conditions": []}
    if len(rows) != 2:
        raise ValueError("must have two condition-specific support rows")
    result = []
    for r in rows:
        n = r["reference_trials"]
        flags = r["flagged_trials"]
        eff = r["effective_unflagged"]
        kept = r["kept"]
        if n is None or flags is None or eff is None or n < 0 or flags < 0 or flags > n or kept < 0:
            raise ValueError(f"invalid trial-support evidence: {r}")
        if not isinstance(eff, (int, float)) or not math.isfinite(eff) or eff < 0:
            raise ValueError(f"nonfinite effective nonblink support: {r}")
        result.append({
            "condition": r["condition"],
            "reference_trials": n,
            "retained_trials": kept,
            "flagged_reference_trials": flags,
            "observed_unflagged": n - flags,
            "flagged_reference_fraction": flags / n if n else None,
            "effective_unflagged": eff,
            "meets_frozen_floor": eff >= FLOOR
        })
    return {"available": bool(data.get("available")), "support_recorded": True,
            "conditions": result,
            "support_valid": all(r["meets_frozen_floor"] for r in result)}


def inspect(records):
    entries = []
    for s in records:
        arms = {arm: analyze_arm(s[arm]) for arm in ARMS}
        numerically_available = all(arms[a]["available"] for a in ARMS)
        valid_support = all(arms[a].get("support_valid", False) for a in ARMS)
        notes = []
        for arm in ARMS:
            for cond in arms[arm]["conditions"]:
                if cond["effective_unflagged"] < FLOOR:
                    notes.append(f'{arm}/{cond["condition"]}: effective={cond["effective_unflagged"]:.2f} < {FLOOR}')
                if cond["flagged_reference_fraction"] is not None and cond["flagged_reference_fraction"] >= 0.90:
                    notes.append(f'{arm}/{cond["condition"]}: >=90% reference epochs flagged (not optical ground truth)')
        entries.append({
            "subject": s["subject"], "original_primary_diagnostic": s["primary_diagnostic"],
            "numerically_available_arms": numerically_available,
            "supports_complete_comparison": numerically_available and valid_support,
            "arms": arms,
            "observed_veog_diagnostic": s.get("veog_blink_diagnostic"),
            "notes": notes
        })
    return entries


def markdown(data):
    text = ["# Original ERP CORE P3 per-condition effective nonblink support",
            "", "Evidence: frozen original MATLAB result markers from Actions run 42809789847.",
            "This audit does not reprocess EEG. A blink flag is not independent optical ground truth.",
            "Frozen minimum effective unflagged support per condition/arm: 5; unchanged.",
            "", "| Subject | Standard target/standard | Minimal target/standard | Common target/standard | Numeric arms | Support floor |",
            "|---:|---:|---:|---:|---|---|"]
    for s in data:
        def fmt(arm):
            rows = s["arms"][arm]["conditions"]
            return "/".join(f'{z["effective_unflagged"]:.2f}' for z in rows) if rows else "unavailable"
        text.append(f'| {s["subject"]} | {fmt("standard")} | {fmt("minimal")} | {fmt("common")} | '
                    f'{"all" if s["numerically_available_arms"] else "incomplete"} | '
                    f'{"yes" if s["supports_complete_comparison"] else "no"} |')
    text += ["", "## Cases requiring trial-level investigation", ""]
    for s in data:
        if s["subject"] not in (4, 8, 17, 18, 23):
            continue
        text.append(f'### Subject {s["subject"]}')
        for note in s["notes"]:
            text.append("- " + note)
        if not s["notes"]:
            text.append("- Source not assessable at Standard arm; inspect filter/evaluation reason")
        text.append("")
    text += [
        "## Required next measurements (NOT calculated by this snapshot)",
        "",
        "- Per-trial rejection reason with epoch IDs and temporal overlap",
        "- VEOG amplitude, baseline drift and EOG/EEG saturation under pinned preprocessing",
        "- ICA removal and regression decisions before/after filtering",
        "- Match suspicious HMM flags to genuinely independent ocular evidence, not a fixed VEOG heuristic",
        "- Recompute effective support after any *prospectively specified* fix; retain every negative record",
        "- Do not silently switch to a longer excerpt or relax the original 5-trial floor",
    ]
    return "\n".join(text) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-dir", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    opts = ap.parse_args()
    records = read_original_snapshots(opts.source_dir)
    entries = inspect(records)
    by_subject = {s["subject"]: s for s in entries}
    assert by_subject[4]["arms"]["standard"]["conditions"][0]["effective_unflagged"] == 4
    assert by_subject[4]["arms"]["standard"]["conditions"][1]["effective_unflagged"] == 3
    assert by_subject[8]["arms"]["minimal"]["conditions"][0]["effective_unflagged"] == 4
    for i in (18, 23):
        assert by_subject[i]["arms"]["minimal"]["conditions"][0]["effective_unflagged"] == 0
    valid = [s["subject"] for s in entries if s["supports_complete_comparison"]]
    assert len(valid) == 15, f"frozen cohort support accounting mismatch, got {len(valid)}"
    opts.out.mkdir(parents=True, exist_ok=True)
    payload = {"source": "original frozen MATLAB result markers; NOT reprocessed EEG",
               "run_id": RUN, "prespecified_floor": FLOOR,
               "valid_subjects": valid, "n_valid": len(valid), "records": entries}
    (opts.out / "support_audit.json").write_text(json.dumps(payload, indent=2) + "\n")
    (opts.out / "support_audit.md").write_text(markdown(entries))
    print("PR51_ORIGINAL_SUPPORT_AUDIT",
          json.dumps({"planned":len(entries),"valid":len(valid),
                      "excluded_or_unavailable":[s["subject"] for s in entries
                                                if not s["supports_complete_comparison"]]}))


if __name__ == "__main__":
    main()
