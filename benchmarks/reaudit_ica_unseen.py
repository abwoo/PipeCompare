#!/usr/bin/env python3
"""Post-run scientific audit on archived six-seed ERP CORE development evidence.

No EEG re-fitting, no seed re-selection, and no new participants. The raw
GitHub run is immutable; this checks the original artifacts against the
pre-outcome independent validation amendment A1.
"""
import argparse
import hashlib
import json
from pathlib import Path
from validate_ica_unseen_seeds import inspect

EXPECTED = {11, 15, 21}
RAW_RUN = "43169562528"
RAW_SHA = "5c9388562974c8305da95cb59e7249982788d589"


def reaudit(evidence: Path):
    files = sorted(evidence.rglob("result.json"))
    rows, seen, failures = [], set(), []
    if len(files) != len(EXPECTED):
        failures.append(f"expected 3 archived EEG result.json files, found {len(files)}")
    for file in files:
        result = json.loads(file.read_text())
        subject = result.get("subject")
        if subject not in EXPECTED or subject in seen:
            failures.append(f"unexpected/duplicate subject: {subject}")
            continue
        seen.add(subject)
        meta_file = file.with_name("exploratory_identity.json")
        if not meta_file.is_file():
            failures.append(f"missing source identity for subject {subject}")
            continue
        meta = json.loads(meta_file.read_text())
        if str(meta.get("workflow_run")) != RAW_RUN or meta.get("implementation_commit") != RAW_SHA:
            failures.append(f"source run/commit mismatch for subject {subject}")
            continue
        try:
            outcome = inspect(result)
        except (ValueError, TypeError, KeyError, IndexError) as exc:
            failures.append(f"subject {subject}: {type(exc).__name__}: {exc}")
            continue
        old_file = file.with_name("ica_unseen.json")
        if not old_file.is_file():
            failures.append(f"missing original scientific interpretation for subject {subject}")
            continue
        prior = json.loads(old_file.read_text())
        if prior.get("subject") != subject:
            failures.append(f"original analysis subject mismatch for {subject}")
            continue
        for which in ("originalBaseline", "shadowUnseen"):
            if prior[which]["candidateId"] != outcome[which]["candidateId"]:
                failures.append(f"candidate ID drift in {which} for subject {subject}")
        baseline = outcome["originalBaseline"]
        shadow = outcome["shadowUnseen"]
        old_range = baseline.get("unseenP3Range")
        new_range = shadow.get("unseenP3Range")
        rows.append({
            "subject": subject,
            "rawResultSha256": hashlib.sha256(file.read_bytes()).hexdigest(),
            "originalCandidateId": baseline["candidateId"],
            "shadowCandidateId": shadow["candidateId"],
            "originalUnseenValid": baseline["validAllThree"],
            "shadowUnseenValid": shadow["validAllThree"],
            "oldP3Range": old_range,
            "shadowP3Range": new_range,
            "rangeChange": new_range - old_range
                if old_range is not None and new_range is not None else None,
            "oldValidBeforeCorrection": prior["originalBaseline"]["validAllThree"],
            "shadowValidBeforeCorrection": prior["shadowUnseen"]["validAllThree"],
            "trainingSelectedOnlyFromSeeds012": outcome["selectionSeeds"] == [0, 1, 2],
            "oldUnseenSeedChecks": baseline.get("seeds", []),
            "shadowUnseenSeedChecks": shadow.get("seeds", []),
        })
    missing = EXPECTED - seen
    if missing:
        failures.append(f"missing subjects {sorted(missing)}")
    rows.sort(key=lambda x: x["subject"])
    return {
        "analysis": "ICA independent validation A1, reuse archived raw EEG",
        "rawWorkflowRun": RAW_RUN, "rawWorkflowHead": RAW_SHA,
        "newMatlabComputation": False,
        "newParticipantsEvaluated": False,
        "pr51Merged": False,
        "original20SubjectVerdict": "PARTIAL", "original20SubjectValidPairs": 15,
        "completeSourceEvidence": len(rows) == 3 and not failures,
        "failures": failures, "rows": rows,
        "allShadowUnseenFeasible": (
            all(row["shadowUnseenValid"] for row in rows)
            if len(rows) == 3 and not failures else None),
        "anyShadowP3RangeImprovement": any(
            row["rangeChange"] is not None and row["rangeChange"] < 0
            for row in rows),
        "productionPromotionAuthorized": False,
        "freeze14SubjectValidationAuthorized": False,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    answer = reaudit(args.evidence)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(answer, indent=2, allow_nan=True))
    for row in answer["rows"]:
        print("PR51_INDEPENDENT_ICA_A1", json.dumps({
            k: row[k] for k in (
                "subject", "originalCandidateId", "shadowCandidateId",
                "originalUnseenValid", "shadowUnseenValid",
                "oldP3Range", "shadowP3Range", "rangeChange",
                "oldValidBeforeCorrection", "shadowValidBeforeCorrection")}))
    print("PR51_INDEPENDENT_ICA_A1_SUMMARY", json.dumps({
        k: answer[k] for k in (
            "completeSourceEvidence", "failures", "allShadowUnseenFeasible",
            "anyShadowP3RangeImprovement", "productionPromotionAuthorized")}))
    if not answer["completeSourceEvidence"]:
        raise SystemExit("Archived scientific evidence incomplete; fail closed")


if __name__ == "__main__":
    main()
