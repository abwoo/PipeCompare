#!/usr/bin/env python3
"""Read-only forensic audit of frozen ERP CORE GitHub Actions result logs.

No EEG fitting, production edits, threshold changes, or scientific reclassification.
Uses the exact PIPECOMPARE_ERPCORE_RESULT records emitted by the frozen MATLAB run.
All 20 planned subjects must be present and all identities must agree.
"""
import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import re
import urllib.request
import zipfile
from io import BytesIO

REPO = "abwoo/PipeCompare"
MARKER = "PIPECOMPARE_ERPCORE_RESULT "
RUN_ID = 42809789847


def get_bytes(url, token):
    req = urllib.request.Request(url, headers={
        "Authorization": f"Bearer {token}",
        "Accept": "application/vnd.github+json",
        "X-GitHub-Api-Version": "2022-11-28",
        "User-Agent": "PipeCompare-forensics"})
    with urllib.request.urlopen(req, timeout=60) as response:
        return response.read()


def fetch_log(job_id, token):
    data = get_bytes(f"https://api.github.com/repos/{REPO}/actions/jobs/{job_id}/logs", token)
    if data.startswith(b"PK\x03\x04"):
        with zipfile.ZipFile(BytesIO(data)) as archive:
            names = sorted(x for x in archive.namelist() if x.endswith((".txt", ".log")))
            if not names:
                raise ValueError(f"no text logs in job {job_id}")
            data = b"\n".join(archive.read(name) for name in names)
    return data.decode("utf-8", "replace")


def fetch_jobs(run_id, token):
    url = f"https://api.github.com/repos/{REPO}/actions/runs/{run_id}/jobs?per_page=100"
    jobs = json.loads(get_bytes(url, token))["jobs"]
    if len(jobs) >= 100:
        raise ValueError("job list may be paginated; do not silently omit subjects")
    return jobs


def result_from_log(text):
    found = []
    for line in text.splitlines():
        if MARKER in line:
            source = line.split(MARKER, 1)[1].strip()
            found.append(json.loads(source))
    if len(found) != 1:
        raise ValueError(f"expected exactly one frozen MATLAB result; got {len(found)}")
    return found[0]


def failure_counts(record):
    reasons = record.get("noRecommendationReasons") or []
    if record.get("erp", {}).get("available"):
        return {}
    return dict(sorted(Counter(str(x) for x in reasons).items(),
                       key=lambda p: (-p[1], p[0])))


def support(record, arm):
    rows = record.get("selectedSupport") if arm == "erp" else record.get(arm, {}).get("support")
    out = []
    for row in rows or []:
        out.append({
            "condition": str(row.get("condition")),
            "kept": int(row.get("kept", 0)),
            "reference_trials": len(row.get("referenceIds") or []),
            "flagged_trials": sum(bool(x) for x in (row.get("flagged") or [])),
            "effective_unflagged": row.get("effectiveUnflaggedN")
        })
    return out


def seed_summary(record):
    obj = record.get("icaSeedCheck")
    if not isinstance(obj, dict):
        return {"available": False}
    picks = obj.get("recommended") or []
    ids = obj.get("ids") or []
    matrix = obj.get("objective") or []
    vals = []
    for seed_index, rec_id in enumerate(picks):
        if rec_id in ids:
            row = matrix[ids.index(rec_id)]
            vals.append(row[seed_index] if seed_index < len(row) else None)
        else:
            vals.append(None)
    return {"available": bool(picks), "recommended_candidate_ids": picks,
            "recommended_objective": vals, "same": obj.get("same")}


def classify(record):
    reasons = failure_counts(record)
    if not record.get("erp", {}).get("available"):
        if any("Unable to perform assignment" in key for key in reasons):
            return "icremove-execution-error-and-other-gates"
        if any("peak latency shifted" in key for key in reasons):
            return "signal-shape-gate"
        return "other-unavailable-standard"
    if not record.get("minimal", {}).get("available") or not record.get("common", {}).get("available"):
        return "fixed-comparator-unavailable"
    for arm in ("erp", "minimal", "common"):
        if any(row["effective_unflagged"] is not None
               and row["effective_unflagged"] < 5 for row in support(record, arm)):
            return "effective-unflagged-support-below-frozen-floor"
    if not seed_summary(record)["same"]:
        return "ica-seed-recommendation-instability"
    return "available"


def summarize(record, log_digest):
    return {
        "subject": record["subject"], "completed": record.get("completed"),
        "algorithm_commit": record.get("algorithmCommit"),
        "execution_commit": record.get("executionCommit"),
        "workflow_run": record.get("workflowRun"),
        "protocol_sha256": record.get("protocolSha256"),
        "log_sha256": log_digest,
        "candidate_count": record.get("candidates"), "feasible_count": record.get("feasible"),
        "primary_diagnostic": classify(record),
        "reasons": failure_counts(record),
        "standard": {
            "available": record.get("erp", {}).get("available"),
            "kept": record.get("erp", {}).get("kept"),
            "reported_difference_se": record.get("erp", {}).get("reportedDifferenceSE"),
            "support": support(record, "erp")},
        "minimal": {
            "available": record.get("minimal", {}).get("available"),
            "support": support(record, "minimal")},
        "common": {
            "available": record.get("common", {}).get("available"),
            "support": support(record, "common")},
        "ica_seeds": seed_summary(record),
        "veog_blink_diagnostic": record.get("hmmVeog"),
        "errors": record.get("errors") or []
    }


def build_report(rows, run_id):
    lines = [
        "# ERP CORE / PR51 frozen-run root cause audit",
        "",
        f"Original GitHub Actions run: https://github.com/{REPO}/actions/runs/{run_id}",
        "",
        "Read-only diagnostic from original, frozen per-subject MATLAB result logs.",
        "Does not rerun an experiment, amend a protocol, or establish improvement.",
        "",
        "| Subject | Standard feasible | Primary finding | ICA seed picks |",
        "|---:|---:|---|---|",
    ]
    for r in rows:
        seed = r["ica_seeds"]
        picks = "/".join(str(x) for x in seed.get("recommended_candidate_ids", [])) or "-"
        lines.append(f'| {r["subject"]} | {r["feasible_count"]}/{r["candidate_count"]} '
                     f'| {r["primary_diagnostic"]} | {picks} |')
    lines += ["", "## Exact zero-feasible reason groups", ""]
    for r in rows:
        if r["feasible_count"] != 0:
            continue
        lines.append(f'### Subject {r["subject"]}')
        for reason, n in r["reasons"].items():
            lines.append(f"- {n} candidates: {reason}")
        lines.append("")
    lines += ["## Prespecified constraints and guardrails", "",
              "- Original cohort and all planned records retained; no subject replacement.",
              "- The original 20/20 complete-pairs requirement remains **unchanged**.",
              "- The 18/20 ICA-stability requirement remains **unchanged**.",
              "- Missing ICA-seed records are not counted as a stable or an unstable ICA fit.",
              "- Compare effective unflagged support, not just retained trial counts.",
              "- Never reinterpret the existing 15-pair pilot as a complete 20-pair confirmation.",
              "- Candidate exceptions must be distinguished from physiological exclusions.",
              "- Compare amplitude stability and pipeline identity stability separately.",
              ""]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-id", type=int, default=RUN_ID)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--offline-dir", type=Path)
    args = parser.parse_args()
    if args.offline_dir is None:
        token = os.environ.get("GH_TOKEN")
        if not token:
            parser.error("GH_TOKEN required when fetching GitHub Actions logs")
        jobs = fetch_jobs(args.run_id, token)
        pairs = []
        for j in jobs:
            match = re.search(r"ERP CORE P3 participant (\d+)", j["name"])
            if match:
                pairs.append((int(match.group(1)), fetch_log(j["id"], token)))
    else:
        pairs = [(int(p.stem), p.read_text()) for p in sorted(args.offline_dir.glob("*.log"))]
    by_subject = {}
    for subject, content in pairs:
        record = result_from_log(content)
        if record["subject"] != subject:
            raise ValueError(f"subject mismatch in job {subject}")
        if subject in by_subject:
            raise ValueError(f"duplicate subject job {subject}; no silent replacement")
        by_subject[subject] = summarize(record, hashlib.sha256(content.encode()).hexdigest())
    expected = [1, 2, 3, 4, 5, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]
    if sorted(by_subject) != expected:
        raise ValueError(f"coverage mismatch; expected {expected}, found {sorted(by_subject)}")
    rows = [by_subject[s] for s in expected]
    for key in ("algorithm_commit", "execution_commit", "workflow_run", "protocol_sha256"):
        if len({r[key] for r in rows}) != 1:
            raise ValueError(f"mixed frozen-run identity: {key}")
    if any(not r["completed"] or r["errors"] for r in rows):
        raise ValueError("at least one frozen run failed; do not treat as complete")
    args.out.mkdir(parents=True, exist_ok=True)
    output = {"run_id": args.run_id, "source": "original GitHub Actions job logs",
              "subjects": rows}
    (args.out / "forensics.json").write_text(json.dumps(output, indent=2, allow_nan=True) + "\n")
    (args.out / "forensics.md").write_text(build_report(rows, args.run_id) + "\n")
    print("PR51_FORENSICS", json.dumps({
        "run_id": args.run_id, "subjects": len(rows),
        "findings": dict(Counter(r["primary_diagnostic"] for r in rows)),
        "unavailable_standard": [r["subject"] for r in rows if not r["standard"]["available"]],
        "unstable_ica": [r["subject"] for r in rows
                         if r["ica_seeds"]["available"] and not r["ica_seeds"]["same"]]
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
