#!/usr/bin/env python3
"""Explain real ICA recommendation stability from original six-start artifacts.

Read-only: no EEG fits, no threshold edits, no outcome-based selection.
The MAT archive describes IC removal ONLY in seed-zero fitting; it does not
contain all seed-specific IC maps or waveforms. Missing fields never
masquerade as an unchanged ICA decomposition.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.io import loadmat
from validate_ica_unseen_seeds import grid, inspect, valid

SUBJECTS = {11, 15, 21}
SOURCE_RUN = "43169562528"
SOURCE_SHA = "5c9388562974c8305da95cb59e7249982788d589"


def elements(obj):
    if isinstance(obj, np.ndarray):
        return list(obj.flat)
    if isinstance(obj, (list, tuple)):
        return list(obj)
    return [obj] if obj is not None else []


def scalar(obj):
    return obj.item() if isinstance(obj, np.ndarray) and obj.size == 1 else obj


def removal_snapshot(candidate):
    """No component-index alignment is claimed beyond one fixed seed."""
    steps = elements(getattr(candidate, "steps", None))
    kinds = [str(scalar(getattr(s, "type", ""))) for s in steps]
    ica = "ica" in kinds
    removed = []
    recorded = False
    for s in steps:
        if str(scalar(getattr(s, "type", ""))) == "icremove":
            n = scalar(getattr(s, "icsRemoved", float("nan")))
            if isinstance(n, (float, int, np.number)) and math.isfinite(float(n)):
                removed.append(int(n))
                recorded = True
    raw = scalar(getattr(candidate, "icsRemoved", float("nan")))
    if not removed and isinstance(raw, (float, int, np.number)) and math.isfinite(float(raw)):
        removed = [int(raw)]
        recorded = True
    if not ica:
        mechanism = "no-ica-step"
    elif not recorded:
        mechanism = "ica-removal-unreported"
    elif sum(removed) == 0:
        mechanism = "ica-fitted-no-components-removed-at-seed0"
    else:
        mechanism = "ica-fitted-components-removed-at-seed0"
    return {"seed0IcaPresent": ica,
            "seed0RemovedComponents": sum(removed) if recorded else None,
            "seed0RemovalRecorded": recorded, "seed0Mechanism": mechanism,
            "seed0IcRemovePresent": "icremove" in kinds}


def condition_summary(m):
    if not isinstance(m, dict):
        return {"available": False}
    eff = m.get("effectiveUnflagged")
    p3 = m.get("correctedDifferenceEstimate")
    return {"available": bool(m.get("hasMeasurement")),
            "p3GainCorrected": p3 if valid(p3) else None,
            "effectiveUnflagged": eff if isinstance(eff, list) and len(eff) == 2 else None,
            "bothConditionsAtLeast5": (
                all(valid(v) and v + 1e-9 >= 5 for v in eff)
                if isinstance(eff, list) and len(eff) == 2 else False)}


def describe_subject(source, mat_path, source_path):
    subject = source["subject"]
    if subject not in SUBJECTS:
        raise ValueError(f"Unexpected subject {subject}")
    if not source["completed"] or source["errors"]:
        raise ValueError(f"Unfinished EEG result {subject}")
    identity = json.loads(source_path.with_name("exploratory_identity.json").read_text())
    if str(identity.get("workflow_run")) != SOURCE_RUN or identity.get("implementation_commit") != SOURCE_SHA:
        raise ValueError(f"Unpinned original artifact {subject}")
    archive = loadmat(str(mat_path), squeeze_me=True, struct_as_record=False)["audit"]
    candidates = {int(scalar(c.id)): c for c in elements(archive.cands)}
    if len(candidates) != len(elements(archive.cands)):
        raise ValueError(f"Duplicate archived baseline candidate ID {subject}")
    reviewed = inspect(source)
    c = source["icaSeedCheck"]
    ids = c["ids"]
    n = len(ids)
    measure = grid(c["candidateMeasurements"], n, "candidateMeasurements")
    status = grid(c["rankStatus"], n, "rankStatus")
    strata = grid(c["strata"], n, "strata")
    originals = c["recommended"]
    original = originals[0]
    shadow = reviewed["shadowUnseen"]["candidateId"]
    rows = []
    for label, ident in (("originalSeed0Recommended", original),
                         ("train012Minimax", shadow)):
        if ident is None:
            rows.append({"role": label, "candidateId": None, "abstained": True})
            continue
        if ident not in ids:
            raise ValueError(f"Unrecognized ICA recipe {ident} for {subject}")
        z = candidates.get(int(ident))
        # The MAT archive has only baseline-feasible candidate records.
        # If absent, make removal status unknown, not an implicit zero.
        mechanism = removal_snapshot(z) if z is not None else {
            "seed0IcaPresent": None, "seed0RemovedComponents": None,
            "seed0RemovalRecorded": False, "seed0Mechanism": "seed0-candidate-not-archived",
            "seed0IcRemovePresent": None}
        at = ids.index(ident)
        observations = []
        for seed in range(6):
            observations.append({"seed": seed,
                 "rankStatus": status[at][seed],
                 "stratum": strata[at][seed],
                 "objective": c["objective"][at][seed],
                 **condition_summary(measure[at][seed])})
        vals = [a["p3GainCorrected"] for a in observations]
        rows.append({"role": label, "candidateId": int(ident),
            "abstained": False, **mechanism,
            "sixSeedFiniteP3Count": sum(v is not None for v in vals),
            "allSixP3Range": max(vals)-min(vals) if all(v is not None for v in vals) else None,
            "minEffectiveUnflaggedAcrossSix": [
                min(a["effectiveUnflagged"][i] for a in observations)
                if all(a["effectiveUnflagged"] is not None for a in observations)
                else None for i in range(2)],
            "allSixSeedScientificStatusFeasible": all(
                a["rankStatus"]=="feasible" for a in observations),
            "unseenSeedObservations": observations[3:]})
    return {"subject": subject, "sourceSHA256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
            "sourceRun": SOURCE_RUN, "trainRecommendations": originals[:3],
            "unseenRecommendations": originals[3:],
            "recommendedIDsStableAcrossSix": len(set(originals)) == 1,
            "originalSourceVerdict": "PARTIAL", "rows": rows,
            "icaMapMatchingAvailable": False, "waveformAcrossSeedsAvailable": False,
            "unseenScienceCorrected": {
                "originalValid": reviewed["originalBaseline"]["validAllThree"],
                "minimaxValid": reviewed["shadowUnseen"]["validAllThree"]}}


def analyze(folder):
    result_files = sorted(folder.rglob("result.json"))
    if len(result_files) != 3:
        raise ValueError(f"Expected three full real EEG result sets, got {len(result_files)}")
    records = [describe_subject(json.loads(p.read_text()),
                p.with_name("scoring_audit.mat"), p) for p in result_files]
    if {r["subject"] for r in records} != SUBJECTS:
        raise ValueError("Unexpected or duplicate development subjects")
    records.sort(key=lambda x: x["subject"])
    return {"method": "archived six-ICA-start mechanistic audit",
            "newEEGFits": 0, "newSubjects": 0, "changedProduction": False,
            "originalFrozenVerdict": "PARTIAL", "validOriginalPairs": 15,
            "componentMapComparisonsPerformed": False,
            "waveformShapeComparisonsPerformed": False,
            "subjects": records, "promotionAuthorized": False}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--evidence",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    report=analyze(a.evidence)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(report,indent=2,allow_nan=False))
    for subj in report["subjects"]:
        for record in subj["rows"]:
            print("PR51_ICA_MECHANISM",json.dumps({
                "subject":subj["subject"],"role":record["role"],
                "id":record.get("candidateId"),
                "removedComponentsSeed0":record.get("seed0RemovedComponents"),
                "mechanism":record.get("seed0Mechanism"),
                "minUnflagged":record.get("minEffectiveUnflaggedAcrossSix"),
                "sixSeedP3Range":record.get("allSixP3Range"),
                "unseenFeasible":subj["unseenScienceCorrected"][
                    "originalValid" if record["role"]=="originalSeed0Recommended" else "minimaxValid"]}))
if __name__ == "__main__":
    main()
