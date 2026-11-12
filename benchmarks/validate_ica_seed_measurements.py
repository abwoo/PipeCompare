#!/usr/bin/env python3
"""Validate and summarize per-seed *gain-corrected* ERP CORE P300 measurement.

Read-only exploratory evidence. Candidate ID and raw P3 amplitude are not
interchangeable with calibrated P3 amplitude; missing data are never imputed.
"""
import argparse
import json
import math
from pathlib import Path


def valid_number(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def validate(subject_result, atol=1e-8):
    if not subject_result.get("completed") or subject_result.get("errors"):
        raise ValueError("uncompleted P3 run cannot be reported as seed reproducibility")
    s = subject_result.get("icaSeedCheck")
    if not isinstance(s, dict):
        raise ValueError("ICA seed check is unavailable")
    rec = s.get("recommended")
    measures = s.get("measurements")
    if not (isinstance(rec, list) and isinstance(measures, list)
            and len(rec) == len(measures) == 3):
        raise ValueError("expected baseline + two pinned ICA starts")
    out = []
    for j, (idx, m) in enumerate(zip(rec, measures)):
        if int(m["seed"]) != j:
            raise ValueError("unexpected ICA seed index")
        if idx is None:
            if m["candidateId"] is not None:
                raise ValueError("cannot assign a selected candidate to no recommendation")
            out.append({"seed": j, "available": False, "candidateId": None})
            continue
        if int(m["candidateId"]) != int(idx):
            raise ValueError("seed-selected candidate ID does not match ranking")
        estimates = m.get("conditionEstimates")
        raw = m.get("differenceEstimate")
        gain = m.get("gain")
        corrected = m.get("correctedDifferenceEstimate")
        sme = m.get("sme")
        if not m.get("hasMeasurement"):
            out.append({"seed": j, "available": False, "candidateId": idx})
            continue
        if not (isinstance(estimates, list) and len(estimates) == 2
                and all(valid_number(v) for v in estimates)
                and valid_number(raw)):
            raise ValueError("invalid two-condition P300 point estimate")
        if not math.isclose(estimates[0]-estimates[1], raw, rel_tol=0, abs_tol=atol):
            raise ValueError("raw P300 difference was not derived from condition estimates")
        if not valid_number(gain) or gain <= 0:
            raise ValueError("valid P300 estimates require positive per-seed calibration gain")
        if not valid_number(corrected):
            raise ValueError("gain-corrected P300 difference is missing")
        if not math.isclose(raw/gain, corrected, rel_tol=0, abs_tol=atol):
            raise ValueError("P300 gain-corrected difference is inconsistent")
        effective = m.get("effectiveUnflagged")
        source = m.get("sourceUnflagged")
        if not (isinstance(effective,list) and isinstance(source,list)
                and len(effective)==len(source)==2):
            raise ValueError("per-condition trial support is not aligned")
        for upper, n in zip(source, effective):
            if valid_number(upper) and valid_number(n):
                if n < -atol or n > upper + atol:
                    raise ValueError("effective support must not exceed unflagged source trials")
        out.append({"seed":j, "available":True, "candidateId":idx,
                    "rawP3Difference":raw, "correctedP3Difference":corrected,
                    "gain":gain, "sme":sme, "effectiveUnflagged":effective})
    baseline = out[0]
    erp = subject_result.get("erp") or {}
    if erp.get("available"):
        if not baseline.get("available"):
            raise ValueError("available ERP CORE readout lacks baseline ICA seed estimate")
        if not valid_number(erp.get("differenceAmplitude")) or not math.isclose(
                erp["differenceAmplitude"],baseline["correctedP3Difference"],
                rel_tol=0,abs_tol=atol):
            raise ValueError("baseline gain-corrected P300 differs from ERP CORE report")
    corrected = [r.get("correctedP3Difference") for r in out]
    delta = [v-corrected[0] if valid_number(v) and valid_number(corrected[0])
             else None for v in corrected]
    return {"subject": int(subject_result["subject"]),
            "candidateIds":rec,"stableCandidateId":bool(s["same"]),
            "bySeed":out,"correctedDeltaFromFirst":delta,
            "noScientificGatesChanged":True,
            "trialShapeAndRemovedICStillNeedSeparateAudit":True}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--inputs", type=Path, nargs="+", required=True)
    p.add_argument("--out", type=Path, required=True)
    a=p.parse_args()
    results = [validate(json.loads(path.read_text())) for path in a.inputs]
    subjects=[x["subject"] for x in results]
    if len(set(subjects)) != len(subjects):
        raise ValueError("duplicated ERP CORE subject")
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps({
        "exploratory":True,"confirmatory":False,
        "reportedScale":"gain-corrected P3 mean amplitude difference in microvolts",
        "subjects":results},indent=2))
    for record in results:
        print("PR51_ICA_MEASUREMENT_COMPARE",json.dumps(record))


if __name__=="__main__":
    main()
