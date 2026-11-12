#!/usr/bin/env python3
"""Independent VEOG/HEOG cross-check of ERP CORE trial-support records.

Diagnostic only: never modifies EEG, ICA, rank, epoch masks, or frozen
thresholds. Uses source event row IDs (1-based MATLAB reference IDs), and
explicitly normalizes '11'/'31' as target/standard to avoid the pilot's
old target-label mismatch. A 50 uV proxy is NOT a clinical blink ground truth.
"""
import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
from scipy.io import loadmat

EYE_PROXY_UV = 50.0
ARMS = ("erp", "minimal", "common")
COND_NAMES = {"11": "target", "31": "standard",
              "target": "target", "standard": "standard"}


def canonical(value):
    name = str(value)
    if name not in COND_NAMES:
        raise ValueError(f"unexpected event condition label: {value!r}")
    return COND_NAMES[name]


def effective(weights):
    w = [float(z) for z in weights if math.isfinite(float(z)) and z > 0]
    s = sum(w)
    ss = sum(x*x for x in w)
    return s*s/ss if ss > 0 else 0.0


def eye_max(signal, latency, sample_rate):
    """MATLAB: epochSamp=-round(.2*fs):round(.8*fs); baseline epochSamp<0."""
    left = int(np.floor(.2*sample_rate+.5))
    right = int(np.floor(.8*sample_rate+.5))
    offsets = np.arange(-left, right+1, dtype=int)
    idx = int(latency)-1 + offsets
    if idx[0] < 0 or idx[-1] >= len(signal):
        raise ValueError("support trial has out-of-range EOG epoch")
    segment = signal[idx].astype(float)
    segment -= np.mean(segment[offsets < 0])
    return float(np.max(np.abs(segment)))


def crosscheck(record, mat):
    """Return two-condition arm summaries and aligned per-trial proxy evidence."""
    if not record.get("completed") or record.get("notRun", 0) != 0:
        raise ValueError("cannot inspect incomplete real-data run")
    n = int(record["subject"])
    fs = float(np.asarray(mat["fs"]).item())
    latencies = np.asarray(mat["eventLatency"]).reshape(-1)
    event_types = [str(x).strip() for x in np.asarray(mat["eventTypes"]).reshape(-1)]
    veog = np.asarray(mat["veog"]).reshape(-1).astype(float)
    heog = np.asarray(mat["heog"]).reshape(-1).astype(float)
    if len(event_types) != len(latencies):
        raise ValueError("event latency/type mismatch")
    if fs <= 0 or len(heog) != len(veog):
        raise ValueError("invalid EOG source identity")
    rows, trials = [], []
    for arm in ARMS:
        source = record.get("selectedSupport") if arm == "erp" else (record.get(arm) or {}).get("support")
        if not source:
            for condition in ("target", "standard"):
                rows.append({"subject": n, "arm": arm, "condition": condition,
                             "available": False, "reason": "support-unavailable"})
            continue
        seen = set()
        for condition_data in source:
            name = canonical(condition_data["condition"])
            if name in seen:
                raise ValueError(f"duplicate support condition {arm}/{name}")
            seen.add(name)
            ids = condition_data["referenceIds"]
            weights = condition_data["weights"]
            flagged = condition_data["flagged"]
            scores = condition_data["scores"]
            if not all(len(seq) == len(ids) for seq in (weights, flagged, scores)):
                raise ValueError(f"misaligned support vectors {arm}/{name}")
            ef_hmm, ef_veog, proxy_positive, hmm_positive, matching = [], [], 0, 0, 0
            proxy_on_hmm = 0
            for ident, score, w, flag in zip(ids, scores, weights, flagged):
                j = int(ident)
                if j != ident or not 1 <= j <= len(latencies):
                    raise ValueError("invalid 1-based source event reference")
                if canonical(event_types[j-1]) != name:
                    raise ValueError(f"source event/condition mismatch subject={n} id={j}: {name}")
                eye = eye_max(veog, latencies[j-1], fs)
                horizontal = eye_max(heog, latencies[j-1], fs)
                proxy = eye > EYE_PROXY_UV
                sc = float(score) if score is not None else float("nan")
                weight = float(w) if w is not None else float("nan")
                good = math.isfinite(sc) and math.isfinite(weight) and weight > 0
                positive_w = weight if good else 0.0
                proxy_positive += int(proxy)
                hmm_positive += int(bool(flag))
                matching += int(bool(flag) == proxy)
                proxy_on_hmm += int(proxy and bool(flag))
                ef_hmm.append(positive_w if not flag else 0.0)
                ef_veog.append(positive_w if not proxy else 0.0)
                trials.append({"subject": n, "arm": arm, "condition": name,
                               "source_event_id": j, "score": score,
                               "weight": positive_w, "hmm_blink_proxy": bool(flag),
                               "veog_max_uv": eye, "heog_max_uv": horizontal,
                               "veog_above_50uv": proxy})
            frozen_eff = float(condition_data["effectiveUnflaggedN"])
            if abs(effective(ef_hmm)-frozen_eff) > 1e-7:
                raise ValueError(f"frozen HMM support mismatch {n}/{arm}/{name}")
            rows.append({"subject": n, "arm": arm, "condition": name, "available": True,
                         "referenceN": int(condition_data["referenceN"]),
                         "kept": int(condition_data["kept"]),
                         "hmm_flagged": hmm_positive, "veog_proxy_positive": proxy_positive,
                         "both_positive": proxy_on_hmm, "proxy_agreement": matching,
                         "hmm_effective_unflagged": frozen_eff,
                         "veog_proxy_effective_unflagged": effective(ef_veog),
                         "above_original_floor_hmm": frozen_eff >= 5,
                         "reason": "diagnostic-only-no-gate-change"})
        if seen != {"target", "standard"}:
            raise ValueError(f"expected both condition identities for {n}/{arm}; got {seen}")
    return rows, trials


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--result", required=True, type=Path)
    ap.add_argument("--mat", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()
    result = json.loads(a.result.read_text())
    mat = loadmat(a.mat, simplify_cells=True)
    # scipy simplify_cells may return eventTypes as a flat string array.
    rows, trial_rows = crosscheck(result, mat)
    a.out.mkdir(parents=True, exist_ok=True)
    for name, data in (("veog_summary.csv", rows), ("veog_trials.csv", trial_rows)):
        if data:
            keys = list(dict.fromkeys(key for row in data for key in row))
            with (a.out/name).open("w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=keys)
                w.writeheader()
                w.writerows(data)
    (a.out/"veog_summary.json").write_text(json.dumps({
        "scope": "real-P3-exploratory-condition-matched-eye-proxy-diagnostic",
        "source_subject": result["subject"], "proxy_threshold_uv": EYE_PROXY_UV,
        "proxy_is_ground_truth": False, "scientific_gates_changed": False,
        "conditions": rows, "trial_count": len(trial_rows)}, indent=2))
    print("PR51_VEOG_CROSSCHECK", json.dumps({"subject": result["subject"], "rows": rows}))


if __name__ == "__main__":
    main()
