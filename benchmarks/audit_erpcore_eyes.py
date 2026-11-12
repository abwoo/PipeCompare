#!/usr/bin/env python3
"""Non-prescriptive trial-level eye/EEG amplitude audit on real ERP CORE P3.

No modification of frozen candidate pipeline, events, acceptance, training,
ICA or artifact flags. VEOG is a signal proxy, NOT independent blink truth.
Waveforms are baseline centered only; these summaries are NOT comparable
to the filtered ±100-uV MATLAB rejection rule or algorithm HMM.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.io import loadmat

EXPECTED_PROTOCOL = "bb65c9a522e9f45210e5ee941a669630b75587a5255664aa0f1402b572e89c6b"

def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for data in iter(lambda: f.read(2 ** 20), b""):
            h.update(data)
    return h.hexdigest()


def one_text(v):
    while isinstance(v, np.ndarray) and v.size == 1:
        v = v.item()
    if isinstance(v, bytes):
        return v.decode()
    if isinstance(v, np.ndarray):
        return "".join(str(x) for x in v.flat)
    return str(v)


def quantile(x):
    arr = np.asarray(x, float)
    return {"min": float(np.min(arr)), "p25": float(np.quantile(arr, .25)),
            "median": float(np.median(arr)), "p75": float(np.quantile(arr, .75)),
            "p90": float(np.quantile(arr, .90)), "max": float(np.max(arr))}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()

    source = json.loads(args.source.read_text())
    assert source["protocol_sha256"] == EXPECTED_PROTOCOL, "wrong source/protocol"
    assert source["dataset"] == "ERP CORE" and source["task"] == "P3"
    d = loadmat(args.input, squeeze_me=True, struct_as_record=False)
    fs = int(float(d["fs"]))
    assert fs == 256
    labels = [one_text(x) for x in np.atleast_1d(d["labels"]).ravel()]
    assert labels.count("Pz") == 1 and len(labels) == 16
    eeg = np.asarray(d["data"], float)
    assert eeg.shape[0] == len(labels)
    pz = eeg[labels.index("Pz"), :]
    veog = np.asarray(d["veog"], float).reshape(-1)
    heog = np.asarray(d["heog"], float).reshape(-1)
    assert pz.shape == veog.shape == heog.shape and pz.size == 300 * fs
    codes = [one_text(x) for x in np.atleast_1d(d["eventTypes"]).reshape(-1)]
    latencies = np.asarray(d["eventLatency"], float).reshape(-1)
    assert len(codes) == len(latencies)
    assert all(x in ("11", "31") for x in codes)
    assert sum(x == "11" for x in codes) == source["n_valid_targets"]
    assert sum(x == "31" for x in codes) == source["n_valid_standards"]

    sample_offsets = np.arange(-round(.2*fs), round(.8*fs)+1, dtype=int)
    baseline_mask = sample_offsets < 0
    window_mask = (sample_offsets >= round(.3*fs)) & (sample_offsets <= round(.6*fs))
    trials = []
    for index, (code, lat) in enumerate(zip(codes, latencies), 1):
        epoch = int(round(lat)) - 1 + sample_offsets
        if epoch[0] < 0 or epoch[-1] >= len(veog):
            raise ValueError(f"incomplete event epoch {index}")
        metrics = {"subject":source["subject"],"trial_index":index,
                   "condition":"target" if code=="11" else "standard",
                   "event_sample_1based":int(round(lat)), "event_time_s":float((lat-1)/fs)}
        for name, trace in (("veog",veog), ("heog",heog),("pz",pz)):
            e = trace[epoch]
            e = e - np.mean(e[baseline_mask])
            metrics[name + "_epoch_maxabs_uv"] = float(np.max(np.abs(e)))
            metrics[name + "_epoch_p2p_uv"] = float(np.ptp(e))
            metrics[name + "_window_mean_uv"] = float(np.mean(e[window_mask]))
            metrics[name + "_baseline_sd_uv"] = float(np.std(e[baseline_mask]))
        # Whole-EEG trial rejection is driven by the maximum across EEG
        # channels, not solely by Pz. These raw, epoch-baseline-corrected
        # peaks are descriptive and are NOT the filtered MATLAB thresholds.
        ch_epoch = eeg[:, epoch] - np.mean(eeg[:, epoch][:, baseline_mask], axis=1, keepdims=True)
        ch_abs = np.max(np.abs(ch_epoch), axis=1)
        metrics["eeg_max_channel_uv"] = float(np.max(ch_abs))
        metrics["eeg_max_channel_label"] = labels[int(np.argmax(ch_abs))]
        metrics["eeg_channels_gt100_uv"] = "|".join([labels[k] for k in np.flatnonzero(ch_abs > 100)])
        metrics["eeg_any_channel_gt100_uv"] = int(np.any(ch_abs > 100))
        metrics["eeg_n_channels_gt100_uv"] = int(np.count_nonzero(ch_abs > 100))
        metrics["fp1_maxabs_uv"] = float(ch_abs[labels.index("Fp1")])
        metrics["fp2_maxabs_uv"] = float(ch_abs[labels.index("Fp2")])
        trials.append(metrics)

    groups = {}
    for label in ("target","standard"):
        rows = [r for r in trials if r["condition"]==label]
        groups[label] = {
            "n":len(rows),
            "veog_epoch_maxabs_uv":quantile([r["veog_epoch_maxabs_uv"] for r in rows]),
            "heog_epoch_maxabs_uv":quantile([r["heog_epoch_maxabs_uv"] for r in rows]),
            "pz_epoch_maxabs_uv":quantile([r["pz_epoch_maxabs_uv"] for r in rows]),
            "veog_above_50uv_fraction":sum(r["veog_epoch_maxabs_uv"]>50 for r in rows)/len(rows),
            "veog_above_100uv_fraction":sum(r["veog_epoch_maxabs_uv"]>100 for r in rows)/len(rows),
            "pz_above_100uv_fraction":sum(r["pz_epoch_maxabs_uv"]>100 for r in rows)/len(rows),
            "eeg_max_channel_uv":quantile([r["eeg_max_channel_uv"] for r in rows]),
            "eeg_any_channel_above_100uv_fraction":sum(r["eeg_any_channel_gt100_uv"] for r in rows)/len(rows),
            "fp1_above_100uv_fraction":sum(r["fp1_maxabs_uv"]>100 for r in rows)/len(rows),
            "fp2_above_100uv_fraction":sum(r["fp2_maxabs_uv"]>100 for r in rows)/len(rows),
            "eeg_channel_gt100_counts":{label:sum(label in r["eeg_channels_gt100_uv"].split("|") for r in rows)
                                         for label in labels},
        }
    args.out.mkdir(parents=True,exist_ok=True)
    with (args.out/"trial_amplitude_audit.csv").open("w",newline="") as fp:
        wr = csv.DictWriter(fp,fieldnames=list(trials[0]))
        wr.writeheader()
        wr.writerows(trials)
    result={
        "subject":source["subject"],
        "source_protocol_sha256":source["protocol_sha256"],
        "source_set_sha256":source["source_set_sha256"],
        "source_fdt_sha256":source["source_fdt_sha256"],
        "source_events_sha256":source["source_events_sha256"],
        "input_mat_sha256":sha(args.input),
        "purpose":"Descriptive VEOG/HEOG/16-EEG-electrode distributions; no optical ground truth",
        "algorithm_changed":False,
        "filter_applied":"NONE; per-epoch prestimulus baseline subtraction ONLY",
        "not_equivalent_to_frozen_MATLAB_rejection":True,
        "erpcore_300s_fs":fs,
        "n_trials":len(trials),
        "by_condition":groups
    }
    (args.out/"audit.json").write_text(json.dumps(result,indent=2)+"\n")
    print("PR51_REAL_EYE_SIGNAL_AUDIT",json.dumps({"subject":source["subject"],"groups":groups}))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
