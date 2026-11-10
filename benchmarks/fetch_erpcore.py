#!/usr/bin/env python3
"""Fetch one ERP CORE P3 participant and prepare a PipeCompare input.

Per-subject download + SHA256 verification from the frozen protocol pins,
then: parse the EEGLAB .set/.fdt, downsample 1024 -> 256 Hz, select the
16 PR51 channels plus VEOG/HEOG/P9/P10, cut the frozen excerpt, keep only
correct trials with RT 200-1000 ms, relabel stimulus events to '11'
(target) / '31' (standard) so the Standard search runs on exactly the
PR51 pilot's code path, and write input.mat + source.json.

Relabeling rationale (frozen protocol section 6): pop_pipecompare 'events'
accepts one condition per event type and cannot express "target = 5 codes,
standard = 20 codes". The original codes, correctness and RTs are preserved
in source.json and in input.mat's trial table, so the relabeling is fully
auditable. No other trial filtering is applied here; artifact handling is
left to the pipeline under test.
"""
import argparse
import csv
import hashlib
import json
import os
import sys
import urllib.request
from pathlib import Path

import numpy as np
from scipy.io import loadmat
from scipy.signal import resample_poly

CHANNELS_16 = ["Fp1", "Fp2", "F7", "F3", "Fz", "F4", "F8",
               "C3", "Cz", "C4", "P3", "Pz", "P4", "O1", "Oz", "O2"]
# ERP CORE channel labels are uppercase in the .set file.
VEOG_LABEL = "VEOG_lower"
HEOG_LEFT = "HEOG_left"
HEOG_RIGHT = "HEOG_right"
REF_P9 = "P9"
REF_P10 = "P10"


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url, dest, expect_sha256):
    req = urllib.request.Request(url, headers={"User-Agent": "PipeCompare-erp-core-fetch/1.0"})
    with urllib.request.urlopen(req, timeout=120) as r, open(dest, "wb") as f:
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            f.write(chunk)
    got = sha256_file(dest)
    if got != expect_sha256:
        raise SystemExit(f"SHA256 mismatch for {dest.name}: got {got}, expected {expect_sha256}")
    return got


def parse_events_tsv(path):
    """Return list of dicts: sample (int, 1-based EEGLAB), onset_s, value, trial_type."""
    rows = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            rows.append({
                "onset": float(r["onset"]),
                "sample": int(r["sample"]),
                "trial_type": r["trial_type"],
                "value": r["value"],
                "stim_file": r.get("stim_file", ""),
            })
    return rows


def trial_filter(events, excerpt_seconds, rt_lo=0.2, rt_hi=1.0):
    """Keep correct trials with RT in [rt_lo, rt_hi] inside the excerpt.

    ERP CORE P3: stimulus events have two-digit codes; tens digit = block
    target letter (1=A..5=E), ones digit = the shown letter. Target when
    tens == ones (11, 22, 33, 44, 55); otherwise standard. A response event
    follows each stimulus: value 201 = correct, 202 = error. RT = response
    onset - stimulus onset (seconds).
    Returns (kept, dropped_counts) where kept is a list of dicts with
    original code, is_target, rt, onset_s, sample (1-based at 1024 Hz).
    """
    dropped = {"no_response": 0, "incorrect": 0, "rt_outside": 0,
               "outside_excerpt": 0, "epoch_incomplete": 0}
    kept = []
    n = len(events)
    for i, ev in enumerate(events):
        if ev["trial_type"] != "stimulus":
            continue
        code = ev["value"]
        if len(code) != 2 or not code.isdigit():
            continue
        onset = ev["onset"]
        if onset < 0.2 or onset > excerpt_seconds - 0.8:
            dropped["outside_excerpt" if onset > excerpt_seconds - 0.8 else "epoch_incomplete"] += 1
            continue
        # First response event after this stimulus.
        resp = None
        for j in range(i + 1, n):
            if events[j]["trial_type"] == "response":
                resp = events[j]
                break
            if events[j]["trial_type"] == "stimulus":
                break
        if resp is None:
            dropped["no_response"] += 1
            continue
        if resp["value"] != "201":
            dropped["incorrect"] += 1
            continue
        rt = resp["onset"] - onset
        if not (rt_lo <= rt <= rt_hi):
            dropped["rt_outside"] += 1
            continue
        kept.append({
            "code": code,
            "is_target": code[0] == code[1],
            "rt": rt,
            "onset": onset,
            "sample1024": ev["sample"],
        })
    return kept, dropped


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subject", type=int, required=True)
    ap.add_argument("--out", required=True, help="output directory for input.mat + source.json")
    ap.add_argument("--protocol", default=None,
                    help="path to erpcore_protocol.json (default: benchmarks/erpcore_protocol.json next to this script)")
    ap.add_argument("--cache", default=None, help="optional directory with pre-downloaded source files")
    args = ap.parse_args()

    script_dir = Path(__file__).resolve().parent
    protocol_path = Path(args.protocol) if args.protocol else script_dir / "erpcore_protocol.json"
    protocol_sha256 = hashlib.sha256(protocol_path.read_bytes()).hexdigest()
    p = json.loads(protocol_path.read_text())
    entry = next((f for f in p["source"]["files"] if f["subject"] == args.subject), None)
    if entry is None:
        raise SystemExit(f"Subject {args.subject} is not in the frozen protocol cohort.")

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    raw = out / "raw"
    raw.mkdir(exist_ok=True)

    files = {}
    for kind in ("set", "fdt", "events"):
        info = entry[kind]
        dest = raw / info["name"]
        if args.cache:
            cached = Path(args.cache) / f"sub-{args.subject:03d}" / info["name"]
            if cached.exists():
                import shutil
                shutil.copyfile(cached, dest)
                got = sha256_file(dest)
                if got != info["sha256"]:
                    raise SystemExit(f"SHA256 mismatch for cached {dest.name}")
                files[kind] = {"path": str(dest), "sha256": got, "size": dest.stat().st_size}
                continue
        got = download(info["url"], dest, info["sha256"])
        files[kind] = {"path": str(dest), "sha256": got, "size": dest.stat().st_size}
        print(f"{info['name']}: OK sha256={got[:16]}...")

    # ---- Parse the EEGLAB .set ----
    mat = loadmat(files["set"]["path"], squeeze_me=True, struct_as_record=False)
    eeg = mat["EEG"]
    data_file = str(eeg.data)  # .fdt filename recorded in the .set
    fdt_path = raw / data_file
    if fdt_path.name != Path(files["fdt"]["path"]).name:
        raise SystemExit(f".set points at {data_file} but protocol pins {files['fdt']['path']}")
    nbchan = int(eeg.nbchan)
    pnts = int(eeg.pnts)
    srate = int(eeg.srate)
    assert srate == p["source"]["sampling_rate_hz"], f"expected 1024 Hz, got {srate}"
    labels = [str(c.labels) for c in np.atleast_1d(eeg.chanlocs)]
    fdt = np.fromfile(fdt_path, dtype="<f4")
    assert fdt.size == nbchan * pnts, f".fdt size {fdt.size} != {nbchan}x{pnts}"
    # EEGLAB writes .fdt column-major (time-major on disk); a C-order
    # reshape scrambles channels into white noise (verified on sub-001:
    # flat spectrum above 30 Hz, no power below). Fortran order restores
    # real EEG (1/f spectrum, neighbour correlation ~0.98).
    data1024 = fdt.reshape((nbchan, pnts), order="F").astype(np.float64)
    # Cross-check event count against events.tsv.
    n_eeg_events = len(np.atleast_1d(eeg.event))
    events = parse_events_tsv(files["events"]["path"])
    assert n_eeg_events == len(events), f"EEG events {n_eeg_events} != events.tsv rows {len(events)}"

    # ---- Channel selection ----
    lower = {lab.lower(): i for i, lab in enumerate(labels)}
    def idx(name):
        if name.lower() not in lower:
            raise SystemExit(f"Channel {name} missing from recording")
        return lower[name.lower()]
    sel16 = [idx(c) for c in CHANNELS_16]
    veog_i = idx(VEOG_LABEL)
    heog_i = [idx(HEOG_LEFT), idx(HEOG_RIGHT)]
    ref_i = [idx(REF_P9), idx(REF_P10)]

    # ---- Downsample 1024 -> 256 Hz on the full recording, then cut ----
    fs_out = p["excerpt"]["sampling_rate"]
    assert srate % fs_out == 0
    factor = srate // fs_out
    data256 = resample_poly(data1024, 1, factor, axis=1)
    n256 = data256.shape[1]
    excerpt_s = p["excerpt"]["duration_seconds"]
    n_excerpt = int(excerpt_s * fs_out)
    assert n256 >= n_excerpt, f"recording too short: {n256/fs_out:.1f}s < {excerpt_s}s"
    excerpt = data256[:, :n_excerpt].copy()
    # DC offset removal over the excerpt (as in the PR51 pilot's inject()).
    excerpt -= excerpt.mean(axis=1, keepdims=True)

    def trace(rows):
        return excerpt[rows, :]

    data16 = trace(sel16)
    veog = trace([veog_i])[0]
    heog = trace(heog_i)[0] - trace(heog_i)[1]  # bipolar, paper definition
    ref_p9 = trace([ref_i[0]])[0]
    ref_p10 = trace([ref_i[1]])[0]

    # ---- Trial selection ----
    kept, dropped = trial_filter(events, excerpt_s,
                                 p["trial_rules"]["rt_min_seconds"],
                                 p["trial_rules"]["rt_max_seconds"])
    n_target = sum(1 for t in kept if t["is_target"])
    n_standard = len(kept) - n_target
    if n_target < p["trial_rules"]["min_valid_targets"]:
        raise SystemExit(f"Only {n_target} valid targets (< {p['trial_rules']['min_valid_targets']})")
    if n_standard < p["trial_rules"]["min_valid_standards"]:
        raise SystemExit(f"Only {n_standard} valid standards (< {p['trial_rules']['min_valid_standards']})")

    # Event latencies at 256 Hz, 1-based EEGLAB samples.
    lat256 = [int(round((t["sample1024"] - 1) / factor)) + 1 for t in kept]
    types = ["11" if t["is_target"] else "31" for t in kept]

    # ---- Write input.mat ----
    from scipy.io import savemat
    trial_table = {
        "originalCode": np.array([t["code"] for t in kept], dtype=object),
        "isTarget": np.array([t["is_target"] for t in kept], dtype=bool),
        "rtSeconds": np.array([t["rt"] for t in kept]),
        "onsetSeconds": np.array([t["onset"] for t in kept]),
        "sample1024": np.array([t["sample1024"] for t in kept]),
    }
    savemat(out / "input.mat", {
        "data": data16,                      # 16 x (duration*256), uV, demeaned
        "fs": float(fs_out),
        "labels": np.array(CHANNELS_16, dtype=object),
        "veog": veog,                        # VEOG_lower, demeaned
        "heog": heog,                        # HEOG_left - HEOG_right, demeaned
        "refP9": ref_p9,                     # P9, demeaned (for the reference computation)
        "refP10": ref_p10,                   # P10, demeaned
        "eventTypes": np.array(types, dtype=object),   # '11' target / '31' standard
        "eventLatency": np.array(lat256),     # 1-based samples at 256 Hz
        "eventOnsetSeconds": np.array([t["onset"] for t in kept]),
        "trials": trial_table,
    })
    source = {
        "dataset": "ERP CORE",
        "task": "P3",
        "subject": args.subject,
        "source_set_sha256": files["set"]["sha256"],
        "source_fdt_sha256": files["fdt"]["sha256"],
        "source_events_sha256": files["events"]["sha256"],
        "source_set_url": entry["set"]["url"],
        "source_fdt_url": entry["fdt"]["url"],
        "source_events_url": entry["events"]["url"],
        "license": p["source"]["license"],
        "sampling_rate": fs_out,
        "samples": n_excerpt,
        "offset_seconds": p["excerpt"]["offset_seconds"],
        "duration_seconds": excerpt_s,
        "channels": CHANNELS_16,
        "eog_channels": {"veog": VEOG_LABEL, "heog": "HEOG_left - HEOG_right"},
        "reference_channels": [REF_P9, REF_P10],
        "protocol_sha256": protocol_sha256,
        "n_events_total": len(events),
        "n_stimuli_in_excerpt": sum(1 for e in events
                                     if e["trial_type"] == "stimulus"
                                     and 0.2 <= e["onset"] <= excerpt_s - 0.8),
        "n_valid_targets": n_target,
        "n_valid_standards": n_standard,
        "dropped": dropped,
        "relabeling": "target codes 11/22/33/44/55 -> '11'; all other two-digit stimulus codes -> '31'; response events dropped after correctness/RT filtering; original codes in input.mat trials.originalCode",
    }
    (out / "source.json").write_text(json.dumps(source, indent=2))
    print(json.dumps({"subject": args.subject, "n_valid_targets": n_target,
                      "n_valid_standards": n_standard, "dropped": dropped}, indent=2))


if __name__ == "__main__":
    main()
