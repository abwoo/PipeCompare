#!/usr/bin/env python3
"""Fetch the frozen public validation subset, verify bytes, prepare MATLAB inputs.

No credentials or subject contact information are needed. Raw data are not
committed or uploaded with result artifacts. See open_data_protocol.json for
source versions, attribution, subject split and SHA-256 values.
"""
import argparse
import csv
import hashlib
import io
import json
import time
import urllib.request
import zipfile
from pathlib import Path

import numpy as np
from scipy.io import savemat


def fetch(url, path, expected):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and hashlib.sha256(path.read_bytes()).hexdigest() == expected:
        return
    for attempt in range(3):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "PipeCompare-public-validation/1"})
            with urllib.request.urlopen(request, timeout=180) as response:
                data = response.read()
            actual = hashlib.sha256(data).hexdigest()
            if actual != expected:
                raise ValueError(f"checksum mismatch for {path.name}: {actual}")
            path.write_bytes(data)
            return
        except Exception:
            if attempt == 2:
                raise
            time.sleep(2 ** attempt)


def read_edf(path):
    """Read EDF calibrated EEG; annotations are deliberately not task labels here."""
    with path.open("rb") as stream:
        header = stream.read(256)
        nbytes, nrecords = int(header[184:192]), int(header[236:244])
        duration, signals = float(header[244:252]), int(header[252:256])
        fields = stream.read(nbytes - 256)
        offset = 0

        def field(width):
            nonlocal offset
            result = [fields[offset + k * width:offset + (k + 1) * width].decode("ascii").strip()
                      for k in range(signals)]
            offset += width * signals
            return result

        labels = field(16)
        field(80)  # transducers
        units = field(8)
        pmin, pmax = np.array(field(8), float), np.array(field(8), float)
        dmin, dmax = np.array(field(8), float), np.array(field(8), float)
        field(80)  # prefilter
        sizes = np.array(field(8), int)
        field(32)
        assert nrecords > 0 and duration > 0 and offset == len(fields)
        raw = np.frombuffer(stream.read(), dtype="<i2")
    assert raw.size == nrecords * sum(sizes), "truncated or unexpected EDF length"
    raw = raw.reshape(nrecords, sum(sizes))
    starts = np.r_[0, np.cumsum(sizes)]
    keep = [i for i, name in enumerate(labels) if "Annotations" not in name]
    assert len(keep) == 64 and len(set(sizes[keep])) == 1
    rows = []
    for i in keep:
        assert units[i].lower() in ("uv", "µv"), f"unexpected EDF unit: {units[i]}"
        digital = raw[:, starts[i]:starts[i + 1]].reshape(-1).astype(float)
        rows.append((digital - dmin[i]) * (pmax[i] - pmin[i]) / (dmax[i] - dmin[i]) + pmin[i])
    return np.stack(rows), float(sizes[keep[0]] / duration), [labels[i].rstrip(".") for i in keep]


def blink_record(archive, subject, destination):
    name = f"S{subject:02d}"
    rows = list(csv.reader(io.StringIO(archive.read(name + "_data.csv").decode()), delimiter=";"))
    values = np.asarray([[float(v) for v in r[:3]] for r in rows[1:] if r], dtype=float)
    t, x = values[:, 0], values[:, 1:3].T
    dt = np.diff(t)
    assert np.all(dt > 0), "blink timestamps must be strictly increasing"
    # Timestamps are quantized on a 1024-Hz clock: median(dt) suggests
    # 256 Hz although the complete recording agrees with 250 Hz. Validate
    # the full time grid against metadata rather than using that median.
    fs = 250.0
    effective_fs = (t.size - 1) / (t[-1] - t[0])
    jitter = np.max(np.abs(t - t[0] - np.arange(t.size) / fs))
    assert jitter < 0.5 / fs, "blink timestamps disagree with the declared sampling grid"
    labels = list(csv.reader(io.StringIO(archive.read(name + "_labels.csv").decode())))
    ncorrupt = int(labels[0][1])
    corrupt = np.asarray([[float(v) for v in row[:2]] for row in labels[1:1 + ncorrupt]]).reshape(-1, 2)
    assert labels[1 + ncorrupt][0].strip() == "blinks"
    events = np.asarray([[float(v) for v in row[:2]] for row in labels[2 + ncorrupt:] if row]).reshape(-1, 2)
    # Annotation intervals, including open-ended intervals, are unavailable
    # observations. They are not quiet examples or false-positive exposure.
    valid = np.ones(t.size, bool)
    for a, b in corrupt:
        valid &= ~((t >= a) & (t <= (t[-1] if b < 0 else b)))
    indices = np.rint((events[:, 0] - t[0]) * fs).astype(int)
    include = (indices >= 0) & (indices < t.size)
    include[include] &= valid[indices[include]]
    events = events[include]
    savemat(destination / (name + ".mat"), dict(data=x, fs=fs, time=t, valid=valid,
            annotations=events, declared_fs=250.0, effective_fs=effective_fs,
            timestamp_jitter_seconds=jitter), do_compression=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--shard", type=int, required=True)
    parser.add_argument("--shards", type=int, default=6)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--components-only", action="store_true")
    args = parser.parse_args()
    assert 0 <= args.shard < args.shards
    protocol = json.loads(Path(__file__).with_name("open_data_protocol.json").read_text())
    args.out.mkdir(parents=True, exist_ok=True)
    provenance = []
    for entry in ([] if args.components_only else protocol["physionet"]["files"][args.shard::args.shards]):
        path = args.out / Path(entry["path"]).name
        url = "https://physionet-open.s3.amazonaws.com/eegmmidb/1.0.0/" + entry["path"]
        fetch(url, path, entry["sha256"])
        data, fs, labels = read_edf(path)
        savemat(path.with_suffix(".mat"), dict(data=data, fs=fs, labels=np.array(labels, dtype=object)), do_compression=True)
        provenance.append(dict(source="EEGMMIDB", subject=entry["subject"], url=url, sha256=entry["sha256"]))
        print(f"Verified EEGMMIDB subject {entry['subject']}: {data.shape}, {fs:g} Hz", flush=True)
    blink = protocol["blink"]
    url = ("https://raw.githubusercontent.com/felixshier/Blink-Detection/" + blink["mirror_commit"]
           + "/data/EEG-IO/EEG-IO.zip")
    archive = args.out / "EEG-IO.zip"
    fetch(url, archive, blink["archive_sha256"])
    with zipfile.ZipFile(archive) as z:
        for subject in blink["test_subjects"][args.shard::args.shards]:
            blink_record(z, subject, args.out)
            provenance.append(dict(source="EEG-IO", subject=subject, url=url, sha256=blink["archive_sha256"]))
    for entry in protocol["erp_core"]["files"]:
        if not (args.components_only and args.shards == 1) and f"sub-{11 + args.shard:03d}/" not in entry["path"]:
            continue
        assert entry["checksum_algorithm"] == "sha256"
        path = args.out / "erp_core" / Path(entry["path"]).name
        fetch(entry["url"], path, entry["checksum"])
        provenance.append(dict(source="ERP CORE", subject=11 + args.shard, url=entry["url"], sha256=entry["checksum"]))
    (args.out / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")


if __name__ == "__main__":
    main()
