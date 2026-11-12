#!/usr/bin/env python3
"""Independent per-event VEOG amplitude probe for frozen ERP CORE P3 excerpts.

Diagnostic only: a baseline-corrected, absolute 50-uV VEOG proxy is NOT a
clinical blink label. Does not change ICA, preprocessing or scientific gates.
Original selectedSupport is joined by excerpt event-row ID when available.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from pathlib import Path
import numpy as np
from scipy.io import loadmat

PROTOCOL="bb65c9a522e9f45210e5ee941a669630b75587a5255664aa0f1402b572e89c6b"
FROZEN_RUN="42809789847"
FIXED_VEOG_PROXY_UV=50.0

def as_labels(raw):
    return [str(x) for x in np.atleast_1d(raw).reshape(-1)]

def read_selected(path, source, subject):
    if path is None:
        return {}, None
    obj=json.loads(Path(path).read_text())
    if (obj["subject"] != subject or str(obj["source_run_id"]) != FROZEN_RUN
        or obj["protocol_sha256"] != PROTOCOL):
        raise ValueError("selectedSupport identity mismatch")
    sha=obj["source_files_sha256"]
    if (sha["set"] != source["source_set_sha256"] or
        sha["fdt"] != source["source_fdt_sha256"] or
        sha["events"] != source["source_events_sha256"]):
        raise ValueError("selectedSupport EEG/event source hash mismatch")
    choices={}
    for group in obj["selectedSupport"]:
        cond=str(group["condition"])
        fields=("referenceIds","flagged","weights","scores")
        lengths=[len(group[k]) for k in fields]
        if len(set(lengths))!=1 or lengths[0]!=group["referenceN"]:
            raise ValueError(f"misaligned saved trial vectors for {cond}")
        for event_id,flag,weight,score in zip(*(group[k] for k in fields)):
            key=(int(event_id),cond)
            if key in choices:
                raise ValueError("duplicate frozen selectedSupport event index")
            choices[key]={"hmm_flag":bool(flag),"selected_weight":float(weight),
                          "selected_p3_score":float(score) if score is not None and math.isfinite(float(score)) else None}
    return choices, obj["source_job_id"]

def analyze(input_mat: Path, source_json: Path, subject: int, selected_support=None):
    source=json.loads(source_json.read_text())
    if source["subject"]!=subject or source["protocol_sha256"]!=PROTOCOL:
        raise ValueError("ERP CORE source identity mismatch")
    data=loadmat(input_mat,squeeze_me=True,struct_as_record=False)
    fs=float(data["fs"])
    if fs!=256 or source["sampling_rate"]!=256:
        raise ValueError("wrong pinned 256-Hz excerpt rate")
    veog=np.asarray(data["veog"],dtype=float).ravel()
    latency=np.asarray(data["eventLatency"],dtype=float).ravel()
    labels=as_labels(data["eventTypes"])
    if len(latency)!=len(labels) or not all(t in ("11","31") for t in labels):
        raise ValueError("invalid event table")
    if labels.count("11")!=source["n_valid_targets"] or labels.count("31")!=source["n_valid_standards"]:
        raise ValueError("source event counts differ from frozen input")
    lookup,source_job=read_selected(selected_support,source,subject)
    offsets=np.arange(-round(0.2*fs),round(0.8*fs)+1,dtype=int)
    mask=offsets<0
    rows=[]
    joined=set()
    for i,(loc,condition) in enumerate(zip(latency,labels),start=1):
        if not math.isfinite(loc) or loc!=int(loc):
            raise ValueError("non-integral excerpt latency")
        positions=(int(loc)-1)+offsets # MATLAB 1-based to Python 0-based.
        if positions.min()<0 or positions.max()>=veog.size:
            raise ValueError(f"incomplete requested epoch for event {i}")
        seg=veog[positions]
        base=float(np.mean(seg[mask]))
        vmax=float(np.max(np.abs(seg-base)))
        if not math.isfinite(vmax):
            raise ValueError(f"non-finite VEOG proxy for event {i}")
        ref=lookup.get((i,condition))
        if ref is not None: joined.add((i,condition))
        rows.append({"subject":subject,"event_index_1based":i,"condition":condition,
                    "veog_peak_abs_baseline_corrected_uv":round(vmax,9),
                    "veog_exceeds_frozen_diagnostic_50uv":vmax>FIXED_VEOG_PROXY_UV,
                    "selected_standard_hmm_flag":ref["hmm_flag"] if ref is not None else None,
                    "selected_standard_weight":ref["selected_weight"] if ref is not None else None,
                    "selected_standard_p3_score":ref["selected_p3_score"] if ref is not None else None})
    if len(joined)!=len(lookup):
        raise ValueError("reference event mismatch between raw VEOG and frozen selectedSupport")
    summary=[]
    for cond in ("11","31"):
        r=[z for z in rows if z["condition"]==cond]
        labelled=[z for z in r if z["selected_standard_hmm_flag"] is not None]
        hf=lambda f:sum(1 for v in labelled if v["selected_standard_hmm_flag"] is f)
        cross=lambda flag,amp:sum(1 for v in labelled if v["selected_standard_hmm_flag"] is flag and v["veog_exceeds_frozen_diagnostic_50uv"] is amp)
        summary.append({"condition":cond,"events":len(r),"veog_proxy_above50":sum(v["veog_exceeds_frozen_diagnostic_50uv"] for v in r),
                        "veog_proxy_median_uv":float(np.median([v["veog_peak_abs_baseline_corrected_uv"] for v in r])),
                        "selected_standard_labelled":len(labelled),
                        "hmm_flagged":hf(True),"hmm_unflagged":hf(False),
                        "hmm_flagged_and_veog_above50":cross(True,True),
                        "hmm_flagged_and_veog_below_equal50":cross(True,False),
                        "hmm_unflagged_and_veog_above50":cross(False,True),
                        "hmm_unflagged_and_veog_below_equal50":cross(False,False)})
    output={"study":"PR51 frozen ERP CORE P3 per-event VEOG proxy; purely diagnostic",
            "subject":subject,"protocol_sha256":PROTOCOL,"source_sha256":{"set":source["source_set_sha256"],
              "fdt":source["source_fdt_sha256"],"events":source["source_events_sha256"]},
            "original_selected_support_run":FROZEN_RUN if source_job else None,
            "original_selected_support_job":source_job,
            "external_proxy_threshold_uv":FIXED_VEOG_PROXY_UV,
            "proxy_is_not_ground_truth":True,"n_events":len(rows),
            "conditions":summary,"events":rows}
    return output

def save(result, out):
    out.mkdir(parents=True,exist_ok=True)
    filename=f"subject_{result['subject']:03d}"
    (out/f"{filename}.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    fields=["subject","event_index_1based","condition","veog_peak_abs_baseline_corrected_uv",
            "veog_exceeds_frozen_diagnostic_50uv","selected_standard_hmm_flag",
            "selected_standard_weight","selected_standard_p3_score"]
    with (out/f"{filename}.csv").open("w",newline="") as file:
        writer=csv.DictWriter(file,fieldnames=fields);writer.writeheader();writer.writerows(result["events"])

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--input",type=Path,required=True)
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--subject",type=int,required=True)
    p.add_argument("--selected-support",type=Path)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    data=analyze(a.input,a.source,a.subject,a.selected_support)
    save(data,a.out)
    print("PR51_VEOG_PROBE",json.dumps({"subject":a.subject,"events":data["n_events"],
        "cond":data["conditions"],"source_job":data["original_selected_support_job"],
        "not_blink_ground_truth":True},separators=(",",":")))
if __name__=="__main__":
    main()
