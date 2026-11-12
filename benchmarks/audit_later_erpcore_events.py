#!/usr/bin/env python3
"""Read-only ERP CORE trial availability beyond the original frozen 300 s.

Downloads only SHA256-pinned .set metadata and events.tsv, NOT raw EEG .fdt.
Counts correct task events; makes no claim that extra events are clean or
usable under the unchanged HMM, ICLabel and effective-n science constraints.
"""
import argparse
import hashlib
import json
import math
import urllib.request
from pathlib import Path

from scipy.io import loadmat
from fetch_erpcore import parse_events_tsv, trial_filter

SUBJECTS=(4,8)
WINDOWS=(300,360,420,480)


def checked_download(url,expected,path):
    req=urllib.request.Request(url,headers={"User-Agent":"PipeCompare-trial-support-only/1.0"})
    with urllib.request.urlopen(req,timeout=90) as resp, path.open("wb") as f:
        for block in iter(lambda:resp.read(1<<20),b""):
            f.write(block)
    digest=hashlib.sha256(path.read_bytes()).hexdigest()
    if digest!=expected:
        raise ValueError(f"Source SHA256 mismatch for {path.name}")
    return digest


def count_windows(events,seconds,windows=WINDOWS,rt=(.2,1.0)):
    """Separate task-event availability from EEG artifact eligibility."""
    if not math.isfinite(seconds) or seconds<=300:
        raise ValueError("Full EEG must extend beyond frozen 300-second excerpt")
    planned=[x for x in windows if x<=seconds]
    full=math.floor(seconds)
    if full>300 and full not in planned:
        planned.append(full)
    planned=sorted(set(planned))
    if 300 not in planned:
        raise ValueError("Frozen original 300-second reference not available")
    out=[]
    for win in planned:
        kept,dropped=trial_filter(events,win,*rt)
        target=sum(bool(x["is_target"]) for x in kept)
        standard=len(kept)-target
        out.append({"windowSeconds":win,"validTargetsBeforeEEGQC":target,
            "validStandardsBeforeEEGQC":standard,
            "correctRTTrialN":len(kept),"droppedByTaskRules":dropped})
    base=out[0]
    for row in out:
        row["extraValidTargetsVs300"]=row["validTargetsBeforeEEGQC"]-base["validTargetsBeforeEEGQC"]
        row["extraValidStandardsVs300"]=row["validStandardsBeforeEEGQC"]-base["validStandardsBeforeEEGQC"]
    return out


def analyze(protocol,outdir):
    p=json.loads(protocol.read_text())
    if p["excerpt"]["duration_seconds"]!=300 or p["source"]["sampling_rate_hz"]!=1024:
        raise ValueError("Original frozen duration or sampling rate has changed")
    if p["trial_rules"]["rt_min_seconds"]!=.2 or p["trial_rules"]["rt_max_seconds"]!=1.0:
        raise ValueError("Original RT rules unexpectedly changed")
    rows=[]
    for sid in SUBJECTS:
        src=next((z for z in p["source"]["files"] if z["subject"]==sid),None)
        if src is None: raise ValueError(f"Missing pinned source subject {sid}")
        downloads={}
        for kind in ("set","events"):
            info=src[kind]
            dest=outdir/info["name"]
            downloads[kind]=checked_download(info["url"],info["sha256"],dest)
        meta=loadmat(outdir/src["set"]["name"],squeeze_me=True,struct_as_record=False)["EEG"]
        channel_count=int(meta.nbchan)
        samples=int(meta.pnts)
        sr=int(meta.srate)
        if sr!=1024 or channel_count!=p["source"]["channels_total"]:
            raise ValueError(f"EEG metadata mismatch for subject {sid}")
        if samples*channel_count*4!=src["fdt"]["size"]:
            raise ValueError(f"FDT size inconsistent with source metadata for {sid}")
        events=parse_events_tsv(outdir/src["events"]["name"])
        if len(events)!=len(meta.event) if isinstance(meta.event,list) else False:
            # MATLAB event may be numpy ndarray: source validation in the
            # original full EEG loader remains authoritative.
            raise ValueError(f"Event count mismatch for subject {sid}")
        seconds=samples/sr
        windows=count_windows(events,seconds,rt=(
            p["trial_rules"]["rt_min_seconds"],p["trial_rules"]["rt_max_seconds"]))
        rows.append({"subject":sid,"fullDurationSeconds":seconds,
            "sourceSetSha256":downloads["set"],"sourceEventsSha256":downloads["events"],
            "sourceFdtSha256PinnedNotDownloaded":src["fdt"]["sha256"],
            "sourceEEGDownloaded":False,"trialSelectionUnchanged":True,
            "windows":windows,"hmmUnflaggedSupport":None,
            "eligiblePairedP300":None,"original20SubjectVerdict":"PARTIAL"})
    return {"study":"ERP CORE P3 exploratory later-event support, subject 4 and 8",
        "originalFrozenWindowSeconds":300,
        "noEEGDownloaded":True,"noHMMRelabeling":True,"originalEffectiveFloor":5,
        "newParticipantsEvaluated":False,"promotionAuthorized":False,
        "freshFullWindowScientificValidationPerformed":False,
        "rows":rows}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--protocol",type=Path,default=Path("benchmarks/erpcore_protocol.json"))
    p.add_argument("--work",type=Path,default=Path("download_event_meta"))
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    args.work.mkdir(parents=True,exist_ok=True)
    result=analyze(args.protocol,args.work)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2))
    for item in result["rows"]:
        print("PR51_EXTRA_EVENT_SUPPORT",json.dumps({
            "subject":item["subject"],"fullSeconds":item["fullDurationSeconds"],
            "windows":[{k:q[k] for k in ("windowSeconds","validTargetsBeforeEEGQC",
                      "validStandardsBeforeEEGQC","extraValidTargetsVs300",
                      "extraValidStandardsVs300")} for q in item["windows"]]}))
if __name__=="__main__":
    main()
