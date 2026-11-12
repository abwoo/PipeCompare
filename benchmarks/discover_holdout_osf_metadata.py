#!/usr/bin/env python3
"""Read *public OSF file metadata only* for the predeclared P3 holdout.

Never accesses trial-level EEG or generates scientific results. This is a
source identity discovery step; a human-reviewed SHA-pinned manifest must be
committed before any held-out EEG download or scientific evaluation.
"""
import argparse
from collections import deque
import json
from pathlib import Path
import re
from urllib.parse import urlparse
from urllib.request import Request, urlopen


def endpoint(uri):
    host = urlparse(uri).hostname
    if host != "api.osf.io" or urlparse(uri).scheme != "https":
        raise ValueError(f"untrusted OSF API metadata endpoint: {host}")
    return uri


def fetch_json(url):
    req = Request(endpoint(url), headers={"Accept":"application/vnd.api+json",
                                          "User-Agent":"PipeCompare-Holdout-Metadata/1.0"})
    with urlopen(req, timeout=30) as response:
        return json.load(response)


def related_files(entry):
    rel=(entry.get("relationships") or {}).get("files") or {}
    val=(rel.get("links") or {}).get("related")
    if isinstance(val,dict):
        return val.get("href")
    if isinstance(val,str):
        return val
    return None


def discover(component, subjects, limit=5000):
    root=f"https://api.osf.io/v2/nodes/{component}/files/osfstorage/"
    pending=deque([root]);visited=set();files={};discovered=0
    subjects=set(subjects)
    expected={
        f"sub-{s:03d}_task-P3_{part}":(s,part)
        for s in subjects for part in ("eeg.set","eeg.fdt","events.tsv")
    }
    while pending:
        url=pending.popleft()
        if url in visited:continue
        visited.add(url)
        if len(visited)>limit:
            raise RuntimeError("OSF tree traversal exceeded safe metadata page budget")
        response=fetch_json(url)
        items=response.get("data")
        if isinstance(items,dict):items=[items]
        if not isinstance(items,list):
            raise RuntimeError("Unexpected OSF response lacks file metadata list")
        for record in items:
            kind=record.get("attributes",{}).get("kind")
            name=record.get("attributes",{}).get("name")
            if kind=="folder":
                child=related_files(record)
                if child:pending.append(endpoint(child))
            elif kind=="file":
                discovered+=1
                if name not in expected:continue
                if name in files:raise ValueError(f"duplicate heldout P3 file name {name}")
                attrs=record["attributes"]
                hashes=(attrs.get("extra") or {}).get("hashes") or {}
                sha=str(hashes.get("sha256") or "")
                if not re.fullmatch(r"[a-f0-9]{64}",sha):
                    raise ValueError(f"missing expected SHA256 in OSF metadata {name}")
                download=(record.get("links") or {}).get("download")
                if not download or urlparse(download).scheme!="https" or (
                        urlparse(download).hostname not in {"osf.io","files.osf.io"}):
                    raise ValueError(f"unknown public OSF download link {name}")
                s,part=expected[name]
                files[name]={"subject":s,"type":part,
                    "name":name,"sha256":sha,"md5":hashes.get("md5"),
                    "size":attrs.get("size"),"osf_id":record.get("id"),
                    "download":download}
        nxt=(response.get("links") or {}).get("next")
        if nxt:pending.append(endpoint(nxt))
    missing=sorted(set(expected)-set(files))
    if missing:
        raise ValueError(f"metadata-only crawl incomplete: missing {len(missing)}/{len(expected)} files; examples: {missing[:8]}; pages={len(visited)} files_scanned={discovered}")
    return {"metadata_only":True,"eeg_bytes_downloaded":0,
            "osf_node":component,"subject_ids":sorted(subjects),
            "files":sorted(files.values(),key=lambda v:(v["subject"],v["type"])),
            "pages_visited":len(visited),"files_scanned":discovered,
            "status":"CANDIDATE_SOURCE_MANIFEST_NOT_YET_FROZEN"}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--holdout",type=Path,required=True)
    p.add_argument("--original",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    plan=json.loads(args.holdout.read_text())
    original=json.loads(args.original.read_text())
    if plan["status"]!="PRE-FREEZE-DESIGN-NOT_EXECUTED":
        raise ValueError("No held-out source discovery after opening outcome data")
    subjects=plan["previously_unseen_14"]
    if set(subjects)&{x["subject"] for x in original["source"]["files"]}:
        raise ValueError("Development EEG exposure contamination")
    result=discover(original["source"]["osf_p3_component"],subjects)
    if len(result["files"])!=42:
        raise ValueError("Need exactly 3 files for each of 14 unexposed participants")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2))
    print("PR51_OSF_METADATA_DISCOVERY",json.dumps({
        "subjects":subjects,"expected_files":42,
        "found":len(result["files"]),"pages":result["pages_visited"],
        "no_eeg_download":True}))


if __name__=="__main__":
    main()
