"""Fetch exactly the frozen 10 x 60-second public pilot; no subject data is committed."""
import argparse
import hashlib
import json
import urllib.request
from pathlib import Path
import numpy as np
from scipy.io import savemat

PROTOCOL = Path(__file__).with_name("public_pilot_protocol.json")

def read_edf(path):
    with Path(path).open("rb") as f:
        h=f.read(256)
        if len(h)!=256: raise ValueError("Truncated EDF header")
        nbytes,nrec=int(h[184:192]),int(h[236:244])
        duration,ns=float(h[244:252]),int(h[252:256])
        fields=f.read(nbytes-256); offset=0
        def field(w):
            nonlocal offset
            vals=[fields[offset+i*w:offset+(i+1)*w].decode("ascii").strip() for i in range(ns)]
            offset+=w*ns
            return vals
        labels=field(16);field(80);units=field(8)
        pmin,pmax=np.array(field(8),float),np.array(field(8),float)
        dmin,dmax=np.array(field(8),float),np.array(field(8),float)
        field(80);sizes=np.array(field(8),int);field(32)
        raw=np.frombuffer(f.read(),dtype="<i2")
    if offset!=len(fields) or nrec<=0 or duration<=0 or raw.size!=nrec*sum(sizes):
        raise ValueError("Incomplete or inconsistent EDF")
    raw=raw.reshape(nrec,sum(sizes)); starts=np.r_[0,np.cumsum(sizes)]
    keep=[i for i,s in enumerate(labels) if "Annotations" not in s]
    if len(keep)!=64 or len(set(sizes[keep]))!=1: raise ValueError("Unexpected source layout")
    rows=[]
    for i in keep:
        if units[i].lower() not in ("uv","µv") or dmax[i]<=dmin[i]: raise ValueError("Invalid EDF calibration")
        a=raw[:,starts[i]:starts[i+1]].reshape(-1).astype(float)
        rows.append((a-dmin[i])*(pmax[i]-pmin[i])/(dmax[i]-dmin[i])+pmin[i])
    return np.stack(rows),float(sizes[keep[0]]/duration),[labels[i].rstrip(".") for i in keep]

def prepare(path, out, subject, p):
    entry=next(f for f in p["source"]["files"] if f["subject"]==subject)
    digest=hashlib.sha256(Path(path).read_bytes()).hexdigest()
    if digest!=entry["sha256"]: raise ValueError("EDF SHA256 mismatch")
    x,fs,labels=read_edf(path)
    if fs!=p["excerpt"]["sampling_rate"] or not np.isfinite(x).all(): raise ValueError("Invalid source samples/rate")
    start=round(p["excerpt"]["offset_seconds"]*fs); n=round(p["excerpt"]["duration_seconds"]*fs)
    if x.shape[1]<start+n: raise ValueError("Source shorter than prespecified excerpt")
    lower=[s.lower() for s in labels]
    wanted=p["excerpt"]["channels"]
    loc=[lower.index(s.lower()) for s in wanted]
    data=x[loc,start:start+n]
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    savemat(out/"input.mat",dict(data=data,fs=fs,labels=np.array(wanted,dtype=object)),do_compression=True)
    evidence=dict(subject=subject,source_url=p["source"]["base_url"]+entry["path"],source_sha256=digest,
                  samples=n,sampling_rate=fs,channels=wanted,offset_seconds=p["excerpt"]["offset_seconds"],
                  duration_seconds=p["excerpt"]["duration_seconds"],protocol_sha256=hashlib.sha256(PROTOCOL.read_bytes()).hexdigest())
    (out/"source.json").write_text(json.dumps(evidence,indent=2)+"\n")
    return evidence

def main():
    a=argparse.ArgumentParser();a.add_argument("--subject",type=int,required=True);a.add_argument("--out",type=Path,required=True)
    args=a.parse_args();p=json.loads(PROTOCOL.read_text())
    entry=next(f for f in p["source"]["files"] if f["subject"]==args.subject)
    args.out.mkdir(parents=True,exist_ok=True);path=args.out/"source.edf"
    url=p["source"]["base_url"]+entry["path"]
    with urllib.request.urlopen(urllib.request.Request(url,headers={"User-Agent":"PipeCompare-public-pilot/1"}),timeout=90) as r:
        path.write_bytes(r.read())
    print(json.dumps(prepare(path,args.out,args.subject,p)),flush=True)
if __name__=="__main__":main()
