"""Fetch one prespecified participant with a pinned official source manifest."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.io import savemat
from fetch_open_data import fetch, read_edf

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--stage',choices=('development','confirmation'),required=True)
    parser.add_argument('--subject',type=int,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    protocol=json.loads(Path(__file__).with_name('optimization_protocol.json').read_text())
    people=protocol[args.stage+'_subjects'];assert args.subject in people
    args.out.mkdir(parents=True,exist_ok=True)
    manifest=args.out/'SHA256SUMS.txt'
    fetch('https://physionet.org/files/eegmmidb/1.0.0/SHA256SUMS.txt',manifest,protocol['source_manifest_sha256'])
    checks={}
    for line in manifest.read_text().splitlines():
        if not line.strip():
            continue
        digest,path=line.split(maxsplit=1)
        checks[path.lstrip('*./')]=digest
    entries=[]
    for subject in people:
        path=f'S{subject:03d}/S{subject:03d}R01.edf'
        entries.append(dict(subject=subject,path=path,sha256=checks[path]))
    frozen=dict(stage=args.stage,protocol=protocol,files=entries)
    # All input identities and hashes are written before fitting any data.
    (args.out/'frozen_inputs.json').write_text(json.dumps(frozen,indent=2)+'\n')
    entry=next(e for e in entries if e['subject']==args.subject)
    path=args.out/Path(entry['path']).name
    fetch('https://physionet-open.s3.amazonaws.com/eegmmidb/1.0.0/'+entry['path'],path,entry['sha256'])
    data,fs,labels=read_edf(path)
    savemat(path.with_suffix('.mat'),dict(data=data,fs=fs,labels=np.asarray(labels,dtype=object)),do_compression=True)
    print(f'FROZEN_OPTIMIZATION_INPUT stage={args.stage} participant={args.subject} sha256={entry["sha256"]}')

if __name__=='__main__':
    main()
