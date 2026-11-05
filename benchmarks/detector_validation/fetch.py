"""Fetch only complete EDF records covering two fixed short windows per subject."""
import csv, io, json, hashlib, math, urllib.request, pathlib, time
import numpy as np
from scipy.io import savemat
ROOT='https://s3.amazonaws.com/openneuro.org/ds006386/'
OUT=pathlib.Path('detector_results'); OUT.mkdir(exist_ok=True)
def get(path, span=None):
    for attempt in range(3):
        try:
            req=urllib.request.Request(ROOT+path,headers={'Range':f'bytes={span[0]}-{span[1]}'} if span else {})
            with urllib.request.urlopen(req,timeout=90) as r:
                if span: assert r.status==206, 'Server did not honor range'
                b=r.read()
            if span: assert len(b)==span[1]-span[0]+1
            return b
        except Exception:
            if attempt==2: raise
            time.sleep(2)
def digest(b):return hashlib.sha256(b).hexdigest()
for sub in range(1,11):
    p=f'sub-{sub}/eeg/sub-{sub}_task-artifact_run-01_'
    annotations=get(f'derivatives/Manual_Annotations/sub{sub}_run01.csv')
    events=get(p+'events.tsv')
    ev=list(csv.DictReader(io.StringIO(events.decode('utf-8-sig')),delimiter='\t'))
    # First instructed blink block; first eye-open baseline follows eye-closed.
    blink=float(next(r['onset'] for r in ev if r['trial_type']=='4'))
    ann=list(csv.DictReader(io.StringIO(annotations.decode('utf-8-sig'))))
    opened=min(float(r['start_time']) for r in ann if r['label']=='open_base')+2
    url=p+'eeg.edf'; h0=get(url,(0,255)); hb=int(h0[184:192]); h=get(url,(0,hb-1))
    ns=int(h[252:256]); dur=float(h[244:252]); nr=int(h[236:244]); pos=256
    fields=[]
    for width in [16,80,8,8,8,8,8,80,8,32]:
        fields.append([h[pos+i*width:pos+(i+1)*width].decode().strip() for i in range(ns)]);pos+=ns*width
    labels,_,units,pmin,pmax,dmin,dmax,_,samples,_=fields
    samples=np.array(samples,int); rb=int(samples.sum()*2); starts=np.r_[0,np.cumsum(samples)]
    ch=[labels.index('Fp1'),labels.index('Fp2')]; fs=samples[ch[0]]/dur
    assert fs==1000 and all(samples[c]/dur==fs for c in ch)
    provenance={'subject':sub,'source':ROOT+url,'headerSHA256':digest(h),'annotationsSHA256':digest(annotations),'eventsSHA256':digest(events),'windows':[]}
    (OUT/f'S{sub:03d}_annotations.csv').write_bytes(annotations)
    for name,start in [('open',opened),('blink',blink)]:
        a=math.floor(start/dur); z=math.ceil((start+60)/dur); assert z<=nr
        span=(hb+a*rb,hb+z*rb-1); b=get(url,span)
        data=np.frombuffer(b,dtype='<i2').reshape(z-a,-1); signals=[]
        for c in ch:
            x=data[:,starts[c]:starts[c+1]].reshape(-1).astype(float)
            x=(x-float(dmin[c]))*(float(pmax[c])-float(pmin[c]))/(float(dmax[c])-float(dmin[c]))+float(pmin[c])
            scale={'uV':1,'µV':1,'V':1e6}[units[c]]; signals.append(x*scale)
        offset=round((start-a*dur)*fs); x=np.array(signals)[:,offset:offset+int(60*fs)]
        assert x.shape==(2,60000) and np.isfinite(x).all()
        # Task interval overlap is descriptive only, not sample-level blink truth.
        relevant=[r for r in ann if float(r['stop_time'])>start and float(r['start_time'])<start+60]
        savemat(OUT/f'S{sub:03d}_{name}.mat',{'X':x,'fs':fs,'start':start})
        provenance['windows'].append({'name':name,'start':start,'duration':60,'range':span,'rangeSHA256':digest(b),'units': [units[c] for c in ch],'overlappingAnnotations':relevant})
    (OUT/f'S{sub:03d}_provenance.json').write_text(json.dumps(provenance,indent=2))
    print('INPUT_READY',sub,flush=True)
