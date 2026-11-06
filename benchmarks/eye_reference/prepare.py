import pathlib, urllib.request, zipfile, hashlib, re, json
import numpy as np
from scipy.io import loadmat,savemat
P=pathlib.Path('eye_reference_results');P.mkdir(exist_ok=True)
url='https://www.eyetracking-eeg.org/testdata/reading.zip'
b=urllib.request.urlopen(url,timeout=90).read();digest=hashlib.sha256(b).hexdigest()
assert digest=='54dd72106c5d019f7bd95d18a3a613a18f5be72f929f9aa87da5036e90544651'
(P/'source.zip').write_bytes(b)
with zipfile.ZipFile(P/'source.zip') as z:
 for name in ['reading_eeg.set','reading_eeg.fdt','reading_eyelink.asc']:(P/name).write_bytes(z.read(name))
E=loadmat(P/'reading_eeg.set',simplify_cells=True)['EEG'];fs=float(E['srate']);n=int(E['pnts']);nc=int(E['nbchan'])
x=np.fromfile(P/'reading_eeg.fdt',dtype='<f4');assert len(x)==n*nc
x=x.reshape((nc,n),order='F').astype(float);assert np.isfinite(x).all()
labels=[c['labels'] for c in E['chanlocs']];txt=(P/'reading_eyelink.asc').read_text()
a=np.array([int(v) for v in re.findall(r'MSG\s+(\d+)\s+MYKEYWORD\s+3(?:\s+\(\d+\))?\s*$',txt,re.M)])
b=np.array([(e['latency']-1)/fs for e in E['event'] if str(e['type'])=='3'])
assert len(a)==len(b)>=5 and np.all(np.diff(a)>0) and np.all(np.diff(b)>0)
slope=(b[-1]-b[0])/(a[-1]-a[0]);convert=lambda t:b[0]+(np.asarray(t)-a[0])*slope
residual=convert(a)-b;assert np.max(np.abs(residual[1:-1]))<=.010,'Independent sync residual >10ms'
assert not any(str(e['type'])=='boundary' and b[0]<(e['latency']-1)/fs<b[-1] for e in E['event'])
start=int(round(b[0]*fs));end=int(round(b[-1]*fs));t=np.arange(start,end+1)/fs
# Binocular overlapping optical blink events. Unilateral events remain unknown.
blinks={eye:[] for eye in ['L','R']}
for eye,on,off in re.findall(r'^EBLINK\s+([LR])\s+(\d+)\s+(\d+)\s+\d+',txt,re.M):
 on,off=convert([int(on),int(off)])
 if on>=t[0] and off<=t[-1]:blinks[eye].append([float(on),float(off)])
refs=[]
for l in blinks['L']:
 matches=[r for r in blinks['R'] if max(l[0],r[0])<=min(l[1],r[1])]
 if matches:refs.append([min([l[0]]+[r[0] for r in matches]),max([l[1]]+[r[1] for r in matches])])
merged=[]
for lo,hi in sorted(refs):
 if merged and lo<=merged[-1][1]:merged[-1][1]=max(hi,merged[-1][1])
 else:merged.append([lo,hi])
refs=merged;assert refs,'No binocular reference blinks'
# Optical pupil availability; unknown samples do not become certified non-blinks.
samples=[]
for line in txt.splitlines():
 if re.match(r'^\d+\s',line):
  fields=line.split()
  if len(fields)>=7:
   try:samples.append([int(fields[0]),float(fields[3]),float(fields[6])])
   except ValueError:pass
S=np.array(samples);st=convert(S[:,0]);ix=np.searchsorted(st,t);ix=np.clip(ix,0,len(st)-1)
valid=(np.abs(st[ix]-t)<=.003)&(S[ix,1]>0)&(S[ix,2]>0)
reference=np.zeros(len(t),bool);unknown=np.zeros(len(t),bool)
for lo,hi in refs:reference|=(t>=lo)&(t<=hi)
for intervals in blinks.values():
 for lo,hi in intervals:unknown|=(t>=lo-.1)&(t<=hi+.1)
valid=valid&~unknown;valid|=reference
X=np.stack([x[labels.index('IO1')]-x[labels.index('Fp1')],x[labels.index('IO2')]-x[labels.index('Fp2')]])[:,start:end+1]
F=x[labels.index('Fp1'),start:end+1][None,:]
savemat(P/'input.mat',{'X':X,'F':F,'fs':fs,'time':t})
np.savetxt(P/'optical_reference.csv',np.c_[t,valid,reference],delimiter=',',header='seconds,referenceEligible,binocularBlink',comments='')
meta={'source':url,'zipSHA256':digest,'fs':fs,'startSeconds':float(t[0]),'endSeconds':float(t[-1]),'samples':len(t),'syncCount':len(a),'maxInteriorSyncResidualMs':float(max(abs(residual[1:-1]))*1000),'syncResidualMs':(residual*1000).tolist(),'binocularIntervals':refs,'leftIntervals':blinks['L'],'rightIntervals':blinks['R'],'note':'One participant; optical algorithm reference, not adjudicated video truth; no event timing fitted to EEG amplitude.'}
(P/'reference.json').write_text(json.dumps(meta,indent=2));print('REFERENCE_READY',json.dumps(meta))
