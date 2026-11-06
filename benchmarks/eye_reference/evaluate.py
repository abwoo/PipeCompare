import pathlib,json
import numpy as np
from scipy.optimize import linear_sum_assignment
P=pathlib.Path('eye_reference_results');meta=json.loads((P/'reference.json').read_text());fs=meta['fs']
M=np.genfromtxt(P/'masks.csv',delimiter=',',names=True);R=np.genfromtxt(P/'optical_reference.csv',delimiter=',',names=True)
assert len(M)==len(R)==meta['samples'] and np.max(abs(M['seconds']-R['seconds']))<1e-8
eligible=R['referenceEligible']>0;background=eligible & ~(R['binocularBlink']>0);background_seconds=float(background.sum()/fs);t=M['seconds'];refs=np.array(meta['binocularIntervals']);out={}
for field in ['hmmRaw','hmm','robust','union','fp1Hmm']:
 mask=M[field]>0;d=np.diff(np.r_[0,mask,0].astype(int));starts=np.where(d==1)[0];ends=np.where(d==-1)[0]-1
 pred=np.c_[t[starts],t[ends]];valid_fraction=np.array([eligible[a:b+1].mean() for a,b in zip(starts,ends)])
 # Maximum cardinality one-to-one event match at frozen 100ms tolerance.
 if len(pred):
  overlaps=(pred[:,None,0]<=refs[None,:,1]+.1)&(pred[:,None,1]>=refs[None,:,0]-.1)
  pi,ri=linear_sum_assignment(-overlaps.astype(int));pairs=[(int(p),int(r)) for p,r in zip(pi,ri) if overlaps[p,r]]
 else:pairs=[]
 matched={p for p,r in pairs};unmatched=[i for i in range(len(pred)) if i not in matched]
 comparable=[i for i in unmatched if valid_fraction[i]>=.8];indeterminate=[i for i in unmatched if valid_fraction[i]<.8]
 saccades=np.array(meta['opticalSaccadeIntervals'])
 saccade_overlap=[i for i in comparable if len(saccades) and np.any((pred[i,0]<=saccades[:,1])&(pred[i,1]>=saccades[:,0]))]
 all_saccade_overlap=[i for i in unmatched if len(saccades) and np.any((pred[i,0]<=saccades[:,1])&(pred[i,1]>=saccades[:,0]))]
 rec={'unmatchedOverlappingSaccade':len(all_saccade_overlap),'unmatchedEligibleOverlappingSaccade':len(saccade_overlap),'unmatchedEligibleWithoutSaccade':len(comparable)-len(saccade_overlap),'referenceEvents':len(refs),'predictedEvents':len(pred),'matchedEvents':len(pairs),'missedReferenceEvents':len(refs)-len(pairs),'referenceRecall':len(pairs)/len(refs),'unmatchedReferenceEligible':len(comparable),'unmatchedIndeterminate':len(indeterminate),'eligibleSeconds':float(eligible.sum()/fs),'eligibleBackgroundSeconds':background_seconds,'unmatchedPerEligibleMinute':len(comparable)/(eligible.sum()/fs/60) if background_seconds>0 else None,'markedFraction':float(mask.mean()),'segmentsOver2Seconds':int(sum((pred[:,1]-pred[:,0])>2)) if len(pred) else 0,'pairs':pairs,'predictedIntervals':pred.tolist(),'referenceEligibleFractions':valid_fraction.tolist(),'unmatchedEligibleIndices':comparable,'unmatchedIndeterminateIndices':indeterminate}
 assert rec['matchedEvents']+len(unmatched)==len(pred)
 out[field]=rec
summary={'scope':'One participant, independent optical algorithm reference; not final natural-blink truth or P3 accuracy','maxSyncResidualMs':meta['maxInteriorSyncResidualMs'],'durationSeconds':len(t)/fs,'results':out}
(P/'summary.json').write_text(json.dumps(summary,indent=2))
print('EYE_REFERENCE_RESULTS',json.dumps(summary));print('EYE_REFERENCE_EVALUATION_PASSED')
