"""Recompute diagnostic facts from exported samples; render scientific figures."""
import csv,json,math,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=True)
h=json.loads((p/'hmm_audit.json').read_text());events=json.loads((p/'events.json').read_text());r=json.loads((p/'rejection_audit.json').read_text())
t=np.genfromtxt(p/'hmm_trace.csv',delimiter=',',names=True)
fig,axs=plt.subplots(3,2,figsize=(12,9),layout='constrained')
for i,q in enumerate(h['cases']):
 fs=q['sampleRate'];lat=events['latencySamples'][q['trial']-1];rel=t['timeSeconds']-(lat-1)/fs
 bl=np.arange(lat+round(-.2*fs),lat+1)-1;win=np.arange(lat+round(.3*fs),lat+round(.6*fs)+1)-1
 for label,ix in [('baseline',bl),('window',win)]:
  assert int(np.sum(t['hmmRaw'][ix]))==q[label+'RawHmmSamples']
  assert int(np.sum(t['hmmPadded'][ix]))==q[label+'PaddedHmmSamples']
 both=np.union1d(bl,win)
 assert bool(np.any(t['hmmRaw'][both]))==q['rawHmmOverlap']
 assert bool(np.any(t['controlHmmRaw'][both]))==q['controlRawOverlap']
 assert float(np.max(np.abs(t['eogUV'][both]-t['noInjectedEyeEogUV'][both])))==q['maxEogDifferenceInWindows']
 assert not q['dilationOnly'] and q['rawHmmOverlap'] and q['controlRawOverlap']
 use=(rel>=-.4)&(rel<=.8)
 for ax in axs[i]:
  ax.axvspan(-.2,0,color='#5d9cc9',alpha=.12,label='Baseline')
  ax.axvspan(.3,.6,color='#70b798',alpha=.15,label='P3 window')
  ax.set_xlim(-.4,.8);ax.set_xlabel('Time from event (s)');ax.grid(alpha=.18)
 axs[i,0].plot(rel[use],t['eogUV'][use],color='#272e3b',label='Original EOG')
 axs[i,0].plot(rel[use],t['noInjectedEyeEogUV'][use],color='#007f9e',ls='--',label='No injected eye signal')
 axs[i,0].set_ylabel('EOG (microvolts)');axs[i,0].set_title(f'Trial {q["trial"]}: observed background waveform')
 axs[i,1].plot(rel[use],t['hmmRobustZ'][use],color='#7c3a8c',label='Smoothed HMM feature')
 axs[i,1].axhline(3.5,color='#a34532',ls=':',label='Segment peak criterion (3.5)')
 # The criterion is for the maximum of a state segment, not a pointwise classifier.
 active=t['hmmRaw'][use].astype(bool);axs[i,1].fill_between(rel[use],0,t['hmmRobustZ'][use],where=active,color='#db9e40',alpha=.32,label='Raw HMM active')
 axs[i,1].set_ylabel('Robust z');axs[i,1].set_title('Active state overlaps without 50 ms padding')
 if i==0:
  axs[i,0].legend(fontsize=8);axs[i,1].legend(fontsize=8)
fig.suptitle('S005: three additional HMM flags — diagnostic, not natural-artifact ground truth',fontsize=13)
fig.savefig(out/'hmm_windows.png',dpi=180);fig.savefig(out/'hmm_windows.svg');plt.close(fig)
fig,axs=plt.subplots(2,1,figsize=(11,7),layout='constrained')
summary=[]
for i,q in enumerate(r):
 arm=q['arm'];uv=q['thresholdUV'];channel=q['crossingChannels'][0]
 curves={}
 for stage in ['before_eog','after_eog','before_rejection','without_eog_same_other_decisions']:
  rr=list(csv.DictReader((p/f'trial23_{arm}_{stage}.csv').open()))
  expected=next(z for z in q['stagePeaks'] if z['stage']==stage)
  maxima={ch:max(abs(float(z['valueUV'])) for z in rr if z['channel']==ch) for ch in expected['channels']}
  for ch,mx in zip(expected['channels'],expected['perChannelMaxAbsUV']):assert math.isclose(maxima[ch],mx,abs_tol=1e-9)
  zz=[z for z in rr if z['channel']==channel];curves[stage]=(np.array([float(z['timeSeconds']) for z in zz]),np.array([float(z['valueUV']) for z in zz]))
  if stage=='before_rejection':
   assert {ch for ch,mx in maxima.items() if mx>uv}==set(q['crossingChannels'])
   for ch in q['channelPeaks']:
    zzch=[z for z in rr if z['channel']==ch['channel']];peak=max(zzch,key=lambda z:abs(float(z['valueUV'])))
    assert math.isclose(float(peak['valueUV']),ch['signedUV'],abs_tol=1e-9)
    assert math.isclose(float(peak['timeSeconds']),ch['timeSeconds'],abs_tol=1e-12)
 ax=axs[i];ax.axvspan(-.2,0,color='#5d9cc9',alpha=.12,label='Baseline');ax.axvspan(.3,.6,color='#70b798',alpha=.15,label='P3 window')
 for stage,color,label in [('before_eog','#007f9e','Before EOG regression'),('after_eog','#a34532','After EOG regression')]:ax.plot(*curves[stage],color=color,label=label)
 for sign in [-1,1]:ax.axhline(sign*uv,color='#555',ls='--',lw=1,label=f'Frozen threshold ±{uv:g} µV' if sign==1 else None)
 peak=next(x for x in q['channelPeaks'] if x['channel']==channel);ax.plot(peak['timeSeconds'],peak['signedUV'],'o',color='#a34532')
 name='Single EOG, 0.1–20 Hz' if arm=='baseline' else 'Expanded EOG, 0.1–40 Hz'
 ax.set_title(f'{name}: {channel}, peak {peak["maxAbsUV"]:.3f} µV at {1000*peak["timeSeconds"]:.2f} ms')
 ax.set_xlim(-.2,.8);ax.set_ylabel('Baseline-corrected amplitude (µV)');ax.set_xlabel('Time from event (s)');ax.grid(alpha=.18);ax.legend(fontsize=8,loc='lower left')
 summary.append({'arm':arm,'channel':channel,'thresholdUV':uv,'peakUV':peak['maxAbsUV'],'timeMs':1000*peak['timeSeconds'],'beforeEogPeakUV':q['stagePeaks'][0]['maxAbsUV'],'replayMaxError':q['replayMaxError']})
 assert q['replayMaxError']==0 and q['rejectionExactlyReproduced']
fig.suptitle('Trial 23: both crossings already exist before EOG regression at the final frozen thresholds',fontsize=12)
fig.savefig(out/'trial23_rejection.png',dpi=180);fig.savefig(out/'trial23_rejection.svg');plt.close(fig)
result={'hmmCases':3,'selectedRecipes':2,'sourceSampleVerification':'passed','scope':'seed-zero diagnostic; fixed decisions, not new adaptive search','rejection':summary}
(out/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print('WAVEFORM_SAMPLE_VERIFICATION_PASSED '+json.dumps(result))
