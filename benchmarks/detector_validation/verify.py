import csv,json,pathlib,math
p=pathlib.Path('detector_results');rows=json.loads((p/'external_summary.json').read_text());controls=json.loads((p/'known_truth_summary.json').read_text())
assert len(rows)==40 and len(controls)==40
assert {(r['subject'],r['block'],r['channels']) for r in rows}=={(s,b,c) for s in range(1,11) for b in ['open','blink'] for c in [1,2]}
for r in rows:
    f=p/f"S{r['subject']:03d}_{r['block']}_ch{r['channels']}.csv"
    mask=list(csv.DictReader(f.open()));assert len(mask)==60000
    h=[int(x['hmm']) for x in mask];hp=[int(x['hmmPadded']) for x in mask];rp=[int(x['robustPadded']) for x in mask]
    values={'hmmFraction':sum(h)/len(h),'hmmPaddedFraction':sum(hp)/len(h),'robustPaddedFraction':sum(rp)/len(h),'unionFraction':sum(a or b for a,b in zip(hp,rp))/len(h),'hmmOnlyFraction':sum(a and not b for a,b in zip(hp,rp))/len(h)}
    for key,value in values.items():assert abs(r[key]-value)<1e-12
    assert r['hmmSegments']==sum(v and not prev for prev,v in zip([0]+h[:-1],h))
for r in controls:
    assert r['totalEvents']==(0 if r['amplitude']==0 else 10)
    assert 0<=r['hitEvents']<=r['totalEvents']
    for key in ['outsideSupportFraction','markedFraction']:assert math.isfinite(r[key]) and 0<=r[key]<=1
summary={'executionChecks':'passed','externalInvalidFits':sum(not r['hmmValid'] for r in rows),'syntheticInvalidFits':sum(not r['hmmValid'] for r in controls),'external':{},'synthetic':{}}
for c in [1,2]:
    summary['external'][str(c)]={b:{k:sum(r[k] for r in rows if r['channels']==c and r['block']==b)/10 for k in ['hmmPaddedFraction','robustPaddedFraction','hmmOnlyFraction']} for b in ['open','blink']}
for amp in [0,4,8,16]:
    q=[r for r in controls if r['amplitude']==amp];assert len(q)==10
    summary['synthetic'][str(amp)]={'hits':sum(r['hitEvents'] for r in q),'events':sum(r['totalEvents'] for r in q),'outsideSupportFraction':sum(r['outsideSupportFraction'] for r in q)/10}
(p/'verified_summary.json').write_text(json.dumps(summary,indent=2));print('DETECTOR_INDEPENDENT_CHECKS_PASSED',json.dumps(summary))
