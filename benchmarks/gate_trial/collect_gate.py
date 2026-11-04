import json,math,sys,statistics as st
from pathlib import Path
p=Path(sys.argv[1]);r=json.loads((p/'records.json').read_text())
assert len(r)==10 and {(x['subject'],x['seed']) for x in r}=={(s,k) for s in [4,5] for k in range(5)}
reference={x['subject']:x for x in json.loads((Path(__file__).resolve().parents[1]/'results/within_trial_20261104/records.json').read_text())}
for x in r:
 for arm in ['baseline','expanded','gated']:
  a=x[arm];t=x['populationTruth'];assert a['available']
  assert abs(st.mean((v/a['gain']-z)**2 for v,z in zip(a['estimate'],t))-a['mse'])<1e-7
  assert abs(st.mean((v-z)**2 for v,z in zip(a['estimate'],t))-a['rawMse'])<1e-7
 assert x['gated']['pipeline']==(x['expanded']['pipeline'] if x['gate']['switched'] else x['baseline']['pipeline'])
nonworse=all(x['gated'][k]<=x['baseline'][k]+1e-7 for x in r for k in ['mse','rawMse'])
improved=any(x['gated']['mse']<x['baseline']['mse']-1e-7 for x in r)
recovered=all(x['gated'][k]<=reference[5]['baseline']['erp'][k]+1e-7 for x in r if x['subject']==5 for k in ['mse','rawMse'])
s={'historicalPairs':10,'participants':2,'independentConfirmation':False,'switches':sum(x['gate']['switched'] for x in r),'noStartRegression':nonworse,'someImprovement':improved,'S005Recovered':recovered,'decision':'proceed-to-frozen-unseen-starts-5-to-9' if nonworse and improved and recovered else 'reject-gate','rows':[{'subject':x['subject'],'seed':x['seed'],'baselineMSE':x['baseline']['mse'],'expandedMSE':x['expanded']['mse'],'gatedMSE':x['gated']['mse'],'switched':x['gate']['switched']} for x in r]}
(p/'summary.json').write_text(json.dumps(s,indent=2)+'\n');print('REPLACEMENT_GATE_SUMMARY '+json.dumps(s))
