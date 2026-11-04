import json,sys,hashlib
from pathlib import Path
from collect_focus import read,close,outcome,stage_check,STAGES
ROOT=Path(__file__).parent

def pair(folder):
 b=read(folder/'baseline/result.json');v=read(folder/'variant/result.json');s=b['subject'];seed=b['seed']
 assert s in [4,5] and seed in range(5)
 ref=next(x for x in read(ROOT/'eog_reference.json') if x['subject']==s and x['seed']==seed)
 for key in ['subject','seed','source','populationTruth','realizedTruth','blinkPrevalence','searchOptions','eeglabSha','executionCommit','workflowRun']:assert b[key]==v[key],key
 for key in ['source','populationTruth','realizedTruth','searchOptions','eeglabSha']:assert b[key]==ref['within'][key],key
 for key in ['mse','rawMse']:close(b['erp'][key],ref['within']['erp'][key])
 assert b['erp']['pipeline']==ref['within']['erp']['pipeline']
 traces={}
 for arm,d in [('baseline',b),('variant',v)]:
  assert d['completed'] and not d['errors'] and d['arm']==arm
  outcome(d['erp'],d['populationTruth']);close(d['replayMaxError'],0,1e-6)
  traces[arm]={}
  for stage in STAGES:
   t=read(folder/arm/(stage+'.json'));pick=stage_check(t,d,s,seed,stage,arm)
   traces[arm][stage]=t
   if stage=='after_eog':close(pick['readout']['mse'],d['erp']['mse']);assert pick['pipeline']==d['erp']['pipeline']
 # Instrumentation and expanded branch cannot change pre-EOG outputs.
 for stage in ['initial','after_cutoffs']:
  bb=traces['baseline'][stage];vv=traces['variant'][stage]
  assert bb['recommended']==vv['recommended']
  assert bb['rows']==vv['rows'],'Pre-EOG change outside authorized variant'
 return {'subject':s,'seed':seed,'baseline':b,'variant':v,'original':ref['original'],'traces':traces}

def summarize(records):
 assert len(records)==10 and {(r['subject'],r['seed']) for r in records}=={(s,k) for s in [4,5] for k in range(5)}
 assert len({r['baseline']['executionCommit'] for r in records})==1
 rows=[]
 for r in records:
  b=r['baseline']['erp'];v=r['variant']['erp'];o=r['original']
  rows.append({'subject':r['subject'],'seed':r['seed'],'baselineMSE':b['mse'],'variantMSE':v['mse'],'difference':v['mse']-b['mse'],'baselineRawMSE':b['rawMse'],'variantRawMSE':v['rawMse'],'originalMSE':o['mse'],'originalRawMSE':o['rawMse']})
 nonworse=all(x['difference']<=1e-7 and x['variantRawMSE']<=x['baselineRawMSE']+1e-7 for x in rows)
 improve=any(x['difference'] < -1e-7 for x in rows)
 recovered=all(x['variantMSE']<=x['originalMSE']+1e-7 and x['variantRawMSE']<=x['originalRawMSE']+1e-7 for x in rows if x['subject']==5)
 return {'pairs':10,'participants':2,'independentConfirmation':False,'rows':rows,'noStartRegressions':nonworse,'someImprovement':improve,'S005Recovered':recovered,'decision':'proceed-to-ten-person-guardrail' if nonworse and improve and recovered else 'reject-expansion','productionActivated':False}

if __name__=='__main__':
 for f,h in read(ROOT/'eog_plan.json')['hashes'].items():assert hashlib.sha256((ROOT/f).read_bytes()).hexdigest()==h
 if sys.argv[1]=='--pair':
  r=pair(Path(sys.argv[2]));print('EOG_PAIR_VALIDATED',r['subject'],r['seed'])
 else:
  out=Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=True)
  try:
   rr=[pair(p) for p in sorted(Path(sys.argv[1]).glob('eog-pair-*'))];s=summarize(rr)
   (out/'records.json').write_text(json.dumps(rr,separators=(',',':'))+'\n');(out/'summary.json').write_text(json.dumps(s,indent=2)+'\n')
   print('PIPECOMPARE_EOG_SUMMARY '+json.dumps(s))
  except Exception as ex:
   (out/'failure.txt').write_text(repr(ex));raise
