"""Independently reconcile every retained-trial count with the identity CSVs."""
import csv,json,sys
from pathlib import Path
p=Path(sys.argv[1]);refs=json.loads((p/'reference.json').read_text());data=json.loads((p/'retention.json').read_text())
flags={s:[{k:int(v) for k,v in r.items()} for r in csv.DictReader((p/f'S{s:03d}_trial_flags.csv').open())] for s in [4,5]}
def seq(x):return x if isinstance(x,list) else [x]
for q in refs:
 rr=flags[q['subject']];assert len(rr)==q['trials']
 assert sum(x['injectedBlink'] for x in rr)==q['injectedBlinks']
 assert sum(x['combinedFlag'] for x in rr)==q['flagged']
 for x in rr:
  assert x['combinedFlag']==int(any(x[k] for k in ['finalBlink','finalHorizontal','finalMuscle']))
  assert x['finalBlink']==int(any(x[k] for k in ['epochBlink','sampleRobustBlink','sampleHmmBlink']))
  assert x['finalHorizontal']==int(any(x[k] for k in ['epochHorizontal','sampleRobustHorizontal','sampleHmmHorizontal']))
  assert x['epochMuscle']==x['finalMuscle']
 for c in q['conditions']:
  zz=[x for x in rr if x['condition']==c['condition']]
  assert len(zz)==c['referenceTrials']
  assert sum(not x['combinedFlag'] for x in zz)==c['unflagged']
  assert sum(x['injectedBlink'] for x in zz)==c['injectedBlinks']
selected=[]
for d in data:
 rr=flags[d['subject']];absent=set(seq(d['absentIds']))
 for c in d['conditions']:
  zz=[x for x in rr if x['condition']==c['condition']]
  kept=[x for x in zz if x['trialId'] not in absent];gone=[x for x in zz if x['trialId'] in absent]
  assert len(kept)==c['retained']
  assert sum(not x['combinedFlag'] for x in kept)==c['retainedUnflagged']
  assert set(seq(c['retainedUnflaggedIds']))=={x['trialId'] for x in kept if not x['combinedFlag']}
  assert set(seq(c['removedUnflaggedIds']))=={x['trialId'] for x in gone if not x['combinedFlag']}
  for pool,prefix in [(kept,'retained'),(gone,'removed')]:
   assert sum(x['injectedBlink'] for x in pool)==c[prefix+'InjectedBlink']
   assert sum(not x['injectedBlink'] for x in pool)==c[prefix+'NoInjectedBlink']
 if d['selected']:
  assert d['rejectionIdsExplainAbsence']
  rejection={v for s in d['rejectionSteps'] if 'ureventIds' in s for v in seq(s['ureventIds'])}
  assert rejection==absent
  selected.append(d)
assert len(data)==1270 and len(selected)==60
for arm in ['baseline','variant']:
 for stage in ['initial','after_cutoffs','after_eog']:
  subset=[x for x in selected if x['subject']==5 and x['arm']==arm and x['stage']==stage]
  assert {x['seed'] for x in subset}==set(range(5))
  for x in subset:
   assert [c['retainedUnflagged'] for c in x['conditions']]==([1,3] if stage=='after_eog' else [2,3])
extra=[x for x in flags[5] if not x['injectedBlink'] and x['combinedFlag']]
assert {x['trialId'] for x in extra}=={3,31,38}
for x in extra:
 assert x['sampleHmmBlink']==1 and sum(x[k] for k in ['epochBlink','epochHorizontal','epochMuscle','sampleRobustBlink','sampleRobustHorizontal','sampleHmmHorizontal'])==0
summary={'candidateSnapshots':len(data),'selectedSnapshots':len(selected),'extraS005Flags':[x['trialId'] for x in extra],'extraFlagSource':'continuous EOG HMM only','s005SelectedSupportBeforeEOG':[2,3],'s005SelectedSupportAfterEOG':[1,3],'lostUnflaggedTrial':23,'allFiveStartsAgree':True,'productionChanged':False}
(p/'verification.json').write_text(json.dumps(summary,indent=2)+'\n')
print('INDEPENDENT_SUPPORT_VERIFICATION_PASSED '+json.dumps(summary))
