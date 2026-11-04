import csv,json,math,collections
from pathlib import Path
p=Path('input')
data=json.loads((p/'diagnosis.json').read_text())
rows=list(csv.DictReader((p/'invalid_draws.csv').open()))
def seq(x): return x if isinstance(x,list) else [x]
by=collections.defaultdict(list)
for r in rows: by[int(r['subject']),int(r['seed'])].append(r)
output=[]
for d in data:
    rr=by[d['subject'],d['seed']]
    expected=set()
    detail=[]
    for c in d['candidateDiagnostics']:
        for cond in c['conditions']:
            reps=seq(cond['badReplicates'])
            causes=seq(cond['causes'])
            assert len(reps)==len(causes)
            for rep,cause in zip(reps,causes):
                expected.add((c['id'],cond['condition'],rep,cause))
        if c['badScoreCount']:
            detail.append({'id':c['id'],'pipeline':c['pipeline'],'badScoreCount':c['badScoreCount'],'conditions':[v for v in c['conditions'] if v['badEyeCount'] or v['badSMECount']]})
    observed={(int(r['candidate']),int(r['condition']),int(r['replicate']),r['cause']) for r in rr}
    assert expected==observed and len(observed)==len(rr)
    scorekeys={(int(r['candidate']),int(r['replicate'])) for r in rr}
    assert len(scorekeys)==d['badScoreCells']
    assert len({k[1] for k in scorekeys})==d['badBootstrapColumns']
    causes=collections.Counter(r['cause'] for r in rr)
    q={k:v for k,v in d.items() if k!='candidateDiagnostics'}
    q['causeCounts']=dict(causes)
    q['oldFormulaInvalidRows']=sum(not math.isfinite(float(r['oldFormulaEyeSq'])) for r in rr)
    q['badSMERows']=sum(not math.isfinite(float(r['sme'])) for r in rr)
    q['badEyeRows']=sum(not math.isfinite(float(r['eyeSq'])) for r in rr)
    q['invalidRowCount']=len(rr)
    q['badCandidateIds']=[v['id'] for v in detail]
    q['selectedDiagnostics']=[]
    for c in d['candidateDiagnostics']:
        if c['id'] in [d['incumbentId'],d['challengerId']]:
            short={'id':c['id'],'pipeline':c['pipeline'],'badScoreCount':c['badScoreCount'],'conditions':[]}
            for v in c['conditions']:
                short['conditions'].append({k:w for k,w in v.items() if k not in ['badReplicates','causes']})
            q['selectedDiagnostics'].append(short)
    q['badCandidateSupport']=[{'id':v['id'],'conditions':[{k:w for k,w in z.items() if k not in ['badReplicates','causes']} for z in v['conditions']]} for v in detail]
    output.append(q)
    print('VERIFIED_DIAGNOSIS '+json.dumps(q,separators=(',',':')))
assert len(data)==10 and len(rows)==101496
Path('summary.json').write_text(json.dumps(output,indent=2))
print('INDEPENDENT_DIAGNOSTIC_COUNTS_PASSED records=10 rows='+str(len(rows)))
