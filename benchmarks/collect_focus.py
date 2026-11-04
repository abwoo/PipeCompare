"""Validate every diagnostic pair; ICA starts are repetitions, not participants."""
import json,sys,math,hashlib,statistics as st
from pathlib import Path
ROOT=Path(__file__).parent
STAGES=['initial','after_cutoffs','after_eog']
def read(p):return json.loads(Path(p).read_text())
def close(a,b,tol=1e-7):
    assert isinstance(a,(int,float)) and isinstance(b,(int,float)) and math.isfinite(a) and math.isfinite(b)
    assert abs(a-b)<=tol,(a,b)
def outcome(e,truth):
    assert e['available'] and len(e['estimate'])==len(truth)==2
    g=e['gain'];assert isinstance(g,(int,float)) and math.isfinite(g) and g>0
    corrected=[x/g for x in e['estimate']]
    close(st.mean((x-t)**2 for x,t in zip(corrected,truth)),e['mse'])
    close(st.mean((x-t)**2 for x,t in zip(e['estimate'],truth)),e['rawMse'])
    for x,y in zip(corrected,e['correctedEstimate']):close(x,y)
def stage_check(d,source,subject,seed,stage,arm=None):
    assert d['subject']==subject and d['seed']==seed and d['stage']==stage
    assert d['populationTruth']==source['populationTruth'] and d['realizedTruth']==source['realizedTruth']
    if arm is not None:assert d['arm']==arm
    rows=d['rows'];assert len(rows)==len({q['pipeline'] for q in rows})
    assert len(rows)==len({q['id'] for q in rows})
    rec=[q for q in rows if q['selected']];assert len(rec)==1
    assert rec[0]['id']==d['recommended'] and rec[0]['ranking']['status']=='feasible'
    for q in rows:
        if q['available']:outcome(q['readout'],d['populationTruth'])
    return rec[0]
def pair(folder):
    arms={a:read(folder/a/'result.json') for a in ['baseline','variant']}
    b,v=arms.values();subject=b['subject'];seed=b['seed'];assert subject in [4,5] and seed in range(5)
    assert v['subject']==subject and v['seed']==seed
    prior={p['subject']:p for p in read(ROOT/'results/within_trial_20261104/records.json')}[subject]
    for key in ['source','populationTruth','realizedTruth','blinkPrevalence','searchOptions','eeglabSha','executionCommit','workflowRun']:
        assert b[key]==v[key],key
    assert b['source']==prior['baseline']['source']
    assert b['populationTruth']==prior['baseline']['populationTruth'] and b['realizedTruth']==prior['baseline']['realizedTruth']
    assert b['eeglabSha']==prior['baseline']['eeglabSha']
    expected=dict(prior['baseline']['searchOptions']);expected['icaSeeds']=0
    assert b['searchOptions']==expected,'Unexpected option change'
    traces={};counter={}
    for a,d in arms.items():
        assert d['arm']==a and d['completed'] and not d['errors']
        close(d['replayMaxError'],0,1e-6);outcome(d['erp'],d['populationTruth'])
        if seed==0:
            for k in ['mse','rawMse']:close(d['erp'][k],prior[a]['erp'][k])
            assert d['erp']['pipeline']==prior[a]['erp']['pipeline'],'Seed-zero selected path drift'
        traces[a]={}
        for stage in STAGES:
            t=read(folder/a/(stage+'.json'));pick=stage_check(t,d,subject,seed,stage,a)
            if stage=='after_eog':
                close(pick['readout']['mse'],d['erp']['mse']);assert pick['pipeline']==d['erp']['pipeline']
            traces[a][stage]=t
        for score in ['baseline','variant']:
            name=a+'_scored_'+score;counter[name]={}
            for stage in STAGES:
                t=read(folder/name/(stage+'.json'));assert t['sourceArm']==a
                stage_check(t,d,subject,seed,stage)
                old=traces[a][stage]
                assert [q['pipeline'] for q in old['rows']]==[q['pipeline'] for q in t['rows']]
                for before,after in zip(old['rows'],t['rows']):
                    assert before['available']==after['available']
                    if before['available']:
                        assert before['readout']['estimate']==after['readout']['estimate']
                        close(before['readout']['mse'],after['readout']['mse'])
                    if a==score:
                        assert before['selected']==after['selected'] and before['ranking']['status']==after['ranking']['status']
                        for k in ['objective','shrunk']:
                            if before['ranking'][k] is not None:close(before['ranking'][k],after['ranking'][k],1e-9)
                counter[name][stage]=t
    return {'subject':subject,'seed':seed,'baseline':b,'variant':v,'traces':traces,'counterfactuals':counter}
def summarize(records):
    assert len(records)==10 and {(r['subject'],r['seed']) for r in records}=={(s,k) for s in [4,5] for k in range(5)},'Missing or duplicate pairs'
    assert len({r['baseline']['executionCommit'] for r in records})==1
    out={'pairs':10,'independentParticipants':2,'independentConfirmation':False,'subjects':[]}
    for subject in [4,5]:
        rr=sorted([r for r in records if r['subject']==subject],key=lambda r:r['seed'])
        rows=[{'seed':r['seed'],'baselineMSE':r['baseline']['erp']['mse'],'variantMSE':r['variant']['erp']['mse'],
            'difference':r['variant']['erp']['mse']-r['baseline']['erp']['mse'],
            'baselineRawMSE':r['baseline']['erp']['rawMse'],'variantRawMSE':r['variant']['erp']['rawMse'],
            'baselinePipeline':r['baseline']['erp']['pipeline'],'variantPipeline':r['variant']['erp']['pipeline']} for r in rr]
        d=[r['difference'] for r in rows];raw=st.mean(r['variantRawMSE']-r['baselineRawMSE'] for r in rows)
        out['subjects'].append({'subject':subject,'rows':rows,'meanDifference':st.mean(d),'medianDifference':st.median(d),
            'differenceRange':[min(d),max(d)],'wins':sum(x<-1e-9 for x in d),'losses':sum(x>1e-9 for x in d),'ties':sum(abs(x)<=1e-9 for x in d),
            'meanRawDifference':raw,'benefitStableAcrossTestedStarts':all(x<-1e-9 for x in d) and raw<=1e-9})
    return out
if __name__=='__main__':
    plan=read(ROOT/'focus_plan.json')
    for path,key in [('within_variant/Measure.m','variant_measure_sha256'),('focus_harness/Executor.m','executor_harness_sha256')]:
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==plan[key]
    if sys.argv[1]=='--pair':
        r=pair(Path(sys.argv[2]));print('FOCUS_PAIR_VALIDATED',r['subject'],r['seed'])
    else:
        out=Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=True)
        try:
            records=[pair(p) for p in sorted(Path(sys.argv[1]).glob('focus-pair-*'))]
            s=summarize(records)
            (out/'records.json').write_text(json.dumps(records,separators=(',',':'))+'\n')
            (out/'summary.json').write_text(json.dumps(s,indent=2)+'\n')
            print('PIPECOMPARE_FOCUS_SUMMARY '+json.dumps(s))
        except Exception as ex:
            (out/'failure.txt').write_text(repr(ex));raise
