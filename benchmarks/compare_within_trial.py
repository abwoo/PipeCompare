"""Frozen paired development comparison; never selects a variant using truth per person."""
import sys,json,math,hashlib
from pathlib import Path
import numpy as np
from summarize_public_pilot import validate_record,interval
ROOT=Path(__file__).parent

def compare(folder):
    protocol=json.loads((ROOT/'public_pilot_protocol.json').read_text())
    digest=hashlib.sha256((ROOT/'public_pilot_protocol.json').read_bytes()).hexdigest()
    plan=json.loads((ROOT/'within_trial_plan.json').read_text())
    assert hashlib.sha256((ROOT/'within_variant/Measure.m').read_bytes()).hexdigest()==plan['variant_measure_sha256']
    original={r['subject']:r for r in json.loads((ROOT/'results/public_pilot_20261103/records.json').read_text())}
    rows=[];records=[];execution=set()
    for subject in range(1,11):
        b=json.loads((folder/f'ranking-pair-{subject}/baseline/result.json').read_text())
        v=json.loads((folder/f'ranking-pair-{subject}/variant/result.json').read_text())
        assert b['subject']==v['subject']==subject
        validate_record(b,protocol,digest)
        vp=dict(protocol,algorithm_commit=v['executionCommit']);validate_record(v,vp,digest)
        assert b['completed'] and v['completed'] and b['erp']['available'] and v['erp']['available']
        assert b['executionCommit']==v['executionCommit'];execution.add(v['executionCommit'])
        assert b['source']==v['source'] and b['searchOptions']==v['searchOptions']
        assert b['populationTruth']==v['populationTruth'] and b['realizedTruth']==v['realizedTruth']
        assert b['eeglabSha']==v['eeglabSha']
        for arm in ['minimal','common']:
            assert np.allclose(b[arm]['estimate'],v[arm]['estimate'],atol=1e-7,rtol=0), 'Fixed control changed'
        assert math.isclose(b['erp']['mse'],original[subject]['erp']['mse'],abs_tol=1e-7,rel_tol=0),'Baseline drift'
        rows.append(dict(subject=subject,baselineMSE=b['erp']['mse'],variantMSE=v['erp']['mse'],difference=v['erp']['mse']-b['erp']['mse'],baselineRawMSE=b['erp']['rawMse'],variantRawMSE=v['erp']['rawMse'],baselineStable=b['icaSeedCheck']['same'],variantStable=v['icaSeedCheck']['same']))
        records.append(dict(subject=subject,baseline=b,variant=v))
    assert len(execution)==1
    d=np.array([r['difference'] for r in rows]);idx=np.random.RandomState(510003).randint(0,10,(10000,10))
    bm=np.mean([r['baselineMSE'] for r in rows]);vm=np.mean([r['variantMSE'] for r in rows]);br=np.mean([r['baselineRawMSE'] for r in rows]);vr=np.mean([r['variantRawMSE'] for r in rows])
    retain=bool(vm<bm-1e-12 and vr<=br+1e-12)
    summary=dict(pairs=10,executionCommit=next(iter(execution)),baselineMSE=float(bm),variantMSE=float(vm),difference=float(d.mean()),ci95=interval(d[idx].mean(1)),baselineRawMSE=float(br),variantRawMSE=float(vr),wins=int((d<-1e-9).sum()),losses=int((d>1e-9).sum()),ties=int((abs(d)<=1e-9).sum()),decision='retain-exploratory-variant' if retain else 'restore-original',independentConfirmation=False,rows=rows)
    return summary,records

if __name__=='__main__':
    out=Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=True)
    try:
        summary,records=compare(Path(sys.argv[1]));(out/'records.json').write_text(json.dumps(records,indent=2)+'\n');(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print('PIPECOMPARE_RANKING_COMPARISON '+json.dumps(summary))
    except Exception as ex:
        (out/'failure.txt').write_text(repr(ex));raise
