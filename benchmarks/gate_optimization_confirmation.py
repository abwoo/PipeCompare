"""Admit a frozen confirmation study only after development execution/controls.
This engineering screen is not a significance test or scientific success claim.
"""
import argparse
import json
import math
import os
from pathlib import Path

def gate(public, synthetic, criteria):
    assert public['stage']==synthetic['stage']=='development', 'Wrong source stage'
    reasons=[];metrics={}
    for name, report in [('public',public),('synthetic',synthetic)]:
        if report['paired']<criteria['minimum_paired_development_records']:
            reasons.append(name+': insufficient paired development records')
        if name=='public':
            old,new=report['standardAvailable'],report['balancedAvailable']
        else:
            old=report['planned']-report['standard']['unavailable']
            new=report['planned']-report['balanced']['unavailable']
        if new<old:
            reasons.append(name+': recommendation coverage declined')
        before=report['standardPairedRMSE']**2;after=report['balancedPairedRMSE']**2
        if not all(math.isfinite(v) and v>=0 for v in [before,after]):
            reasons.append(name+': nonfinite paired loss');relative=None
        elif before==0:
            relative=0 if after==0 else None
            if after>0:reasons.append(name+': worsened a zero-loss control')
        else:
            relative=(after-before)/before
            if relative>criteria['maximum_relative_paired_mse_increase']:
                reasons.append(name+': paired MSE exceeds the development tolerance')
        metrics[name+'_relative_paired_mse_change']=relative
    short={};outside={}
    for mode in ('long','multiscale'):
        selected=[r for r in public['stretches'] if r['mode']==mode and r['seconds'] in (1,5)]
        assert len(selected)==2, 'Missing short-fault development controls'
        values=[r['meanRecall'] for r in selected]
        assert all(v is not None and math.isfinite(v) for v in values)
        short[mode]=sum(values)/len(values)
        values=[r['outsideInjectedShare'] for r in public['stretches'] if r['mode']==mode]
        outside[mode]=sum(values)/len(values)
    gain=short['multiscale']-short['long'];metrics['short_fault_recall_gain']=gain
    excess=outside['multiscale']-outside['long'];metrics['additional_outside_injection_marking']=excess
    if excess>criteria['maximum_additional_outside_injection_marking']:
        reasons.append('short detector: excessive extra marking outside injections')
    accuracy_gain=any(metrics[n+'_relative_paired_mse_change'] is not None and
                      metrics[n+'_relative_paired_mse_change']<0 for n in ('public','synthetic'))
    if not accuracy_gain and gain<criteria['minimum_short_fault_recall_gain_if_no_accuracy_gain']:
        reasons.append('no development accuracy gain or substantial short-fault detection gain')
    return dict(ready_for_confirmation=not reasons,reasons=reasons,metrics=metrics,criteria=criteria,
        scope='Development engineering screen only; held-out outcomes remain unknown.')

def main():
    parser=argparse.ArgumentParser();parser.add_argument('root',type=Path)
    parser.add_argument('out',type=Path);args=parser.parse_args()
    criteria=json.loads(Path(__file__).with_name('optimization_confirmation_gate.json').read_text())
    public=json.loads((args.root/'optimization-public-summary'/'summary.json').read_text())
    synthetic=json.loads((args.root/'optimization-synthetic-summary'/'summary.json').read_text())
    result=gate(public,synthetic,criteria);args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    if os.environ.get('GITHUB_OUTPUT'):
        with open(os.environ['GITHUB_OUTPUT'],'a') as f:
            f.write('ready='+str(result['ready_for_confirmation']).lower()+'\n')
    print('PIPECOMPARE_CONFIRMATION_ADMISSION '+json.dumps(result,allow_nan=False))

if __name__=='__main__':main()
