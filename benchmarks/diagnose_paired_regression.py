"""Attribute the completed paired development regression without fitting EEG."""
import csv,json,math,sys
from pathlib import Path
from statistics import mean
from diagnose_selection import candidates

root,out=map(Path,sys.argv[1:])
analyses={}
for population in ('public','synthetic'):
    folder=root/('development-optimization-'+population+'-summary')
    arms={}
    for recipe in ('standard','balanced'):
        with (folder/(recipe+'_calibration.csv')).open(newline='') as f:
            raw=list(csv.DictReader(f))
        groups=candidates(raw)
        arms[recipe]={}
        for record,cs in groups.items():
            picks=[c for c in cs if c['selected']]
            assert len(picks)==1, 'Expected one assessable recommendation'
            pick=picks[0]
            score=min(cs,key=lambda c:c['reported_mse'])
            oracle=min(cs,key=lambda c:c['mse'])
            association=min(cs,key=lambda c:c['association_mse'])
            arms[recipe][record]=dict(selected=pick,score_minimum=score,oracle=oracle,association_minimum=association,
                candidates=len(cs),best_actual=sorted(cs,key=lambda c:c['mse'])[:4],
                best_reported=sorted(cs,key=lambda c:c['reported_mse'])[:4])
    assert set(arms['standard'])==set(arms['balanced'])
    records=[]
    for record in arms['standard']:
        a,b=(arms[k][record] for k in ('standard','balanced'))
        records.append(dict(recording=record,difference=b['selected']['mse']-a['selected']['mse'],standard=a,balanced=b))
    analyses[population]=dict(records=sorted(records,key=lambda r:r['difference'],reverse=True),
        standard_mse=mean(r['standard']['selected']['mse'] for r in records),
        balanced_mse=mean(r['balanced']['selected']['mse'] for r in records),
        balanced_score_minimum_mse=mean(r['balanced']['score_minimum']['mse'] for r in records),
        balanced_oracle_mse=mean(r['balanced']['oracle']['mse'] for r in records),
        balanced_association_minimum_mse=mean(r['balanced']['association_minimum']['mse'] for r in records))
report=dict(source_run=40282756637,source_revision='4116ea11f3d5d39e4d3d282894dcf487172bd50d',
    scope='Post-hoc development attribution only. Truth never enters an operational selector.',analyses=analyses)
out.parent.mkdir(parents=True,exist_ok=True)
out.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
for population,data in analyses.items():
    print('ATTRIBUTION_TOTAL '+json.dumps(dict(population=population,**{k:v for k,v in data.items() if k!='records'})))
    for row in data['records']:
        print('ATTRIBUTION_RECORD '+json.dumps(dict(population=population,**row),allow_nan=False))
