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
    unavailable=[]
    for recipe in ('standard','balanced'):
        with (folder/(recipe+'_calibration.csv')).open(newline='') as f:
            raw=list(csv.DictReader(f))
        groups=candidates(raw)
        arms[recipe]={}
        for record,cs in groups.items():
            picks=[c for c in cs if c['selected']]
            if len(picks)!=1:
                unavailable.append(dict(recording=record,arm=recipe,reason='No unique assessable recommendation'))
                continue
            pick=picks[0]
            score=min(cs,key=lambda c:c['reported_mse'])
            oracle=min(cs,key=lambda c:c['mse'])
            association=min(cs,key=lambda c:c['association_mse'])
            arms[recipe][record]=dict(selected=pick,score_minimum=score,oracle=oracle,association_minimum=association,
                candidates=len(cs),best_actual=sorted(cs,key=lambda c:c['mse'])[:4],
                best_reported=sorted(cs,key=lambda c:c['reported_mse'])[:4])
    paired=sorted(set(arms['standard'])&set(arms['balanced']))
    for record in sorted(set(arms['standard'])^set(arms['balanced'])):
        arm='standard' if record in arms['standard'] else 'balanced'
        unavailable.append(dict(recording=record,arm=arm,reason='Assessable in one arm only; paired comparison unavailable'))
    records=[]
    for record in paired:
        a,b=(arms[k][record] for k in ('standard','balanced'))
        records.append(dict(recording=record,difference=b['selected']['mse']-a['selected']['mse'],standard=a,balanced=b))
    avg=lambda get: mean(map(get,records)) if records else None
    analyses[population]=dict(records=sorted(records,key=lambda r:r['difference'],reverse=True),
        assessable_paired=len(records),excluded_saved_records=sorted(unavailable,key=lambda u:(u['recording'],u['arm'])),
        standard_mse=avg(lambda r:r['standard']['selected']['mse']),
        balanced_mse=avg(lambda r:r['balanced']['selected']['mse']),
        balanced_score_minimum_mse=avg(lambda r:r['balanced']['score_minimum']['mse']),
        balanced_oracle_mse=avg(lambda r:r['balanced']['oracle']['mse']),
        balanced_association_minimum_mse=avg(lambda r:r['balanced']['association_minimum']['mse']))
report=dict(source_run=40282756637,source_revision='4116ea11f3d5d39e4d3d282894dcf487172bd50d',
    scope='Post-hoc development attribution only. Truth never enters an operational selector.',analyses=analyses)
out.parent.mkdir(parents=True,exist_ok=True)
out.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
for population,data in analyses.items():
    print('ATTRIBUTION_TOTAL '+json.dumps(dict(population=population,**{k:v for k,v in data.items() if k not in ('records','excluded_saved_records')})))
    for row in data['records']:
        print('ATTRIBUTION_RECORD '+json.dumps(dict(population=population,**row),allow_nan=False))
    for row in data['excluded_saved_records']:
        print('ATTRIBUTION_UNAVAILABLE '+json.dumps(dict(population=population,**row)))
