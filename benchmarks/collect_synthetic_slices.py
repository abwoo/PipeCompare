"""Require both recipes for every frozen seed before producing summary inputs."""
import argparse
import json
import shutil
from pathlib import Path
from prepare_synthetic_slice import prepare

def collect(root,out,stage,revision,protocol):
    expected={(seed,recipe) for seed in protocol[stage+'_synthetic_seeds'] for recipe in protocol['recipes']}
    found=set()
    for snapshot in root.glob('*/standard/snapshot.json'):
        found.add(snapshot)
    for snapshot in root.glob('*/balanced/snapshot.json'):
        found.add(snapshot)
    seen=set()
    for snapshot in sorted(found):
        token=json.loads(snapshot.read_text());key=(token['seed'],token['recipe'])
        assert key in expected and key not in seen, 'Unexpected or duplicate simulation arm'
        assert token['stage']==stage and token['revision']==revision, 'Wrong source revision or stage'
        assert prepare(snapshot.parent.parent,None,*key,stage,revision,protocol), f'Incomplete simulation fit: {key}'
        seen.add(key)
        target=out/f'optimization-synthetic-{key[0]}'/key[1];target.mkdir(parents=True,exist_ok=True)
        for name in ('execution.json','benchmark.csv','calibration.csv','reused_control.json','snapshot.json'):
            f=snapshot.parent/name
            if f.exists():shutil.copy2(f,target/name)
    assert seen==expected, f'Missing planned simulation arms: {expected-seen}'
    print(f'Complete: {len(seen)} arms from {len({key[0] for key in expected})} frozen seeds')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--stage',required=True);p.add_argument('--revision',required=True);p.add_argument('--protocol',type=Path,required=True)
    a=p.parse_args();collect(a.root,a.out,a.stage,a.revision,json.loads(a.protocol.read_text()))
