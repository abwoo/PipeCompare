"""Restore only complete, pinned CPU-MEX controls; validate every resumed arm."""
import argparse
import csv
import hashlib
import json
import os
import shutil
import tempfile
from pathlib import Path

SOURCE_RUN = 40159157180
SOURCE_HEAD = '5e22e7534be739dea0d7a41ef4b95c92038f51fc'
BACKEND = '092ecb8c04031ddcb1d9a7ba6ebe76a250c4672c'
PATCH = 'cpu-gcc12-unused-gpu-include-and-explicit-limits-header'
RESULTS = ('execution.json', 'benchmark.csv', 'calibration.csv')
# Match the writers in pipecompare_benchmark and Calibration.empty. An empty
# calibration table is valid unavailable evidence; a missing schema is not.
CSV_FIELDS = {
    'benchmark.csv': set('recording pipeline steps trueTarget trueStandard target standard errTarget errStandard errDifference kept badChannel lockedBlinkShare'.split()),
    'calibration.csv': set('recording steps condition actualError reportedError sme eyeBias gain selected conditionalError shapeSpread totalBias rawError rawConditionalError rawEstimate correctedEstimate populationTruth realizedTruth associationSensitivityError artifactBiasSE'.split()),
}

def completed(folder, seed, recipe):
    f = folder / 'execution.json'
    if not f.exists():
        return False
    e = json.loads(f.read_text())
    assert e['seed'] == seed and e['recipe'] == recipe, 'Wrong execution identity'
    if e['completed'] is not True:
        return False
    assert e['artifactVariance'] == ('legacy-total' if recipe == 'standard' else 'pooled')
    b = e['nativeBackend']
    assert b['matconvnetCommit'] == BACKEND and b['cpuBuildPatch'] == PATCH, 'Wrong backend'
    assert b['sameWinningClasses'] is True and 0 <= b['maxProbabilityDifference'] < 1e-5
    assert b['mex'].endswith('vl_nnconv.mexa64'), 'Native MEX evidence missing'
    assert '2026' in b['matlab'], 'Wrong MATLAB release'
    for name in ('benchmark.csv', 'calibration.csv'):
        with (folder / name).open(newline='', encoding='utf-8-sig') as f:
            reader = csv.DictReader(f)
            fields = reader.fieldnames
            assert fields and len(fields) == len(set(fields)) and CSV_FIELDS[name] <= set(fields), 'Missing or duplicate result columns: ' + name
            rows = list(reader)
        assert all(None not in r and all(v is not None for v in r.values()) for r in rows), 'Malformed result row: ' + name
        assert all(float(r['recording']) == seed for r in rows), 'Wrong result seed'
        if name == 'benchmark.csv':
            assert len(rows) == 3 and {r['pipeline'] for r in rows} == {'PipeCompare', 'minimal', 'common'}
    return True

def prepare(state, source, seed, recipe, stage, revision, protocol):
    assert stage in ('development', 'confirmation') and recipe in protocol['recipes']
    assert seed in protocol[stage + '_synthetic_seeds'] and revision
    folder = state / recipe
    folder.mkdir(parents=True, exist_ok=True)
    token = dict(seed=seed, recipe=recipe, stage=stage, revision=revision,
                 protocol_sha256=hashlib.sha256(json.dumps(protocol, sort_keys=True).encode()).hexdigest())
    snapshot = folder / 'snapshot.json'
    if snapshot.exists():
        assert json.loads(snapshot.read_text()) == token, 'Different checkpoint experiment'
    else:
        assert not any(folder.iterdir()), 'Unidentified checkpoint contents'
        # Publish the snapshot and imported evidence together; interrupted copies
        # must remain retryable and must never look like a completed current arm.
        with tempfile.TemporaryDirectory(prefix='.import-', dir=state) as temp:
            pending = Path(temp) / recipe
            pending.mkdir()
            if stage == 'development' and source is not None and completed(source / recipe, seed, recipe):
                execution = json.loads((source / recipe / 'execution.json').read_text())
                assert execution.get('sourceCommit', SOURCE_HEAD) == SOURCE_HEAD and execution.get('stage', stage) == stage, 'Wrong imported execution source'
                for name in RESULTS + ('benchmark_summary.txt',):
                    f = source / recipe / name
                    if f.exists():
                        shutil.copy2(f, pending / name)
                provenance = dict(sourceRun=SOURCE_RUN, sourceCommit=SOURCE_HEAD, targetRevision=revision,
                                  reason='Completed frozen CPU-MEX arm; execution-only recovery after parity tests',
                                  sha256={name: hashlib.sha256((pending/name).read_bytes()).hexdigest() for name in RESULTS})
                (pending / 'reused_control.json').write_text(json.dumps(provenance, indent=2))
            (pending / 'snapshot.json').write_text(json.dumps(token))
            folder.rmdir()  # Refuse unexpected concurrent writes.
            pending.rename(folder)
    proof = folder / 'reused_control.json'
    if proof.exists():
        assert stage == 'development', 'Confirmation cannot reuse development results'
        p = json.loads(proof.read_text())
        assert set(p['sha256']) == set(RESULTS), 'Incomplete imported digest manifest'
        assert p['sourceRun'] == SOURCE_RUN and p['sourceCommit'] == SOURCE_HEAD and p['targetRevision'] == revision
        assert all(hashlib.sha256((folder/name).read_bytes()).hexdigest() == digest for name, digest in p['sha256'].items()), 'Imported result modified'
    execution_file = folder / 'execution.json'
    if execution_file.exists():
        execution = json.loads(execution_file.read_text())
        if proof.exists():
            # The pinned 5e22 source predates embedded revision/stage fields.
            # Its run identity is checked by the workflow; the full digest proof
            # binds those unchanged legacy bytes to this development snapshot.
            assert execution.get('sourceCommit', SOURCE_HEAD) == SOURCE_HEAD and execution.get('stage', stage) == stage, 'Wrong imported execution source'
        else:
            assert execution.get('sourceCommit') == revision and execution.get('stage') == stage, 'Wrong execution source or missing provenance'
    return completed(folder, seed, recipe)

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--state', type=Path, required=True);p.add_argument('--source', type=Path)
    p.add_argument('--seed', type=int, required=True);p.add_argument('--recipe', required=True)
    p.add_argument('--stage', required=True);p.add_argument('--revision', required=True)
    p.add_argument('--protocol', type=Path, required=True)
    p.add_argument('--compact-completed', action='store_true')
    a = p.parse_args()
    done = prepare(a.state,a.source,a.seed,a.recipe,a.stage,a.revision,json.loads(a.protocol.read_text()))
    if done and a.compact_completed:
        shutil.rmtree(a.state/a.recipe/'checkpoints',ignore_errors=True)
    if os.environ.get('GITHUB_OUTPUT'):
        with open(os.environ['GITHUB_OUTPUT'],'a') as f:f.write('complete='+str(done).lower()+'\n')
    print(json.dumps(dict(seed=a.seed,recipe=a.recipe,completed=done)))

if __name__ == '__main__':main()

