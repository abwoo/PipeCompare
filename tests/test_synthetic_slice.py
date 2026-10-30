import importlib.util
import csv
import json
import tempfile
import unittest
from pathlib import Path
spec=importlib.util.spec_from_file_location('slices',Path(__file__).parents[1]/'benchmarks/prepare_synthetic_slice.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

class SliceTest(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name);self.source=self.root/'source';self.folder=self.source/'standard';self.folder.mkdir(parents=True)
        self.protocol={'recipes':['standard','balanced'],'development_synthetic_seeds':[40101], 'confirmation_synthetic_seeds':[41001]}
        self.e=dict(seed=40101,recipe='standard',completed=True,artifactVariance='legacy-total',nativeBackend=dict(
            matconvnetCommit=m.BACKEND,cpuBuildPatch=m.PATCH,sameWinningClasses=True,maxProbabilityDifference=1e-7,mex='/matlab/mex/vl_nnconv.mexa64',matlab='R2026a'))
        self.write()
        self.benchmark_fields='recording pipeline steps trueTarget trueStandard target standard errTarget errStandard errDifference kept badChannel lockedBlinkShare'.split()
        self.calibration_fields='recording steps condition actualError reportedError sme eyeBias gain selected conditionalError shapeSpread totalBias rawError rawConditionalError rawEstimate correctedEstimate populationTruth realizedTruth associationSensitivityError artifactBiasSE'.split()
        self.csv('benchmark.csv',self.benchmark_fields,[dict(recording=40101,pipeline=name) for name in ('PipeCompare','minimal','common')])
        self.csv('calibration.csv',self.calibration_fields,[dict(recording=40101,actualError=0.2)])
    def csv(self,name,fields,rows):
        with (self.folder/name).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
            for row in rows:w.writerow({key:row.get(key,0) for key in fields})
    def write(self):(self.folder/'execution.json').write_text(json.dumps(self.e))
    def prepare(self,**kw):
        args=dict(state=self.root/'state',source=self.source,seed=40101,recipe='standard',stage='development',revision='abc',protocol=self.protocol);args.update(kw)
        return m.prepare(**args)
    def test_completed_arm_is_reused_and_bound(self):
        self.assertTrue(self.prepare());self.assertTrue(self.prepare())
        with self.assertRaises(AssertionError):self.prepare(revision='changed')
    def test_partial_arm_is_not_reused(self):
        self.e['completed']=False;self.write();self.assertFalse(self.prepare())
        self.assertFalse((self.root/'state/standard/benchmark.csv').exists())
    def test_wrong_backend_and_identity_are_rejected(self):
        self.e['nativeBackend']['cpuBuildPatch']='fallback';self.write()
        with self.assertRaises(AssertionError):self.prepare()
        self.e['nativeBackend']['cpuBuildPatch']=m.PATCH;self.e['seed']=999;self.write()
        with self.assertRaises(AssertionError):m.completed(self.folder,40101,'standard')
    def test_incomplete_or_tampered_result_is_rejected(self):
        (self.folder/'calibration.csv').unlink()
        with self.assertRaises(FileNotFoundError):self.prepare()
        self.csv('calibration.csv',self.calibration_fields,[dict(recording=40101)])
        self.assertTrue(self.prepare(state=self.root/'fresh'))
        (self.root/'fresh/standard/benchmark.csv').write_text('changed')
        with self.assertRaises(AssertionError):self.prepare(state=self.root/'fresh')
    def test_confirmation_cannot_import_development(self):
        self.assertFalse(self.prepare(seed=41001,stage='confirmation'))
    def test_unidentified_checkpoint_is_refused(self):
        d=self.root/'state/standard';d.mkdir(parents=True);(d/'unexpected.mat').write_text('old')
        with self.assertRaises(AssertionError):self.prepare()
    def test_empty_calibration_is_valid_unavailable_not_missing(self):
        self.csv('calibration.csv',self.calibration_fields,[])
        self.assertTrue(self.prepare())

    def test_empty_digest_manifest_is_rejected(self):
        self.prepare()
        f=self.root/'state/standard/reused_control.json'
        proof=json.loads(f.read_text());proof['sha256']={};f.write_text(json.dumps(proof))
        with self.assertRaises(AssertionError):self.prepare()

    def test_partial_digest_manifest_is_rejected(self):
        self.prepare()
        f=self.root/'state/standard/reused_control.json'
        proof=json.loads(f.read_text());del proof['sha256']['calibration.csv'];f.write_text(json.dumps(proof))
        with self.assertRaises(AssertionError):self.prepare()

    def test_missing_provenance_is_rejected(self):
        self.prepare();(self.root/'state/standard/reused_control.json').unlink()
        with self.assertRaises(AssertionError):self.prepare()

    def test_wrong_source_commit_is_rejected_before_import(self):
        self.e['sourceCommit']='wrong';self.write()
        with self.assertRaises(AssertionError):self.prepare()

    def test_current_result_requires_matching_stage_and_commit(self):
        self.assertFalse(self.prepare(source=None))
        import shutil
        dst=self.root/'state/standard'
        self.e['sourceCommit']='abc';self.e['stage']='development';self.write()
        for name in m.RESULTS:shutil.copy2(self.folder/name,dst/name)
        self.assertTrue(self.prepare(source=None))
        self.e['stage']='confirmation';self.write();shutil.copy2(self.folder/'execution.json',dst/'execution.json')
        with self.assertRaises(AssertionError):self.prepare(source=None)

    def test_partial_current_execution_with_wrong_revision_is_rejected(self):
        self.assertFalse(self.prepare(source=None))
        self.e.update(completed=False,sourceCommit='wrong',stage='development')
        (self.root/'state/standard/execution.json').write_text(json.dumps(self.e))
        with self.assertRaises(AssertionError):self.prepare(source=None)

    def test_copy_failure_after_first_file_remains_retryable(self):
        from unittest.mock import patch
        copy=m.shutil.copy2
        calls=[]
        def interrupted(source,target):
            calls.append(source)
            if len(calls)==2:raise OSError('interrupted after execution')
            return copy(source,target)
        with patch.object(m.shutil,'copy2',side_effect=interrupted):
            with self.assertRaises(OSError):self.prepare()
        self.assertFalse((self.root/'state/standard/execution.json').exists())
        self.assertTrue(self.prepare())

    def test_failed_import_can_retry_without_partial_state(self):
        from unittest.mock import patch
        with patch.object(m.shutil,'copy2',side_effect=OSError('interrupted copy')):
            with self.assertRaises(OSError):self.prepare()
        self.assertTrue(self.prepare())

    def test_incomplete_csv_schema_is_rejected(self):
        (self.folder/'calibration.csv').write_text('recording,actualError\n40101,0.2\n')
        with self.assertRaises(AssertionError):self.prepare()

    def test_duplicate_header_and_truncated_row_are_rejected(self):
        fields=self.calibration_fields+[self.calibration_fields[-1]]
        self.csv('calibration.csv',fields,[])
        with self.assertRaises(AssertionError):self.prepare()
        self.csv('calibration.csv',self.calibration_fields,[])
        with (self.folder/'calibration.csv').open('a') as f:f.write('40101,only-two-fields\n')
        with self.assertRaises(AssertionError):self.prepare()

    def test_zero_byte_csv_is_rejected_with_clear_error(self):
        (self.folder/'calibration.csv').write_text('')
        with self.assertRaises(AssertionError):self.prepare()


    def test_changed_ranking_cannot_import_old_balanced_recommendations(self):
        import shutil
        self.folder=self.source/'balanced'
        shutil.copytree(self.source/'standard',self.folder)
        self.e.update(recipe='balanced',artifactVariance='pooled');self.write()
        self.protocol['reference_cleanup_guard']=True
        self.assertFalse(self.prepare(recipe='balanced'))
        self.assertFalse((self.root/'state/balanced/reused_control.json').exists())
        self.assertFalse((self.root/'state/balanced/benchmark.csv').exists())
        self.assertTrue(self.prepare(recipe='standard',state=self.root/'legacy'))


class CollectorTest(unittest.TestCase):
    prepare=SliceTest.prepare
    write=SliceTest.write
    csv=SliceTest.csv
    def setUp(self):
        SliceTest.setUp(self)
        import sys
        sys.path.insert(0,str(Path(__file__).parents[1]/'benchmarks'))
        from collect_synthetic_slices import collect
        self.collect=collect
        self.input=self.root/'slices';self.output=self.root/'out'
        self.protocol['recipes']=['standard']
        self.prepare(state=self.input/'development-synthetic-40101-standard-slice-6')
    def run_collect(self):self.collect(self.input,self.output,'development','abc',self.protocol)
    def test_missing_duplicate_and_partial_are_refused(self):
        import shutil
        self.run_collect()
        d=next(self.input.iterdir());shutil.copytree(d,self.input/'duplicate')
        with self.assertRaises(AssertionError):self.run_collect()
        shutil.rmtree(self.input/'duplicate')
        f=d/'standard/execution.json';e=json.loads(f.read_text());e['completed']=False;f.write_text(json.dumps(e))
        with self.assertRaises(AssertionError):self.run_collect()
        shutil.rmtree(d)
        with self.assertRaises(AssertionError):self.run_collect()
    def test_revision_and_seed_are_refused(self):
        d=next(self.input.iterdir());f=d/'standard/snapshot.json';t=json.loads(f.read_text());t['revision']='wrong';f.write_text(json.dumps(t))
        with self.assertRaises(AssertionError):self.run_collect()
        t['revision']='abc';t['seed']=999;f.write_text(json.dumps(t))
        with self.assertRaises(AssertionError):self.run_collect()

    def test_missing_arm_does_not_publish_partial_summary(self):
        self.protocol['development_synthetic_seeds'].append(40102)
        # Refresh the existing arm token for the expanded frozen plan.
        d=next(self.input.iterdir());f=d/'standard/snapshot.json'
        t=json.loads(f.read_text())
        import hashlib
        t['protocol_sha256']=hashlib.sha256(json.dumps(self.protocol,sort_keys=True).encode()).hexdigest()
        f.write_text(json.dumps(t))
        with self.assertRaises(AssertionError):self.run_collect()
        self.assertFalse(self.output.exists())

    def test_nonempty_output_is_not_mixed_with_new_results(self):
        self.output.mkdir();old=self.output/'old.csv';old.write_text('preserve')
        with self.assertRaises(AssertionError):self.run_collect()
        self.assertEqual(old.read_text(),'preserve')
        self.assertEqual(list(self.output.iterdir()),[old])

    def test_copy_failure_does_not_publish_partial_summary(self):
        from unittest.mock import patch
        with patch('collect_synthetic_slices.shutil.copy2',side_effect=OSError('copy interrupted')):
            with self.assertRaises(OSError):self.run_collect()
        self.assertFalse(self.output.exists())
        self.run_collect()
        self.assertTrue((self.output/'optimization-synthetic-40101/standard/benchmark.csv').exists())

if __name__=='__main__':unittest.main()

