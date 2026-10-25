import importlib.util
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
        (self.folder/'benchmark.csv').write_text('recording,pipeline\n40101,PipeCompare\n40101,minimal\n40101,common\n')
        (self.folder/'calibration.csv').write_text('recording,actualError\n40101,0.2\n')
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
        (self.folder/'calibration.csv').write_text('recording\n40101\n')
        self.assertTrue(self.prepare(state=self.root/'fresh'))
        (self.root/'fresh/standard/benchmark.csv').write_text('changed')
        with self.assertRaises(AssertionError):self.prepare(state=self.root/'fresh')
    def test_confirmation_cannot_import_development(self):
        self.assertFalse(self.prepare(seed=41001,stage='confirmation'))
    def test_unidentified_checkpoint_is_refused(self):
        d=self.root/'state/standard';d.mkdir(parents=True);(d/'unexpected.mat').write_text('old')
        with self.assertRaises(AssertionError):self.prepare()
    def test_empty_calibration_is_valid_unavailable_not_missing(self):
        (self.folder/'calibration.csv').write_text('recording,actualError\n')
        self.assertTrue(self.prepare())


class CollectorTest(unittest.TestCase):
    prepare=SliceTest.prepare
    write=SliceTest.write
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

if __name__=='__main__':unittest.main()
