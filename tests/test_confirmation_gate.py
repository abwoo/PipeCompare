import json
import sys
import unittest
from copy import deepcopy
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'benchmarks'))
from gate_optimization_confirmation import gate

class GateTests(unittest.TestCase):
    def setUp(self):
        self.c=json.loads((Path(__file__).resolve().parents[1]/'benchmarks/optimization_confirmation_gate.json').read_text())
        self.p=dict(stage='development',paired=4,standardAvailable=4,balancedAvailable=4,
            standardPairedRMSE=2,balancedPairedRMSE=1.8,stretches=[])
        for mode in ('long','multiscale'):
            for seconds in (1,5):
                self.p['stretches'].append(dict(mode=mode,seconds=seconds,meanRecall=0 if mode=='long' else 1,outsideInjectedShare=0))
        self.s=dict(stage='development',paired=4,planned=4,standard=dict(unavailable=0),balanced=dict(unavailable=0),
            standardPairedRMSE=1,balancedPairedRMSE=1)
    def test_accuracy_cannot_hide_coverage_loss(self):
        p=deepcopy(self.p);p['balancedAvailable']=3
        self.assertFalse(gate(p,self.s,self.c)['ready_for_confirmation'])
    def test_admission_is_not_granted_to_a_worse_search(self):
        p=deepcopy(self.p);p['balancedPairedRMSE']=3
        self.assertFalse(gate(p,self.s,self.c)['ready_for_confirmation'])
    def test_valid_development_can_only_admit_confirmation(self):
        result=gate(self.p,self.s,self.c)
        self.assertTrue(result['ready_for_confirmation'])
        self.assertIn('outcomes remain unknown',result['scope'])
        p=deepcopy(self.p);p['stage']='confirmation'
        with self.assertRaises(AssertionError):gate(p,self.s,self.c)

if __name__=='__main__':unittest.main()
