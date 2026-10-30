"""Keep development inputs and locked confirmation inputs disjoint."""
import json
import unittest
from pathlib import Path

class FrozenProtocolTests(unittest.TestCase):
    def setUp(self):
        root=Path(__file__).resolve().parents[1]
        self.new=json.loads((root/'benchmarks/optimization_protocol.json').read_text())
        self.old=json.loads((root/'benchmarks/open_data_protocol.json').read_text())

    def test_confirmation_has_no_previously_used_participant(self):
        seen={r['subject'] for r in self.old['physionet']['files']}
        confirm=self.new['confirmation_subjects']
        self.assertEqual(len(confirm),len(set(confirm)))
        self.assertFalse(seen.intersection(confirm))
        self.assertFalse(set(self.new['development_subjects']).intersection(confirm))

    def test_training_and_test_partition_the_new_cohort(self):
        train=set(self.new['confirmation_training_subjects'])
        test=set(self.new['confirmation_test_subjects'])
        self.assertFalse(train.intersection(test))
        self.assertEqual(train|test,set(self.new['confirmation_subjects']))
        self.assertEqual(len(test),20)

    def test_new_simulations_are_disjoint_and_not_human_units(self):
        train=self.new['development_synthetic_seeds'];test=self.new['confirmation_synthetic_seeds']
        self.assertFalse(set(train).intersection(test))
        self.assertEqual(len(test),30);self.assertEqual(len(set(test)),30)
        self.assertTrue(all(seed>20060 for seed in train+test))
        self.assertIn('never candidate rows',self.new['evaluation']['units'])

if __name__=='__main__':
    unittest.main()
