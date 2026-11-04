import unittest,copy,json
from pathlib import Path
from collect_focus import outcome,summarize
class FocusChecks(unittest.TestCase):
 def setUp(self):
  self.prior=json.loads((Path(__file__).parent/'results/within_trial_20261104/records.json').read_text())
  self.records=[{'subject':s,'seed':k,'baseline':self.prior[s-1]['baseline'],'variant':self.prior[s-1]['variant']} for s in [4,5] for k in range(5)]
 def test_known_outcomes(self):
  for r in self.prior:
   for arm in ['baseline','variant']:outcome(r[arm]['erp'],r[arm]['populationTruth'])
 def test_reject_tampered_error(self):
  e=copy.deepcopy(self.prior[3]['variant']['erp']);e['mse']+=.01
  with self.assertRaises(AssertionError):outcome(e,self.prior[3]['variant']['populationTruth'])
 def test_reject_undefined(self):
  e=copy.deepcopy(self.prior[3]['variant']['erp']);e['gain']=float('nan')
  with self.assertRaises(AssertionError):outcome(e,self.prior[3]['variant']['populationTruth'])
 def test_missing_pair(self):
  with self.assertRaises(AssertionError):summarize(self.records[:-1])
 def test_duplicate_pair(self):
  with self.assertRaises(AssertionError):summarize(self.records[:-1]+[self.records[0]])
 def test_one_reversal_is_not_stable(self):
  rr=copy.deepcopy(self.records);rr[1]['variant']['erp']['mse']=rr[1]['baseline']['erp']['mse']+1
  self.assertFalse(summarize(rr)['subjects'][0]['benefitStableAcrossTestedStarts'])
if __name__=='__main__':unittest.main()
