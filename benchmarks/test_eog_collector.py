import unittest,copy
from collect_eog import summarize,read,ROOT
class EogGate(unittest.TestCase):
 def setUp(self):
  self.rows=[{'subject':x['subject'],'seed':x['seed'],'baseline':x['within'],'variant':copy.deepcopy(x['within']),'original':x['original']} for x in read(ROOT/'eog_reference.json')]
 def test_unchanged_rejects(self):self.assertEqual(summarize(self.rows)['decision'],'reject-expansion')
 def test_missing_rejects(self):
  with self.assertRaises(AssertionError):summarize(self.rows[:-1])
 def test_recovered_passes(self):
  for x in self.rows:
   for k in ['mse','rawMse']:x['variant']['erp'][k]=.1
  self.assertEqual(summarize(self.rows)['decision'],'proceed-to-ten-person-guardrail')
 def test_one_bad_start_rejects(self):
  for x in self.rows:
   for k in ['mse','rawMse']:x['variant']['erp'][k]=.1
  self.rows[0]['variant']['erp']['mse']=99
  self.assertEqual(summarize(self.rows)['decision'],'reject-expansion')
if __name__=='__main__':unittest.main()
