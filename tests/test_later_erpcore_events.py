"""Task-event-only exploration never confuses trial count with clean EEG."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"benchmarks"))
from audit_later_erpcore_events import count_windows

def events():
    result=[]
    for sec,code in [(10,"11"),(100,"12"),(330,"22"),(415,"31")]:
        result.extend([
            {"onset":float(sec),"value":code,"trial_type":"stimulus"},
            {"onset":float(sec)+.5,"value":"201","trial_type":"response"}])
    return result

class ExtendedEventSupport(unittest.TestCase):
    def test_frozen_300s_counts_preserved(self):
        z=count_windows(events(),470)
        self.assertEqual(z[0]["windowSeconds"],300)
        self.assertEqual(z[0]["validTargetsBeforeEEGQC"],1)
        self.assertEqual(z[0]["validStandardsBeforeEEGQC"],1)
        self.assertEqual(z[0]["extraValidTargetsVs300"],0)

    def test_extra_events_not_labeled_artifact_free(self):
        z=count_windows(events(),470)
        by={r["windowSeconds"]:r for r in z}
        self.assertEqual(by[360]["extraValidTargetsVs300"],1)
        self.assertEqual(by[360]["extraValidStandardsVs300"],0)
        self.assertEqual(by[420]["extraValidTargetsVs300"],1)
        self.assertEqual(by[420]["extraValidStandardsVs300"],1)
        self.assertNotIn("hmmUnflaggedSupport",by[420])

    def test_unavailable_longer_duration_fails_closed(self):
        with self.assertRaises(ValueError):
            count_windows(events(),300)

    def test_windows_beyond_source_record_not_invented(self):
        z=count_windows(events(),375)
        self.assertEqual([q["windowSeconds"] for q in z],[300,360,375])

if __name__=="__main__":
    unittest.main()
