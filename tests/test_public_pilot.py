import copy
import hashlib
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"benchmarks"))
from fetch_public_pilot import read_edf
from summarize_public_pilot import paired_summary,calibration,summarize

PFILE=Path(__file__).resolve().parents[1]/"benchmarks"/"public_pilot_protocol.json"
P=json.loads(PFILE.read_text())
DIGEST=hashlib.sha256(PFILE.read_bytes()).hexdigest()
def arm(error):
    v=[float(error),float(error)]
    return dict(available=True,estimate=v,correctedEstimate=v,actualError=v,rawError=v,conditionalError=v,
                reportedError=[1.,1.],mse=error**2,rawMse=error**2,conditionalMse=error**2,reportedMse=1.,gain=1.)
def record(subject=1,error=1):
    source=next(x for x in P["source"]["files"] if x["subject"]==subject)
    return dict(subject=subject,algorithmCommit=P["algorithm_commit"],protocolSha256=DIGEST,
                executionCommit="one-frozen-execution",source=dict(subject=subject,source_sha256=source["sha256"],sampling_rate=160,samples=9600),
                completed=True,notRun=0,errors=[],populationTruth=[0.,0.],realizedTruth=[0.,0.],
                searchOptions=dict(icaSeeds=2,refineCutoffs=True,eogCheck=True,signalTemplate="data"),
                erp=arm(error),minimal=arm(2),common=arm(3),replayMaxError=0.,exportMaxError=0.,
                alpha=dict(completed=True,recommended=True))
class PilotTests(unittest.TestCase):
    def test_distinct_frozen_short_sources(self):
        self.assertEqual([x["subject"] for x in P["source"]["files"]],list(range(1,11)))
        self.assertEqual(len({x["sha256"] for x in P["source"]["files"]}),10)
        self.assertEqual(P["excerpt"]["duration_seconds"],60)
        self.assertNotIn("balanced",P["evaluation"]["comparators"])
    def test_missing_records_remain_in_coverage(self):
        s=summarize([record()],P,DIGEST)
        self.assertEqual((s["expected"],s["present"],s["completed"]),(10,1,1))
        self.assertFalse(s["executionComplete"]);self.assertFalse(s["pilotUsabilityScreen"])
        self.assertEqual(sum(not x["present"] for x in s["coverage"]),9)
    def test_duplicate_source_is_not_an_extra_record(self):
        with self.assertRaises(AssertionError):summarize([record(),record()],P,DIGEST)
    def test_tampered_mse_is_rejected(self):
        r=record();r["erp"]["mse"]=99
        with self.assertRaises(AssertionError):summarize([r],P,DIGEST)
    def test_changed_options_or_provenance_are_rejected(self):
        for key in ("options","source","protocol","algorithm"):
            r=record()
            if key=="options":r["searchOptions"]["icaSeeds"]=0
            elif key=="source":r["source"]["source_sha256"]="wrong"
            elif key=="protocol":r["protocolSha256"]="wrong"
            else:r["algorithmCommit"]="wrong"
            with self.assertRaises(AssertionError):summarize([r],P,DIGEST)
    def test_only_complete_matched_pairs_enter_comparison(self):
        a,b,c=record(1),record(2,10),record(3,20)
        b["minimal"]["available"]=False;c["completed"]=False
        s=paired_summary([a,b,c],"minimal")
        self.assertEqual(s["paired"],1);self.assertEqual(s["meanDifference"],-3)
        self.assertEqual(s["verdict"],"not-established")
    def test_ten_cases_never_establish_stable_calibration(self):
        s=calibration([record(i) for i in range(1,11)],P)
        self.assertEqual(s["assessable"],10);self.assertFalse(s["meetsPilotPrecision"])
        self.assertFalse(s["stableCalibrationEstablished"])
    def test_undefined_bootstrap_mass_widens_calibration_interval(self):
        a,b=record(1,0),record(2,1)
        a["erp"]["reportedMse"]=0
        s=calibration([a,b],P)
        self.assertEqual(s["ci95"][0],0);self.assertTrue(math.isinf(s["ci95"][1]))
        self.assertGreater(s["undefinedBootstrapShare"],0)
    def test_truncated_edf_is_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            f=Path(d)/"bad.edf";f.write_bytes(b"0"*100)
            with self.assertRaises(ValueError):read_edf(f)
if __name__=="__main__":unittest.main()
