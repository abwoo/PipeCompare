#!/usr/bin/env python3
"""Descriptive scientific ICA readout across saved development artifacts.

Never selects a new recipe or changes the original P3 study. Works on
successful per-subject scientific validation JSON files already in GitHub.
"""
import argparse
import json
import math
from pathlib import Path
from statistics import median

COHORT={11,15,21}


def finite(x):
    return isinstance(x,(int,float)) and math.isfinite(x)


def summarize(files):
    if len(files)!=3:
        raise ValueError(f"expected 3 separate subject science results, got {len(files)}")
    by_subject={}
    for path in files:
        x=json.loads(path.read_text())
        sid=x["subject"]
        if sid not in COHORT or sid in by_subject:
            raise ValueError(f"duplicate/unplanned participant {sid}")
        if x["identityErrors"] or x["scienceFailures"] or not x["unalteredFrozenRankAndEstimates"]:
            raise ValueError(f"subject {sid} failed source parity/science constraints")
        choices=[]
        for cand in x["predeclaredCandidateScience"]:
            starts=cand["startReports"]
            fits=[s for s in starts if s["scienceValid"]]
            valid_boot=[s["bootstrapSE"] for s in fits if finite(s.get("bootstrapSE"))]
            if len(starts)!=6 or len(fits)!=6 or len(valid_boot)!=6:
                raise ValueError(f"missing complete six-start scientific evidence for {sid}/{cand['candidateId']}")
            matches=cand["componentMatchesVsSeed0"]
            cor=[a["medianAbsoluteTopographicCorrelation"] for a in matches
                 if a.get("available") and finite(a.get("medianAbsoluteTopographicCorrelation"))]
            rem=[a["removedComponentMatchedJaccard"] for a in matches
                 if a.get("available") and finite(a.get("removedComponentMatchedJaccard"))]
            wave=[a["fullWaveformRmsDifference"] for a in matches
                  if a.get("available") and finite(a.get("fullWaveformRmsDifference"))]
            p3=[s["P3"] for s in fits]
            choice={"candidateId":cand["candidateId"],"seedCount":6,
                    "meanBootstrapSE":sum(valid_boot)/len(valid_boot),
                    "minBootstrapSE":min(valid_boot),"maxBootstrapSE":max(valid_boot),
                    "P300StartRange":max(p3)-min(p3),
                    "medianAbsTopographicCorr":median(cor) if cor else None,
                    "minimumAbsTopographicCorr":min(cor) if cor else None,
                    "medianRemovedComponentJaccard":median(rem) if rem else None,
                    "medianFullWaveformRmsDiff":median(wave) if wave else None,
                    "topographyComparisonsAvailable":len(cor),
                    "removedICComparisonsAvailable":len(rem),
                    "waveformComparisonsAvailable":len(wave)}
            choices.append(choice)
        by_subject[sid]={"subject":sid,"choices":choices}
    if set(by_subject)!=COHORT:
        raise ValueError("Missing frozen development subjects")
    return {"analysis":"descriptive ICA science discovery, six existing development fits",
            "sourceRun":43184910782,
            "originalSixSeedControlRun":43169562528,
            "originalPR51Verdict":"PARTIAL",
            "validOriginalMatchedPairs":15,
            "newSubjects":False,
            "changeRanking":False,
            "holdoutFrozenExperimentAuthorized":False,
            "participants":[by_subject[x] for x in sorted(COHORT)]}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--inputs",type=Path,nargs="+",required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    result=summarize(a.inputs)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2))
    for s in result["participants"]:
        for c in s["choices"]:
            print("PR51_ICA_SCIENCE_DISCOVERY",json.dumps({
                "subject":s["subject"],**c}))
    print("PR51_ICA_SCIENCE_DECISION",
          "development-only; no accuracy advantage or clinical validity established;" +
          " do not merge PR51; do not open untouched 14-subject holdout")


if __name__=="__main__":
    main()
