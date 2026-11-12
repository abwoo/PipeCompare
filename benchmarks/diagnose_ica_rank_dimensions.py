#!/usr/bin/env python3
"""Read-only ICA component-rank dimension diagnosis on archived EEG science."""
import argparse
import json
from pathlib import Path
from validate_ica_unseen_seeds import grid

SUBJECTS={11,15,21}
PREDECLARED={11:[12,25],15:[5,2,6],21:[25,4]}


def inspect(study):
    c=study["icaSeedCheck"];n=len(c["ids"])
    m=grid(c["candidateMeasurements"],n,"candidateMeasurements")
    status=grid(c["rankStatus"],n,"rankStatus")
    details=[]
    for ident in PREDECLARED[study["subject"]]:
        i=c["ids"].index(ident)
        rows=[]
        for seed in range(6):
            z=m[i][seed]
            science=z.get("science",{})
            decomp=science.get("decompositions",[])
            matrices=decomp if isinstance(decomp,list) else [decomp]
            dims=[d.get("nComponents") for d in matrices]
            channels=[d.get("nChannels") for d in matrices]
            rows.append({"seed":seed,"rankStatus":status[i][seed],
                         "recordStatus":science.get("status"),
                         "components":dims,"channels":channels,
                         "removedCount":sum(len(d["removed"] if isinstance(d["removed"],list)
                                                 else [d["removed"]]) for d in matrices)})
        distinct={tuple(r["components"]) for r in rows}
        details.append({"candidateId":ident,"perSeed":rows,
                        "componentDimensionsStable":len(distinct)==1,
                        "observedComponentDimensions":[list(v) for v in sorted(distinct)]})
    return {"subject":study["subject"],"choices":details,
            "unchangedRank":True,"noNewEEGFits":True,"scientificPromotionAllowed":False}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--evidence",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    files=sorted(a.evidence.rglob("result.json"))
    if len(files)!=3:raise SystemExit("Expected exactly three preserved source records")
    rows=[inspect(json.loads(f.read_text())) for f in files]
    if {r["subject"] for r in rows}!=SUBJECTS:raise SystemExit("Wrong source roster")
    rows.sort(key=lambda x:x["subject"])
    out={"purpose":"read-only component dimension rank diagnosis",
         "sourceRun":43184910782,"newSubjects":False,
         "originalPR51Verdict":"PARTIAL","rows":rows}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,indent=2))
    for r in rows:
        for ch in r["choices"]:
            print("PR51_ICA_RANK_DIMENSIONS",json.dumps({
                "subject":r["subject"],"candidateId":ch["candidateId"],
                "stable":ch["componentDimensionsStable"],
                "perSeed":[{"seed":x["seed"],"components":x["components"],
                            "channels":x["channels"],"removedCount":x["removedCount"]}
                           for x in ch["perSeed"]]}))


if __name__=="__main__":
    main()
