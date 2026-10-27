"""Inspect existing development error components. Never refit EEG."""
import csv,json,math
from pathlib import Path
from statistics import mean

def blocks(rows):
    result=[]
    for row in rows:
        if not result or row["recording"]!=result[-1][0]["recording"] or int(row["condition"])<=int(result[-1][-1]["condition"]):
            result.append([])
        result[-1].append(row)
    return result
def feature(b):
    score=[]
    for r in b:
        g=float(r["gain"]);e=float(r["eyeBias"])/g;s=float(r["artifactBiasSE"])
        shape=abs(float(r["correctedEstimate"]))*float(r["shapeSpread"])
        score.append((float(r["sme"])/g)**2+(math.sqrt(max(e*e-s*s,0))+shape)**2)
    return dict(label=b[0]["steps"],mse=mean(float(r["actualError"])**2 for r in b),
        old_risk=mean(float(r["reportedError"])**2 for r in b),debiased_risk=mean(score),
        selected=b[0]["selected"].lower() in ("1","true"),
        components=[{k:float(r[k]) for k in ("actualError","sme","eyeBias","gain","artifactBiasSE","totalBias","shapeSpread","rawEstimate","populationTruth","reportedError")} for r in b])
totals={}
for population in ("public","synthetic"):
    src=Path("source_results")/("development-optimization-"+population+"-summary")
    values=[]
    for recipe in ("standard","balanced"):
        with (src/(recipe+"_calibration.csv")).open() as f:bs=blocks(list(csv.DictReader(f)))
        groups={}
        for b in bs:groups.setdefault(b[0]["recording"],[]).append(feature(b))
        for id,cs in groups.items():
            chosen=next(c for c in cs if c["selected"])
            risk=min(cs,key=lambda c:c["debiased_risk"])
            old=min(cs,key=lambda c:c["old_risk"])
            print("COMPONENTS "+json.dumps(dict(population=population,recipe=recipe,recording=id,selected=chosen,debiased_min=risk,old_min=old),allow_nan=False))
            if recipe=="balanced":values.append(risk["mse"])
    totals[population]=mean(values)
print("DEBIASED_MIN_TOTAL "+json.dumps(totals))
