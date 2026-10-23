#!/usr/bin/env python3
"""Recover derived CSV metrics only; do not refit or alter any pipeline."""
import csv
import json
import math
import random
import statistics
from pathlib import Path

root=Path("recovered")
out=Path("recovered_summary")
out.mkdir(exist_ok=True)
def rows(name):
    result=[]
    for folder in sorted(root.glob("public-validation-[0-5]")):
        file=folder/(name+".csv")
        if file.exists():
            with file.open(encoding="utf-8-sig",newline="") as stream:
                result.extend(csv.DictReader(stream))
    return result
def number(row,key):
    try:return float(row[key])
    except (ValueError,TypeError,KeyError):return float("nan")
def avg(values):
    v=[x for x in values if math.isfinite(x)]
    return statistics.fmean(v) if v else None
def interval(values,seed=20261023):
    v=[x for x in values if math.isfinite(x)]
    if len(v)<2:return [None,None]
    r=random.Random(seed)
    q=sorted(statistics.fmean(r.choices(v,k=len(v))) for _ in range(1999))
    return [q[49],q[1949]]
def calibration(data,reported="reportedError"):
    blocks={}
    for row in data:blocks.setdefault(row["recording"],[]).append(row)
    a=[];p=[]
    for values in blocks.values():
        x=[number(v,"actualError") for v in values]
        y=[number(v,reported) for v in values]
        if not all(math.isfinite(v) for v in x+y) or any(v<0 for v in y):continue
        a.append(statistics.fmean(v*v for v in x))
        p.append(statistics.fmean(v*v for v in y))
    result={"assessableRecords":len(a),"rows":len(data),"ratio":None}
    if not a or statistics.fmean(p)<=0:return result
    result.update(actualRMS=math.sqrt(statistics.fmean(a)),reportedRMS=math.sqrt(statistics.fmean(p)))
    result["ratio"]=result["actualRMS"]/result["reportedRMS"]
    r=random.Random(7407);q=[]
    for _ in range(1999):
        indices=r.choices(range(len(a)),k=len(a))
        denominator=statistics.fmean(p[i] for i in indices)
        if denominator>0:q.append(math.sqrt(statistics.fmean(a[i] for i in indices)/denominator))
    q.sort()
    if len(q)>=20:result["independentPythonRecordBootstrap95"]=[q[int(.025*len(q))],q[min(len(q)-1,int(.975*len(q)))]]
    return result
def accuracy(data,key):
    groups={}
    for row in data:
        x=number(row,key)
        if math.isfinite(x):groups.setdefault(row["recording"],[]).append(x)
    mse=[statistics.fmean(x*x for x in values) for values in groups.values()]
    bias=[statistics.fmean(values) for values in groups.values()]
    return {"records":len(groups),"rmse":math.sqrt(statistics.fmean(mse)) if mse else None,"bias":avg(bias)}
blink=rows("blink");spectral=rows("spectral");stretches=rows("stretches");erp=rows("real_erp");C=rows("calibration")
selected=[r for r in C if number(r,"selected")==1]
result={"sourceRun":39938907558,"scope":"Recovered CSV metrics; source score-checkpoint errors prevent group/replay validation. No threshold or pipeline was refitted or tuned.",
    "calibration":calibration(C),"selectedCalibration":calibration(selected),
    "selectedAccuracy":{key:accuracy(selected,key) for key in ["actualError","rawError","conditionalError","rawConditionalError"]}}
result["blink"]={"records":len(blink),"annotations":sum(number(r,"annotations") for r in blink),
    "macroPrecision":avg([number(r,"precision") for r in blink]),"macroRecall":avg([number(r,"recall") for r in blink]),
    "macroF1":avg([number(r,"f1") for r in blink]),"macroF1RecordBootstrap95":interval([number(r,"f1") for r in blink]),
    "f1Range":[min(number(r,"f1") for r in blink),max(number(r,"f1") for r in blink)],
    "unassessed":sum(number(r,"assessed")!=1 for r in blink),
    "meanFalseBlinkPredictionsPerMinute":avg([number(r,"falseEventsPerMinute") for r in blink])}
result["spectral"]={"participants":len({r["subject"] for r in spectral}),"rows":len(spectral),
    "maxMatchedReferenceAbsError":max(number(r,"referenceMaxError") for r in spectral),
    "unresolvedRows":sum(number(r,"resolved")!=1 for r in spectral)}
cases=[]
for duration in [1,5,10,20]:
    for adjacent in [0,1]:
        D=[r for r in stretches if number(r,"faultSeconds")==duration and number(r,"adjacent")==adjacent]
        cases.append({"seconds":duration,"adjacent":bool(adjacent),"participants":len(D),
            "assessed":sum(number(r,"assessed")==1 for r in D),
            "sampleRecallIncludingUnassessed":avg([number(r,"sampleRecall") for r in D]),
            "falsePositiveShare":avg([number(r,"falsePositiveShare") for r in D]),
            "corruptedRMSE":avg([number(r,"corruptedRMSE") for r in D]),
            "oracleRepairRMSE":avg([number(r,"oracleRepairRMSE") for r in D]),
            "maximumBoundaryJumpUv":max(number(r,"boundaryJumpUv") for r in D),
            "meanBiasUv":avg([number(r,"meanBiasUv") for r in D]),
            "meanAlphaChangeLog10PSD":avg([number(r,"alphaChangeLog10PSD") for r in D])})
result["stretches"]={"participants":len({r["subject"] for r in stretches}),"cases":cases}
result["realERP"]={"participants":len({r["recording"] for r in erp}),
    "selectedParticipants":len({r["recording"] for r in erp if number(r,"selected")==1}),
    "meanSelectedRawRareMinusFrequentUv":avg([number(r,"rawRareMinusFrequentUv") for r in erp if number(r,"selected")==1])}
prevalence=[]
for p in [0,.3,.7,.95,1]:
    people=[n for n in range(21,51) if [0,.3,.7,.95,1][(n-21)%5]==p]
    D=[r for r in selected if int(float(r["recording"])) in people]
    prevalence.append({"prevalence":p,"planned":len(people),"selected":calibration(D)})
result["prevalence"]=prevalence
sensitivity=[]
for row in selected:
    d=row.copy()
    noise=number(row,"sme")/number(row,"gain")
    association=(abs(number(row,"eyeBias"))+abs(number(row,"rawEstimate"))*number(row,"shapeSpread"))/number(row,"gain")
    d["unshrunkAssociationReportedError"]=str(math.hypot(noise,association))
    sensitivity.append(d)
result["exploratoryReadoutOnlySensitivity"]={"unshrunkAssociation":calibration(sensitivity,"unshrunkAssociationReportedError"),
    "scope":"Post-hoc diagnostic at unchanged selected pipelines; association can include real neural activity and is not causal artifact bias. Not a validated replacement scoring rule."}
assert len(blink)==15 and len({r["recording"] for r in blink})==15
assert len({r["subject"] for r in spectral})==30 and len(spectral)==810
assert len(stretches)==240 and len({r["recording"] for r in erp})==6
assert abs(result["selectedCalibration"]["ratio"]-3.70)<.01, "Independent implementation differs from the MATLAB pilot summary"
(out/"recovered_metrics.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
print("PIPECOMPARE_RECOVERED_JSON "+json.dumps(result,allow_nan=False))
