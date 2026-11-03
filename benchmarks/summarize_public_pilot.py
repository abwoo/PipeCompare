"""Independently recompute the frozen pilot, retaining every planned participant."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import numpy as np
PROTOCOL=Path(__file__).with_name("public_pilot_protocol.json")

def finite(x):
    return isinstance(x,(float,int)) and not isinstance(x,bool) and math.isfinite(x)

def interval(x,qs=(.025,.975)):
    x=np.sort(np.asarray(x,float))
    if not len(x): return [None,None]
    return [float(x[max(0,min(len(x)-1,math.ceil((len(x)+1)*q)-1))]) for q in qs]

def paired_summary(records,arm,draws=10000,seed=510003):
    pairs=[(r["erp"]["mse"],r[arm]["mse"]) for r in records if r.get("completed") and
           r["erp"].get("available") and r[arm].get("available") and
           finite(r["erp"].get("mse")) and finite(r[arm].get("mse"))]
    if not pairs: return dict(paired=0,standardMSE=None,controlMSE=None,meanDifference=None,ci95=[None,None],verdict="not-established")
    a=np.asarray(pairs,float);delta=a[:,0]-a[:,1];ci=[None,None]
    if len(a)>=2:
        idx=np.random.RandomState(seed).randint(0,len(a),(draws,len(a)))
        ci=interval(delta[idx].mean(1))
    verdict="not-established"
    if ci[1] is not None and ci[1]<0: verdict="pilot-improvement"
    elif ci[0] is not None and ci[0]>0: verdict="pilot-worsening"
    return dict(paired=len(a),standardMSE=float(a[:,0].mean()),controlMSE=float(a[:,1].mean()),
                meanDifference=float(delta.mean()),ci95=ci,verdict=verdict)

def calibration(records,p):
    rows=[r["erp"] for r in records if r.get("completed") and r["erp"].get("available") and
          finite(r["erp"].get("mse")) and finite(r["erp"].get("reportedMse")) and r["erp"]["reportedMse"]>=0]
    out=dict(assessable=len(rows),planned=p["evaluation"]["expected_records"],ratio=None,ci95=[None,None],
             actualRMS=None,reportedRMS=None,exploratory=True,stableCalibrationEstablished=False,
             relativeHalfWidth=None,meetsPilotPrecision=False,agreement="not-established")
    if not rows:return out
    a=np.array([r["mse"] for r in rows]);b=np.array([r["reportedMse"] for r in rows])
    out["actualRMS"]=float(np.sqrt(a.mean()));out["reportedRMS"]=float(np.sqrt(b.mean()))
    if b.mean()>0:out["ratio"]=float(np.sqrt(a.mean()/b.mean()))
    elif a.mean()>0:out["ratio"]=float("inf")
    if len(rows)<2:return out
    idx=np.random.RandomState(p["evaluation"]["bootstrap_seed"]).randint(0,len(rows),(p["evaluation"]["bootstrap_draws"],len(rows)))
    aa=a[idx].mean(1);bb=b[idx].mean(1)
    ratio=np.full(len(aa),np.nan);pos=bb>0
    ratio[pos]=np.sqrt(aa[pos]/bb[pos]);ratio[(bb==0)&(aa>0)]=np.inf
    lo=ratio.copy();hi=ratio.copy();lo[np.isnan(lo)]=0;hi[np.isnan(hi)]=np.inf
    out["ci95"]=[interval(lo)[0],interval(hi)[1]]
    out["undefinedBootstrapShare"]=float(np.isnan(ratio).mean())
    out["infiniteBootstrapShare"]=float(np.isinf(ratio).mean())
    if finite(out["ratio"]) and out["ratio"]>0 and all(finite(x) for x in out["ci95"]):
        out["relativeHalfWidth"]=(out["ci95"][1]-out["ci95"][0])/(2*out["ratio"])
        out["meetsPilotPrecision"]=len(rows)>=p["evaluation"]["precision_floor_records"] and out["relativeHalfWidth"]<=p["evaluation"]["relative_half_width_target"]
    lower,upper=p["evaluation"]["calibration_agreement_band"];l,u=out["ci95"]
    if l is not None and u is not None:
        if l>=lower and u<=upper:out["agreement"]="preliminary-agreement"
        elif u<lower or l>upper:out["agreement"]="adverse"
    return out

def validate_record(r,p,digest):
    assert r["algorithmCommit"]==p["algorithm_commit"],"Algorithm revision mismatch"
    assert r["protocolSha256"]==digest,"Protocol mismatch"
    entry=next(x for x in p["source"]["files"] if x["subject"]==r["subject"])
    source=r["source"]
    assert source["subject"]==r["subject"] and source["source_sha256"]==entry["sha256"],"Source identity mismatch"
    assert source["sampling_rate"]==p["excerpt"]["sampling_rate"] and source["samples"]==p["excerpt"]["duration_seconds"]*p["excerpt"]["sampling_rate"],"Excerpt mismatch"
    if not r["completed"]: return
    assert not r.get("errors"),"Completed record contains errors"
    assert r.get("notRun")==0,"Incomplete candidate family"
    opts=r["searchOptions"]
    assert opts["icaSeeds"]==2 and opts["refineCutoffs"] and opts["eogCheck"] and opts["signalTemplate"]=="data","Standard options changed"
    for arm in ("erp","minimal","common"):
        x=r[arm]
        if not x["available"]:continue
        estimate=np.asarray(x["estimate"],float);pop=np.asarray(r["populationTruth"],float);real=np.asarray(r["realizedTruth"],float)
        assert estimate.shape==pop.shape==real.shape==(2,),"Unexpected condition shape"
        corrected=estimate/x["gain"]
        assert np.allclose(corrected,x["correctedEstimate"],atol=1e-9,rtol=1e-9)
        for key,want in (("actualError",corrected-pop),("rawError",estimate-pop),("conditionalError",corrected-real)):
            assert np.allclose(want,x[key],atol=1e-8,rtol=1e-8),"Error field mismatch: "+key
        for key,field in (("mse","actualError"),("rawMse","rawError"),("conditionalMse","conditionalError"),("reportedMse","reportedError")):
            want=float(np.mean(np.asarray(x[field],float)**2))
            assert finite(x[key]) and np.isclose(want,x[key],atol=1e-8,rtol=1e-8),"MSE mismatch: "+key
    if r["erp"]["available"]:
        assert finite(r["replayMaxError"]) and r["replayMaxError"]<p["evaluation"]["replay_tolerance_uv"],"Replay failed"
        assert finite(r["exportMaxError"]) and r["exportMaxError"]<p["evaluation"]["replay_tolerance_uv"],"Export failed"

def summarize(records,p,digest):
    expected=[x["subject"] for x in p["source"]["files"]];by={}
    for r in records:
        assert r["subject"] in expected and r["subject"] not in by,"Unexpected or duplicate participant"
        validate_record(r,p,digest);by[r["subject"]]=r
    revisions={r["executionCommit"] for r in records}
    assert len(revisions)<=1,"Mixed execution revisions"
    complete=[r for r in records if r.get("completed")]
    selected=[r for r in complete if r["erp"]["available"]]
    rows=[]
    for subject in expected:
        r=by.get(subject)
        rows.append(dict(subject=subject,present=r is not None,completed=bool(r and r["completed"]),
                         recommended=bool(r and r["erp"]["available"]),
                         errors=[] if r is None else r.get("errors",[])))
    result=dict(scope=p["scope"],algorithmCommit=p["algorithm_commit"],executionCommits=sorted(revisions),protocolSha256=digest,
                expected=len(expected),present=len(by),completed=len(complete),assessableRecommendations=len(selected),
                executionComplete=len(complete)==len(expected),pilotUsabilityScreen=False,
                independentConfirmation=False,coverage=rows,calibration=calibration(records,p),
                comparisons={arm:paired_summary(records,arm,p["evaluation"]["bootstrap_draws"],p["evaluation"]["bootstrap_seed"]) for arm in ("minimal","common")},
                alphaCompleted=sum(bool(r.get("alpha",{}).get("completed")) for r in records),
                alphaRecommendations=sum(bool(r.get("alpha",{}).get("recommended")) for r in records))
    result["pilotUsabilityScreen"]=result["executionComplete"] and len(selected)>=p["evaluation"]["minimum_recommendations_for_pilot_execution"]
    return result

def safe_json(v):
    if isinstance(v,float) and not math.isfinite(v):return "Infinity" if math.isinf(v) and v>0 else ("-Infinity" if math.isinf(v) else None)
    if isinstance(v,dict):return {k:safe_json(x) for k,x in v.items()}
    if isinstance(v,list):return [safe_json(x) for x in v]
    return v

def main():
    a=argparse.ArgumentParser();a.add_argument("root",type=Path);a.add_argument("out",type=Path);args=a.parse_args()
    p=json.loads(PROTOCOL.read_text());digest=hashlib.sha256(PROTOCOL.read_bytes()).hexdigest()
    records=[json.loads(f.read_text()) for f in sorted(args.root.rglob("result.json"))]
    result=summarize(records,p,digest);args.out.mkdir(parents=True,exist_ok=True)
    (args.out/"summary.json").write_text(json.dumps(safe_json(result),indent=2,allow_nan=False)+"\n")
    with (args.out/"coverage.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=["subject","present","completed","recommended","errors"]);w.writeheader();w.writerows(result["coverage"])
    lines=["# Ten-record public EEG pilot","",p["scope"],"",
           f"Execution: {result['completed']}/{result['expected']}; recommendations: {result['assessableRecommendations']}/{result['expected']}.",
           f"Unaltered alpha searches: {result['alphaCompleted']}/10; recommendations: {result['alphaRecommendations']}/10.",
           "", "| Comparator | Complete pairs | Standard MSE | Control MSE | Difference [95% CI] | Pilot verdict |",
           "|---|---:|---:|---:|---|---|"]
    def fmt(x):return "unavailable" if x is None else f"{x:.6g}"
    for arm,c in result["comparisons"].items():
        lines.append(f"| {arm} | {c['paired']} | {fmt(c['standardMSE'])} | {fmt(c['controlMSE'])} | {fmt(c['meanDifference'])} [{fmt(c['ci95'][0])}, {fmt(c['ci95'][1])}] | {c['verdict']} |")
    c=result["calibration"]
    lines+=["",f"Selected actual/reported RMS ratio: {fmt(c['ratio'])}, 95% CI [{fmt(c['ci95'][0])}, {fmt(c['ci95'][1])}], {c['assessable']} assessable records. Agreement: {c['agreement']}.",
            "", "Ten records are exploratory and cannot establish stable calibration. Natural resting EEG is a background source, not a real P3 task or clean neural ground truth.",
            "The controls and whole Standard procedure were fixed before outcomes. This run does not open the PR50 confirmation cohort.",
            "", "| Subject | Present | Completed | Recommendation | Errors |", "|---|---|---|---|---|"]
    for r in result["coverage"]:
        errs="; ".join(x["stage"]+": "+x["message"] for x in r["errors"])
        lines.append(f"| {r['subject']} | {r['present']} | {r['completed']} | {r['recommended']} | {errs.replace('|','/')} |")
    (args.out/"summary.md").write_text("\n".join(lines)+"\n")
    print("PIPECOMPARE_PILOT_SUMMARY "+json.dumps(safe_json(result),allow_nan=False),flush=True)
    if not result["executionComplete"]:raise SystemExit("Some planned records failed or are missing; inspect coverage and artifacts.")
if __name__=="__main__":main()
