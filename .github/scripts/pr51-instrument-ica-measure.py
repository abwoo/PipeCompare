#!/usr/bin/env python3
"""Temporarily enrich cloud-run ICA seed comparison with actual P3 readouts.

This edits only the disposable GitHub Actions runner checkout of Executor.m
AFTER the workflow verifies exactly the three allowed production source
deltas. There is no committed change to Executor.m and no altered selection,
gate, seed, trial inclusion or ranking. Outputs remain exploratory evidence.
"""
from __future__ import annotations
from pathlib import Path
import argparse
import subprocess

START = "obj = [T.objective(ids) nan(numel(ids), n)]; rec = [k nan(1, n)];"
PERSEED = "R2 = pipecompare.eval.Rank.run(c2, ref, ro, levelsOf(env.leaves));"
ENDING = "out = struct('seeds', 1:n, 'ids', ids(:)', 'objective', obj, 'recommended', rec, 'same', all(rec == k));"
BEFORE_HELPER = "end\n\nfunction c = rerunWithIca"

MEASURE = r"""
function a = pr51SeedMeasure(cands, R, ref, id, seed)
% Read-only diagnostic; never changes ICA runs, candidate objectives, or ranking.
a = struct('seed',seed,'candidateId',NaN,'available',false, ...
 'estimate',[],'gain',NaN,'correctedDifferenceAmplitude',NaN, ...
 'differenceSEAnalytic',NaN,'kept',[],'nICRemoved',NaN, ...
 'hasICA',false,'nICADataPoints',NaN,'eyeBiasAssessed',[], ...
 'scores',{{}},'weights',{{}},'diagnosticError','');
try
 if isempty(id) || ~isnumeric(id) || ~isscalar(id) || ~isfinite(id) || ...
         id < 1 || id > numel(cands) || id ~= round(id)
    a.diagnosticError = 'No valid recommendation for this seed'; return;
 end
 a.candidateId=id;
 m=cands(id).m;
 if isempty(m) || ~isfield(m,'objectives') || isempty(m.objectives)
    a.diagnosticError = 'Recommended candidate measurement missing'; return;
 end
 o=m.objectives(1);
 if numel(o.estimate)~=2 || numel(o.sme)~=2
    a.diagnosticError = 'Expected exactly two P3 conditions'; return;
 end
 gcol=matlab.lang.makeValidName(['gain_' ref.objectives{1}]);
 if ~ismember(gcol,R.table.Properties.VariableNames)
    a.diagnosticError = 'Reference gain column missing'; return;
 end
 gain=R.table.(gcol)(id);
 if ~isfinite(gain) || gain <= 0 || ~all(isfinite(o.estimate))
    a.diagnosticError = 'Non-finite gain or P3 estimate'; return;
 end
 a.estimate=double(o.estimate(:)');
 a.gain=gain;
 a.correctedDifferenceAmplitude=(a.estimate(1)-a.estimate(2))/gain;
 a.differenceSEAnalytic=sqrt(sum(o.sme.^2))/gain;
 a.kept=m.kept;
 if isfield(o,'eyeBiasAssessed'), a.eyeBiasAssessed=o.eyeBiasAssessed; end
 a.nICRemoved=0;
 for j=1:numel(cands(id).steps)
    v=cands(id).steps{j};
    if strcmp(v.type,'ica')
        a.hasICA=true;
        if isfield(v,'icaPoints'), a.nICADataPoints=v.icaPoints; end
    end
    if strcmp(v.type,'icremove') && isfield(v,'icsRemoved')
        a.nICRemoved=a.nICRemoved+v.icsRemoved;
        if isfield(v,'icaPoints'), a.nICADataPoints=v.icaPoints; end
    end
 end
 if isfield(o,'X') && isfield(o,'w') && numel(o.X)==2 && numel(o.w)==2
    a.scores=cellfun(@(v) double(v(:,1)'),o.X,'UniformOutput',false);
    a.weights=cellfun(@(v) double(v(:)'),o.w,'UniformOutput',false);
 end
 a.available=true;
catch ME
 a.diagnosticError=ME.message;
end
end
"""

def patch(src: str) -> str:
    for marker in (START, PERSEED, ENDING, BEFORE_HELPER):
        if src.count(marker) != 1:
            raise ValueError(f"Expected one unique frozen ICA audit marker: {marker[:65]}")
    src = src.replace(START, START + "\nseedMeasures=cell(1,n+1);\nseedMeasures{1}=pr51SeedMeasure(allc,R,ref,k,0);", 1)
    src = src.replace(PERSEED, PERSEED + "\n    seedMeasures{seed+1}=pr51SeedMeasure(c2,R2,ref,R2.recommended,seed);", 1)
    src = src.replace(ENDING, ENDING + "\nout.seedMeasures=seedMeasures;", 1)
    src = src.replace(BEFORE_HELPER, "end\n\n" + MEASURE.strip() + "\n\nfunction c = rerunWithIca", 1)
    return src

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("repo",type=Path)
    args=parser.parse_args()
    path=args.repo / "+pipecompare/+run/Executor.m"
    original=path.read_text()
    replacement=patch(original)
    path.write_text(replacement)
    diff=subprocess.run(["git","diff","--","+pipecompare/+run/Executor.m"],
                        cwd=args.repo,capture_output=True,text=True,check=True).stdout
    if "seedMeasures{seed+1}=pr51SeedMeasure" not in diff or diff.count("@@") < 3:
        raise ValueError("Unexpected temporary diagnostic patch; refuse to run")
    out=args.repo.parent / "result"
    out.mkdir(exist_ok=True)
    (out / "ica_measure_instrumentation.diff").write_text(diff)
    print("PR51 ICA per-seed measurement audit instrumented ONLY on disposable runner checkout")

if __name__=="__main__":
    main()
