function replay_reference_guard(root,previous,outDir)
% Apply a decision gate to saved bootstrap bounds and shrunk scores.
% Never refit EEG, ICA, the ranking bootstrap or its covariance model.
if ~isfolder(outDir),mkdir(outDir);end
reports={};files=dir(fullfile(root,'**','score.mat'));
for k=1:numel(files)
    f=fullfile(files(k).folder,files(k).name);
    if ~contains(f,[filesep 'balanced' filesep]),continue;end
    d=load(f,'item');id=str2double(regexp(d.item.name,'\d+','match','once'));
    reports{end+1}=apply(d.item.cands,id,'public',root,previous); %#ok<AGROW>
end
files=dir(fullfile(root,'**','manifest.mat'));
for k=1:numel(files)
    if ~endsWith(files(k).folder,[filesep 'search']),continue;end
    M=load(fullfile(files(k).folder,files(k).name),'identity','result');
    leaves=dir(fullfile(files(k).folder,'leaf_*.mat'));cands=[];ids=[];
    for j=1:numel(leaves)
        d=load(fullfile(leaves(j).folder,leaves(j).name),'cand','identity','id','key');
        assert(strcmp(d.identity,M.identity)&&strcmp(d.key,M.result.leaves(d.id).key));
        if isempty(cands),cands=d.cand;else,cands(end+1)=d.cand;end %#ok<AGROW>
        ids(end+1)=d.id; %#ok<AGROW>
    end
    [~,order]=sort(ids);cands=cands(order);
    token=regexp(files(k).folder,'record_(\d+)','tokens','once');assert(~isempty(token));
    reports{end+1}=apply(cands,str2double(token{1}),'synthetic',root,previous); %#ok<AGROW>
end
assert(numel(reports)==7);
save(fullfile(outDir,'reference_guard_replay.mat'),'reports','-v7');
fid=fopen(fullfile(outDir,'reference_guard_replay.json'),'w');cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(reports,PrettyPrint=true));
fprintf('REFERENCE_GUARD_RECORDS=%d; reused fitted candidates, bootstrap bounds and shrinkage\n',numel(reports));
end

function report=apply(cands,id,population,root,previous)
d=load(fullfile(previous,sprintf('%s_%d.mat',population,id)),'R','report');R=d.R{1};old=d.report;
assert(height(R.table)==numel(cands)&&~isempty(R.recommended));
[repaired,needed]=pipecompare.eval.Rank.referenceProtection(cands);
feasible=strcmp(R.table.status,'feasible');repaired=repaired&feasible;
best=R.best;
% Production critical value, recovered from the saved bestSet intervals.
% Rank.bestSet (+pipecompare/+eval/Rank.m, bestSet) sets hi = Dhat(:,ib) + c,
% so c = diffHi - (objective - objective(best)) for every feasible
% non-best candidate; the saved table already stores bestSet's output, and
% the replay deliberately does not refit the bootstrap, so the production
% critical is recovered here rather than recomputed from replicates
% (bestSet itself is a file-local function in Rank.m and cannot be called
% from this script). All feasible non-best candidates must agree:
% disagreement means the saved artifact is not a bestSet table. (The best
% candidate itself is excluded: bestSet writes hi(ib) = 0, not Dhat + c.)
% An unbounded critical (Inf, from incomplete bootstrap draws) is kept as
% Inf, which the guard records as 'unassessed' rather than as evidence
% for either side.
others=find(feasible & (1:numel(cands))'~=best);
if isempty(others)
    critical=0;
else
    recovered=R.table.diffHi(others)-(R.table.objective(others)-R.table.objective(best));
    if all(isinf(recovered))
        critical=Inf;
    else
        finiteRecovered=recovered(isfinite(recovered));
        assert(numel(finiteRecovered)==numel(recovered) && ...
            max(finiteRecovered)-min(finiteRecovered)<1e-9*max(1,max(abs(finiteRecovered))), ...
            'Saved diffHi does not encode a single bestSet critical value');
        critical=max(0,recovered(1));
    end
end
proposed=R.recommended;
[pick,guard]=pipecompare.eval.Rank.guardReference(R.table.objective,R.table.shrunk, ...
    R.table.bias>pipecompare.eval.Rank.BiasShare*R.table.noise,repaired,needed,proposed,critical);
% Truth is read AFTER the operational choice; it never enters the gate.
C=readtable(fullfile(root,['development-optimization-' population '-summary'],'balanced_calibration.csv'), ...
    'TextType','string','VariableNamingRule','preserve');
C=C(C.recording==id,:);assert(~isempty(C));n=numel(cands(proposed).m.objectives(1).estimate);truth=zeros(1,n);
for q=1:n,v=unique(C.populationTruth(C.condition==q));assert(isscalar(v));truth(q)=v;end
oldMse=mse(cands,proposed,truth);newMse=mse(cands,pick,truth);
assert(abs(oldMse-old.pooledReplay.mse)<1e-9,'Candidate ordering does not reproduce the source artifact');
report=struct('record',id,'population',population,'scope',old.scope,'originalMSE',oldMse, ...
    'guardedMSE',newMse,'proposed',proposed,'selected',pick,'guard',guard, ...
    'productionGuardCritical',critical, ...
    'originalKey',cands(proposed).key,'selectedKey',cands(pick).key, ...
    'sourceRun',40282756637,'scoringRun',40430435867,'fullyHeldOut',false,'preprocessingRefitted',false);
fprintf('REFERENCE_REPLAY %s\n',jsonencode(report));
end

function x=mse(cands,k,truth)
o=cands(k).m.objectives(1);g=pipecompare.utils.fieldOr(cands(k).signal,'gain',1);
x=mean((o.estimate/g(1)-truth).^2);
end
