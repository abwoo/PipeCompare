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
best=R.best;others=find(feasible & (1:numel(cands))'~=best);critical=0;
if ~isempty(others)
    at=others(1);critical=max(0,R.table.diffHi(at)-(R.table.objective(at)-R.table.objective(best)));
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
    'originalKey',cands(proposed).key,'selectedKey',cands(pick).key, ...
    'sourceRun',40282756637,'scoringRun',40430435867,'fullyHeldOut',false,'preprocessingRefitted',false);
fprintf('REFERENCE_REPLAY %s\n',jsonencode(report));
end

function x=mse(cands,k,truth)
o=cands(k).m.objectives(1);g=pipecompare.utils.fieldOr(cands(k).signal,'gain',1);
x=mean((o.estimate/g(1)-truth).^2);
end
