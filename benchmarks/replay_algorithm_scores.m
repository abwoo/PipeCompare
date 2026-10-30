function replay_algorithm_scores(root,outDir)
% Re-score saved development candidates only: never execute preprocessing.
if ~isfolder(outDir),mkdir(outDir);end
records={};files=dir(fullfile(root,'**','score.mat'));
for k=1:numel(files)
    f=fullfile(files(k).folder,files(k).name);
    if ~contains(f,[filesep 'balanced' filesep]),continue;end
    d=load(f,'item');S=d.item;
    levels=recoverLevels(S);
    records{end+1}=inspect(S.cands,S.ref,S.ranking.options,levels,S.ranking.byStratum.recommended, ...
        S.name,'public',root,outDir); %#ok<AGROW>
end
files=dir(fullfile(root,'**','manifest.mat'));
for k=1:numel(files)
    f=fullfile(files(k).folder,files(k).name);
    if ~endsWith(files(k).folder,[filesep 'search']),continue;end
    M=load(f,'result','identity');r=M.result;
    leaves=dir(fullfile(files(k).folder,'leaf_*.mat'));cands=[];kept=[];
    for j=1:numel(leaves)
        d=load(fullfile(leaves(j).folder,leaves(j).name),'cand','id','key','identity');
        assert(strcmp(d.identity,M.identity)&&strcmp(d.key,r.leaves(d.id).key),'Checkpoint identity mismatch');
        if isempty(cands),cands=d.cand;else,cands(end+1)=d.cand;end %#ok<AGROW>
        kept(end+1)=d.id; %#ok<AGROW>
    end
    [kept,order]=sort(kept);cands=cands(order);levels=cell(numel(cands),1);
    for j=1:numel(cands)
        p=r.leaves(kept(j)).path;
        levels{j}=[cellfun(@(v)v.slot,p(:),'UniformOutput',false),cellfun(@(v)v.key,p(:),'UniformOutput',false)];
    end
    token=regexp(f,'record_(\d+)','tokens','once');
    assert(~isempty(token),'Missing recording identity');
    records{end+1}=inspect(cands,r.ref,r.options,levels,[],token{1},'synthetic',root,outDir); %#ok<AGROW>
end
assert(numel(records)==7,'Expected four public summaries and three cached simulation trees');
save(fullfile(outDir,'replay.mat'),'records','-v7');
fprintf('SAVED_SCORE_REPLAY_RECORDS=%d; no preprocessing or ICA executed\n',numel(records));
end

function L=recoverLevels(S)
if isfield(S,'levels'),L=S.levels;return;end
% Historical archives omit factor levels. For these recipes factor names
% are step types; preserve the saved model's exclusion of refined leaves.
L=cell(numel(S.cands),1);feas=find(strcmp(S.table.status,'feasible'));
F=S.ranking.byStratum.factorial;
if isempty(F),return;end
for i=feas(F.use(:))'
    parts=strsplit(S.labels{i},' > ');slot=cell(size(parts));
    for q=1:numel(parts),slot{q}=regexp(parts{q},'^[^(]+','match','once');end
    L{i}=[slot(:),parts(:)];
end
end

function report=inspect(cands,ref,opts,levels,original,name,population,root,outDir)
opts.nSplits=0;opts.upgradeBaseline=[];opts.nBoot=1999;
opts.protectReferenceCleanup=false; % isolate variance-model diagnostics from the later decision guard
for j=1:numel(cands)
    if ~isfield(cands,'message')||isempty(cands(j).message),cands(j).message='saved candidate not eligible';end
end
calib=readtable(fullfile(root,['development-optimization-' population '-summary'],'balanced_calibration.csv'), ...
    'TextType','string','VariableNamingRule','preserve');
id=str2double(regexp(name,'\d+','match','once'));C=calib(calib.recording==id,:);
assert(~isempty(C),'No matching frozen calibration record');
models={'pooled','heteroscedastic'};R=cell(1,2);cc=cell(1,2);
for mode=1:2
    cc{mode}=cands;
    for j=1:numel(cands)
        if ~strcmp(cands(j).status,'ok'),continue;end
        obs=cands(j).m.objectives;
        for o=1:numel(obs)
            obs(o).eyeBiasSEMode=models{mode};
            if ~strcmp(obs(o).kind,'scalar')||~isfield(ref,'blink')||numel(ref.blink)<o||isempty(ref.blink{o}),continue;end
            for q=1:numel(ref.ids)
                x=obs(o).X{q}(:,1);w=ones(size(x));
                if isfield(obs,'w')&&~isempty(obs(o).w),w=obs(o).w{q};w(~isfinite(w))=0;end
                known=~isfield(obs,'eyeBiasAssessed')||obs(o).eyeBiasAssessed(q);
                [a,b,ok,se]=pipecompare.eval.Measure.eyeBias(x,ref.blink{o}{q},w,ones(size(w)),models{mode});
                if ~known,a=NaN;b=NaN;ok=false;se=NaN;end
                obs(o).eyeBias(q)=a;obs(o).eyeBiasSq(q)=b;
                obs(o).eyeBiasAssessed(q)=ok;obs(o).eyeBiasSE(q)=se;
            end
        end
        cc{mode}(j).m.objectives=obs;
    end
    R{mode}=pipecompare.eval.Rank.run(cc{mode},ref,opts,levels);
end
if isempty(original),original=R{1}.recommended;scope='cached main-tree subset; later adaptive leaves absent';
else,scope='completed candidate family; fitted choices fixed';end
truth=zeros(1,numel(ref.ids));
for q=1:numel(truth),v=unique(C.populationTruth(C.condition==q));assert(isscalar(v));truth(q)=v;end
old=metrics(cands,original,truth);base=metrics(cc{1},R{1}.recommended,truth);new=metrics(cc{2},R{2}.recommended,truth);
guard=[];protected=[];
if ~isempty(original)&&strcmp(R{2}.table.status{original},'feasible')&&~isempty(R{2}.recommended)
    feasible=find(strcmp(R{2}.table.status,'feasible'));best=R{2}.best;other=feasible(feasible~=best);
    critical=0;
    if ~isempty(other)
        at=other(1);critical=max(0,R{2}.table.diffHi(at)-(R{2}.table.objective(at)-R{2}.table.objective(best)));
    end
    [pick,guard]=pipecompare.eval.Rank.guardUpgrade(R{2}.table.objective,R{2}.recommended,original,critical);
    protected=metrics(cc{2},pick,truth);
end
report=struct('record',id,'population',population,'scope',scope,'original',old,'pooledReplay',base, ...
    'unequalVariance',new,'protected',protected,'guard',guard, ...
    'sourceRun',40282756637,'sourceRevision','4116ea11f3d5d39e4d3d282894dcf487172bd50d', ...
    'fullyHeldOut',false,'preprocessingRefitted',false);
fprintf('SCORE_REPLAY %s\n',jsonencode(report));
save(fullfile(outDir,sprintf('%s_%d.mat',population,id)),'report','R','-v7');
end

function v=metrics(cands,k,truth)
v=struct('id',[],'mse',NaN,'key','');
if isempty(k),return;end
o=cands(k).m.objectives(1);g=pipecompare.utils.fieldOr(cands(k).signal,'gain',1);
v=struct('id',k,'mse',mean((o.estimate/g(1)-truth).^2),'key',cands(k).key);
end
