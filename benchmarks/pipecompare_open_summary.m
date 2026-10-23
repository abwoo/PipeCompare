function pipecompare_open_summary(root,outDir)
% Record-balanced summaries and training-only group priors. No EEG needed.
if ~isfolder(outDir),mkdir(outDir);end
folders=dir(fullfile(root,'public-validation-*')); S=struct([]); C=table();
metrics={'blink','spectral','stretches','recordings','real_erp'};
for k=1:numel(metrics)
    combined=table();
    for j=1:numel(folders)
        path=fullfile(folders(j).folder,folders(j).name,[metrics{k} '.csv']);
        if ~isfile(path),continue;end
        T=readtable(path,'TextType','string');
        if isempty(combined),combined=T;else,combined=[combined;T];end %#ok<AGROW>
    end
    if ~isempty(combined),writetable(combined,fullfile(outDir,[metrics{k} '.csv']));end
end
for j=1:numel(folders)
    folder=fullfile(folders(j).folder,folders(j).name);
    if isfile(fullfile(folder,'derived_scores.mat'))
        D=load(fullfile(folder,'derived_scores.mat'));S=[S D.S]; %#ok<AGROW>
    end
    if isfile(fullfile(folder,'calibration.csv'))
        D=readtable(fullfile(folder,'calibration.csv'),'TextType','char');
        if isempty(C),C=D;else,C=[C;D];end %#ok<AGROW>
    end
end
result=struct('commit',getenv('GITHUB_SHA'),'recordsWithScores',numel(S), ...
    'plannedRecords',30,'scope','Independent real-background injection; population injected signal, not certified clean human neural truth');
if ~isempty(C)
    writetable(C,fullfile(outDir,'calibration.csv'));
    result.all=pipecompare.eval.Calibration.summarize(C);
    result.selected=pipecompare.eval.Calibration.summarize(C(C.selected,:));
    fprintf('%s\n',pipecompare.eval.Calibration.describe(C,'all candidates'));
    fprintf('%s\n',pipecompare.eval.Calibration.describe(C(C.selected,:),'selected candidates'));
    % Raw errors are separate accuracy targets, not a relabeling of the
    % gain-corrected reported uncertainty as a calibrated raw-data SE.
    R=C;R.actualError=C.rawError;result.rawAccuracy=pipecompare.eval.Calibration.summarize(R);
    R=C;R.actualError=C.conditionalError;result.realizedTarget=pipecompare.eval.Calibration.summarize(R);
end
rows={};
if numel(S)>=3
    for held=1:numel(S)
        train=setdiff(1:numel(S),held);keys=S(train(1)).keys;
        % Family availability is determined from TRAINING records only.
        % A new person's adaptive extra candidates never tune that family.
        for j=train(2:end),keys=intersect(keys,S(j).keys,'stable');end
        [present,loc]=ismember(keys,S(held).keys);keep=find(present);
        if isempty(keep),continue;end
        [Y,V,covariance]=groupInputs(S(train),keys);
        [~,~,~,prior]=pipecompare.Group.twoLevel(Y,V,covariance);
        [y,v,c]=groupInputs(S(held),keys);
        [theta,info]=pipecompare.Group.applyPrior(prior,y,v,c);
        valid=find(isfinite(y)&y>0);if isempty(valid),continue;end
        [~,a]=min(y(valid));a=valid(a);[~,b]=min(theta(valid));b=valid(b);
        T=S(held).table;[~,f]=min(T.shrunk(loc(valid)));f=valid(f);
        errors=nan(1,numel(keys));
        for p=valid
            name=S(held).labels{loc(p)};
            subject=sscanf(S(held).name,'S%dR01.mat');
            at=C.recording==subject & strcmp(C.steps,name);
            if nnz(at)==2,errors(p)=mean(C.actualError(at).^2);end
        end
        oracle=min(errors,[],'omitnan');
        rows(end+1,:)={S(held).name,numel(train),numel(valid),info.assessed, ...
            errors(a),errors(f),errors(b),oracle,errors(b)-errors(a),a,b}; %#ok<AGROW>
    end
end
if ~isempty(rows)
    T=cell2table(rows,'VariableNames',{'recording','trainingParticipants','eligibleCandidates','priorAssessed', ...
        'rawMinimumMSE','factorialMinimumMSE','frozenGroupMSE','oracleMSE','groupMinusIndependentMSE','rawChoice','groupChoice'});
    writetable(T,fullfile(outDir,'participant_held_out.csv'));
    result.groupTransferRecords=height(T);result.groupTransferAssessed=nnz(T.priorAssessed);
    result.meanGroupMinusIndependentMSE=mean(T.groupMinusIndependentMSE,'omitnan');
end
fid=fopen(fullfile(outDir,'summary.json'),'w');fprintf(fid,'%s\n',jsonencode(result));fclose(fid);
disp(result);
assert(numel(S)==30,'PipeCompare:Validation','Only %d of 30 planned real-background records completed. Inspect shard failures.',numel(S));
end

function [Y,V,C]=groupInputs(S,keys)
n=numel(S);p=numel(keys);Y=nan(n,p);V=nan(n,p);C=cell(n,1);
for j=1:n
    [present,loc]=ismember(keys,S(j).keys);
    valid=false(1,p);valid(present)=strcmp(S(j).table.status(loc(present)),'feasible');
    C{j}=nan(p);if ~any(valid),continue;end
    Y(j,valid)=S(j).table.objective(loc(valid));
    options=S(j).ranking.options;
    options.seed=pipecompare.Group.recordingSeeds({S(j).name},options.seed);
    boot=pipecompare.eval.Rank.bootstrap(S(j).cands(loc(valid)),S(j).ref,options,{S(j).ranking.objective});
    C{j}(valid,valid)=cov(boot{1}');V(j,valid)=diag(C{j}(valid,valid));
end
end
