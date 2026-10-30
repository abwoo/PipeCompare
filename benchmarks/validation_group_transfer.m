function result=validation_group_transfer(S,C,outDir,trainingSubjects,testSubjects)
% The participant split comes from frozen source metadata. Intervals are
% conditional on that fitted training prior, not on new training cohorts.
for cohort={trainingSubjects,testSubjects}
    v=cohort{1};
    assert(isnumeric(v)&&isreal(v)&&(isempty(v)||isvector(v))&& ...
        all(isfinite(v(:))&v(:)>0&v(:)==fix(v(:)))&&numel(unique(v))==numel(v), ...
        'PipeCompare:Validation','Cohort IDs must be unique positive integers.');
end
assert(isempty(intersect(trainingSubjects,testSubjects)), ...
    'PipeCompare:Validation','Training and test participants must be disjoint.');
id=zeros(1,numel(S));
for j=1:numel(S)
    token=regexp(S(j).name,'^S([0-9]+)R01\.mat$','tokens','once');
    assert(~isempty(token),'PipeCompare:Validation','Malformed participant recording name.');
    id(j)=str2double(token{1});
end
assert(all(isfinite(id)&id>0&id==fix(id))&&numel(unique(id))==numel(id), ...
    'PipeCompare:Validation','Each participant must have exactly one uniquely identified recording.');
train=find(ismember(id,trainingSubjects));test=find(ismember(id,testSubjects));
train=train(:)';test=test(:)';
result=struct('trainingRequested',numel(trainingSubjects),'testRequested',numel(testSubjects),'trainingObserved',0, ...
    'testCovered',0,'scope','Independent test participants conditional on a frozen training prior');
if isempty(train)||isempty(test),return;end
keys=S(train(1)).keys;
for j=train(2:end),keys=intersect(keys,S(j).keys,'stable');end
if isempty(keys),return;end
[Y,V,Cov]=groupInputs(S(train),keys);
[~,~,~,prior]=pipecompare.Group.twoLevel(Y,V,Cov);
result.trainingWithPointScores=nnz(any(isfinite(Y)&Y>0,2));
result.trainingObserved=nnz(any(isfinite(Y)&Y>0&isfinite(V)&V>=0,2));
rows={};coverageRows={};
for j=test(:)'
    [y,v,covariance]=groupInputs(S(j),keys);
    [theta,info]=pipecompare.Group.applyPrior(prior,y,v,covariance);
    [a,b,valid]=validation_group_choices(y,v,theta);
    [recommended,recommendedContrast,inFamily]=validation_recommended_loss(C,S(j),keys);
    coverageRows(end+1,:)={id(j),nnz(isfinite(y)&y>0), ...
        nnz(isfinite(y)&y>0&isfinite(v)&v>=0),numel(valid),info.assessed, ...
        isfinite(recommended),info.reason{1}}; %#ok<AGROW>
    if isempty(valid),continue;end
    assert(info.assessed||a==b,'PipeCompare:Validation','Unassessed prior changed a matched-support raw choice.');
    [~,loc]=ismember(keys,S(j).keys);
    [raw,rawContrast]=candidateLoss(C,S(j),loc(a));
    [group,groupContrast]=candidateLoss(C,S(j),loc(b));
    rows(end+1,:)={id(j),info.assessed,numel(valid),a,b,raw,group,group-raw, ...
        rawContrast,groupContrast,recommended,recommendedContrast,inFamily,group-recommended}; %#ok<AGROW>
end
if ~isempty(coverageRows)
    coverageTable=cell2table(coverageRows,'VariableNames',{'subject','pointEligibleCandidates', ...
        'finiteCovarianceCandidates','comparableCandidates','priorAssessed','actualRecommendationAvailable','reason'});
    writetable(coverageTable,fullfile(outDir,'frozen_prior_independent_test_coverage.csv'));
    result.testWithPointScores=nnz(coverageTable.pointEligibleCandidates>0);
    result.testNoComparableScore=nnz(coverageTable.comparableCandidates==0);
    result.actualRecommendationsAvailableOverall=nnz(coverageTable.actualRecommendationAvailable);
end
if isempty(rows),return;end
T=cell2table(rows,'VariableNames',{'subject','priorAssessed','candidates','rawChoice','groupChoice', ...
    'rawMinimumMSE','frozenGroupMSE','groupMinusRawMinimumMSE','rawMinimumContrastError','frozenGroupContrastError', ...
    'actualRecommendationMSE','actualRecommendationContrastError','recommendationInTrainingFamily','groupMinusRecommendationMSE'});
writetable(T,fullfile(outDir,'frozen_prior_independent_test.csv'));
result.testCovered=height(T);result.priorAssessed=nnz(T.priorAssessed);
result.meanGroupMinusRawMinimumMSE=mean(T.groupMinusRawMinimumMSE);
result.rawMinimumComparisonConditional95=recordInterval(T.groupMinusRawMinimumMSE);
result.meanGroupMinusRecommendationMSE=mean(T.groupMinusRecommendationMSE,'omitnan');
result.recommendationComparisonConditional95=recordInterval(T.groupMinusRecommendationMSE);
result.actualRecommendationsAvailable=nnz(isfinite(T.actualRecommendationMSE));
result.actualRecommendationsInTrainingFamily=nnz(T.recommendationInTrainingFamily);
result.changedChoicesVersusRawMinimum=nnz(T.rawChoice~=T.groupChoice);
result.harmfulChangesVersusRawMinimum=nnz(T.groupMinusRawMinimumMSE>0&T.rawChoice~=T.groupChoice);
result.harmfulChangesVersusRecommendation=nnz(T.groupMinusRecommendationMSE>0);
result.assessedOnlyMeanGroupMinusRawMinimumMSE=mean(T.groupMinusRawMinimumMSE(T.priorAssessed==1),'omitnan');
result.assessedOnlyMeanGroupMinusRecommendationMSE=mean(T.groupMinusRecommendationMSE(T.priorAssessed==1),'omitnan');
result.comparisonScope='Raw-minimum comparison uses the training-supported family; actual recommendations use the original full per-record family.';
a=mod(T.subject,2)==1;b=~a;
if any(a)&&any(b)
    result.rawMinimumGroupContrastError=mean(T.rawMinimumContrastError(a))-mean(T.rawMinimumContrastError(b));
    result.frozenPriorGroupContrastError=mean(T.frozenGroupContrastError(a))-mean(T.frozenGroupContrastError(b));
    result.actualRecommendationGroupContrastError=mean(T.actualRecommendationContrastError(a),'omitnan')- ...
        mean(T.actualRecommendationContrastError(b),'omitnan');
end
end

function [Y,V,C]=groupInputs(S,keys)
n=numel(S);p=numel(keys);Y=nan(n,p);V=nan(n,p);C=cell(n,1);
for j=1:n
    [present,loc]=ismember(keys,S(j).keys);
    valid=false(1,p);valid(present)=strcmp(S(j).table.status(loc(present)),'feasible');
    C{j}=nan(p);if ~any(valid),continue;end
    Y(j,valid)=S(j).table.objective(loc(valid))';
    options=S(j).ranking.options;
    options.seed=pipecompare.Group.recordingSeeds({S(j).name},options.seed);
    boot=pipecompare.eval.Rank.bootstrap(S(j).cands(loc(valid)),S(j).ref,options,{S(j).ranking.objective});
    C{j}(valid,valid)=cov(boot{1}');V(j,valid)=diag(C{j}(valid,valid))';
end
end

function [loss,contrast]=candidateLoss(C,S,k)
subject=sscanf(S.name,'S%dR01.mat');o=S.cands(k).m.objectives(1);
names=S.table.Properties.VariableNames;gains=names(startsWith(names,'gain_'));
gain=S.table.(gains{1})(k);e=nan(1,numel(o.estimate));
for c=1:numel(e)
    at=C.recording==subject & strcmp(C.steps,S.labels{k}) & C.condition==c & ...
        abs(C.rawEstimate-o.estimate(c))<=1e-9*max(1,abs(o.estimate(c))) & ...
        abs(C.gain-gain)<=1e-9*max(1,abs(gain));
    values=C.actualError(at);
    assert(~isempty(values) && max(values)-min(values)<1e-8,'PipeCompare:Validation', ...
        'Cannot uniquely match a candidate readout to its truth rows.');
    e(c)=mean(values);
end
loss=mean(e.^2);contrast=e(1)-e(2);
end

function ci=recordInterval(x)
x=x(isfinite(x));ci=[NaN NaN];if numel(x)<2,return;end
rs=RandStream('mt19937ar','Seed',20261023);B=1999;
I=randi(rs,numel(x),numel(x),B);q=sort(mean(x(I),1));
at=max(1,min(B,ceil((B+1)*[0.025 0.975])));ci=q(at);
end
