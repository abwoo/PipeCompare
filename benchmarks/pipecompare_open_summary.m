function pipecompare_open_summary(root,outDir)
% Source-aware, record-balanced summaries; no raw EEG is loaded here.
if ~isfolder(outDir),mkdir(outDir);end
folders=dir(fullfile(root,'public-validation-*'));
folders=folders(~cellfun('isempty',regexp({folders.name},'^public-validation-[0-5]$','once')));
S=pipecompare_record_summary();C=table();sourceHeads={};sourceTimers={};expectedSubjects=[];trainingSubjects=[];testSubjects=[];
metrics={'blink','spectral','stretches','recordings','real_erp'};
components=struct();sourceFailures=0;coverage=table();
for k=1:numel(metrics)
    combined=table();
    for j=1:numel(folders)
        file=fullfile(folders(j).folder,folders(j).name,[metrics{k} '.csv']);
        if ~isfile(file),continue;end
        opts=detectImportOptions(file);
        textVars=intersect(opts.VariableNames,{'recording','steps','reason'});
        if ~isempty(textVars),opts=setvartype(opts,textVars,'string');end
        T=readtable(file,opts);
        if isempty(combined),combined=T;else,combined=[combined;T];end %#ok<AGROW>
    end
    if ~isempty(combined)
        writetable(combined,fullfile(outDir,[metrics{k} '.csv']));
        components.(metrics{k})=componentSummary(metrics{k},combined);
        if strcmp(metrics{k},'recordings'),coverage=combined;end
    end
end
for j=1:numel(folders)
    folder=fullfile(folders(j).folder,folders(j).name);
    if isfile(fullfile(folder,'metadata.json'))
        meta=jsondecode(fileread(fullfile(folder,'metadata.json')));
        sourceHeads{end+1}=meta.commit; %#ok<AGROW>
        sourceTimers{end+1}=pipecompare.utils.fieldOr(meta,'fitTimer','implicit'); %#ok<AGROW>
        subjects=sort([meta.protocol.physionet.files.subject]);
        if isempty(expectedSubjects)
            expectedSubjects=subjects;
            if isfield(meta.protocol,'group_independent_test')
                trainingSubjects=meta.protocol.group_independent_test.training_subjects(:)';
                testSubjects=meta.protocol.group_independent_test.test_subjects(:)';
            else
                trainingSubjects=subjects(1:10);testSubjects=subjects(11:end);
            end
        else
            assert(isequal(expectedSubjects,subjects),'PipeCompare:Validation','Mixed source participant cohorts.');
        end
    end
    if isfile(fullfile(folder,'execution.json'))
        execution=jsondecode(fileread(fullfile(folder,'execution.json')));
        sourceFailures=sourceFailures+~execution.completed;
    else,sourceFailures=sourceFailures+1;end
    part=pipecompare_record_summary();
    if isfile(fullfile(folder,'derived_scores.mat'))
        D=load(fullfile(folder,'derived_scores.mat'));part=D.S;
    end
    checkpoints=dir(fullfile(folder,'score_*.mat'));
    for q=1:numel(checkpoints)
        D=load(fullfile(checkpoints(q).folder,checkpoints(q).name));
        at=find(strcmp({part.name},D.item.name),1);
        if isempty(at),part(end+1)=D.item;else,part(at)=D.item;end %#ok<AGROW>
    end
    S=[S part]; %#ok<AGROW>
    if isfile(fullfile(folder,'calibration.csv'))
        D=readtable(fullfile(folder,'calibration.csv'),'TextType','char');
        if isempty(C),C=D;else,C=[C;D];end %#ok<AGROW>
    end
end
if isfield(components,'recordings')
    components.recordings.fitTimingAssessed=numel(sourceTimers)==numel(folders)&& ...
        all(strcmp(sourceTimers,'explicit-tic-handle'));
    if ~components.recordings.fitTimingAssessed
        components.recordings.legacyReportedFitSeconds=components.recordings.totalSeconds;
        components.recordings.totalSeconds=NaN;
    end
end
result=struct('commit',getenv('GITHUB_SHA'),'recordsWithScores',numel(S), ...
    'plannedRecords',30,'sourceShards',numel(folders),'noCalibratableRecordings',30, ...
    'noSelectedCalibration',30,'sourceFailedShards',sourceFailures,'components',components, ...
    'scope','Independent real-background injection; not certified clean human neural truth');
if ~isempty(C)
    writetable(C,fullfile(outDir,'calibration.csv'));
    result.all=pipecompare.eval.Calibration.summarize(C);
    result.selected=pipecompare.eval.Calibration.summarize(C(C.selected==1,:));
    result.selectedRatioPrecision=validation_precision(result.selected);
    fprintf('%s\n',pipecompare.eval.Calibration.describe(C,'all candidates'));
    fprintf('%s\n',pipecompare.eval.Calibration.describe(C(C.selected==1,:),'selected candidates'));
    result.rawPopulationAccuracy=validation_accuracy(C,'rawError');
    result.correctedRealizedAccuracy=validation_accuracy(C,'conditionalError');
    result.rawRealizedAccuracy=validation_accuracy(C,'rawConditionalError');
    result.selectedRawPopulationAccuracy=validation_accuracy(C(C.selected==1,:),'rawError');
    result.selectedCorrectedPopulationAccuracy=validation_accuracy(C(C.selected==1,:),'actualError');
    result.selectedCorrectedRealizedAccuracy=validation_accuracy(C(C.selected==1,:),'conditionalError');
    result.selectedRawRealizedAccuracy=validation_accuracy(C(C.selected==1,:),'rawConditionalError');
    result.noCalibratableRecordings=30-numel(unique(C.recording));
    result.noSelectedCalibration=30-numel(unique(C.recording(C.selected==1)));
    result.calibrationByPrevalence=prevalenceCalibration(C,coverage);
end
if ~isempty(expectedSubjects)
    assert(numel(trainingSubjects)==10&&numel(testSubjects)==20&& ...
        isempty(intersect(trainingSubjects,testSubjects))&& ...
        isequal(sort([trainingSubjects testSubjects]),expectedSubjects), ...
        'PipeCompare:Validation','Training/test participants must be disjoint and match the frozen cohort.');
end
rows={};groupCoverage={};
if numel(S)>=3
    for held=1:numel(S)
        train=setdiff(1:numel(S),held);keys=S(train(1)).keys;
        for j=train(2:end),keys=intersect(keys,S(j).keys,'stable');end
        [present,loc]=ismember(keys,S(held).keys);keep=find(present);
        if isempty(keep),continue;end
        [Y,V,covariance]=groupInputs(S(train),keys);
        [~,~,~,prior]=pipecompare.Group.twoLevel(Y,V,covariance);
        [y,v,c]=groupInputs(S(held),keys);
        [theta,info]=pipecompare.Group.applyPrior(prior,y,v,c);
        [a,b,valid]=validation_group_choices(y,v,theta);
        groupCoverage(end+1,:)={S(held).name,nnz(isfinite(y)&y>0), ...
            nnz(isfinite(y)&y>0&isfinite(v)&v>=0),numel(valid),info.assessed,info.reason{1}}; %#ok<AGROW>
        if isempty(valid),continue;end
        assert(info.assessed||a==b,'PipeCompare:Validation','Unassessed prior changed a matched-support raw choice.');
        T=S(held).table;[~,f]=min(T.shrunk(loc(valid)));f=valid(f);
        errors=nan(1,numel(keys));
        for p=valid
            errors(p)=candidateLoss(C,S(held),loc(p));
        end
        oracle=min(errors,[],'omitnan');
        [recommended,~,inFamily]=validation_recommended_loss(C,S(held),keys);
        rows(end+1,:)={S(held).name,numel(train),numel(valid),info.assessed, ...
            errors(a),errors(f),errors(b),oracle,recommended,inFamily, ...
            errors(b)-recommended,errors(b)-errors(a),a,b}; %#ok<AGROW>
    end
end
if ~isempty(groupCoverage)
    coverageTable=cell2table(groupCoverage,'VariableNames',{'recording','pointEligibleCandidates', ...
        'finiteCovarianceCandidates','comparableCandidates','priorAssessed','reason'});
    writetable(coverageTable,fullfile(outDir,'participant_held_out_coverage.csv'));
    result.groupNoComparableScore=nnz(coverageTable.comparableCandidates==0);
    result.groupPriorAssessedRecords=nnz(coverageTable.priorAssessed);
end
if ~isempty(rows)
    T=cell2table(rows,'VariableNames',{'recording','trainingParticipants','eligibleCandidates','priorAssessed', ...
        'rawMinimumMSE','factorialMinimumMSE','frozenGroupMSE','familyOracleMSE', ...
        'actualRecommendationMSE','recommendationInTrainingFamily','groupMinusRecommendationMSE', ...
        'groupMinusRawMinimumMSE','rawChoice','groupChoice'});
    writetable(T,fullfile(outDir,'participant_held_out.csv'));
    result.groupTransferRecords=height(T);result.groupTransferAssessed=nnz(T.priorAssessed);
    result.meanGroupMinusRawMinimumMSE=mean(T.groupMinusRawMinimumMSE,'omitnan');
    result.meanGroupMinusRecommendationMSE=mean(T.groupMinusRecommendationMSE,'omitnan');
    result.actualRecommendationsInTrainingFamily=nnz(T.recommendationInTrainingFamily);
end
result.sourceRun=getenv('VALIDATION_SOURCE_RUN');result.sourceCommits=unique(sourceHeads);
result.sourceParticipantIds=expectedSubjects;result.trainingSubjectIds=trainingSubjects;result.testSubjectIds=testSubjects;
result.groupFamilyScope='Raw minima and oracle use training-supported keys; actual recommendations may use the full per-record candidate family.';
result.fixedPriorTest=fixedPriorEvaluation(S,C,outDir,trainingSubjects,testSubjects);
fid=fopen(fullfile(outDir,'summary.json'),'w');fprintf(fid,'%s\n',jsonencode(result));fclose(fid);
fprintf('PIPECOMPARE_VALIDATION_JSON %s\n',jsonencode(result));
assert(numel(S)==30 && numel(unique({S.name}))==30,'PipeCompare:Validation', ...
    'Only %d of 30 planned records have saved scores.',numel(S));
assert(sourceFailures==0,'PipeCompare:Validation','%d source shards recorded execution/invariant failures.',sourceFailures);
ids=arrayfun(@(s)sscanf(s.name,'S%dR01.mat'),S);
assert(numel(expectedSubjects)==30&&numel(unique(expectedSubjects))==30&&isequal(sort(ids(:))',expectedSubjects), ...
    'PipeCompare:Validation','Source participants differ from the frozen metadata cohort.');
expected={'blink','records',15;'spectral','participants',30;'stretches','participants',30; ...
    'recordings','participants',30;'real_erp','participants',6};
for k=1:size(expected,1)
    name=expected{k,1};field=expected{k,2};count=expected{k,3};
    assert(isfield(components,name)&&components.(name).(field)==count,'PipeCompare:Validation', ...
        'The %s evaluation did not cover all %d planned participants.',name,count);
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

function S=componentSummary(name,T)
S=struct('rows',height(T));
switch name
    case 'blink'
        S.records=height(T);S.annotations=sum(T.annotations);
        S.macroPrecision=mean(T.precision);S.macroRecall=mean(T.recall);
        S.macroF1=mean(T.f1);S.f1Range=[min(T.f1) max(T.f1)];
        S.macroF1RecordBootstrap95=recordInterval(T.f1);
        S.meanFalseBlinkPredictionsPerMinute=mean(T.falseEventsPerMinute);
        S.unassessed=nnz(T.assessed~=1);
        S.scope='Blink annotations only; unmatched detections can include other ocular activity.';
    case 'spectral'
        S.participants=numel(unique(T.subject));
        S.maxReferenceAbsError=max(T.referenceMaxError);
        S.unresolvedRows=nnz(T.resolved~=1);
        S.scope='Matched pmtm/periodogram estimator agreement; not latent-spectrum unbiasedness.';
    case 'stretches'
        S.participants=numel(unique(T.subject));rows={};
        for duration=[1 5 10 20]
            for adjacent=[0 1]
                D=T(T.faultSeconds==duration&T.adjacent==adjacent,:);
                if isempty(D),continue;end
                rows(end+1,:)={duration,adjacent,height(D),nnz(D.assessed==1), ...
                    mean(D.sampleRecall),mean(D.falsePositiveShare), ...
                    mean(D.corruptedRMSE),mean(D.oracleRepairRMSE), ...
                    max(D.boundaryJumpUv),mean(D.meanBiasUv),mean(D.alphaChangeLog10PSD)}; %#ok<AGROW>
            end
        end
        S.cases=table2struct(cell2table(rows,'VariableNames',{'faultSeconds','adjacent','participants', ...
            'assessed','meanSampleRecallIncludingUnassessed','meanFalsePositiveShare', ...
            'corruptedRMSE','oracleRepairRMSE','maximumBoundaryJumpUv','meanBiasUv','meanAlphaChangeLog10PSD'}));
        S.scope='Oracle reconstruction is separate from detection; target is the original observed EEG.';
    case 'recordings'
        S.participants=numel(unique(T.subject));S.recommendations=nnz(T.recommended==1);
        S.noRecommendation=height(T)-S.recommendations;
        S.maximumReplayError=max(T.replayError,[],'omitnan');
        S.totalSeconds=sum(T.seconds);
    case 'real_erp'
        S.participants=numel(unique(T.recording));
        S.selectedRecords=numel(unique(T.recording(T.selected==1)));
        S.meanSelectedRawRareMinusFrequent=mean(T.rawRareMinusFrequentUv(T.selected==1),'omitnan');
        S.scope='Real-task execution and precision; clean neural truth is unknown.';
end
end

function ci=recordInterval(x)
x=x(isfinite(x));ci=[NaN NaN];if numel(x)<2,return;end
rs=RandStream('mt19937ar','Seed',20261023);B=1999;
I=randi(rs,numel(x),numel(x),B);q=sort(mean(x(I),1));
at=max(1,min(B,ceil((B+1)*[0.025 0.975])));ci=q(at);
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

function result=fixedPriorEvaluation(S,C,outDir,trainingSubjects,testSubjects)
% The participant split comes from frozen source metadata. Intervals are
% conditional on that fitted training prior, not on new training cohorts.
id=arrayfun(@(s)sscanf(s.name,'S%dR01.mat'),S);
train=find(ismember(id,trainingSubjects));test=find(ismember(id,testSubjects));
train=train(:)';test=test(:)';
result=struct('trainingRequested',10,'testRequested',20,'trainingObserved',0, ...
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

function S=prevalenceCalibration(C,coverage)
S=struct([]);
if isempty(coverage),return;end
for p=unique(coverage.blinkPrevalence)'
    people=coverage.subject(coverage.blinkPrevalence==p);
    D=C(ismember(C.recording,people),:);
    all=pipecompare.eval.Calibration.summarize(D);
    selected=pipecompare.eval.Calibration.summarize(D(D.selected==1,:));
    entry=struct('prevalence',p,'plannedParticipants',numel(people),'all',all,'selected',selected);
    if isempty(S),S=entry;else,S(end+1)=entry;end %#ok<AGROW>
end
end
