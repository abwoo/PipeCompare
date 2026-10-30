function pipecompare_optimization_summary(root,outDir,stage)
% Record-balanced paired comparisons, with all missing recommendations kept.
if ~isfolder(outDir),mkdir(outDir);end
protocol=jsondecode(fileread(fullfile(fileparts(mfilename('fullpath')),'optimization_protocol.json')));
expected=protocol.([stage '_subjects']);expected=expected(:)';
folders=dir(fullfile(root,'optimization-public-*'));assert(numel(folders)==numel(expected), ...
    'PipeCompare:Validation','Missing planned public-record artifacts.');
allC=struct('standard',table(),'balanced',table());
S=struct('standard',pipecompare_record_summary(),'balanced',pipecompare_record_summary());
rows={};seen=[];spectra=table();stretches=table();sourceCommits={};
for k=1:numel(folders)
    folder=fullfile(folders(k).folder,folders(k).name);
    meta=jsondecode(fileread(fullfile(folder,'metadata.json')));seen(end+1)=meta.subject; %#ok<AGROW>
    assert(strcmp(meta.stage,stage),'PipeCompare:Validation','Mixed development and confirmation data.');
    sourceCommits{end+1}=meta.commit; %#ok<AGROW>
    D=readtable(fullfile(folder,'spectral.csv'));spectra=joinRows(spectra,D);
    D=readtable(fullfile(folder,'stretches.csv'),'TextType','string');stretches=joinRows(stretches,D);
    score=nan(1,2);raw=nan(1,2);available=false(1,2);sensitivity=nan(1,2);
    recipes={'standard','balanced'};
    for j=1:2
        recipe=recipes{j};part=fullfile(folder,recipe);
        execution=jsondecode(fileread(fullfile(part,'execution.json')));
        assert(execution.completed,'PipeCompare:Validation','Incomplete native recipe fit.');
        assert(~execution.recommended||execution.replayMaxError<1e-6,'PipeCompare:Validation','Native replay failed.');
        C=readtable(fullfile(part,'calibration.csv'),'TextType','char');
        if height(C)==0,C=pipecompare.eval.Calibration.empty();end
        allC.(recipe)=joinRows(allC.(recipe),C);
        D=load(fullfile(part,'score.mat'));S.(recipe)(end+1)=D.item;
        selected=C;if height(C)>0,selected=C(C.selected==1,:);end
        if ~isempty(selected)&&all(isfinite(selected.actualError))
            score(j)=mean(selected.actualError.^2);raw(j)=mean(selected.rawError.^2);
            available(j)=true;sensitivity(j)=mean(selected.associationSensitivityError.^2);
        end
    end
    rows(end+1,:)={meta.subject,meta.prevalence,available(1),available(2),score(1),score(2), ...
        score(2)-score(1),raw(1),raw(2),sensitivity(1),sensitivity(2)}; %#ok<AGROW>
end
assert(isequal(sort(seen),sort(expected)),'PipeCompare:Validation','Participant identity mismatch.');
assert(numel(unique(sourceCommits))==1,'PipeCompare:Validation','Mixed algorithm revisions.');
T=cell2table(rows,'VariableNames',{'subject','prevalence','standardAvailable','balancedAvailable', ...
    'standardMSE','balancedMSE','balancedMinusStandardMSE','standardRawMSE','balancedRawMSE', ...
    'standardAssociationRiskSquared','balancedAssociationRiskSquared'});
writetable(T,fullfile(outDir,'paired_public.csv'));
writetable(spectra,fullfile(outDir,'spectral.csv'));writetable(stretches,fullfile(outDir,'stretches.csv'));
paired=T.standardAvailable&T.balancedAvailable;
result=struct('commit',getenv('GITHUB_SHA'),'stage',stage,'planned',numel(expected), ...
    'sourceCommits',{unique(sourceCommits)},'paired',nnz(paired), ...
    'standardAvailable',nnz(T.standardAvailable),'balancedAvailable',nnz(T.balancedAvailable), ...
    'standardOnly',nnz(T.standardAvailable&~T.balancedAvailable), ...
    'balancedOnly',nnz(T.balancedAvailable&~T.standardAvailable), ...
    'neitherAvailable',nnz(~T.standardAvailable&~T.balancedAvailable), ...
    'meanPairedMSEDifference',mean(T.balancedMinusStandardMSE(paired)), ...
    'pairedDifferenceConditional95',recordInterval(T.balancedMinusStandardMSE(paired)), ...
    'standardPairedRMSE',sqrt(mean(T.standardMSE(paired))), ...
    'balancedPairedRMSE',sqrt(mean(T.balancedMSE(paired))), ...
    'scope','Paired program comparison on assessable participants; coverage reported separately; real-background injected-signal truth.');
for recipe={'standard','balanced'}
    name=recipe{1};C=allC.(name);writetable(C,fullfile(outDir,[name '_calibration.csv']));
    selected=C(C.selected==1,:);
    result.(name).selected=pipecompare.eval.Calibration.summarize(selected);
    result.(name).precision=validation_precision(result.(name).selected);
    if ~isempty(selected)
        D=selected;D.reportedError=D.associationSensitivityError;
        result.(name).associationSensitivity=pipecompare.eval.Calibration.summarize(D);
    end
    if strcmp(stage,'confirmation')
        groupOut=fullfile(outDir,[name '_group']);mkdir(groupOut);
        result.(name).groupTransfer=validation_group_transfer(S.(name),C,groupOut, ...
            protocol.confirmation_training_subjects,protocol.confirmation_test_subjects);
    end
end
% Report every duration, including the no-injection condition. The latter
% measures background marking, not false-positive truth on natural EEG.
rows={};
for mode={'long','multiscale'}
    for duration=unique(stretches.faultSeconds)'
        D=stretches(stretches.mode==string(mode{1})&stretches.faultSeconds==duration,:);
        rows(end+1,:)={mode{1},duration,height(D),nnz(D.assessed),mean(D.sampleRecall,'omitnan'), ...
            mean(D.outsideInjectedShare),mean(D.detectorRepairFaultRMSE,'omitnan'), ...
            mean(D.outsideInjectedRepairRMSE)}; %#ok<AGROW>
    end
end
T=cell2table(rows,'VariableNames',{'mode','seconds','cases','assessed','meanRecall', ...
    'outsideInjectedShare','faultRepairRMSE','outsideInjectedRepairRMSE'});
writetable(T,fullfile(outDir,'stretch_summary.csv'));result.stretches=table2struct(T);
result.spectralMaxReferenceError=max(spectra.referenceMaxError);
result.spectralUnresolvedRows=nnz(spectra.resolved~=1);
fid=fopen(fullfile(outDir,'summary.json'),'w');fprintf(fid,'%s\n',jsonencode(result));fclose(fid);
fprintf('PIPECOMPARE_OPTIMIZATION_JSON %s\n',jsonencode(result));
end

function T=joinRows(T,D)
if isempty(T),T=D;else,T=[T;D];end
end

function ci=recordInterval(x)
x=x(isfinite(x));ci=[NaN NaN];if numel(x)<2,return;end
rs=RandStream('mt19937ar','Seed',20261025);B=1999;
I=randi(rs,numel(x),numel(x),B);q=sort(mean(x(I),1));
at=max(1,min(B,ceil((B+1)*[0.025 0.975])));ci=q(at);
end
