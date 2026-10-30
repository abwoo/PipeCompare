function pipecompare_optimization_synthetic_summary(root,outDir,stage)
if ~isfolder(outDir),mkdir(outDir);end
protocol=jsondecode(fileread(fullfile(fileparts(mfilename('fullpath')),'optimization_protocol.json')));
expected=protocol.([stage '_synthetic_seeds']);expected=expected(:)';
folders=dir(fullfile(root,'optimization-synthetic-*'));
assert(numel(folders)==numel(expected),'PipeCompare:Validation','Missing planned simulation artifacts.');
result=struct('commit',getenv('GITHUB_SHA'),'stage',stage,'planned',numel(expected));
allC=struct();
for recipe={'standard','balanced'}
    name=recipe{1};C=table();bench=table();
    for k=1:numel(folders)
        path=fullfile(folders(k).folder,folders(k).name,name);
        execution=jsondecode(fileread(fullfile(path,'execution.json')));
        assert(execution.completed,'PipeCompare:Validation','Incomplete simulation fit.');
        D=readtable(fullfile(path,'calibration.csv'),'TextType','char');
        if height(D)==0,D=pipecompare.eval.Calibration.empty();end
        if isempty(C),C=D;else,C=[C;D];end %#ok<AGROW>
        D=readtable(fullfile(path,'benchmark.csv'),'TextType','char');
        if isempty(bench),bench=D;else,bench=[bench;D];end %#ok<AGROW>
    end
    assert(isequal(sort(unique(bench.recording))',sort(expected)), ...
        'PipeCompare:Validation','Completed simulation seeds differ from protocol.');
    allC.(name)=C;writetable(C,fullfile(outDir,[name '_calibration.csv']));
    selected=C(C.selected==1,:);result.(name).selected=pipecompare.eval.Calibration.summarize(selected);
    result.(name).precision=validation_precision(result.(name).selected);
    D=selected;D.reportedError=D.associationSensitivityError;
    result.(name).associationSensitivity=pipecompare.eval.Calibration.summarize(D);
    result.(name).unavailable=numel(expected)-numel(unique(selected.recording));
    result.(name).sourceBenchRecords=numel(unique(bench.recording));
end
rows={};
for seed=expected
    a=allC.standard(allC.standard.recording==seed&allC.standard.selected==1,:);
    b=allC.balanced(allC.balanced.recording==seed&allC.balanced.selected==1,:);
    x=mean(a.actualError.^2);y=mean(b.actualError.^2);
    rows(end+1,:)={seed,~isempty(a),~isempty(b),x,y,y-x}; %#ok<AGROW>
end
T=cell2table(rows,'VariableNames',{'seed','standardAvailable','balancedAvailable','standardMSE','balancedMSE','difference'});
writetable(T,fullfile(outDir,'paired_synthetic.csv'));paired=T.standardAvailable&T.balancedAvailable;
result.paired=nnz(paired);result.meanPairedMSEDifference=mean(T.difference(paired));
result.standardPairedRMSE=sqrt(mean(T.standardMSE(paired)));result.balancedPairedRMSE=sqrt(mean(T.balancedMSE(paired)));
result.pairedDifferenceConditional95=recordInterval(T.difference(paired));
fid=fopen(fullfile(outDir,'summary.json'),'w');fprintf(fid,'%s\n',jsonencode(result));fclose(fid);
fprintf('PIPECOMPARE_OPTIMIZATION_SYNTHETIC_JSON %s\n',jsonencode(result));
end

function ci=recordInterval(x)
x=x(isfinite(x));ci=[NaN NaN];if numel(x)<2,return;end
rs=RandStream('mt19937ar','Seed',20261025);B=1999;
I=randi(rs,numel(x),numel(x),B);q=sort(mean(x(I),1));
at=max(1,min(B,ceil((B+1)*[0.025 0.975])));ci=q(at);
end
