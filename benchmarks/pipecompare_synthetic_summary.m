function pipecompare_synthetic_summary(root,outDir,planned)
% Combine disjoint seed shards; resample recordings, never candidate rows.
if ~isfolder(outDir),mkdir(outDir);end
files=dir(fullfile(root,'benchmark-*','calibration.csv'));C=table();
for j=1:numel(files)
    T=readtable(fullfile(files(j).folder,files(j).name),'TextType','char');
    if isempty(C),C=T;else,C=[C;T];end %#ok<AGROW>
end
assert(~isempty(C),'PipeCompare:Validation','No completed simulation calibration output.');
writetable(C,fullfile(outDir,'calibration.csv'));
bench=table();files=dir(fullfile(root,'benchmark-*','benchmark.csv'));
for j=1:numel(files)
    T=readtable(fullfile(files(j).folder,files(j).name),'TextType','char');
    if isempty(bench),bench=T;else,bench=[bench;T];end %#ok<AGROW>
end
if ~isempty(bench),writetable(bench,fullfile(outDir,'benchmark.csv'));end
result=struct('commit',getenv('GITHUB_SHA'),'planned',planned, ...
    'completed',numel(unique(bench.recording)),'calibratableRecords',numel(unique(C.recording)),'all',pipecompare.eval.Calibration.summarize(C), ...
    'selected',pipecompare.eval.Calibration.summarize(C(C.selected==1,:)));
result.selectedRatioPrecision=validation_precision(result.selected);
result.rawPopulationAccuracy=validation_accuracy(C,'rawError');
result.correctedRealizedAccuracy=validation_accuracy(C,'conditionalError');
result.selectedRawPopulationAccuracy=validation_accuracy(C(C.selected==1,:),'rawError');
result.selectedUnavailable=planned-numel(unique(C.recording(C.selected==1)));
result.rawDenoisingComparators=struct();
for name={'PipeCompare','minimal','common'}
    D=bench(strcmp(bench.pipeline,name{1}),:);valid=isfinite(D.errDifference);
    stats=struct('records',nnz(valid),'unavailable',height(D)-nnz(valid), ...
        'differenceRMSE',sqrt(mean(D.errDifference(valid).^2)), ...
        'differenceBias',mean(D.errDifference(valid)));
    result.rawDenoisingComparators.(name{1})=stats;
end
result.sourceRun=getenv('VALIDATION_SOURCE_RUN');
fprintf('PIPECOMPARE_SIMULATION_JSON %s\n',jsonencode(result));
fid=fopen(fullfile(outDir,'summary.json'),'w');fprintf(fid,'%s\n',jsonencode(result));fclose(fid);
fprintf('%s\n',pipecompare.eval.Calibration.describe(C,'all candidates'));
fprintf('%s\n',pipecompare.eval.Calibration.describe(C(C.selected==1,:),'selected candidates'));
assert(result.completed==planned,'PipeCompare:Validation','Only %d of %d planned simulation records completed.',result.completed,planned);
end
