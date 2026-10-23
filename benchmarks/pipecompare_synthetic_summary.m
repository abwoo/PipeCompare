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
result=struct('commit',getenv('GITHUB_SHA'),'planned',planned, ...
    'completed',numel(unique(C.recording)),'all',pipecompare.eval.Calibration.summarize(C), ...
    'selected',pipecompare.eval.Calibration.summarize(C(C.selected,:)));
fid=fopen(fullfile(outDir,'summary.json'),'w');fprintf(fid,'%s\n',jsonencode(result));fclose(fid);
fprintf('%s\n',pipecompare.eval.Calibration.describe(C,'all candidates'));
fprintf('%s\n',pipecompare.eval.Calibration.describe(C(C.selected,:),'selected candidates'));
assert(result.completed==planned,'PipeCompare:Validation','Only %d of %d planned simulation records completed.',result.completed,planned);
end
