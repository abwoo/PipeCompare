function pipecompare_optimization_synthetic(outDir,seed,stage,legacyDir)
% Same seed and full native search; completed legacy controls can be reused.
if nargin<4,legacyDir='';end
protocol=jsondecode(fileread(fullfile(fileparts(mfilename('fullpath')),'optimization_protocol.json')));
assert(ismember(seed,protocol.([stage '_synthetic_seeds'])),'PipeCompare:Validation','Seed outside frozen stage.');
for recipe={'standard','balanced'}
    folder=fullfile(outDir,recipe{1});if ~isfolder(folder),mkdir(folder);end
    varianceMode='pooled';
    if strcmp(recipe{1},'standard')
        varianceMode='legacy-total';
        if strcmp(stage,'development')&&validation_reuse_legacy(legacyDir,folder,'synthetic',seed,[],[]),continue;end
    end
    fid=fopen(fullfile(folder,'execution.json'),'w');
    fprintf(fid,'%s\n',jsonencode(struct('completed',false,'seed',seed,'recipe',recipe{1})));fclose(fid);
    diary(fullfile(folder,'native.log'));closeLog=onCleanup(@()diary('off'));fitTimer=tic;
    pipecompare_benchmark(1,folder,seed,recipe{1},varianceMode);
    seconds=toc(fitTimer);diary('off');clear closeLog
    if ~isfile(fullfile(folder,'calibration.csv'))
        writetable(pipecompare.eval.Calibration.empty(),fullfile(folder,'calibration.csv'));
    end
    fid=fopen(fullfile(folder,'execution.json'),'w');
    execution=struct('completed',true,'seed',seed,'recipe',recipe{1},'artifactVariance',varianceMode,'legacyReused',false,'seconds',seconds);
    backendFile=getenv('PIPECOMPARE_NATIVE_BACKEND_REPORT');
    if isfile(backendFile),execution.nativeBackend=jsondecode(fileread(backendFile));end
    fprintf(fid,'%s\n',jsonencode(execution));
    fclose(fid);
end
end
