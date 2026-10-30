function pipecompare_optimization_synthetic(outDir,seed,stage,legacyDir,recipe,maxSeconds)
% One recipe per worker, with exact checkpoints across bounded time slices.
if nargin<4,legacyDir='';end %#ok<NASGU> % old fallback controls are never imported here
if nargin<5
    for item={'standard','balanced'}
        pipecompare_optimization_synthetic(outDir,seed,stage,'',item{1},Inf);
    end
    return;
end
if nargin<6,maxSeconds=7200;end
protocolFile=fullfile(fileparts(mfilename('fullpath')),'optimization_protocol.json');
protocol=jsondecode(fileread(protocolFile));
assert(ismember(seed,protocol.([stage '_synthetic_seeds'])),'PipeCompare:Validation','Seed outside frozen stage.');
assert(ismember(recipe,{'standard','balanced'}),'PipeCompare:Validation','Unknown recipe.');
folder=fullfile(outDir,recipe);if ~isfolder(folder),mkdir(folder);end
varianceMode=protocol.artifact_variance.(recipe);
identity=struct('seed',seed,'stage',stage,'recipe',recipe,'artifactVariance',varianceMode, ...
    'revision',getenv('PIPECOMPARE_VALIDATION_REVISION'),'protocol',fileread(protocolFile));
assert(~isempty(identity.revision),'PipeCompare:Validation','Validation revision must be pinned.');
identityFile=fullfile(folder,'resume_identity.json');
if isfile(identityFile)
    assert(isequal(jsondecode(fileread(identityFile)),identity),'PipeCompare:Validation','Checkpoint belongs to a different experiment.');
else
    assert(~isfolder(fullfile(folder,'checkpoints')),'PipeCompare:Validation','Unidentified checkpoint state.');
    writeJson(identityFile,identity);
end
executionFile=fullfile(folder,'execution.json');
if isfile(executionFile)
    prior=jsondecode(fileread(executionFile));
    assert(prior.seed==seed&&strcmp(prior.recipe,recipe),'PipeCompare:Validation','Wrong execution identity.');
    if prior.completed,return;end
end
execution=struct('completed',false,'seed',seed,'recipe',recipe,'stage',stage, ...
    'artifactVariance',varianceMode,'sourceCommit',identity.revision,'legacyReused',false);
writeJson(executionFile,execution);
diary(fullfile(folder,'native.log'));closeLog=onCleanup(@()diary('off'));fitTimer=tic; %#ok<NASGU>
try
    pipecompare_benchmark(1,folder,seed,recipe,varianceMode,fullfile(folder,'checkpoints'),maxSeconds);
catch ME
    if ~strcmp(ME.identifier,'PipeCompare:Interrupted'),rethrow(ME);end
    execution.sliceSeconds=toc(fitTimer);execution.reason=ME.message;writeJson(executionFile,execution);
    fprintf('PIPECOMPARE_SLICE_INCOMPLETE seed=%d recipe=%s; checkpoints preserved, not a completed validation\n',seed,recipe);
    return;
end
if ~isfile(fullfile(folder,'calibration.csv'))
    writetable(pipecompare.eval.Calibration.empty(),fullfile(folder,'calibration.csv'));
end
execution.completed=true;execution.seconds=toc(fitTimer);
backendFile=getenv('PIPECOMPARE_NATIVE_BACKEND_REPORT');
assert(isfile(backendFile),'PipeCompare:Validation','Missing native backend evidence.');
execution.nativeBackend=jsondecode(fileread(backendFile));
writeJson(executionFile,execution);
fprintf('PIPECOMPARE_SLICE_COMPLETE seed=%d recipe=%s\n',seed,recipe);
end

function writeJson(file,value)
fid=fopen([file '.tmp'],'w');assert(fid>=0,'PipeCompare:Validation','Cannot save execution state.');
fprintf(fid,'%s\n',jsonencode(value));fclose(fid);movefile([file '.tmp'],file,'f');
end
