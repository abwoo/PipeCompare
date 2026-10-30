function [result, state] = pipecompare_checkpoint_search(folder, passport, startSearch, maxSeconds, mode)
% One full-family search slice, with strict identity and durable progress.
% A completed slice chain is a separate execution, not evidence added to an
% old partial run. Adaptive leaves may be re-derived by Executor.resume.
if nargin < 5, mode = 'start'; end
assert(ischar(folder) && ~isempty(folder), 'PipeCompare:Slice', 'Checkpoint folder is required.');
validateattributes(maxSeconds, {'numeric'}, {'scalar','real','finite','positive','<=',18000});
assert(ismember(mode, {'start','resume'}), 'PipeCompare:Slice', 'Mode must be start or resume.');
required = {'commit','seed','arm','artifactVariance','eeglabSha'};
for k = 1:numel(required)
    assert(isfield(passport, required{k}) && ~isempty(passport.(required{k})), ...
        'PipeCompare:Slice', 'Missing checkpoint passport field %s.', required{k});
end
passport.matlabVersion = version();
assert(~isempty(which('eeg_getversion')), 'PipeCompare:Slice', 'EEGLAB must be loaded.');
passport.eeglabVersion = eeg_getversion();
assert(~isfield(passport, 'checkpointDir'), 'PipeCompare:Slice', 'Passport must be portable across workers.');
if ~isfolder(folder), mkdir(folder); end
meta = fullfile(folder, 'checkpoint_passport.json');
manifest = fullfile(folder, 'manifest.mat');
if isfile(meta)
    prior = jsondecode(fileread(meta));
    assert(isequaln(prior, passport), 'PipeCompare:Slice', ...
        'Checkpoint passport mismatch: code, seed, arm, variance and toolchain must match.');
elseif strcmp(mode, 'resume') || isfile(manifest)
    error('PipeCompare:Slice', 'Refusing an unbound checkpoint: checkpoint_passport.json is missing.');
else
    fid = fopen([meta '.tmp'],'w');
    assert(fid >= 0,'PipeCompare:Slice','Cannot write checkpoint passport.');
    clean = onCleanup(@()fclose(fid));
    fprintf(fid,'%s\n',jsonencode(passport)); clear clean;
    movefile([meta '.tmp'],meta,'f');
end
if strcmp(mode, 'resume')
    assert(isfile(manifest), 'PipeCompare:Slice', 'Resume requires manifest.mat.');
    bound = load(manifest,'result');
    assert(isfield(bound.result.options,'artifactVariance') && ...
        strcmp(bound.result.options.artifactVariance,passport.artifactVariance), ...
        'PipeCompare:Slice', 'Stored manifest variance conflicts with its passport.');
else
    assert(~isfile(manifest), 'PipeCompare:Slice', 'Checkpoint already started; use resume.');
end
state = struct('status','running','execution','checkpointed-full-family', ...
    'complete',false,'partialSearch',true,'exactResume',false, ...
    'adaptiveLeavesMayBeRederived',true,'passport',passport, ...
    'leafFilesBefore',numel(dir(fullfile(folder,'leaf_*.mat'))), ...
    'icaFilesBefore',numel(dir(fullfile(folder,'ica_cache','*.mat'))), ...
    'wallSeconds',0,'reason','');
result = [];
finished = fullfile(folder, 'completed_result.mat');
if strcmp(mode, 'resume') && isfile(finished)
    saved = load(finished, 'result', 'passport');
    assert(isequaln(saved.passport, passport), 'PipeCompare:Slice', 'Completed result passport mismatch.');
    result = saved.result; result.options.checkpoint = folder;
    state.status = 'complete'; state.complete = true; state.partialSearch = false;
    state.reason = 'Reused completed result; no search executed.';
    writeState(folder, state); return;
end
writeState(folder, state); clock = tic;
try
    if strcmp(mode, 'resume')
        result = pipecompare.run.Executor.resume(folder, struct('maxSeconds',maxSeconds,'strictCheckpoint',true));
    else
        result = startSearch();
    end
    assert(isfield(result,'notRun') && result.notRun == 0, 'PipeCompare:Slice', ...
        'Search returned an incomplete candidate family.');
    assert(isfield(result.options,'artifactVariance') && strcmp(result.options.artifactVariance,passport.artifactVariance), ...
        'PipeCompare:Slice', 'Result variance does not match its passport.');
    state.status = 'complete'; state.complete = true; state.partialSearch = false;
    tmp = fullfile(folder, 'completed_result.tmp.mat');
    save(tmp, 'result', 'passport', '-v7.3'); movefile(tmp, finished, 'f');
catch ME
    state.reason = ME.message;
    if strcmp(ME.identifier, 'PipeCompare:Interrupted')
        state.status = 'partial';
    else
        state.status = 'failed'; state.wallSeconds = toc(clock);
        writeState(folder, state); rethrow(ME);
    end
end
state.wallSeconds = toc(clock);
state.leafFilesAfter = numel(dir(fullfile(folder,'leaf_*.mat')));
state.icaFilesAfter = numel(dir(fullfile(folder,'ica_cache','*.mat')));
writeState(folder, state);
end

function writeState(folder, state)
file = fullfile(folder, 'slice_progress.json');
fid = fopen([file '.tmp'], 'w');
assert(fid >= 0, 'PipeCompare:Slice', 'Cannot write slice progress.');
clean = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', jsonencode(state)); clear clean;
movefile([file '.tmp'], file, 'f');
end
