function out = pipecompare_shadow_selection(outDir, firstSeed, nSeeds, varargin)
% Fresh-record fitted-family comparison and saved-output aggregation.
% search: unsliced full-family search; merge: original unsliced fragments.
% slice/resume: one full-family checkpointed search, explicitly labelled as
% a separate execution and excluded from the original unsliced paired mean.
% micro: design-seed diagnosis only. test: known fake candidate statistics.
% The current rule and revision-specific evidence are in docs/PR50_STATUS.md.
% No mode changes Rank, fits on confirmation seeds, or promotes a result.
%
% ---- argument parsing (manual Name-Value pairs, no toolbox needed) ----
if nargin < 2 || isempty(firstSeed), firstSeed = 20031; end
if nargin < 3 || isempty(nSeeds), nSeeds = 8; end
mode = 'search'; headSha = ''; workflowRunId = ''; maxSeconds = Inf;
checkpointDir = ''; sliceSeconds = Inf; seedList = [];
nStrata = []; icaSeeds = []; refineCutoffs = []; eogCheck = [];
reducedFamily = ''; reducedTrials = ''; sliceDiagnostic = false; artifactVariance = '';
armMode = 'all';
testArms = []; fragmentsDir = ''; eeglabVersion = getenv('PIPECOMPARE_EEGLAB_VERSION');
if isempty(eeglabVersion), eeglabVersion = 'unknown'; end
k = 1;
while k <= numel(varargin)
    name = varargin{k};
    assert(ischar(name) && k < numel(varargin), 'PipeCompare:Shadow', 'Name-Value arguments must come in pairs.');
    value = varargin{k + 1};
    switch lower(name)
        case 'mode', mode = value;
        case 'arm', armMode = value;
        case 'headsha', headSha = value;
        case 'workflowrunid', workflowRunId = value;
        case 'maxseconds', maxSeconds = value;
        case 'checkpointdir', checkpointDir = value;
        case 'seedlist', seedList = value;
        case 'sliceseconds', sliceSeconds = value;
        case 'nstrata', nStrata = value;
        case 'icaseeds', icaSeeds = value;
        case 'refinecutoffs', refineCutoffs = value;
        case 'eogcheck', eogCheck = value;
        case 'reducedfamily', reducedFamily = value;
        case 'reducedtrials', reducedTrials = value;
        case 'slicediagnostic', sliceDiagnostic = value;
        case 'artifactvariance', artifactVariance = value;
        case 'testarms', testArms = value;
        case 'fragmentsdir', fragmentsDir = value;
        case 'eeglabversion', eeglabVersion = value;
        otherwise, error('PipeCompare:Shadow', 'Unknown option %s.', name);
    end
    k = k + 2;
end
assert(ischar(outDir) && ~isempty(outDir), 'PipeCompare:Shadow', 'outDir must be a non-empty path.');
assert(isscalar(firstSeed) && firstSeed == floor(firstSeed) && firstSeed > 0, 'PipeCompare:Shadow', 'firstSeed must be a positive integer.');
assert(isscalar(nSeeds) && nSeeds == floor(nSeeds) && nSeeds >= 1 && nSeeds <= 30, 'PipeCompare:Shadow', 'nSeeds must be an integer in 1..30.');
assert(ismember(mode, {'search', 'merge', 'test', 'resume', 'slice', 'micro', 'resume-fast'}), 'PipeCompare:Shadow', 'mode must be search, merge, test, slice, resume, micro or resume-fast.');
% Checkpointed full-family executions remain separate from unsliced evidence.
if ismember(mode, {'slice','resume'})
    assert(nSeeds == 1 && ismember(armMode, {'standard','balanced'}), 'PipeCompare:Shadow', 'A slice runs one seed and one search.');
    assert(~isempty(checkpointDir) && isfinite(maxSeconds) && maxSeconds > 0 && maxSeconds <= 18000, ...
        'PipeCompare:Shadow', 'A slice requires a checkpoint directory and finite budget <= 18000 seconds.');
    assert(~isempty(getenv('EEGLAB_SHA')), 'PipeCompare:Shadow', 'EEGLAB_SHA is required for checkpoint identity.');
    if strcmp(mode, 'resume')
        assert(isfile(fullfile(checkpointDir,'manifest.mat')), 'PipeCompare:Shadow', 'Resume requires an existing checkpoint.');
    else
        assert(~isfile(fullfile(checkpointDir,'manifest.mat')), 'PipeCompare:Shadow', 'Existing checkpoint requires resume mode.');
    end
end
% Round 10 fast modes. Both validate BEFORE any search, so a bad
% fast recipe fails in seconds, never after hours of searching.
if strcmp(mode, 'micro')
    % A diagnostic slice runs one design seed and one search.
    % Scientific slice/resume executions use separate explicit modes.
    assert(nSeeds == 1, 'PipeCompare:Shadow', 'micro mode requires exactly 1 seed.');
    assert(ismember(armMode, {'standard', 'balanced'}), 'PipeCompare:Shadow', 'micro mode requires one arm (standard or balanced), not ''all''.');
    assert(isscalar(firstSeed) && ismember(firstSeed, 40101:40104), 'PipeCompare:Shadow', 'micro mode accepts design seeds 40101-40104 only; a design seed is diagnostic, never held-out evidence.');
    assert(isempty(seedList), 'PipeCompare:Shadow', 'micro mode does not accept a seedList; it runs exactly one design seed.');
    assert(~isempty(checkpointDir), 'PipeCompare:Shadow', 'micro mode requires a non-empty checkpointDir; empty checkpoint plus Inf time is forbidden.');
    assert(isfinite(maxSeconds) && maxSeconds > 0 && maxSeconds <= 2400, 'PipeCompare:Shadow', 'micro mode requires finite maxSeconds in (0, 2400].');
    assert(sliceSeconds == 2400, 'PipeCompare:Shadow', 'micro mode requires sliceSeconds == 2400 (40 minutes); changing it is a new recipe, not micro.');
    assert(isequal(nStrata, 1), 'PipeCompare:Shadow', 'micro mode requires nStrata == 1 (one bench_synth record).');
    assert(isequal(icaSeeds, 1) && isequal(refineCutoffs, false) && isequal(eogCheck, false), 'PipeCompare:Shadow', 'micro recipe requires icaSeeds==1, refineCutoffs==false, eogCheck==false (diagnostic only).');
    assert(ischar(reducedFamily) && ~isempty(reducedFamily) && ischar(reducedTrials) && ~isempty(reducedTrials), 'PipeCompare:Shadow', 'micro recipe requires non-empty reducedFamily and reducedTrials rule text.');
    if ~isequal(sliceDiagnostic, true)
        error('PipeCompare:Shadow', 'Reduced-family micro execution is not wired in Round 10 (pop_pipecompare has no icaSeeds/refineCutoffs/eogCheck options and Presets.searchOptions() hardcodes them); pass slicediagnostic=true to run the slice-diagnostic run on the production family instead, labelled slice-diagnostic, never as a reduced-family result.');
    end
end
if strcmp(mode, 'resume-fast')
    % Guard-only: validate the slice manifest, write a marker, and
    % stop. This mode never calls Executor.resume and never writes
    % a held-out fragment; actual resume execution is a later step
    % after this guard is proven in CI.
    assert(~isempty(checkpointDir), 'PipeCompare:Shadow', 'resume-fast mode requires a non-empty checkpointDir.');
    assert(ischar(artifactVariance) && ~isempty(artifactVariance), 'PipeCompare:Shadow', 'resume-fast mode requires an explicit artifactvariance.');
    manifestFile = fullfile(checkpointDir, 'slice_manifest.json');
    assert(isfile(manifestFile), 'PipeCompare:Shadow', 'resume-fast mode requires slice_manifest.json in the checkpoint directory.');
    manifest = jsondecode(fileread(manifestFile));
    guard = validate_resume_checkpoint(manifest, version(), eeglabVersion, artifactVariance);
    marker = struct('mode', 'resume-fast', 'validated', guard.ok, 'exactResume', guard.exactResume, ...
        'adaptiveLeavesRederived', guard.adaptiveLeavesRederived, 'note', guard.note, ...
        'manifest', manifest);
    if ~isfolder(outDir), mkdir(outDir); end
    writeJson(fullfile(outDir, 'slice_validation_marker.json'), marker);
    out = struct('mode', mode, 'seeds', firstSeed, 'perSeed', struct('seed', {}, 'arms', {}, 'guardChanged', {}, 'guardDecision', {}), 'runJson', '');
    return;
end
assert(ismember(armMode, {'all', 'standard', 'balanced'}), 'PipeCompare:Shadow', 'arm must be all, standard or balanced.');
if ~strcmp(mode, 'test')
    assert(ischar(headSha) && numel(headSha) >= 7, 'PipeCompare:Shadow', 'headSha is required for search and merge modes (pass the run commit SHA).');
end
if strcmp(mode, 'test')
    assert(~isempty(testArms), 'PipeCompare:Shadow', 'test mode requires testArms.');
end
if strcmp(mode, 'merge')
    assert(ischar(fragmentsDir) && isfolder(fragmentsDir), 'PipeCompare:Shadow', 'merge mode requires an existing fragmentsDir.');
end
if strcmp(mode, 'search') && ~isempty(checkpointDir)
    error('PipeCompare:Shadow', 'Full shadow searches must not use a checkpoint directory; use micro slice-diagnostic mode for checkpointed runs.');
end
if ~isfolder(outDir), mkdir(outDir); end

% Round 8 seed_list: a non-contiguous dispatch list overrides the
% contiguous range. NaN entries (from an empty env string) are ignored.
if ~isempty(seedList)
    seedList = seedList(isfinite(seedList));
end
if ~isempty(seedList)
    seeds = sort(unique(seedList(:)'));
else
    seeds = firstSeed:firstSeed + nSeeds - 1;
end
% Fresh-seed assertion (design set + locked confirmation ranges).
% Micro mode is exempt by design: it runs design seeds only, and its
% fragments are labelled design-diagnostic and rejected by merge.
if ~strcmp(mode, 'micro')
    assert(all(seeds >= 20031), 'PipeCompare:Shadow', 'Shadow seeds must use the fresh range (>= 20031).');
    assert(~any(ismember(seeds, 40101:40104)), 'PipeCompare:Shadow', 'Shadow seeds must not reuse design seeds 40101-40104.');
    assert(~any(ismember(seeds, 41001:41030)), 'PipeCompare:Shadow', 'Shadow seeds must not touch locked confirmation seeds 41001-41030.');
end

headShort = headSha(1:min(7, numel(headSha)));
if strcmp(mode, 'merge')
    [perSeed, presentSeeds] = readFragmentsOrArms(fragmentsDir, seeds, headSha, workflowRunId);
    assert(~isempty(perSeed), 'PipeCompare:Shadow', 'No seed fragments found in %s.', fragmentsDir);
    runJson = writeRunOutputs(outDir, headSha, headShort, presentSeeds, perSeed, workflowRunId, eeglabVersion, seeds);
    out = struct('mode', mode, 'seeds', presentSeeds, 'perSeed', perSeed, 'runJson', runJson);
    return;
end

if strcmp(mode, 'test')
    assert(numel(testArms) == nSeeds && isequal([testArms.seed], seeds), ...
        'PipeCompare:Shadow', 'testArms seeds must equal firstSeed:firstSeed+nSeeds-1 in order.');
    perSeed = armsToPerSeed(testArms);
elseif strcmp(armMode, 'standard') || strcmp(armMode, 'balanced')
    % Per-search matrix jobs (Round 6): one seed, one search, one arm
    % fragment. A matrix job always runs exactly one seed.
    assert(nSeeds == 1, 'PipeCompare:Shadow', 'arm mode requires exactly one seed per job.');
    s = seeds(1);
    if strcmp(mode, 'micro')
        fastMode = struct('mode', 'slice-diagnostic', 'partialSearch', false, 'slicedValidation', false, 'designSeedDiagnostic', true, 'reducedFamily', reducedFamily, 'reducedTrials', reducedTrials);
        if strcmp(armMode, 'standard')
            armFrag = runArmStandard(s, maxSeconds, headSha, workflowRunId, eeglabVersion, checkpointDir, fastMode);
        else
            armFrag = runArmBalanced(s, maxSeconds, headSha, workflowRunId, eeglabVersion, checkpointDir, fastMode);
        end
        % A slice-diagnostic fragment is never a held-out fragment:
        % its labels must not claim fullyHeldOut.
        armFrag.labels.fullyHeldOut = false;
    elseif ismember(mode, {'slice','resume'})
        fastMode = struct('mode','sliced-validation','partialSearch',false, ...
            'slicedValidation',true,'designSeedDiagnostic',false,'resume',strcmp(mode,'resume'));
        if strcmp(armMode,'standard')
            armFrag = runArmStandard(s,maxSeconds,headSha,workflowRunId,eeglabVersion,checkpointDir,fastMode);
        else
            armFrag = runArmBalanced(s,maxSeconds,headSha,workflowRunId,eeglabVersion,checkpointDir,fastMode);
        end
        armFrag.labels.fullyHeldOut = false;
        armFrag.execution = 'checkpointed-full-family; separate execution, excluded from unsliced paired means';
    elseif strcmp(armMode, 'standard')
        armFrag = runArmStandard(s, maxSeconds, headSha, workflowRunId, eeglabVersion);
    else
        armFrag = runArmBalanced(s, maxSeconds, headSha, workflowRunId, eeglabVersion);
    end
    fragFile = fullfile(outDir, sprintf('shadow_selection_seed_%d_%s.json', s, armMode));
    writeJson(fragFile, armFrag);
    fprintf('[shadow] seed %d arm %s fragment written: %s\n', s, armMode, fragFile);
    perSeed = struct('seed', {}, 'arms', {}, 'guardChanged', {}, 'guardDecision', {});
    out = struct('mode', mode, 'seeds', seeds, 'perSeed', perSeed, 'runJson', '');
    return;
else
    perSeed = struct('seed', {}, 'arms', {}, 'guardChanged', {}, 'guardDecision', {});
    for s = seeds
        entry = runSeed(s, maxSeconds);
        perSeed(end + 1) = entry; %#ok<AGROW>
        fragFile = fullfile(outDir, sprintf('shadow_selection_seed_%d.json', s));
        writeJson(fragFile, fragmentStruct(entry, headSha, workflowRunId, eeglabVersion));
        fprintf('[shadow] seed %d fragment written: %s\n', s, fragFile);
    end
end

runJson = '';
if nSeeds > 1 || strcmp(mode, 'test')
    runJson = writeRunOutputs(outDir, headSha, headShort, seeds, perSeed, workflowRunId, eeglabVersion);
end
out = struct('mode', mode, 'seeds', seeds, 'perSeed', perSeed, 'runJson', runJson);
end

function entry = runSeed(seed, maxSeconds)
% Run one fresh seed: Standard search + Balanced search, then derive the
% guarded pick (stored production recommendation) and the unguarded pick
% (Rank.run re-execution on the same fitted Balanced family).
entry = struct('seed', seed, 'arms', [], 'guardChanged', [], 'guardDecision', '');
arms = struct('name', {}, 'pick', {}, 'selectedKey', {}, 'mseGainCorrected', {}, 'mseRaw', {}, 'assessable', {}, 'reason', {});

[EEG, truth] = bench_synth(seed);
setBase(EEG);

% Standard arm: its own search (legacy-total). A failed or partial search
% makes only this arm unavailable; it never discards the Balanced arms.
[rStd, stdErr] = runSearch(EEG, seed, 'standard', 'legacy-total', maxSeconds);
% Balanced arm: one search (pooled). Both balanced decisions share it.
[rBal, balErr] = runSearch(EEG, seed, 'balanced', 'pooled', maxSeconds);

% Population truth comes from the Balanced contract; the Standard
% contract must agree with it, otherwise the Standard arm is not
% comparable and is marked unavailable (see contractMatch).
truthVec = [];
contractsAgree = true;
if isempty(balErr)
    truthVec = populationP3(truth, rBal.contract, EEG.srate);
    if ~isempty(stdErr)
        % Standard search failed: nothing to compare contracts against;
        % the Standard arm is unavailable for its own reason below.
        contractsAgree = false;
    elseif ~contractMatch(rStd.contract, rBal.contract)
        contractsAgree = false;
    end
end

% ---- Standard arm ----
stdArm = emptyArm('standard');
if ~isempty(stdErr)
    stdArm.reason = stdErr;
elseif ~isempty(truthVec) && ~contractsAgree
    stdArm.reason = 'contract mismatch with Balanced search';
elseif isempty(rBal) || isempty(truthVec)
    % Balanced search failed: no common truth, so no arm is assessable.
    stdArm.reason = 'no common population truth (Balanced search unavailable)';
else
    pick = storedPick(rStd);
    [mseGC, mseRaw, ok] = armMse(rStd.cands, pick, truthVec);
    stdArm = fillArm(stdArm, pick, rStd, mseGC, mseRaw, ok, 'no recommendation');
end
arms(end + 1) = stdArm; %#ok<AGROW>

% ---- Balanced arms (shared family) ----
[guardArm, unguardArm, entry.guardChanged, entry.guardDecision] = balancedArms(rBal, balErr, truthVec, seed);
arms(end + 1) = guardArm; %#ok<AGROW>
arms(end + 1) = unguardArm; %#ok<AGROW>

entry.arms = arms;
fprintf('[shadow] seed %d: standard pick %s, guarded pick %s, unguarded pick %s\n', ...
    seed, pickText(stdArm.pick), pickText(guardArm.pick), pickText(unguardArm.pick));
end

function [guardArm, unguardArm, guardChanged, guardDecision] = balancedArms(rBal, balErr, truthVec, seed)
% Balanced arms from ONE fitted family (shared by runSeed and by the
% Round 6 Balanced-only job). Logic is exactly the former runSeed
% balanced section, moved here so the two paths cannot drift apart.
guardArm = emptyArm('balancedGuarded');
unguardArm = emptyArm('balancedUnguarded');
guardDecision = ''; guardChanged = [];
if ~isempty(balErr)
    guardArm.reason = balErr;
    unguardArm.reason = balErr;
elseif isempty(truthVec)
    guardArm.reason = 'no population truth';
    unguardArm.reason = 'no population truth';
else
    cands = rBal.cands;
    baseOpts = rBal.ranking.options;
    ub = [];
    if isfield(baseOpts, 'upgradeBaseline'), ub = baseOpts.upgradeBaseline; end
    if ~isempty(ub)
        guardArm.reason = 'explicit upgradeBaseline present; guard contrast undefined';
        unguardArm.reason = 'explicit upgradeBaseline present; guard contrast undefined';
    else
        stored = storedPick(rBal);
        if isempty(stored)
            guardArm.reason = 'no recommendation';
            unguardArm.reason = 'no recommendation (shared family has no guarded pick)';
        else
            optsG = baseOpts; optsG.protectReferenceCleanup = true;
            % Round 8 B1: a reproduction failure must leave an
            % unavailable arm with a reason, not crash the whole job
            % and lose the fragment (and the failure reason with it).
            reproOK = true;
            try
                RG = pipecompare.eval.Rank.run(cands, rBal.ref, optsG, levelsOf(rBal.leaves));
                assert(isscalar(RG.recommended) && RG.recommended == stored, ...
                    'PipeCompare:Shadow', 'Guarded re-execution does not reproduce the stored recommendation (seed %d).', seed);
            catch ME
                reproOK = false;
                pipecompare.utils.log('Shadow: guarded reproduction failed (seed %d): %s', seed, ME.message);
            end
            if reproOK
            optsU = baseOpts; optsU.protectReferenceCleanup = false;
            RU = pipecompare.eval.Rank.run(cands, rBal.ref, optsU, levelsOf(rBal.leaves));
            unguardedPick = RU.recommended;
            if isempty(unguardedPick) || ~isscalar(unguardedPick)
                unguardedPick = [];
            end
            guardDecision = guardDecisionOf(rBal);
            [mseGC, mseRaw, okG] = armMse(cands, stored, truthVec);
            guardArm = fillArm(guardArm, stored, rBal, mseGC, mseRaw, okG, 'no recommendation');
            [mseGC2, mseRaw2, okU] = armMse(cands, unguardedPick, truthVec);
            unguardArm = fillArm(unguardArm, unguardedPick, rBal, mseGC2, mseRaw2, okU, 'no recommendation');
            if isempty(unguardedPick)
                guardChanged = [];
            else
                guardChanged = (stored ~= unguardedPick);
            end
            if ~isempty(unguardedPick) && stored == unguardedPick
                % Round 8 B1: a negative-control violation also becomes
                % an unavailable arm, not a job crash.
                try
                    assert(abs(guardArm.mseGainCorrected - unguardArm.mseGainCorrected) < 1e-12 && ...
                        abs(guardArm.mseRaw - unguardArm.mseRaw) < 1e-12, ...
                        'PipeCompare:Shadow', 'Negative control violated (seed %d): identical picks give different MSEs.', seed);
                catch ME
                    pipecompare.utils.log('Shadow: negative control violated (seed %d): %s', seed, ME.message);
                    guardArm = emptyArm('balancedGuarded'); guardArm.reason = 'negative control violated';
                    unguardArm = emptyArm('balancedUnguarded'); unguardArm.reason = 'negative control violated';
                end
            end
            else
                guardArm.reason = 'guarded reproduction failed';
                unguardArm.reason = 'guarded reproduction failed';
            end
        end
    end
end
end

function frag = runArmStandard(seed, maxSeconds, headSha, workflowRunId, eeglabVersion, checkpointDir, fastMode)
if nargin < 6 || isempty(checkpointDir), checkpointDir = ''; end
if nargin < 7 || isempty(fastMode), fastMode = struct('mode', 'full', 'partialSearch', false, 'slicedValidation', false, 'designSeedDiagnostic', false); end
% Round 6 Standard-only job: one search, one arm fragment. The
% population truth is computed from the Standard contract itself; the
% merge step checks that its fingerprint matches the Balanced
% fingerprint, and only then keeps the Standard arm assessable, so the
% truth vector is identical to the Balanced-contract truth whenever
% the arm is kept (the fingerprints agree on roi, window, baseline
% and condition events, which are exactly the inputs to populationP3).
[EEG, truth] = bench_synth(seed);
setBase(EEG);
[rStd, stdErr, stdTiming] = runSearch(EEG, seed, 'standard', 'legacy-total', maxSeconds, checkpointDir, headSha, fastMode);
stdArm = emptyArm('standard');
fp = struct();
if ~isempty(stdErr)
    stdArm.reason = stdErr;
elseif isempty(rStd)
    stdArm.reason = 'search unavailable';
else
    fp = contractFingerprint(rStd.contract);
    truthVec = populationP3(truth, rStd.contract, EEG.srate);
    pick = storedPick(rStd);
    [mseGC, mseRaw, ok] = armMse(rStd.cands, pick, truthVec);
    stdArm = fillArm(stdArm, pick, rStd, mseGC, mseRaw, ok, 'no recommendation');
end
if ~isempty(stdErr) && contains(stdErr, 'slice interrupted'), fastMode.partialSearch = true; end
frag = armFragmentStruct(seed, stdArm, [], '', fp, headSha, workflowRunId, eeglabVersion, stdTiming, fastMode);
fprintf('[shadow] seed %d: standard pick %s\n', seed, pickText(stdArm.pick));
end

function frag = runArmBalanced(seed, maxSeconds, headSha, workflowRunId, eeglabVersion, checkpointDir, fastMode)
if nargin < 6 || isempty(checkpointDir), checkpointDir = ''; end
if nargin < 7 || isempty(fastMode), fastMode = struct('mode', 'full', 'partialSearch', false, 'slicedValidation', false, 'designSeedDiagnostic', false); end
% Round 6 Balanced-only job: one search, both Balanced arms from that
% one fitted family, plus the contract fingerprint for the merge step.
[EEG, truth] = bench_synth(seed);
setBase(EEG);
[rBal, balErr, balTiming] = runSearch(EEG, seed, 'balanced', 'pooled', maxSeconds, checkpointDir, headSha, fastMode);
fp = struct();
if ~isempty(balErr) || isempty(rBal)
    truthVec = [];
else
    fp = contractFingerprint(rBal.contract);
    truthVec = populationP3(truth, rBal.contract, EEG.srate);
end
[guardArm, unguardArm, guardChanged, guardDecision] = balancedArms(rBal, balErr, truthVec, seed);
arms = [guardArm unguardArm];
if ~isempty(balErr) && contains(balErr, 'slice interrupted'), fastMode.partialSearch = true; end
frag = armFragmentStruct(seed, arms, guardChanged, guardDecision, fp, headSha, workflowRunId, eeglabVersion, balTiming, fastMode);
fprintf('[shadow] seed %d: guarded pick %s, unguarded pick %s\n', seed, pickText(guardArm.pick), pickText(unguardArm.pick));
end

function fp = contractFingerprint(c)
% The contract fields that populationP3 and contractMatch compare:
% first component roi and window, baseline, and condition events.
% Stored in arm fragments so the merge step can check Standard vs
% Balanced comparability without re-running searches. Uses arrayfun
% (not cellfun): c.conditions is a struct array, and cellfun on a
% struct array errors in MATLAB.
fp = struct('roi', {{}}, 'window', [], 'baseline', [], 'conditionEvents', {{}});
try
    s = c.components(1);
    fp.roi = cellstr(s.roi);
    fp.window = s.window;
    fp.baseline = c.baseline;
    fp.conditionEvents = arrayfun(@(x) x.events, c.conditions, 'UniformOutput', false);
catch
    fp = struct();
end
end

function frag = armFragmentStruct(seed, arms, guardChanged, guardDecision, fp, headSha, workflowRunId, eeglabVersion, timing, fastMode)
% One per-search arm fragment. 'arms' is a single standard arm or the
% two balanced arms; 'contract' is the search's fingerprint (empty
% struct when the search itself failed, in which case the arm already
% carries its failure reason).
% Round 9 timing (report-only): frag.timing is the search timing from
% searchTiming (wallSeconds + per-step-type summary). It never changes
% pick or MSE, and old fragments without a timing field are still
% accepted by the merge step and the Python checker (both ignore it).
% The per-seed runSeed path is NOT wired here: it runs two searches
% for one seed, so a single wallSeconds would be ambiguous; only the
% per-search arm fragments (the CI matrix path) carry timing.
if nargin < 9 || isempty(timing), timing = timingEmpty(NaN); end
if nargin < 10 || isempty(fastMode), fastMode = struct('mode', 'full', 'partialSearch', false, 'slicedValidation', false, 'designSeedDiagnostic', false); end
frag = struct('seed', seed, 'arms', arms, 'guardChanged', guardChanged, ...
    'guardDecision', guardDecision, 'contract', fp, 'timing', timing, 'fastMode', fastMode, 'labels', variantLabels(), ...
    'provenance', struct('script', 'benchmarks/pipecompare_shadow_selection.m', 'commit', headSha, ...
    'workflowRunId', workflowRunId, 'matlabVersion', version(), 'eeglabVersion', eeglabVersion, ...
    'eeglabSha', getenv('EEGLAB_SHA'), ...
    'inputArtifactRunId', [], 'inputArtifactNote', 'fresh synthetic seeds; no frozen input artifact'));
end

function [r, err, timing] = runSearch(EEG, seed, recipe, varianceMode, maxSeconds, checkpointDir, headSha, fastMode)
if nargin < 6 || isempty(checkpointDir), checkpointDir = ''; end
if nargin < 7, headSha = ''; end
if nargin < 8, fastMode = struct(); end
% Full scientific definitions are unchanged; a finite budget uses the
% checkpoint helper. Cooperative interruption saves completed leaves/ICA.
% The budget is checked between operations, not inside an individual fit.
% A failed/partial search never supplies a scientific recommendation.
r = []; err = ''; timing = timingEmpty(NaN);
t0 = tic;
try
    if isempty(checkpointDir)
        r = initialSearch(EEG, recipe, varianceMode, checkpointDir, maxSeconds);
    else
        passport = struct('commit',headSha,'seed',seed,'arm',recipe, ...
            'artifactVariance',varianceMode,'eeglabSha',getenv('EEGLAB_SHA'));
        sliceMode = 'start';
        if isfield(fastMode,'resume') && fastMode.resume, sliceMode = 'resume'; end
        [r, sliceState] = pipecompare_checkpoint_search(checkpointDir, passport, ...
            @()initialSearch(EEG,recipe,varianceMode,checkpointDir,maxSeconds), maxSeconds, sliceMode);
        if ~sliceState.complete
            err = 'slice interrupted (partialSearch)';
            timing = timingEmpty(toc(t0));
            return;
        end
    end
    if ~isfield(r, 'notRun') || r.notRun > 0
        err = 'partial search';
        r = [];
    elseif ~isfield(r, 'ranking') || ~isfield(r, 'cands') || ~isfield(r, 'ref') || ~isfield(r, 'contract')
        err = 'search result missing ranking/cands/ref/contract';
        r = [];
    end
catch ME
    if strcmp(ME.identifier,'PipeCompare:Slice') || startsWith(ME.identifier,'PipeCompare:Checkpoint'), rethrow(ME); end
    if strcmp(ME.identifier, 'PipeCompare:Interrupted')
        err = 'slice interrupted (partialSearch)';
    else
        err = ['search failed: ' ME.message];
    end
    r = [];
end
wall = toc(t0);
if isempty(err)
    timing = searchTiming(r, wall);
else
    % A failed/partial search still reports how long it ran, but no
    % per-type summary: there is no complete candidate family to
    % summarise, and inventing one from a partial family would mislead.
    timing = timingEmpty(wall);
end
end

function r = initialSearch(EEG, recipe, varianceMode, checkpointDir, maxSeconds)
[~, ~, r] = pop_pipecompare(EEG, 'measure', 'P3', 'events', {'11','31'}, ...
    'recipe',recipe,'artifactvariance',varianceMode,'reference','average', ...
    'show','off','checkpoint',checkpointDir,'maxseconds',maxSeconds);
end

function t = timingEmpty(wall)
% Timing placeholder: wallSeconds may be NaN (caller unknown) or a
% measured value (failed search); byType is empty and the flag is
% false, so no reader can mistake the absence of timing for zero time.
t = struct('wallSeconds', wall, 'stepSecondsAvailable', false, ...
    'byType', struct('type', {}, 'medianSeconds', {}, 'maxSeconds', {}, 'count', {}));
end

function t = searchTiming(r, wall)
% Per-step-type summary over ALL candidates of one completed search:
% for each step type, the median and maximum of the finite per-step
% seconds recorded by Executor.stepFacts (Round 8), plus how many
% candidate-steps contributed. Median/max are robust to the long tail
% of a few slow candidates; a sum is deliberately NOT reported, because
% shared tree prefixes would be counted once per descendant candidate.
% stepSecondsAvailable is false when no candidate carries a finite
% step seconds field (old fragments / old Executor), never silently 0.
t = timingEmpty(wall);
if isempty(r) || ~isfield(r, 'cands') || isempty(r.cands), return; end
types = {}; values = {};
for ci = 1:numel(r.cands)
    cand = r.cands(ci);
    if ~isfield(cand, 'steps') || isempty(cand.steps), continue; end
    steps = cand.steps;
    if iscell(steps), list = steps; else, list = num2cell(steps); end
    for si = 1:numel(list)
        st = list{si};
        if ~isstruct(st) || ~isfield(st, 'type') || ~isfield(st, 'seconds'), continue; end
        sec = st.seconds;
        if ~isscalar(sec) || ~isfinite(sec), continue; end
        k = find(strcmp(types, st.type), 1);
        if isempty(k), types{end + 1} = st.type; values{end + 1} = sec; %#ok<AGROW>
        else, values{k}(end + 1) = sec; %#ok<AGROW>
        end
    end
end
if isempty(types), return; end
byType = struct('type', {}, 'medianSeconds', {}, 'maxSeconds', {}, 'count', {});
for k = 1:numel(types)
    v = values{k};
    byType(k) = struct('type', types{k}, 'medianSeconds', median(v), 'maxSeconds', max(v), 'count', numel(v)); %#ok<AGROW>
end
t.stepSecondsAvailable = true;
t.byType = byType;
end

function pick = storedPick(r)
% The search's own recommendation: scalar, or [] when there is none.
pick = [];
if isfield(r, 'ranking') && isfield(r.ranking, 'recommended') && isscalar(r.ranking.recommended)
    pick = r.ranking.recommended;
end
end

function decision = guardDecisionOf(rBal)
% The reference-guard decision string of the single-stratum ranking, or
% '' when the guard is empty/absent. Used for reporting only; the
% guardChanged flag is computed from picks, not from this string.
decision = '';
try
    by = rBal.ranking.byStratum;
    if isscalar(by) && isfield(by, 'upgradeGuard') && ~isempty(by.upgradeGuard) && isfield(by.upgradeGuard, 'decision')
        decision = by.upgradeGuard.decision;
    end
catch
    decision = '';
end
end

function a = emptyArm(name)
a = struct('name', name, 'pick', [], 'selectedKey', '', 'mseGainCorrected', NaN, 'mseRaw', NaN, 'assessable', false, 'reason', '');
end

function a = fillArm(a, pick, r, mseGC, mseRaw, ok, noPickReason)
if isempty(pick) || ~ok
    a.reason = noPickReason;
    if ~isempty(pick) && ~ok
        a.reason = 'selected candidate has no finite estimate';
    end
    return;
end
a.pick = pick;
a.selectedKey = r.cands(pick).key;
a.mseGainCorrected = mseGC;
a.mseRaw = mseRaw;
a.assessable = true;
end

function t = pickText(pick)
if isempty(pick), t = 'none'; else, t = num2str(pick); end
end

function perSeed = armsToPerSeed(testArms)
% Test hook: convert caller-supplied fake arms into per-seed entries using
% the SAME arm-MSE, pairing and output code as production. testArms is a
% struct array with fields: seed, truth, standardCands, standardPick,
% balancedCands, balancedGuardedPick, balancedUnguardedPick, guardDecision.
% The two balanced picks share balancedCands, exactly like production.
perSeed = struct('seed', {}, 'arms', {}, 'guardChanged', {}, 'guardDecision', {});
for i = 1:numel(testArms)
    t = testArms(i);
    arms = struct('name', {}, 'pick', {}, 'selectedKey', {}, 'mseGainCorrected', {}, 'mseRaw', {}, 'assessable', {}, 'reason', {});
    stdArm = emptyArm('standard');
    [mseGC, mseRaw, ok] = armMse(t.standardCands, t.standardPick, t.truth);
    stdArm = fillArmTest(stdArm, t.standardCands, t.standardPick, mseGC, mseRaw, ok);
    arms(end + 1) = stdArm; %#ok<AGROW>
    guardArm = emptyArm('balancedGuarded');
    [mseGC, mseRaw, okG] = armMse(t.balancedCands, t.balancedGuardedPick, t.truth);
    guardArm = fillArmTest(guardArm, t.balancedCands, t.balancedGuardedPick, mseGC, mseRaw, okG);
    arms(end + 1) = guardArm; %#ok<AGROW>
    unguardArm = emptyArm('balancedUnguarded');
    [mseGC2, mseRaw2, okU] = armMse(t.balancedCands, t.balancedUnguardedPick, t.truth);
    unguardArm = fillArmTest(unguardArm, t.balancedCands, t.balancedUnguardedPick, mseGC2, mseRaw2, okU);
    arms(end + 1) = unguardArm; %#ok<AGROW>
    entry = struct('seed', t.seed, 'arms', arms, 'guardChanged', [], 'guardDecision', t.guardDecision);
    if ~isempty(t.balancedGuardedPick) && ~isempty(t.balancedUnguardedPick)
        entry.guardChanged = (t.balancedGuardedPick ~= t.balancedUnguardedPick);
        if t.balancedGuardedPick == t.balancedUnguardedPick
            assert(abs(guardArm.mseGainCorrected - unguardArm.mseGainCorrected) < 1e-12 && ...
                abs(guardArm.mseRaw - unguardArm.mseRaw) < 1e-12, ...
                'PipeCompare:Shadow', 'Negative control violated (test seed %d).', t.seed);
        end
    end
    perSeed(end + 1) = entry; %#ok<AGROW>
end
end

function a = fillArmTest(a, cands, pick, mseGC, mseRaw, ok)
if isempty(pick) || ~ok
    if ~isempty(pick) && ~ok
        a.reason = 'selected candidate has no finite estimate';
    else
        a.reason = 'no recommendation';
    end
    return;
end
a.pick = pick;
if isfield(cands(pick), 'key'), a.selectedKey = cands(pick).key; end
a.mseGainCorrected = mseGC;
a.mseRaw = mseRaw;
a.assessable = true;
end

function [mseGC, mseRaw, ok] = armMse(cands, pick, truth)
% Gain-corrected (primary) and raw (secondary) selected MSE for one pick.
% Mirrors benchmarks/replay_reference_guard.m mse(): objectives(1) is the
% P3 measure objective; gain defaults to 1 when absent.
ok = false; mseGC = NaN; mseRaw = NaN;
if isempty(pick) || ~isscalar(pick) || pick < 1 || pick > numel(cands), return; end
c = cands(pick);
if ~isfield(c, 'm') || isempty(c.m) || ~isfield(c.m, 'objectives') || isempty(c.m.objectives) || numel(c.m.objectives) ~= 1, return; end
est = c.m.objectives(1).estimate(:);
truth = truth(:);
if numel(est) ~= numel(truth) || ~all(isfinite(est)) || ~all(isfinite(truth)), return; end
if isfield(c, 'signal') && ~isempty(c.signal)
    g = pipecompare.utils.fieldOr(c.signal, 'gain', 1);
else
    g = 1;
end
g = g(1);
% Round 8 B1: negative gain is also undefined in Rank (gainOf);
% only a finite positive gain gives a meaningful corrected MSE.
if ~isfinite(g) || g <= 0, return; end
mseGC = mean(((est / g) - truth).^2);
mseRaw = mean((est - truth).^2);
ok = isfinite(mseGC) && isfinite(mseRaw);
end

function tf = contractMatch(cStd, cBal)
% The two searches must measure the same thing, or their MSEs are not
% comparable. Compares the first component (roi, window), baseline and
% condition events only; other contract fields may legitimately differ.
% Round 6 fix: the original code used cellfun on c.conditions, but
% conditions is a struct array, and MATLAB cellfun requires a cell
% array input, so it errored on every call and the try/catch below
% returned false every time (the Round 5 Standard arm was therefore
% unavailable on every seed). arrayfun is the correct iterator for a
% struct array.
tf = false;
try
    sStd = cStd.components(1); sBal = cBal.components(1);
    tf = isequal(cellstr(sStd.roi), cellstr(sBal.roi)) && isequal(sStd.window, sBal.window) && ...
        isequal(cStd.baseline, cBal.baseline) && ...
        isequal(arrayfun(@(x) x.events, cStd.conditions, 'UniformOutput', false), ...
        arrayfun(@(x) x.events, cBal.conditions, 'UniformOutput', false));
catch
    tf = false;
end
end

function want = populationP3(truth, c, fs)
% Copied from pipecompare_benchmark.m populationP3, logic unchanged:
% the P3 actually in the data against the population (not realized)
% trial mean, for each condition.
comp = c.components(1); P = truth.population;
roi = ismember(upper(truth.labels), upper(cellstr(comp.roi)));
eeg = ~startsWith(upper(truth.labels), 'EOG'); topo = truth.topo - mean(truth.topo(eeg));
tw = (ceil(comp.window(1) * fs - 1e-6):floor(comp.window(2) * fs + 1e-6)) / fs;
sigma = sqrt(truth.sigma^2 + P.latencySd^2);
shape = mean(truth.sigma / sigma * exp(-0.5 * ((tw - P.latencyMean) / sigma).^2));
if ~isempty(c.baseline)
    tb = (ceil(c.baseline(1) * fs - 1e-6):floor(c.baseline(2) * fs + 1e-6)) / fs;
    shape = shape - mean(truth.sigma / sigma * exp(-0.5 * ((tb - P.latencyMean) / sigma).^2));
end
want = zeros(1, numel(c.conditions));
for k = 1:numel(want)
    assert(numel(c.conditions(k).events) == 1, 'PipeCompare:Benchmark', 'Population truth needs one event type per condition.');
    assert(any(strcmp(c.conditions(k).events{1}, {'11', '31'})), 'PipeCompare:Benchmark', 'Unknown population event type.');
    if strcmp(c.conditions(k).events{1}, '11'), amp = P.targetAmp; else, amp = P.standardAmp; end
    want(k) = mean(topo(roi)) * amp * shape;
end
end

function setBase(EEG)
% Copied from pipecompare_benchmark.m setBase, logic unchanged: the
% simulated EEG becomes the current EEGLAB dataset before the search.
assignin('base', 'BENCH_TMP', EEG);
evalin('base', ['global ALLCOM; ALLEEG = []; [ALLEEG, EEG, CURRENTSET] = eeg_store([], BENCH_TMP, 0); ' ...
    'EEG.history = BENCH_TMP.history; ALLEEG(CURRENTSET).history = BENCH_TMP.history; clear BENCH_TMP;']);
end

function L = levelsOf(leaves)
% Copied from pipecompare.run.Executor levelsOf, logic unchanged: per
% pipeline, its settings as one {slot, key} row per step it runs, for the
% ranking's factorial model. The shadow re-runs must use the same levels
% as the search, or the factorial shrinkage (and hence the ranking)
% would silently differ from production.
L = cell(numel(leaves), 1);
for i = 1:numel(leaves)
    p = leaves(i).path;
    L{i} = [cellfun(@(in) in.slot, p(:), 'UniformOutput', false), cellfun(@(in) in.key, p(:), 'UniformOutput', false)];
end
end

function f = fragmentStruct(entry, headSha, workflowRunId, eeglabVersion)
% Per-seed fragment: the unit the collect job merges. Labels are the
% variant A constants; the merge step asserts them on every fragment.
f = struct('seed', entry.seed, 'arms', entry.arms, 'guardChanged', entry.guardChanged, ...
    'guardDecision', entry.guardDecision, 'labels', variantLabels(), ...
    'provenance', struct('script', 'benchmarks/pipecompare_shadow_selection.m', 'commit', headSha, ...
    'workflowRunId', workflowRunId, 'matlabVersion', version(), 'eeglabVersion', eeglabVersion, ...
    'eeglabSha', getenv('EEGLAB_SHA'), ...
    'inputArtifactRunId', [], 'inputArtifactNote', 'fresh synthetic seeds; no frozen input artifact'));
end

function labels = variantLabels()
labels = struct('fullyHeldOut', true, 'preprocessingRefitted', false, 'adaptiveLeavesIncluded', false);
end

function [perSeed, presentSeeds] = readFragmentsOrArms(fragmentsDir, seeds, headSha, workflowRunId)
% Merge dispatcher: complete per-seed fragments use the original
% readFragments path (which requires every planned seed); per-search
% arm fragments use combineArmFragments (which allows partial
% coverage). When any arm fragment is present, only arm fragments
% are combined (the Round 6 workflow produces only arm fragments, so
% a mixed run should not occur in production; a complete fragment
% in a mixed directory is ignored rather than silently merged).
files = dir(fullfile(fragmentsDir, '**', 'shadow_selection_seed_*.json'));
hasArm = false;
for i = 1:numel(files)
    nm = files(i).name;
    if endsWith(nm, '_standard.json') || endsWith(nm, '_balanced.json'), hasArm = true; break; end
end
if ~hasArm
    perSeed = readFragments(fragmentsDir, seeds, headSha, workflowRunId);
    presentSeeds = seeds;
    return;
end
perSeed = combineArmFragments(fragmentsDir, seeds, headSha, workflowRunId);
presentSeeds = [perSeed.seed];
end

function assertNotFastFragment(f, nm)
% A micro/sliced/design fragment must never enter a held-out merge.
% This is defense in depth: the labels check below would also reject
% a design fragment (fullyHeldOut false), but an explicit fastMode
% check gives the exact reason instead of a generic label failure.
if isfield(f, 'fastMode') && isfield(f.fastMode, 'mode') && ~strcmp(f.fastMode.mode, 'full')
    error('PipeCompare:Shadow', 'Fragment %s is a fast-mode fragment (mode %s) and cannot be merged into a held-out run.', nm, f.fastMode.mode);
end
if isfield(f, 'fastMode') && isfield(f.fastMode, 'designSeedDiagnostic') && isequal(f.fastMode.designSeedDiagnostic, true)
    error('PipeCompare:Shadow', 'Fragment %s is a design-seed diagnostic and cannot be merged into a held-out run.', nm);
end
end

function perSeed = combineArmFragments(fragmentsDir, seeds, headSha, workflowRunId)
% Combine Round 6 per-search arm fragments into per-seed entries.
% A seed with at least one arm fragment appears in the output; a
% missing arm is unavailable with reason 'arm fragment missing'. If
% both arm fragments are present but their contract fingerprints
% differ (or either fingerprint is empty), the Standard arm is
% unavailable with reason 'contract mismatch with Balanced search'
% and its pick/MSE are cleared, matching the original contractMatch
% semantics. Labels, commit and workflow run are asserted on every
% fragment, exactly as in readFragments.
files = dir(fullfile(fragmentsDir, '**', 'shadow_selection_seed_*.json'));
bySeed = struct();
for i = 1:numel(files)
    nm = files(i).name;
    isStd = endsWith(nm, '_standard.json');
    isBal = endsWith(nm, '_balanced.json');
    if ~isStd && ~isBal, continue; end
    f = jsondecode(fileread(fullfile(files(i).folder, nm)));
    assertNotFastFragment(f, nm);
    % Round 8 B1: a fragment whose filename seed disagrees with its
    % content seed must fail, not be silently relabelled by content.
    tok = regexp(nm, '^shadow_selection_seed_(\d+)_(standard|balanced)\.json$', 'tokens');
    assert(~isempty(tok) && str2double(tok{1}{1}) == f.seed, 'PipeCompare:Shadow', 'Fragment %s filename seed disagrees with content seed %d.', nm, f.seed);
    assert(isfield(f, 'seed') && isfield(f, 'arms') && isfield(f, 'labels'), 'PipeCompare:Shadow', 'Fragment %s is missing seed/arms/labels.', nm);
    assert(isfield(f.labels, 'fullyHeldOut') && isfield(f.labels, 'preprocessingRefitted') && isfield(f.labels, 'adaptiveLeavesIncluded'), ...
        'PipeCompare:Shadow', 'Fragment %s has incomplete labels.', nm);
    assert(isequal(f.labels.fullyHeldOut, true) && isequal(f.labels.preprocessingRefitted, false) && isequal(f.labels.adaptiveLeavesIncluded, false), ...
        'PipeCompare:Shadow', 'Fragment %s does not carry variant A labels.', nm);
    if isfield(f, 'provenance') && isfield(f.provenance, 'commit') && ~isempty(f.provenance.commit)
        assert(strcmp(f.provenance.commit, headSha), 'PipeCompare:Shadow', 'Fragment %s is from a different commit.', nm);
    end
    if ~isempty(workflowRunId) && isfield(f, 'provenance') && isfield(f.provenance, 'workflowRunId') && ~isempty(f.provenance.workflowRunId)
        assert(strcmp(f.provenance.workflowRunId, workflowRunId), 'PipeCompare:Shadow', 'Fragment %s is from a different workflow run.', nm);
    end
    key = sprintf('s%d', f.seed);
    if ~isfield(bySeed, key), bySeed.(key) = struct('standard', [], 'balanced', []); end
    % Round 8 B1: duplicate fragments for one seed-arm must fail,
    % not silently overwrite each other.
    if isStd
        assert(isempty(bySeed.(key).standard), 'PipeCompare:Shadow', 'Duplicate standard fragment for seed %d.', f.seed);
        bySeed.(key).standard = f;
    else
        assert(isempty(bySeed.(key).balanced), 'PipeCompare:Shadow', 'Duplicate balanced fragment for seed %d.', f.seed);
        bySeed.(key).balanced = f;
    end
end
perSeed = struct('seed', {}, 'arms', {}, 'guardChanged', {}, 'guardDecision', {});
names = {'standard', 'balancedGuarded', 'balancedUnguarded'};
for s = seeds
    key = sprintf('s%d', s);
    if ~isfield(bySeed, key), continue; end
    entry = bySeed.(key);
    arms = struct('name', {}, 'pick', {}, 'selectedKey', {}, 'mseGainCorrected', {}, 'mseRaw', {}, 'assessable', {}, 'reason', {});
    stdArm = emptyArm('standard');
    if isempty(entry.standard)
        stdArm.reason = 'arm fragment missing';
    else
        found = armByNameInFragment(entry.standard, 'standard');
        stdArm = armFromFragment(stdArm, found);
    end
    guardArm = emptyArm('balancedGuarded');
    unguardArm = emptyArm('balancedUnguarded');
    gc = []; gd = '';
    if isempty(entry.balanced)
        guardArm.reason = 'arm fragment missing';
        unguardArm.reason = 'arm fragment missing';
    else
        guardArm = armFromFragment(guardArm, armByNameInFragment(entry.balanced, 'balancedGuarded'));
        unguardArm = armFromFragment(unguardArm, armByNameInFragment(entry.balanced, 'balancedUnguarded'));
        if isfield(entry.balanced, 'guardChanged'), gc = entry.balanced.guardChanged; end
        if isfield(entry.balanced, 'guardDecision'), gd = entry.balanced.guardDecision; end
    end
    if ~isempty(entry.standard) && ~isempty(entry.balanced)
        if ~fingerprintsMatch(entry.standard, entry.balanced)
            stdArm = emptyArm('standard');
            stdArm.reason = 'contract mismatch with Balanced search';
        end
    end
    arms(end + 1) = stdArm; %#ok<AGROW>
    arms(end + 1) = guardArm; %#ok<AGROW>
    arms(end + 1) = unguardArm; %#ok<AGROW>
    perSeed(end + 1) = struct('seed', s, 'arms', arms, 'guardChanged', gc, 'guardDecision', gd); %#ok<AGROW>
end
end

function arm = armFromFragment(arm, found)
% Copy one fragment arm into a canonical per-seed arm. An empty
% 'found' (arm absent from its fragment) leaves the arm unavailable
% with reason 'arm fragment missing'.
if isempty(found)
    arm.reason = 'arm fragment missing';
    return;
end
arm.pick = found.pick;
if isfield(found, 'selectedKey'), arm.selectedKey = found.selectedKey; end
arm.mseGainCorrected = numOrNan(found.mseGainCorrected);
arm.mseRaw = numOrNan(found.mseRaw);
arm.assessable = logical(found.assessable);
if isfield(found, 'reason'), arm.reason = found.reason; end
end

function found = armByNameInFragment(frag, name)
% Find an arm by name in a fragment, handling the JSON shape where a
% single arm encodes as a scalar struct rather than a 1-element array.
found = [];
arms = frag.arms;
for j = 1:numel(arms)
    if strcmp(arms(j).name, name), found = arms(j); return; end
end
end

function tf = fingerprintsMatch(stdFrag, balFrag)
% True only when both fragments carry a non-empty contract
% fingerprint and the fingerprints are identical on roi, window,
% baseline and condition events. Any missing or differing field is a
% mismatch, matching the conservative original semantics.
tf = false;
try
    if ~isfield(stdFrag, 'contract') || ~isfield(balFrag, 'contract'), return; end
    cs = stdFrag.contract; cb = balFrag.contract;
    if isempty(cs) || isempty(cb), return; end
    if ~isfield(cs, 'roi') || ~isfield(cb, 'roi') || ~isfield(cs, 'window') || ~isfield(cb, 'window') ...
            || ~isfield(cs, 'baseline') || ~isfield(cb, 'baseline') ...
            || ~isfield(cs, 'conditionEvents') || ~isfield(cb, 'conditionEvents'), return; end
    tf = isequal(cellstr(cs.roi), cellstr(cb.roi)) && isequal(cs.window, cb.window) && ...
        isequal(cs.baseline, cb.baseline) && isequal(cs.conditionEvents, cb.conditionEvents);
catch
    tf = false;
end
end

function perSeed = readFragments(fragmentsDir, seeds, headSha, workflowRunId)
% Read per-seed fragment JSONs and rebuild per-seed entries in seed
% order. Fails on any mismatch: wrong count, wrong seed list, labels
% that are not the variant A constants, fragments from a different
% commit or workflow run, or arms that do not match the three canonical
% names. JSON null decodes to [], so empty MSE fields are normalized
% back to NaN here; assessable flags are authoritative.
files = dir(fullfile(fragmentsDir, '**', 'shadow_selection_seed_*.json'));
assert(numel(files) == numel(seeds), 'PipeCompare:Shadow', 'Fragment count %d does not match planned seeds %d.', numel(files), numel(seeds));
bySeed = struct();
for i = 1:numel(files)
    f = jsondecode(fileread(fullfile(files(i).folder, files(i).name)));
    assertNotFastFragment(f, files(i).name);
    assert(isfield(f, 'seed') && isfield(f, 'arms') && isfield(f, 'labels'), 'PipeCompare:Shadow', 'Fragment %s is missing seed/arms/labels.', files(i).name);
    assert(isfield(f.labels, 'fullyHeldOut') && isfield(f.labels, 'preprocessingRefitted') && isfield(f.labels, 'adaptiveLeavesIncluded'), ...
        'PipeCompare:Shadow', 'Fragment %s has incomplete labels.', files(i).name);
    assert(isequal(f.labels.fullyHeldOut, true) && isequal(f.labels.preprocessingRefitted, false) && isequal(f.labels.adaptiveLeavesIncluded, false), ...
        'PipeCompare:Shadow', 'Fragment %s does not carry variant A labels.', files(i).name);
    if isfield(f, 'provenance') && isfield(f.provenance, 'commit') && ~isempty(f.provenance.commit)
        assert(strcmp(f.provenance.commit, headSha), 'PipeCompare:Shadow', 'Fragment %s is from a different commit.', files(i).name);
    end
    if ~isempty(workflowRunId) && isfield(f, 'provenance') && isfield(f.provenance, 'workflowRunId') && ~isempty(f.provenance.workflowRunId)
        assert(strcmp(f.provenance.workflowRunId, workflowRunId), 'PipeCompare:Shadow', 'Fragment %s is from a different workflow run.', files(i).name);
    end
    key = sprintf('s%d', f.seed);
    bySeed.(key) = f;
end
perSeed = struct('seed', {}, 'arms', {}, 'guardChanged', {}, 'guardDecision', {});
names = {'standard', 'balancedGuarded', 'balancedUnguarded'};
for s = seeds
    key = sprintf('s%d', s);
    assert(isfield(bySeed, key), 'PipeCompare:Shadow', 'Missing fragment for seed %d.', s);
    f = bySeed.(key);
    arms = struct('name', {}, 'pick', {}, 'selectedKey', {}, 'mseGainCorrected', {}, 'mseRaw', {}, 'assessable', {}, 'reason', {});
    for a = 1:3
        found = [];
        for j = 1:numel(f.arms)
            if strcmp(f.arms(j).name, names{a}), found = f.arms(j); break; end
        end
        assert(~isempty(found), 'PipeCompare:Shadow', 'Fragment for seed %d is missing arm %s.', s, names{a});
        arm = emptyArm(names{a});
        arm.pick = found.pick;
        if isfield(found, 'selectedKey'), arm.selectedKey = found.selectedKey; end
        arm.mseGainCorrected = numOrNan(found.mseGainCorrected);
        arm.mseRaw = numOrNan(found.mseRaw);
        arm.assessable = logical(found.assessable);
        if isfield(found, 'reason'), arm.reason = found.reason; end
        arms(end + 1) = arm; %#ok<AGROW>
    end
    gc = [];
    if isfield(f, 'guardChanged'), gc = f.guardChanged; end
    gd = '';
    if isfield(f, 'guardDecision'), gd = f.guardDecision; end
    perSeed(end + 1) = struct('seed', s, 'arms', arms, 'guardChanged', gc, 'guardDecision', gd); %#ok<AGROW>
end
end

function v = numOrNan(x)
% JSON null decodes to []; normalize it back to NaN for MSE fields.
if isempty(x), v = NaN; else, v = x; end
end

function runJson = writeRunOutputs(outDir, headSha, headShort, seeds, perSeed, workflowRunId, eeglabVersion, plannedSeeds)
% Shared statistics + output path for multi-seed search runs, test mode
% and merge mode: one implementation, so the contract cannot drift
% between the MATLAB path and the CI collect path. 'seeds' is the
% list actually present in perSeed; 'plannedSeeds' (merge mode only)
% is the planned list, used for the output filename and for
% 'nSeedsPlanned' when a partial run is published. When omitted,
% planned equals present, which is the test-mode and full-run case.
if nargin < 8 || isempty(plannedSeeds), plannedSeeds = seeds; end
stats = pairedStats(perSeed, numel(plannedSeeds));
firstSeed = plannedSeeds(1); nSeeds = numel(plannedSeeds);
base = sprintf('shadow_selection_%s_%d_%d', headShort, firstSeed, nSeeds);
jsonPath = fullfile(outDir, [base '.json']);
recipes = struct('standard', struct('recipe', 'standard', 'artifactVariance', 'legacy-total', 'protectReferenceCleanup', true), ...
    'balancedGuarded', struct('recipe', 'balanced', 'artifactVariance', 'pooled', 'protectReferenceCleanup', true), ...
    'balancedUnguarded', struct('recipe', 'balanced', 'artifactVariance', 'pooled', 'protectReferenceCleanup', false, ...
    'note', 'same fitted family as balancedGuarded; decision-only contrast'));
run = struct('commit', headSha, 'seeds', seeds, 'firstSeed', firstSeed, 'nSeedsPlanned', nSeeds, ...
    'recipes', recipes, 'perSeed', perSeed, 'pairedDifferences', stats.pairedDifferences, ...
    'pairedMeanDifference', stats.pairedMeanDifference, 'recordBootstrapCI95', stats.recordBootstrapCI95, ...
    'coverage', stats.coverage, 'unavailable', stats.unavailable, ...
    'provenance', struct('script', 'benchmarks/pipecompare_shadow_selection.m', 'commit', headSha, ...
    'workflowRunId', workflowRunId, 'matlabVersion', version(), 'eeglabVersion', eeglabVersion, ...
    'eeglabSha', getenv('EEGLAB_SHA'), ...
    'inputArtifactRunId', [], 'inputArtifactNote', 'fresh synthetic seeds; no frozen input artifact', ...
    'outputPath', jsonPath), ...
    'labels', variantLabels(), ...
    'labelDefinitions', struct('fullyHeldOut', 'every seed is fresh (seed >= 20031, not in design 40101-40104, not in locked confirmation 41001-41030)', ...
    'preprocessingRefitted', 'guarded and unguarded decisions share one fitted Balanced family; nothing is refitted per arm', ...
    'adaptiveLeavesIncluded', 'adaptive leaves are NOT re-derived per arm; candidates generated during the single search run sit in the shared family for both decisions, and a per-arm adaptive-leaf effect is NOT measured here (that is variant B)'), ...
    'frozenCriteria', struct('text', ['Non-inferiority (frozen, docs/VALIDATION.md fitted-family shadow summary): the paired gain-corrected ', ...
    'population-target selected MSE difference (guarded Balanced minus Standard) has point estimate <= 0 AND the 95% ', ...
    'record-bootstrap interval upper bound <= +5% of the Standard mean. Point <= 0 with the upper bound above +5%, or ', ...
    'point > 0 with the interval covering 0, is NOT ESTABLISHED (not a failure). Interval lower bound above +5% is ', ...
    'FAILED. No paired seeds or incomplete output is FAILED. Variant A verdicts are conditional ', ...
    '(adaptiveLeavesIncluded = false) and cannot close P0-5.'], 'toleranceFraction', 0.05, ...
    'source', 'docs/VALIDATION.md fitted-family shadow summary'), ...
    'primaryVerdict', stats.verdicts.guardedMinusStandard.verdict, ...
    'verdicts', stats.verdicts);
writeJson(jsonPath, run);
% Per-seed CSV (design section 4.4 columns).
% Round 8 B1: CSV readers must be able to see which pipeline won.
rows = cell(0, 7);
for i = 1:numel(perSeed)
    for a = 1:3
        arm = perSeed(i).arms(a);
        rows(end + 1, :) = {perSeed(i).seed, arm.name, arm.mseGainCorrected, arm.mseRaw, double(arm.assessable), arm.pick, arm.selectedKey}; %#ok<AGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'seed', 'arm', 'selectedMSE_gainCorrected', 'selectedMSE_raw', 'assessable', 'pick', 'selectedKey'});
writetable(T, fullfile(outDir, [base '.csv']));
% Summary CSV (the design's second option: a separate summary file).
pairs = {'guardedMinusStandard', 'guardedMinusUnguarded'};
srows = cell(0, 8);
for p = 1:2
    pr = pairs{p};
    srows(end + 1, :) = {pr, stats.coverage.(['paired_' pr]), stats.pairedMeanDifference.(pr), ...
        stats.recordBootstrapCI95.(pr).lower, stats.recordBootstrapCI95.(pr).upper, ...
        stats.verdicts.(pr).referenceMean, stats.verdicts.(pr).tolerance, stats.verdicts.(pr).verdict}; %#ok<AGROW>
end
S = cell2table(srows, 'VariableNames', {'pair', 'nPaired', 'meanDifference', 'ciLower', 'ciUpper', 'referenceMean', 'tolerance', 'verdict'});
writetable(S, fullfile(outDir, [base '_summary.csv']));
fprintf('[shadow] run JSON: %s\n', jsonPath);
fprintf('[shadow] primary pair guardedMinusStandard: mean %.6g, CI [%.6g, %.6g], verdict %s\n', ...
    stats.pairedMeanDifference.guardedMinusStandard, stats.recordBootstrapCI95.guardedMinusStandard.lower, ...
    stats.recordBootstrapCI95.guardedMinusStandard.upper, stats.verdicts.guardedMinusStandard.verdict);
runJson = jsonPath;
end

function stats = pairedStats(perSeed, plannedCount)
% Round 8 B1: coverage.planned is the PLANNED count, not the present
% count; coverage.present is the present count. The old code set
% planned = numel(perSeed), which made the checker test tautological.
if nargin < 2 || isempty(plannedCount), plannedCount = numel(perSeed); end
% Paired statistics over seeds. A pair uses only seeds assessable in
% BOTH arms (the intersection); anything else is listed unavailable.
% Bootstrap: 10000 resamples of the paired seeds, rng seeded with 1 in
% this statistics phase (after all searches), percentile interval.
names = {'standard', 'balancedGuarded', 'balancedUnguarded'};
pairDefs = struct('pair', {'guardedMinusStandard', 'guardedMinusUnguarded'}, ...
    'first', {'balancedGuarded', 'balancedGuarded'}, 'second', {'standard', 'balancedUnguarded'});
pairedDifferences = struct(); pairedMeanDifference = struct(); recordBootstrapCI95 = struct();
verdicts = struct(); coverage = struct('planned', plannedCount, 'present', numel(perSeed));
unavailable = struct('seed', {}, 'arm', {}, 'reason', {});
for a = 1:3
    n = 0;
    for i = 1:numel(perSeed)
        arm = armByName(perSeed(i), names{a});
        if arm.assessable, n = n + 1; else, unavailable(end + 1) = struct('seed', perSeed(i).seed, 'arm', names{a}, 'reason', arm.reason); end %#ok<AGROW>
    end
    coverage.(names{a}) = n;
end
for p = 1:numel(pairDefs)
    d = pairDefs(p);
    pairedSeeds = []; diffsGC = []; diffsRaw = [];
    for i = 1:numel(perSeed)
        a1 = armByName(perSeed(i), d.first); a2 = armByName(perSeed(i), d.second);
        if a1.assessable && a2.assessable
            pairedSeeds(end + 1) = perSeed(i).seed; %#ok<AGROW>
            diffsGC(end + 1) = a1.mseGainCorrected - a2.mseGainCorrected; %#ok<AGROW>
            diffsRaw(end + 1) = a1.mseRaw - a2.mseRaw; %#ok<AGROW>
        end
    end
    pairedDifferences.(d.pair) = struct('seeds', pairedSeeds, 'gainCorrected', diffsGC, 'raw', diffsRaw);
    if isempty(diffsGC)
        meanDiff = NaN; ci = [NaN NaN]; reps = [];
        refMean = NaN;
    else
        meanDiff = mean(diffsGC);
        refVals = [];
        for i = 1:numel(perSeed)
            a2 = armByName(perSeed(i), d.second);
            if armByName(perSeed(i), d.first).assessable && a2.assessable, refVals(end + 1) = a2.mseGainCorrected; end %#ok<AGROW>
        end
        refMean = mean(refVals);
        rng(1);
        reps = zeros(1, 10000);
        n = numel(diffsGC);
        for b = 1:10000
            reps(b) = mean(diffsGC(randi(n, n, 1)));
        end
        ci = [shadowPercentile(reps, 2.5), shadowPercentile(reps, 97.5)];
    end
    pairedMeanDifference.(d.pair) = meanDiff;
    recordBootstrapCI95.(d.pair) = struct('lower', ci(1), 'upper', ci(2), 'draws', 10000, 'seed', 1, 'replicates', reps);
    tol = 0.05 * refMean;
    if isempty(diffsGC)
        verdict = 'failed'; reason = 'no paired seeds';
    elseif ci(1) > tol
        verdict = 'failed'; reason = 'interval lower bound above +5% tolerance';
    elseif meanDiff <= 0 && ci(2) <= tol
        verdict = 'established'; reason = '';
    else
        verdict = 'not-established'; reason = 'point or interval does not meet the frozen rule';
    end
    verdicts.(d.pair) = struct('verdict', verdict, 'referenceArm', d.second, 'referenceMean', refMean, 'tolerance', tol, 'reason', reason);
    coverage.(['paired_' d.pair]) = numel(pairedSeeds);
end
stats = struct('pairedDifferences', pairedDifferences, 'pairedMeanDifference', pairedMeanDifference, ...
    'recordBootstrapCI95', recordBootstrapCI95, 'coverage', coverage, 'unavailable', unavailable, 'verdicts', verdicts);
end

function arm = armByName(entry, name)
arm = [];
for a = 1:numel(entry.arms)
    if strcmp(entry.arms(a).name, name), arm = entry.arms(a); return; end
end
error('PipeCompare:Shadow', 'Seed %d is missing arm %s.', entry.seed, name);
end

function v = shadowPercentile(x, p)
% Percentile by linear interpolation, the same rule as MATLAB prctile's
% default, implemented here so the script needs no Statistics Toolbox.
% (The Python checker implements this exact formula.)
x = sort(x(:));
n = numel(x);
if n == 0, v = NaN; return; end
if n == 1, v = x(1); return; end
pos = (n - 1) * p / 100; % 0-based position
lo = floor(pos); hi = ceil(pos);
if lo == hi
    v = x(lo + 1);
else
    v = x(lo + 1) + (pos - lo) * (x(hi + 1) - x(lo + 1));
end
end

function writeJson(file, value)
% Atomic-ish write: temp file then move, mirroring
% benchmarks/pipecompare_optimization_synthetic.m writeJson.
fid = fopen([file '.tmp'], 'w');
assert(fid >= 0, 'PipeCompare:Shadow', 'Cannot write %s.', file);
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', jsonencode(value));
clear cleanup;
movefile([file '.tmp'], file, 'f');
end
