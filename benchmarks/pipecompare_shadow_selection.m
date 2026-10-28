function out = pipecompare_shadow_selection(outDir, firstSeed, nSeeds, varargin)
%PIPECOMPARE_SHADOW_SELECTION Held-out shadow comparison, fitted-family variant A.
%
% For each fresh synthetic seed, run two searches (Standard and Balanced)
% and compare three arms on the SAME held-out seed:
%   standard         Standard recipe, artifactVariance legacy-total,
%                    production guard (protectReferenceCleanup = true).
%   balancedGuarded  Balanced recipe, artifactVariance pooled,
%                    production guard (protectReferenceCleanup = true).
%   balancedUnguarded Balanced recipe, artifactVariance pooled, guard
%                    disabled (protectReferenceCleanup = false), decided
%                    on the SAME fitted Balanced family as balancedGuarded.
%
% Variant A fits each family ONCE per seed. The guarded pick is the
% Balanced search's stored production recommendation. The unguarded pick
% comes from re-running pipecompare.eval.Rank.run on the same fitted
% candidates, reference and options, with only protectReferenceCleanup
% flipped to false. Rank.bootstrap is seeded from opts.seed, and the
% ranking options are taken verbatim from the stored ranking, so the
% guarded re-run must reproduce the stored recommendation exactly; if it
% does not, the script errors instead of silently reporting drifted
% decisions. No EEG, ICA or candidate is refitted per arm.
%
% Labels (variant A; see docs/METHODS.md shadow section for the full
% definitions):
%   fullyHeldOut = true: every seed is fresh: seed >= 20031, not in the
%       design set 40101-40104, not in the locked confirmation range
%       41001-41030. The script asserts this; the workflow prep job
%       asserts it again before any job starts.
%   preprocessingRefitted = false: guarded and unguarded decisions share
%       one fitted Balanced family; nothing is refitted per arm.
%   adaptiveLeavesIncluded = false: adaptive leaves are NOT re-derived
%       per arm. If the single search run generated adaptive candidates,
%       they sit in the shared family for BOTH decisions; a per-arm
%       adaptive-leaf effect is NOT measured here. That is variant B.
%       A variant A verdict is therefore conditional and cannot close
%       P0-5 (see frozenCriteria).
%
% Modes (Name-Value 'mode'):
%   'search' (default): run fresh searches. Writes one per-seed fragment
%       JSON per seed immediately after that seed finishes, and, when
%       nSeeds > 1, the run-level JSON + CSV described below. A matrix
%       job runs one seed (nSeeds = 1) and writes only its fragment; a
%       collect job merges fragments with mode 'merge'.
%   'search' with Name-Value 'arm': run only one search for the seed.
%       'arm', 'standard' runs only the Standard search and writes
%       shadow_selection_seed_<seed>_standard.json; 'arm', 'balanced'
%       runs only the Balanced search and writes
%       shadow_selection_seed_<seed>_balanced.json (both Balanced arms,
%       derived from that one fitted family). 'arm', 'all' (default) is
%       the original two-search runSeed path. The CI workflow uses the
%       per-search arms (see Round 6 amendment below); the per-seed
%       path is kept for local runs and backward compatibility.
%   'merge': read fragment JSONs from 'fragmentsDir' and write the
%       run-level JSON + CSV. Accepts complete per-seed fragments
%       (shadow_selection_seed_<seed>.json) and/or per-search arm
%       fragments (shadow_selection_seed_<seed>_standard.json and
%       shadow_selection_seed_<seed>_balanced.json). Arm fragments for
%       one seed are combined into one per-seed entry; a missing arm
%       is unavailable with reason 'arm fragment missing'. If both arm
%       fragments are present but their contract fingerprints differ,
%       the Standard arm is unavailable with reason 'contract mismatch
%       with Balanced search' (the original contractMatch semantics).
%       A seed with no fragment at all is absent from the run output;
%       the run JSON's 'seeds' field lists only the seeds actually
%       present, while 'nSeedsPlanned' records the planned count, so a
%       partial run can still be published and checked instead of
%       being discarded. The collect job must therefore run even when
%       search jobs fail (workflow 'if: always()').
%   'test': unit-test hook. No EEG, no search: the caller supplies
%       'testArms' (see below) and the script runs only the arm-MSE,
%       pairing, bootstrap, negative-control and output code. The CI
%       workflow never sets this mode.
%
% Round 6 amendment (2026-10-28, supersedes the execution description
% below; statistics, seeds, arm definitions and frozen criteria are
% UNCHANGED). Round 5 (head 23f3010) failed for two independent
% reasons: (a) one matrix job ran both searches for a seed and 7 of 8
% seeds exceeded the 330-minute step timeout, so no fragment was
% produced; (b) contractMatch used cellfun on a struct array
% (c.conditions), which MATLAB rejects with 'Input expected to be a
% cell array, was struct instead', and its try/catch swallowed the
% error and returned false every time, so the Standard arm was
% unavailable on every seed even when both searches succeeded (the
% only completed seed, 20033, shows exactly that reason). Round 6
% fixes (b) by using arrayfun in contractMatch, and fixes (a) by
% splitting each seed into two independent matrix jobs (Standard-only
% and Balanced-only, each with a 360-minute job timeout and a
% 350-minute step timeout, fragment upload with 'if: always()').
% Raising 330 to 360 alone was rejected: it adds only 30 minutes
% while seed 20031 was still searching at 330 minutes with an unknown
% remainder. Shrinking the search scope was rejected: it would change
% the arm definitions and break comparability with production search.
% The Round 5 run is void as paired evidence for the new run: its
% 20033 fragment is diagnostic only and is never counted in the new
% run's paired evidence set. If a single search still exceeds
% 360 minutes, GitHub-hosted runners offer no higher job timeout, so
% the next step is search slicing or checkpointing, not another
% timeout raise.
%
% Unavailable semantics (same mode as benchmarks/diagnose_paired_regression.py):
% a seed arm with no recommendation, no finite estimate, or an ambiguous
% target is assessable = false and is excluded from every pair that needs
% it; it is listed in 'unavailable' with a reason, never silently folded
% into a mean. A paired comparison uses only seeds assessable in BOTH
% arms (the intersection). If a pair has no paired seeds, its verdict is
% 'failed' (incomplete output), not 'not-established'.
%
% Negative control: when the guarded and unguarded picks are identical,
% their MSE difference must be exactly 0 (same pick, same family, same
% truth). A non-zero difference means the scoring code is broken, so the
% script errors instead of reporting it.
%
% Selected MSE definitions (per condition, then averaged):
%   gain-corrected (primary): mean((estimate/gain(1) - populationTruth).^2)
%   raw (secondary):          mean((estimate - populationTruth).^2)
% where estimate = cands(pick).m.objectives(1).estimate and gain comes
% from cands(pick).signal.gain (default 1). objectives(1) is the P3
% measure objective of these single-measure searches; if a candidate
% carries anything other than exactly one objective, the arm is
% unavailable ('multiple objectives; target ambiguous').
% populationTruth is the population P3 from bench_synth truth, computed
% with the same formula as pipecompare_benchmark.m populationP3 (copied
% below, logic unchanged). The Standard and Balanced contracts must agree
% on the first component (roi, window), baseline and condition events;
% otherwise the Standard arm is unavailable ('contract mismatch').
%
% Frozen criteria (three-layer-validation-plan.md section 3.1, quoted in
% the output JSON under frozenCriteria): non-inferiority means the paired
% gain-corrected mean difference (guarded Balanced minus Standard) has
% point estimate <= 0 AND the 95% record-bootstrap interval upper bound
% <= +5% of the Standard mean. Point <= 0 with the upper bound above +5%,
% or point > 0 with the interval covering 0, is NOT ESTABLISHED (not a
% failure). Interval lower bound above +5% is FAILED. Variant A verdicts
% are conditional (adaptiveLeavesIncluded = false).
%
% Bootstrap: 10000 resamples of the paired seeds, MATLAB rng seeded with
% 1 in the statistics phase (after all searches, so it cannot disturb
% search randomness), percentile interval by linear interpolation (the
% same rule as MATLAB prctile's default, implemented here so the script
% needs no Statistics Toolbox). The bootstrap replicates are stored in
% the JSON so an independent checker can recompute the interval exactly
% without re-running any random draw.
%
% Checkpoints: none. Like the benchmark matrix pattern (one record per
% job, job timeout as the bound), searches run without checkpoint dirs,
% so no checkpoint identity question arises in variant A. (Whether a
% guard switch changes ICA checkpoint identity matters only for variant
% B, which runs two searches per seed; it must be settled there.)
%
% Workflow defaults (see .github/workflows/shadow-selection.yml):
% first_seed 20031, 8 seeds, MATLAB R2026a + EEGLAB 2026.0.0. Round 5
% used one per-seed job (340 minutes, MATLAB step 330 minutes). Round 6
% uses one per-search job (360 minutes, MATLAB step 350 minutes):
% Standard-only and Balanced-only jobs per seed, each uploading its
% arm fragment with 'if: always()'.
%
% This script never changes production ranking code: it only calls the
% production search (pop_pipecompare) and the production ranker
% (Rank.run). It never writes to the confirmation ranges and never
% touches the gate or legacy artifacts.

% ---- argument parsing (manual Name-Value pairs, no toolbox needed) ----
if nargin < 2 || isempty(firstSeed), firstSeed = 20031; end
if nargin < 3 || isempty(nSeeds), nSeeds = 8; end
mode = 'search'; headSha = ''; workflowRunId = ''; maxSeconds = Inf;
checkpointDir = ''; sliceSeconds = Inf; seedList = [];
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
assert(ismember(mode, {'search', 'merge', 'test', 'resume'}), 'PipeCompare:Shadow', 'mode must be search, merge, test or resume.');
% Round 8 slicing: 'resume' is a design-level mode. It is NOT a silent
% continuation of a timed-out search: the caller must pass a checkpoint
% directory created by a sliced run, and Executor.resume re-derives
% adaptive leaves (see Executor.resume comments). Until a sliced run
% is proven in MATLAB, resume must not be used to complete Round 6
% evidence; see round8-complete/MANIFEST.md.
if strcmp(mode, 'resume')
    error('PipeCompare:Shadow', 'Resume mode is design-only in Round 8: checkpoint slicing has not yet been proven in MATLAB, so it must not complete held-out evidence.');
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
assert(all(seeds >= 20031), 'PipeCompare:Shadow', 'Shadow seeds must use the fresh range (>= 20031).');
assert(~any(ismember(seeds, 40101:40104)), 'PipeCompare:Shadow', 'Shadow seeds must not reuse design seeds 40101-40104.');
assert(~any(ismember(seeds, 41001:41030)), 'PipeCompare:Shadow', 'Shadow seeds must not touch locked confirmation seeds 41001-41030.');

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
    if strcmp(armMode, 'standard')
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

function frag = runArmStandard(seed, maxSeconds, headSha, workflowRunId, eeglabVersion)
% Round 6 Standard-only job: one search, one arm fragment. The
% population truth is computed from the Standard contract itself; the
% merge step checks that its fingerprint matches the Balanced
% fingerprint, and only then keeps the Standard arm assessable, so the
% truth vector is identical to the Balanced-contract truth whenever
% the arm is kept (the fingerprints agree on roi, window, baseline
% and condition events, which are exactly the inputs to populationP3).
[EEG, truth] = bench_synth(seed);
setBase(EEG);
[rStd, stdErr] = runSearch(EEG, seed, 'standard', 'legacy-total', maxSeconds);
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
frag = armFragmentStruct(seed, stdArm, [], '', fp, headSha, workflowRunId, eeglabVersion);
fprintf('[shadow] seed %d: standard pick %s\n', seed, pickText(stdArm.pick));
end

function frag = runArmBalanced(seed, maxSeconds, headSha, workflowRunId, eeglabVersion)
% Round 6 Balanced-only job: one search, both Balanced arms from that
% one fitted family, plus the contract fingerprint for the merge step.
[EEG, truth] = bench_synth(seed);
setBase(EEG);
[rBal, balErr] = runSearch(EEG, seed, 'balanced', 'pooled', maxSeconds);
fp = struct();
if ~isempty(balErr) || isempty(rBal)
    truthVec = [];
else
    fp = contractFingerprint(rBal.contract);
    truthVec = populationP3(truth, rBal.contract, EEG.srate);
end
[guardArm, unguardArm, guardChanged, guardDecision] = balancedArms(rBal, balErr, truthVec, seed);
arms = [guardArm unguardArm];
frag = armFragmentStruct(seed, arms, guardChanged, guardDecision, fp, headSha, workflowRunId, eeglabVersion);
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

function frag = armFragmentStruct(seed, arms, guardChanged, guardDecision, fp, headSha, workflowRunId, eeglabVersion)
% One per-search arm fragment. 'arms' is a single standard arm or the
% two balanced arms; 'contract' is the search's fingerprint (empty
% struct when the search itself failed, in which case the arm already
% carries its failure reason).
% Round 8 timing/diffusion/calibration (report-only): timing and
% diffusion fields are added by the caller via frag.timing /
% frag.diffusion when available; they never change pick or MSE.
frag = struct('seed', seed, 'arms', arms, 'guardChanged', guardChanged, ...
    'guardDecision', guardDecision, 'contract', fp, 'labels', variantLabels(), ...
    'provenance', struct('script', 'benchmarks/pipecompare_shadow_selection.m', 'commit', headSha, ...
    'workflowRunId', workflowRunId, 'matlabVersion', version(), 'eeglabVersion', eeglabVersion, ...
    'eeglabSha', getenv('EEGLAB_SHA'), ...
    'inputArtifactRunId', [], 'inputArtifactNote', 'fresh synthetic seeds; no frozen input artifact'));
end

function [r, err] = runSearch(EEG, seed, recipe, varianceMode, maxSeconds)
% One production search, no checkpoints (benchmark matrix pattern).
% Round 8 slicing (design-only, see above): a sliced run would pass a
% checkpoint directory and a finite maxSeconds here instead of ''/Inf,
% with parallel=false (Executor asserts time slices require serial
% prefix sharing). That change is NOT wired in this bundle because
% Executor.resume re-derives adaptive leaves and has not been proven
% in MATLAB for the shadow arms.
% Any exception makes only the affected recipe unavailable for this seed.
r = []; err = '';
try
    [~, ~, r] = pop_pipecompare(EEG, 'measure', 'P3', 'events', {'11', '31'}, ...
        'recipe', recipe, 'artifactvariance', varianceMode, 'reference', 'average', ...
        'show', 'off', 'checkpoint', '', 'maxseconds', maxSeconds);
    if ~isfield(r, 'notRun') || r.notRun > 0
        err = 'partial search';
        r = [];
    elseif ~isfield(r, 'ranking') || ~isfield(r, 'cands') || ~isfield(r, 'ref') || ~isfield(r, 'contract')
        err = 'search result missing ranking/cands/ref/contract';
        r = [];
    end
catch ME
    err = ['search failed: ' ME.message];
    r = [];
end
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
    'frozenCriteria', struct('text', ['Non-inferiority (frozen, three-layer-validation-plan.md section 3.1): the paired gain-corrected ', ...
    'population-target selected MSE difference (guarded Balanced minus Standard) has point estimate <= 0 AND the 95% ', ...
    'record-bootstrap interval upper bound <= +5% of the Standard mean. Point <= 0 with the upper bound above +5%, or ', ...
    'point > 0 with the interval covering 0, is NOT ESTABLISHED (not a failure). Interval lower bound above +5% is ', ...
    'FAILED. No paired seeds or incomplete output is FAILED. Variant A verdicts are conditional ', ...
    '(adaptiveLeavesIncluded = false) and cannot close P0-5.'], 'toleranceFraction', 0.05, ...
    'source', 'three-layer-validation-plan.md section 3.1'), ...
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
