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
%   'merge': read per-seed fragment JSONs from 'fragmentsDir' and write
%       the run-level JSON + CSV. Fails if the fragment count or the
%       seed list does not match firstSeed/firstSeed+nSeeds-1 exactly.
%   'test': unit-test hook. No EEG, no search: the caller supplies
%       'testArms' (see below) and the script runs only the arm-MSE,
%       pairing, bootstrap, negative-control and output code. The CI
%       workflow never sets this mode.
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
% first_seed 20031, 8 seeds, MATLAB R2026a + EEGLAB 2026.0.0, per-seed
% job timeout 340 minutes (benchmark.yml per-record job bound), MATLAB
% run-step timeout 330 minutes.
%
% This script never changes production ranking code: it only calls the
% production search (pop_pipecompare) and the production ranker
% (Rank.run). It never writes to the confirmation ranges and never
% touches the gate or legacy artifacts.

% ---- argument parsing (manual Name-Value pairs, no toolbox needed) ----
if nargin < 2 || isempty(firstSeed), firstSeed = 20031; end
if nargin < 3 || isempty(nSeeds), nSeeds = 8; end
mode = 'search'; headSha = ''; workflowRunId = ''; maxSeconds = Inf;
testArms = []; fragmentsDir = ''; eeglabVersion = getenv('PIPECOMPARE_EEGLAB_VERSION');
if isempty(eeglabVersion), eeglabVersion = 'unknown'; end
k = 1;
while k <= numel(varargin)
    name = varargin{k};
    assert(ischar(name) && k < numel(varargin), 'PipeCompare:Shadow', 'Name-Value arguments must come in pairs.');
    value = varargin{k + 1};
    switch lower(name)
        case 'mode', mode = value;
        case 'headsha', headSha = value;
        case 'workflowrunid', workflowRunId = value;
        case 'maxseconds', maxSeconds = value;
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
assert(ismember(mode, {'search', 'merge', 'test'}), 'PipeCompare:Shadow', 'mode must be search, merge or test.');
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

seeds = firstSeed:firstSeed + nSeeds - 1;
% Fresh-seed assertion (design set + locked confirmation ranges).
assert(all(seeds >= 20031), 'PipeCompare:Shadow', 'Shadow seeds must use the fresh range (>= 20031).');
assert(~any(ismember(seeds, 40101:40104)), 'PipeCompare:Shadow', 'Shadow seeds must not reuse design seeds 40101-40104.');
assert(~any(ismember(seeds, 41001:41030)), 'PipeCompare:Shadow', 'Shadow seeds must not touch locked confirmation seeds 41001-41030.');

headShort = headSha(1:min(7, numel(headSha)));
if strcmp(mode, 'merge')
    perSeed = readFragments(fragmentsDir, seeds, headSha, workflowRunId);
    runJson = writeRunOutputs(outDir, headSha, headShort, seeds, perSeed, workflowRunId, eeglabVersion);
    out = struct('mode', mode, 'seeds', seeds, 'perSeed', perSeed, 'runJson', runJson);
    return;
end

if strcmp(mode, 'test')
    assert(numel(testArms) == nSeeds && isequal([testArms.seed], seeds), ...
        'PipeCompare:Shadow', 'testArms seeds must equal firstSeed:firstSeed+nSeeds-1 in order.');
    perSeed = armsToPerSeed(testArms);
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
guardArm = emptyArm('balancedGuarded');
unguardArm = emptyArm('balancedUnguarded');
guardDecision = '';
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
        % An explicit upgrade baseline replaces the reference guard with
        % guardUpgrade, so a protect=false re-run would NOT be an
        % unguarded decision. Refuse to invent semantics: both balanced
        % arms are unavailable for this seed. (The benchmark search path
        % never sets upgradeBaseline; this branch is defensive.)
        guardArm.reason = 'explicit upgradeBaseline present; guard contrast undefined';
        unguardArm.reason = 'explicit upgradeBaseline present; guard contrast undefined';
    else
        stored = storedPick(rBal);
        if isempty(stored)
            guardArm.reason = 'no recommendation';
            unguardArm.reason = 'no recommendation (shared family has no guarded pick)';
        else
            % Fidelity re-run: the guarded re-execution must reproduce
            % the stored production recommendation exactly. Any drift
            % means the re-execution path is wrong, so fail loudly
            % instead of reporting drifted decisions.
            optsG = baseOpts; optsG.protectReferenceCleanup = true;
            RG = pipecompare.eval.Rank.run(cands, rBal.ref, optsG, levelsOf(rBal.leaves));
            assert(isscalar(RG.recommended) && RG.recommended == stored, ...
                'PipeCompare:Shadow', 'Guarded re-execution does not reproduce the stored recommendation (seed %d).', seed);
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
                entry.guardChanged = [];
            else
                entry.guardChanged = (stored ~= unguardedPick);
            end
            entry.guardDecision = guardDecision;
            % Negative control: identical picks must give identical MSEs.
            if ~isempty(unguardedPick) && stored == unguardedPick
                assert(abs(guardArm.mseGainCorrected - unguardArm.mseGainCorrected) < 1e-12 && ...
                    abs(guardArm.mseRaw - unguardArm.mseRaw) < 1e-12, ...
                    'PipeCompare:Shadow', 'Negative control violated (seed %d): identical picks give different MSEs.', seed);
            end
        end
    end
end
arms(end + 1) = guardArm; %#ok<AGROW>
arms(end + 1) = unguardArm; %#ok<AGROW>

entry.arms = arms;
fprintf('[shadow] seed %d: standard pick %s, guarded pick %s, unguarded pick %s\n', ...
    seed, pickText(stdArm.pick), pickText(guardArm.pick), pickText(unguardArm.pick));
end

function [r, err] = runSearch(EEG, seed, recipe, varianceMode, maxSeconds)
% One production search, no checkpoints (benchmark matrix pattern).
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
if ~isfinite(g) || g == 0, return; end
mseGC = mean(((est / g) - truth).^2);
mseRaw = mean((est - truth).^2);
ok = isfinite(mseGC) && isfinite(mseRaw);
end

function tf = contractMatch(cStd, cBal)
% The two searches must measure the same thing, or their MSEs are not
% comparable. Compares the first component (roi, window), baseline and
% condition events only; other contract fields may legitimately differ.
tf = false;
try
    sStd = cStd.components(1); sBal = cBal.components(1);
    tf = isequal(cellstr(sStd.roi), cellstr(sBal.roi)) && isequal(sStd.window, sBal.window) && ...
        isequal(cStd.baseline, cBal.baseline) && ...
        isequal(cellfun(@(x) x.events, cStd.conditions, 'UniformOutput', false), ...
        cellfun(@(x) x.events, cBal.conditions, 'UniformOutput', false));
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
    'inputArtifactRunId', [], 'inputArtifactNote', 'fresh synthetic seeds; no frozen input artifact'));
end

function labels = variantLabels()
labels = struct('fullyHeldOut', true, 'preprocessingRefitted', false, 'adaptiveLeavesIncluded', false);
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

function runJson = writeRunOutputs(outDir, headSha, headShort, seeds, perSeed, workflowRunId, eeglabVersion)
% Shared statistics + output path for multi-seed search runs, test mode
% and merge mode: one implementation, so the contract cannot drift
% between the MATLAB path and the CI collect path.
stats = pairedStats(perSeed);
firstSeed = seeds(1); nSeeds = numel(seeds);
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
rows = cell(0, 5);
for i = 1:numel(perSeed)
    for a = 1:3
        arm = perSeed(i).arms(a);
        rows(end + 1, :) = {perSeed(i).seed, arm.name, arm.mseGainCorrected, arm.mseRaw, double(arm.assessable)}; %#ok<AGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'seed', 'arm', 'selectedMSE_gainCorrected', 'selectedMSE_raw', 'assessable'});
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

function stats = pairedStats(perSeed)
% Paired statistics over seeds. A pair uses only seeds assessable in
% BOTH arms (the intersection); anything else is listed unavailable.
% Bootstrap: 10000 resamples of the paired seeds, rng seeded with 1 in
% this statistics phase (after all searches), percentile interval.
names = {'standard', 'balancedGuarded', 'balancedUnguarded'};
pairDefs = struct('pair', {'guardedMinusStandard', 'guardedMinusUnguarded'}, ...
    'first', {'balancedGuarded', 'balancedGuarded'}, 'second', {'standard', 'balancedUnguarded'});
pairedDifferences = struct(); pairedMeanDifference = struct(); recordBootstrapCI95 = struct();
verdicts = struct(); coverage = struct('planned', numel(perSeed));
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
