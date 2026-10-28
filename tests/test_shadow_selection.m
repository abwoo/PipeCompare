function tests = test_shadow_selection
% Tests for benchmarks/pipecompare_shadow_selection.m, variant A contract:
% arm MSE definitions (gain-corrected primary, raw secondary), paired
% differences on arm intersections only, unavailable handling, labels,
% required JSON fields, frozen-criteria reference, and CSV/JSON agreement.
%
% Scope note: these tests exercise the script's statistics and output
% code through its documented test hook ('mode', 'test' with caller
% supplied fake arms). They do NOT run EEG searches and do NOT test the
% guard decision logic itself (covered by test_algorithm_regressions and
% test_methodology). The two balanced picks share one balanced candidate
% array in every test, exactly like production.
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(here); addpath(fullfile(here, '..')); addpath(fullfile(here, '..', 'benchmarks'));
end

function testArmMseGainCorrectedAndRawAndLabels(tc)
% Gain division, raw MSE, pairing, labels and required fields.
truth = [0.5 1.0];
stdCands = mkCands({[1 2], [9 9]}, {[2 2], [1 1]});
balCands = mkCands({[1 2], [3 1]}, {[2 2], [1 1]});
arms = struct('seed', 20031, 'truth', truth, ...
    'standardCands', stdCands, 'standardPick', 1, ...
    'balancedCands', balCands, 'balancedGuardedPick', 2, 'balancedUnguardedPick', 1, ...
    'guardDecision', 'maintain');
outDir = mkOutDir(tc, 'mse');
out = pipecompare_shadow_selection(outDir, 20031, 1, 'mode', 'test', 'testArms', arms, ...
    'headSha', 'testhead1234567890', 'workflowRunId', 'test');
verifyFalse(tc, isempty(out.runJson));
data = jsondecode(fileread(out.runJson));
% Standard pick 1: gain-corrected mean([(1/2-0.5)^2, (2/2-1)^2]) = 0;
% raw mean([(1-0.5)^2, (2-1)^2]) = 0.625.
% Guarded pick 2: gain-corrected mean([(3-0.5)^2, (1-1)^2]) = 3.125; raw 3.125.
% Unguarded pick 1 (shared family): gain-corrected 0; raw 0.625.
ps = data.perSeed(1);
verifyEqual(tc, ps.arms(1).mseGainCorrected, 0, 'AbsTol', 1e-12);
verifyEqual(tc, ps.arms(1).mseRaw, 0.625, 'AbsTol', 1e-12);
verifyEqual(tc, ps.arms(2).mseGainCorrected, 3.125, 'AbsTol', 1e-12);
verifyEqual(tc, ps.arms(2).mseRaw, 3.125, 'AbsTol', 1e-12);
verifyEqual(tc, ps.arms(3).mseGainCorrected, 0, 'AbsTol', 1e-12);
verifyTrue(tc, ps.guardChanged);
verifyEqual(tc, data.pairedDifferences.guardedMinusStandard.gainCorrected, 3.125, 'AbsTol', 1e-12);
verifyEqual(tc, data.pairedDifferences.guardedMinusUnguarded.gainCorrected, 3.125, 'AbsTol', 1e-12);
verifyEqual(tc, data.pairedDifferences.guardedMinusStandard.raw, 2.5, 'AbsTol', 1e-12);
verifyEqual(tc, data.pairedMeanDifference.guardedMinusStandard, 3.125, 'AbsTol', 1e-12);
% Labels are exactly the variant A constants.
verifyTrue(tc, data.labels.fullyHeldOut);
verifyFalse(tc, data.labels.preprocessingRefitted);
verifyFalse(tc, data.labels.adaptiveLeavesIncluded);
% Single-seed bootstrap CI is a point; with a zero Standard mean the
% +5% tolerance is 0, so a positive difference fails by the frozen rule.
% (This documents the rule edge; it is not a verdict on the guard.)
verifyEqual(tc, data.verdicts.guardedMinusStandard.verdict, 'failed');
verifyEqual(tc, data.coverage.planned, 1);
verifyEqual(tc, data.coverage.standard, 1);
verifyEqual(tc, data.coverage.paired_guardedMinusStandard, 1);
end

function testPairedIntersectionUnavailableAndVerdicts(tc)
% Pairing uses intersections only; single-arm-only seeds are listed
% unavailable, never folded into a mean. Verdicts follow the frozen rule.
truth = [0 0];
% Seed 20031: standard 0.5, guarded 0, unguarded 0.5.
a31 = struct('seed', 20031, 'truth', truth, ...
    'standardCands', mkCands({[1 0]}, {[1 1]}), 'standardPick', 1, ...
    'balancedCands', mkCands({[1 0], [0 0]}, {[1 1], [1 1]}), ...
    'balancedGuardedPick', 2, 'balancedUnguardedPick', 1, 'guardDecision', 'maintain');
% Seed 20032: standard 1.0, guarded 4.0, unguarded unavailable.
a32 = struct('seed', 20032, 'truth', truth, ...
    'standardCands', mkCands({[1 1]}, {[1 1]}), 'standardPick', 1, ...
    'balancedCands', mkCands({[2 2]}, {[1 1]}), ...
    'balancedGuardedPick', 1, 'balancedUnguardedPick', [], 'guardDecision', 'maintain');
% Seed 20033: standard unavailable, guarded 0.5, unguarded 0.5.
a33 = struct('seed', 20033, 'truth', truth, ...
    'standardCands', mkCands({[1 0]}, {[1 1]}), 'standardPick', [], ...
    'balancedCands', mkCands({[1 0]}, {[1 1]}), ...
    'balancedGuardedPick', 1, 'balancedUnguardedPick', 1, 'guardDecision', 'maintain');
arms = [a31 a32 a33];
outDir = mkOutDir(tc, 'paired');
out = pipecompare_shadow_selection(outDir, 20031, 3, 'mode', 'test', 'testArms', arms, ...
    'headSha', 'testhead1234567890', 'workflowRunId', 'test');
data = jsondecode(fileread(out.runJson));
% Guarded-minus-Standard pairs seeds {20031, 20032}: diffs [-0.5, 3.0].
% (:)' normalizes jsondecode column vectors to rows before comparing.
verifyEqual(tc, data.pairedDifferences.guardedMinusStandard.seeds(:)', [20031 20032]);
verifyEqual(tc, data.pairedDifferences.guardedMinusStandard.gainCorrected(:)', [-0.5 3.0], 'AbsTol', 1e-12);
verifyEqual(tc, data.pairedMeanDifference.guardedMinusStandard, 1.25, 'AbsTol', 1e-12);
% Point > 0 but the bootstrap interval covers 0: not-established.
verifyEqual(tc, data.verdicts.guardedMinusStandard.verdict, 'not-established');
verifyTrue(tc, data.recordBootstrapCI95.guardedMinusStandard.lower < 0);
verifyTrue(tc, data.recordBootstrapCI95.guardedMinusStandard.upper > 0);
% Guarded-minus-Unguarded pairs seeds {20031, 20033}: diffs [-0.5, 0].
verifyEqual(tc, data.pairedDifferences.guardedMinusUnguarded.seeds(:)', [20031 20033]);
verifyEqual(tc, data.pairedMeanDifference.guardedMinusUnguarded, -0.25, 'AbsTol', 1e-12);
verifyEqual(tc, data.verdicts.guardedMinusUnguarded.verdict, 'established');
% Coverage and unavailable listing.
verifyEqual(tc, data.coverage.standard, 2);
verifyEqual(tc, data.coverage.balancedGuarded, 3);
verifyEqual(tc, data.coverage.balancedUnguarded, 2);
verifyEqual(tc, data.coverage.paired_guardedMinusStandard, 2);
verifyEqual(tc, data.coverage.paired_guardedMinusUnguarded, 2);
unav = data.unavailable;
verifyEqual(tc, numel(unav), 2);
pairs = strcat({unav.arm}, '@', arrayfun(@num2str, [unav.seed], 'UniformOutput', false));
verifyTrue(tc, any(strcmp(pairs, 'balancedUnguarded@20032')));
verifyTrue(tc, any(strcmp(pairs, 'standard@20033')));
end

function testNegativeControlAndDefaultGain(tc)
% Identical guarded/unguarded picks give a zero difference, and a
% candidate without a signal field uses gain 1.
truth = [0 0];
balCands = mkCandsNoSignal({[2 0], [9 9]});
arms = struct('seed', 20031, 'truth', truth, ...
    'standardCands', mkCands({[1 0]}, {[1 1]}), 'standardPick', 1, ...
    'balancedCands', balCands, ...
    'balancedGuardedPick', 1, 'balancedUnguardedPick', 1, 'guardDecision', 'maintain');
outDir = mkOutDir(tc, 'control');
out = pipecompare_shadow_selection(outDir, 20031, 1, 'mode', 'test', 'testArms', arms, ...
    'headSha', 'testhead1234567890', 'workflowRunId', 'test');
data = jsondecode(fileread(out.runJson));
% Gain defaults to 1: mean([2^2, 0^2]) = 2.
verifyEqual(tc, data.perSeed(1).arms(2).mseGainCorrected, 2, 'AbsTol', 1e-12);
verifyEqual(tc, data.perSeed(1).arms(2).mseRaw, 2, 'AbsTol', 1e-12);
verifyFalse(tc, data.perSeed(1).guardChanged);
verifyEqual(tc, data.pairedDifferences.guardedMinusUnguarded.gainCorrected, 0, 'AbsTol', 1e-12);
end

function testCsvMatchesJsonAndRequiredFields(tc)
% The CSVs agree with the JSON, and every section 4.4 field is present.
truth = [0.5 1.0];
arms = struct('seed', 20031, 'truth', truth, ...
    'standardCands', mkCands({[1 2]}, {[2 2]}), 'standardPick', 1, ...
    'balancedCands', mkCands({[1 2], [3 1]}, {[2 2], [1 1]}), ...
    'balancedGuardedPick', 2, 'balancedUnguardedPick', 1, 'guardDecision', 'maintain');
outDir = mkOutDir(tc, 'csv');
out = pipecompare_shadow_selection(outDir, 20031, 1, 'mode', 'test', 'testArms', arms, ...
    'headSha', 'testhead1234567890', 'workflowRunId', 'test');
data = jsondecode(fileread(out.runJson));
for f = {'commit', 'seeds', 'recipes', 'perSeed', 'pairedDifferences', 'pairedMeanDifference', ...
        'recordBootstrapCI95', 'coverage', 'provenance', 'labels', 'frozenCriteria'}
    verifyTrue(tc, isfield(data, f{1}), sprintf('missing field %s', f{1}));
end
verifyEqual(tc, data.recipes.standard.artifactVariance, 'legacy-total');
verifyEqual(tc, data.recipes.balancedGuarded.artifactVariance, 'pooled');
verifyTrue(tc, data.recipes.balancedGuarded.protectReferenceCleanup);
verifyFalse(tc, data.recipes.balancedUnguarded.protectReferenceCleanup);
verifyTrue(tc, contains(data.frozenCriteria.text, '+5%'));
verifyTrue(tc, isempty(data.provenance.inputArtifactRunId));
base = 'shadow_selection_testhea_20031_1';
T = readtable(fullfile(outDir, [base '.csv']));
verifyEqual(tc, height(T), 3);
verifyEqual(tc, string(T.arm), ["standard"; "balancedGuarded"; "balancedUnguarded"]);
verifyEqual(tc, T.selectedMSE_gainCorrected, [0; 3.125; 0], 'AbsTol', 1e-12);
verifyEqual(tc, T.selectedMSE_raw, [0.625; 3.125; 0.625], 'AbsTol', 1e-12);
verifyEqual(tc, T.assessable, [1; 1; 1]);
S = readtable(fullfile(outDir, [base '_summary.csv']));
verifyEqual(tc, height(S), 2);
verifyEqual(tc, S.meanDifference(1), 3.125, 'AbsTol', 1e-12);
end

function testMergeArmFragmentsMatchingFingerprint(tc)
% Round 6 merge: Standard + Balanced arm fragments with identical
% contract fingerprints combine into one per-seed entry, and the
% paired statistics are computed from the combined arms.
fragDir = mkOutDir(tc, 'merge_match');
stdFrag = mkArmFrag(20031, mkArm('standard', 1, 1.0), [], '', mkFingerprint([0.3 0.6]));
balFrag = mkArmFrag(20031, [mkArm('balancedGuarded', 1, 0.5) mkArm('balancedUnguarded', 1, 0.5)], false, 'maintain', mkFingerprint([0.3 0.6]));
writeFrag(fragDir, 'shadow_selection_seed_20031_standard.json', stdFrag);
writeFrag(fragDir, 'shadow_selection_seed_20031_balanced.json', balFrag);
outDir = mkOutDir(tc, 'merge_match_out');
out = pipecompare_shadow_selection(outDir, 20031, 1, 'mode', 'merge', 'fragmentsDir', fragDir, ...
    'headSha', 'testhead1234567890', 'workflowRunId', 'test');
verifyFalse(tc, isempty(out.runJson));
data = jsondecode(fileread(out.runJson));
verifyEqual(tc, data.seeds, 20031);
verifyEqual(tc, data.nSeedsPlanned, 1);
ps = data.perSeed(1);
verifyEqual(tc, ps.arms(1).mseGainCorrected, 1.0, 'AbsTol', 1e-12);
verifyEqual(tc, ps.arms(2).mseGainCorrected, 0.5, 'AbsTol', 1e-12);
verifyTrue(tc, ps.arms(1).assessable);
verifyEqual(tc, data.pairedMeanDifference.guardedMinusStandard, -0.5, 'AbsTol', 1e-12);
verifyEqual(tc, data.verdicts.guardedMinusStandard.verdict, 'established');
end

function testMergeArmFragmentsMismatchedFingerprint(tc)
% Round 6 merge: differing contract fingerprints make the Standard
% arm unavailable with the original contract-mismatch reason, while
% the Balanced pair (guarded minus unguarded) is still computed.
fragDir = mkOutDir(tc, 'merge_mismatch');
stdFrag = mkArmFrag(20031, mkArm('standard', 1, 1.0), [], '', mkFingerprint([0.3 0.6]));
balFrag = mkArmFrag(20031, [mkArm('balancedGuarded', 1, 0.5) mkArm('balancedUnguarded', 1, 0.5)], false, 'maintain', mkFingerprint([0.3 0.7]));
writeFrag(fragDir, 'shadow_selection_seed_20031_standard.json', stdFrag);
writeFrag(fragDir, 'shadow_selection_seed_20031_balanced.json', balFrag);
outDir = mkOutDir(tc, 'merge_mismatch_out');
out = pipecompare_shadow_selection(outDir, 20031, 1, 'mode', 'merge', 'fragmentsDir', fragDir, ...
    'headSha', 'testhead1234567890', 'workflowRunId', 'test');
data = jsondecode(fileread(out.runJson));
ps = data.perSeed(1);
verifyFalse(tc, ps.arms(1).assessable);
verifyTrue(tc, strcmp(ps.arms(1).reason, 'contract mismatch with Balanced search'));
verifyEqual(tc, data.coverage.paired_guardedMinusStandard, 0);
verifyEqual(tc, data.verdicts.guardedMinusStandard.verdict, 'failed');
verifyEqual(tc, data.pairedMeanDifference.guardedMinusUnguarded, 0, 'AbsTol', 1e-12);
verifyEqual(tc, data.verdicts.guardedMinusUnguarded.verdict, 'established');
end

function testMergeArmFragmentsPartialOnlyBalanced(tc)
% Round 6 partial publish: only the Balanced fragment for seed
% 20031 is present and seed 20032 has no fragment at all. The run
% JSON lists only the present seed, but its filename and
% nSeedsPlanned record the planned count (2), and the Standard arm
% is unavailable with reason 'arm fragment missing'.
fragDir = mkOutDir(tc, 'merge_partial');
balFrag = mkArmFrag(20031, [mkArm('balancedGuarded', 1, 0.5) mkArm('balancedUnguarded', 1, 0.5)], false, 'maintain', mkFingerprint([0.3 0.6]));
writeFrag(fragDir, 'shadow_selection_seed_20031_balanced.json', balFrag);
outDir = mkOutDir(tc, 'merge_partial_out');
out = pipecompare_shadow_selection(outDir, 20031, 2, 'mode', 'merge', 'fragmentsDir', fragDir, ...
    'headSha', 'testhead1234567890', 'workflowRunId', 'test');
data = jsondecode(fileread(out.runJson));
verifyEqual(tc, data.seeds, 20031);
verifyEqual(tc, data.nSeedsPlanned, 2);
verifyEqual(tc, numel(data.perSeed), 1);
ps = data.perSeed(1);
verifyFalse(tc, ps.arms(1).assessable);
verifyTrue(tc, strcmp(ps.arms(1).reason, 'arm fragment missing'));
verifyEqual(tc, data.coverage.paired_guardedMinusStandard, 0);
verifyEqual(tc, data.pairedMeanDifference.guardedMinusUnguarded, 0, 'AbsTol', 1e-12);
end

function testMergeArmFragmentsFilenameSeedMismatch(tc)
% Round 9: a fragment whose filename seed disagrees with its content
% seed must fail the merge, not be silently relabelled by content.
% (Production merge already asserts this since Round 8 B1; this test
% is the MATLAB cover that was missing.)
fragDir = mkOutDir(tc, 'merge_filename_mismatch');
stdFrag = mkArmFrag(20031, mkArm('standard', 1, 1.0), [], '', mkFingerprint([0.3 0.6]));
balFrag = mkArmFrag(20031, [mkArm('balancedGuarded', 1, 0.5) mkArm('balancedUnguarded', 1, 0.5)], false, 'maintain', mkFingerprint([0.3 0.6]));
writeFrag(fragDir, 'shadow_selection_seed_20032_standard.json', stdFrag);
writeFrag(fragDir, 'shadow_selection_seed_20031_balanced.json', balFrag);
outDir = mkOutDir(tc, 'merge_filename_mismatch_out');
verifyError(tc, @() pipecompare_shadow_selection(outDir, 20031, 1, 'mode', 'merge', 'fragmentsDir', fragDir, ...
    'headSha', 'testhead1234567890', 'workflowRunId', 'test'), 'PipeCompare:Shadow');
end

function testMergeArmFragmentsDuplicateSeedArm(tc)
% Round 9: two fragments for the same seed-arm (here in two
% subdirectories, as the collect job downloads per-job artifacts)
% must fail the merge, not silently overwrite each other.
% (Production merge already asserts this since Round 8 B1; this test
% is the MATLAB cover that was missing.)
fragDir = mkOutDir(tc, 'merge_duplicate');
subDir = fullfile(fragDir, 'job2');
mkdir(subDir);
stdFrag = mkArmFrag(20031, mkArm('standard', 1, 1.0), [], '', mkFingerprint([0.3 0.6]));
balFrag = mkArmFrag(20031, [mkArm('balancedGuarded', 1, 0.5) mkArm('balancedUnguarded', 1, 0.5)], false, 'maintain', mkFingerprint([0.3 0.6]));
writeFrag(fragDir, 'shadow_selection_seed_20031_standard.json', stdFrag);
writeFrag(subDir, 'shadow_selection_seed_20031_standard.json', stdFrag);
writeFrag(fragDir, 'shadow_selection_seed_20031_balanced.json', balFrag);
outDir = mkOutDir(tc, 'merge_duplicate_out');
verifyError(tc, @() pipecompare_shadow_selection(outDir, 20031, 1, 'mode', 'merge', 'fragmentsDir', fragDir, ...
    'headSha', 'testhead1234567890', 'workflowRunId', 'test'), 'PipeCompare:Shadow');
end

function testMergeCsvPickAndSelectedKeyMatchFragments(tc)
% Round 9: in MERGE mode (not just test mode), the per-seed CSV must
% carry the same pick and selectedKey as the merged fragments, with
% distinct keys per arm. The Python checker requires these columns;
% before this test no MATLAB test covered them on the merge path.
fragDir = mkOutDir(tc, 'merge_csv_pick');
stdArm = mkArm('standard', 1, 1.0); stdArm.selectedKey = 'std-key';
guardArm = mkArm('balancedGuarded', 2, 0.5); guardArm.selectedKey = 'guard-key';
unguardArm = mkArm('balancedUnguarded', 2, 0.5); unguardArm.selectedKey = 'guard-key';
stdFrag = mkArmFrag(20031, stdArm, [], '', mkFingerprint([0.3 0.6]));
balFrag = mkArmFrag(20031, [guardArm unguardArm], false, 'maintain', mkFingerprint([0.3 0.6]));
writeFrag(fragDir, 'shadow_selection_seed_20031_standard.json', stdFrag);
writeFrag(fragDir, 'shadow_selection_seed_20031_balanced.json', balFrag);
outDir = mkOutDir(tc, 'merge_csv_pick_out');
out = pipecompare_shadow_selection(outDir, 20031, 1, 'mode', 'merge', 'fragmentsDir', fragDir, ...
    'headSha', 'testhead1234567890', 'workflowRunId', 'test');
verifyFalse(tc, isempty(out.runJson));
data = jsondecode(fileread(out.runJson));
ps = data.perSeed(1);
base = 'shadow_selection_testhea_20031_1';
T = readtable(fullfile(outDir, [base '.csv']));
verifyEqual(tc, height(T), 3);
verifyEqual(tc, string(T.arm), ["standard"; "balancedGuarded"; "balancedUnguarded"]);
verifyEqual(tc, T.pick, [1; 2; 2]);
verifyEqual(tc, string(T.selectedKey), ["std-key"; "guard-key"; "guard-key"]);
% CSV must agree with the merged JSON, arm by arm, not just exist.
for a = 1:3
    verifyEqual(tc, T.pick(a), ps.arms(a).pick);
    verifyTrue(tc, strcmp(char(T.selectedKey(a)), ps.arms(a).selectedKey));
end
end

function a = mkArm(name, pick, mse)
a = struct('name', name, 'pick', pick, 'selectedKey', 'k', 'mseGainCorrected', mse, 'mseRaw', mse, 'assessable', true, 'reason', '');
end

function fp = mkFingerprint(window)
fp = struct('roi', {{'Pz'}}, 'window', window, 'baseline', [-0.2 0], 'conditionEvents', {{ {'11'}, {'31'} }});
end

function frag = mkArmFrag(seed, arms, guardChanged, guardDecision, fp)
frag = struct('seed', seed, 'arms', arms, 'guardChanged', guardChanged, 'guardDecision', guardDecision, ...
    'contract', fp, 'labels', struct('fullyHeldOut', true, 'preprocessingRefitted', false, 'adaptiveLeavesIncluded', false), ...
    'provenance', struct('script', 'benchmarks/pipecompare_shadow_selection.m', 'commit', 'testhead1234567890', ...
    'workflowRunId', 'test', 'matlabVersion', 'test', 'eeglabVersion', 'test', ...
    'inputArtifactRunId', [], 'inputArtifactNote', 'fresh synthetic seeds; no frozen input artifact'));
end

function writeFrag(dir, filename, frag)
fid = fopen(fullfile(dir, filename), 'w');
assert(fid >= 0, 'Cannot write fragment %s.', filename);
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', jsonencode(frag));
clear cleanup;
end

function cands = mkCands(estimates, gains)
% One candidate per cell entry: m.objectives(1).estimate + signal.gain.
n = numel(estimates);
cands = struct('m', {}, 'signal', {});
for i = 1:n
    cands(i).m = struct('objectives', struct('estimate', estimates{i}));
    cands(i).signal = struct('gain', gains{i});
end
end

function cands = mkCandsNoSignal(estimates)
% Candidates with no signal field at all (gain defaults to 1).
n = numel(estimates);
cands = struct('m', {});
for i = 1:n
    cands(i).m = struct('objectives', struct('estimate', estimates{i}));
end
end

function dir = mkOutDir(tc, tag) %#ok<INUSD>
% A fresh output directory per test (removed by the caller is not
% required: tempdir is scratch, and the run files document the test).
dir = fullfile(tempdir, ['shadow_test_' tag]);
if isfolder(dir), rmdir(dir, 's'); end
mkdir(dir);
end
