function tests = test_engine
%TEST_ENGINE End-to-end searches on synthetic EEG with known truth, and
%   engine invariants (prefix sharing, resume, parallel, units, rules).
tests = functiontests(localfunctions);
end

function setupOnce(~)
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
end

function testFilterCheckExcludesOnlyWhatTheFullCheckExcludes(tc)
% The filter check runs nothing for pipelines whose filters alone break a
% limit; without it, the same pipelines are run and excluded by the full
% signal check, and every other result is the same.
[EEG, truth] = nqc_synth();
nqc_setBase(EEG);
c = pipecompare.eval.Contract('conditions', {'target', {'11'}; 'standard', {'31'}}, ...
    'epoch', [-0.2 1.0], 'baseline', [-0.2 0], 'components', {'P3', [0.30 0.50], truth.roi});
p = pipecompare.plan.Plan();
p = p.add('highpass', 'cutoff', {0.1, 1, 2});
p = p.add('lowpass', 'cutoff', {20, 30});
p = p.add('epoch'); p = p.add('baseline');
a = pipecompare.PipeCompare.optimize(p, c);
b = pipecompare.PipeCompare.optimize(p, c, struct('filterCheck', false));
ex = strcmp({a.cands.status}, 'excluded');
verifyTrue(tc, any(ex));
verifyTrue(tc, all(strcmp(b.ranking.table.status(ex), 'rejected')));    % the full check excludes them too
verifyEqual(tc, a.ranking.table.status, b.ranking.table.status);
verifyEqual(tc, a.ranking.table.objective(~ex), b.ranking.table.objective(~ex), 'RelTol', 1e-12);
verifyEqual(tc, a.ranking.recommended, b.ranking.recommended);
end

function testEndToEndRecoversSensibleChoice(tc)
[EEG, truth] = nqc_synth();
nqc_setBase(EEG);
c = pipecompare.eval.Contract('conditions', {'target', {'11'}; 'standard', {'31'}}, ...
    'epoch', [-0.2 1.0], 'baseline', [-0.2 0], 'components', {'P3', [0.30 0.50], truth.roi});
p = pipecompare.plan.Plan();
p = p.add('highpass', 'cutoff', {0.1, 2});
p = p.add('lowpass', 'cutoff', 30);
p = p.add('epoch'); p = p.add('baseline');
p = p.add('reject_threshold', 'uv', {75, 1000});
r = pipecompare.PipeCompare.optimize(p, c);
T = r.ranking.table;
verifyEqual(tc, numel(r.leaves), 4);
verifyEqual(tc, r.report.nNodes, 2 + 2 + 2 + 2 + 4);
% 2 Hz high-pass distorts the P3 -> rejected by the signal check
hp2 = find(contains(r.labels, 'highpass(cutoff=2)'));
verifyTrue(tc, all(strcmp(T.status(hp2), 'rejected')));
% ... before anything runs: its filters alone already do (filter check)
verifyTrue(tc, all(strcmp({r.cands(hp2).status}, 'excluded')));
verifyTrue(tc, all(startsWith(T.reason(hp2), 'filters alone: ')));
% with movement artifacts on ~15% of trials, rejecting at 75 uV beats no rejection
rec = r.ranking.recommended;
verifyTrue(tc, contains(r.labels{rec}, 'highpass(cutoff=0.1)'));
verifyTrue(tc, contains(r.labels{rec}, 'uv=75'));
% every candidate carries its full command history
verifyTrue(tc, any(contains(r.cands(rec).coms, 'pop_eegfiltnew')));
verifyTrue(tc, any(contains(r.cands(rec).coms, 'pop_rejepoch')));
% the live dataset was not modified
cur = evalin('base', 'EEG');
verifyEqual(tc, cur.history, EEG.history);
% adopt: new dataset, full history, identical SME on replay
nBefore = evalin('base', 'numel(ALLEEG)');
pipecompare.PipeCompare.adopt(r);
verifyEqual(tc, evalin('base', 'numel(ALLEEG)'), nBefore + 1);
adopted = evalin('base', 'EEG');
verifyTrue(tc, contains(adopted.history, 'pop_eegthresh'));
verifyTrue(tc, contains(adopted.history, 'pop_rejepoch'));
m = pipecompare.eval.Measure.candidate(adopted, c, r.ref);
verifyEqual(tc, m.composite, r.cands(rec).m.composite, 'RelTol', 1e-9);
end

function testIcaSharedAcrossThresholds(tc)
EEG = nqc_synth(struct('seconds', 240, 'nPerCond', 50));
nqc_setBase(EEG);
c = pipecompare.eval.Contract('conditions', {'target', {'11'}; 'standard', {'31'}}, ...
    'epoch', [-0.2 1.0], 'baseline', [-0.2 0], 'components', {'P3', [0.30 0.50], {'Pz','P3','P4','POz'}});
p = pipecompare.plan.Plan();
p = p.add('highpass', 'cutoff', 0.1); p = p.add('lowpass', 'cutoff', 30);
p = p.add('reref', 'mode', 'average');
p = p.add('ica', 'fitHighpass', 1); p = p.add('icremove', 'threshold', {0.8, 0.9});
p = p.add('epoch'); p = p.add('baseline');
r = pipecompare.PipeCompare.optimize(p, c);
verifyEqual(tc, numel(r.leaves), 2);
verifyEqual(tc, r.report.nNodes, 4 + 2 * 3); % one ICA for both thresholds
verifyTrue(tc, all(strcmp({r.cands.status}, 'ok')));
verifyTrue(tc, any(contains(r.cands(1).coms, 'EEGica = pop_eegfiltnew')));
verifyTrue(tc, any(contains(r.cands(1).coms, 'pop_icflag')));
end

function testNativeCommandStep(tc)
EEG = nqc_synth(struct('seconds', 120, 'nPerCond', 30, 'artifactTrials', 0));
nqc_setBase(EEG);
c = pipecompare.eval.Contract('conditions', {'target', {'11'}; 'standard', {'31'}}, ...
    'epoch', [-0.2 1.0], 'baseline', [-0.2 0], 'components', {'P3', [0.30 0.50], {'Pz','P3','P4'}});
p = pipecompare.plan.Plan();
p = p.addNative('EEG = pop_eegfiltnew(EEG, ''locutoff'',0.1,''hicutoff'',30,''plotfreqz'',1);');
p = p.add('epoch'); p = p.add('baseline');
r = pipecompare.PipeCompare.optimize(p, c);
verifyEqual(tc, r.cands(1).status, 'ok');
verifyTrue(tc, contains(r.cands(1).coms{1}, '''plotfreqz'',0'));
verifyEqual(tc, r.cands(1).signal.source, 'injection');
verifyGreaterThan(tc, r.cands(1).signal.amplitudeError, 0); % the filter's effect on the known signal is measured
end

% ------------------------------------------------- numerical validation

function testPrefixSharingEqualsIndependentRunsWithIca(tc)
% Shared-prefix execution must give exactly what each pipeline gives when
% run alone from the starting dataset - with real EEGLAB processing,
% including ICA (runica with rndreset 'no' is deterministic).
EEG = nqc_synth(struct('seconds', 150, 'nPerCond', 30));
nqc_setBase(EEG);
c = pipecompare.eval.Contract('conditions', {'t', {'11'}; 's', {'31'}}, 'epoch', [-0.2 1], ...
    'baseline', [-0.2 0], 'components', {'P3', [0.3 0.5], {'Pz','P3','P4'}});
p = pipecompare.plan.Plan();
p = p.add('highpass', 'cutoff', {0.1, 0.5}); p = p.add('lowpass', 'cutoff', 30);
p = p.add('reref', 'mode', 'average'); p = p.add('ica', 'fitHighpass', 1);
p = p.add('icremove', 'threshold', 0.8);
p = p.add('epoch'); p = p.add('baseline'); p = p.add('reject_threshold', 'uv', {80, 1000});
r = pipecompare.PipeCompare.optimize(p, c);
[r.cands.ica] = deal({});   % replay fits ICA again instead of reusing the search's
for k = 1:numel(r.leaves)
    E = pipecompare.run.Executor.replay(r, k);
    m = pipecompare.eval.Measure.candidate(E, c, r.ref);
    verifyEqual(tc, m.composite, r.cands(k).m.composite, 'RelTol', 1e-10, r.labels{k});
    h = regexp(E.history, '[^\n]+', 'match');
    verifyEqual(tc, h(end-numel(r.cands(k).coms)+1:end), r.cands(k).coms);
end
end

function testSourceDatasetIsNeverModified(tc)
EEG = nqc_synth(struct('seconds', 90, 'nPerCond', 20));
nqc_setBase(EEG);
before = evalin('base', 'EEG');
c = nqc_c();
p = pipecompare.plan.Plan(); p = p.add('highpass', 'cutoff', {0.1, 0.5}); p = p.add('epoch'); p = p.add('baseline');
pipecompare.PipeCompare.optimize(p, c);
after = evalin('base', 'EEG');
verifyEqual(tc, after.data, before.data);
verifyEqual(tc, after.history, before.history);
verifyEqual(tc, numel(after.event), numel(before.event));
A = evalin('base', 'ALLEEG');
verifyEqual(tc, numel(A), 1);
end

function testResumeEqualsUninterruptedSearch(tc)
EEG = nqc_synth(struct('seconds', 120, 'nPerCond', 25));
nqc_setBase(EEG);
c = nqc_c();
p = pipecompare.plan.Plan();
p = p.add('highpass', 'cutoff', {0.1, 0.5}); p = p.add('lowpass', 'cutoff', {20, 30});
p = p.add('epoch'); p = p.add('baseline'); p = p.add('reject_threshold', 'uv', {80, 1000});
full = pipecompare.PipeCompare.optimize(p, c);
d = tempname; cleanup = onCleanup(@() rmdir(d, 's'));
verifyError(tc, @() pipecompare.PipeCompare.optimize(p, c, struct('checkpoint', d, 'stopAfter', 3)), 'PipeCompare:Interrupted');
verifyEqual(tc, numel(dir(fullfile(d, 'leaf_*.mat'))), 3);
evalin('base', 'EEG.setname = ''changed after the interruption'';');   % the live data may change meanwhile
res = pipecompare.PipeCompare.resume(d);
verifyEqual(tc, {res.cands.status}, {full.cands.status});
for k = 1:numel(full.cands)
    verifyEqual(tc, res.cands(k).m.composite, full.cands(k).m.composite, 'RelTol', 1e-12);
    verifyEqual(tc, res.cands(k).coms, full.cands(k).coms);
end
verifyEqual(tc, res.ranking.recommended, full.ranking.recommended);
verifyEqual(tc, res.rootFingerprint, full.rootFingerprint);   % identity of the starting dataset kept
% resuming with the relative folder name still stores the absolute path
% (it is used to read the starting dataset back from any folder later)
[parent, name] = fileparts(d);
here = pwd; cd(parent); c3 = onCleanup(@() cd(here)); %#ok<NASGU>
res2 = pipecompare.PipeCompare.resume(name);
cd(here);
w = what(d); verifyEqual(tc, res2.options.checkpoint, w(1).path);
res2.root = [];
E = pipecompare.run.Executor.rootOf(res2);
verifyEqual(tc, E.nbchan, EEG.nbchan);
end

function testParallelEqualsSerial(tc)
assumeTrue(tc, license('test', 'Distrib_Computing_Toolbox') && ~isempty(ver('parallel')), 'Parallel Computing Toolbox not available');
EEG = nqc_synth(struct('seconds', 120, 'nPerCond', 25));
nqc_setBase(EEG);
c = nqc_c();
p = pipecompare.plan.Plan();
p = p.add('highpass', 'cutoff', {0.1, 0.5}); p = p.add('lowpass', 'cutoff', {20, 30});
p = p.add('epoch'); p = p.add('baseline'); p = p.add('reject_threshold', 'uv', {80, 1000});
s = pipecompare.PipeCompare.optimize(p, c);
q = pipecompare.PipeCompare.optimize(p, c, struct('parallel', true));
for k = 1:numel(s.cands)
    verifyEqual(tc, q.cands(k).m.composite, s.cands(k).m.composite, 'RelTol', 1e-12);
    verifyEqual(tc, q.cands(k).coms, s.cands(k).coms);
end
end

function testVoltDataAreConvertedNotMisread(tc)
EEG = nqc_synth(struct('seconds', 90, 'nPerCond', 20));
V = EEG; V.data = V.data * 1e-6;
c = nqc_c();
p = pipecompare.plan.Plan(); p = p.add('highpass', 'cutoff', 0.1); p = p.add('epoch'); p = p.add('baseline');
p = p.add('reject_threshold', 'uv', {80, 1000});
nqc_setBase(EEG); a = pipecompare.PipeCompare.optimize(p, c);
nqc_setBase(V);   b = pipecompare.PipeCompare.optimize(p, c, struct('dataUnit', 'V'));
verifyEqual(tc, [b.cands.status], [a.cands.status]);
verifyEqual(tc, arrayfun(@(x) x.m.composite, b.cands), arrayfun(@(x) x.m.composite, a.cands), 'RelTol', 1e-4);
verifyTrue(tc, any(contains(b.rootComs, '1e6')));
s = pipecompare.live.DataState.fromEEG(V);
verifyTrue(tc, any(contains(s.warnings, 'volts')));
end

function testTrialRuleExcludesPracticeBlock(tc)
EEG = nqc_synth(struct('seconds', 120, 'nPerCond', 30));
nqc_setBase(EEG);
c0 = nqc_c();
c1 = pipecompare.eval.Contract('conditions', {'t', {'11'}; 's', {'31'}}, 'epoch', [-0.2 1], ...
    'baseline', [-0.2 0], 'components', {'P3', [0.3 0.5], {'Pz','P3','P4'}}, ...
    'trials', struct('mode', 'time_ranges', 'ranges', [30 Inf]));
p = pipecompare.plan.Plan(); p = p.add('highpass', 'cutoff', 0.1); p = p.add('epoch'); p = p.add('baseline');
a = pipecompare.PipeCompare.optimize(p, c0); b = pipecompare.PipeCompare.optimize(p, c1);
lat = [EEG.event.latency] / EEG.srate;
expected = sum(lat >= 30 & ismember({EEG.event.type}, {'11'}));
verifyEqual(tc, b.ref.n(1), expected);
verifyLessThan(tc, sum(b.ref.n), sum(a.ref.n));
end

function testInjectionUnaffectedByTrialRule(tc)
% Trials excluded by a trial rule still carry the injected signal, so the
% injection check does not read the exclusion as amplitude loss.
EEG = nqc_synth(struct('seconds', 120, 'nPerCond', 30, 'artifactTrials', 0));
nqc_setBase(EEG);
c = pipecompare.eval.Contract('conditions', {'t', {'11'}; 's', {'31'}}, 'epoch', [-0.2 1], ...
    'baseline', [-0.2 0], 'components', {'P3', [0.3 0.5], {'Pz','P3','P4'}}, ...
    'trials', struct('mode', 'time_ranges', 'ranges', [30 Inf]));
p = pipecompare.plan.Plan(); p = p.add('highpass', 'cutoff', 0.1); p = p.add('epoch'); p = p.add('baseline');
p = p.add('reject_threshold', 'uv', 1000);
r = pipecompare.PipeCompare.optimize(p, c);
verifyLessThan(tc, r.cands(1).signal.amplitudeError, 0.05);
verifyEqual(tc, r.ranking.table.status{1}, 'feasible');
end

function testNoStepPipelineIsEvaluated(tc)
% "Do nothing" is a legitimate option: when every slot chooses 'none' the
% starting copy itself is a candidate.
EEG = nqc_synth(struct('seconds', 90, 'nPerCond', 20));
nqc_setBase(nqc_epoched(EEG));
p = pipecompare.plan.Plan(); p = p.addChoice('hp', {'reject_threshold', 'uv', 100}, 'none');
r = pipecompare.PipeCompare.optimize(p, nqc_contract());
verifyEqual(tc, numel(r.leaves), 2);
verifyTrue(tc, all(strcmp({r.cands.status}, 'ok')));
end

function testEveryCatalogStepRunsThroughEeglab(tc)
% Each catalog step executes as a native EEGLAB call, is recorded in the
% candidate's history, and its decisions replay on the injected copy.
EEG = nqc_synth(struct('seconds', 150, 'nPerCond', 30, 'noisyChannels', {{'T7'}}, 'noisyUv', 200));
nqc_setBase(EEG);
p = pipecompare.plan.Plan();
p = p.add('linenoise', 'freq', 50);
p = p.add('highpass', 'cutoff', 0.5);
p = p.add('badchannels', 'measure', 'kurt', 'threshold', 5, 'action', 'remove');
p = p.add('restore');
p = p.add('reref', 'mode', 'average');
p = p.add('epoch'); p = p.add('baseline');
p = p.addChoice('rej', {'reject_jointprob', 'sd', 5}, {'reject_kurtosis', 'sd', 5});
r = pipecompare.PipeCompare.optimize(p, nqc_c());
verifyEqual(tc, {r.cands.status}, {'ok', 'ok'});
coms = strjoin(r.cands(1).coms, newline);
for f = {'pop_eegfiltnew', 'revfilt', 'bad-channel test', 'pop_interp', 'pop_reref', 'pop_epoch', 'pop_rmbase', 'pop_jointprob', 'pop_rejepoch'}
    verifyTrue(tc, contains(coms, f{1}), f{1});
end
verifyTrue(tc, contains(strjoin(r.cands(2).coms, newline), 'pop_rejkurt'));
verifyEmpty(tc, [r.cands.unmatched]);          % every decision was replayed, none re-run
end

function testAsrDecisionsAreReplayedOnTheSignalCopy(tc)
% ASR's window-by-window reconstructions on the real data are recorded
% (verified against EEGLAB's own output) and applied to the injected copy:
% decision-matched, like every other data-driven step.
assumeTrue(tc, exist('pop_clean_rawdata', 'file') == 2, 'clean_rawdata plugin not installed');
assumeTrue(tc, exist('yulewalk', 'file') == 2, 'ASR at 250 Hz needs the Signal Processing Toolbox');
EEG = nqc_synth(struct('seconds', 120, 'nPerCond', 20));
nqc_setBase(EEG);
p = pipecompare.plan.Plan();
p = p.add('highpass', 'cutoff', 1); p = p.add('asr', 'cutoff', {10, 20});
p = p.add('epoch'); p = p.add('baseline');
r = pipecompare.PipeCompare.optimize(p, nqc_c(), struct('filterCheck', false));   % both run (a 1 Hz high-pass alone fails the P3 check)
for k = 1:2
    verifyEqual(tc, r.cands(k).status, 'ok');
    verifyEmpty(tc, r.cands(k).unmatched);                     % decision-matched
    verifyTrue(tc, any(contains(r.cands(k).coms, '''MaxMem'',64')));
    verifyTrue(tc, isfinite(r.cands(k).signal.amplitudeError));
end
% the recorded decisions are ASR's: replaying them on the real data gives
% exactly EEGLAB's output, and they are linear in the data
[~, E] = evalc('pop_eegfiltnew(EEG, ''locutoff'', 1, ''plotfreqz'', 0)');
[~, E1] = evalc('pop_clean_rawdata(E, ''FlatlineCriterion'',''off'',''ChannelCriterion'',''off'',''LineNoiseCriterion'',''off'',''Highpass'',''off'',''BurstCriterion'',10,''WindowCriterion'',''off'',''BurstRejection'',''off'',''Distance'',''Euclidian'',''MaxMem'',64)');
rec = pipecompare.run.AsrRecord.record(E, 10);
verifyEqual(tc, pipecompare.run.AsrRecord.apply(rec, E.data), double(E1.data), 'AbsTol', 1e-9);
verifyGreaterThan(tc, sum(~[rec.updates.trivial]), 0);              % ASR did change windows
A = randn(size(E.data)); B = randn(size(E.data));
verifyEqual(tc, pipecompare.run.AsrRecord.apply(rec, A + 2 * B), ...
    pipecompare.run.AsrRecord.apply(rec, A) + 2 * pipecompare.run.AsrRecord.apply(rec, B), 'AbsTol', 1e-8);
end

function testPeakLatencyObjective(tc)
EEG = nqc_synth(struct('seconds', 200, 'nPerCond', 40, 'p3Jitter', 0.03));
nqc_setBase(EEG);
c = pipecompare.eval.Contract('conditions', {'t', {'11'}}, 'epoch', [-0.2 1], 'baseline', [-0.2 0], ...
    'components', {'P3', [0.25 0.6], {'Pz','P3','P4'}, 'peakLatency'});
p = pipecompare.plan.Plan(); p = p.add('highpass', 'cutoff', 0.1); p = p.add('lowpass', 'cutoff', {10, 30});
p = p.add('epoch'); p = p.add('baseline'); p = p.add('reject_threshold', 'uv', {75, 1000});
r = pipecompare.PipeCompare.optimize(p, c, struct('objective', {{'P3.peakLatency'}}, 'nBootPeakOuter', 100, 'nBootPeak', 300));
verifyEqual(tc, r.ref.units{1}, 'ms');
ok = strcmp({r.cands.status}, 'ok');
verifyTrue(tc, all(ok));
est = arrayfun(@(x) x.m.objectives(1).estimate, r.cands);
verifyLessThan(tc, abs(est - 400), 60);   % true peak at 400 ms
verifyNotEmpty(tc, r.ranking.recommended);
end

function testNamedChannelRepairOptions(tc)
% A broken electrode is either kept, removed or interpolated - only as the
% user allows; the result reports which choice keeps more trials.
EEG = nqc_synth(struct('seconds', 150, 'nPerCond', 30, 'noisyChannels', {{'O1','O2'}}, 'noisyUv', 80, 'artifactTrials', 0, 'blinkRate', 0));
nqc_setBase(EEG);
c = nqc_c();
p = pipecompare.plan.Plan();
p = p.add('highpass', 'cutoff', 0.1);
p = p.addChoice('o1o2', {'channels', 'labels', {'O1','O2'}, 'action', 'interpolate'}, ...
    {'channels', 'labels', {'O1','O2'}, 'action', 'remove'}, 'none');
p = p.add('epoch'); p = p.add('baseline'); p = p.add('reject_threshold', 'uv', 100);
r = pipecompare.PipeCompare.optimize(p, c);
T = r.ranking.table;
keepAll = find(~contains(r.labels, 'channels('));
verifyEqual(tc, T.status{keepAll}, 'rejected');       % noisy O1/O2: the threshold removes every epoch
verifyEqual(tc, T.minRetention(keepAll), 0);
verifyTrue(tc, contains(T.reason{keepAll}, 'retention 0%'));
fixed = find(contains(r.labels, 'channels('));
verifyTrue(tc, all(T.minRetention(fixed) > 0.9));
verifyTrue(tc, all(contains(r.labels([r.ranking.byStratum.recommended]), 'channels(')));
end

function testVerboseNormalShowsCommandsNotEeglabChatter(tc)
EEG = nqc_synth(struct('seconds', 60, 'nPerCond', 10, 'artifactTrials', 0));
nqc_setBase(EEG);
p = pipecompare.plan.Plan(); p = p.add('highpass', 'cutoff', 0.5); p = p.add('epoch'); p = p.add('baseline');
txt = evalc('pipecompare.PipeCompare.optimize(p, nqc_c());');
verifyTrue(tc, contains(txt, 'EEG = pop_eegfiltnew(EEG, ''locutoff'',0.5'));
verifyFalse(tc, contains(txt, 'pop_eegfiltnew() - performing'));   % EEGLAB's own progress lines
verifyTrue(tc, contains(txt, 'Recommended'));
end

function testTrialRuleRestrictsTheOutputToo(tc)
% Audit: with a trial rule only the eligible trials were scored, but the
% adopted dataset still held every epoch. Adopted data and the exported
% script must hold exactly the scored trials.
EEG = nqc_synth(struct('seconds', 120, 'nPerCond', 30, 'artifactTrials', 0));
nqc_setBase(EEG);
c = pipecompare.eval.Contract('conditions', {'t', {'11'}; 's', {'31'}}, 'epoch', [-0.2 1], ...
    'baseline', [-0.2 0], 'components', {'P3', [0.3 0.5], {'Pz','P3','P4'}}, ...
    'trials', struct('mode', 'time_ranges', 'ranges', [30 Inf]));
p = pipecompare.plan.Plan(); p = p.add('highpass', 'cutoff', 0.1); p = p.add('epoch'); p = p.add('baseline');
r = pipecompare.PipeCompare.optimize(p, c);
pipecompare.PipeCompare.adopt(r, 1);
ad = evalin('base', 'EEG');
verifyEqual(tc, ad.trials, sum(r.ref.n));
verifyTrue(tc, contains(ad.history, 'notrial'));
d = tempname; mkdir(d); cleanup = onCleanup(@() rmdir(d, 's')); %#ok<NASGU>
pipecompare.PipeCompare.writeScript(r, 1, fullfile(d, 'nqc_rule_script.m'));
addpath(d); c2 = onCleanup(@() rmpath(d)); %#ok<NASGU>
out = nqc_rule_script(EEG);
verifyEqual(tc, out.trials, sum(r.ref.n));
end

function testEarlierTrialRuleDoesNotPersist(tc)
% Audit: after switching the rule back to 'all', the earlier restriction
% was still applied (it travelled in EEG.etc of the dataset).
EEG = nqc_synth(struct('seconds', 90, 'nPerCond', 20, 'artifactTrials', 0));
p = pipecompare.plan.Plan(); p = p.add('highpass', 'cutoff', 0.1); p = p.add('epoch'); p = p.add('baseline');
nqc_setBase(EEG); a = pipecompare.PipeCompare.optimize(p, nqc_c());
E2 = EEG; E2.etc.pipecompare.eligibleUrevents = [EEG.event(1:4).urevent];   % left by an earlier search
nqc_setBase(E2); b = pipecompare.PipeCompare.optimize(p, nqc_c());
verifyEqual(tc, b.ref.n, a.ref.n);
verifyTrue(tc, any(contains(b.rootComs, 'rmfield(EEG.etc.pipecompare, ''eligibleUrevents'')')));
end

function testVoltConversionScriptScalesIcaLikeTheSearch(tc)
% Audit: the search scaled ICA weights with the data, the script did not.
EEG = nqc_synth(struct('seconds', 90, 'nPerCond', 20));
V = EEG; V.data = V.data * 1e-6;
rng(3); V.icaweights = randn(V.nbchan); V.icasphere = eye(V.nbchan); V.icachansind = 1:V.nbchan;
V.icawinv = []; V.icaact = []; V = eeg_checkset(V);
nqc_setBase(V);
% no epoch step: pop_epoch renormalises the components and would hide the scale
p = pipecompare.plan.Plan(); p = p.add('highpass', 'cutoff', 0.1);
r = pipecompare.PipeCompare.optimize(p, nqc_c(), struct('dataUnit', 'V'));
d = tempname; mkdir(d); cleanup = onCleanup(@() rmdir(d, 's')); %#ok<NASGU>
pipecompare.PipeCompare.writeScript(r, 1, fullfile(d, 'nqc_volt_script.m'));
addpath(d); c2 = onCleanup(@() rmpath(d)); %#ok<NASGU>
out = nqc_volt_script(V);
ref = pipecompare.run.Executor.replay(r, 1);           % the search's own in-memory path
% equal up to the per-component scale EEGLAB may set (option_scaleicarms:
% eeg_checkset rescales each row to RMS microvolts in some EEGLAB versions);
% a missing volts -> microvolts conversion would be a factor of 1e6
ratio = out.icaweights ./ ref.icaweights;
verifyEqual(tc, ratio, repmat(ratio(:, 1), 1, size(ratio, 2)), 'RelTol', 1e-6);
verifyLessThan(tc, max(abs(log10(abs(ratio(:, 1))))), 2);
m = pipecompare.eval.Measure.candidate(out, r.contract, r.ref);
verifyEqual(tc, [m.objectives.agg], [r.cands(1).m.objectives.agg], 'RelTol', 1e-6);
end

function testCheckpointFolderBelongsToOneSearch(tc)
% Audit: an existing root.mat was kept while manifest.mat was replaced, so
% a second search in the same folder could resume on another dataset.
EEG = nqc_synth(struct('seconds', 90, 'nPerCond', 20));
p = pipecompare.plan.Plan(); p = p.add('highpass', 'cutoff', {0.1, 0.5}); p = p.add('epoch'); p = p.add('baseline');
d = tempname; cleanup = onCleanup(@() rmdir(d, 's')); %#ok<NASGU>
nqc_setBase(EEG);
verifyError(tc, @() pipecompare.PipeCompare.optimize(p, nqc_c(), struct('checkpoint', d, 'stopAfter', 1)), 'PipeCompare:Interrupted');
other = nqc_synth(struct('seconds', 90, 'nPerCond', 20, 'seed', 99));
nqc_setBase(other);
verifyError(tc, @() pipecompare.PipeCompare.optimize(p, nqc_c(), struct('checkpoint', d)), 'PipeCompare:Checkpoint');
c3 = nqc_c(); c3.baseline = [-0.1 0];
nqc_setBase(EEG);
verifyError(tc, @() pipecompare.PipeCompare.optimize(p, c3, struct('checkpoint', d)), 'PipeCompare:Checkpoint');
% the folder still resumes its own search, and the same search may start over
res = pipecompare.PipeCompare.resume(d);
verifyTrue(tc, all(strcmp({res.cands.status}, 'ok')));
again = pipecompare.PipeCompare.optimize(p, nqc_c(), struct('checkpoint', d));
verifyEqual(tc, again.identity, res.identity);
% a starting dataset from another search is refused on resume
d2 = tempname; c4 = onCleanup(@() rmdir(d2, 's')); %#ok<NASGU>
nqc_setBase(other); pipecompare.PipeCompare.optimize(p, nqc_c(), struct('checkpoint', d2));
copyfile(fullfile(d2, 'root.mat'), fullfile(d, 'root.mat'));
verifyError(tc, @() pipecompare.PipeCompare.resume(d), 'PipeCompare:Checkpoint');
end

function testAdoptRefusesInfeasibleAndNonReproducingCandidates(tc)
EEG = nqc_synth(struct('seconds', 90, 'nPerCond', 20));
nqc_setBase(EEG);
p = pipecompare.plan.Plan(); p = p.add('highpass', 'cutoff', 0.1); p = p.add('epoch'); p = p.add('baseline');
r = pipecompare.PipeCompare.optimize(p, nqc_c(), struct('minTrials', 1000));
n0 = evalin('base', 'numel(ALLEEG)');
verifyError(tc, @() pipecompare.PipeCompare.adopt(r, 1), 'PipeCompare:Adopt');          % rejected: not adopted silently
verifyEqual(tc, evalin('base', 'numel(ALLEEG)'), n0);
pipecompare.PipeCompare.adopt(r, 1, true);                                          % explicit override
verifyEqual(tc, evalin('base', 'numel(ALLEEG)'), n0 + 1);
nqc_setBase(EEG);
r = pipecompare.PipeCompare.optimize(p, nqc_c());
r.cands(1).m.objectives(1).agg = 2 * r.cands(1).m.objectives(1).agg;      % evaluation no longer matches the data
n0 = evalin('base', 'numel(ALLEEG)');
verifyError(tc, @() pipecompare.PipeCompare.adopt(r, 1), 'PipeCompare:ReplayMismatch');
verifyEqual(tc, evalin('base', 'numel(ALLEEG)'), n0);
end

function testRestoreChannelsRemovedBeforePipeCompare(tc)
% Audit: channels removed in EEGLAB before entering could not be restored
% by the plan (restore demanded a removal inside the plan).
EEG = nqc_synth(struct('seconds', 90, 'nPerCond', 20));
n = EEG.nbchan;
EEG = pop_select(EEG, 'nochannel', {'O1', 'O2'});
assert(numel(EEG.chaninfo.removedchans) >= 2);
nqc_setBase(EEG);
p = pipecompare.plan.Plan(); p = p.add('highpass', 'cutoff', 0.1); p = p.add('restore'); p = p.add('epoch'); p = p.add('baseline');
r = pipecompare.PipeCompare.optimize(p, nqc_c());
verifyEqual(tc, r.cands(1).status, 'ok');
verifyEqual(tc, r.cands(1).interpolatedFraction, 2 / n, 'AbsTol', 1e-12);
verifyTrue(tc, any(contains(r.cands(1).coms, 'preRemoved')));
pipecompare.PipeCompare.adopt(r, 1);
ad = evalin('base', 'EEG');
verifyEqual(tc, ad.nbchan, n);
verifyTrue(tc, all(ismember({'O1', 'O2'}, {ad.chanlocs.labels})));
end

function testNativeRejectionWorkflowKeepsDialogSettingsAndDecisions(tc)
% Mark -> reject captured from EEGLAB's dialogs is one step: every dialog
% setting is kept (asymmetric limits here), both commands are in the
% history, and the signal check applies the epochs the real data rejected
% (decision-matched) instead of re-deciding on the injected copy.
EEG = nqc_synth(struct('seconds', 120, 'nPerCond', 30));
nqc_setBase(EEG);
wf = sprintf(['EEG = pop_eegthresh(EEG,1,[1:32],-60,120,-0.2,0.996,0,0);\n', ...
    'EEG = pop_rejepoch(EEG, EEG.reject.rejthresh, 0);']);
p = pipecompare.plan.Plan(); p = p.add('highpass', 'cutoff', 0.1); p = p.add('epoch'); p = p.add('baseline');
p = p.addNative(wf, 'reject');
r = pipecompare.PipeCompare.optimize(p, nqc_c());
c = r.cands(1);
verifyEqual(tc, c.status, 'ok');
verifyGreaterThan(tc, c.rejectedEpochs, 0);
verifyEmpty(tc, c.unmatched);                                   % decision-matched
verifyLessThan(tc, c.signal.amplitudeError, 0.05);
verifyTrue(tc, any(contains(c.coms, '-60,120')) && any(contains(c.coms, 'pop_rejepoch')));
% the decision equals what the asymmetric limits mark on the processed data
E = pipecompare.run.Executor.replay(r, 1);
verifyEqual(tc, E.trials, 60 - c.rejectedEpochs);
% legality: a rejection workflow on continuous data is refused with a reason
q = pipecompare.plan.Plan(); q = q.addNative(wf, 'reject'); q = q.add('epoch');
verifyError(tc, @() q.enumerate(pipecompare.live.DataState.fromEEG(EEG), nqc_c()), 'PipeCompare:NoLegalPipeline');
end

function testReplayReusesTheSearchIca(tc)
% Adopting does not fit ICA again: the search's decomposition is applied.
EEG = nqc_synth(struct('seconds', 90, 'nPerCond', 20));
nqc_setBase(EEG);
p = pipecompare.plan.Plan(); p = p.add('highpass', 'cutoff', 0.1); p = p.add('ica');
r = pipecompare.PipeCompare.optimize(p, nqc_c());
d = r.cands(1).ica;
verifyNumElements(tc, d, 1);
verifyTrue(tc, contains(d{1}.coms{1}, 'EEGica = pop_runica(EEGica, ''icatype'''));
verifyTrue(tc, contains(d{1}.coms{1}, 'EEGica = pipecompare.run.Steps.cleanForIca(EEGica);'));
W = d{1}.icaweights([2 1 3:end], :);   % not what runica gives: proves it is applied, not refitted
r.cands(1).ica{1}.icaweights = W;
E = pipecompare.run.Executor.replay(r, 1);
verifyEqual(tc, E.icaweights, W);
verifyTrue(tc, contains(E.history, d{1}.coms{1}));
end

function testIclabelThresholdsShareOneClassification(tc)
% Sibling ICLabel steps reuse one classification; each gives what a fresh
% classification gives.
assumeTrue(tc, exist('pop_iclabel', 'file') == 2, 'ICLabel not installed');
EEG = nqc_synth(struct('seconds', 150, 'nPerCond', 30));
nqc_setBase(EEG);
p = pipecompare.plan.Plan(); p = p.add('highpass', 'cutoff', 0.1); p = p.add('ica', 'fitHighpass', 1);
p = p.add('icremove', 'threshold', {0.5, 0.95}); p = p.add('epoch'); p = p.add('baseline');
r = pipecompare.PipeCompare.optimize(p, nqc_c());
[r.cands.ica] = deal({});
for k = 1:numel(r.leaves)
    verifyEqual(tc, r.cands(k).status, 'ok');
    E = pipecompare.run.Executor.replay(r, k);   % classifies again
    m = pipecompare.eval.Measure.candidate(E, r.contract, r.ref);
    verifyEqual(tc, m.composite, r.cands(k).m.composite, 'RelTol', 1e-10, r.labels{k});
end
verifyGreaterThanOrEqual(tc, r.cands(1).icsRemoved + r.cands(2).icsRemoved, 1);
% a given classification is used as it is: component 2 is Eye
E = EEG; rng(4); E.icaweights = randn(E.nbchan); E.icasphere = eye(E.nbchan); E.icachansind = 1:E.nbchan;
E.icawinv = []; E.icaact = []; E = eeg_checkset(E);
% pop_icflag compares with strict inequalities, so no probability is exactly 0 or 1
cls = repmat([0.94 0.01 0.01 0.01 0.01 0.01 0.01], E.nbchan, 1); cls(2, :) = [0.01 0.01 0.94 0.01 0.01 0.01 0.01];
lab = struct('classification', struct('ICLabel', struct('classes', {{'Brain','Muscle','Eye','Heart', ...
    'Line Noise','Channel Noise','Other'}}, 'classifications', cls, 'version', 'default')), ...
    'com', 'EEG = pop_iclabel(EEG, ''default'');');
[E2, coms, info] = pipecompare.run.Steps.icRemove(E, struct('classes', {{'Eye'}}, 'threshold', 0.8), ...
    struct('iclabel', lab));
verifyEqual(tc, info.comps, 2);
verifyEqual(tc, size(E2.icaweights, 1), E.nbchan - 1);
verifyEqual(tc, coms{1}, lab.com);
end

function testNativeIcaWorkflowRemovesTheFlaggedComponents(tc)
assumeTrue(tc, exist('pop_iclabel', 'file') == 2, 'ICLabel not installed');
EEG = nqc_synth(struct('seconds', 150, 'nPerCond', 30));
nqc_setBase(EEG);
wf = sprintf(['EEG = pop_iclabel(EEG, ''default'');\n', ...
    'EEG = pop_icflag(EEG, [NaN NaN;0.5 1;0.5 1;NaN NaN;NaN NaN;NaN NaN;NaN NaN]);\n', ...
    'EEG = pop_subcomp(EEG, [], 0);']);
p = pipecompare.plan.Plan(); p = p.add('highpass', 'cutoff', 0.1); p = p.add('ica', 'fitHighpass', 1);
p = p.addNative(wf, 'icclean'); p = p.add('epoch'); p = p.add('baseline');
r = pipecompare.PipeCompare.optimize(p, nqc_c());
c = r.cands(1);
verifyEqual(tc, c.status, 'ok');
verifyGreaterThan(tc, c.icsRemoved, 0);                        % blinks are flagged as Eye
verifyEmpty(tc, c.unmatched);
verifyEqual(tc, sum(contains(c.coms, {'pop_iclabel', 'pop_icflag', 'pop_subcomp'})), 3);
% without any ICA before it, the workflow is not a legal step
q = pipecompare.plan.Plan(); q = q.addNative(wf, 'icclean');
verifyError(tc, @() q.enumerate(pipecompare.live.DataState.fromEEG(EEG), nqc_c()), 'PipeCompare:NoLegalPipeline');
end

function testEeglabConfigurationsAreSearchedAsCandidates(tc)
% Several configurations of one step from EEGLAB dialogs, plus "skip",
% are the searched alternatives of that step; each stays intact.
EEG = nqc_synth(struct('seconds', 90, 'nPerCond', 20));
nqc_setBase(EEG);
p = pipecompare.plan.Plan();
p = p.addNative('EEG = pop_eegfiltnew(EEG, ''locutoff'',0.1,''hicutoff'',30,''filtorder'',3300,''plotfreqz'',0);', 'filter');
p = p.addAlternative('filter', pipecompare.plan.Plan.nativeAlt('EEG = pop_eegfiltnew(EEG, ''locutoff'',0.1,''hicutoff'',20,''plotfreqz'',0);'));
p = p.setSkippable('filter', true);
p = p.add('epoch'); p = p.add('baseline');
r = pipecompare.PipeCompare.optimize(p, nqc_c());
verifyEqual(tc, numel(r.leaves), 3);
verifyEqual(tc, sum(contains(r.labels, 'filtorder=3300')), 1);     % the full dialog setting is kept (and shown)
verifyEqual(tc, sum(startsWith(r.labels, 'epoch')), 1);            % the skipped variant
verifyTrue(tc, any(strcmp(r.marginal.parameter, 'filter.(alternative)')));
p = p.setSkippable('filter', false);
verifyEqual(tc, numel(p.Slots(1).alternatives), 2);
end

function testEeglabReferencesAreStrataAndCheckedCorrectly(tc)
% Audit: a reference set in EEGLAB's pop_reref dialog was neither
% stratified nor modelled in the signal check (its channels were ignored).
EEG = nqc_synth(struct('seconds', 90, 'nPerCond', 20, 'artifactTrials', 0));
nqc_setBase(EEG);
p = pipecompare.plan.Plan();
p = p.addEeglab('EEG = pop_reref(EEG, []);', 'ref', 'arg2', {[], {'P7', 'P8'}});
p = p.add('epoch'); p = p.add('baseline');
r = pipecompare.PipeCompare.optimize(p, nqc_c());
verifyEqual(tc, numel(unique({r.leaves.stratum})), 2);          % never ranked against each other
verifyEqual(tc, numel(r.ranking.byStratum), 2);
for k = 1:2
    verifyLessThan(tc, r.cands(k).signal.amplitudeError, 0.01);  % a change of reference is not distortion
    verifyGreaterThan(tc, r.cands(k).signal.waveformCorr, 0.999);
end
end

function testExcludedChannelsAcrossSteps(tc)
% Real-data case: bad-channel detection removed A1/A2, then the average
% reference listed them in 'exclude' and every candidate failed. Channels
% an earlier step removed are skipped; unknown labels are still errors.
% badchannels 'exclude' keeps EOG out of the test, and the channel it
% removes is the right one (pop_rejchan indexes the tested channels).
EEG = nqc_synth(struct('seconds', 90, 'nPerCond', 20));
rng(1); k = find(strcmp({EEG.chanlocs.labels}, 'PO8'));
EEG.data(k, :) = EEG.data(k, :) + 2000 * (rand(1, EEG.pnts) > 0.999);   % one spiky channel
nqc_setBase(EEG);
c = pipecompare.eval.Contract('conditions', {'t', {'11'}; 's', {'31'}}, 'epoch', [-0.2 1], ...
    'baseline', [-0.2 0], 'components', {'P3', [0.3 0.5], {'Pz','P3','P4'}});
p = pipecompare.plan.Plan();
p = p.add('badchannels', 'measure', 'prob', 'threshold', 5, 'action', 'remove', 'exclude', {{'EOG1', 'EOG2'}});
p = p.add('channels', 'labels', {'EOG1'}, 'action', 'remove');
p = p.add('reref', 'mode', 'average', 'exclude', {{'EOG1', 'EOG2'}});
p = p.add('epoch'); p = p.add('baseline');
r = pipecompare.PipeCompare.optimize(p, c);
verifyEqual(tc, r.cands(1).status, 'ok');
E = pipecompare.run.Executor.replay(r, 1);
verifyFalse(tc, any(strcmp({E.chanlocs.labels}, 'PO8')));            % the spiky channel, not a neighbour
verifyTrue(tc, any(strcmp({E.chanlocs.labels}, 'EOG2')));            % EOG was not tested
verifyTrue(tc, any(contains(r.cands(1).coms, 'pop_reref')));
q = pipecompare.plan.Plan(); q = q.add('reref', 'mode', 'average', 'exclude', {{'EOGX'}}); q = q.add('epoch'); q = q.add('baseline');
r = pipecompare.PipeCompare.optimize(q, c);
verifyEqual(tc, r.cands(1).status, 'failed');
verifyTrue(tc, contains(r.cands(1).message, 'not in the dataset: EOGX'));
end

function testBadChannelsDetectedOnAHighPassedCopy(tc)
% detectHighpass: pop_rejchan runs on a high-passed copy (slow drifts
% distort kurtosis and probability); the channels are interpolated or
% removed in the data as they are, whose other channels stay unfiltered.
EEG = nqc_synth(struct('seconds', 90, 'nPerCond', 20));
rng(1); k = find(strcmp({EEG.chanlocs.labels}, 'PO8'));
EEG.data(k, :) = EEG.data(k, :) + 2000 * (rand(1, EEG.pnts) > 0.999);   % one spiky channel
pz = @(E) double(E.data(strcmp({E.chanlocs.labels}, 'Pz'), :));
for action = {'interpolate', 'remove'}
    p = struct('measure', 'prob', 'threshold', 5, 'exclude', {{'EOG1', 'EOG2'}}, ...
        'detectHighpass', 1, 'action', action{1});
    [E, coms, info] = pipecompare.run.Steps.badChannels(EEG, p, struct('highpass', 0));
    verifyTrue(tc, ismember('PO8', info.badChannels));
    verifyFalse(tc, any(ismember({'EOG1', 'EOG2'}, info.badChannels)));
    verifyTrue(tc, contains(coms{1}, 'on a 1 Hz high-passed copy'));
    verifyEqual(tc, pz(E), pz(EEG));                                   % not filtered
    verifyEqual(tc, any(strcmp({E.chanlocs.labels}, 'PO8')), strcmp(action{1}, 'interpolate'));
end
[~, coms] = pipecompare.run.Steps.badChannels(EEG, p, struct('highpass', 1));   % already high-passed: no copy
verifyFalse(tc, any(contains(coms, 'copy')));
end

function testBadChannelMeasuresCombine(tc)
% Standard detects with 'kurt+prob+resid': kurtosis finds a spiky
% channel, resid a noisy one (kurtosis and joint probability flagged it
% only for being less heavy-tailed than channels with blinks, which they
% no longer see), and a channel is bad when any measure flags it.
EEG = nqc_synth(struct('seconds', 90, 'nPerCond', 20, 'noisyChannels', {{'T7'}}, 'noisyUv', 200));
rng(1); k = find(strcmp({EEG.chanlocs.labels}, 'PO8'));
EEG.data(k, :) = EEG.data(k, :) + 2000 * (rand(1, EEG.pnts) > 0.999);   % one spiky channel
p = struct('measure', 'kurt+prob+resid', 'threshold', 5, 'exclude', {{'EOG1', 'EOG2'}}, ...
    'detectHighpass', 1, 'action', 'interpolate');
[~, coms, info] = pipecompare.run.Steps.badChannels(EEG, p, struct('highpass', 0));
verifyTrue(tc, all(ismember({'PO8', 'T7'}, info.badChannels)), strjoin(info.badChannels, ' '));
verifyFalse(tc, ismember('FPz', info.badChannels), strjoin(info.badChannels, ' '));   % blinks do not make it bad
verifyTrue(tc, contains(coms{1}, 'measure kurt+prob+resid'));
p.action = 'remove'; p.detectHighpass = 0;                          % one copy-free path, several measures
[E, ~, info] = pipecompare.run.Steps.badChannels(EEG, p, struct('highpass', 0));
verifyFalse(tc, any(ismember(info.badChannels, {E.chanlocs.labels})));
verifyEqual(tc, E.nbchan, EEG.nbchan - numel(info.badChannels));
end

function testAChannelThatIsJustNoisierIsBad(tc)
% resid: a channel with noise of its own (P8 with several times the
% noise, no spikes) is bad, as the other channels cannot predict it
% (kurtosis and joint probability depend only on the shape of a channel's
% values, not their size); no good channel is flagged. Data whose channels
% are linearly dependent (average-referenced) are not tested, which is
% said.
EEG0 = nqc_synth(struct('seconds', 90, 'nPerCond', 20));
EEG = EEG0; k = find(strcmp({EEG.chanlocs.labels}, 'P8'));
rng(2); EEG.data(k, :) = EEG.data(k, :) + single(25 * randn(1, EEG.pnts));
p = struct('measure', 'resid', 'threshold', 5, 'exclude', {{'EOG1', 'EOG2'}}, ...
    'detectHighpass', 1, 'detectLowpass', 40, 'action', 'interpolate');
[E, coms, info] = pipecompare.run.Steps.badChannels(EEG, p, struct('highpass', 0));
verifyEqual(tc, info.badChannels, {'P8'});
verifyTrue(tc, contains(coms{1}, 'on a 1 Hz high-passed, 40 Hz low-passed copy'), coms{1});
verifyEqual(tc, info.interpolated, {'P8'});
verifyNotEqual(tc, E.data(k, :), EEG.data(k, :));
[~, ~, info] = pipecompare.run.Steps.badChannels(EEG0, p, struct('highpass', 0));
verifyEmpty(tc, info.badChannels);
eog = find(ismember({EEG0.chanlocs.labels}, {'EOG1', 'EOG2'}));
[~, A] = evalc('pop_reref(EEG0, [], ''exclude'', eog)');
[txt, ~, ~, info] = evalc('pipecompare.run.Steps.badChannels(A, p, struct(''highpass'', 0))');
verifyEmpty(tc, info.badChannels);
verifyTrue(tc, contains(txt, 'linearly dependent'), txt);
end

function testChannelBadForPartOfTheRecording(tc)
% stretches: P8 gets noise of its own (30 uV) in the last 20 s of 150 s
% only. As a whole it is not bad, but the two-state HMM of its noise the
% other channels do not explain finds that stretch, and only there is it
% interpolated (the rest of P8 and every other channel as they were); the
% same recording without it has no stretch. The signal-check copy gets the
% same interpolation, and the saved command does it again.
EEG0 = nqc_synth(struct('seconds', 150, 'nPerCond', 20));
EEG = EEG0; lab = {EEG.chanlocs.labels}; k = find(strcmp(lab, 'P8'));
from = round(130 * EEG.srate) + 1;
rng(4); EEG.data(k, from:end) = EEG.data(k, from:end) + single(30 * randn(1, EEG.pnts - from + 1));
p = struct('measure', 'resid', 'threshold', 5, 'exclude', {{'EOG1', 'EOG2'}}, 'detectHighpass', 1, ...
    'detectLowpass', 40, 'action', 'interpolate', 'stretches', true);
[E, coms, info] = pipecompare.run.Steps.badChannels(EEG, p, struct('highpass', 0));
verifyEmpty(tc, info.badChannels);
verifyNumElements(tc, info.stretches, 1);
verifyEqual(tc, info.stretches.channel, 'P8');
verifyEqual(tc, info.stretches.samples(1) / EEG.srate, 130, 'AbsTol', 3);
verifyEqual(tc, info.stretches.samples(2), EEG.pnts);
in = info.stretches.samples(1):info.stretches.samples(2);
verifyEqual(tc, E.data(k, 1:in(1) - 1), EEG.data(k, 1:in(1) - 1));
verifyLessThan(tc, std(diff(double(E.data(k, in)))), 0.5 * std(diff(double(EEG.data(k, in)))));   % its own noise gone
verifyEqual(tc, E.data(setdiff(1:EEG.nbchan, k), :), EEG.data(setdiff(1:EEG.nbchan, k), :));
verifyTrue(tc, contains(coms{end}, 'pipecompare.run.Steps.interpStretches(EEG, struct(''channel'', {''P8''}'), coms{end});
verifyEqual(tc, pipecompare.run.Steps.interpStretches(EEG, info.stretches).data, E.data);
step = struct('type', 'badchannels', 'params', p, 'key', 'badchannels', 'slot', 'badchannels', 'label', 'bad');
S = EEG; S.data = single(randn(size(EEG.data)));
info.badIdx = [];
S2 = pipecompare.run.Steps.replayDecision(step, S, info, struct());
verifyEqual(tc, S2.data, pipecompare.run.Steps.interpStretches(S, info.stretches).data);
verifyNotEqual(tc, S2.data(k, in), S.data(k, in));
[~, ~, info0] = pipecompare.run.Steps.badChannels(EEG0, p, struct('highpass', 0));
verifyEmpty(tc, info0.stretches);
end

function testBadChannelsInterpolatedByTheirEffectOnTheError(tc)
% repair 'error': a bad channel is left as it is only when interpolating
% it would measurably raise the total error. P8, just noisier, is
% interpolated; Pz, the measured electrode, also noisier but with very
% noisy neighbours (not tested), is left as it is: interpolating it would
% bring their noise in and lose the P3 there. The interpolation the
% choice computes with is EEGLAB's own.
EEG = nqc_synth(struct('seconds', 150, 'nPerCond', 40, 'artifactTrials', 0));
lab = {EEG.chanlocs.labels};
near = {'P3', 'P4', 'CP1', 'CP2', 'POz', 'PO3', 'PO4'};
rng(3);
for k = find(ismember(lab, near)), EEG.data(k, :) = EEG.data(k, :) + single(80 * randn(1, EEG.pnts)); end
for k = find(ismember(lab, {'P8', 'Pz'})), EEG.data(k, :) = EEG.data(k, :) + single(25 * randn(1, EEG.pnts)); end
b = find(ismember(lab, {'P8', 'Pz'}));
M = pipecompare.eval.ChannelChoice.interpMatrix(EEG, b);
[~, I] = evalc('pop_interp(EEG, b, ''spherical'')');
verifyEqual(tc, M * double(EEG.data(:, 1:500)), double(I.data(:, 1:500)), 'AbsTol', 1e-3);
nqc_setBase(EEG);
c = pipecompare.eval.Contract('conditions', {'target', {'11'}; 'standard', {'31'}}, ...
    'epoch', [-0.2 1.0], 'baseline', [-0.2 0], 'components', {'P3', [0.30 0.50], {'Pz'}});
ex = [{'EOG1', 'EOG2'} near];
p = pipecompare.plan.Plan();
p = p.add('badchannels', 'measure', 'resid', 'threshold', 5, 'exclude', {ex}, 'detectHighpass', 1, ...
    'detectLowpass', 40, 'repair', 'error');
p = p.add('reref', 'mode', 'average', 'exclude', {{'EOG1', 'EOG2'}});
p = p.add('highpass', 'cutoff', 0.1);
p = p.add('epoch'); p = p.add('baseline');
r = pipecompare.PipeCompare.optimize(p, c);
f = r.cands(1).steps{1};
verifyEqual(tc, sort(cellstr(f.interpolated)), {'P8'});
verifyEqual(tc, cellstr(f.leftAsIs), {'Pz'});
verifyTrue(tc, contains(strjoin(r.cands(1).coms, ' '), 'left as they are'));
verifyTrue(tc, any(contains(pipecompare.simple.Presets.stepsText(r, 1), 'Pz left as they are')));
% adopting it replays the same choice
[~, E] = evalc('pipecompare.run.Executor.replay(r, 1)');
verifyTrue(tc, contains(E.history, sprintf('pop_interp(EEG, %d, ''spherical'')', find(strcmp(lab, 'P8')))), E.history);
verifyFalse(tc, contains(E.history, sprintf('pop_interp(EEG, %s', mat2str(sort(b)))));
end

function testLineNoiseAloneDoesNotMakeAChannelBad(tc)
% detectLowpass: line noise on one channel (P4), which a later low-pass
% removes, makes it bad on a copy that keeps it, not on a copy low-passed
% to the band the analysis keeps.
EEG = nqc_synth(struct('seconds', 90, 'nPerCond', 20));
k = find(strcmp({EEG.chanlocs.labels}, 'P4'));
EEG.data(k, :) = EEG.data(k, :) + single(40 * sin(2 * pi * 50 * (0:EEG.pnts - 1) / EEG.srate));
p = struct('measure', 'kurt+prob+resid', 'threshold', 5, 'exclude', {{'EOG1', 'EOG2'}}, ...
    'detectHighpass', 1, 'action', 'interpolate');
[~, ~, info] = pipecompare.run.Steps.badChannels(EEG, p, struct('highpass', 0));
verifyTrue(tc, ismember('P4', info.badChannels), strjoin(info.badChannels, ' '));
p.detectLowpass = 40;
[~, coms, info] = pipecompare.run.Steps.badChannels(EEG, p, struct('highpass', 0));
verifyFalse(tc, ismember('P4', info.badChannels), strjoin(info.badChannels, ' '));
verifyTrue(tc, contains(coms{1}, '40 Hz low-passed copy'), coms{1});
end

function testBlinksDoNotMakeFrontalChannelsBad(tc)
% Kurtosis and joint probability are measured without the stretches with
% blinks (found in the eye channels): frequent large blinks make the
% frontal channels' values heavy-tailed, yet those channels are good.
EEG = nqc_synth(struct('seconds', 90, 'nPerCond', 20, 'blinkRate', 0.4, 'artifactTrials', 0));
p = struct('measure', 'kurt+prob', 'threshold', 5, 'exclude', {{'EOG1', 'EOG2'}}, ...
    'detectHighpass', 1, 'action', 'interpolate');
[~, coms, info] = pipecompare.run.Steps.badChannels(EEG, p, struct('highpass', 0));
verifyFalse(tc, any(ismember({'FPz', 'F3', 'Fz', 'F4'}, info.badChannels)), strjoin(info.badChannels, ' '));
verifyTrue(tc, contains(coms{1}, 'with blinks left out for kurtosis and probability'), coms{1});
[~, H] = evalc('pop_eegfiltnew(EEG, ''locutoff'', 1, ''plotfreqz'', 0)');   % as the test sees the data
[T, share] = pipecompare.run.Steps.withoutBlinks(H);
verifyGreaterThan(tc, share, 0.05); verifyLessThan(tc, share, 0.5);
verifyEqual(tc, T.pnts, H.pnts - round(share * H.pnts), 'AbsTol', 1);
end

function testFlatChannelsAreBad(tc)
% A channel with no signal (a constant, or almost nothing) is bad whatever
% the measures say, and does not keep the others from being detected; the
% data line names it before the run.
EEG = nqc_synth(struct('seconds', 90, 'nPerCond', 20, 'noisyChannels', {{'T7'}}, 'noisyUv', 200));
f1 = find(strcmp({EEG.chanlocs.labels}, 'Cz')); f2 = find(strcmp({EEG.chanlocs.labels}, 'O1'));
EEG.data(f1, :) = 0;                                                % disconnected: a constant
EEG.data(f2, :) = 1e-3 * randn(1, EEG.pnts);                        % almost nothing
p = struct('measure', 'kurt+prob+resid', 'threshold', 5, 'exclude', {{'EOG1', 'EOG2'}}, ...
    'detectHighpass', 1, 'action', 'interpolate');
[E, coms, info] = pipecompare.run.Steps.badChannels(EEG, p, struct('highpass', 0));
verifyTrue(tc, all(ismember({'Cz', 'O1', 'T7'}, info.badChannels)), strjoin(info.badChannels, ' '));
verifyTrue(tc, contains(coms{1}, 'flat channels [Cz O1]'), coms{1});
verifyGreaterThan(tc, std(double(E.data(f1, :))), 1);              % interpolated from its neighbours
p.action = 'remove'; p.detectHighpass = 0; p.measure = 'kurt';      % the copy-free path
[E, ~, info] = pipecompare.run.Steps.badChannels(EEG, p, struct('highpass', 0));
verifyTrue(tc, all(ismember({'Cz', 'O1'}, info.badChannels)), strjoin(info.badChannels, ' '));
verifyFalse(tc, any(ismember({'Cz', 'O1'}, {E.chanlocs.labels})));
s = pipecompare.live.DataState.fromEEG(EEG);
w = s.warnings(contains(s.warnings, 'Flat channel'));
verifyNumElements(tc, w, 1);
verifyTrue(tc, contains(w{1}, 'Cz, O1') && ~contains(w{1}, 'T7'), w{1});
s = pipecompare.live.DataState.fromEEG(nqc_synth(struct('seconds', 30, 'nPerCond', 5)));
verifyFalse(tc, any(contains(s.warnings, 'Flat channel')));        % none on normal data
end

function testFixedEeglabCommandsAreDecisionMatched(tc)
% A fixed EEGLAB transform re-run on the injected copy is the same
% operation (not flagged); a data-driven one (pop_rejchan) is flagged.
EEG = nqc_synth(struct('seconds', 90, 'nPerCond', 20));
nqc_setBase(EEG);
p = pipecompare.plan.Plan();
p = p.addEeglab('EEG = pop_eegfiltnew(EEG, ''locutoff'',0.1,''hicutoff'',30,''plotfreqz'',0);', 'filter');
p = p.addNative('EEG = pop_reref(EEG, []);', 'ref');
p = p.add('epoch'); p = p.add('baseline'); p = p.add('reject_threshold', 'uv', 150);
r = pipecompare.PipeCompare.optimize(p, nqc_c());
verifyEmpty(tc, r.cands(1).unmatched);
q = pipecompare.plan.Plan();
q = q.addNative('EEG = pop_rejchan(EEG, ''elec'',[1:32],''threshold'',5,''norm'',''on'',''measure'',''kurt'');', 'bad');
q = q.add('epoch'); q = q.add('baseline');
r = pipecompare.PipeCompare.optimize(q, nqc_c());
verifyEqual(tc, r.cands(1).unmatched, {'native'});
end

function testAdoptAndInspectAfterACheckpointedSearch(tc)
% A checkpointed search does not keep the starting dataset in the result;
% adopt / replay / inspect read it back from the checkpoint folder, also
% after changing folder (the path is stored absolute).
EEG = nqc_synth(struct('seconds', 90, 'nPerCond', 20));
nqc_setBase(EEG);
p = pipecompare.plan.Plan(); p = p.add('highpass', 'cutoff', {0.1, 0.5}); p = p.add('epoch'); p = p.add('baseline');
base = tempname; mkdir(base); c0 = onCleanup(@() rmdir(base, 's')); %#ok<NASGU>
here = pwd; cd(base); c1 = onCleanup(@() cd(here)); %#ok<NASGU>
r = pipecompare.PipeCompare.optimize(p, nqc_c(), struct('checkpoint', 'nqc_run1'));   % relative, as in the README
cd(here);
verifyEmpty(tc, r.root);                                          % not kept in memory
verifyTrue(tc, isfile(fullfile(r.options.checkpoint, 'root.mat')));
n0 = evalin('base', 'numel(ALLEEG)');
pipecompare.PipeCompare.adopt(r);
verifyEqual(tc, evalin('base', 'numel(ALLEEG)'), n0 + 1);
E = pipecompare.run.Executor.rootOf(r);
verifyEqual(tc, double(E.data), double(EEG.data), 'AbsTol', 1e-4);
r2 = r; r2.identity = 'another search';
verifyError(tc, @() pipecompare.run.Executor.rootOf(r2), 'PipeCompare:Checkpoint');
end

function testRemovedEpochsAreReadNotGuessed(tc)
% pop_rejepoch's epochs are read from its argument before it runs, so an
% epoch that has no events is not mistaken for a removed one; a removal
% that cannot be told apart is re-run (flagged) instead of guessed.
EEG = nqc_synth(struct('seconds', 60, 'nPerCond', 10));
[~, E] = evalc('pop_epoch(EEG, {''11'', ''31''}, [-0.2 0.8])');
E.event([E.event.epoch] == 5) = []; E = eeg_checkset(E, 'eventconsistency');   % epoch 5 has no event
[E2, ~, info] = pipecompare.run.Steps.native(E, 'EEG = pop_rejepoch(EEG, [2 3], 0);');
verifyEqual(tc, info.decisions.idx, [2 3]);
verifyEqual(tc, E2.trials, E.trials - 2);
[E3, ~, info] = pipecompare.run.Steps.native(E, 'EEG = pop_eegthresh(EEG,1,[1:32],-150,150,-0.2,0.796,0,1);');
verifyLessThan(tc, E3.trials, E.trials);                          % it did remove epochs
verifyEmpty(tc, info.decisions);                                  % not decision-matched: re-run on the copy
end

function testStrataDoNotDependOnStepOrder(tc)
% The same measure-defining choices in another order are one stratum.
p = pipecompare.plan.Plan();
p = p.add('reref', 'mode', 'average');
p = p.addEeglab('EEG = pop_reref(EEG, {''P7'',''P8''});', 'ref2');
p.OrderMode = 'search';
leaves = p.enumerate(nqc_fakeState(false, 500), nqc_contract());
verifyEqual(tc, numel(leaves), 2);
verifyEqual(tc, numel(unique({leaves.stratum})), 1);
q = pipecompare.plan.Plan(); q = q.add('reref', 'mode', 'channels');          % no reference channels
try
    q.enumerate(nqc_fakeState(false, 500), nqc_contract());
    verifyFail(tc, 'reref channels without channels must be refused at plan time');
catch ME
    verifyEqual(tc, ME.identifier, 'PipeCompare:NoLegalPipeline');
    verifyTrue(tc, contains(ME.message, 'needs the reference channels'));
end
end

function testSeveralStrataAskTheUserToChoose(tc)
EEG = nqc_synth(struct('seconds', 90, 'nPerCond', 20, 'artifactTrials', 0));
nqc_setBase(EEG);
p = pipecompare.plan.Plan();
p = p.addChoice('ref', {'reref', 'mode', 'average'}, {'reref', 'mode', 'channels', 'channels', {'P7', 'P8'}});
p = p.add('epoch'); p = p.add('baseline');
r = pipecompare.PipeCompare.optimize(p, nqc_c());
verifyEqual(tc, numel(r.ranking.byStratum), 2);
verifyEmpty(tc, r.ranking.recommended);
try
    pipecompare.PipeCompare.adopt(r);
    verifyFail(tc, 'adopt without a choice must explain the strata');
catch ME
    verifyEqual(tc, ME.identifier, 'PipeCompare:Adopt');
    verifyTrue(tc, contains(ME.message, 'strata'));
end
pipecompare.PipeCompare.adopt(r, r.ranking.byStratum(1).recommended);   % a chosen stratum's recommendation
end

function testRoiOnAChannelRestoredByThePlan(tc)
% Pz was removed before PipeCompare; the plan restores it and measures there.
% The ROI check must not refuse this, the restored Pz carries the known
% signal, and a pipeline without the restore fails with that reason.
EEG = nqc_synth(struct('seconds', 90, 'nPerCond', 20, 'artifactTrials', 0));
E = pop_select(EEG, 'nochannel', {'Pz'});
nqc_setBase(E);
c = pipecompare.eval.Contract('conditions', {'t', {'11'}; 's', {'31'}}, 'epoch', [-0.2 1], ...
    'baseline', [-0.2 0], 'components', {'P3', [0.3 0.5], {'Pz'}});
p = pipecompare.plan.Plan(); p = p.addChoice('rest', {'restore'}, 'none'); p = p.add('epoch'); p = p.add('baseline');
r = pipecompare.PipeCompare.optimize(p, c);
k = find(contains(r.labels, 'restore'));
verifyEqual(tc, r.cands(k).status, 'ok');
% the restored Pz is checked against the true field there: interpolating
% the peak channel of a focal component loses ~28% of it, and that is what
% the check reports (test_signal/testInterpolatingAnRoiChannelIsMeasured)
verifyGreaterThan(tc, r.cands(k).signal.amplitudeError, 0.15);
verifyLessThan(tc, r.cands(k).signal.amplitudeError, 0.45);
verifyGreaterThan(tc, r.cands(k).signal.waveformCorr, 0.95);
j = setdiff(1:2, k);
verifyEqual(tc, r.cands(j).status, 'failed');
verifyTrue(tc, contains(r.cands(j).message, 'ROI channel(s) missing'));
c2 = pipecompare.eval.Contract('conditions', {'t', {'11'}}, 'epoch', [-0.2 1], 'baseline', [-0.2 0], ...
    'components', {'X', [0.3 0.5], {'NoSuchChannel'}});
verifyError(tc, @() pipecompare.PipeCompare.optimize(p, c2), 'PipeCompare:Contract');   % unknown channels still refused
end

% ---------------------------------------------------------------- helpers
function testEpochsRepairedByInterpolationWithinEpochs(tc)
% An epoch with at most n channels over the limit keeps them, interpolated
% within that epoch exactly as EEGLAB's eeg_interp does; an epoch with more
% is rejected; the signal copy gets the same channels in the same epochs.
EEG = nqc_synth(struct('seconds', 150, 'nPerCond', 30, 'artifactTrials', 0));
[~, E] = evalc('nqc_epoched(EEG)');
ch = find(strcmpi({E.chanlocs.labels}, 'T7')); ch2 = find(strcmpi({E.chanlocs.labels}, 'P8'));
E.data(ch, :, [3 7 12]) = E.data(ch, :, [3 7 12]) + 500;      % one channel far over the limit
E.data(ch2, :, 7) = E.data(ch2, :, 7) - 500;                   % two in epoch 7
E.data(:, :, 20) = E.data(:, :, 20) + 500;                     % every channel: rejected
pz = find(strcmpi({E.chanlocs.labels}, 'Pz'));
E.data(pz, :, 25) = E.data(pz, :, 25) + 500;                   % a measured electrode: rejected, not repaired
in = struct('type', 'reject_threshold', 'params', struct('uv', 400, 'exclude', {{}}, 'interpolate', 3), ...
    'key', 'rej', 'slot', 'rej', 'label', 'rej');
x = struct('contract', nqc_c(), 'highpass', 0);
[E2, coms, info] = pipecompare.run.Steps.run(in, E, x);
verifyEqual(tc, info.epochsInterpolated, 3);
verifyEqual(tc, info.epochsNotRepaired, 1);
verifyEqual(tc, info.rejIdx, [20 25]);
verifyEqual(tc, E2.trials, E.trials - 2);
verifyTrue(tc, any(contains(coms, 'pipecompare.run.Steps.interpolateEpochs')));
[~, R] = evalc('pop_select(E, ''trial'', [3 12])');
[~, R] = evalc('eeg_interp(R, ch, ''spherical'')');
verifyEqual(tc, double(E2.data(ch, :, [3 12])), double(R.data(ch, :, :)), 'AbsTol', 1e-3);
[~, R] = evalc('pop_select(E, ''trial'', 7)');
[~, R] = evalc('eeg_interp(R, [ch ch2], ''spherical'')');
verifyEqual(tc, double(E2.data([ch ch2], :, 7)), double(R.data([ch ch2], :, :)), 'AbsTol', 1e-3);
verifyEqual(tc, E2.data(:, :, 1), E.data(:, :, 1));            % other epochs untouched
S = pipecompare.run.Steps.replayDecision(in, E, info, x);       % the same decisions on a copy
verifyEqual(tc, double(S.data), double(E2.data), 'AbsTol', 1e-4);
in.params.interpolate = 0;                                      % off: rejected
[~, ~, info0] = pipecompare.run.Steps.run(in, E, x);
verifyEqual(tc, sort(info0.rejIdx), [3 7 12 20 25]);
in.params.interpolate = 1;                                      % epoch 7 has two channels over: rejected
[~, ~, info1] = pipecompare.run.Steps.run(in, E, x);
verifyEqual(tc, sort(info1.rejIdx), [7 20 25]);
end

function testRejectionLimitFromTheDataIsTheExactOptimum(tc)
% uv 'auto': every limit that keeps a different set of epochs is scored;
% the best one equals a brute-force search, the one chosen is the highest
% the data show equivalent to it, and pop_eegthresh rejects exactly
% the epochs over it.
EEG = nqc_synth(struct('seconds', 150, 'nPerCond', 30));
[~, E] = evalc('nqc_epoched(EEG)');
c = nqc_c();
ref = pipecompare.eval.Measure.reference(E, c);
ref.blink = {};   % the SME part alone (what blinks add: testBlinksCountInTheTotalError)
in = struct('type', 'reject_threshold', 'params', struct('uv', 'auto', 'exclude', {{}}), 'key', 'rej', 'slot', 'rej', 'label', 'rej');
x = struct('contract', c, 'highpass', 0, 'ref', ref, 'rank', struct());
[E2, coms, info] = pipecompare.run.Steps.run(in, E, x);
u = info.uvChosen; th = info.threshold;
m = reshape(max(max(abs(double(E.data)), [], 1), [], 2), 1, []);
verifyEqual(tc, sort(info.rejIdx(:)'), find(m > u));              % exactly the epochs over the limit
verifyTrue(tc, contains(coms{1}, 'pop_eegthresh'));
verifyGreaterThanOrEqual(tc, u, th.uvBest);
verifyLessThanOrEqual(tc, th.objectiveChosen, 1.05 * th.objectiveBest);   % equivalent to the best, not only not distinguished
verifyGreaterThan(tc, info.rejected, 0);                          % the movement artifacts go
% brute force over every limit (one per value of m)
T = pipecompare.eval.Measure.trials(E, c); mt = m(T.epoch);
need = max(10, ceil(0.5 * ref.n));
best = Inf; uBest = NaN;
for v = unique(mt(:)')
    s2 = zeros(1, 2);
    for k = 1:2
        y = T.data{1}(T.cond == k & mt(:) <= v);
        if numel(y) < need(k), s2(k) = Inf; else, s2(k) = var(y) / numel(y); end
    end
    o = sqrt(mean(s2));
    if o <= best, best = o; uBest = v; end
end
verifyEqual(tc, th.objectiveBest, best, 'RelTol', 1e-9);
verifyEqual(tc, find(m > th.uvBest), find(m > uBest));            % the same epochs kept
% the objective at the chosen limit is the ranking's own (SME of the data it leaves)
mc = pipecompare.eval.Measure.candidate(E2, c, ref);
verifyEqual(tc, mc.composite, th.objectiveChosen, 'RelTol', 1e-9);
% peak-to-peak: the same, on the moving-window peak-to-peak values
in.params.method = 'peaktopeak'; in.params.window = 200;
[~, ~, info] = pipecompare.run.Steps.run(in, E, x);
pp = max(pipecompare.run.Steps.peakToPeak(E, 1:E.nbchan, 200), [], 1);
verifyEqual(tc, sort(info.rejIdx(:)'), find(pp > info.uvChosen));
% peak measures have no score per trial: a clear error
x.contract = pipecompare.eval.Contract('conditions', {'t', {'11'}; 's', {'31'}}, 'epoch', [-0.2 1], ...
    'baseline', [-0.2 0], 'components', {'P3', [0.3 0.5], {'Pz'}, {'peakAmplitude', 'positive'}});
verifyError(tc, @() pipecompare.run.Steps.run(in, E, x), 'PipeCompare:Threshold');
end

function testRejectionLimitFromTheDataEndToEnd(tc)
% In a search: with movement artifacts, the limit chosen from the data
% beats rejecting almost nothing; adopt rebuilds the same pipeline (same
% limit, same SME) and the saved script keeps 'auto'.
[EEG, truth] = nqc_synth();
nqc_setBase(EEG);
c = pipecompare.eval.Contract('conditions', {'target', {'11'}; 'standard', {'31'}}, ...
    'epoch', [-0.2 1.0], 'baseline', [-0.2 0], 'components', {'P3', [0.30 0.50], truth.roi});
p = pipecompare.plan.Plan();
p = p.add('highpass', 'cutoff', 0.1); p = p.add('lowpass', 'cutoff', 30);
p = p.add('epoch'); p = p.add('baseline');
p = p.add('reject_threshold', 'uv', {'auto', 1000});
r = pipecompare.PipeCompare.optimize(p, c);
rec = r.ranking.recommended;
verifyTrue(tc, contains(r.labels{rec}, 'auto'), r.labels{rec});
st = r.cands(rec).steps{end};
verifyTrue(tc, isfinite(st.uvChosen) && st.uvChosen < 1000);
verifyEqual(tc, r.ranking.byStratum.crossfit.nested, 1);         % its limit chosen again in the selection check
L = pipecompare.simple.Presets.stepsText(r, rec);
verifyTrue(tc, any(contains(L, sprintf('+/-%g uV', st.uvChosen)) & contains(L, 'limit chosen from these data')), strjoin(L, newline));
pipecompare.PipeCompare.adopt(r);
adopted = evalin('base', 'EEG');
verifyTrue(tc, contains(adopted.history, 'pop_eegthresh'));
m = pipecompare.eval.Measure.candidate(adopted, c, r.ref);
verifyEqual(tc, m.composite, r.cands(rec).m.composite, 'RelTol', 1e-9);
f = fullfile(tempdir, 'pc_auto_limit_script.m'); cl = onCleanup(@() delete(f)); %#ok<NASGU>
pipecompare.PipeCompare.writeScript(r, rec, f);
verifyTrue(tc, contains(fileread(f), 'auto'), fileread(f));
end

function testIclabelThresholdFromTheDataIsTheExactOptimum(tc)
% threshold 'auto': every threshold that removes a different set of
% components is scored with the rejection that follows at its best limit;
% each set's objective equals what EEGLAB's own processing gives
% (pop_subcomp, epoching, baseline, a brute-force search over the limits,
% the gain from the signal copy), the one chosen removes no more
% components than the best, and pop_icflag removes exactly that set.
EEG = nqc_synth(struct('seconds', 150, 'nPerCond', 30));
c = nqc_c();
n = EEG.nbchan; rng(4);
EEG.icaweights = eye(n) + 0.3 * randn(n); EEG.icasphere = eye(n); EEG.icachansind = 1:n;
EEG.icawinv = []; EEG.icaact = []; EEG = eeg_checkset(EEG);
q = 0.2 * ones(1, n);
q(1:6) = [0.97 0.9 0.8 0.65 0.55 0.45];      % 5 components over 0.5: 6 sets with removing nothing
q(strcmpi({EEG.chanlocs.labels}, 'Pz')) = 0.6;
q = single(q);
ref = pipecompare.eval.Measure.reference(EEG, c);
ref.blink = {};   % the SME part alone (what blinks add: testBlinksCountInTheTotalError)
[S, truth] = pipecompare.eval.Injection.prepare(EEG, c, ref);
ins = @(type, params) struct('type', type, 'params', params, 'key', type, 'slot', type, 'label', type);
after = {ins('epoch', struct()), ins('baseline', struct()), ins('reject_threshold', struct('uv', 'auto', 'exclude', {{}}))};
x = struct('contract', c, 'highpass', 0, 'ref', ref, 'rank', struct(), 'after', {after}, ...
    'signal', struct('S', S, 'truth', truth));
[t, info] = pipecompare.eval.IcThreshold.choose(EEG, q, x);
vals = sort(unique(double(q(q > 0.5))), 'descend');
sets = [{[]} arrayfun(@(v) find(q >= v), vals, 'UniformOutput', false)];
verifyEqual(tc, info.sets, numel(sets));
need = max(10, ceil(0.5 * ref.n));
for k = 1:numel(sets)
    E = EEG; Sk = S;
    if ~isempty(sets{k}), E = pop_subcomp(E, sets{k}, 0); Sk = pop_subcomp(Sk, sets{k}, 0); end
    [~, E] = evalc('nqc_epoched(E)'); [~, Sk] = evalc('nqc_epoched(Sk)');
    r = pipecompare.eval.Injection.compare(Sk, c, truth, {});
    m = reshape(max(max(abs(double(E.data)), [], 1), [], 2), 1, []);
    T = pipecompare.eval.Measure.trials(E, c); mt = m(T.epoch);
    best = Inf;
    for v = unique(mt(:)')
        s2 = zeros(1, 2);
        for cc = 1:2
            y = T.data{1}(T.cond == cc & mt(:) <= v) / r.gain;
            if numel(y) < need(cc), s2(cc) = Inf; else, s2(cc) = var(y) / numel(y); end
        end
        best = min(best, sqrt(mean(s2)));
    end
    if isempty(pipecompare.eval.Rank.signalReasons(r, struct())) && isfinite(best)
        verifyEqual(tc, info.objectives(k), best, 'RelTol', 1e-4, sprintf('set %d', k));
    else
        verifyTrue(tc, isnan(info.objectives(k)), sprintf('set %d: %s', k, info.why{k}));
    end
end
verifyEqual(tc, info.objectiveBest, min(info.objectives));
verifyLessThanOrEqual(tc, info.removedChosen, info.removedBest);
verifyLessThanOrEqual(tc, info.objectiveChosen, 1.05 * info.objectiveBest);   % equivalent to the best
verifyEqual(tc, sum(q > t), info.removedChosen);                 % pop_icflag flags q > t
verifyTrue(tc, any(isnan(info.objectives)), mat2str(info.objectives));   % removing Pz loses the signal
% in the step: the same threshold, and pop_icflag removes exactly that set
if exist('pop_iclabel', 'file') == 2
    cls = zeros(n, 7, 'single'); cls(:, 3) = q; cls(:, 1) = 1 - q;
    x.iclabel = struct('classification', struct('ICLabel', struct('classes', {{'Brain','Muscle','Eye','Heart', ...
        'Line Noise','Channel Noise','Other'}}, 'classifications', cls, 'version', 'default')), ...
        'com', 'EEG = pop_iclabel(EEG, ''default'');');
    [~, ~, si] = pipecompare.run.Steps.icRemove(EEG, struct('classes', {{'Eye'}}, 'threshold', 'auto'), x);
    verifyEqual(tc, si.thresholdChosen, t);
    verifyEqual(tc, si.comps, find(q > t));
end
% peak measures have no score per trial: a clear error
x.contract = pipecompare.eval.Contract('conditions', {'t', {'11'}; 's', {'31'}}, 'epoch', [-0.2 1], ...
    'baseline', [-0.2 0], 'components', {'P3', [0.3 0.5], {'Pz'}, {'peakAmplitude', 'positive'}});
verifyError(tc, @() pipecompare.eval.IcThreshold.choose(EEG, q, x), 'PipeCompare:IcThreshold');
end

function testIcComponentsChosenBySme(tc)
% threshold 'sme': components are added one at a time, in the order of
% their own variance at the score over their share of the known signal
% (the artifact components first); each set on the path is scored
% exactly as EEGLAB's own processing gives (pop_subcomp, epoching,
% baseline, the best limit, the gain from the signal copy); a component
% that may not be removed (ICLabel takes it for brain activity) never is;
% the one chosen removes no more components than the best and is
% equivalent to it; in the step, exactly those components are removed.
% Three components carry a strong artifact into the measured electrodes
% (and almost none of the known signal): they are the first ones removed.
EEG = nqc_synth(struct('seconds', 150, 'nPerCond', 30, 'blinkRate', 0, 'artifactTrials', 0));
c = nqc_c();
n = EEG.nbchan; rng(4);
roi = find(ismember(lower({EEG.chanlocs.labels}), lower(c.components(1).roi)));
xyz = [[EEG.chanlocs.X]' [EEG.chanlocs.Y]' [EEG.chanlocs.Z]']; xyz = xyz ./ vecnorm(xyz, 2, 2);
[~, far] = sort(xyz * mean(xyz(roi, :), 1)');                  % farthest from the measured electrodes first
art = sort(far(1:3))'; evk = far(4);
A = eye(n); A(roi, [art evk]) = 1;                               % component k: its channel and the ROI
src = filter(1, [1 -0.9], randn(numel(art), EEG.pnts), [], 2);
EEG.data = EEG.data + single(A(:, art) * (30 * src / std(src(:))));
% one more carries a time-locked deflection into the measured window: brain
% activity ICLabel did not recognise; it may not be removed
lat = [EEG.event(ismember({EEG.event.type}, {'11', '31'})).latency];
p3 = zeros(1, EEG.pnts); tt = (0:EEG.pnts-1);
for L = lat, p3 = p3 + 6 * exp(-0.5 * ((tt - L - 0.4 * EEG.srate) / (0.07 * EEG.srate)) .^ 2); end
EEG.data = EEG.data + single(A(:, evk) * p3);
EEG.icawinv = A; EEG.icaweights = inv(A); EEG.icasphere = eye(n); EEG.icachansind = 1:n;
EEG.icaact = []; EEG = eeg_checkset(EEG);
ref = pipecompare.eval.Measure.reference(EEG, c);
ref.blink = {};   % the SME part alone (what blinks add: testBlinksCountInTheTotalError)
[S, truth] = pipecompare.eval.Injection.prepare(EEG, c, ref);
ins = @(type, params) struct('type', type, 'params', params, 'key', type, 'slot', type, 'label', type);
after = {ins('epoch', struct()), ins('baseline', struct()), ins('reject_threshold', struct('uv', 'auto', 'exclude', {{}}))};
x = struct('contract', c, 'highpass', 0, 'ref', ref, 'rank', struct(), 'after', {after}, ...
    'signal', struct('S', S, 'truth', truth));
brain = roi(:)';                                                % may not be removed
cand = setdiff(1:n, brain);
[comps, info] = pipecompare.eval.IcThreshold.chooseBySme(EEG, cand, x);
verifyFalse(tc, any(ismember(info.order, brain)), mat2str(info.order));
verifyEqual(tc, info.protected, evk);                            % time-locked: kept
verifyFalse(tc, ismember(evk, info.order));
verifyEqual(tc, numel(unique(info.order)), numel(info.order));
verifyEqual(tc, info.sets, numel(info.order) + 1);
verifyGreaterThan(tc, info.sets, 3);
verifyEqual(tc, sort(info.order(1:3)), art);                     % the artifact components first
need = max(10, ceil(0.5 * ref.n));
for k = 1:info.sets
    R = info.order(1:k-1);
    E = EEG; Sk = S;
    if ~isempty(R), E = pop_subcomp(E, R, 0); Sk = pop_subcomp(Sk, R, 0); end
    [~, E] = evalc('nqc_epoched(E)'); [~, Sk] = evalc('nqc_epoched(Sk)');
    r = pipecompare.eval.Injection.compare(Sk, c, truth, {});
    m = reshape(max(max(abs(double(E.data)), [], 1), [], 2), 1, []);
    T = pipecompare.eval.Measure.trials(E, c); mt = m(T.epoch);
    best = Inf;
    for v = unique(mt(:)')
        s2 = zeros(1, 2);
        for cc = 1:2
            y = T.data{1}(T.cond == cc & mt(:) <= v) / r.gain;
            if numel(y) < need(cc), s2(cc) = Inf; else, s2(cc) = var(y) / numel(y); end
        end
        best = min(best, sqrt(mean(s2)));
    end
    if isempty(pipecompare.eval.Rank.signalReasons(r, struct())) && isfinite(best)
        verifyEqual(tc, info.objectives(k), best, 'RelTol', 1e-4, sprintf('set %d', k));
    else
        verifyTrue(tc, isnan(info.objectives(k)), sprintf('set %d: %s', k, info.why{k}));
    end
end
f = info.objectives;
verifyLessThan(tc, min(f), f(1));                                % removing components helps here
verifyEqual(tc, info.objectiveBest, min(f));
verifyEqual(tc, comps, sort(info.order(1:info.removedChosen)));
verifyLessThanOrEqual(tc, info.removedChosen, info.removedBest);
verifyLessThanOrEqual(tc, info.objectiveChosen, 1.05 * info.objectiveBest);   % equivalent to the best
% in the step: exactly those components removed, none ICLabel takes for brain activity
if exist('pop_iclabel', 'file') == 2
    % the others Other, as on recordings ICLabel recognises little (an Eye
    % component would not be kept for being time-locked)
    cls = repmat(single([0.1 0 0 0 0 0 0.9]), n, 1); cls(brain, :) = repmat(single([0.9 0 0 0 0 0 0.1]), numel(brain), 1);
    x.iclabel = struct('classification', struct('ICLabel', struct('classes', {{'Brain','Muscle','Eye','Heart', ...
        'Line Noise','Channel Noise','Other'}}, 'classifications', cls, 'version', 'default')), ...
        'com', 'EEG = pop_iclabel(EEG, ''default'');');
    [E2, coms, si] = pipecompare.run.Steps.icRemove(EEG, struct('classes', {{'Eye'}}, 'threshold', 'sme'), x);
    verifyEqual(tc, si.comps, comps);
    verifyEqual(tc, si.icsConsidered, numel(cand));
    verifyTrue(tc, any(contains(coms, 'gcompreject')), strjoin(coms, newline));
    E3 = EEG; if ~isempty(comps), E3 = pop_subcomp(EEG, comps, 0); end
    verifyEqual(tc, double(E2.data), double(E3.data), 'AbsTol', 1e-3);
end
% peak measures have no score per trial: a clear error
x.contract = pipecompare.eval.Contract('conditions', {'t', {'11'}; 's', {'31'}}, 'epoch', [-0.2 1], ...
    'baseline', [-0.2 0], 'components', {'P3', [0.3 0.5], {'Pz'}, {'peakAmplitude', 'positive'}});
verifyError(tc, @() pipecompare.eval.IcThreshold.chooseBySme(EEG, cand, x), 'PipeCompare:IcThreshold');
end

function testBlinkComponentICLabelCallsOtherIsNotKept(tc)
% threshold 'sme': a component that carries blinks time-locked to the
% stimulus into the measured window is not kept as time-locked brain
% activity when its activity follows an eye channel, though ICLabel
% would call it Other; a component with time-locked activity of its own
% (F3's here) still is.
EEG = nqc_synth(struct('seconds', 150, 'nPerCond', 30, 'blinkRate', 0, 'artifactTrials', 0));
c = nqc_c();
n = EEG.nbchan; lab = {EEG.chanlocs.labels};
roi = find(ismember(lower(lab), lower(c.components(1).roi)));
eog = find(strcmp(lab, 'EOG1')); evk = find(strcmp(lab, 'F3'));
A = eye(n); A(roi, eog) = 0.5; A(roi, evk) = 1;
lat = [EEG.event(ismember({EEG.event.type}, {'11', '31'})).latency];
rng(5); tt = (0:EEG.pnts-1);
blink = zeros(1, EEG.pnts); p3 = zeros(1, EEG.pnts);
for L = lat
    if rand < 0.5, blink = blink + 100 * exp(-0.5 * ((tt - L - 0.4 * EEG.srate) / (0.05 * EEG.srate)) .^ 2); end
    p3 = p3 + 6 * exp(-0.5 * ((tt - L - 0.4 * EEG.srate) / (0.07 * EEG.srate)) .^ 2);
end
EEG.data = EEG.data + single(A(:, eog) * blink + A(:, evk) * p3);
EEG.icawinv = A; EEG.icaweights = inv(A); EEG.icasphere = eye(n); EEG.icachansind = 1:n;
EEG.icaact = []; EEG = eeg_checkset(EEG);
ref = pipecompare.eval.Measure.reference(EEG, c); ref.blink = {};
[S, truth] = pipecompare.eval.Injection.prepare(EEG, c, ref);
ins = @(type, params) struct('type', type, 'params', params, 'key', type, 'slot', type, 'label', type);
x = struct('contract', c, 'highpass', 0, 'ref', ref, 'rank', struct(), ...
    'after', {{ins('epoch', struct()), ins('baseline', struct())}}, 'signal', struct('S', S, 'truth', truth));
[~, info] = pipecompare.eval.IcThreshold.chooseBySme(EEG, setdiff(1:n, roi), x);
k = find(setdiff(1:n, roi) == eog);
verifyTrue(tc, info.evoked.eye(k));
verifyFalse(tc, ismember(eog, info.protected), mat2str(info.protected));
verifyTrue(tc, ismember(evk, info.protected), mat2str(info.protected));
verifyTrue(tc, ismember(eog, info.order(1:min(3, end))), mat2str(info.order));   % the blinks: removed early
end

function testFineTunedCutoffsAreExactWithTheDecisionsFixed(tc)
% With a pipeline's decisions held fixed (components and epochs removed),
% the objective of other cutoffs computed on the measured electrodes alone
% equals filtering every channel with EEGLAB and removing the same
% components and epochs; at the pipeline's own cutoffs it is the
% pipeline's objective. In the search, better cutoffs run as one more
% pipeline, marked as fine-tuned.
assumeTrue(tc, exist('pop_iclabel', 'file') == 2, 'ICLabel not installed');
EEG = nqc_synth(struct('seconds', 150, 'nPerCond', 30));
nqc_setBase(EEG);
c = nqc_c();
p = pipecompare.plan.Plan();
p = p.add('reref', 'mode', 'average'); p = p.add('ica', 'fitHighpass', 1);
p = p.add('highpass', 'cutoff', {0.1, 0.5}); p = p.add('lowpass', 'cutoff', {20, 40});
p = p.add('icremove', 'threshold', 0.8);
p = p.add('epoch'); p = p.add('baseline'); p = p.add('reject_threshold', 'uv', 100);
r = pipecompare.PipeCompare.optimize(p, c, struct('refineCutoffs', true));
verifyEqual(tc, r.report.nLeaves, 4);
k = r.ranking.recommended;
verifyNotEmpty(tc, k);
cand = r.cands(k); path = r.leaves(k).path;
types = cellfun(@(i) i.type, path, 'UniformOutput', false);
dec = struct('comps', cand.steps{strcmp(types, 'icremove')}.comps, 'rejIdx', cand.steps{strcmp(types, 'reject_threshold')}.rejIdx);
% the data before the first filter: the reference, then the search's ICA
ctx = struct('contract', c);
[~, D] = evalc('pipecompare.run.Steps.run(path{1}, r.root, ctx)'); D = withIca(D, cand.ica{1});
[S, truth] = pipecompare.eval.Injection.prepare(r.root, c, r.ref);
f = find(ismember(types, {'highpass', 'lowpass'}));
env = struct('contract', c, 'ref', r.ref, 'rank', r.options, 'S', S, 'truth', truth, 'grid', {{[0.1 0.3 0.5], [20 30 40]}});
[~, info] = pipecompare.eval.Cutoff.search(D, path(f(1):end), dec, struct('m', cand.m, 'gain', cand.signal.gain), env);
own = cellfun(@(i) i.params.cutoff, path(f));   % the recommended may itself be a fine-tuned pipeline
verifyEqual(tc, info.pairs, numel(unique([0.1 0.3 0.5 own(1)])) * numel(unique([20 30 40 own(2)])));
verifyEqual(tc, info.reproduced, info.objectiveBase, 'RelTol', 1e-4);
% 0.3 Hz and 30 Hz: every channel filtered by EEGLAB, the same decisions,
% the gain from the known-signal copy processed the same way
pathX = path;
pathX{f(1)} = pipecompare.plan.Plan.reparam(path{f(1)}, 'cutoff', 0.3);
pathX{f(2)} = pipecompare.plan.Plan.reparam(path{f(2)}, 'cutoff', 30);
E = D;
[~, Sx] = evalc('pipecompare.run.Steps.replayDecision(path{1}, S, struct(), ctx)'); Sx = withIca(Sx, cand.ica{1});
for q = [f(:)' find(strcmp(types, 'icremove')) find(strcmp(types, 'epoch')) find(strcmp(types, 'baseline'))]
    in = pathX{q};
    if strcmp(in.type, 'icremove')
        if ~isempty(dec.comps), E = pop_subcomp(E, dec.comps, 0); Sx = pop_subcomp(Sx, dec.comps, 0); end
        continue;
    end
    [~, E] = evalc('pipecompare.run.Steps.run(in, E, ctx)');
    [~, Sx] = evalc('pipecompare.run.Steps.run(in, Sx, ctx)');
end
if ~isempty(dec.rejIdx), E = pop_rejepoch(E, dec.rejIdx, 0); Sx = pop_rejepoch(Sx, dec.rejIdx, 0); end
m = pipecompare.eval.Measure.candidate(E, c, r.ref);
g = pipecompare.eval.Injection.compare(Sx, c, truth, pathX).gain;
j = find(all(abs(info.grid - [0.3 30]) < 1e-9, 2));
verifyEqual(tc, info.objectives(j), pipecompare.eval.Rank.pointObjective(m, r.ref, r.options, g), 'RelTol', 2e-3);
% in the search: one more pipeline when other cutoffs are better
verifyTrue(tc, any(numel(r.cands) == [4 5]));
if numel(r.cands) == 5
    verifyTrue(tc, any(cellfun(@(x) x.cutoffRefined, r.cands(5).steps)));
    verifyTrue(tc, any(contains(pipecompare.simple.Presets.stepsText(r, 5), 'fine-tuned from these data')));
    verifyTrue(tc, strcmp(r.cands(5).status, 'ok'));
end
% the grid: round steps between the compared cutoffs, those included
verifyEqual(tc, pipecompare.eval.Cutoff.grid([0.1 0.3 0.5 1]), 0.1:0.05:1, 'AbsTol', 1e-12);
verifyEqual(tc, pipecompare.eval.Cutoff.grid([20 30 40]), 20:2:40, 'AbsTol', 1e-12);
end

function testWeightedEpochsCountByTheirNoise(tc)
% Weighting the epochs by their noise keeps every epoch and, on data with
% large artifacts in some trials, measures more precisely than equal
% weights. The weights are stored with the epochs, the scores are their
% weighted mean (checked by hand), the signal check gives nearly the same
% gain (the known signal is the same in every epoch), and the command recorded
% for the step gives the same weights again.
EEG = nqc_synth(struct('seconds', 150, 'nPerCond', 30));
nqc_setBase(EEG);
c = nqc_c();
p = pipecompare.plan.Plan();
p = p.add('highpass', 'cutoff', 0.5); p = p.add('epoch'); p = p.add('baseline');
p = p.addChoice('weight', {'weight_trials'}, 'none');
r = pipecompare.PipeCompare.optimize(p, c);
verifyEqual(tc, {r.cands.status}, {'ok', 'ok'});
weighted = arrayfun(@(l) any(cellfun(@(x) strcmp(x.type, 'weight_trials'), l.path)), r.leaves);
kw = find(weighted); k0 = find(~weighted);
verifyEqual(tc, r.cands(kw).m.retention, [1 1]);
verifyLessThan(tc, r.ranking.table.objective(kw), r.ranking.table.objective(k0));
verifyEqual(tc, r.cands(kw).signal.gain, r.cands(k0).signal.gain, 'RelTol', 0.02);   % (filter tails of neighbouring epochs differ a little)
% the scores by hand
ctx = struct('contract', c);
E = r.root;
for q = 1:numel(r.leaves(kw).path) - 1
    [~, E] = evalc('pipecompare.run.Steps.run(r.leaves(kw).path{q}, E, ctx)');
end
before = E;
[~, E] = evalc('pipecompare.run.Steps.run(r.leaves(kw).path{end}, E, ctx)');
w = pipecompare.eval.Measure.epochWeights(E);
verifyNumElements(tc, w, E.trials);
verifyLessThan(tc, min(w), 0.5);                                % the artifact trials count less
times = E.xmin + (0:E.pnts - 1) / E.srate; win = times >= 0.3 - 1e-9 & times <= 0.5 + 1e-9;
roi = ismember(lower({E.chanlocs.labels}), lower({'Pz', 'P3', 'P4'}));
score = reshape(mean(mean(double(E.data(roi, win, :)), 1), 2), 1, []);
for k = 1:2
    sel = pipecompare.eval.Measure.lockingEvents(E, c.conditions(k).events) > 0;
    verifyEqual(tc, r.cands(kw).m.objectives(1).estimate(k), sum(w(sel) .* score(sel)) / sum(w(sel)), 'RelTol', 1e-6);
end
% the weighted average of each condition gives the same estimates
A = pipecompare.eval.Measure.weightedAverages(E, c);
verifyEqual(tc, size(A.data), [E.nbchan E.pnts 2]);
for k = 1:2
    verifyEqual(tc, mean(mean(A.data(roi, win, k), 1), 2), r.cands(kw).m.objectives(1).estimate(k), 'RelTol', 1e-6);
end
% adopting it stores them; the results say so, and name the unweighted twin
pipecompare.PipeCompare.adopt(r, kw, true);
adopted = evalin('base', 'EEG');
verifyEqual(tc, adopted.etc.pipecompare.weightedERP.data, A.data, 'AbsTol', 1e-4);
r.ranking.recommended = kw;
t = pipecompare.simple.Presets.weightsText(r);
verifyTrue(tc, contains(t, 'weightedERP') && contains(t, sprintf('pipeline %d', k0)), t);
r.ranking.recommended = k0;
verifyEmpty(tc, pipecompare.simple.Presets.weightsText(r));
% the recorded command weights the epochs again
com = r.cands(kw).coms{find(contains(r.cands(kw).coms, 'weightTrials'), 1)};
EEG = before; eval(com);
verifyEqual(tc, pipecompare.eval.Measure.epochWeights(EEG), w, 'RelTol', 1e-12);
% band power is not weighted
b = pipecompare.eval.Contract('analysis', 'bandpower', 'segment', 2, 'bands', {'alpha', [8 13], {'Pz'}});
verifyError(tc, @() pipecompare.run.Steps.run(r.leaves(kw).path{end}, before, struct('contract', b)), 'PipeCompare:Weights');
end

function testIclabelThresholdFromTheDataEndToEnd(tc)
% In a search with ICA: the threshold is chosen and listed in the steps;
% adopt rebuilds the same pipeline (same SME), and the saved script's
% steps, run again on the starting dataset, choose the same threshold.
assumeTrue(tc, exist('pop_iclabel', 'file') == 2, 'ICLabel not installed');
[EEG, truth] = nqc_synth(struct('seconds', 150, 'nPerCond', 30));
nqc_setBase(EEG);
c = pipecompare.eval.Contract('conditions', {'target', {'11'}; 'standard', {'31'}}, ...
    'epoch', [-0.2 1.0], 'baseline', [-0.2 0], 'components', {'P3', [0.30 0.50], truth.roi});
p = pipecompare.plan.Plan();
p = p.add('highpass', 'cutoff', 0.1); p = p.add('ica', 'fitHighpass', 1);
p = p.add('icremove', 'threshold', 'auto'); p = p.add('epoch'); p = p.add('baseline');
p = p.add('reject_threshold', 'uv', 'auto');
r = pipecompare.PipeCompare.optimize(p, c);
verifyEqual(tc, r.cands(1).status, 'ok', r.cands(1).message);
st = r.cands(1).steps{3};
verifyEqual(tc, st.type, 'icremove');
verifyTrue(tc, isfinite(st.thresholdChosen) && st.thresholdChosen > 0 && st.thresholdChosen <= 1);
verifyEqual(tc, r.ranking.byStratum.crossfit.nested, 1);         % its threshold chosen again in the selection check
L = pipecompare.simple.Presets.stepsText(r, 1);
verifyTrue(tc, any(contains(L, 'threshold chosen from these data')), strjoin(L, newline));
pipecompare.PipeCompare.adopt(r, 1, true);
adopted = evalin('base', 'EEG');
m = pipecompare.eval.Measure.candidate(adopted, c, r.ref);
verifyEqual(tc, m.composite, r.cands(1).m.composite, 'RelTol', 1e-9);
path = r.leaves(1).path;
steps = [cellfun(@(x) x.type, path, 'UniformOutput', false)' cellfun(@(x) x.params, path, 'UniformOutput', false)'];
[~, E3] = evalc('pipecompare.PipeCompare.apply(EEG, steps, c)');
verifyTrue(tc, contains(E3.history, 'pop_icflag'));
m3 = pipecompare.eval.Measure.candidate(E3, c, r.ref);
verifyEqual(tc, m3.composite, r.cands(1).m.composite, 'RelTol', 1e-6);
end

function testIcComponentsChosenBySmeEndToEnd(tc)
% In a search with ICA: the components are chosen by their effect on the
% SME and listed in the steps; adopt rebuilds the same pipeline (same
% SME) with the components the search chose.
assumeTrue(tc, exist('pop_iclabel', 'file') == 2, 'ICLabel not installed');
[EEG, truth] = nqc_synth(struct('seconds', 150, 'nPerCond', 30));
nqc_setBase(EEG);
c = pipecompare.eval.Contract('conditions', {'target', {'11'}; 'standard', {'31'}}, ...
    'epoch', [-0.2 1.0], 'baseline', [-0.2 0], 'components', {'P3', [0.30 0.50], truth.roi});
p = pipecompare.plan.Plan();
p = p.add('highpass', 'cutoff', 0.1); p = p.add('ica', 'fitHighpass', 1);
p = p.add('icremove', 'threshold', 'sme'); p = p.add('epoch'); p = p.add('baseline');
p = p.add('reject_threshold', 'uv', 'auto');
r = pipecompare.PipeCompare.optimize(p, c);
verifyEqual(tc, r.cands(1).status, 'ok', r.cands(1).message);
st = r.cands(1).steps{3};
verifyEqual(tc, st.type, 'icremove');
verifyTrue(tc, isfinite(st.icsConsidered));
verifyEqual(tc, r.ranking.byStratum.crossfit.nested, 1);         % its components chosen again in the selection check
L = pipecompare.simple.Presets.stepsText(r, 1);
verifyTrue(tc, any(contains(L, 'effect on the noise')), strjoin(L, newline));
pipecompare.PipeCompare.adopt(r, 1, true);
adopted = evalin('base', 'EEG');
m = pipecompare.eval.Measure.candidate(adopted, c, r.ref);
verifyEqual(tc, m.composite, r.cands(1).m.composite, 'RelTol', 1e-9);
verifyTrue(tc, contains(adopted.history, 'gcompreject'));
end

function testPeakToPeakRejection(tc)
% Peak-to-peak in moving 200 ms windows: a slow drift that crosses an
% absolute limit is kept, a fast deflection is rejected; the marks are
% per channel, so epochs can be repaired the same way.
EEG = nqc_synth(struct('seconds', 90, 'nPerCond', 20, 'artifactTrials', 0, 'blinkRate', 0));
[~, E] = evalc('nqc_epoched(EEG)');
ch = find(strcmpi({E.chanlocs.labels}, 'T7'));
E.data(ch, :, 5) = E.data(ch, :, 5) + linspace(0, 300, E.pnts);          % slow drift
t = (1:E.pnts) - round(0.5 * E.srate);
E.data(ch, :, 9) = E.data(ch, :, 9) + 150 * exp(-0.5 * (t / (0.03 * E.srate)) .^ 2);   % fast deflection
in = struct('type', 'reject_threshold', 'params', struct('uv', 100, 'exclude', {{}}), 'key', 'rej', 'slot', 'rej', 'label', 'rej');
x = struct('contract', nqc_c(), 'highpass', 0);
[~, ~, a] = pipecompare.run.Steps.run(in, E, x);
in.params.method = 'peaktopeak'; in.params.window = 200;
[~, coms, b] = pipecompare.run.Steps.run(in, E, x);
verifyTrue(tc, ismember(5, a.rejIdx) && ~ismember(5, b.rejIdx), mat2str(b.rejIdx));
verifyTrue(tc, ismember(9, b.rejIdx), mat2str(b.rejIdx));
verifyTrue(tc, contains(coms{1}, 'pipecompare.run.Steps.markPeakToPeak'));
in.params.interpolate = 3;                                       % repaired like any other test
[~, ~, c] = pipecompare.run.Steps.run(in, E, x);
verifyFalse(tc, ismember(9, c.rejIdx));
verifyGreaterThan(tc, c.epochsInterpolated, 0);
end

function testVariableP3DoesNotMisleadTheChoice(tc)
% Simulation: the P3 differs from trial to trial in size (SD 50%) and in
% latency (SD 50 ms), real activity that the SME counts as noise. Every
% pipeline that passes still measures each condition's P3 without bias:
% its estimate, divided by its signal gain, is within 3 standard errors
% of the P3 actually in the data (the mean, over all trials, of each
% trial's P3 at the measured electrodes in the window). Weighting the
% epochs included: their weights do not count the trials with a large
% P3 less. The recommended pipeline is one of them.
[EEG, truth] = nqc_synth(struct('seconds', 300, 'nPerCond', 60, 'p3AmpSd', 0.5, 'p3Jitter', 0.05));
nqc_setBase(EEG);
c = pipecompare.eval.Contract('conditions', {'target', {'11'}; 'standard', {'31'}}, ...
    'epoch', [-0.2 1.0], 'baseline', [-0.2 0], 'components', {'P3', [0.30 0.50], truth.roi});
p = pipecompare.plan.Plan();
p = p.add('highpass', 'cutoff', 0.1); p = p.add('lowpass', 'cutoff', {20, 40});
p = p.add('epoch'); p = p.add('baseline');
p = p.add('reject_threshold', 'uv', 'auto');
p = p.addChoice('weight', {'weight_trials'}, 'none');
r = pipecompare.PipeCompare.optimize(p, c);
% the P3 in the data, per trial: its mean over the ROI and the window
fs = EEG.srate; tt = (-round(0.2 * fs):round(1.0 * fs)) / fs;
win = tt >= 0.30 - 1e-9 & tt <= 0.50 + 1e-9;
roi = ismember(upper({EEG.chanlocs.labels}), upper(truth.roi));
g = arrayfun(@(k) mean(exp(-0.5 * ((tt(win) - truth.p3Latency - truth.p3Shift(k)) / truth.p3Sigma) .^ 2)), 1:numel(truth.p3Amp));
p3 = mean(truth.p3Topo(roi)) * truth.p3Amp .* g;
want = [mean(p3(strcmp(truth.types, '11'))) mean(p3(strcmp(truth.types, '31')))];
T = r.ranking.table;
feas = find(strcmp(T.status, 'feasible'))';
verifyNotEmpty(tc, feas);
verifyTrue(tc, ismember(r.ranking.recommended, feas));
for k = feas
    o = r.cands(k).m.objectives(1); gk = r.cands(k).signal.gain(1);
    err = abs(o.estimate / gk - want) ./ (o.sme / gk);
    fprintf('variable P3, pipeline %d (%s): estimate %s uV, P3 in the data %s uV, %s standard errors off\n', ...
        k, r.labels{k}, mat2str(o.estimate / gk, 3), mat2str(want, 3), mat2str(err, 2));
    verifyLessThan(tc, err, 3, r.labels{k});
end
end

function testBlinksCountInTheTotalError(tc)
% Simulation: in 30% of the target trials a blink time-locked to the
% stimulus falls in the window, measured at Cz (a blink at FPz reaches Pz
% with well under 1% of its size, Cz with a few %). It moves the target mean much more than
% it raises the SME. The estimate of what the blinks add (eyeBias) is
% close to what they add in fact (the same pipeline on the same data
% without those blinks: everything else is drawn the same), a pipeline
% that rejects the blink trials adds almost nothing, and it is
% recommended.
o = struct('seconds', 300, 'nPerCond', 80, 'lockedBlinks', 0.3, 'lockedBlinkUv', 500, 'artifactTrials', 0, 'blinkRate', 0);
EEG = nqc_synth(o);
o.lockedBlinks = 0; EEG0 = nqc_synth(o);
nqc_setBase(EEG);
c = pipecompare.eval.Contract('conditions', {'target', {'11'}; 'standard', {'31'}}, ...
    'epoch', [-0.2 1.0], 'baseline', [-0.2 0], 'components', {'P3', [0.30 0.50], {'Cz'}});
p = pipecompare.plan.Plan();
p = p.add('epoch'); p = p.add('baseline');
p = p.add('reject_threshold', 'uv', {1000, 100}, 'exclude', {{'EOG1', 'EOG2'}});
r = pipecompare.PipeCompare.optimize(p, c);
verifyEqual(tc, sort(r.ref.blinkInfo.channels), {'EOG1', 'EOG2'});
T = r.ranking.table;
kAll = find(contains(r.labels, '1000')); kRej = setdiff(1:2, kAll);
verifyTrue(tc, all(strcmp(T.status, 'feasible')), strjoin(T.reason, '; '));
% what the blinks add in fact, to the target mean of the pipeline keeping every trial
[~, E1] = evalc('nqc_epoched(EEG)'); [~, E0] = evalc('nqc_epoched(EEG0)');
m1 = pipecompare.eval.Measure.candidate(E1, c, r.ref); m0 = pipecompare.eval.Measure.candidate(E0, c, r.ref);
actual = m1.objectives(1).estimate(1) - m0.objectives(1).estimate(1);
o1 = r.cands(kAll).m.objectives(1);
eb = T.(matlab.lang.makeValidName(['eyeBias_' r.ref.objectives{1}]));   % the table column of the P3's eye bias
x = o1.X{1}(:, 1); b = r.ref.blink{1}{1};
in = ~isnan(x); x = x(in); b = b(in); pb = mean(b);            % the trials it scored
verifyGreaterThan(tc, sum(b), 10);                             % the blink trials are found
se = pb * std(x) * sqrt(1 / sum(b) + 1 / sum(~b));
fprintf(['blinks: they add %.2f uV to the target mean; estimated %.2f uV (SE %.2f); rejecting them: %.2f uV; ', ...
    'noise %.3f keeping them, %.3f rejecting them\n'], actual, eb(kAll), se, eb(kRej), T.noise(kAll), T.noise(kRej));
verifyGreaterThan(tc, actual, 2 * se);                            % a real bias
verifyLessThan(tc, abs(eb(kAll) - actual), 3 * se);
verifyLessThan(tc, abs(eb(kRej)), abs(eb(kAll)) / 4);
verifyEqual(tc, r.ranking.recommended, kRej);
r2 = r; r2.ranking.recommended = kAll;                            % the results' words for the other one
verifyTrue(tc, contains(pipecompare.simple.Presets.biasText(r2), 'muscle activity in the trials it keeps'), pipecompare.simple.Presets.biasText(r2));
% the total error: noise and bias together
verifyEqual(tc, T.objective, sqrt(T.noise .^ 2 + T.bias .^ 2), 'RelTol', 1e-9);
end

function testEyeMovementsAndMuscleCountInTheBias(tc)
% Simulation: in 30% of the target trials the eyes move to the left
% during the window (left minus right EOG), and in 15% of all trials a
% muscle burst falls in it, measured at F3. The eye movements are found
% in the difference of the left and right EOG channels (where a blink
% cancels), the muscle bursts from the fast activity at F3; what the
% trials with them add to the target mean is estimated close to what the
% eye movements add in fact (the same data without them).
o = struct('seconds', 300, 'nPerCond', 80, 'lockedSaccades', 0.3, 'lockedSaccadeUv', 100, ...
    'lockedEmg', 0.15, 'lockedEmgUv', 30, 'artifactTrials', 0, 'blinkRate', 0);
[EEG, truth] = nqc_synth(o);
o.lockedSaccades = 0; EEG0 = nqc_synth(o);
nqc_setBase(EEG);
c = pipecompare.eval.Contract('conditions', {'target', {'11'}; 'standard', {'31'}}, ...
    'epoch', [-0.2 1.0], 'baseline', [-0.2 0], 'components', {'P3', [0.30 0.50], {'F3'}});
p = pipecompare.plan.Plan();
p = p.add('epoch'); p = p.add('baseline');
p = p.add('reject_threshold', 'uv', {1000, 150}, 'exclude', {{'EOG1', 'EOG2'}});
r = pipecompare.PipeCompare.optimize(p, c);
bi = r.ref.blinkInfo;
verifyEqual(tc, bi.horizontal, 'EOG2-EOG1');
verifyEqual(tc, bi.kinds, {'blinks', 'eye movements', 'muscle'});
% the trials found, against the trials simulated (both conditions in event order)
tgt = strcmp(truth.types, '11');
b = false(1, numel(tgt)); b(tgt) = r.ref.blink{1}{1}; b(~tgt) = r.ref.blink{1}{2};
verifyGreaterThan(tc, mean(b(truth.lockedSaccade)), 0.8);
verifyGreaterThan(tc, mean(b(truth.lockedEmg)), 0.8);
verifyLessThan(tc, mean(b(~truth.lockedSaccade & ~truth.lockedEmg)), 0.1);
verifyGreaterThan(tc, bi.kindShare(1, 2), 0.1);
verifyGreaterThan(tc, bi.kindShare(1, 3), 0.1);
% what they add to the target mean of the pipeline keeping every trial
T = r.ranking.table;
kAll = find(contains(r.labels, '1000'));
[~, E1] = evalc('nqc_epoched(EEG)'); [~, E0] = evalc('nqc_epoched(EEG0)');
m1 = pipecompare.eval.Measure.candidate(E1, c, r.ref); m0 = pipecompare.eval.Measure.candidate(E0, c, r.ref);
actual = m1.objectives(1).estimate(1) - m0.objectives(1).estimate(1);
eb = T.(matlab.lang.makeValidName(['eyeBias_' r.ref.objectives{1}]));
x = r.cands(kAll).m.objectives(1).X{1}(:, 1); bb = r.ref.blink{1}{1};
in = ~isnan(x); x = x(in); bb = bb(in);
se = mean(bb) * std(x) * sqrt(1 / sum(bb) + 1 / sum(~bb));
fprintf('eye movements: they add %.2f uV to the target mean at F3; estimated %.2f uV (SE %.2f)\n', actual, eb(kAll), se);
verifyGreaterThan(tc, abs(actual), 2 * se);
verifyLessThan(tc, abs(eb(kAll) - actual), 3 * se);
end

function testShapeSpreadGrowsWithTheHighPass(tc)
% How much a pipeline's filters would change the gain for another signal
% shape: none without a filter, small for a 0.1 Hz high-pass, larger at
% 1 Hz; a low-pass at 30 Hz changes the broad P3 shapes little.
EEG = nqc_synth(struct('seconds', 120, 'nPerCond', 20));
c = nqc_c();
ref = pipecompare.eval.Measure.reference(EEG, c);
S = pipecompare.eval.Injection.prepare(EEG, c, ref);
ins = @(type, v) struct('type', type, 'params', struct('cutoff', v), 'key', sprintf('%s%g', type, v), 'slot', type, 'label', type);
leaves = struct('path', {{ins('epoch', 0)}, {ins('highpass', 0.1)}, {ins('highpass', 1)}, {ins('lowpass', 30)}});
sp = pipecompare.eval.Injection.shapeSpreads(S, c, struct('kind', 'erp'), leaves);
sp = [sp{:}];
fprintf('shape spread: none %.3f, HP 0.1 %.3f, HP 1 %.3f, LP 30 %.3f\n', sp);
verifyEqual(tc, sp(1), 0);
verifyLessThan(tc, sp(2), 0.02);
verifyGreaterThan(tc, sp(3), 0.05);
verifyGreaterThan(tc, sp(3), 5 * sp(2));
verifyLessThan(tc, sp(4), 0.05);
end

function testRejectingEpochsForTheirLargeComponentIsReported(tc)
% A limit low enough to reject epochs for their large P3 at the measured
% electrodes keeps epochs with less P3 than all of them: the rejection
% reports it (keptShift) and the results warn. Rejecting the same number
% of epochs for a deflection at a far channel does not.
EEG = nqc_synth(struct('seconds', 90, 'nPerCond', 20, 'artifactTrials', 0, 'blinkRate', 0));
[~, E] = evalc('nqc_epoched(EEG)');
c = nqc_c(); x = struct('contract', c, 'highpass', 0);
labs = lower({E.chanlocs.labels});
roi = ismember(labs, lower(c.components(1).roi));
eeg = ~pipecompare.simple.Presets.isNonEeg(E.chanlocs(:)');
uv = ceil(max(abs(double(E.data(:))))) + 1;                  % over the limit: none yet
times = E.xmin + (0:E.pnts - 1) / E.srate; win = times >= 0.3 & times <= 0.5;
T = pipecompare.eval.Measure.trials(E, c, true);
ep = T.epoch(T.cond == 1); ep = ep(1:6);
in = struct('type', 'reject_threshold', 'params', struct('uv', uv, 'exclude', {{}}), 'key', 'rej', 'slot', 'rej', 'label', 'rej');
A = E; A.data(roi, win, ep) = A.data(roi, win, ep) + uv + 20;    % a large P3 in 6 target epochs
[~, ~, a] = pipecompare.run.Steps.run(in, A, x);
verifyEqual(tc, sort(a.rejIdx(:))', sort(ep(:))');
verifyLessThan(tc, a.keptShift.z, -pipecompare.simple.Presets.KeptShiftZ);
verifyEqual(tc, a.keptShift.condition, 't');
r = struct('ranking', struct('recommended', 1), 'cands', struct('id', 1, 'steps', {{a}}), 'ref', struct('units', {{'uV'}}));
t = pipecompare.simple.Presets.keptShiftText(r);
verifyTrue(tc, contains(t, 'had more P3 than those it kept'), t);
far = find(eeg & ~roi, 1, 'last');
B = E; B.data(far, win, ep) = B.data(far, win, ep) + uv + 20;    % the same epochs, for a far channel
[~, ~, b] = pipecompare.run.Steps.run(in, B, x);
verifyEqual(tc, sort(b.rejIdx(:))', sort(ep(:))');
verifyGreaterThan(tc, b.keptShift.z, -pipecompare.simple.Presets.KeptShiftZ);
r.cands.steps = {b};
verifyEmpty(tc, pipecompare.simple.Presets.keptShiftText(r));
end

function testIcaIsFittedWithoutExtremeStretches(tc)
% The copy ICA is fitted on leaves out windows far noisier than the rest
% (at most 20%), and nothing on clean data.
EEG = nqc_synth(struct('seconds', 90, 'nPerCond', 20, 'artifactTrials', 0, 'blinkRate', 0));
[T, left] = pipecompare.run.Steps.cleanForIca(EEG);
verifyLessThan(tc, left, 0.05);
fs = EEG.srate; eeg = ~pipecompare.simple.Presets.isNonEeg(EEG.chanlocs(:)');
rng(2);
for k = [10 30 50 70 80]
    EEG.data(eeg, (k - 1) * fs + (1:fs)) = EEG.data(eeg, (k - 1) * fs + (1:fs)) + 300 * randn(sum(eeg), fs);
end
[T, left] = pipecompare.run.Steps.cleanForIca(EEG);
verifyGreaterThanOrEqual(tc, left, 5 / 90 - 1e-9);
verifyLessThanOrEqual(tc, T.pnts, EEG.pnts - 5 * fs);
verifyLessThan(tc, max(abs(double(T.data(eeg, :))), [], 'all'), 1000);   % the bursts are gone
end

function testRepairedEpochsInASearch(tc)
% In a search, repairing epochs is compared like any setting; its
% decisions are replayed on the signal copy and listed in the steps.
EEG = nqc_synth(struct('seconds', 150, 'nPerCond', 30));
% T7 loses contact after every 4th event: far over the limit in those
% epochs only (a channel noisy all the time would fail every epoch)
ch = strcmpi({EEG.chanlocs.labels}, 'T7');
for e = 1:4:numel(EEG.event)
    t = round(EEG.event(e).latency) + (0:round(0.5 * EEG.srate));
    EEG.data(ch, t) = EEG.data(ch, t) + 400;
end
nqc_setBase(EEG);
p = pipecompare.plan.Plan();
p = p.add('highpass', 'cutoff', 0.5);
p = p.add('epoch'); p = p.add('baseline');
p = p.add('reject_threshold', 'uv', 100, 'interpolate', {0, 3});
r = pipecompare.PipeCompare.optimize(p, nqc_c());
verifyEqual(tc, {r.cands.status}, {'ok', 'ok'});
verifyEmpty(tc, [r.cands.unmatched]);                          % decisions replayed, not re-run
k = find(contains(r.labels, 'interpolate=3'));
f = r.cands(k).steps{end};
verifyGreaterThan(tc, f.epochsInterpolated, 0);
verifyLessThan(tc, r.cands(k).rejectedEpochs, r.cands(3 - k).rejectedEpochs);
L = pipecompare.simple.Presets.stepsText(r, k);
verifyTrue(tc, any(contains(L, 'kept by interpolating up to 3 channel(s)')), strjoin(L, newline));
end

function E = withIca(E, d)
% the decomposition d (as the search stores it) on dataset E
E.icaweights = d.icaweights; E.icasphere = d.icasphere; E.icachansind = d.icachansind;
E.icawinv = []; E.icaact = [];
E = eeg_checkset(E);
end

function testEyeChannelsRegressedOut(tc)
% Steps.eogRegress: blinks locked to the P3 (only they differ between the
% two recordings) are taken out of the scalp channels by the coefficients
% fitted on the data with them; the eye channels stay as they are and the
% signal-check copy gets the same coefficients. In the search (eogCheck)
% the recommended pipeline also runs with the eye channels regressed out,
% which here, measured where the blinks reach, lowers its total error.
o = struct('seconds', 150, 'nPerCond', 30, 'artifactTrials', 0, 'lockedBlinkUv', 600);
E1 = nqc_synth(setfield(o, 'lockedBlinks', 0.6)); E0 = nqc_synth(setfield(o, 'lockedBlinks', 0)); %#ok<SFLD>
c = nqc_c();
[R1, info] = pipecompare.run.Steps.eogRegress(E1, {}, {c.conditions.events}, c.epoch);
verifyEqual(tc, sort(info.eyes), {'EOG1', 'EOG2'});
labs = {E1.chanlocs.labels}; ei = ismember(labs, info.eyes); ti = ismember(labs, info.targets);
verifyEqual(tc, R1.data(ei, :), E1.data(ei, :));                  % the eye channels as they are
D = double(E1.data - E0.data);                                    % the locked blinks alone
left = D(ti, :) - info.B * D(ei, :);
fz = strcmp(labs(ti), 'Fz');
verifyLessThan(tc, max(abs(left(fz, :))), 0.1 * max(abs(D(strcmp(labs, 'Fz'), :))));
in = struct('type', 'eogregress', 'params', struct('channels', {{}}), 'key', 'eogregress', 'slot', 'eog', 'label', 'eog');
S = E1; S.data = single(randn(size(E1.data)));
S2 = pipecompare.run.Steps.replayDecision(in, S, struct('eog', info), struct('contract', c));
verifyEqual(tc, double(S2.data(ti, :)), double(S.data(ti, :)) - info.B * double(S.data(ei, :)), 'AbsTol', 1e-4);
[~, ~, com] = evalc('pipecompare.run.Steps.run(in, E1, struct(''contract'', c))');
verifyTrue(tc, contains(com{1}, 'pipecompare.run.Steps.eogRegress(EEG'));
% in the search, measured at Cz, which the blinks reach (Pz they hardly do)
c = pipecompare.eval.Contract('conditions', {'t', {'11'}; 's', {'31'}}, 'epoch', [-0.2 1], ...
    'baseline', [-0.2 0], 'components', {'P3', [0.3 0.5], {'Cz'}});
nqc_setBase(E1);
p = pipecompare.plan.Plan();
p = p.add('highpass', 'cutoff', {0.1, 0.5}); p = p.add('lowpass', 'cutoff', 30);
p = p.add('epoch'); p = p.add('baseline');
r = pipecompare.PipeCompare.optimize(p, c, struct('eogCheck', true));
k = find(arrayfun(@(l) any(cellfun(@(i) strcmp(i.type, 'eogregress'), l.path)), r.leaves));
verifyNumElements(tc, k, 1);
verifyEqual(tc, r.cands(k).status, 'ok');
base = sscanf(char(regexp(strjoin(r.cands(k).notes, ' '), 'pipeline (\d+) with the eye', 'tokens', 'once')), '%d');
verifyLessThan(tc, r.ranking.table.objective(k), r.ranking.table.objective(base));
verifyTrue(tc, any(contains(pipecompare.simple.Presets.stepsText(r, k), 'Eye channels EOG')));
r0 = pipecompare.PipeCompare.optimize(p, c);
verifyEqual(tc, numel(r0.cands), 2);                              % not by default
end

function c = nqc_c()
c = pipecompare.eval.Contract('conditions', {'t', {'11'}; 's', {'31'}}, 'epoch', [-0.2 1], ...
    'baseline', [-0.2 0], 'components', {'P3', [0.3 0.5], {'Pz','P3','P4'}});
end

% ------------------------------------------------------------ band power
function testBandPowerOfContinuousDataEndToEnd(tc)
% Continuous alpha power: a 9 Hz low-pass destroys the band and is
% rejected by the signal check; 30 Hz keeps it. Segments are marked with
% eeg_regepochs and paired across candidates by urevent.
EEG = nqc_synth(struct('seconds', 150, 'nPerCond', 10, 'alphaUv', 10, 'artifactTrials', 0));
nqc_setBase(EEG);
c = pipecompare.eval.Contract('analysis', 'bandpower', 'segment', 2, 'bands', {'alpha', [8 12], {'Oz', 'O1', 'O2'}});
p = pipecompare.plan.Plan(); p = p.add('highpass', 'cutoff', 1); p = p.add('lowpass', 'cutoff', {9, 30});
p = p.add('epoch');
r = pipecompare.PipeCompare.optimize(p, c);
lp9 = contains(r.labels, 'cutoff=9'); lp30 = contains(r.labels, 'cutoff=30');
verifyEqual(tc, r.ranking.table.status{lp9}, 'rejected');
verifyTrue(tc, contains(r.ranking.table.reason{lp9}, 'amplitude'));
verifyEqual(tc, r.ranking.table.status{lp30}, 'feasible');
verifyGreaterThanOrEqual(tc, r.ref.n, 70);                        % 150 s / 2 s, minus edges
verifyEqual(tc, r.ref.units{1}, 'log10(uV^2)');
verifyTrue(tc, r.ref.segmented);
verifyTrue(tc, any(contains(r.rootComs, 'eeg_regepochs')));       % pure EEGLAB, in the script
end

function testBandPowerContractIsValidated(tc)
EEG = nqc_synth(struct('seconds', 60, 'nPerCond', 10));
st = pipecompare.live.DataState.fromEEG(EEG);
short = pipecompare.eval.Contract('analysis', 'bandpower', 'segment', 0.1, 'bands', {'theta', [4 8], {'Fz'}});
verifyError(tc, @() short.validate(st), 'PipeCompare:Contract');      % fewer than two 4 Hz cycles
high = pipecompare.eval.Contract('analysis', 'bandpower', 'segment', 2, 'bands', {'gamma', [100 140], {'Fz'}});
verifyError(tc, @() high.validate(st), 'PipeCompare:Contract');       % above Nyquist (125 Hz)
ok = pipecompare.eval.Contract('analysis', 'bandpower', 'segment', 2, 'bands', {'theta', [4 8], {'Fz'}});
ok.validate(st);
[~, Ep] = evalc('pop_epoch(EEG, {''11''}, [-0.2 0.8])');
verifyError(tc, @() ok.validate(pipecompare.live.DataState.fromEEG(Ep)), 'PipeCompare:Contract');   % segments need continuous data
% event-related band power: a baseline outside the epoch is caught before the run
ev = {'analysis', 'bandpower', 'conditions', {'t', {'11'}}, 'epoch', [-0.5 1.5], 'bands', {'theta', [4 8], {'Fz'}}};
noBase = pipecompare.eval.Contract(ev{:}); noBase.validate(st);
inside = pipecompare.eval.Contract(ev{:}, 'baseline', [-0.5 0]); inside.validate(st);
outside = pipecompare.eval.Contract(ev{:}, 'baseline', [-1 0]);
verifyError(tc, @() outside.validate(st), 'PipeCompare:Contract');
reversed = pipecompare.eval.Contract(ev{:}, 'baseline', [0 -0.2]);
verifyError(tc, @() reversed.validate(st), 'PipeCompare:Contract');
% segment together with conditions: one of them would be ignored, so it is an error
verifyError(tc, @() pipecompare.eval.Contract(ev{:}, 'segment', 2), 'PipeCompare:Contract');
end

function testEventRelatedBandPower(tc)
% With conditions and epoch instead of segments, each epoch's band power
% is scored; trials are independent (no block resampling).
EEG = nqc_synth(struct('seconds', 120, 'nPerCond', 30, 'alphaUv', 10));
nqc_setBase(EEG);
c = pipecompare.eval.Contract('analysis', 'bandpower', 'conditions', {'t', {'11'}; 's', {'31'}}, ...
    'epoch', [-0.5 1.5], 'bands', {'alpha', [8 12], {'Oz', 'O1', 'O2'}});
p = pipecompare.plan.Plan(); p = p.add('highpass', 'cutoff', {0.5, 1}); p = p.add('epoch');
r = pipecompare.PipeCompare.optimize(p, c);
verifyFalse(tc, r.ref.segmented);
verifyEqual(tc, r.ref.n, [30 30]);
verifyTrue(tc, all(strcmp(r.ranking.table.status, 'feasible')));
end

function testBandPowerScriptKeepsTheContract(tc)
% Audit: writeScript wrote the default (ERP) contract for band power, so
% the saved script did not segment the data and its epoch step failed.
EEG = nqc_synth(struct('seconds', 90, 'nPerCond', 20, 'alphaUv', 10, 'artifactTrials', 0));
nqc_setBase(EEG);
cs = {pipecompare.eval.Contract('analysis', 'bandpower', 'segment', 2, 'bands', {'alpha', [8 12], {'Oz', 'O1', 'O2'}}), ...
      pipecompare.eval.Contract('analysis', 'bandpower', 'conditions', {'t', {'11'}; 's', {'31'}}, ...
          'epoch', [-0.5 1.5], 'bands', {'alpha', [8 12], {'Oz', 'O1', 'O2'}})};
d = tempname; mkdir(d); cleanup = onCleanup(@() rmdir(d, 's')); %#ok<NASGU>
addpath(d); c2 = onCleanup(@() rmpath(d)); %#ok<NASGU>
for i = 1:numel(cs)
    p = pipecompare.plan.Plan(); p = p.add('highpass', 'cutoff', 1); p = p.add('epoch');
    r = pipecompare.PipeCompare.optimize(p, cs{i});
    name = sprintf('nqc_band_script%d', i);
    pipecompare.PipeCompare.writeScript(r, 1, fullfile(d, [name '.m']));
    out = feval(name, EEG);
    m = pipecompare.eval.Measure.candidate(out, r.contract, r.ref);
    verifyEqual(tc, [m.objectives.agg], [r.cands(1).m.objectives.agg], 'RelTol', 1e-6);
end
end

function testErpScriptKeepsTheComponents(tc)
% NQC-034: writeScript left the components out of an ERP contract, so the
% contract in the saved script could not be scored again.
EEG = nqc_synth(struct('seconds', 90, 'nPerCond', 20, 'artifactTrials', 0));
nqc_setBase(EEG);
base = {'epoch', [-0.2 1], 'baseline', [-0.2 0]};
cs = {pipecompare.eval.Contract('conditions', {'t', {'11'}; 's', {'31'}}, base{:}, 'components', ...
          {'P3', [0.3 0.5], {'Pz', 'P3', 'P4'}, {'peakAmplitude', 'negative'}; 'N1', [0.08 0.12], {'Cz'}, 'mean'}), ...
      pipecompare.eval.Contract('conditions', {'l', {'11'}; 'r', {'31'}}, base{:}, 'components', ...
          {'LAT', [0.2 0.4], {'P3', 'P4'}, 'mean', {'P4', 'P3'}})};   % contralateral minus ipsilateral
d = tempname; mkdir(d); cleanup = onCleanup(@() rmdir(d, 's')); %#ok<NASGU>
for i = 1:numel(cs)
    p = pipecompare.plan.Plan(); p = p.add('highpass', 'cutoff', 0.1); p = p.add('epoch'); p = p.add('baseline');
    r = pipecompare.PipeCompare.optimize(p, cs{i});
    f = fullfile(d, sprintf('nqc_erp_script%d.m', i));
    pipecompare.PipeCompare.writeScript(r, 1, f);
    L = strtrim(strsplit(fileread(f), newline));
    contract = []; eval(L{startsWith(L, 'contract = ')});
    verifyEqual(tc, contract.toStruct(), r.contract.toStruct());
    contract.validate(r.state);
end
end

function testLongHighPassRunsInTheFrequencyDomain(tc)
% A 0.1 Hz high-pass (an 8251-point FIR at 250 Hz) uses firfilt's fftfilt
% option when the Signal Processing Toolbox is there: the same filter, so
% the same data as EEGLAB's default time-domain filtering.
EEG = nqc_synth(struct('seconds', 60, 'nPerCond', 10));
in = struct('type', 'highpass', 'params', struct('cutoff', 0.1));
[~, E, coms] = evalc('pipecompare.run.Steps.run(in, EEG, struct())');
[~, D] = evalc('pop_eegfiltnew(EEG, ''locutoff'', 0.1, ''plotfreqz'', 0)');
verifyEqual(tc, double(E.data), double(D.data), 'AbsTol', 1e-3);   % uV, single-precision data
verifyEqual(tc, contains(coms{1}, 'usefftfilt'), exist('fftfilt', 'file') == 2);
% a short FIR (1 Hz: 826 points) stays in the time domain
in.params.cutoff = 1;
[~, ~, coms] = evalc('pipecompare.run.Steps.run(in, EEG, struct())');
verifyFalse(tc, contains(coms{1}, 'usefftfilt'));
end

function testCheckpointWithAProgressWindowCanBeResumed(tc)
% The panel passes its progress window's callback with a checkpoint
% folder: it is not stored and not part of the search identity.
EEG = nqc_synth(struct('seconds', 60, 'nPerCond', 12));
nqc_setBase(EEG);
p = pipecompare.plan.Plan();
p = p.add('highpass', 'cutoff', {0.5, 1}); p = p.add('epoch'); p = p.add('baseline');
d = tempname; cleanup = onCleanup(@() rmdir(d, 's')); %#ok<NASGU>
w = pipecompare.gui.Progress(2, false); c2 = onCleanup(@() delete(w)); %#ok<NASGU>
verifyError(tc, @() pipecompare.PipeCompare.optimize(p, nqc_c(), struct('checkpoint', d, 'stopAfter', 1, ...
    'progress', @(n) w.step(n), 'filterCheck', false)), 'PipeCompare:Interrupted');   % both run
M = load(fullfile(d, 'manifest.mat'));
verifyEmpty(tc, M.result.options.progress);
res = pipecompare.PipeCompare.resume(d);
verifyTrue(tc, all(strcmp({res.cands.status}, 'ok')));
end

function testIcaMatchesPopRunica(tc)
% EEGLAB before 2025.1 opens runica's Interrupt window for every
% pop_runica call; there PipeCompare fits ICA itself, as pop_runica does
% (same rank, weights and component order), without that window.
EEG = nqc_synth(struct('seconds', 60, 'nPerCond', 10, 'artifactTrials', 0));
[~, EEG] = evalc('pop_reref(EEG, [])');   % rank deficient: PCA to one dimension less
in = struct('type', 'ica', 'params', struct('extended', 1, 'fitHighpass', 0, 'fitClean', 0), ...
    'key', 'ica', 'slot', 'ica', 'label', 'ica');
[~, A] = evalc('pipecompare.run.Steps.run(in, EEG, struct())');
[~, B] = evalc('pop_runica(EEG, ''icatype'', ''runica'', ''extended'', 1, ''rndreset'', ''no'', ''interrupt'', ''off'')');
verifyEqual(tc, size(A.icaweights), size(B.icaweights));
verifyEqual(tc, A.icaweights, B.icaweights, 'AbsTol', 1e-8);
verifyEqual(tc, A.icasphere, B.icasphere, 'AbsTol', 1e-8);
verifyEqual(tc, double(A.icawinv), double(B.icawinv), 'AbsTol', 1e-6);
end

function testIcaStartedDifferentlyChecksTheRecommendation(tc)
% Another seed: runica sees the samples in another order, so the
% decomposition differs, and the same seed gives the same one. In the
% search, the pipelines deciding the recommendation run again with each
% seed and the result says whether the recommendation stays.
EEG = nqc_synth(struct('seconds', 60, 'nPerCond', 10, 'artifactTrials', 0));
in = struct('type', 'ica', 'params', struct('extended', 1, 'fitHighpass', 0, 'fitClean', 0), ...
    'key', 'ica', 'slot', 'ica', 'label', 'ica');
[~, A] = evalc('pipecompare.run.Steps.run(in, EEG, struct())');
[~, B, coms] = evalc('pipecompare.run.Steps.run(in, EEG, struct(''icaSeed'', 1))');
[~, C] = evalc('pipecompare.run.Steps.run(in, EEG, struct(''icaSeed'', 1))');
verifyGreaterThan(tc, max(abs(A.icaweights(:) - B.icaweights(:))), 1e-6);
verifyEqual(tc, B.icaweights, C.icaweights);
verifyTrue(tc, contains(strjoin(coms, ' '), 'random order (seed 1)'));
assumeTrue(tc, exist('pop_iclabel', 'file') == 2, 'ICLabel not installed');
EEG = nqc_synth(struct('seconds', 150, 'nPerCond', 30));
nqc_setBase(EEG);
p = pipecompare.plan.Plan();
p = p.add('reref', 'mode', 'average'); p = p.add('ica', 'fitHighpass', 1);
p = p.add('highpass', 'cutoff', {0.1, 0.5}); p = p.add('lowpass', 'cutoff', 30);
p = p.add('icremove', 'threshold', 0.8);
p = p.add('epoch'); p = p.add('baseline'); p = p.add('reject_threshold', 'uv', 100);
r = pipecompare.PipeCompare.optimize(p, nqc_c(), struct('icaSeeds', 1));
q = r.icaSeeds;
verifyNotEmpty(tc, q);
verifyEqual(tc, q.seeds, 1);
verifyEqual(tc, size(q.objective), [numel(q.ids) 2]);
verifyEqual(tc, q.recommended(1), r.ranking.recommended);
verifyTrue(tc, all(isfinite(q.objective(:))));
verifyEqual(tc, q.same, q.recommended(2) == q.recommended(1));
verifyTrue(tc, contains(pipecompare.simple.Presets.precisionText(r), 'With ICA started in'));
r0 = pipecompare.PipeCompare.optimize(p, nqc_c());
verifyEmpty(tc, r0.icaSeeds);
end

function testAsrMemoryFitsManyChannels(tc)
% clean_rawdata before 2.8 (EEGLAB 2022) stops with "Not enough memory"
% when MaxMem cannot hold ASR's lookahead buffer: it is raised with the
% number of channels and the sampling rate, and stays 64 MB otherwise.
verifyEqual(tc, pipecompare.run.AsrRecord.maxMemMB(64, 1000), 64);
for cf = [64 1000; 160 1000; 256 2048]'
    C = cf(1); P = round(max(0.5, 1.5 * C / cf(2)) / 2 * cf(2));
    verifyGreaterThan(tc, pipecompare.run.AsrRecord.maxMemMB(C, cf(2)) * 2^20, C * C * P * 24);
end
end
