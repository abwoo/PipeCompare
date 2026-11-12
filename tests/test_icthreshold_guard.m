function tests = test_icthreshold_guard
% Regressions for PR51 real ERP CORE subject 18/23: non-finite SME must
% be a scored infeasible option, not a MATLAB indexed-assignment crash.
tests = functiontests(localfunctions);
end

function setupOnce(~)
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
end

function testNonfiniteObjectiveReturnsScalarInfeasibleScore(tc)
% Start with genuine two-condition epoched trials, then make every Pz
% measurement undefined. Count/retention requirements still pass:
% Threshold.curve can report ok=true even when every obj is NaN.
[E, contract, ref, opts] = fixture();
pz = strcmpi({E.chanlocs.labels}, 'Pz');
verifyEqual(tc, nnz(pz), 1);
E.data(pz, :, :) = NaN;
[f, limit, why] = pipecompare.eval.IcThreshold.score(E, [], [], [], contract, ref, opts, []);
[v, ok, obj] = pipecompare.eval.Threshold.curve(E, zeros(1,E.trials), contract, ref, opts, [], 0);
fprintf('PR51_GUARD_DEBUG f=%s limit=%s why=%s v=%s ok=%s obj=%s\n', ...
    mat2str(f), mat2str(limit), why, mat2str(v), mat2str(ok), mat2str(obj));
verifyTrue(tc, isscalar(f) && isnan(f), 'NaN objective must not produce an empty score.');
verifyTrue(tc, isscalar(limit) && isnan(limit), 'NaN objective must not produce an empty cutoff.');
verifyTrue(tc, contains(why, 'non-finite'), 'Diagnostic must distinguish invalid objective from too few trials.');
end

function testAutoRejectThresholdRejectsNonfiniteObjectiveAsEvidence(tc)
% The default rejection search can retain enough trials while all P3
% scores are invalid. Never fit a rejection threshold to NaN as zero SME.
[E, contract, ref, opts] = fixture();
pz = strcmpi({E.chanlocs.labels}, 'Pz');
E.data(pz, :, :) = NaN;
m = 1:E.trials;
ctx = struct('contract', contract, 'ref', ref, 'rank', opts);
[limit, info] = pipecompare.eval.Threshold.choose(E, m, ctx);
verifyEqual(tc, info.feasible, 0);
verifyGreaterThan(tc, info.limitsWithEnoughTrials, 0);
verifyGreaterThanOrEqual(tc, limit, max(m)); % keep all data, do not optimize missing scores
verifyTrue(tc, contains(info.note, 'non-finite'));
end

function testIcaTextHandlesNoRecommendedPipeline(tc)
% A legitimate no-feasible-candidate outcome has ranking.recommended = [].
% Report formatting must not crash after the search has completed.
step = struct('type', 'icremove', 'icsTotal', 12, 'icsRemoved', 0, ...
    'icsBrain', 3, 'icsOther', 1, 'icaPoints', 25000, ...
    'icsConsidered', 4, 'otherMedian', 0.1);
candidate = struct('id', 1, 'steps', {{step}});
result = struct('cands', candidate, 'ranking', struct('recommended', []));
t = pipecompare.simple.Presets.icaText(result);
verifyTrue(tc, ischar(t));
verifyNotEmpty(tc, t);
end

function testFiniteObjectiveStillChoosesKnownLimit(tc)
[E, contract, ref, opts] = fixture();
[f, limit, why] = pipecompare.eval.IcThreshold.score(E, [], [], [], contract, ref, opts, []);
verifyTrue(tc, isscalar(f) && isfinite(f));
verifyEqual(tc, limit, 0);
verifyEmpty(tc, why);
end

function [E, contract, ref, opts] = fixture()
EEG = nqc_synth(struct('seconds', 75, 'nPerCond', 20, ...
    'artifactTrials', 0, 'blinkRate', 0));
E = nqc_epoched(EEG);
contract = pipecompare.eval.Contract('conditions', {'target', {'11'}; 'standard', {'31'}}, ...
    'epoch', [-0.2 1], 'baseline', [-0.2 0], ...
    'components', {'P3', [0.3 0.5], {'Pz'}});
ref = pipecompare.eval.Measure.reference(E, contract);
ref.blink = {};
opts = pipecompare.eval.Rank.defaults();
opts.nBoot = 20;
end
