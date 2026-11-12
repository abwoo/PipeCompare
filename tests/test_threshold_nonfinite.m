function tests = test_threshold_nonfinite
% Regression for rejection threshold search when retention is sufficient
% but blink-bias cannot be estimated. A finite fallback is not a claim of
% a valid scientific objective or of an acceptable retained sample.
tests = functiontests(localfunctions);
end

function setupOnce(~)
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
end

function testNoFiniteObjectiveFailsClosed(tc)
[E, c, ref, m] = fixture();
ref.blink = {{true(ref.n(1), 1), false(ref.n(2), 1)}};
ctx = struct('contract', c, 'ref', ref, 'rank', struct('nBoot', 40));
[u, info] = pipecompare.eval.Threshold.choose(E, m, ctx);
verifyEqual(tc, u, 100, 'AbsTol', 1e-12);
verifyTrue(tc, isnan(info.objectiveBest));
verifyTrue(tc, isnan(info.objectiveChosen));
verifyTrue(tc, contains(info.note, 'no rejection limit has a finite objective'));
verifyEqual(tc, info.compared, 0);
verifyEqual(tc, info.rejectedChosen, 0);
end

function testFiniteObjectiveStillSelects(tc)
[E, c, ref, m] = fixture();
ref.blink = {{false(ref.n(1), 1), false(ref.n(2), 1)}};
ctx = struct('contract', c, 'ref', ref, 'rank', struct('nBoot', 40));
[u, info] = pipecompare.eval.Threshold.choose(E, m, ctx);
verifyTrue(tc, isfinite(u));
verifyTrue(tc, isfinite(info.objectiveBest));
verifyTrue(tc, isfinite(info.objectiveChosen));
verifyTrue(tc, isempty(info.note));
verifyGreaterThanOrEqual(tc, info.compared, 1);
end

function [E, c, ref, m] = fixture()
raw = nqc_synth(struct('seconds', 120, 'nPerCond', 30, ...
    'artifactTrials', 0, 'blinkRate', 0));
c = pipecompare.eval.Contract('conditions', {'target', {'11'}; 'standard', {'31'}}, ...
    'epoch', [-0.2 1], 'baseline', [-0.2 0], ...
    'components', {'P3', [0.3 0.5], {'Pz'}});
ref = pipecompare.eval.Measure.reference(raw, c);
[~, E] = evalc('pop_epoch(raw, c.allEvents(), c.epoch, ''epochinfo'', ''yes'')');
m = 100 * ones(1, E.trials);
end
