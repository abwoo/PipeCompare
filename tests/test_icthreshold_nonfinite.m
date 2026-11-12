function tests = test_icthreshold_nonfinite
% Regression: all retained target trials are blink-flagged, so eye bias
% is not identifiable. The threshold search must return an unavailable
% objective, not an empty argmin assigned into f(k)/u(k).
tests = functiontests(localfunctions);
end

function setupOnce(~)
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
end

function testAllBlinkTargetsFailClosedWithoutAnException(tc)
[E, c, ref, ctx] = fixture();
ref.blink = {{true(ref.n(1), 1), false(ref.n(2), 1)}};
ctx.ref = ref;
q = single(0.2 * ones(1, E.nbchan));  % no artifact IC; one candidate only
[t, info] = pipecompare.eval.IcThreshold.choose(E, q, ctx);
verifyTrue(tc, isfinite(t));
verifyEqual(tc, info.sets, 1);
verifyEqual(tc, info.feasible, 0);
verifyTrue(tc, isnan(info.objectives(1)));
verifyTrue(tc, contains(info.why{1}, 'no finite objective'));
verifyEqual(tc, info.removedChosen, 0);
end

function testFiniteComparisonStillWorks(tc)
[E, ~, ref, ctx] = fixture();
ref.blink = {{false(ref.n(1), 1), false(ref.n(2), 1)}};
ctx.ref = ref;
q = single(0.2 * ones(1, E.nbchan));
[t, info] = pipecompare.eval.IcThreshold.choose(E, q, ctx);
verifyTrue(tc, isfinite(t));
verifyEqual(tc, info.sets, 1);
verifyEqual(tc, info.feasible, 1);
verifyTrue(tc, isfinite(info.objectives(1)));
verifyEqual(tc, info.removedChosen, 0);
end

function [E, c, ref, ctx] = fixture()
E = nqc_synth(struct('seconds', 120, 'nPerCond', 30, ...
    'artifactTrials', 0, 'blinkRate', 0));
c = pipecompare.eval.Contract('conditions', {'target', {'11'}; 'standard', {'31'}}, ...
    'epoch', [-0.2 1], 'baseline', [-0.2 0], ...
    'components', {'P3', [0.3 0.5], {'Pz'}});
n = E.nbchan;
E.icaweights = eye(n);
E.icasphere = eye(n);
E.icachansind = 1:n;
E.icawinv = [];
E.icaact = [];
E = eeg_checkset(E);
ref = pipecompare.eval.Measure.reference(E, c);
make = @(type) struct('type', type, 'params', struct(), ...
    'key', type, 'slot', type, 'label', type);
after = {make('epoch'), make('baseline')};
ctx = struct('contract', c, 'ref', ref, 'after', {after}, ...
    'rank', struct('nBoot', 40), 'highpass', 0);
end
