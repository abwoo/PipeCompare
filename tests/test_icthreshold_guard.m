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
