function tests = test_statistics
%TEST_STATISTICS Numerical validation of the evaluation statistics against
%   external criteria (definitions, simulation with known truth).
tests = functiontests(localfunctions);
end

function setupOnce(~)
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
end

% ----------------------------------------------------------------- SME
function testSmeAnalytic(tc)
S = [1; 2; 3; 4; NaN];
verifyEqual(tc, pipecompare.eval.Measure.sme(S), std([1 2 3 4]) / 2, 'AbsTol', 1e-12);
end

function testSmeMatchesEmpiricalSdOfTheMean(tc)
% SME must estimate the SD of the averaged score across replications.
rs = RandStream('mt19937ar', 'Seed', 11);
N = 40; sigma = 7; reps = 4000;
m = zeros(reps, 1); sme = zeros(reps, 1);
for r = 1:reps
    x = 3 + sigma * randn(rs, N, 1);
    m(r) = mean(x); sme(r) = pipecompare.eval.Measure.sme(x);
end
verifyEqual(tc, mean(sme), std(m), 'RelTol', 0.05);
verifyEqual(tc, mean(sme), sigma / sqrt(N), 'RelTol', 0.05);
end

function testBootstrappedSmeOfPeakLatencyMatchesReplications(tc)
% bSME of the peak latency of the average must estimate the SD of that
% latency across independent replications (Luck et al., 2021).
rs = RandStream('mt19937ar', 'Seed', 4);
fs = 250; t = (0.2:1/fs:0.6) * 1000; N = 40; reps = 300;
o = struct('kind', 'peakLatency', 'polarity', 'positive', 'times', t);
opts = pipecompare.eval.Measure.defaults(); opts.nBootPeak = 400;
lat = zeros(reps, 1); bsme = zeros(reps, 1);
for r = 1:reps
    X = 5 * exp(-0.5 * ((t - 400) / 60) .^ 2) + 8 * randn(rs, N, numel(t));
    [lat(r), bsme(r)] = pipecompare.eval.Measure.estimate(o, X, opts, r);
end
fprintf('peak latency: empirical SD %.2f ms, mean bSME %.2f ms\n', std(lat), mean(bsme));
verifyEqual(tc, mean(bsme), std(lat), 'RelTol', 0.20);
end

% ------------------------------------------------- comparison and ranking
function testTieCriterionErrorRateAndPower(tc)
% With equal true noise the better-looking candidate must rarely be
% declared distinguishable (post-selection false "worse" rate <= ~5% at
% the default alpha); with 1.4x noise it must usually be.
rs = RandStream('mt19937ar', 'Seed', 5);
N = 100; sims = 150; o = pipecompare.eval.Rank.defaults(); o.nBoot = 800;
ref = nqc_ref(N);
fa = 0; hit = 0;
for k = 1:sims
    common = randn(rs, N, 1);
    a = common + randn(rs, N, 1); b = common + randn(rs, N, 1);
    c = common + 1.71 * randn(rs, N, 1);       % total SD 1.4x that of a
    fa = fa + worse(a, b, ref, o);
    hit = hit + worse(a, c, ref, o);
end
fprintf('not-distinguished test: false-worse rate %.3f, power %.3f\n', fa / sims, hit / sims);
verifyLessThan(tc, fa / sims, 0.06);
verifyGreaterThan(tc, hit / sims, 0.8);
end

function testManyCandidatesKeepTheFamilywiseErrorBounded(tc)
% 24 candidates with the same true noise: the share of searches in which
% ANY of them is declared worse must stay near 5% (per-comparison
% intervals against the selected best did this in 87/100 searches).
rs = RandStream('mt19937ar', 'Seed', 12);
K = 24; sims = 80; o = pipecompare.eval.Rank.defaults(); o.nBoot = 600;
for N = [30 100]
    ref = nqc_ref(N);
    fwer = 0; frac = 0;
    for k = 1:sims
        common = randn(rs, N, 1);
        cs = arrayfun(@(j) nqc_cand(common + randn(rs, N, 1)), 1:K);
        R = pipecompare.eval.Rank.run(cs, ref, o);
        fwer = fwer + any(~R.table.notDistinguished); frac = frac + mean(~R.table.notDistinguished);
    end
    fprintf('%d candidates, %d trials: any false "worse" %.3f, mean share %.4f\n', K, N, fwer / sims, frac / sims);
    verifyLessThanOrEqual(tc, fwer / sims, 0.08);
end
end

function testPowerDoesNotCollapseWithManyCandidates(tc)
% A clearly noisier candidate (1.4x SD) among 23 equal ones must still
% usually be declared worse.
rs = RandStream('mt19937ar', 'Seed', 13);
N = 100; sims = 40; ref = nqc_ref(N); o = pipecompare.eval.Rank.defaults(); o.nBoot = 600;
hit = 0;
for k = 1:sims
    common = randn(rs, N, 1);
    cs = arrayfun(@(j) nqc_cand(common + randn(rs, N, 1)), 1:23);
    cs(24) = nqc_cand(common + 1.71 * randn(rs, N, 1));
    R = pipecompare.eval.Rank.run(cs, ref, o);
    hit = hit + ~R.table.notDistinguished(24);
end
fprintf('power with 24 candidates: %.3f\n', hit / sims);
verifyGreaterThan(tc, hit / sims, 0.8);
end

function testMissingSignalCheckIsNotAPass(tc)
% NaN in an applicable signal metric means the check failed or was not
% computed; it must reject the candidate. A metric the check declares not
% applicable is skipped. A candidate without any signal check is rejected.
N = 40; ref = nqc_ref(N); x = randn(N, 1);
o = pipecompare.eval.Rank.defaults(); o.nBoot = 100;
ok = struct('source', 'injection', 'amplitudeError', 0.01, 'latencyShiftMs', 0, 'artifactPct', 0, ...
    'waveformCorr', 0.99, 'topoCorr', 0.99, 'chain', '', 'notApplicable', {{}});
nanAmp = ok; nanAmp.amplitudeError = NaN; nanAmp.waveformCorr = NaN;
naTopo = ok; naTopo.topoCorr = NaN; naTopo.notApplicable = {'topoCorr'};
cs = [nqc_cand(x, 'signal', ok) nqc_cand(x, 'signal', nanAmp) nqc_cand(x, 'signal', naTopo) nqc_cand(x)];
cs(4).signal = [];
R = pipecompare.eval.Rank.run(cs, ref, o);
verifyEqual(tc, R.table.status, {'feasible'; 'rejected'; 'feasible'; 'rejected'});
verifyTrue(tc, contains(R.table.reason{2}, 'amplitudeError not computed'));
verifyTrue(tc, contains(R.table.reason{2}, 'waveformCorr not computed'));
verifyTrue(tc, contains(R.table.reason{4}, 'signal check missing'));
% a NaN retention / interpolation count is not a pass either
c5 = nqc_cand(x); c5.interpolatedFraction = NaN;
R = pipecompare.eval.Rank.run(c5, ref, o);
verifyEqual(tc, R.table.status, {'rejected'});
end

function testPeakComparisonSpreadMatchesReplications(tc)
% For peak measures the paired resampling must reproduce how much the
% difference between two candidates' bSME varies across independent
% replications. Resampling with replacement around the inner bootstrap
% overstated it by ~50% (no power); half-samples rescaled to n do not.
rs = RandStream('mt19937ar', 'Seed', 80); fs = 250; t = (0.2:1/fs:0.6) * 1000; N = 60;
o = pipecompare.eval.Rank.defaults(); ref = nqc_ref(N, {'P3.lat'}, {'ms'});
sdBoot = zeros(1, 40); dpt = zeros(1, 40);
for s = 1:40
    base = 5 * exp(-0.5 * ((t - 400) / 60) .^ 2) + 4 * randn(rs, N, numel(t));
    c1 = nqc_cand(base + 6 * randn(rs, N, numel(t)), 'kind', {'peakLatency'}, 'times', t);
    c2 = nqc_cand(base + 9 * randn(rs, N, numel(t)), 'kind', {'peakLatency'}, 'times', t);
    B = pipecompare.eval.Rank.bootstrap([c1 c2], ref, o);
    sdBoot(s) = std(B{1}(2, :) - B{1}(1, :));
    dpt(s) = c2.m.composite - c1.m.composite;
end
fprintf('peak latency: SD of the difference across replications %.2f ms, paired resampling %.2f ms\n', std(dpt), mean(sdBoot));
verifyEqual(tc, mean(sdBoot), std(dpt), 'RelTol', 0.25);
end

function testIdenticalCandidatesAreNotDistinguished(tc)
% Two pipelines that give identical data (e.g. removing vs interpolating a
% channel outside the ROI) are both kept; a clearly noisier one is not.
rng(7); N = 50; ref = nqc_ref(N); x = randn(N, 1);
o = pipecompare.eval.Rank.defaults(); o.nBoot = 300;
R = pipecompare.eval.Rank.run([nqc_cand(x) nqc_cand(x) nqc_cand(2 * x)], ref, o);
verifyEqual(tc, R.table.notDistinguished, [true; true; false]);
end

function testEquivalenceNeedsTheWholeIntervalWithinTheMargin(tc)
% Equivalent = the simultaneous interval of the difference from the best
% lies within 5% of the best; not distinguished alone is absence of
% evidence. Identical or 2% noisier data are equivalent; much noisier data
% are neither; two independent equally noisy pipelines on few trials are
% not distinguished but not shown equivalent either.
rng(9); N = 50; ref = nqc_ref(N); x = randn(N, 1);
o = pipecompare.eval.Rank.defaults(); o.nBoot = 300;
R = pipecompare.eval.Rank.run([nqc_cand(x) nqc_cand(x) nqc_cand(1.02 * x) nqc_cand(2 * x)], ref, o);
verifyEqual(tc, R.table.equivalent, [true; true; true; false]);
verifyEqual(tc, R.table.notDistinguished, [true; true; true; false]);
rs = RandStream('mt19937ar', 'Seed', 21); N = 20; ref = nqc_ref(N);
R = pipecompare.eval.Rank.run([nqc_cand(randn(rs, N, 1)) nqc_cand(randn(rs, N, 1))], ref, o);
verifyTrue(tc, all(R.table.notDistinguished));
verifyEqual(tc, sum(R.table.equivalent), 1);                   % only the best itself
end

function testConfidenceSetEquivalentIsWithinTheMargin(tc)
% The automatic thresholds choose among the options shown equivalent to
% the best, not among all those not distinguished from it: the
% equivalent ones are a subset, always hold the best, and none is more
% than the margin worse at the point estimate. Options: the same scores
% scaled (1% and 2% noisier, three times noisier), so the paired
% bootstrap replicates are the same SME scaled.
rs = RandStream('mt19937ar', 'Seed', 41); N = 40; B = 999;
x = randn(rs, N, 1); W = pipecompare.eval.Rank.resampleWeights(rs, N, B, false);
n = sum(W, 1); m1 = sum(x .* W, 1) ./ n;
smeB = sqrt((sum(x .^ 2 .* W, 1) - n .* m1 .^ 2) ./ (n - 1) ./ n);
s = [1 1.01 1.02 3]';
pts = s * (std(x) / sqrt(N)); Bq = s * smeB;
[keep, equiv] = pipecompare.eval.Rank.confidenceSet(pts, Bq, 1, 0.02, 0.05);
verifyEqual(tc, keep(:)', [true true true false]);
verifyEqual(tc, equiv(:)', [true true true false]);
verifyTrue(tc, all(pts(equiv) <= 1.05 * pts(1)));
[~, e0] = pipecompare.eval.Rank.confidenceSet(pts, Bq, 1, 0.02);   % no margin: the best only
verifyEqual(tc, e0(:)', [true false false false]);
end

function testWeightedMeanAndItsSme(tc)
% Equal weights give the plain mean and SD / sqrt(n). With trials of
% unequal noise, inverse-variance weights give a mean with a much smaller
% error, still centred on the true value, and the weighted SME (sandwich)
% matches its spread over replications; the ranking's bootstrap takes the
% weights (each trial drawn once: the weighted SME).
rs = RandStream('mt19937ar', 'Seed', 51);
opts = pipecompare.eval.Measure.defaults();
o = struct('kind', 'scalar', 'polarity', 'positive', 'times', []);
x = randn(rs, 30, 1);
[e1, s1] = pipecompare.eval.Measure.estimate(o, x, opts, 0, ones(30, 1));
verifyEqual(tc, [e1 s1], [mean(x) std(x) / sqrt(30)], 'RelTol', 1e-12);
N = 30; sd = [ones(20, 1); 4 * ones(10, 1)]; w = 1 ./ sd .^ 2;
R = 4000; ep = zeros(R, 1); ew = ep; sw = ep;
for r = 1:R
    y = 2 + sd .* randn(rs, N, 1);
    ep(r) = mean(y);
    [ew(r), sw(r)] = pipecompare.eval.Measure.estimate(o, y, opts, 0, w);
end
fprintf('error of the mean: equal weights %.3f, inverse-variance %.3f (SME %.3f)\n', std(ep), std(ew), mean(sw));
verifyLessThan(tc, std(ew), 0.6 * std(ep));                     % 0.22 against 0.45 in theory
verifyEqual(tc, mean(ew), 2, 'AbsTol', 0.02);
verifyEqual(tc, mean(sw), std(ew), 'RelTol', 0.08);
o.X = {y}; o.w = {w};
verifyEqual(tc, pipecompare.eval.Measure.smeBoot(o, {ones(N, 1)}, opts, 0), sw(end), 'RelTol', 1e-10);
% peak measures: equal weights give the unweighted peak and bSME
p = struct('kind', 'peakAmplitude', 'polarity', 'positive', 'times', 1:20);
Y = randn(rs, 25, 20) + [zeros(1, 9) 3 zeros(1, 10)];
[a0, b0] = pipecompare.eval.Measure.estimate(p, Y, opts, 5);
[a1, b1] = pipecompare.eval.Measure.estimate(p, Y, opts, 5, ones(25, 1));
verifyEqual(tc, [a1 b1], [a0 b0], 'RelTol', 1e-10);
end

function testWhatBlinksAddToTheMean(tc)
% Measure.eyeBias: the share of the blink trials times the difference of
% their mean from the others' (what they add to the mean), its square
% estimated without bias (minus its variance), 0 with no blink trial kept
% or no other trial, and every mean (a column of weights) at once.
x = [1 1 1 3 3 NaN]'; b = logical([0 0 0 1 1 1])';
[bias, sq] = pipecompare.eval.Measure.eyeBias(x, b, ones(6, 1));
verifyEqual(tc, bias, 0.4 * 2, 'AbsTol', 1e-12);
v = var([1 1 1 3 3]);
verifyEqual(tc, sq, 0.8 ^ 2 - 0.4 ^ 2 * v * (1 / 2 + 1 / 3), 'AbsTol', 1e-12);
K = [1 1 1 0 0 0; 0 0 0 1 1 0; 1 1 1 1 1 1]';
[bias, sq] = pipecompare.eval.Measure.eyeBias(x, b, K);
verifyEqual(tc, bias, [0 0 0.8], 'AbsTol', 1e-12);
verifyEqual(tc, sq(1:2), [0 0]);
% unbiased: over many samples without a real bias, the mean of sq is about 0
rng(3); q = zeros(1, 2000);
for k = 1:2000
    [~, q(k)] = pipecompare.eval.Measure.eyeBias(randn(40, 1), rand(40, 1) < 0.3, ones(40, 1));
end
verifyLessThan(tc, mean(q), 0.01);
end

function testSelectionFrequencyOfTheRecommendation(tc)
% A pipeline clearly better than the others is recommended in nearly
% every resample of the trials. Among pipelines with the same true noise
% the one recommended on these trials is recommended in only about half
% of the resamples (simulated in Python: median 54%, 95th percentile 86%,
% with 5 such pipelines and 60 trials).
rs = RandStream('mt19937ar', 'Seed', 61); N = 60;
o = pipecompare.eval.Rank.defaults(); o.nBoot = 499; o.nSplits = 0;
x = randn(rs, N, 1);
R = pipecompare.eval.Rank.run([nqc_cand(x) nqc_cand(3 * x) nqc_cand(3.2 * x)], nqc_ref(N), o);
verifyEqual(tc, R.byStratum.recommended, 1);
verifyGreaterThan(tc, R.byStratum.stability, 0.95);
verifyEqual(tc, sum(R.table.chosenShare), 1, 'AbsTol', 1e-12);
f = zeros(1, 20);
for k = 1:20
    R = pipecompare.eval.Rank.run(arrayfun(@(j) nqc_cand(randn(rs, N, 1)), 1:5), nqc_ref(N), o);
    f(k) = R.byStratum.stability;
    verifyEqual(tc, R.table.chosenShare(R.byStratum.recommended), f(k));
end
fprintf('5 equally good pipelines: the recommended one is recommended again in %.0f%% of resamples (mean)\n', 100 * mean(f));
verifyLessThan(tc, mean(f), 0.75);
end

function testCrossfitRemovesTheOptimismOfTheBest(tc)
% 24 pipelines with the same true noise: the best-looking one's SME
% underestimates the true SME (selection), the cross-fitted value of the
% pipeline chosen on other trials does not. With one pipeline there is
% nothing to select and no overstatement.
rs = RandStream('mt19937ar', 'Seed', 31);
K = 24; N = 60; sims = 30; truth = 1 / sqrt(N);
o = pipecompare.eval.Rank.defaults(); o.nBoot = 200; o.nSplits = 10;
ref = nqc_ref(N);
ap = zeros(1, sims); ho = zeros(1, sims); one = zeros(1, sims);
for k = 1:sims
    cs = arrayfun(@(j) nqc_cand(randn(rs, N, 1)), 1:K);
    R = pipecompare.eval.Rank.run(cs, ref, o);
    ap(k) = R.byStratum.crossfit.apparent; ho(k) = R.byStratum.crossfit.honest;
    R1 = pipecompare.eval.Rank.run(cs(1), ref, o);
    one(k) = R1.byStratum.crossfit.overstatement;
end
fprintf('true SME %.4f: best apparent %.4f, cross-fitted %.4f; one pipeline: overstatement %.3f\n', ...
    truth, mean(ap), mean(ho), mean(one));
verifyLessThan(tc, mean(ap), 0.95 * truth);                     % the selection flatters the best
verifyEqual(tc, mean(ho), truth, 'RelTol', 0.05);               % trials not used to choose it do not
verifyGreaterThan(tc, mean(ho ./ ap - 1), 0.05);
verifyLessThan(tc, abs(mean(one)), 0.04);
end

function testCrossfitHoldsOutTheLimitChosenFromTheData(tc)
% One pipeline whose rejection limit was chosen from the data (uv 'auto')
% on rejection statistics unrelated to the scores: the limit's advantage
% is pure chance. Without its stored scores the selection check holds the
% kept trials fixed and finds no overstatement; with them (nested) it
% chooses the limit again on each half and finds the chance advantage
% (simulated in Python: about 3% against 0%).
rs = RandStream('mt19937ar', 'Seed', 41);
N = 60; sims = 30; ref = nqc_ref(N);
o = pipecompare.eval.Rank.defaults(); o.nBoot = 200; o.nSplits = 10;
ovN = zeros(1, sims); ovP = zeros(1, sims);
for k = 1:sims
    x = randn(rs, N, 1); M = rand(rs, N, 1);
    v = unique(M)'; obj = inf(size(v));
    for j = 1:numel(v)
        kp = M <= v(j);
        if sum(kp) >= 30, obj(j) = std(x(kp)) / sqrt(sum(kp)); end
    end
    u = v(find(obj <= (1 + o.equivalenceMargin) * min(obj), 1, 'last'));   % the most lenient within the margin
    xc = x; xc(M > u) = NaN;
    plain = nqc_cand(xc); plain.steps = {};
    nested = plain;
    C = struct('X', {{x - mean(x)}}, 'M', {{M}}, 'use', 1, 'v', [], 'gained', false);
    nested.steps = {struct('type', 'reject_threshold', 'cutoffRefined', false, ...
        'curves', struct('options', {{C}}, 'ok', true, 'chosen', 1))};
    cP = pipecompare.eval.Rank.crossfit(plain, ref, o, 'composite', plain.m.objectives.agg);
    cN = pipecompare.eval.Rank.crossfit(nested, ref, o, 'composite', nested.m.objectives.agg);
    verifyEqual(tc, [cP.nested cN.nested], [0 1]);
    ovP(k) = cP.overstatement; ovN(k) = cN.overstatement;
end
fprintf('limit chosen from the data: overstatement %.3f nested, %.3f with the kept trials fixed\n', mean(ovN), mean(ovP));
verifyGreaterThan(tc, mean(ovN) - mean(ovP), 0.015);
verifyLessThan(tc, abs(mean(ovP)), 0.02);
% a stored choice that does not reproduce the candidate's scores is not used
bad = nested; bad.steps{1}.curves.options{1}.X{1} = randn(rs, N, 1);
verifyEqual(tc, pipecompare.eval.Rank.crossfit(bad, ref, o, 'composite', 1).nested, 0);
end

function testCrossfitNeedsEnoughTrialsAndCanBeOff(tc)
x = randn(3, 1); ref = nqc_ref(3); o = pipecompare.eval.Rank.defaults(); o.nBoot = 100; o.minTrials = 2;
R = pipecompare.eval.Rank.run([nqc_cand(x) nqc_cand(2 * x)], ref, o);
verifyEmpty(tc, R.byStratum.crossfit);                          % fewer than 4 trials: not checked
o.nSplits = 0; x = randn(40, 1); ref = nqc_ref(40);
R = pipecompare.eval.Rank.run([nqc_cand(x) nqc_cand(2 * x)], ref, o);
verifyEmpty(tc, R.byStratum.crossfit);
end

function testRankingIndependentOfOtherCandidates(tc)
% Adding a candidate must not change how two others compare (0.6 failed
% this, v06_reproductions R4).
rng(2); N = 80; ref = nqc_ref(N);
base = randn(N, 1);
A = nqc_cand(base * 1.0); B = nqc_cand(base * 1.2);
o = pipecompare.eval.Rank.defaults(); o.nBoot = 500;
r1 = pipecompare.eval.Rank.run([A B], ref, o);
r2 = pipecompare.eval.Rank.run([A B nqc_cand(base * 3)], ref, o);
r3 = pipecompare.eval.Rank.run([A B nqc_cand(base * 1.1 + 0.3 * randn(N, 1))], ref, o);
verifyEqual(tc, r1.table.objective(1:2), r2.table.objective(1:2));
verifyEqual(tc, r1.table.objective(1:2), r3.table.objective(1:2));
verifyEqual(tc, r1.recommended, r2.recommended);
end

function testConstraintBoundaries(tc)
N = 20; ref = nqc_ref(N);
x = randn(N, 1);
atLimit = x; atLimit(11:end) = NaN;          % 10 of 20 trials = 50%
below = x; below(10:end) = NaN;              % 9 of 20 = 45%
o = pipecompare.eval.Rank.defaults(); o.nBoot = 100;
R = pipecompare.eval.Rank.run([nqc_cand(atLimit) nqc_cand(below)], ref, o);
verifyEqual(tc, R.table.status, {'feasible'; 'rejected'});
verifyTrue(tc, contains(R.table.reason{2}, 'retention'));
end

function testNoFeasibleCandidateGivesNoRecommendation(tc)
N = 30; ref = nqc_ref(N);
o = pipecompare.eval.Rank.defaults(); o.minTrials = 1000; o.nBoot = 50;
R = pipecompare.eval.Rank.run([nqc_cand(randn(N, 1)) nqc_cand(randn(N, 1))], ref, o);
verifyEmpty(tc, R.recommended);
verifyEmpty(tc, R.byStratum);
txt = evalc('pipecompare.eval.Rank.print(R, {''a'',''b''})');
verifyTrue(tc, contains(txt, 'NO FEASIBLE PIPELINE'));
end

function testStrataAreRankedSeparately(tc)
rng(3); N = 60; ref = nqc_ref(N);
o = pipecompare.eval.Rank.defaults(); o.nBoot = 200;
c1 = nqc_cand(randn(N, 1), 'stratum', 'ref=average');
c2 = nqc_cand(3 * randn(N, 1), 'stratum', 'ref=mastoids');   % much noisier, other measure
R = pipecompare.eval.Rank.run([c1 c2], ref, o);
verifyEqual(tc, numel(R.byStratum), 2);
verifyEmpty(tc, R.recommended);                  % no cross-stratum winner
verifyEqual(tc, sort([R.byStratum.recommended]), [1 2]);
end

function testObjectiveAndUnitSafety(tc)
% Measures in different units cannot be combined; the user chooses the
% one to optimize, the others are reported (sme_ columns).
rng(5); N = 80;
ref = nqc_ref(N, {'obj1', 'obj2'}, {'uV', 'ms'});
A = nqc_cand({{randn(N, 1)}, {5 * randn(N, 1)}}, 'units', {'uV', 'ms'});
B = nqc_cand({{2 * randn(N, 1)}, {1 * randn(N, 1)}}, 'units', {'uV', 'ms'});
o = pipecompare.eval.Rank.defaults(); o.nBoot = 100;
verifyError(tc, @() pipecompare.eval.Rank.run([A B], ref, o), 'PipeCompare:Objective');     % composite of uV and ms
o.objective = 'obj2';
R = pipecompare.eval.Rank.run([A B], ref, o);
verifyEqual(tc, R.best, 2);
verifyEqual(tc, R.objective, 'obj2');
verifyTrue(tc, all(isfinite(R.table.sme_obj1)));                                    % still reported
o.objective = {'obj1', 'obj2'};
verifyError(tc, @() pipecompare.eval.Rank.run([A B], ref, o), 'PipeCompare:Objective');     % one objective, not a list
end

function tf = worse(a, b, ref, o)
R = pipecompare.eval.Rank.run([nqc_cand(a) nqc_cand(b)], ref, o);
tf = sum(R.table.notDistinguished) < 2;
end

% ----------------------------------------------- gain-corrected objective
function testScalingTheDataDoesNotChangeTheObjective(tc)
% score = g*a + e: a candidate that scales signal and noise by 0.92 has a
% raw SME 8% lower but measures the signal no more precisely. SME/g is
% the same for both, so neither is preferred for it, and the least
% aggressive one (no amplitude change) is recommended.
rs = RandStream('mt19937ar', 'Seed', 11); N = 60; ref = nqc_ref(N);
x = 3 + randn(rs, N, 1);
sig = @(g) struct('source', 'test', 'amplitudeError', abs(g - 1), 'latencyShiftMs', 0, 'artifactPct', 0, ...
    'waveformCorr', 1, 'topoCorr', NaN, 'chain', '', 'notApplicable', {{'topoCorr'}}, 'gain', g);
a = nqc_cand(x, 'signal', sig(1)); b = nqc_cand(0.92 * x, 'signal', sig(0.92));
verifyLessThan(tc, b.m.objectives.agg, a.m.objectives.agg);          % raw SME prefers b
o = pipecompare.eval.Rank.defaults(); o.nBoot = 199;
R = pipecompare.eval.Rank.run([a b], ref, o);
verifyEqual(tc, R.table.objective(2), R.table.objective(1), 'RelTol', 1e-12);
verifyTrue(tc, all(R.table.notDistinguished));
verifyEqual(tc, R.recommended, 1);
verifyEqual(tc, R.table.gain_P3_mean(2), 0.92);
end

function testNonPositiveGainLeavesTheObjectiveUndefined(tc)
x = randn(RandStream('mt19937ar', 'Seed', 12), 40, 1);
sig = struct('source', 'test', 'amplitudeError', 0, 'latencyShiftMs', 0, 'artifactPct', 0, ...
    'waveformCorr', 1, 'topoCorr', NaN, 'chain', '', 'notApplicable', {{'topoCorr'}}, 'gain', -0.5);
o = pipecompare.eval.Rank.defaults(); o.nBoot = 99;
R = pipecompare.eval.Rank.run(nqc_cand(x, 'signal', sig), nqc_ref(40), o);
verifyEqual(tc, R.table.status{1}, 'rejected');
verifyTrue(tc, contains(R.table.reason{1}, 'objective undefined'));
end

function testBootstrapSizesGiveAnExactQuantile(tc)
% (B+1)*alpha integer: the critical value is an order statistic of the
% bootstrap maxima, with no interpolation (Davison & Hinkley, 1997).
o = pipecompare.eval.Rank.defaults();
verifyEqual(tc, mod((o.nBoot + 1) * o.alpha, 1), 0, 'AbsTol', 1e-9);
verifyEqual(tc, mod((o.nBootPeakOuter + 1) * o.alpha, 1), 0, 'AbsTol', 1e-9);
end

function testReasonSummaryKeepsEachReasonWhole(tc)
% A reason that contains ';' (a step key with a list) or ends with a full
% stop is counted as one reason, not cut apart.
x = randn(RandStream('mt19937ar', 'Seed', 3), 30, 1);
a = nqc_cand(x); a.status = 'failed';
a.message = 'reref(exclude={EOG1;EOG2}): channel(s) not in the dataset: EOG1, EOG2 (Edit > Channel locations).';
o = pipecompare.eval.Rank.defaults(); o.nBoot = 99;
R = pipecompare.eval.Rank.run([a a], nqc_ref(30), o);
txt = evalc('pipecompare.eval.Rank.print(R, {''a'', ''b''})');
verifyTrue(tc, contains(txt, [a.message ' (x2)']));
end

% ----------------------------------------------- dependent segments
function testBlockBootstrapKeepsTheErrorRateWithDependentSegments(tc)
% Segment scores of one recording are autocorrelated. Two candidates with
% the same true noise (AR(1), phi = 0.6) must rarely be declared different;
% the moving-block bootstrap keeps the false "worse" rate near the ERP one.
rs = RandStream('mt19937ar', 'Seed', 21);
N = 120; sims = 150; phi = 0.6;
o = pipecompare.eval.Rank.defaults(); o.nBoot = 599;
ar = @(n) filter(1, [1 -phi], randn(rs, n, 1)) * sqrt(1 - phi ^ 2);
fa = [0 0];
for k = 1:sims
    common = ar(N);
    a = nqc_cand(common + ar(N)); b = nqc_cand(common + ar(N));
    for mode = 1:2
        ref = nqc_ref(N); ref.segmented = (mode == 1);
        R = pipecompare.eval.Rank.run([a b], ref, o);
        fa(mode) = fa(mode) + any(~R.table.notDistinguished);
    end
end
fprintf('AR(1) segments: false "worse" rate block %.3f, i.i.d. %.3f\n', fa / sims);
verifyLessThanOrEqual(tc, fa(1) / sims, 0.07);
end

function testSameScoresAreRecognised(tc)
% Pipelines that pass the checks with exactly the same score: the choices
% compared make no difference (e.g. no epoch over any rejection threshold).
R.table = table({'feasible'; 'feasible'; 'rejected'}, [0.44437; 0.44437; 0.3], 'VariableNames', {'status', 'objective'});
verifyTrue(tc, pipecompare.eval.Rank.sameScores(R));
R.table.objective(2) = 0.45;
verifyFalse(tc, pipecompare.eval.Rank.sameScores(R));
R.table.status(2) = {'rejected'};
verifyFalse(tc, pipecompare.eval.Rank.sameScores(R));   % one pipeline: nothing compared
end
