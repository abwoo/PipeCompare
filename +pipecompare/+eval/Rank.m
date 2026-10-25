classdef Rank
    %RANK Decide which candidates are acceptable and which is recommended.
    %
    %   R = pipecompare.eval.Rank.run(cands, ref, opts)
    %
    %   1. Execution failures are reported, never ranked.
    %   2. Constraints, each with a readable reason:
    %        trials per condition    >= minTrials
    %        retention per condition >= minRetention
    %        interpolated channels   <= maxInterpolated (fraction)
    %        signal preservation (pipecompare.eval.Injection): amplitude error
    %          <= maxAmplitudeError, peak shift <= maxLatencyShiftMs,
    %          artifactual deflection <= maxArtifactPct, waveform
    %          correlation >= minWaveformCorr, topography correlation >=
    %          minTopoCorr. A metric that was not computed fails.
    %   3. Objective (opts.objective): the gain-corrected SME, SME/g, where
    %      g is the factor by which the candidate scales a known signal in
    %      that measure (pipecompare.eval.Injection: window mean for mean
    %      amplitude, peak for peak amplitude, g = 1 for latency). With a
    %      pipeline acting linearly on signal + noise, score = g*a + e, so
    %      the score divided by g estimates the same quantity a for every
    %      candidate, with standard error SME/g; ranking SME/g is ranking
    %      signal-to-noise g*a/SME (Zhang, Garrett & Luck, 2024, rank filters
    %      by SNR = signal/SME, not by SME, because filters attenuate the
    %      signal too). It is invariant to any overall scaling of the data,
    %      which raw SME is not: a pipeline that shrinks signal and noise by
    %      9% passes a 10% amplitude-error limit and would lower raw SME by
    %      9% without being more precise. 'composite' = RMS over measures
    %      and conditions (when they share a unit), or the name of one
    %      measure (e.g. 'P3.peakLatency'), needed when the units differ.
    %      The other measures are reported, not ranked.
    %   3b. Total error: the objective adds to each condition's SME/g the
    %      bias b of its gain-corrected mean, sqrt((SME/g)^2 + b^2) (the
    %      root mean square error of the estimate; RMS over measures and
    %      conditions as above). b = (e + |mean| s) / g: e, what the trials
    %      with a blink add to the mean (Measure.eyeBias, from the eye
    %      channels of the data as they were; counted only beyond 2 of
    %      its standard errors, so 0 unless the data show it), and s, how much the
    %      pipeline's filters would change its gain if the signal had
    %      another shape (Injection.shapeSpreads: the largest change over
    %      a family of shapes, a minimax). The SME counts what varies from
    %      trial to trial; a blink in a share of the trials also moves the
    %      mean, and much more than it raises the SME. Columns noise
    %      (SME/g alone) and bias (b alone) of the table, eyeBias_ (e, in
    %      the measure's unit, the condition where it is largest) and
    %      shapeSpread_ (s) per measure (docs/METHODS.md, section 15).
    %   4. Uncertainty: a paired bootstrap over trials (the same resampled
    %      trials, by urevent, for every candidate). A candidate is "not
    %      distinguished from the best" unless the data show it is worse
    %      than some other candidate, with intervals simultaneous over all
    %      candidate pairs, conditional on the fitted candidates. A roughly
    %      5% family-wise rate was seen in the tested simulations (alpha
    %      = 0.02 chosen for 2-24 candidates, 20-100 trials,
    %      unequal condition sizes, correlated candidates; see
    %      test_statistics). This is the bootstrap max-statistic of White
    %      (2000) and Romano & Wolf (2005); the set kept is a model
    %      confidence set in the sense of Hansen, Lunde & Nason (2011). The
    %      critical value is the ((B+1)(1-alpha))-th ordered bootstrap
    %      maximum, with B chosen so that (B+1)*alpha is an integer
    %      (Davison & Hinkley, 1997): B = 1999 (40 values beyond it), and
    %      999 for peak measures whose nested bootstrap costs more (20).
    %      Not distinguished = absence of evidence, not equivalence.
    %      Equivalent: the upper bound of the difference from the best lies
    %      within opts.equivalenceMargin (default 5%) of the best objective,
    %      so the data show the difference is negligible (a one-sided test
    %      of equivalence, Schuirmann, 1987; the best is the lowest, so only
    %      the upper side is open). The bound is simultaneous over the
    %      candidates not shown worse (step-down, Romano & Wolf, 2005).
    %      Not distinguished but not equivalent means too few trials to tell.
    %   4b. Selection bias: the best of many objectives looks better than
    %      it is (part of its advantage is chance). crossfit estimates the
    %      objective of the pipeline chosen this way on trials not used to
    %      choose it: the trials of each condition are split in halves,
    %      the pipeline best on one half is scored on the other (rescaled
    %      to the full number of trials), and back; opts.nSplits random
    %      splits (consecutive halves for segments of one recording). A
    %      choice a pipeline makes from the measure itself (a rejection
    %      limit, an ICLabel threshold or components chosen from the data)
    %      is made again on the selecting half too (nested cross-fitting),
    %      so its advantage by chance is held out as well.
    %   5. The objectives are first shrunk through the pipelines'
    %      settings (factorialModel, shrink: each step is a factor and its
    %      settings the levels; the main effects are estimated from every
    %      pipeline that shares a setting, and what is left, interactions
    %      and noise, is shrunk by the positive-part James-Stein factor
    %      from the bootstrap covariance: factorial empirical Bayes). With
    %      30-40 trials nearly every pipeline is "not distinguished from
    %      the best", and the lowest raw objective is mostly chance; in
    %      simulation (18 pipelines, high-pass x low-pass x weights) the
    %      lowest shrunk objective had up to 60% less excess error than the
    %      lowest raw one. Among candidates the data do not distinguish
    %      from the best, a biased one (bias above BiasShare of its noise)
    %      only when no unbiased one is among them; of the others, the one
    %      with the lowest shrunk objective, unless candidates shown
    %      equivalent to the best (step 4) keep more trials: then the least
    %      aggressive of those (highest minimum retention, then smallest
    %      signal distortion, to AmpStep: two pipelines keeping the same
    %      trials, as with and without trial weights, differ in it only by
    %      filter tails of neighbouring epochs). Keeping more trials is a
    %      reason to choose only between pipelines the data show to be
    %      equally good: among those merely not told apart it cost 5-11% of
    %      true error in simulation.
    %   5b. Stability: how often the recommendation would be the same with
    %      other trials. In each paired bootstrap resample of the trials
    %      (those of step 4), the choice of steps 4 and 5 is made again on
    %      the resampled objectives (with the critical value of step 4:
    %      a candidate is not distinguished from the best when its
    %      objective is within it of the lowest); the share of resamples in
    %      which each candidate is recommended is its selection frequency
    %      (the bootstrap selection frequency of model selection; Efron,
    %      2014, JASA). A low frequency of the recommended pipeline says the
    %      data do not single it out.
    %   6. Candidates are compared only within a stratum: candidates that
    %      differ in a parameter that changes the measured quantity (e.g.
    %      the reference) are ranked separately and never against each
    %      other.
    %
    %   No condition difference or p-value is optimized. Data-dependent
    %   preprocessing still needs independent downstream validation.

    properties (Constant)
        BiasShare = 0.5   % a bias above this share of the noise makes a pipeline "biased" (step 5)
        AmpStep = 0.005   % signal distortions closer than this (0.5% of the amplitude) count as the same (step 5)
    end

    methods (Static)
        function o = defaults()
            o = struct('minTrials', 10, 'minRetention', 0.5, 'maxInterpolated', 0.2, ...
                'maxAmplitudeError', 0.10, 'maxLatencyShiftMs', 10, 'maxArtifactPct', 0.05, ...
                'minWaveformCorr', 0.95, 'minTopoCorr', 0.90, ...
                'objective', 'composite', 'nBoot', 1999, 'nBootPeakOuter', 999, 'nBootPeakInner', 100, ...
                'nBootPeak', 1000, 'alpha', 0.02, 'seed', 1, 'equivalenceMargin', 0.05, 'nSplits', 20);
        end

        function R = run(cands, ref, opts, levels)
            % levels (optional): per candidate, its settings as a cell of
            % {slot, key} rows (Executor), for the factorial shrinkage of
            % step 5; without them the objectives are not shrunk.
            if nargin < 3, opts = struct(); end
            if nargin < 4, levels = {}; end
            opts = pipecompare.utils.withDefaults(opts, pipecompare.eval.Rank.defaults());
            n = numel(cands);
            objName = resolveObjective(opts.objective, ref);
            nO = numel(ref.objectives);
            status = repmat({''}, n, 1); reasons = repmat({''}, n, 1); whyList = repmat({{}}, n, 1);
            notes = repmat({''}, n, 1); gains = nan(n, nO);
            stratum = repmat({''}, n, 1);
            primary = nan(n, 1); minRet = nan(n, 1); minKept = nan(n, 1); interp = nan(n, 1);
            ampErr = nan(n, 1); latSh = nan(n, 1); artPct = nan(n, 1); wCorr = nan(n, 1); tCorr = nan(n, 1);
            associationRisk=nan(n,nO);artifactSE=nan(n,nO);
            objAgg = nan(n, nO); noise = nan(n, 1); bias = nan(n, 1); eyeB = nan(n, nO); spreads = nan(n, nO);
            for i = 1:n
                c = cands(i);
                if isfield(c, 'stratum') && ~isempty(c.stratum), stratum{i} = c.stratum; end
                if isfield(c, 'notes') && ~isempty(c.notes), notes{i} = strjoin(c.notes, '; '); end
                if strcmp(c.status, 'excluded')      % excluded before it ran (e.g. its filters alone distort the signal)
                    status{i} = 'rejected'; whyList{i} = {c.message}; reasons{i} = c.message; continue;
                elseif strcmp(c.status, 'rejected')      % e.g. a rejection step removed every epoch
                    status{i} = 'rejected'; minRet(i) = 0; whyList{i} = {sprintf('retention 0%% (%s)', c.message)};
                    reasons{i} = whyList{i}{1}; continue;
                elseif ~strcmp(c.status, 'ok')
                    status{i} = 'failed'; reasons{i} = c.message; whyList{i} = {c.message}; continue;
                end
                m = c.m;
                g = gainOf(c, nO); gains(i, :) = g; sp = spreadOf(c, nO); spreads(i, :) = sp;
                objAgg(i, :) = [m.objectives.agg] ./ g;
                [primary(i), noise(i), bias(i)] = primaryPoint(m, objName, ref, g, sp);
                for k = 1:nO
                    [~,~,~,sensitivity]=pipecompare.eval.Rank.conditionErrors(m.objectives(k),g(k),sp(k));
                    associationRisk(i,k)=sqrt(mean(sensitivity.associationError.^2));
                    artifactSE(i,k)=max(sensitivity.artifactBiasSE);
                    eb = pipecompare.utils.fieldOr(m.objectives(k), 'eyeBias', 0);
                    if all(isfinite(eb))
                        [~, at] = max(abs(eb)); eyeB(i, k) = eb(at);   % the condition where blinks add most
                    end
                end
                minRet(i) = min(m.retention); minKept(i) = min(m.kept);
                interp(i) = c.interpolatedFraction;
                sg = c.signal;
                if ~isempty(sg)
                    ampErr(i) = sg.amplitudeError; latSh(i) = sg.latencyShiftMs; artPct(i) = sg.artifactPct;
                    wCorr(i) = sg.waveformCorr; tCorr(i) = sg.topoCorr;
                end
                why = {};
                assessedUse = 1:nO;
                if ~strcmp(objName,'composite'), assessedUse = find(strcmp(ref.objectives,objName)); end
                for k = assessedUse
                    assessed = pipecompare.utils.fieldOr(m.objectives(k),'eyeBiasAssessed',true);
                    if any(~assessed)
                        why{end+1} = 'artifact bias unassessed: no artifact-free comparison trials'; %#ok<AGROW>
                        break;
                    end
                end
                % every comparison is written so that a missing value (NaN) fails it:
                % a check that was not computed is never a pass
                if ~(minKept(i) >= opts.minTrials), why{end+1} = sprintf('only %d trials in a condition (< %d)', minKept(i), opts.minTrials); end %#ok<AGROW>
                if ~(minRet(i) >= opts.minRetention), why{end+1} = sprintf('retention %.0f%% (< %.0f%%)', 100*minRet(i), 100*opts.minRetention); end %#ok<AGROW>
                if ~(interp(i) <= opts.maxInterpolated), why{end+1} = sprintf('%.0f%% channels interpolated (> %.0f%%)', 100*interp(i), 100*opts.maxInterpolated); end %#ok<AGROW>
                if isempty(sg)
                    why{end+1} = 'signal check missing'; %#ok<AGROW>
                else
                    why = [why pipecompare.eval.Rank.signalReasons(sg, opts)]; %#ok<AGROW>
                end
                if ~isfinite(primary(i)), why{end+1} = 'objective undefined'; end %#ok<AGROW>
                if isempty(why), status{i} = 'feasible'; else, status{i} = 'rejected'; reasons{i} = strjoin(why, '; '); end
                whyList{i} = why;
            end

            lo = nan(n, 1); hi = nan(n, 1); hiKept = nan(n, 1); notDist = false(n, 1); equiv = false(n, 1);
            recommended = []; best = [];
            byStratum = struct('stratum', {}, 'best', {}, 'recommended', {}, 'set', {}, 'crossfit', {}, 'stability', {}, ...
                'factorial', {}, 'bootstrapCompleteShare', {}, 'bootstrapUnavailableShare', {});
            chosen = nan(n, 1); shrunk = nan(n, 1);
            strata = unique(stratum(strcmp(status, 'feasible')), 'stable');
            for s = 1:numel(strata)
                feas = find(strcmp(status, 'feasible') & strcmp(stratum, strata{s}));
                boot = pipecompare.eval.Rank.bootstrap(cands(feas), ref, opts, {objName});
                pts = primary(feas);
                [~, ib] = min(pts);
                [keep, l, h, hEq, crit] = bestSet(pts(:), boot{1}, ib, opts.alpha);
                lo(feas) = l; hi(feas) = h; notDist(feas(keep)) = true;
                % equivalent: the whole interval within the margin of the best
                hEq(~keep) = NaN; hiKept(feas) = hEq;
                equiv(feas) = keep & hEq <= opts.equivalenceMargin * pts(ib);
                t = feas(keep);
                % step 5: the objectives shrunk through the pipelines'
                % settings (factorial empirical Bayes), then the choice
                F = pipecompare.eval.Rank.factorialModel(levelsAt(levels, feas, cands(feas)), boot{1}, pts(:));
                th = pipecompare.eval.Rank.shrink(F, pts(:)); shrunk(feas) = th;
                biased = double(bias > pipecompare.eval.Rank.BiasShare * noise);
                rec = feas(choose(keep, th, biased(feas), equiv(feas), minRet(feas), ampStep(ampErr(feas))));
                cf = pipecompare.eval.Rank.crossfit(cands(feas), ref, opts, objName, pts(ib), F);
                [f, unavailableShare] = selectionFrequency(boot{1}, crit, minRet(feas), ampStep(ampErr(feas)), biased(feas), equiv(feas), F);
                completeShare=mean(all(isfinite(boot{1}),1));
                if completeShare<1
                    for q=feas(:)'
                        note=sprintf('paired bootstrap incomplete in %.1f%% of draws; incomplete maxima are unbounded; covariance fits require complete paired support',100*(1-completeShare));
                        pieces={notes{q},note};notes{q}=strjoin(pieces(~cellfun(@isempty,pieces)),'; ');
                    end
                end
                chosen(feas) = f;
                byStratum(end+1) = struct('stratum', strata{s}, 'best', feas(ib), 'recommended', rec, 'set', t, ...
                    'crossfit', cf, 'stability', f(feas == rec), 'factorial', F, ...
                    'bootstrapCompleteShare',completeShare,'bootstrapUnavailableShare',unavailableShare); %#ok<AGROW>
            end
            if isscalar(byStratum)
                recommended = byStratum.recommended; best = byStratum.best;
            end
            feasAll = find(strcmp(status, 'feasible'));
            order = [sortBy(feasAll, primary); find(strcmp(status, 'rejected')); find(strcmp(status, 'failed'))];
            T = table((1:n)', stratum, status, primary, noise, bias, shrunk, lo, hi, notDist, equiv, hiKept, chosen, minRet, minKept, interp, ...
                ampErr, latSh, artPct, wCorr, tCorr, reasons, notes, ...
                'VariableNames', {'id','stratum','status','objective','noise','bias','shrunk','diffLo','diffHi','notDistinguished','equivalent','diffHiKept', ...
                'chosenShare','minRetention','minTrials','interpolated','ampError', ...
                'latencyShiftMs','artifactPct','waveformCorr','topoCorr','reason','note'});
            for k = 1:nO
                T.(matlab.lang.makeValidName(['associationRisk_' ref.objectives{k}]))=associationRisk(:,k);
                T.(matlab.lang.makeValidName(['artifactBiasSE_' ref.objectives{k}]))=artifactSE(:,k);
                T.(matlab.lang.makeValidName(['sme_' ref.objectives{k}])) = objAgg(:, k);   % gain-corrected
                T.(matlab.lang.makeValidName(['gain_' ref.objectives{k}])) = gains(:, k);
                T.(matlab.lang.makeValidName(['eyeBias_' ref.objectives{k}])) = eyeB(:, k);   % in the measure's unit
                T.(matlab.lang.makeValidName(['shapeSpread_' ref.objectives{k}])) = spreads(:, k);
            end
            R = struct('table', T, 'whyList', {whyList}, 'order', order, 'best', best, 'recommended', recommended, ...
                'byStratum', byStratum, 'objective', objName, 'options', opts, 'units', {ref.units}, ...
                'validationScope','conditional-on-fitted-candidates','fullyHeldOut',false);
        end

        function why = signalReasons(sg, opts)
            % The signal-check limits a result sg (Injection.compare) breaks,
            % each as a readable reason. A metric that was not computed fails;
            % one listed in sg.notApplicable is skipped.
            opts = pipecompare.utils.withDefaults(opts, pipecompare.eval.Rank.defaults());
            why = {};
            na = {}; if isfield(sg, 'notApplicable'), na = cellstr(sg.notApplicable); end
            lim = {'amplitudeError', opts.maxAmplitudeError, 1, 'component amplitude changed by %.0f%%', 100; ...
                'latencyShiftMs', opts.maxLatencyShiftMs, 1, 'peak latency shifted by %.0f ms', 1; ...
                'artifactPct', opts.maxArtifactPct, 1, 'artifactual deflection of %.0f%%', 100; ...
                'waveformCorr', opts.minWaveformCorr, -1, 'recovered waveform r = %.2f', 1; ...
                'topoCorr', opts.minTopoCorr, -1, 'recovered topography r = %.2f', 1};
            for q = 1:size(lim, 1)
                if any(strcmp(lim{q, 1}, na)), continue; end
                v = sg.(lim{q, 1});
                if ~isfinite(v)
                    why{end+1} = sprintf('%s: %s not computed', sg.source, lim{q, 1}); %#ok<AGROW>
                elseif ~(lim{q, 3} * v <= lim{q, 3} * lim{q, 2})
                    why{end+1} = sprintf(['%s: ' lim{q, 4}], sg.source, lim{q, 5} * v); %#ok<AGROW>
                end
            end
        end

        function boot = bootstrap(cands, ref, opts, objNames)
            % Paired bootstrap. Returns one (candidates x B) matrix per
            % objective in objNames ('composite' allowed).
            if nargin < 4, objNames = {'composite'}; end
            opts = pipecompare.utils.withDefaults(opts, pipecompare.eval.Rank.defaults());
            s = RandStream('mt19937ar', 'Seed', opts.seed);
            peak = isfield(cands(1).m, 'objectives') && any(~strcmp({cands(1).m.objectives.kind}, 'scalar'));
            B = opts.nBoot; if peak, B = opts.nBootPeakOuter; end
            nC = numel(ref.ids); W = cell(1, nC); indices = cell(1, nC); scale = ones(1, nC);
            for c = 1:nC
                N = ref.n(c);
                if peak
                    % Peak measures already carry an inner bootstrap (bSME).
                    % Resampling WITH replacement around it overstates how
                    % much bSME varies (simulated: SD of a difference 11 ms
                    % vs 7.4 ms across replications), which left no power.
                    % Half-samples without replacement, rescaled to the full
                    % sample (bSME ~ 1/sqrt(n)), give the right spread
                    % (7.0 ms); the same half-samples for every candidate.
                    m = floor(N / 2);
                    W{c} = zeros(N, B);
                    for b = 1:B, W{c}(randperm(s, N, m), b) = 1; end
                    scale(c) = sqrt(m / N);
                else
                    [W{c}, indices{c}] = pipecompare.eval.Rank.resampleWeights(s, N, B, isfield(ref, 'segmented') && ref.segmented);
                end
            end
            per = perObjective(cands, ref, opts, W, scale, indices);
            boot = cellfun(@(q) combine(per, q, ref), objNames, 'UniformOutput', false);
        end

        function [W, idx] = resampleWeights(s, N, B, segmented)
            % N x B counts of each trial in B bootstrap resamples (random
            % stream s). Consecutive segments of one recording are dependent
            % (slow changes of state and noise), so they are resampled in
            % blocks of l = round(N^(1/3)) consecutive segments, the
            % moving-block bootstrap (Kunsch, 1989) with the block-length
            % order of Hall, Horowitz & Jing (1995); trials are in time order.
            if segmented && N >= 8
                l = max(2, round(N ^ (1 / 3))); nb = ceil(N / l);
                starts = randi(s, N - l + 1, nb, B);
                idx = zeros(nb * l, B);
                for j = 1:l, idx(j:l:end, :) = starts + (j - 1); end
                idx = idx(1:N, :);
            else
                idx = randi(s, N, N, B);
            end
            W = full(sparse(idx, repmat(1:B, N, 1), 1, N, B));
        end

        function [keep, equiv] = confidenceSet(pts, Bq, ib, alpha, margin)
            % Which of K options (point objectives pts, paired bootstrap
            % replicates Bq, K x B; ib the point-best) the data do not show
            % to be worse than another one: the test of step 4. equiv:
            % those also shown equivalent to the best (the whole interval of
            % the difference within margin times the best objective; the
            % best itself always).
            [keep, ~, ~, hEq] = bestSet(pts(:), Bq, ib, alpha);
            if nargin < 5, margin = 0; end
            equiv = keep & hEq <= margin * pts(ib);
            equiv(ib) = true;
        end

        function [chosen, nEq] = fixedSequence(pts, Bq, ib, seq, alpha, margin)
            % The most lenient of nested options shown equivalent to the
            % best, by fixed-sequence testing (Maurer, Hothorn & Lehmacher,
            % 1995; Westfall & Krishen, 2001). pts: point objectives (K x
            % 1); Bq: paired bootstrap replicates (K x B); ib: the best;
            % seq: the options more lenient than the best, nearest first.
            % Each in turn is compared with the best alone: equivalent when
            % the upper bound of the difference (the difference plus the
            % (1-alpha) bootstrap quantile of its deviations) is within
            % margin times the best objective; the first that is not ends
            % the sequence. Each test is made only when all before it
            % passed, so no correction for their number is needed, and the
            % options far from the best no longer widen the bound of those
            % near it, as the bound simultaneous over all pairs did (with it
            % the rule chose the best itself in every simulated run).
            % chosen: the last option that passed (ib when none); nEq:
            % how many are equivalent, the best included.
            chosen = ib; nEq = 1;
            for i = seq(:)'
                d = pts(i) - pts(ib);
                dev = (Bq(i, :) - Bq(ib, :)) - d; dev(~isfinite(dev))=Inf;
                if ~isfinite(d) || isempty(dev) || ~(d + upperQuantile(dev, alpha) <= margin * pts(ib)), break; end
                chosen = i; nEq = nEq + 1;
            end
        end

        function v = pointObjective(m, ref, opts, g, spread)
            % The ranking's objective of one candidate's measures m with
            % the signal gains g and shape spreads (Injection.shapeSpreads;
            % default 0), one per measure: what Rank.run ranks by.
            opts = pipecompare.utils.withDefaults(opts, pipecompare.eval.Rank.defaults());
            if nargin < 5 || isempty(spread), spread = zeros(1, numel(g)); end
            v = primaryPoint(m, resolveObjective(opts.objective, ref), ref, g, spread);
        end

        function name = objectiveName(obj, ref)
            % 'composite' or the one measure the objective option names
            name = resolveObjective(obj, ref);
        end

        function [total, noise, bias, sensitivity] = conditionErrors(o, gain, spread)
            % The same per-condition error used by primaryPoint. Benchmark
            % calibration must include both artifact and shape bias.
            noise = o.sme / gain;
            bias = biasOf(o, gain, spread);
            total = sqrt(noise .^ 2 + bias .^ 2);
            association=pipecompare.utils.fieldOr(o,'eyeBias',nan(size(o.estimate)));
            associationError=sqrt(noise.^2+((abs(association)+abs(o.estimate)*spread)/gain).^2);
            se=pipecompare.utils.fieldOr(o,'eyeBiasSE',nan(size(o.estimate)))/gain;
            sensitivity=struct('associationError',associationError,'artifactBiasSE',se, ...
                'scope','Association sensitivity, not causal bias or a calibrated absolute error bound');
        end

        function cf = crossfit(cands, ref, opts, objName, apparent, F)
            % The objective of the pipeline that is best on one half of the
            % trials, scored on the other half (rescaled to all trials by
            % sqrt(m/N), as the half-samples of bootstrap; SME ~ 1/sqrt(n)),
            % averaged over both directions and opts.nSplits random splits.
            % Conditional selection diagnostic: preprocessing and the
            % candidate family have already been fitted on all trials.
            % This is not an independent validation of the whole pipeline.
            % overstatement = honest / apparent - 1 (honest is a legacy
            % field name). [] when a
            % condition has fewer than 4 trials or nSplits is 0.
            % The choices a pipeline makes from the measure itself (a
            % rejection limit, an ICLabel threshold or components chosen by
            % the SME: the steps that store the scores of every option they
            % chose among, 'curves') are made again on the selecting half
            % (nestedHalf), so their advantage by chance is held out too
            % (nested cross-fitting); nested: how many candidates had one.
            % A candidate is not nested when the stored scores of its own
            % choice do not reproduce its scores (a step after the choice
            % changed them, or its cutoffs were fine-tuned afterwards).
            % The other data-driven decisions (bad channels, the ICA
            % decomposition, fine-tuned cutoffs, trial weights) were made
            % on all trials.
            % F (optional): the factorial model of step 5; the pipeline
            % chosen on a half is then the one with the lowest shrunk
            % objective there (the same shrinkage, fixed), as on all trials.
            cf = [];
            if nargin < 5, apparent = NaN; end
            if nargin < 6, F = []; end
            opts = pipecompare.utils.withDefaults(opts, pipecompare.eval.Rank.defaults());
            if opts.nSplits < 1 || any(ref.n < 4), return; end
            s = RandStream('mt19937ar', 'Seed', opts.seed + 101);
            nC = numel(ref.ids);
            segmented = isfield(ref, 'segmented') && ref.segmented;
            nS = opts.nSplits; if segmented, nS = 1; end
            honest = []; inSample = [];
            inner = arrayfun(@(c) innerOf(c, ref, objName), cands, 'UniformOutput', false);
            nest = find(~cellfun(@isempty, inner));
            for sp = 1:nS
                WA = cell(1, nC); WB = cell(1, nC); scA = ones(1, nC); scB = ones(1, nC);
                for c = 1:nC
                    N = ref.n(c); m = floor(N / 2);
                    % consecutive segments are dependent: earlier vs later half
                    if segmented, p = 1:N; else, p = randperm(s, N); end
                    WA{c} = zeros(N, 1); WA{c}(p(1:m)) = 1; scA(c) = sqrt(m / N);
                    WB{c} = zeros(N, 1); WB{c}(p(m+1:end)) = 1; scB(c) = sqrt((N - m) / N);
                end
                vA = combine(perObjective(cands, ref, opts, WA, scA), objName, ref);
                vB = combine(perObjective(cands, ref, opts, WB, scB), objName, ref);
                % vXY: chosen on half X, scored on half Y
                vAA = vA; vAB = vB; vBB = vB; vBA = vA;
                for k = nest(:)'
                    [vAA(k), vAB(k)] = nestedHalf(inner{k}, WA, WB, scA, scB, ref, opts);
                    [vBB(k), vBA(k)] = nestedHalf(inner{k}, WB, WA, scB, scA, ref, opts);
                end
                vAA(~isfinite(vAA)) = NaN; vBB(~isfinite(vBB)) = NaN;
                if all(isnan(vAA)) || all(isnan(vBB)), continue; end
                [~, ia] = min(pipecompare.eval.Rank.shrink(F, vAA(:))); [~, ib] = min(pipecompare.eval.Rank.shrink(F, vBB(:)));
                honest = [honest vAB(ia) vBA(ib)]; inSample = [inSample vAA(ia) vBB(ib)]; %#ok<AGROW>
            end
            ok = isfinite(honest);
            if ~any(ok), return; end
            h = mean(honest(ok));
            cf = struct('apparent', apparent, 'honest', h, 'overstatement', h / apparent - 1, ...
                'halfInSample', mean(inSample(ok)), 'nSplits', nS, 'nested', numel(nest), ...
                'scope','conditional-selection','fullyHeldOut',false,'totalCandidates',numel(cands), ...
                'fixedDecisions',{{'candidate eligibility','HMM and bad channels','ICA','trial weights', ...
                'fine-tuned cutoffs','signal checks','factorial shrinkage'}});
        end

        function F = factorialModel(levels, Bq, pts)
            % Step 5's model of the K pipelines' objectives through their
            % settings. levels: per pipeline its settings ({slot, key}
            % rows; {} = unknown, or a pipeline left out of the model, as
            % one whose cutoffs were fine-tuned afterwards); Bq: paired
            % bootstrap replicates (K x B); pts: the objectives (K x 1).
            % Each slot (step) is a factor whose levels are its settings
            % ('absent' when skipped); the design D has an intercept and a
            % column per level but the first of each factor (main effects).
            % With H = D D^+ and R = I - H, the objectives are y = H y + r:
            % the main effects, which every pipeline sharing a setting
            % estimates together, and r, the interactions plus noise. With
            % Sigma the bootstrap covariance of y (the noise of the
            % objectives, paired), E|r|^2 = tau^2 tr(R) + tr(R Sigma) when
            % the interactions are independent with variance tau^2, so r is
            % shrunk towards 0 by the positive-part James-Stein factor
            % s = max(0, 1 - (df - 2)/df tr(R Sigma) / |r|^2), df = tr(R)
            % (Efron & Morris, 1973; empirical Bayes): the interactions
            % the data show are kept, the rest is taken for noise. F: the
            % pipelines in the model (use), H, s, df; [] when fewer than 4
            % pipelines are in it, df < 3, or paired covariance is unassessed.
            F = [];
            K = numel(pts);
            if numel(levels) ~= K, return; end
            use = find(~cellfun(@isempty, levels(:)'));
            ok = all(isfinite(Bq(use, :)), 1);
            if numel(use) < 4 || size(Bq,2) < 20 || ~all(ok) || ~all(isfinite(pts(use))), return; end
            slots = {};
            for i = use, slots = [slots levels{i}(:, 1)']; end %#ok<AGROW>
            slots = unique(slots, 'stable');
            D = ones(numel(use), 1);
            for q = 1:numel(slots)
                v = repmat({'absent'}, numel(use), 1);
                for a = 1:numel(use)
                    L = levels{use(a)}; at = find(strcmp(L(:, 1), slots{q}), 1);
                    if ~isempty(at), v{a} = L{at, 2}; end
                end
                [u, ~, g] = unique(v, 'stable');
                for j = 2:numel(u), D(:, end+1) = double(g == j); end %#ok<AGROW>
            end
            H = D * pinv(D); Rm = eye(numel(use)) - H;
            df = trace(Rm);
            if df < 3 - 1e-9, return; end
            Sig = cov(Bq(use, ok)');
            r = Rm * pts(use);
            s = max(0, 1 - (df - 2) / df * trace(Rm * Sig) / max(r' * r, eps));
            F = struct('use', use, 'H', H, 's', s, 'df', round(df), 'factors', numel(slots));
        end

        function th = shrink(F, v)
            % The objectives v (K x 1) through the factorial model F
            % (factorialModel): H v + s (I - H) v for the pipelines in it,
            % v itself for the others, and all of v when F is [] or a value
            % in the model is not finite.
            th = v;
            if isempty(F), return; end
            x = v(F.use);
            if ~all(isfinite(x)), return; end
            fit = F.H * x;
            th(F.use) = fit + F.s * (x - fit);
        end

        function T = marginal(R, leaves, searched)
            % For each searched parameter: how each value fares across all
            % evaluated candidates (median objective, median minimum
            % retention) and among feasible ones (best objective).
            rows = {};
            st = R.table.status; comp = R.table.objective; ret = R.table.minRetention;
            for s = 1:numel(searched)
                slot = searched(s).slot; param = searched(s).param;
                vals = choiceValues(leaves, slot, param);
                u = unique(vals, 'stable');
                for v = 1:numel(u)
                    sel = strcmp(vals, u{v});
                    ev = sel & ~strcmp(st, 'failed');
                    f = sel & strcmp(st, 'feasible');
                    best = min(comp(f)); if isempty(best), best = NaN; end
                    rows(end+1, :) = {sprintf('%s.%s', slot, param), u{v}, sum(sel), sum(f), ...
                        median(comp(ev), 'omitnan'), median(ret(ev), 'omitnan'), best}; %#ok<AGROW>
                end
            end
            if isempty(rows), T = table(); return; end
            T = cell2table(rows, 'VariableNames', {'parameter','value','n','nFeasible', ...
                'medianObjective','medianMinRetention','bestFeasibleObjective'});
        end

        function T = robustness(cands, R, leaves, searched, ref)
            % Multiverse summary: how much each measure's point estimate
            % (per condition) varies over the feasible pipelines, and which
            % searched choice accounts for most of that variation. Every
            % allowed pipeline has been run, so this is the full multiverse
            % of the plan (Steegen et al., 2016; for ERPs, Clayson et al.,
            % 2021), not a sample. Within each stratum only: estimates under
            % different references measure different quantities.
            % share = eta^2 of the choice = sum_v n_v (mean_v - mean)^2 /
            % sum_i (x_i - mean)^2 over the feasible pipelines: the fraction
            % of the spread that differs between that choice's values. It
            % describes sensitivity; it is not used for the ranking.
            rows = {};
            st = R.table.status; strat = R.table.stratum;
            for su = unique(strat(strcmp(st, 'feasible')), 'stable')'
                feas = find(strcmp(st, 'feasible') & strcmp(strat, su{1}));
                if numel(feas) < 2, continue; end
                for o = 1:numel(ref.objectives)
                    for c = 1:numel(ref.names)
                        x = arrayfun(@(i) cands(i).m.objectives(o).estimate(c), feas);
                        ok = isfinite(x); x = x(ok); f = feas(ok);
                        if numel(x) < 2, continue; end
                        top = '-'; share = NaN; ss = sum((x - mean(x)) .^ 2);
                        for q = 1:numel(searched)
                            v = choiceValues(leaves(f), searched(q).slot, searched(q).param);
                            [u, ~, g] = unique(v);
                            if numel(u) < 2 || ss <= 0, continue; end
                            between = sum(accumarray(g, 1) .* (accumarray(g, x(:)) ./ accumarray(g, 1) - mean(x)) .^ 2);
                            if isnan(share) || between / ss > share
                                share = between / ss; top = sprintf('%s.%s', searched(q).slot, searched(q).param);
                            end
                        end
                        rows(end+1, :) = {su{1}, ref.objectives{o}, ref.names{c}, ref.units{o}, numel(x), ...
                            min(x), median(x), max(x), max(x) - min(x), top, share}; %#ok<AGROW>
                    end
                end
            end
            if isempty(rows), T = table(); return; end
            T = cell2table(rows, 'VariableNames', {'stratum','measure','condition','unit','nPipelines', ...
                'min','median','max','range','mostInfluentialChoice','share'});
        end

        function tf = sameScores(R)
            % true when the pipelines that passed the checks all have the
            % same score: the choices compared make no difference on these data
            v = R.table.objective(strcmp(R.table.status, 'feasible'));
            tf = numel(v) > 1 && max(v) - min(v) <= 1e-9 * max(abs(v));
        end

        function c = commonReasons(R, k)
            % The k most frequent reasons why candidates are not feasible, as
            % 'reason (xN)'. Each reason is kept as its own item (a step key
            % may hold ';'), and only free-standing numbers are masked, so
            % reasons that differ only in a value are counted together.
            parts = [R.whyList{~strcmp(R.table.status, 'feasible')}];
            parts = regexprep(parts, '(?<![\w.])\d+(\.\d+)?', '#');
            [u, ~, ic] = unique(parts(~cellfun(@isempty, parts)));
            c = {};
            if isempty(u), return; end
            [cnt, ord] = sort(accumarray(ic(:), 1), 'descend');
            c = arrayfun(@(j) sprintf('%s (x%d)', u{ord(j)}, cnt(j)), 1:min(k, numel(ord)), 'UniformOutput', false);
        end

        function print(R, labels)
            T = R.table; o = R.options;
            pipecompare.utils.log(['Ranking by %s estimated total error (noise: SME / signal gain, and bias: artifact association and signal-shape sensitivity; lower = smaller fitted-model error). diff = difference from the best, ', ...
                '%.0f%% conditional interval simultaneous over the fitted candidates; "nd" = not distinguished from the best by these data ', ...
                '(not equivalence); "eq" = equivalent: the whole interval within %.0f%% of the best.'], R.objective, ...
                100 * (1 - o.alpha), 100 * o.equivalenceMargin);
            fprintf('   %-4s %-9s %9s %7s %19s %3s %6s %6s %6s %5s  %s\n', 'id', 'status', 'objective', 'bias', 'diff vs best [CI]', 'nd', 'minRet', 'interp', 'ampErr', 'art', 'pipeline');
            for k = R.order(:)'
                mark = ' '; if any([R.byStratum.recommended] == k), mark = '*'; end
                ci = ''; if isfinite(T.diffLo(k)), ci = sprintf('[%+.3f %+.3f]', T.diffLo(k), T.diffHi(k)); end
                nd = ''; if T.equivalent(k), nd = 'eq'; elseif T.notDistinguished(k), nd = 'nd'; end
                fprintf('  %s%-4d %-9s %9.4g %7.3g %19s %3s %5.0f%% %5.0f%% %5.0f%% %4.1f%%  %s\n', mark, k, T.status{k}, ...
                    T.objective(k), T.bias(k), ci, nd, 100*T.minRetention(k), 100*T.interpolated(k), ...
                    100*T.ampError(k), 100*T.artifactPct(k), labels{k});
                if ~isempty(T.stratum{k}), fprintf('        stratum: %s\n', T.stratum{k}); end
                if ~isempty(T.reason{k}), fprintf('        -> %s\n', T.reason{k}); end
                if ~isempty(T.note{k}), fprintf('        note: %s\n', T.note{k}); end
            end
            if isempty(R.byStratum)
                pipecompare.utils.log('NO FEASIBLE PIPELINE: no candidate satisfies the constraints (none relaxed).');
                common = pipecompare.eval.Rank.commonReasons(R, 3);
                if ~isempty(common), pipecompare.utils.log('Most common reasons: %s', strjoin(common, '; ')); end
                return;
            end
            for s = 1:numel(R.byStratum)
                b = R.byStratum(s);
                lab = ''; if ~isempty(b.stratum), lab = sprintf(' [stratum %s]', b.stratum); end
                pipecompare.utils.log(['Recommended (*)%s: candidate %d (fewest trials lost, then least distortion, among the %d ', ...
                    'not distinguished from the best); best objective: candidate %d.'], lab, b.recommended, numel(b.set), b.best);
                if isfield(b,'bootstrapCompleteShare')&&isfinite(b.bootstrapCompleteShare)&&b.bootstrapCompleteShare<1
                    pipecompare.utils.log('Paired bootstrap: %.1f%% complete draws; %.1f%% had no assessed recommendation. Incomplete maxima contribute unbounded uncertainty.', ...
                        100*b.bootstrapCompleteShare,100*b.bootstrapUnavailableShare);
                end
                if isfinite(b.stability)
                    pipecompare.utils.log('Stability: in resamples of the trials, candidate %d is recommended in %.0f%% of them.', ...
                        b.recommended, 100 * b.stability);
                end
                if ~isempty(b.crossfit)
                    pipecompare.utils.log(['Selection check: the pipeline best on half of the trials has objective %.4g on the ', ...
                        'other half (rescaled to all trials), against %.4g for the best on all trials: the best looks %.0f%% ', ...
                        'better in this conditional diagnostic (%d split(s), both directions; the limits, thresholds or components chosen ', ...
                        'from the data made again on the selecting half for %d pipeline(s)).'], b.crossfit.honest, b.crossfit.apparent, ...
                        100 * b.crossfit.overstatement, b.crossfit.nSplits, pipecompare.utils.fieldOr(b.crossfit, 'nested', 0));
                    pipecompare.utils.log(['Conditional on preprocessing, candidate eligibility and shrinkage fitted on all trials; ', ...
                        'this is not held-out validation of the complete pipeline or of new recordings.']);
                end
            end
            if numel(R.byStratum) > 1
                pipecompare.utils.log(['Strata differ in what is measured (e.g. the reference); their results are not ', ...
                    'comparable, so there is no overall recommendation. Choose the stratum that fits your analysis ', ...
                    'and adopt its recommendation: adopt(result, id).']);
            end
        end
    end
end

% ---------------------------------------------------------------------
function name = resolveObjective(obj, ref)
% 'composite' (all measures share a unit) or the name of one measure.
if iscell(obj)
    assert(isscalar(obj), 'PipeCompare:Objective', 'Choose one objective (or ''composite''), not a list.');
    obj = obj{1};
end
name = char(obj);
if strcmp(name, 'composite')
    assert(numel(unique(ref.units)) == 1, 'PipeCompare:Objective', ...
        ['The measures have different units (%s); a composite would add incompatible quantities. ', ...
         'Choose the measure to optimize, e.g. ''%s''.'], strjoin(unique(ref.units), ', '), ref.objectives{1});
else
    assert(any(strcmp(name, ref.objectives)), 'PipeCompare:Objective', 'Unknown objective %s; available: composite, %s', ...
        name, strjoin(ref.objectives, ', '));
end
end

function [v, noise, bias] = primaryPoint(m, name, ref, g, spread)
% The total error: the RMS over measures and conditions of
% sqrt((SME/g)^2 + b^2), b the bias of the gain-corrected mean
% (biasOf); noise and bias: the same RMS of SME/g alone and of b alone.
use = 1:numel(g);
if ~strcmp(name, 'composite'), use = find(strcmp(ref.objectives, name)); end
s2 = []; b2 = [];
for k = use
    [~, ns, bs] = pipecompare.eval.Rank.conditionErrors(m.objectives(k), g(k), spread(k));
    s2 = [s2 ns .^ 2]; %#ok<AGROW>
    b2 = [b2 bs .^ 2]; %#ok<AGROW>
end
noise = sqrt(mean(s2)); bias = sqrt(mean(b2));
v = sqrt(mean(s2 + b2));
end

function b = biasOf(o, g, spread, eyeSq)
% The bias of measure o's gain-corrected mean in each condition (1 x
% conditions), in its unit: what blinks add to it (Measure.eyeBias:
% the part of it beyond 2 of its standard errors, squared; eyeSq
% replaces it, e.g. a bootstrap resample's) plus what the filters would
% change if the signal had another shape (|mean| x spread,
% Injection.shapeSpreads), both divided by the gain g.
if nargin < 4, eyeSq = pipecompare.utils.fieldOr(o, 'eyeBiasSq', zeros(size(o.estimate))); end
est = o.estimate; est(~isfinite(est)) = 0;
b = (sqrt(eyeSq) + abs(est) * spread) / g;
end

function s = spreadOf(c, nO)
% the shape spread of each measure (Injection.shapeSpreads), 0 when unknown
s = zeros(1, nO);
if isfield(c, 'signal') && isstruct(c.signal) && isfield(c.signal, 'shapeSpread') && numel(c.signal.shapeSpread) == nO
    s = double(c.signal.shapeSpread(:)');
    s(~isfinite(s)) = 0;
end
end

function g = gainOf(c, nO)
% Signal gain per measure from the signal check (1 when the candidate has
% none, e.g. synthetic records). A gain that is not positive leaves the
% objective undefined: the signal was lost or inverted.
g = ones(1, nO);
if isfield(c, 'signal') && isstruct(c.signal) && isfield(c.signal, 'gain') && numel(c.signal.gain) == nO
    g = double(c.signal.gain(:)');
end
g(~(g > 0)) = NaN;
end

function [keep, lo, hi, hiEq, c] = bestSet(pts, Bq, ib, alpha)
% Which candidates the data do not show to be worse than another one.
% pts: point objectives (K x 1); Bq: paired bootstrap replicates (K x B);
% ib: index of the point-best. lo/hi: interval of each difference from
% the best (0 for the best itself). c: the critical value (a candidate is
% kept when its objective is within c of the lowest).
K = numel(pts); keep = false(K, 1); lo = nan(K, 1); hi = nan(K, 1);
keep(ib) = true;
% Simultaneous over ALL ordered pairs (bootstrap max statistic): with
% probability ~>= 1-alpha no candidate is wrongly declared worse than any
% other, however many candidates are compared and whichever one looks
% best. This is what makes "not distinguished from the best" valid after
% selecting the best among many. The statistic is the raw difference, not
% a t-ratio: with few trials the per-pair bootstrap SD is itself noisy and
% the maximum of t-ratios over many pairs is driven by the pairs whose SD
% happens to be underestimated (simulation: 40% false "worse" at 30
% trials and 24 candidates). Raw differences keep the error rate bounded
% at the cost of power for pairs that differ little.
B = size(Bq, 2); maxD = -inf(1, B); complete=all(isfinite(Bq),1);
Dhat = pts - pts';                             % K x K: i minus j
for i = 1:K
    D = Bq(i, :) - Bq;                         % K x B: candidate i minus each j
    D(~isfinite(D)) = NaN;
    maxD = max(maxD, max(D - Dhat(i, :)', [], 1, 'omitnan'));
end
% Preserve all B draws: an undefined family maximum is unbounded.
% Dropping that draw, or only its unavailable pairs, truncates the tail.
maxD(~complete)=Inf;
if B==0, c=Inf; else, c=upperQuantile(maxD,alpha); end
worse = Dhat > c;
keep = keep | ~any(worse, 2);
lo = Dhat(:, ib) - c; hi = Dhat(:, ib) + c;
lo(ib) = 0; hi(ib) = 0;
% For equivalence, the upper bounds are simultaneous over the candidates
% kept only (step-down, Romano & Wolf, 2005): a candidate already shown
% worse would otherwise widen every interval, and no difference among the
% others could ever be shown negligible.
k = find(keep); maxK = -inf(1, B);
for i = k(:)'
    D = Bq(i, :) - Bq(k, :); D(~isfinite(D)) = NaN;
    maxK = max(maxK, max(D - Dhat(i, k)', [], 1, 'omitnan'));
end
maxK(~all(isfinite(Bq(k,:)),1))=Inf;
if B==0, cK=Inf; else, cK=max(upperQuantile(maxK,alpha),0); end
hiEq = Dhat(:, ib) + cK; hiEq(ib) = 0; hiEq(~keep) = Inf;
end

function a = ampStep(a)
% signal distortions to AmpStep (step 5)
a = round(a / pipecompare.eval.Rank.AmpStep) * pipecompare.eval.Rank.AmpStep;
end

function k = choose(keep, th, biased, eqv, ret, amp)
% Step 5 among K candidates: of those kept (not distinguished from the
% best), the unbiased ones (all when none is); among them the one with the
% lowest shrunk objective th, and with it those shown equivalent to the
% best (eqv); of these, the highest minimum retention, the smallest
% amplitude error (to AmpStep), the lowest shrunk objective.
U = find(keep(:) & ~biased(:));
if isempty(U), U = find(keep(:)); end
[~, j] = min(th(U));
E = unique([U(eqv(U)); U(j)]);
[~, o] = sortrows(round([-ret(E), amp(E), th(E)], 10));
k = E(o(1));
end

function [f, unavailable] = selectionFrequency(Bq, c, ret, amp, biased, eqv, F)
% The share of the bootstrap resamples (columns of Bq, K x B) in which
% each of K candidates is recommended: those within c of the lowest
% objective of the resample, then the rule of step 5 (choose) on the
% resample's objectives shrunk with the same model F; which ones are
% equivalent to the best (eqv) is that of all trials.
K = size(Bq, 1); cnt = zeros(K, 1); used = 0;
for b = 1:size(Bq, 2)
    v = Bq(:, b); ok = isfinite(v);
    if ~any(ok), continue; end
    th = pipecompare.eval.Rank.shrink(F, v);
    k = choose(ok & v - min(v(ok)) <= c, th, biased, eqv, ret, amp);
    cnt(k) = cnt(k) + 1; used = used + 1;
end
B=size(Bq,2);f=cnt/max(B,1);unavailable=(B-used)/max(B,1);
if used==0,f(:)=NaN;end
if B==0,unavailable=NaN;end
end

function v = upperQuantile(z, alpha)
% The (1-alpha) quantile of the bootstrap distribution as the order
% statistic z_((B+1)(1-alpha)) (Davison & Hinkley, 1997): exact when
% (B+1)*alpha is an integer, rounded up (conservative) otherwise.
z = sort(z(:));
r = min(numel(z), ceil((numel(z) + 1) * (1 - alpha) - 1e-9));
v = z(max(r, 1));
end

function L = levelsAt(levels, feas, cands)
% the settings of the candidates feas for factorialModel; {} for one whose
% cutoffs were fine-tuned after the comparison (a setting no other
% pipeline shares, so the model would leave it as it is anyway) and for
% all when none were given
L = repmat({{}}, numel(feas), 1);
if isempty(levels), return; end
for a = 1:numel(feas)
    st = pipecompare.utils.fieldOr(cands(a), 'steps', {});
    if any(cellfun(@(x) isfield(x, 'cutoffRefined') && isequal(x.cutoffRefined, true), st)), continue; end
    L{a} = levels{feas(a)};
end
end

function i = sortBy(idx, val)
[~, o] = sort(val(idx)); i = idx(o); i = i(:);
end

function vals = choiceValues(leaves, slot, param)
% the value each pipeline takes for one searched choice ('absent' when its
% step is skipped; the step type for an alternative)
vals = cell(numel(leaves), 1);
for i = 1:numel(leaves)
    vals{i} = 'absent';
    for q = 1:numel(leaves(i).path)
        in = leaves(i).path{q};
        if ~strcmp(in.slot, slot), continue; end
        if strcmp(param, '(alternative)'), vals{i} = in.type;
        elseif isfield(in.params, param), vals{i} = valText(in.params.(param)); end
    end
end
end

function t = valText(v)
if ischar(v) || isstring(v), t = char(v);
elseif isnumeric(v) || islogical(v), t = mat2str(v);
elseif iscell(v), t = strjoin(cellfun(@valText, v, 'UniformOutput', false), ',');
else, t = class(v);
end
end

function per = perObjective(cands, ref, opts, W, scale, indices)
% Total error of every candidate and objective under the trial weights
% W{c} (N x B per condition): candidates x B x objectives. The SME
% (Measure.smeBoot) and what blinks add (Measure.eyeBias) are those of
% each resample; the shape part of the bias is that of all trials.
nO = numel(ref.objectives); B = size(W{1}, 2);
if nargin < 6, indices = {}; end
per = zeros(numel(cands), B, nO);
blink = pipecompare.utils.fieldOr(ref, 'blink', {});
for k = 1:numel(cands)
    objs = cands(k).m.objectives;
    g = gainOf(cands(k), nO); sp = spreadOf(cands(k), nO);
    for o = 1:nO
        ob = objs(o); ob.segmented = pipecompare.utils.fieldOr(ref, 'segmented', false);
        sm = pipecompare.eval.Measure.smeBoot(ob, W, opts, 7919 * o, scale, indices) / g(o);
        nC = numel(W); b2 = zeros(1, B);
        for c = 1:nC
            eq = zeros(1, B);
            if strcmp(objs(o).kind, 'scalar') && numel(blink) >= o && ~isempty(blink{o})
                K = W{c};
                if isfield(objs(o), 'w') && ~isempty(objs(o).w), wt = objs(o).w{c}; wt(isnan(wt)) = 0; K = K .* wt; end
                [~, eq] = pipecompare.eval.Measure.eyeBias(objs(o).X{c}(:, 1), blink{o}{c}, K, W{c}, ...
                    pipecompare.utils.fieldOr(objs(o),'eyeBiasSEMode','legacy-total'));
            end
            oc = struct('estimate', objs(o).estimate(c));
            b2 = b2 + biasOf(oc, g(o), sp(o), eq) .^ 2;
        end
        per(k, :, o) = sqrt(sm .^ 2 + b2 / nC);
    end
end
end

function v = combine(per, name, ref)
% The objective by its name: 'composite' = RMS over objectives.
if strcmp(name, 'composite')
    v = sqrt(mean(per .^ 2, 3));
else
    v = per(:, :, strcmp(ref.objectives, name));
end
end

function in = innerOf(c, ref, objName)
% The choice candidate c made from the measure itself, for nestedHalf: the
% scores of every option of the last such step (components chosen before
% a rejection limit come first: their options cover the limit too),
% divided by the gain, with zeros where a trial is absent (M: NaN), the
% options that passed the signal check (ok) and c's trial weights (w, {}
% when not weighted; trials c did not keep get the median weight). []
% when c has none, its cutoffs were fine-tuned after it, or the stored
% scores of its own choice do not reproduce c's scores: the same trials
% kept (those up to one limit) and the same scores (cosine of the centred
% scores at least 0.999).
in = [];
if ~strcmp(c.status, 'ok') || ~isfield(c, 'steps') || isempty(c.steps), return; end
ty = cellfun(@(st) st.type, c.steps, 'UniformOutput', false);
has = cellfun(@(st) isfield(st, 'curves') && ~isempty(st.curves), c.steps);
refined = cellfun(@(st) isfield(st, 'cutoffRefined') && isequal(st.cutoffRefined, true), c.steps);
if ~any(has) || any(refined), return; end
k = find(has & strcmp(ty, 'icremove'), 1, 'last');
if isempty(k), k = find(has, 1, 'last'); end
cv = c.steps{k}.curves;
nO = numel(ref.objectives);
use = 1:nO;
if ~strcmp(objName, 'composite'), use = find(strcmp(ref.objectives, objName)); end
O = cv.options{cv.chosen};
objs = c.m.objectives(use);
if ~isequal(O.use(:)', use(:)') || ~all(strcmp({objs.kind}, 'scalar')), return; end
g = gainOf(c, nO); g = g(use);
if ~all(g > 0), return; end
nC = numel(ref.ids);
for q = 1:nC
    kept = ~isnan(objs(1).X{q}(:, 1)); pres = ~isnan(O.M{q});
    if nnz(kept) < 2 || any(kept & ~pres), return; end
    drop = pres & ~kept;
    if any(drop) && ~(max(O.M{q}(kept)) < min(O.M{q}(drop))), return; end
    for j = 1:numel(use)
        a = objs(j).X{q}(kept, 1); b = O.X{q}(kept, j);
        a = a - mean(a); b = b - mean(b);
        if norm(a) > 0 && (a' * b) / max(norm(a) * norm(b), eps) < 0.999, return; end
    end
end
w = {};
if isfield(objs(1), 'w') && ~isempty(objs(1).w)
    w = objs(1).w;
    for q = 1:nC
        w{q}(isnan(w{q})) = median(w{q}(~isnan(w{q})));
    end
end
opt = cv.options;
for j = 1:numel(opt)
    for q = 1:nC
        X = opt{j}.X{q};
        if ~opt{j}.gained, X = X ./ g; end
        X(isnan(X)) = 0; opt{j}.X{q} = X;
    end
end
% the shape part of each measure's bias in each condition (fixed: that of
% all trials) and the blink trials, for the total error
sp = spreadOf(c, nO); sp = sp(use);
shape = zeros(numel(use), nC);
for j = 1:numel(use), shape(j, :) = biasOf(struct('estimate', objs(j).estimate), g(j), sp(j), zeros(1, nC)); end
blink = pipecompare.utils.fieldOr(ref, 'blink', {});
bl = {};
if numel(blink) >= max(use) && all(~cellfun(@isempty, blink(use))), bl = blink(use); end
modes=arrayfun(@(j)pipecompare.utils.fieldOr(objs(j),'eyeBiasSEMode','legacy-total'),use,'UniformOutput',false);
in = struct('options', {opt}, 'ok', logical(cv.ok), 'w', {w}, 'shape', shape, 'blink', {bl},'biasVarianceMode',{modes});
end

function [vH, vO] = nestedHalf(in, H, O, scH, scO, ref, opts)
% The candidate's choice (innerOf) made again on the trials of half H
% (H{c}: 1 for the reference trials in it), as its step makes it on all
% of them, with the margin in place of the bootstrap: per option its best
% limit (a fixed limit: that one), the first option (fewest components)
% within equivalenceMargin of the best option, and in it the most lenient
% limit within equivalenceMargin of its best. Its objective on H (vH) and
% on the other half O (vO), rescaled by scH and scO like perObjective;
% NaN when no option keeps enough trials of H.
vH = NaN; vO = NaN;
nC = numel(ref.ids); nOpt = numel(in.options);
need = max(opts.minTrials, ceil(opts.minRetention * ref.n - 1e-9));
mg = 1 + opts.equivalenceMargin;
objs = cell(1, nOpt); lims = cell(1, nOpt); best = inf(1, nOpt);
for j = find(in.ok(:)')
    P = in.options{j};
    v = P.v;
    if isempty(v)
        v = [];
        for q = 1:nC, v = [v; P.M{q}(H{q} > 0)]; end %#ok<AGROW>
        v = unique(v(~isnan(v)))';
    end
    if isempty(v), continue; end
    v = v(:)';
    ok = true(1, numel(v)); acc = zeros(1, numel(v));
    for q = 1:nC
        K = double(P.M{q}(:) <= v) .* H{q}(:);           % trials x limits
        n = sum(K, 1);
        ok = ok & n >= max(ceil(need(q) * sum(H{q}) / ref.n(q) - 1e-9), 2);
        for o = 1:size(P.X{q}, 2)
            x = P.X{q}(:, o);
            if pipecompare.utils.fieldOr(ref,'segmented',false)
                [~, se] = pipecompare.eval.Measure.meanError(repmat(x,1,numel(v)),K,true);
                acc = acc + se.^2;
            else
                s1 = x' * K; s2 = (x .^ 2)' * K;
                acc = acc + max((s2 - s1 .^ 2 ./ n) ./ (n - 1), 0) ./ n;
            end
            acc = acc + halfBias(in, x, o, q, K) .^ 2;
        end
    end
    obj = sqrt(acc / (nC * size(P.X{1}, 2)));
    obj(~ok | ~isfinite(obj)) = Inf;
    objs{j} = obj; lims{j} = v; best(j) = min(obj);
end
if ~any(isfinite(best)), return; end
j = find(best <= mg * min(best), 1);
u = lims{j}(find(objs{j} <= mg * best(j), 1, 'last'));
P = in.options{j};
vH = scoredAt(in, P, u, H, scH, ref, opts); vO = scoredAt(in, P, u, O, scO, ref, opts);
end

function v = scoredAt(in, P, u, W, sc, ref, opts)
% The total error of option P with the trials kept at limit u, under the
% trial weights W (Measure.smeBoot, rescaled by sc; the bias as in
% halfBias; RMS over the measures)
nC = numel(P.X); per = zeros(1, size(P.X{1}, 2));
for o = 1:numel(per)
    X = cell(1, nC); b2 = 0;
    for q = 1:nC
        X{q} = P.X{q}(:, o); X{q}(~(P.M{q} <= u)) = NaN;
        K = W{q}; if ~isempty(in.w), K = K .* in.w{q}; end
        b2 = b2 + halfBias(in, X{q}, o, q, K) ^ 2;
    end
    ob = struct('kind','scalar','X',{X},'w',{in.w},'segmented',pipecompare.utils.fieldOr(ref,'segmented',false));
    per(o) = sqrt(pipecompare.eval.Measure.smeBoot(ob, W, opts, 0, sc) ^ 2 + b2 / nC);
end
v = sqrt(mean(per .^ 2));
end

function b = halfBias(in, x, o, q, K)
% The bias of measure o (of the options' measures) in condition q with
% weights K (a column per mean): what the blinks among those trials add
% (Measure.eyeBias; x is already divided by the gain) plus the shape part.
eq = zeros(1, size(K, 2));
if ~isempty(in.blink), [~, eq] = pipecompare.eval.Measure.eyeBias(x, in.blink{o}{q}, K,ones(size(K)),in.biasVarianceMode{o}); end
b = sqrt(eq) + in.shape(o, q);
end

