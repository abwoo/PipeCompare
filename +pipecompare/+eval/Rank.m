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
    %   4. Uncertainty: a paired bootstrap over trials (the same resampled
    %      trials, by urevent, for every candidate). A candidate is "not
    %      distinguished from the best" unless the data show it is worse
    %      than some other candidate, with intervals simultaneous over all
    %      candidate pairs: the probability that ANY candidate is wrongly
    %      declared worse stays <= ~5% for any number of candidates (alpha
    %      = 0.02 chosen by simulation for 2-24 candidates, 20-100 trials,
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
    %   5. Among candidates the data do not distinguish from the best, the
    %      recommendation is the least aggressive: highest minimum trial
    %      retention, then smallest signal distortion, then best objective.
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
    %   No experimental effect (condition difference, p-value) is used, so
    %   choosing a pipeline cannot inflate the effect tested later.

    methods (Static)
        function o = defaults()
            o = struct('minTrials', 10, 'minRetention', 0.5, 'maxInterpolated', 0.2, ...
                'maxAmplitudeError', 0.10, 'maxLatencyShiftMs', 10, 'maxArtifactPct', 0.05, ...
                'minWaveformCorr', 0.95, 'minTopoCorr', 0.90, ...
                'objective', 'composite', 'nBoot', 1999, 'nBootPeakOuter', 999, 'nBootPeakInner', 100, ...
                'nBootPeak', 1000, 'alpha', 0.02, 'seed', 1, 'equivalenceMargin', 0.05, 'nSplits', 20);
        end

        function R = run(cands, ref, opts)
            if nargin < 3, opts = struct(); end
            opts = pipecompare.utils.withDefaults(opts, pipecompare.eval.Rank.defaults());
            n = numel(cands);
            objName = resolveObjective(opts.objective, ref);
            nO = numel(ref.objectives);
            status = repmat({''}, n, 1); reasons = repmat({''}, n, 1); whyList = repmat({{}}, n, 1);
            notes = repmat({''}, n, 1); gains = nan(n, nO);
            stratum = repmat({''}, n, 1);
            primary = nan(n, 1); minRet = nan(n, 1); minKept = nan(n, 1); interp = nan(n, 1);
            ampErr = nan(n, 1); latSh = nan(n, 1); artPct = nan(n, 1); wCorr = nan(n, 1); tCorr = nan(n, 1);
            objAgg = nan(n, nO);
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
                g = gainOf(c, nO); gains(i, :) = g;
                objAgg(i, :) = [m.objectives.agg] ./ g;
                primary(i) = primaryPoint(m, objName, ref, g);
                minRet(i) = min(m.retention); minKept(i) = min(m.kept);
                interp(i) = c.interpolatedFraction;
                sg = c.signal;
                if ~isempty(sg)
                    ampErr(i) = sg.amplitudeError; latSh(i) = sg.latencyShiftMs; artPct(i) = sg.artifactPct;
                    wCorr(i) = sg.waveformCorr; tCorr(i) = sg.topoCorr;
                end
                why = {};
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
            byStratum = struct('stratum', {}, 'best', {}, 'recommended', {}, 'set', {}, 'crossfit', {}, 'stability', {});
            chosen = nan(n, 1);
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
                % least aggressive among those not distinguished from the best
                key = [-minRet(t), ampErr(t), primary(t)];
                [~, o] = sortrows(round(key, 10));
                cf = pipecompare.eval.Rank.crossfit(cands(feas), ref, opts, objName, pts(ib));
                f = selectionFrequency(boot{1}, crit, minRet(feas), ampErr(feas));
                chosen(feas) = f;
                byStratum(end+1) = struct('stratum', strata{s}, 'best', feas(ib), 'recommended', t(o(1)), 'set', t, ...
                    'crossfit', cf, 'stability', f(feas == t(o(1)))); %#ok<AGROW>
            end
            if isscalar(byStratum)
                recommended = byStratum.recommended; best = byStratum.best;
            end
            feasAll = find(strcmp(status, 'feasible'));
            order = [sortBy(feasAll, primary); find(strcmp(status, 'rejected')); find(strcmp(status, 'failed'))];
            T = table((1:n)', stratum, status, primary, lo, hi, notDist, equiv, hiKept, chosen, minRet, minKept, interp, ...
                ampErr, latSh, artPct, wCorr, tCorr, reasons, notes, ...
                'VariableNames', {'id','stratum','status','objective','diffLo','diffHi','notDistinguished','equivalent','diffHiKept', ...
                'chosenShare','minRetention','minTrials','interpolated','ampError', ...
                'latencyShiftMs','artifactPct','waveformCorr','topoCorr','reason','note'});
            for k = 1:nO
                T.(matlab.lang.makeValidName(['sme_' ref.objectives{k}])) = objAgg(:, k);   % gain-corrected
                T.(matlab.lang.makeValidName(['gain_' ref.objectives{k}])) = gains(:, k);
            end
            R = struct('table', T, 'whyList', {whyList}, 'order', order, 'best', best, 'recommended', recommended, ...
                'byStratum', byStratum, 'objective', objName, 'options', opts, 'units', {ref.units});
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
            nC = numel(ref.ids); W = cell(1, nC); scale = ones(1, nC);
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
                    W{c} = pipecompare.eval.Rank.resampleWeights(s, N, B, isfield(ref, 'segmented') && ref.segmented);
                end
            end
            per = perObjective(cands, ref, opts, W, scale);
            boot = cellfun(@(q) combine(per, q, ref), objNames, 'UniformOutput', false);
        end

        function W = resampleWeights(s, N, B, segmented)
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

        function v = pointObjective(m, ref, opts, g)
            % The ranking's objective of one candidate's measures m with
            % the signal gains g (one per measure): what Rank.run ranks by.
            opts = pipecompare.utils.withDefaults(opts, pipecompare.eval.Rank.defaults());
            v = primaryPoint(m, resolveObjective(opts.objective, ref), ref, g);
        end

        function name = objectiveName(obj, ref)
            % 'composite' or the one measure the objective option names
            name = resolveObjective(obj, ref);
        end

        function cf = crossfit(cands, ref, opts, objName, apparent)
            % The objective of the pipeline that is best on one half of the
            % trials, scored on the other half (rescaled to all trials by
            % sqrt(m/N), as the half-samples of bootstrap; SME ~ 1/sqrt(n)),
            % averaged over both directions and opts.nSplits random splits.
            % Selection and evaluation use different trials, so the value is
            % not flattered by the selection; apparent (the best objective on
            % all trials) is. overstatement = honest / apparent - 1. A choice
            % made on half the trials is somewhat worse than one made on all,
            % so the overstatement is, if anything, too large. [] when a
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
            cf = [];
            if nargin < 5, apparent = NaN; end
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
                [~, ia] = min(vAA); [~, ib] = min(vBB);
                honest = [honest vAB(ia) vBA(ib)]; inSample = [inSample vAA(ia) vBB(ib)]; %#ok<AGROW>
            end
            ok = isfinite(honest);
            if ~any(ok), return; end
            h = mean(honest(ok));
            cf = struct('apparent', apparent, 'honest', h, 'overstatement', h / apparent - 1, ...
                'halfInSample', mean(inSample(ok)), 'nSplits', nS, 'nested', numel(nest));
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
            pipecompare.utils.log(['Ranking by %s gain-corrected SME (SME / signal gain; lower = more precise). diff = difference from the best, ', ...
                '%.0f%% interval simultaneous over all candidates; "nd" = not distinguished from the best by these data ', ...
                '(not equivalence); "eq" = equivalent: the whole interval within %.0f%% of the best.'], R.objective, ...
                100 * (1 - o.alpha), 100 * o.equivalenceMargin);
            fprintf('   %-4s %-9s %9s %19s %3s %6s %6s %6s %5s  %s\n', 'id', 'status', 'objective', 'diff vs best [CI]', 'nd', 'minRet', 'interp', 'ampErr', 'art', 'pipeline');
            for k = R.order(:)'
                mark = ' '; if any([R.byStratum.recommended] == k), mark = '*'; end
                ci = ''; if isfinite(T.diffLo(k)), ci = sprintf('[%+.3f %+.3f]', T.diffLo(k), T.diffHi(k)); end
                nd = ''; if T.equivalent(k), nd = 'eq'; elseif T.notDistinguished(k), nd = 'nd'; end
                fprintf('  %s%-4d %-9s %9.4g %19s %3s %5.0f%% %5.0f%% %5.0f%% %4.1f%%  %s\n', mark, k, T.status{k}, ...
                    T.objective(k), ci, nd, 100*T.minRetention(k), 100*T.interpolated(k), ...
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
                if isfinite(b.stability)
                    pipecompare.utils.log('Stability: in resamples of the trials, candidate %d is recommended in %.0f%% of them.', ...
                        b.recommended, 100 * b.stability);
                end
                if ~isempty(b.crossfit)
                    pipecompare.utils.log(['Selection check: the pipeline best on half of the trials has objective %.4g on the ', ...
                        'other half (rescaled to all trials), against %.4g for the best on all trials: the best looks %.0f%% ', ...
                        'better than it is (%d split(s), both directions; the limits, thresholds or components chosen ', ...
                        'from the data made again on the selecting half for %d pipeline(s)).'], b.crossfit.honest, b.crossfit.apparent, ...
                        100 * b.crossfit.overstatement, b.crossfit.nSplits, pipecompare.utils.fieldOr(b.crossfit, 'nested', 0));
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

function v = primaryPoint(m, name, ref, g)
% gain-corrected objective: RMS over measures and conditions of SME/g
if strcmp(name, 'composite')
    v = sqrt(mean(cell2mat(arrayfun(@(k) m.objectives(k).sme / g(k), 1:numel(g), 'UniformOutput', false)) .^ 2));
else
    k = strcmp(ref.objectives, name);
    v = m.objectives(k).agg / g(k);
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
B = size(Bq, 2); maxD = -inf(1, B);
Dhat = pts - pts';                             % K x K: i minus j
for i = 1:K
    D = Bq(i, :) - Bq;                         % K x B: candidate i minus each j
    D(~isfinite(D)) = NaN;
    maxD = max(maxD, max(D - Dhat(i, :)', [], 1, 'omitnan'));
end
z = maxD(isfinite(maxD));
if isempty(z), c = Inf; else, c = upperQuantile(z, alpha); end
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
z = maxK(isfinite(maxK));
if isempty(z), cK = 0; else, cK = max(upperQuantile(z, alpha), 0); end
hiEq = Dhat(:, ib) + cK; hiEq(ib) = 0; hiEq(~keep) = Inf;
end

function f = selectionFrequency(Bq, c, ret, amp)
% The share of the bootstrap resamples (columns of Bq, K x B) in which
% each of K candidates is recommended: those within c of the lowest
% objective of the resample, then the rule of step 5 (highest minimum
% retention, smallest amplitude error, lowest objective).
K = size(Bq, 1); cnt = zeros(K, 1); used = 0;
for b = 1:size(Bq, 2)
    v = Bq(:, b); ok = isfinite(v);
    if ~any(ok), continue; end
    t = find(ok & v - min(v(ok)) <= c);
    [~, o] = sortrows(round([-ret(t), amp(t), v(t)], 10));
    cnt(t(o(1))) = cnt(t(o(1))) + 1; used = used + 1;
end
f = cnt / max(used, 1);
if used == 0, f(:) = NaN; end
end

function v = upperQuantile(z, alpha)
% The (1-alpha) quantile of the bootstrap distribution as the order
% statistic z_((B+1)(1-alpha)) (Davison & Hinkley, 1997): exact when
% (B+1)*alpha is an integer, rounded up (conservative) otherwise.
z = sort(z(:));
r = min(numel(z), ceil((numel(z) + 1) * (1 - alpha) - 1e-9));
v = z(max(r, 1));
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

function per = perObjective(cands, ref, opts, W, scale)
% Gain-corrected SME of every candidate and objective under the trial
% weights W{c} (N x B per condition): candidates x B x objectives.
nO = numel(ref.objectives); B = size(W{1}, 2);
per = zeros(numel(cands), B, nO);
for k = 1:numel(cands)
    objs = cands(k).m.objectives;
    g = gainOf(cands(k), nO);
    for o = 1:nO
        per(k, :, o) = pipecompare.eval.Measure.smeBoot(objs(o), W, opts, 7919 * o, scale) / g(o);
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
in = struct('options', {opt}, 'ok', logical(cv.ok), 'w', {w});
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
            s1 = x' * K; s2 = (x .^ 2)' * K;
            acc = acc + max((s2 - s1 .^ 2 ./ n) ./ (n - 1), 0) ./ n;
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
vH = scoredAt(P, u, H, scH, in.w, opts); vO = scoredAt(P, u, O, scO, in.w, opts);
end

function v = scoredAt(P, u, W, sc, w, opts)
% The objective of option P with the trials kept at limit u, under the
% trial weights W (Measure.smeBoot, rescaled by sc; RMS over the measures)
nC = numel(P.X); per = zeros(1, size(P.X{1}, 2));
for o = 1:numel(per)
    X = cell(1, nC);
    for q = 1:nC
        X{q} = P.X{q}(:, o); X{q}(~(P.M{q} <= u)) = NaN;
    end
    per(o) = pipecompare.eval.Measure.smeBoot(struct('kind', 'scalar', 'X', {X}, 'w', {w}), W, opts, 0, sc);
end
v = sqrt(mean(per .^ 2));
end
