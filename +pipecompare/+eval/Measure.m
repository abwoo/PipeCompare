classdef Measure
    %MEASURE Trial-level scores and their measurement error, per objective.
    %
    %   ref = pipecompare.eval.Measure.reference(rootEEG, contract)
    %   m   = pipecompare.eval.Measure.candidate(EEG, contract, ref, opts)
    %
    %   Every candidate is scored the same way. Data still continuous are
    %   epoched with the contract window; the contract baseline is removed
    %   (idempotent when the pipeline already did it). Then one objective is
    %   computed per contract component:
    %
    %     mean          trial score = mean amplitude over ROI and window;
    %                   precision = analytic SME = SD / sqrt(N)  (uV)
    %     peakAmplitude precision of the peak amplitude of the averaged
    %     peakLatency   ROI waveform = bootstrapped SME (bSME: SD of the
    %                   measure over bootstrap averages; uV or ms)
    %
    %   A component scored contralateral minus ipsilateral (N2pc, LRP; see
    %   pipecompare.eval.Contract) has, per trial, the electrode
    %   contralateral to the trial's condition minus the other one as its
    %   waveform: the SME is that of the difference wave.
    %
    %   (Luck, Stewart, Simmons & Rhemtulla, 2021, Psychophysiology.)
    %   SME falls when noise is removed and rises when trials are lost, so it
    %   prices the noise-vs-trial-count trade-off directly. It does NOT tell
    %   whether the signal survived; that is the job of
    %   pipecompare.eval.Injection.
    %
    %   Trials are identified by the urevent index of their time-locking
    %   event, so the same physical trial is paired across candidates.
    %
    %   Epochs weighted by a weight_trials step (pipecompare.run.Steps.
    %   weightTrials) are averaged with their weights w: estimate =
    %   sum(w .* x) / sum(w), and SME = sqrt(n / (n - 1) * sum(w.^2 .* (x -
    %   estimate).^2)) / sum(w), the sandwich (robust) standard error of a
    %   weighted mean conditional on the fitted weights; equal weights give
    %   SD / sqrt(n). Adaptive weights and dependent trials still need
    %   separate validation. Peak measures average with the weights.
    %
    %   w = pipecompare.eval.Measure.epochWeights(EEG)
    %   The weight of each epoch (1 x trials; all 1 without weighting): the
    %   weighted average is sum(EEG.data .* reshape(w, 1, 1, []), 3) / sum(w).
    %
    %   A = pipecompare.eval.Measure.weightedAverages(EEG, contract)
    %   That weighted average of each condition: A.data (channels x samples
    %   x conditions), A.conditions, A.labels, A.times (ms), A.trials and
    %   A.sumWeights per condition. EEGLAB's and ERPLAB's own averages count
    %   every epoch equally, so they do not give the ERP a weighted
    %   pipeline was scored on.

    properties (Constant)
        BlinkZ = 3   % an epoch's eye activity this many robust SDs above the others' is a blink
        BiasZ = 2    % what blinks add counts only beyond this many of its standard errors
    end

    methods (Static)
        function o = defaults()
            o = struct('nBootPeak', 1000, 'seed', 1);
        end

        function ref = reference(EEG, contract)
            % trial identities only: the ROI may be restored later in the plan
            T = pipecompare.eval.Measure.trials(EEG, contract, true);
            conds = contract.conditions;
            ref = struct('ids', {cell(1, numel(conds))}, 'names', {{conds.name}});
            for c = 1:numel(conds)
                ref.ids{c} = unique(T.id(T.cond == c))';
                assert(numel(ref.ids{c}) >= 2, 'PipeCompare:Contract', ...
                    'Condition %s has %d eligible trial(s) in the current dataset.', conds(c).name, numel(ref.ids{c}));
            end
            ref.n = cellfun(@numel, ref.ids);
            ref.objectives = contract.objectiveNames();
            ref.units = contract.objectiveUnits();
            % consecutive segments of one recording are not independent:
            % Rank resamples them in blocks (moving-block bootstrap)
            ref.segmented = contract.isSegmented();
            [ref.blink, ref.blinkInfo] = pipecompare.eval.Measure.blinks(EEG, contract, ref);
        end

        function [B, info] = blinks(EEG, contract, ref)
            % Which reference trials have an artifact that can move the mean
            % of a mean-amplitude component's score: B{j}{c} (ref.n(c) x 1
            % logical) for component j and condition c, {} when it cannot
            % be told (band power, or no channel to tell it from). Three
            % kinds, in the data as they are (before any processing), each
            % found when its robust z (median and MAD over all reference
            % trials) exceeds BlinkZ; a trial with any of them counts:
            %   blinks: the largest, over the eye channels
            %   (Injection.eyeChannels: EOG by type or name, else Fp1, Fp2,
            %   FPz; never a measured electrode), of the absolute window
            %   mean minus the baseline mean (what an eye artifact adds to
            %   the score is carried by that difference);
            %   eye movements: the same of the difference of a left and a
            %   right channel (horizontalPair), in which a blink, the same
            %   on both sides, cancels and a horizontal movement doubles;
            %   muscle: the log of the largest, over the measured
            %   electrodes, SD of the first difference in the window (fast
            %   activity, which a slow ERP hardly has).
            % In continuous data, a trial also counts when a blink or eye
            % movement found sample by sample (sampleHits) falls in its
            % window or baseline: comparing the trials with each other
            % finds few when most trials have a blink (n-back tasks with a
            % response on every trial).
            % info: channels (eye channels' labels), horizontal ('L-R' or
            % ''), share (of the trials, per component), kinds and
            % kindShare (components x kinds).
            kinds = {'blinks', 'eye movements', 'muscle'};
            B = {}; info = struct('channels', {{}}, 'horizontal', '', 'share', [], 'kinds', {kinds}, 'kindShare', [], ...
                'invalidEyeInput', false, 'reason', '', 'available', {{}});
            if ~strcmp(contract.analysis, 'erp') || isempty(contract.components), return; end
            labs = lower({EEG.chanlocs.labels});
            roi = find(ismember(labs, lower(contract.allRoi())));
            eyes = setdiff(pipecompare.eval.Injection.eyeChannels(EEG.chanlocs), roi);
            [hl, hr] = horizontalPair(EEG.chanlocs, eyes, roi);
            use = unique([eyes hl hr roi]);
            if isempty(use), return; end
            info.channels = {EEG.chanlocs(eyes).labels};
            if ~isempty(eyes) && ~any(isfinite(EEG.data(eyes,:,:)), 'all')
                % A declared eye channel with no observations is unavailable
                % evidence, rather than observed absence of ocular events.
                info.invalidEyeInput = true;
                info.reason = 'no finite eye-channel samples';
                return;
            end
            if ~isempty(hl), info.horizontal = sprintf('%s-%s', EEG.chanlocs(hl).labels, EEG.chanlocs(hr).labels); end
            E = EEG; E.icaact = []; E.icaweights = []; E.icasphere = []; E.icawinv = []; E.icachansind = [];
            [~, E] = evalc('pop_select(E, ''channel'', use)');
            H = [];
            if E.trials == 1
                H = sampleHits(E, ismember(use, eyes), ismember(use, hl), ismember(use, hr), contract);
                [~, E] = evalc('pop_epoch(E, contract.allEvents(), contract.epoch, ''epochinfo'', ''yes'')');
            end
            T = pipecompare.eval.Measure.trials(E, contract, true);
            times = E.xmin + (0:E.pnts-1) / E.srate;
            D = double(E.data);
            De = D(ismember(use, eyes), :, :); Dh = [];
            if ~isempty(hl), Dh = D(use == hl, :, :) - D(use == hr, :, :); end
            bsel = false(size(times));
            if ~isempty(contract.baseline)
                bsel = times >= contract.baseline(1) - 1e-9 & times <= contract.baseline(2) + 1e-9;
                if any(bsel)
                    De = De - mean(De(:, bsel, :), 2);
                    if ~isempty(Dh), Dh = Dh - mean(Dh(:, bsel, :), 2); end
                end
            end
            z = @(v) (v - median(v)) / max(1.4826 * median(abs(v - median(v))), eps);
            perTrial = @(v) reshape(v, [], 1);
            nJ = numel(contract.components); nC = numel(ref.ids);
            info.available = cell(1,nJ);
            B = cell(1, nJ); info.share = nan(1, nJ); info.kindShare = nan(nJ, numel(kinds));
            for j = 1:nJ
                B{j} = arrayfun(@(c) false(ref.n(c), 1), 1:nC, 'UniformOutput', false);
                info.available{j} = arrayfun(@(c) false(ref.n(c),1),1:nC,'UniformOutput',false);
                comp = contract.components(j);
                w = times >= comp.window(1) - 1e-9 & times <= comp.window(2) + 1e-9;
                if ~strcmp(comp.measure, 'mean') || ~any(w), continue; end
                available=true(E.trials,1);
                if ~isempty(eyes)
                    available=reshape(all(all(isfinite(D(ismember(use,eyes),w|bsel,:)),1),2),[],1);
                end
                K = false(numel(T.epoch), numel(kinds));
                if ~isempty(De)
                    v = perTrial(max(abs(mean(De(:, w, :), 2)), [], 1)); K(:, 1) = z(v(T.epoch)) > pipecompare.eval.Measure.BlinkZ;
                end
                if ~isempty(Dh)
                    v = perTrial(abs(mean(Dh(:, w, :), 2))); K(:, 2) = z(v(T.epoch)) > pipecompare.eval.Measure.BlinkZ;
                end
                if ~isempty(H)
                    % also the trials with a blink or eye movement found in
                    % the continuous data, in the window or the baseline
                    hit = reshape(any(H.data(:, w | bsel, :), 2), size(H.data, 1), []);
                    K(:, 1) = K(:, 1) | hit(1, T.epoch)';
                    if size(hit, 1) > 1, K(:, 2) = K(:, 2) | hit(2, T.epoch)'; end
                end
                r = ismember(use, roi(ismember(labs(roi), lower(comp.roi))));
                if any(r) && sum(w) > 2
                    v = perTrial(log(max(std(diff(D(r, w, :), 1, 2), 0, 2), [], 1) + eps));
                    K(:, 3) = z(v(T.epoch)) > pipecompare.eval.Measure.BlinkZ;
                end
                isB = any(K, 2);
                for c = 1:nC
                    rr = find(T.cond == c);
                    [tf, loc] = ismember(T.id(rr), ref.ids{c});
                    B{j}{c}(loc(tf)) = isB(rr(tf));
                    info.available{j}{c}(loc(tf)) = available(T.epoch(rr(tf)));
                end
                info.share(j) = mean(isB); info.kindShare(j, :) = mean(K, 1);
            end
        end

        function [bias, sq, assessed] = eyeBias(x, b, K, counts)
            % What the blink trials add to a weighted mean score. x: scores
            % (n x 1, NaN = absent), b: blink trials (n x 1 logical), K:
            % weights (n x L: a column per mean). Optional counts are the
            % physical resampling multiplicities, before trial weights;
            % without them each positive weight is one observation. Pass
            % counts separately for bootstrap draws, so duplicating a trial
            % differs from increasing its reliability weight. The
            % mean is sum(K .* x) / sum(K); with blink share p (of the
            % weight) and blink and other trials' means mB and mN, it is mN
            % + p (mB - mN): bias = p (mB - mN), what the blinks add. The
            % estimate is noisy (se = p s sqrt(qB/swB^2 + qN/swN^2),
            % q the sums of squared weights per physical observation, s
            % the pooled weighted SD, often larger than the SME itself when
            % most trials have a blink), so what counts is only the part
            % beyond what chance makes: sq = max(|bias| - BiasZ se, 0)^2
            % (soft thresholding). Taking only the variance off (bias^2 -
            % se^2, kept at 0 or more) left a bias where there is none: in
            % simulation with 70% blink trials and no real effect it added a
            % third to the total error on average and pushed the choices
            % towards rejecting blink trials; with soft thresholding about
            % 1%. Zero when no kept trial has a detected artifact. NaN
            % when every kept trial has one, or no trial remains: there is
            % no artifact-free comparison and the bias is unidentifiable.
            % assessed says which weighted means have an estimate. 1 x L.
            if nargin < 4, counts = ones(size(K)); end
            assert(all(isfinite(K(:)) & K(:) >= 0), 'PipeCompare:Weights', 'Bias weights must be finite and nonnegative.');
            assert(all(isfinite(counts(:)) & counts(:) >= 0 & counts(:) == fix(counts(:))), ...
                'PipeCompare:Weights', 'Resampling counts must be nonnegative integers.');
            counts = double(counts) + zeros(size(K));
            assert(all(K(counts == 0) == 0), 'PipeCompare:Weights', 'A positive weight needs a physical observation.');
            ok = isfinite(x); x(~ok) = 0; b = double(b(:) & ok(:)); o = double(ok(:)) - b;
            K = double(K) .* double(ok(:));
            Q = K .^ 2 ./ max(counts, 1);
            nB = b' * K; nN = o' * K; n = nB + nN;
            sB = (x .* b)' * K; sN = (x .* o)' * K;
            s2 = (x .^ 2)' * K;
            mu = (sB + sN) ./ n;
            v = max(s2 ./ n - mu .^ 2, 0) ./ max(1 - sum(Q,1) ./ n .^ 2, eps);
            p = nB ./ n;
            bias = p .* (sB ./ nB - sN ./ nN);
            se = p .* sqrt(v .* ((b' * Q) ./ nB .^ 2 + (o' * Q) ./ nN .^ 2));
            sq = max(abs(bias) - pipecompare.eval.Measure.BiasZ * se, 0) .^ 2;
            clean = nB == 0 & nN > 0;
            assessed = nN > 0;
            bias(clean) = 0; sq(clean) = 0;
            bias(~assessed) = NaN; sq(~assessed) = NaN;
        end

        function m = candidate(EEG, contract, ref, opts)
            if nargin < 4, opts = struct(); end
            opts = pipecompare.utils.withDefaults(opts, pipecompare.eval.Measure.defaults());
            T = pipecompare.eval.Measure.trials(EEG, contract);
            nC = numel(ref.ids); nO = numel(ref.objectives);
            m = struct('kept', zeros(1, nC), 'retention', zeros(1, nC), 'extraTrials', 0, ...
                'objectives', [], 'composite', NaN);
            % align rows of every condition to the reference trial ids
            rowOf = cell(1, nC);
            for c = 1:nC
                rows = find(T.cond == c);
                [tf, loc] = ismember(T.id(rows), ref.ids{c});
                m.extraTrials = m.extraTrials + sum(~tf);
                rows = rows(tf); loc = loc(tf);
                [loc, first] = unique(loc, 'stable'); rows = rows(first);
                rowOf{c} = struct('rows', rows, 'loc', loc);
                m.kept(c) = numel(loc);
                m.retention(c) = m.kept(c) / ref.n(c);
            end
            objs = repmat(struct('name', '', 'unit', '', 'kind', '', 'polarity', '', 'times', [], ...
                'X', {{}}, 'w', {{}}, 'estimate', [], 'sme', [], 'agg', NaN, 'eyeBias', [], 'eyeBiasSq', [], ...
                'eyeBiasAssessed', [], 'segmented', false), 1, nO);
            blink = pipecompare.utils.fieldOr(ref, 'blink', {});
            % the epochs' weights (weight_trials), by the trial's urevent
            wRow = [];
            tw = pipecompare.eval.Measure.storedWeights(EEG);
            if ~isempty(tw)
                [ok, at] = ismember(T.id, tw.urevent);
                assert(all(ok), 'PipeCompare:Weights', 'Some epochs have no weight (weight_trials).');
                wRow = tw.weight(at); wRow = wRow(:);
            end
            for k = 1:nO
                o = objs(k);
                o.name = ref.objectives{k}; o.unit = ref.units{k};
                [o.kind, o.polarity] = objKind(contract, k);
                o.segmented = pipecompare.utils.fieldOr(ref, 'segmented', false);
                o.times = T.times{k};
                o.X = cell(1, nC); o.estimate = nan(1, nC); o.sme = nan(1, nC);
                % what blinks add to each condition's mean (Measure.eyeBias;
                % 0 without eye channels, for peak measures and band power)
                o.eyeBias = zeros(1, nC); o.eyeBiasSq = zeros(1, nC);
                if strcmp(o.kind, 'scalar') && numel(blink) >= k && ~isempty(blink{k})
                    o.eyeBiasAssessed = false(1, nC);
                end
                if ~isempty(wRow), o.w = cell(1, nC); end
                for c = 1:nC
                    X = nan(ref.n(c), size(T.data{k}, 2));
                    X(rowOf{c}.loc, :) = T.data{k}(rowOf{c}.rows, :);
                    o.X{c} = X;
                    w = [];
                    if ~isempty(wRow)
                        w = nan(ref.n(c), 1); w(rowOf{c}.loc) = wRow(rowOf{c}.rows); o.w{c} = w;
                    end
                    [o.estimate(c), o.sme(c)] = pipecompare.eval.Measure.estimate(o, X, opts, k * 1000 + c, w);
                    if strcmp(o.kind, 'scalar') && numel(blink) >= k && ~isempty(blink{k})
                        K = ones(ref.n(c), 1); if ~isempty(w), K = w; K(isnan(K)) = 0; end
                        [o.eyeBias(c), o.eyeBiasSq(c), o.eyeBiasAssessed(c)] = pipecompare.eval.Measure.eyeBias(X(:, 1), blink{k}{c}, K);
                    end
                    if strcmp(o.kind, 'scalar') && isfield(ref, 'blinkInfo')
                        invalid=pipecompare.utils.fieldOr(ref.blinkInfo,'invalidEyeInput',false);
                        if isfield(ref.blinkInfo,'available') && numel(ref.blinkInfo.available)>=k && ...
                                ~isempty(ref.blinkInfo.available{k})
                            kept=isfinite(X(:,1));if ~isempty(w),kept=kept&w>0;end
                            invalid=invalid||any(kept&~ref.blinkInfo.available{k}{c}(:));
                        end
                        if invalid
                            o.eyeBias(c)=NaN;o.eyeBiasSq(c)=NaN;o.eyeBiasAssessed(c)=false;
                        end
                    end
                end
                o.agg = sqrt(mean(o.sme .^ 2));
                objs(k) = o;
            end
            m.objectives = objs;
            if numel(unique(ref.units)) == 1
                m.composite = sqrt(mean([objs.sme] .^ 2));
            end
        end

        function [s, resolution] = bandScores(X, fs, band)
            % X: ROI channels x samples x trials. Per trial: log10 of the mean
            % power in the band (multitaper: the one-sided power spectral
            % density averaged over the tapers of Measure.tapers, the band's
            % bins and the ROI channels).
            [~, n, nt] = size(X);
            assert(isscalar(fs) && isfinite(fs) && fs>0 && numel(band)==2 && all(isfinite(band)) && ...
                band(1)>=0 && band(1)<band(2) && band(2)<=fs/2, ...
                'PipeCompare:Measure', 'A one-sided PSD band must lie between zero and Nyquist.');
            [H, resolution] = pipecompare.eval.Measure.tapers(n, fs, band);
            f = (0:n-1) * fs / n;
            % Locate the Nyquist bin by index. Fractional sample rates can
            % round its computed frequency above/below fs/2 by one ulp.
            if mod(n,2)==0, f(n/2+1)=fs/2; end
            sel = f >= band(1) & f <= band(2);
            assert(any(sel), 'PipeCompare:Measure', 'The window is too short to resolve the band [%g %g] Hz.', band(1), band(2));
            Xp = permute(X, [2 1 3]);                              % samples x ROI x trials
            P = 0;
            sided = 2 * ones(nnz(sel), 1);
            bins=find(sel);
            sided(bins == 1 | (mod(n,2) == 0 & bins == n/2+1)) = 1;
            for k = 1:size(H, 2)
                F = fft(Xp .* H(:, k), [], 1);
                P = P + sided .* abs(F(sel, :, :)) .^ 2 / (fs * sum(H(:, k) .^ 2));
            end
            s = reshape(log10(mean(mean(P / size(H, 2), 1), 2)), nt, 1);
        end

        function [H, resolution] = tapers(n, fs, band)
            % Discrete prolate spheroidal (Slepian) tapers for n samples at fs
            % Hz: the orthonormal windows whose power is most concentrated
            % within +-W Hz (eigenvectors of Percival and Walden's
            % tridiagonal matrix, the largest eigenvalues first). Time-
            % bandwidth NW = W T = 2 (3 tapers), lowered for a band narrower
            % than 4 / T so that the smoothing (2 W) stays within the band:
            % NW = 1.5 gives 2 tapers, NW = 1 a single one. Cached per size.
            persistent cache
            if isempty(cache), cache = containers.Map('KeyType', 'char', 'ValueType', 'any'); end
            NW = max(1, min(2, (band(2) - band(1)) * n / fs / 2));
            K = max(1, floor(2 * NW + 1e-9) - 1);
            resolution = struct('NW',NW,'tapers',K,'binWidth',fs/n,'smoothingWidth',2*NW*fs/n, ...
                'bandWidth',diff(band),'resolved',diff(band)>=2*NW*fs/n-1e-12,'unit','log10(uV^2/Hz)');
            key = sprintf('%d_%.17g', n, NW);
            if isKey(cache, key), H = cache(key); return; end
            k = (0:n-1)';
            d = ((n - 1 - 2 * k) / 2) .^ 2 * cos(2 * pi * NW / n);
            e = k(2:end) .* (n - k(2:end)) / 2;
            M = spdiags([[e; 0] d [0; e]], -1:1, n, n);
            % each eigenvalue exactly by bisection on the Sturm count, its
            % vector by inverse iteration (as LAPACK's dstebz and dstein);
            % Krylov solvers (eigs) do not converge for long segments, the
            % gaps between the largest eigenvalues being small relative to them
            e2 = e .^ 2; lo = min(d) - 2 * max([e; 0]); hi = max(d) + 2 * max([e; 0]);
            H = zeros(n, K);
            for j = 1:K
                a = lo; b = hi;                       % the j-th largest: n - j eigenvalues below it
                for it = 1:200
                    m = (a + b) / 2;
                    if countBelow(d, e2, m) <= n - j, a = m; else, b = m; end
                    if b - a <= 1e-14 * max(abs(a), abs(b)), break; end
                end
                lam = (a + b) / 2;
                A = M - (lam + 1e-8 * abs(lam)) * speye(n);
                v = 1 + (0:n-1)' / n;
                for it = 1:4
                    v = A \ v;
                    v = v - H(:, 1:j-1) * (H(:, 1:j-1)' * v);
                    v = v / norm(v);
                end
                H(:, j) = v;
            end
            cache(key) = H;
        end

        function v = sme(S)
            % Analytic SME of a mean score per column, ignoring missing (NaN) trials.
            v = nan(1, size(S, 2));
            for j = 1:size(S, 2)
                x = S(~isnan(S(:, j)), j);
                if numel(x) >= 2, v(j) = std(x) / sqrt(numel(x)); end
            end
        end

        function [est, sme] = estimate(o, X, opts, salt, w)
            % Point estimate and its SME for one condition (w: the trials'
            % weights, [] = equal).
            if nargin < 5, w = []; end
            if strcmp(o.kind, 'scalar') && pipecompare.utils.fieldOr(o, 'segmented', false)
                [est, sme] = pipecompare.eval.Measure.meanError(X, w, true); return;
            end
            keep = ~isnan(X(:, 1)); Y = X(keep, :); n = size(Y, 1);
            est = NaN; sme = NaN;
            if n < 2, return; end
            if nargin > 4 && ~isempty(w)
                w = w(keep); w = w(:); sw = sum(w);
                switch o.kind
                    case 'scalar'
                        est = (w' * Y) / sw;
                        sme = sqrt(n / (n - 1) * sum(w .^ 2 .* (Y - est) .^ 2)) / sw;
                    otherwise
                        est = peakOf((w' * Y) / sw, o);
                        s = RandStream('mt19937ar', 'Seed', opts.seed + salt);
                        idx = randi(s, n, n, opts.nBootPeak);
                        W = full(sparse(idx, repmat(1:opts.nBootPeak, n, 1), 1, n, opts.nBootPeak));
                        avgs = (W' * (w .* Y)) ./ (W' * w);
                        sme = std(peakOf(avgs, o));
                end
                return;
            end
            switch o.kind
                case 'scalar'
                    est = mean(Y); sme = std(Y) / sqrt(n);
                otherwise
                    est = peakOf(mean(Y, 1), o);
                    s = RandStream('mt19937ar', 'Seed', opts.seed + salt);
                    idx = randi(s, n, n, opts.nBootPeak);
                    W = full(sparse(idx, repmat(1:opts.nBootPeak, n, 1), 1, n, opts.nBootPeak));
                    avgs = (W' * Y) / n;
                    sme = std(peakOf(avgs, o));
            end
        end

        function [mu, se] = meanError(X, w, dependent)
            % Sandwich SE of a scalar mean, one series per column. For
            % consecutive segments use a Bartlett HAC estimate with
            % floor(sqrt(T)) lags. Missing/rejected positions retain their
            % place in time (zero influence), rather than closing gaps.
            if nargin < 2 || isempty(w), w = ones(size(X)); end
            if nargin < 3, dependent = false; end
            if isvector(w) && size(X, 2) > 1, w = repmat(w(:), 1, size(X, 2)); end
            w = double(w); X = double(X);
            w(~isfinite(X) | ~isfinite(w) | w <= 0) = 0; X(~isfinite(X)) = 0;
            sw = sum(w, 1); n = sum(w > 0, 1);
            mu = sum(w .* X, 1) ./ sw;
            z = w .* (X - mu); z(w == 0) = 0;
            v = sum(z .^ 2, 1);
            if dependent
                lag = min(size(X, 1) - 1, max(1, floor(sqrt(size(X, 1)))));
                for k = 1:lag
                    v = v + 2 * (1 - k / (lag + 1)) * sum(z(1:end-k, :) .* z(1+k:end, :), 1);
                end
            end
            se = sqrt(max(v, 0) .* n ./ max(n - 1, 1)) ./ sw;
            mu(n < 2) = NaN; se(n < 2) = NaN;
        end

        function v = smeBoot(o, Wc, opts, salt, scale, indices)
            % SME of objective o under each outer bootstrap draw (paired:
            % Wc{c} holds the same trial counts for every candidate).
            % scale(c) rescales condition c (half-samples -> full sample).
            % Returns 1 x B, RMS over conditions.
            if nargin < 5, scale = ones(1, numel(Wc)); end
            if nargin < 6, indices = {}; end
            B = size(Wc{1}, 2); acc = zeros(1, B);
            weighted = isfield(o, 'w') && ~isempty(o.w);
            for c = 1:numel(Wc)
                X = o.X{c}; M = double(~isnan(X(:, 1))); X0 = X; X0(isnan(X0)) = 0;
                Nb = M' * Wc{c};
                if weighted
                    wt = o.w{c}; wt(~M) = 0;
                end
                if strcmp(o.kind, 'scalar') && pipecompare.utils.fieldOr(o, 'segmented', false)
                    wt = ones(size(X, 1), 1); if weighted, wt = o.w{c}; end
                    if ~isempty(indices) && ~isempty(indices{c})
                        % The order of a block resample matters to a HAC
                        % SE. Counts alone cannot reconstruct its adjacency.
                        I = indices{c}; Z = reshape(X(I), size(I));
                        wb = reshape(wt(I), size(I));
                    else
                        Z = repmat(X, 1, B); wb = Wc{c} .* wt;
                    end
                    [~, sm] = pipecompare.eval.Measure.meanError(Z, wb, true); sm2 = sm .^ 2;
                elseif weighted && strcmp(o.kind, 'scalar')
                    % the weighted SME (see above) of each resample: a trial
                    % drawn k times counts k times with its weight
                    s0 = wt' * Wc{c}; s1 = (wt .* X0)' * Wc{c};
                    q0 = (wt .^ 2)' * Wc{c}; q1 = (wt .^ 2 .* X0)' * Wc{c}; q2 = (wt .^ 2 .* X0 .^ 2)' * Wc{c};
                    mu = s1 ./ s0;
                    sm2 = max(q2 - 2 * mu .* q1 + mu .^ 2 .* q0, 0) .* Nb ./ (Nb - 1) ./ s0 .^ 2;
                elseif strcmp(o.kind, 'scalar')
                    s1 = X0' * Wc{c}; s2 = (X0 .^ 2)' * Wc{c};
                    sm2 = max((s2 - s1 .^ 2 ./ Nb) ./ (Nb - 1), 0) ./ Nb;
                else
                    % nested bootstrap: bSME within each outer draw
                    sm2 = zeros(1, B);
                    s = RandStream('mt19937ar', 'Seed', opts.seed + salt + c);
                    Bin = opts.nBootPeakInner; keep = find(M);
                    for b = 1:B
                        w = Wc{c}(keep, b); n = sum(w);
                        if n < 2, sm2(b) = Inf; continue; end
                        nz = w > 0; kk = keep(nz); w = w(nz);
                        edges = [0; cumsum(w(:)) / n]; edges(end) = 1;
                        pick = discretize(rand(s, n, Bin), edges);   % index into kk
                        Win = full(sparse(pick, repmat(1:Bin, n, 1), 1, numel(kk), Bin));
                        if weighted
                            avgs = (Win' * (wt(kk) .* X(kk, :))) ./ (Win' * wt(kk));
                        else
                            avgs = (Win' * X(kk, :)) / n;
                        end
                        sm2(b) = var(peakOf(avgs, o));
                    end
                end
                sm2(:, Nb < 2) = Inf;
                acc = acc + sm2 * scale(c) ^ 2;
            end
            v = sqrt(acc / numel(Wc));
        end

        function tw = storedWeights(EEG)
            % The weights a weight_trials step stored ([] when none):
            % urevent (of each epoch's time-locking event), weight, events
            % (the event types epochs are locked to).
            tw = [];
            if isfield(EEG, 'etc') && isstruct(EEG.etc) && isfield(EEG.etc, 'pipecompare') && ...
                    isfield(EEG.etc.pipecompare, 'trialWeights')
                tw = EEG.etc.pipecompare.trialWeights;
            end
        end

        function A = weightedAverages(EEG, contract)
            T = pipecompare.eval.Measure.trials(EEG, contract, true);
            w = pipecompare.eval.Measure.epochWeights(EEG);
            nC = numel(contract.conditions);
            A = struct('conditions', {{contract.conditions.name}}, 'labels', {{EEG.chanlocs.labels}}, ...
                'times', 1000 * (EEG.xmin + (0:EEG.pnts-1) / EEG.srate), 'data', nan(EEG.nbchan, EEG.pnts, nC), ...
                'trials', zeros(1, nC), 'sumWeights', zeros(1, nC));
            for c = 1:nC
                ep = T.epoch(T.cond == c); wc = w(ep);
                A.trials(c) = numel(ep); A.sumWeights(c) = sum(wc);
                if isempty(ep), continue; end
                A.data(:, :, c) = sum(double(EEG.data(:, :, ep)) .* reshape(wc, 1, 1, []), 3) / sum(wc);
            end
        end

        function w = epochWeights(EEG)
            % The weight of each epoch (1 when the epochs are not weighted,
            % or for an epoch not locked to an event of a condition).
            w = ones(1, EEG.trials);
            tw = pipecompare.eval.Measure.storedWeights(EEG);
            if isempty(tw) || EEG.trials < 2, return; end
            lock = pipecompare.eval.Measure.lockingEvents(EEG, tw.events);
            for k = find(lock > 0)
                at = find(tw.urevent == double(EEG.event(lock(k)).urevent), 1);
                if ~isempty(at), w(k) = tw.weight(at); end
            end
        end

        function lock = lockingEvents(EEG, codes)
            % For each epoch, the index in EEG.event of the event it is
            % time-locked to: the first event of one of the codes within
            % half a sample of time 0 (0 = none). The one definition used
            % for scoring, for the trial rule and for the signal check.
            halfSample = 1000 / EEG.srate / 2 + 1e-6;   % eventlatency is in ms
            lock = zeros(1, EEG.trials);
            for k = 1:EEG.trials
                ep = EEG.epoch(k);
                lat = cellify(ep.eventlatency); ev = cellify(ep.event); ty = ep.eventtype;
                if ischar(ty) || isstring(ty), ty = {char(ty)}; elseif ~iscell(ty), ty = num2cell(ty); end
                for q = 1:numel(lat)
                    if abs(double(lat{q})) > halfSample, continue; end
                    if ~any(strcmp(strtrim(char(string(ty{q}))), codes)), continue; end
                    lock(k) = double(ev{q});
                    break;
                end
            end
        end

        function T = trials(EEG, contract, idsOnly)
            % One row per epoch: urevent id, condition index, the epoch's
            % index in EEG, per-objective data (idsOnly: identities and
            % conditions only).
            if nargin < 3, idsOnly = false; end
            conds = contract.conditions;
            codes = contract.allEvents();
            win = contract.epoch;
            if EEG.trials == 1 && (~isfield(EEG, 'epoch') || isempty(EEG.epoch))
                [~, EEG] = evalc('pop_epoch(EEG, codes, win, ''epochinfo'', ''yes'')');
            end
            % (EEGLAB turns a single epoch into continuous data, so its times are lost)
            assert(EEG.trials > 1, 'PipeCompare:Contract', ...
                'There is only one epoch; each condition needs at least two trials.');
            tol = 1.5 / EEG.srate;
            assert(EEG.xmin <= win(1) + tol && EEG.xmax >= win(2) - tol, ...
                'PipeCompare:Measure', 'Epochs [%g %g] s do not cover the contract epoch [%g %g] s.', EEG.xmin, EEG.xmax, win(1), win(2));
            times = EEG.xmin + (0:EEG.pnts-1) / EEG.srate;
            labels = lower({EEG.chanlocs.labels});
            nT = EEG.trials;
            id = nan(nT, 1); cond = zeros(nT, 1);
            eligible = [];
            if isfield(EEG, 'etc') && isfield(EEG.etc, 'pipecompare') && isfield(EEG.etc.pipecompare, 'eligibleUrevents')
                eligible = EEG.etc.pipecompare.eligibleUrevents;
            end
            lock = pipecompare.eval.Measure.lockingEvents(EEG, codes);
            for k = find(lock > 0)
                e = EEG.event(lock(k));
                t = strtrim(char(string(e.type)));
                assert(isfield(e, 'urevent') && ~isempty(e.urevent), 'PipeCompare:Urevent', ...
                    'Time-locking events have no urevent index.');
                if ~isempty(eligible) && ~ismember(double(e.urevent), eligible), continue; end
                id(k) = double(e.urevent);
                cond(k) = find(arrayfun(@(cd) any(strcmp(cd.events, t)), conds), 1);
            end
            keep = cond > 0;
            T = struct('id', id(keep), 'cond', cond(keep), 'epoch', find(keep), 'data', {{}}, 'times', {{}}, ...
                'labels', {labels});
            if idsOnly, return; end
            data = double(EEG.data(:, :, keep));
            bl = contract.baseline;
            if ~isempty(bl)
                bsel = times >= bl(1) - 1e-9 & times <= bl(2) + 1e-9;
                assert(any(bsel), 'PipeCompare:Measure', 'No samples in the baseline window.');
                data = data - mean(data(:, bsel, :), 2);
            end
            if strcmp(contract.analysis, 'bandpower')
                w = times >= win(1) - 1e-9 & times <= win(2) + 1e-9;
                for b = 1:numel(contract.bands)
                    band = contract.bands(b);
                    T.data{b} = pipecompare.eval.Measure.bandScores(data(roiIndex(band.roi, labels), w, :), EEG.srate, band.freq);
                    T.times{b} = [];
                end
                return;
            end
            for j = 1:numel(contract.components)
                comp = contract.components(j);
                roi = roiIndex(comp.roi, labels);
                w = times >= comp.window(1) - 1e-9 & times <= comp.window(2) + 1e-9;
                if contract.isLateral(j)
                    % contralateral minus ipsilateral: the first ROI electrode
                    % minus the second, sign-flipped where the second is the
                    % contralateral one of the trial's condition
                    [~, ci] = ismember(lower(comp.contra), lower(comp.roi));
                    sgn = 3 - 2 * ci(T.cond);
                    wave = permute(data(roi(1), w, :) - data(roi(2), w, :), [3 2 1]) .* sgn(:);
                else
                    wave = permute(mean(data(roi, w, :), 1), [3 2 1]);   % trials x samples
                end
                if size(wave, 2) ~= sum(w), wave = reshape(wave, [], sum(w)); end
                if strcmp(comp.measure, 'mean')
                    T.data{j} = mean(wave, 2);
                else
                    T.data{j} = wave;
                end
                T.times{j} = 1000 * times(w);
            end
        end
    end
end

% ---------------------------------------------------------------------
function c = cellify(x)
if iscell(x), c = x; else, c = num2cell(x); end
end

function idx = roiIndex(roi, labels)
[ok, idx] = ismember(lower(roi), labels);
if ~all(ok)
    error('PipeCompare:RoiMissing', 'ROI channel(s) missing in this candidate: %s', strjoin(roi(~ok), ', '));
end
end

function [kind, pol] = objKind(contract, k)
pol = '';
if strcmp(contract.analysis, 'bandpower'), kind = 'scalar'; return; end
comp = contract.components(k);
pol = comp.polarity;
switch comp.measure
    case 'mean', kind = 'scalar';
    otherwise, kind = comp.measure;
end
end

function c = countBelow(d, e2, x)
% The number of eigenvalues below x of the symmetric tridiagonal matrix with
% diagonal d and squared off-diagonal e2 (Sturm sequence).
c = 0; q = d(1) - x;
if q < 0, c = 1; end
for i = 2:numel(d)
    if q == 0, q = realmin; end
    q = d(i) - x - e2(i-1) / q;
    if q < 0, c = c + 1; end
end
end

function v = peakOf(avgs, o)
% avgs: draws x samples. Peak of each averaged waveform.
if strcmp(o.polarity, 'negative'), [a, i] = min(avgs, [], 2); else, [a, i] = max(avgs, [], 2); end
if strcmp(o.kind, 'peakLatency'), v = reshape(o.times(i), [], 1); else, v = a; end
end

function H = sampleHits(E, isEye, isL, isR, contract)
% The blinks and eye movements of the continuous data E, found sample by
% sample (the trial rule compares the trials with each other, so when
% most trials have a blink the median trial has one, and few are found):
% where an eye channel, or the left minus the right channel, minus its
% moving median over 1 s (slow drifts out) is beyond 5 robust SDs (median
% and MAD over the samples: blinks fill far fewer samples than trials),
% or in the active state of a two-state HMM fitted to it
% (pipecompare.eval.Hmm.events: weak blinks too), with 50 ms on either
% side. Epoched like the data: one row per kind (blinks; eye movements
% when there is a pair). [] without eye channels.
H = [];
if ~any(isEye), return; end
rows = {double(E.data(isEye, :))};
if any(isL) && any(isR), rows{2} = double(E.data(isL, :)) - double(E.data(isR, :)); end
w = round(0.05 * E.srate); M = zeros(numel(rows), E.pnts);
for k = 1:numel(rows)
    Y = rows{k} - movmedian(rows{k}, round(E.srate), 2);
    Y = abs(Y - median(Y, 2));
    on = any(Y > 5 * max(1.4826 * median(Y, 2), eps), 1) | pipecompare.eval.Hmm.events(rows{k}, E.srate);
    M(k, :) = conv(double(on), ones(1, 2 * w + 1), 'same') > 0;
end
H = E; H.data = single(M); H.nbchan = size(M, 1);
H.chanlocs = struct('labels', arrayfun(@(k) sprintf('hit%d', k), 1:H.nbchan, 'UniformOutput', false));
[~, H] = evalc('pop_epoch(H, contract.allEvents(), contract.epoch, ''epochinfo'', ''yes'')');
H.data = H.data > 0.5;
end

function [l, r] = horizontalPair(chanlocs, eyes, roi)
% The left and the right channel whose difference shows horizontal eye
% movements: of the EOG channels, the leftmost and the rightmost by
% position (EEGLAB's Y axis points to the left ear) when they lie on
% either side, else those whose names end in L and R (e.g. EYEL and EYER);
% without EOG channels, or when those give no pair, F7 and F8. Empty when
% there is none or one of them is measured.
l = []; r = [];
eog = [];
if numel(eyes) >= 2, eog = eyes(pipecompare.simple.Presets.isNonEeg(chanlocs(eyes))); end
if numel(eog) >= 2
    y = nan(1, numel(eog));
    for k = 1:numel(eog)
        c = chanlocs(eog(k));
        if isfield(c, 'Y') && ~isempty(c.X) && ~isempty(c.Y) && ~isempty(c.Z)
            y(k) = c.Y / max(norm([c.X c.Y c.Z]), eps);
        end
    end
    [ym, a] = max(y); [yn, b] = min(y);
    if ym > 0.2 && yn < -0.2, l = eog(a); r = eog(b); return; end
    labs = upper(strtrim({chanlocs(eog).labels}));
    a = find(~cellfun(@isempty, regexp(labs, 'L$', 'once'))); b = find(~cellfun(@isempty, regexp(labs, 'R$', 'once')));
    if isscalar(a) && isscalar(b), l = eog(a); r = eog(b); return; end
end
labs = upper({chanlocs.labels});
a = find(strcmp(labs, 'F7')); b = find(strcmp(labs, 'F8'));
if isscalar(a) && isscalar(b) && ~any(ismember([a b], roi)), l = a; r = b; end
end

