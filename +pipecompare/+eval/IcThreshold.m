classdef IcThreshold
    %ICTHRESHOLD The ICLabel threshold of IC removal chosen from the data
    %   (icremove with threshold = 'auto').
    %
    %   [t, info] = pipecompare.eval.IcThreshold.choose(EEG, q, ctx)
    %
    %   q: per component, the largest ICLabel probability among the artifact
    %   classes of the step (below 1: pop_icflag flags a class probability p
    %   when t < p < 1). A component is removed when q > t, so lowering t
    %   only adds components: the removed sets are nested and there are only
    %   as many different ones as different values of q. Every one of them
    %   that removes only components with q > MinProbability (more likely
    %   that artifact than anything else) is compared, and removing nothing:
    %   the exact best threshold, not the best of a few thresholds tried.
    %
    %   Removing components is linear (the data minus the removed
    %   components' projections, as pop_subcomp computes it), and the steps
    %   up to the epoch rejection that follows (ctx.after: epoching and
    %   baseline removal) commute with it, so each set is applied to the
    %   epoched data directly. For each set, the rejection that follows is
    %   applied as it would run: a fixed limit, or with uv 'auto' the best of
    %   all its limits (pipecompare.eval.Threshold). The objective is the
    %   ranking's (gain-corrected SME), the gain and the signal check coming
    %   from the known-signal copy (ctx.signal) with the same components
    %   removed; a set that fails the signal check or keeps too few trials is
    %   not a candidate. Among the sets shown equivalent to the best one
    %   (paired bootstrap, simultaneous intervals: the whole interval of the
    %   difference within equivalenceMargin of the best objective), the one
    %   removing the fewest components is chosen: no component is removed
    %   unless keeping it would make the measure measurably less precise. A
    %   set only not distinguished from the best (too few trials to tell) is
    %   not enough. t is the roundest number among the thresholds that remove
    %   that set.
    %
    %   The signal check of each set is that of all epochs (the rejected
    %   ones change it only where epochs overlap); the candidate's own check
    %   after the whole pipeline is the one the ranking uses.
    %   Needs per-trial scores: mean amplitude or band power.

    %
    %   [comps, info] = pipecompare.eval.IcThreshold.chooseBySme(EEG, cand, ctx)
    %
    %   (icremove with threshold = 'sme') The components to remove chosen by
    %   their effect on the objective instead of their ICLabel probability.
    %   cand: the components that may be removed (those whose likeliest
    %   ICLabel class is not Brain: a component ICLabel takes for brain
    %   activity is never removed). Starting from removing nothing, the
    %   component whose removal gives the lowest objective is added, one at
    %   a time (forward selection), each set evaluated exactly as above:
    %   the data without those components (linear, so computed from the
    %   activations directly), the rejection that follows at its best
    %   limit, the gain and the signal check from the known-signal copy (a
    %   set that fails it cannot be chosen). This gives nested sets ordered
    %   by what each component does to the measure, not by what it looks
    %   like. The search stops after MaxRemoved components, when no further
    %   set passes, or when the objective has risen more than twice
    %   equivalenceMargin above the best on the path. Among the sets on the
    %   path, the one with the fewest components among those shown
    %   equivalent to the best is chosen, as above. info.order: the
    %   components in the order they were added.
    %   Components that carry time-locked activity into a measured window
    %   are not candidates (info.protected): removing one alone would move
    %   a condition's mean score by more than EvokedShare times its SME,
    %   with |t| > EvokedT against the component's own trial-to-trial
    %   spread. The SME counts real trial-to-trial brain activity as noise,
    %   so without this a brain component that ICLabel does not recognise
    %   could be removed for lowering it. Components whose likeliest
    %   ICLabel class is Eye are exempt (blinks can be time-locked).
    %
    %   Both: info.curves holds the scores of every set (fewest components
    %   first), with which Rank.crossfit makes the choice again on half of
    %   the trials.

    properties (Constant)
        MinProbability = 0.5   % lowest threshold compared
        MaxRemoved = 30        % 'sme': at most this many components removed
        EvokedT = 3            % 'sme': a component this clearly time-locked ...
        EvokedShare = 0.5      % ... and moving a condition's mean by this many SMEs is kept
    end

    methods (Static)
        function [comps, info] = chooseBySme(EEG, cand, ctx)
            contract = ctx.contract;
            ref = pipecompare.utils.fieldOr(ctx, 'ref', []);
            if isempty(ref), ref = pipecompare.eval.Measure.reference(EEG, contract); end
            opts = pipecompare.utils.withDefaults(pipecompare.utils.fieldOr(ctx, 'rank', struct()), ...
                pipecompare.eval.Rank.defaults());
            assert(pipecompare.eval.Threshold.applies(contract), 'PipeCompare:IcThreshold', ...
                ['Components chosen by their effect on the SME (threshold = ''sme'') need mean amplitude or band power ', ...
                 '(peak measures have no score per trial); use ICLabel thresholds instead (e.g. 0.7 | 0.8 | 0.9).']);
            [before, rej] = splitAfter(pipecompare.utils.fieldOr(ctx, 'after', {}));
            Ep = runBefore(EEG, before, ctx);
            sig = pipecompare.utils.fieldOr(ctx, 'signal', []);
            Sp = []; PS = [];
            if ~isempty(sig) && ~isempty(sig.S), Sp = runBefore(sig.S, before, ctx); PS = projection(EEG, Sp.data); end
            P = projection(EEG, Ep.data);
            chans = [];
            if ~isempty(rej), chans = pipecompare.run.Steps.testedChannels(Ep, rej.params); end
            % components that carry time-locked activity into a measured
            % window are not candidates (the SME counts real trial-to-trial
            % brain activity as noise; ICLabel's Eye components excepted)
            cand = cand(:)';
            ev = evokedGuard(EEG, Ep, P, contract, cand);
            protected = cand(ev.protect);
            cand = cand(~ev.protect);
            % the path: nothing removed, then one component more at a time
            [f, u, why, X, M, C] = scoreSet(Ep, Sp, rej, chans, contract, ref, opts, sig);
            sets = {zeros(1, 0)}; why = {why}; X = {X}; M = {M}; C = {C};
            cur = zeros(1, 0); base = P.PX; baseS = []; if ~isempty(PS), baseS = PS.PX; end
            for step = 1:min(numel(cand), pipecompare.eval.IcThreshold.MaxRemoved)
                left = setdiff(cand, cur, 'stable');
                fb = Inf; kb = 0; out = {};
                for c = left
                    Ek = withData(Ep, P, base - P.A(:, c) * P.S(c, :));
                    Sk = [];
                    if ~isempty(PS), Sk = withData(Sp, PS, baseS - PS.A(:, c) * PS.S(c, :)); end
                    [fc, uc, wc, Xc, Mc, Cc] = scoreSet(Ek, Sk, rej, chans, contract, ref, opts, sig);
                    if fc < fb, fb = fc; kb = c; out = {uc, wc, Xc, Mc, Cc}; end
                end
                if kb == 0, break; end                       % no further set passes
                cur(end+1) = kb; %#ok<AGROW>
                base = base - P.A(:, kb) * P.S(kb, :);
                if ~isempty(PS), baseS = baseS - PS.A(:, kb) * PS.S(kb, :); end
                sets{end+1} = cur; f(end+1) = fb; u(end+1) = out{1}; %#ok<AGROW>
                why{end+1} = out{2}; X{end+1} = out{3}; M{end+1} = out{4}; C{end+1} = out{5}; %#ok<AGROW>
                if fb > (1 + 2 * opts.equivalenceMargin) * min(f), break; end   % well past the best
            end
            K = numel(sets);
            info = struct('sets', K, 'considered', numel(cand), 'feasible', sum(isfinite(f)), ...
                'notDistinguished', 0, 'equallyGood', 0, 'objectiveBest', NaN, 'objectiveChosen', NaN, ...
                'removedBest', 0, 'removedChosen', 0, 'objectives', f, 'order', cur, 'why', {why}, 'note', '', ...
                'protected', protected, 'evoked', ev);
            F = find(isfinite(f));
            if isempty(F)
                chosen = 1; best = 1;
                info.note = sprintf('removing no component does not meet the limits (%s); none removed', why{1});
            else
                best = F(find(f(F) == min(f(F)), 1));     % ties: fewer components
                Bq = zeros(numel(F), opts.nBoot);
                for i = 1:numel(F)
                    Bq(i, :) = pipecompare.eval.Threshold.bootstrap(X{F(i)}, M{F(i)}, u(F(i)), opts, ref, C{F(i)});
                end
                [keep, equiv] = pipecompare.eval.Rank.confidenceSet(f(F)', Bq, find(F == best), ...
                    opts.alpha, opts.equivalenceMargin);
                chosen = min(F(equiv));
                info.notDistinguished = sum(keep); info.equallyGood = sum(equiv);
                info.objectiveBest = f(best); info.objectiveChosen = f(chosen);
            end
            comps = sort(sets{chosen});
            info.curves = struct('options', {C}, 'ok', isfinite(f), 'chosen', chosen);   % for Rank.crossfit
            info.removedBest = numel(sets{best}); info.removedChosen = numel(comps);
            pipecompare.utils.log(['ICA components chosen by their effect on the SME: %d removed (best %d; %d of %d ', ...
                'components could be removed, %d kept because they carry time-locked activity into the measured ', ...
                'window, %d sets on the path, %d equally good).'], numel(comps), numel(sets{best}), ...
                numel(cand), size(EEG.icaweights, 1), numel(protected), K, info.equallyGood);
        end

        function [t, info] = choose(EEG, q, ctx)
            contract = ctx.contract;
            ref = pipecompare.utils.fieldOr(ctx, 'ref', []);
            if isempty(ref), ref = pipecompare.eval.Measure.reference(EEG, contract); end
            opts = pipecompare.utils.withDefaults(pipecompare.utils.fieldOr(ctx, 'rank', struct()), ...
                pipecompare.eval.Rank.defaults());
            assert(pipecompare.eval.Threshold.applies(contract), 'PipeCompare:IcThreshold', ...
                ['An ICLabel threshold chosen from the data (threshold = ''auto'') needs mean amplitude or band power ', ...
                 '(peak measures have no score per trial); set the thresholds instead (e.g. 0.7 | 0.8 | 0.9).']);
            [before, rej] = splitAfter(pipecompare.utils.fieldOr(ctx, 'after', {}));
            Ep = runBefore(EEG, before, ctx);
            sig = pipecompare.utils.fieldOr(ctx, 'signal', []);
            Sp = [];
            if ~isempty(sig) && ~isempty(sig.S), Sp = runBefore(sig.S, before, ctx); end
            % the removed sets: none, then down the distinct values of q
            qr = q(:)';                          % as ICLabel stores them (pop_icflag compares these)
            q = double(qr);
            c = sort(unique(q(q > pipecompare.eval.IcThreshold.MinProbability)), 'descend');
            K = numel(c) + 1;                    % set 1 removes nothing, set k those with q >= c(k-1)
            hi = [1 c]; lo = [c max([0 q(q <= pipecompare.eval.IcThreshold.MinProbability)])];
            R = [{false(size(qr))} arrayfun(@(x) qr >= x, c, 'UniformOutput', false)];   % a threshold in [lo(k), hi(k)) removes R{k} (q > t)
            chans = [];
            if ~isempty(rej), chans = pipecompare.run.Steps.testedChannels(Ep, rej.params); end
            f = nan(1, K); u = nan(1, K); why = repmat({''}, 1, K); X = cell(1, K); M = cell(1, K); C = cell(1, K);
            for k = 1:K
                comps = find(R{k});
                Ek = Ep; Ek.data = removeComponents(EEG, Ep.data, comps);
                Sk = [];
                if ~isempty(Sp), Sk = Sp; Sk.data = removeComponents(EEG, Sp.data, comps); end
                [f(k), u(k), why{k}, X{k}, M{k}, C{k}] = scoreSet(Ek, Sk, rej, chans, contract, ref, opts, sig);
            end
            % per set (fewest components first): its objective at its best
            % limit (NaN: not a candidate, why says why) and the number removed
            info = struct('thresholdBest', NaN, 'thresholdChosen', NaN, 'sets', K, 'feasible', sum(isfinite(f)), ...
                'notDistinguished', 0, 'equallyGood', 0, 'objectiveBest', NaN, 'objectiveChosen', NaN, 'removedBest', 0, 'removedChosen', 0, ...
                'objectives', f, 'removed', cellfun(@sum, R), 'why', {why}, 'note', '');
            F = find(isfinite(f));
            if isempty(F)
                % no set meets the limits: remove the fewest (the ranking
                % then reports what is wrong)
                chosen = 1; best = 1;
                info.note = sprintf('no threshold meets the limits (%s); the fewest components removed', why{1});
            else
                best = F(find(f(F) == min(f(F)), 1));     % ties: fewer components
                Bq = zeros(numel(F), opts.nBoot);
                for i = 1:numel(F)
                    Bq(i, :) = pipecompare.eval.Threshold.bootstrap(X{F(i)}, M{F(i)}, u(F(i)), opts, ref, C{F(i)});
                end
                [keep, equiv] = pipecompare.eval.Rank.confidenceSet(f(F)', Bq, find(F == best), ...
                    opts.alpha, opts.equivalenceMargin);
                chosen = min(F(equiv));
                info.notDistinguished = sum(keep); info.equallyGood = sum(equiv);
                info.objectiveBest = f(best); info.objectiveChosen = f(chosen);
            end
            info.curves = struct('options', {C}, 'ok', isfinite(f), 'chosen', chosen);   % for Rank.crossfit
            t = roundIn(lo(chosen), hi(chosen), qr, R{chosen});
            info.thresholdBest = roundIn(lo(best), hi(best), qr, R{best}); info.thresholdChosen = t;
            info.removedBest = sum(R{best}); info.removedChosen = sum(R{chosen});
            pipecompare.utils.log(['ICLabel threshold chosen from the data: %g (best %g; %d different component sets, ', ...
                '%d equally good; %d of %d components removed).'], t, info.thresholdBest, K, info.equallyGood, ...
                info.removedChosen, numel(q));
        end
    end
end

function ev = evokedGuard(EEG, Ep, P, contract, cand)
% For each candidate component and each condition: what removing it alone
% does to the condition's mean score of each mean-amplitude component (the
% component's part of the score, linear: its ROI weight times its
% activation's window mean, baseline removed), as a t value against its
% own trial-to-trial spread and as a share of the score's standard error
% (SME). Protected (not removed by 'sme'): |t| > EvokedT and shift >
% EvokedShare SME in some condition, i.e. clearly time-locked and large
% enough to change the measure; components whose likeliest ICLabel class
% is Eye are not protected (blinks can be time-locked to the stimulus).
% Band power and contralateral minus ipsilateral components are not
% checked.
nK = numel(cand);
ev = struct('protect', false(1, nK), 't', zeros(1, nK), 'shift', zeros(1, nK));
if ~strcmp(contract.analysis, 'erp') || isempty(contract.components), return; end
T = pipecompare.eval.Measure.trials(Ep, contract);
times = Ep.xmin + (0:Ep.pnts-1) / Ep.srate;
labs = lower({Ep.chanlocs.labels});
[~, rowOf] = ismember(1:Ep.nbchan, P.ch);                     % data channel -> row of A
isEye = false(1, nK);
if isfield(EEG.etc, 'ic_classification') && isfield(EEG.etc.ic_classification, 'ICLabel')
    L = EEG.etc.ic_classification.ICLabel;
    [~, top] = max(L.classifications, [], 2);
    isEye = strcmp(cellstr(L.classes(top(cand))), 'Eye'); isEye = isEye(:)';
end
bl = [];
if ~isempty(contract.baseline)
    bl = times >= contract.baseline(1) - 1e-9 & times <= contract.baseline(2) + 1e-9;
end
for j = 1:numel(contract.components)
    comp = contract.components(j);
    if contract.isLateral(j) || ~strcmp(comp.measure, 'mean'), continue; end
    roi = find(ismember(labs, lower(comp.roi)));
    roi = roi(rowOf(roi) > 0);
    if isempty(roi), continue; end
    w = times >= comp.window(1) - 1e-9 & times <= comp.window(2) + 1e-9;
    sc = T.data{j};                                              % the scores (trials)
    for i = 1:nK
        k = cand(i);
        Sk = reshape(P.S(k, :), Ep.pnts, Ep.trials);
        a = mean(Sk(w, :), 1);
        if ~isempty(bl) && any(bl), a = a - mean(Sk(bl, :), 1); end
        d = mean(P.A(rowOf(roi), k)) * a(T.epoch);               % the component's part of each trial's score
        for c = unique(T.cond(:))'
            r = T.cond == c; n = sum(r);
            if n < 3, continue; end
            dm = mean(d(r)); sd = std(d(r)); sme = std(sc(r)) / sqrt(n);
            t = dm / max(sd / sqrt(n), eps); sh = abs(dm) / max(sme, eps);
            if abs(t) > abs(ev.t(i)), ev.t(i) = t; end
            ev.shift(i) = max(ev.shift(i), sh);
            if abs(t) > pipecompare.eval.IcThreshold.EvokedT && sh > pipecompare.eval.IcThreshold.EvokedShare && ~isEye(i)
                ev.protect(i) = true;
            end
        end
    end
end
end

function P = projection(EEG, X)
% The data X (channels x samples x epochs) as pop_subcomp rebuilds them
% from all components (PX, on the ICA channels ch), and the activations S:
% removing the set R gives PX - A(:, R) * S(R, :) on those channels.
W = double(EEG.icaweights) * double(EEG.icasphere);
A = double(EEG.icawinv); if isempty(A), A = pinv(W); end
sz = size(X); X2 = reshape(X, sz(1), []);
ch = EEG.icachansind;
S = W * double(X2(ch, :));
P = struct('PX', A * S, 'A', A, 'S', S, 'ch', ch, 'sz', sz);
end

function E = withData(E, P, Y)
% E with its ICA channels replaced by Y (as P.PX), the others unchanged
X2 = reshape(E.data, P.sz(1), []);
X2(P.ch, :) = Y;
E.data = reshape(X2, P.sz);
end

function [f, u, why, X, M, C] = scoreSet(Ek, Sk, rej, chans, contract, ref, opts, sig)
% One set of components removed (Ek: the data without them; Sk: the
% known-signal copy without them, [] = no signal check): its objective at
% the best limit of the rejection that follows (NaN when it fails the
% signal check or keeps too few trials; why says why), that limit, and
% the scores and limit statistics of Threshold.curve for the bootstrap
% (C: aligned with the reference trials, for Rank.crossfit).
f = NaN; u = NaN; why = ''; g = [];
if ~isempty(Sk)
    r = pipecompare.eval.Injection.compare(Sk, contract, sig.truth, {});
    g = double(r.gain(:)');
    w = pipecompare.eval.Rank.signalReasons(r, opts);
    if ~all(g > 0), w{end+1} = 'signal lost or inverted'; end
    if ~isempty(w), why = strjoin(w, '; '); end
end
if isempty(rej)
    m = zeros(1, Ek.trials); v = 0;              % nothing rejected
else
    m = pipecompare.run.Steps.limitStat(Ek, chans, rej.params);
    v = [];                                     % every limit (uv 'auto')
    if ~ischar(rej.params.uv), v = rej.params.uv; end
end
[v, ok, obj, X, M, C] = pipecompare.eval.Threshold.curve(Ek, m, contract, ref, opts, g, v);
if ~any(ok)
    if isempty(why), why = 'too few trials kept'; end
    return;
end
F = find(ok);
b = F(find(obj(F) == min(obj(F)), 1, 'last'));
u = v(b);
if isempty(why), f = obj(b); end
end

function [before, rej] = splitAfter(after)
% the steps between IC removal and the epoch rejection, and that rejection
% ([] when none follows)
rej = [];
k = find(cellfun(@(in) strcmp(in.type, 'reject_threshold'), after), 1);
if isempty(k), before = after; else, before = after(1:k-1); rej = after{k}; end
for i = 1:numel(before)
    assert(any(strcmp(before{i}.type, {'epoch', 'baseline'})), 'PipeCompare:IcThreshold', ...
        'An ICLabel threshold chosen from the data needs only epoching and baseline removal before the epoch rejection (found %s).', ...
        before{i}.type);
end
end

function E = runBefore(E, before, ctx)
% the steps between IC removal and the rejection, then epochs (as the
% measure takes them) when the data are still continuous
for i = 1:numel(before)
    [~, E] = evalc('pipecompare.run.Steps.run(before{i}, E, ctx)');
end
if E.trials == 1
    [~, E] = evalc('pop_epoch(E, ctx.contract.allEvents(), ctx.contract.epoch, ''epochinfo'', ''yes'')');
end
end

function X = removeComponents(EEG, X, comps)
% the data without components comps, as pop_subcomp computes it: the
% back-projection of the components kept (nothing changes when none is
% removed: pop_subcomp is not run)
if isempty(comps), return; end
W = double(EEG.icaweights) * double(EEG.icasphere);
A = double(EEG.icawinv); if isempty(A), A = pinv(W); end
keep = setdiff(1:size(W, 1), comps);
sz = size(X); X2 = reshape(X, sz(1), []);
ch = EEG.icachansind;
X2(ch, :) = A(:, keep) * (W(keep, :) * double(X2(ch, :)));
X = reshape(X2, sz);
end

function t = roundIn(lo, hi, q, target)
% the roundest threshold in [lo, hi) (hi = 1 included: nothing has q > 1)
% that removes the components target, checked with the comparison
% pop_icflag makes (q > t)
for d = 0:12
    t = ceil(lo * 10 ^ d) / 10 ^ d;
    if (t < hi || t == 1) && isequal(q > t, target), return; end
end
t = lo;
end
