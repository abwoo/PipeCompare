classdef Threshold
    %THRESHOLD The amplitude limit of epoch rejection chosen from the data
    %   (reject_threshold with uv = 'auto').
    %
    %   [u, info] = pipecompare.eval.Threshold.choose(EEG, m, ctx)
    %
    %   m: per epoch, the value the limit is compared with: the largest
    %   absolute value over the tested channels and the whole epoch, or
    %   with method peaktopeak the largest moving-window peak-to-peak
    %   value. An epoch is rejected when m > u. Raising u only adds epochs,
    %   so the kept sets are nested and there are only as many different
    %   ones as different values of m: every limit between two neighbouring
    %   values keeps the same epochs. The objective of the ranking (the SME
    %   of each condition's mean score, RMS over measures and conditions, or
    %   the one measure named by the objective option) is computed for
    %   every one of these sets: the exact best limit, not the best of a few
    %   limits tried. The signal gain is left out: it is the same at every
    %   limit (the known signal is the same in every epoch).
    %
    %   Only limits that keep at least minTrials trials and minRetention of
    %   the trials of every condition (the ranking's own limits) are
    %   considered. Among them, the limits shown equivalent to the best one
    %   (the paired bootstrap and simultaneous intervals of the ranking, over
    %   at most MaxCompared limits spread over the range, the best included:
    %   the whole interval of the difference within equivalenceMargin of the
    %   best objective) are equally good, and the highest of them is chosen:
    %   no epoch is rejected unless keeping it would make the measure
    %   measurably less precise. A limit only not distinguished from the
    %   best (too few trials to tell) is not enough: the lenient limits have
    %   the widest intervals, so that set reaches limits that reject almost
    %   nothing even when rejecting clearly helps. u is the roundest number
    %   among the limits that keep the same epochs.
    %
    %   The scores are those of the data at this step (the steps after it
    %   are not run): exact when nothing after it changes the scores other
    %   than the baseline, which the scores already have removed. Epochs
    %   weighted after it (weight_trials) are not modelled: the limit is
    %   the one for equal weights.
    %   Needs per-trial scores: mean amplitude or band power.

    properties (Constant)
        MaxCompared = 60   % limits compared by the bootstrap (spread over the range, the best included)
    end

    methods (Static)
        function [u, info] = choose(EEG, m, ctx)
            contract = ctx.contract;
            % the trials of the starting data (minRetention is a share of them);
            % without them (a preview), those at this step
            ref = pipecompare.utils.fieldOr(ctx, 'ref', []);
            if isempty(ref), ref = pipecompare.eval.Measure.reference(EEG, contract,pipecompare.utils.fieldOr(ctx,'rank',struct())); end
            opts = pipecompare.utils.withDefaults(pipecompare.utils.fieldOr(ctx, 'rank', struct()), ...
                pipecompare.eval.Rank.defaults());
            assert(pipecompare.eval.Threshold.applies(contract), 'PipeCompare:Threshold', ...
                ['A rejection limit chosen from the data (uv = ''auto'') needs mean amplitude or band power ', ...
                 '(peak measures have no score per trial); set the limits instead (e.g. 75 | 100 | 150).']);
            [v, ok, obj, X, M, C] = pipecompare.eval.Threshold.curve(EEG, m, contract, ref, opts, [], []);
            U = numel(v);
            % curves: the choice as Rank.crossfit makes it again on half of
            % the trials (one option: these data; its limit chosen there)
            info = struct('uvBest', NaN, 'uvChosen', NaN, 'limits', U, 'feasible', sum(ok), 'compared', 0, ...
                'notDistinguished', 0, 'equallyGood', 0, 'objectiveBest', NaN, 'objectiveChosen', NaN, 'rejectedBest', 0, 'rejectedChosen', 0, ...
                'note', '', 'curves', struct('options', {{C}}, 'ok', true, 'chosen', 1));
            nAll = numel([M{:}]);
            if ~any(ok)
                % no limit keeps enough trials: reject nothing (the ranking
                % then reports the trials missing)
                u = roundIn(v(end), 2 * v(end)); info.uvChosen = u; info.uvBest = u;
                info.note = 'no limit keeps enough trials; nothing rejected';
                return;
            end
            F = find(ok);
            best = F(find(obj(F) == min(obj(F)), 1, 'last'));   % ties: the higher limit
            pick = F;
            if numel(F) > pipecompare.eval.Threshold.MaxCompared
                pick = unique([F(round(linspace(1, numel(F), pipecompare.eval.Threshold.MaxCompared - 1))) best]);
            end
            Bq = pipecompare.eval.Threshold.bootstrap(X, M, v(pick), opts, ref, C);
            ib = find(pick == best);
            keep = pipecompare.eval.Rank.confidenceSet(obj(pick)', Bq, ib, opts.alpha);
            % the highest limit shown equivalent to the best, going up from it
            [c, nEq] = pipecompare.eval.Rank.fixedSequence(obj(pick)', Bq, ib, ib+1:numel(pick), ...
                opts.alpha, opts.equivalenceMargin);
            chosen = pick(c);
            u = roundIn(v(chosen), nextAbove(v, chosen));
            info.uvBest = roundIn(v(best), nextAbove(v, best)); info.uvChosen = u;
            info.compared = numel(pick); info.notDistinguished = sum(keep); info.equallyGood = nEq;
            info.objectiveBest = obj(best); info.objectiveChosen = obj(chosen);
            info.rejectedBest = sum([M{:}] > v(best)); info.rejectedChosen = sum([M{:}] > v(chosen));
            pipecompare.utils.log(['Rejection limit chosen from the data: %g uV (best %g uV; %d different limits, ', ...
                '%d equally good; %d/%d trials rejected).'], u, info.uvBest, U, info.equallyGood, info.rejectedChosen, nAll);
        end

        function [v, ok, obj, X, M, C] = curve(EEG, m, contract, ref, opts, g, v)
            % The objective at each limit v (empty: every value of m, one
            % limit per kept set), and whether the limit keeps enough
            % trials. X{c}: the centred scores of condition c's reference
            % trials (time order) divided by the gain g of each measure
            % (empty: 1), so that their SD is the gain-corrected one; M{c}:
            % their values of m. C: the same scores and values of m at the
            % rows of all ref.n(c) reference trials (NaN: not present at
            % this step), the measures they belong to (use), the limits
            % given (v: [] every one) and whether the gain is divided out
            % (gained): what the ranking's selection check needs to make
            % this choice again on half of the trials (Rank.crossfit).
            vIn = v;
            name = pipecompare.eval.Rank.objectiveName(opts.objective, ref);
            use = 1:numel(ref.objectives);
            if ~strcmp(name, 'composite'), use = find(strcmp(ref.objectives, name)); end
            gained = ~isempty(g);
            if isempty(g), g = ones(1, numel(ref.objectives)); end
            T = pipecompare.eval.Measure.trials(EEG, contract);
            m = double(m(:)');
            nC = numel(ref.ids); X = cell(1, nC); M = cell(1, nC);
            C = struct('X', {cell(1, nC)}, 'M', {cell(1, nC)}, 'use', use, 'v', vIn, 'gained', gained, 'loc', {cell(1, nC)});
            for c = 1:nC
                % the reference trials of condition c, in time order
                rows = find(T.cond == c);
                [tf, loc] = ismember(T.id(rows), ref.ids{c});
                rows = rows(tf); loc = loc(tf);
                [~, first] = unique(loc); rows = rows(first);     % sorted by loc: time order
                Xc = [T.data{use}] ./ double(g(use));
                X{c} = Xc(rows, :) - mean(Xc(rows, :), 1);       % centred: the SD does not change
                M{c} = m(T.epoch(rows));
                C.X{c} = nan(ref.n(c), numel(use)); C.X{c}(unique(loc), :) = X{c};
                C.M{c} = nan(ref.n(c), 1); C.M{c}(unique(loc)) = M{c};
                C.loc{c} = unique(loc);                          % the reference row of each row of X{c}
            end
            if isempty(v), v = unique([M{:}]); end   % one limit per kept set (keep m <= v(j))
            v = v(:)';
            assert(~isempty(v), 'PipeCompare:Threshold', 'No trial of the conditions is left at this step.');
            need = max(opts.minTrials, ceil(opts.minRetention * ref.n - 1e-9));
            segmented = pipecompare.utils.fieldOr(ref,'segmented',false);
            ok = true(1, numel(v)); acc = zeros(1, numel(v));
            for c = 1:nC
                K = double(M{c}(:) <= v);                     % trials x limits
                n = sum(K, 1);
                ok = ok & n >= max(need(c), 2);
                for o = 1:size(X{c}, 2)
                    if segmented
                        fullK = double(C.M{c}(:) <= v);
                        [~, se] = pipecompare.eval.Measure.meanError(repmat(C.X{c}(:,o),1,numel(v)),fullK,true);
                        acc = acc + se.^2;
                    else
                        s1 = X{c}(:, o)' * K; s2 = (X{c}(:, o) .^ 2)' * K;
                        acc = acc + max((s2 - s1 .^ 2 ./ n) ./ (n - 1), 0) ./ n;
                    end
                    acc = acc + eyeSq(X{c}(:, o), ref, use(o), c, C.loc{c}, K);   % what blinks add (total error)
                end
            end
            obj = sqrt(acc / (nC * size(X{1}, 2)));
            obj(~ok | ~isfinite(obj)) = NaN;
        end

        function tf = applies(contract)
            % per-trial scores exist: band power, or only mean amplitudes
            tf = strcmp(contract.analysis, 'bandpower') || ...
                (~isempty(contract.components) && all(strcmp({contract.components.measure}, 'mean')));
        end

        function Bq = bootstrap(X, M, v, opts, ref, C)
            % The objective at each limit v under the paired bootstrap of the
            % ranking (same resampled trials for every limit): limits x B.
            % C (curve's): with it, what blinks add is resampled too.
            s = RandStream('mt19937ar', 'Seed', opts.seed + 202);
            segmented = isfield(ref, 'segmented') && ref.segmented;
            nC = numel(X); B = opts.nBoot;
            acc = zeros(numel(v), B);
            for c = 1:nC
                if segmented
                    loc = 1:size(X{c},1); if nargin > 5, loc = C.loc{c}; end
                    fullM = nan(ref.n(c),1); fullM(loc) = M{c};
                    [~, I] = pipecompare.eval.Rank.resampleWeights(s,ref.n(c),B,true);
                    for o = 1:size(X{c},2)
                        fullX = nan(ref.n(c),1); fullX(loc) = X{c}(:,o);
                        Z = reshape(fullX(I),size(I)); marks = reshape(fullM(I),size(I));
                        for j = 1:numel(v)
                            [~, se] = pipecompare.eval.Measure.meanError(Z,double(marks<=v(j)),true);
                            acc(j,:) = acc(j,:) + se.^2;
                        end
                    end
                    continue;
                end
                W = pipecompare.eval.Rank.resampleWeights(s, size(X{c}, 1), B, segmented);
                K = double(M{c}(:) <= v(:)');                 % trials x limits
                Nb = K' * W;
                for o = 1:size(X{c}, 2)
                    x = X{c}(:, o);
                    s1 = K' * (x .* W); s2 = K' * (x .^ 2 .* W);
                    sm2 = max((s2 - s1 .^ 2 ./ Nb) ./ (Nb - 1), 0) ./ Nb;
                    sm2(Nb < 2) = Inf;
                    if nargin > 5
                        for j = 1:numel(v)
                            sm2(j, :) = sm2(j, :) + eyeSq(x, ref, C.use(o), c, C.loc{c}, K(:, j) .* W, K(:, j) .* W);
                        end
                    end
                    acc = acc + sm2;
                end
            end
            Bq = sqrt(acc / (nC * size(X{1}, 2)));
        end
    end
end

function q = eyeSq(x, ref, k, c, loc, K, counts)
% The squared bias blinks add to measure k's mean in condition c
% (Measure.eyeBias) for the rows loc of the reference trials (scores x,
% weights K, a column per mean); 0 when the blinks are not known.
q = zeros(1, size(K, 2));
blink = pipecompare.utils.fieldOr(ref, 'blink', {});
if numel(blink) < k || isempty(blink{k}), return; end
mode=pipecompare.utils.fieldOr(ref,'artifactVariance','pooled');
if nargin<7||strcmp(mode,'legacy-total'),counts=ones(size(K));end
[~, q] = pipecompare.eval.Measure.eyeBias(x, blink{k}{c}(loc), K,counts,mode);
end

function b = nextAbove(v, j)
% the next value of m above the limit (the same epochs are kept up to it)
if j < numel(v), b = v(j + 1); else, b = 2 * v(j); end   % above all: any limit keeps every epoch
end

function u = roundIn(a, b)
% the roundest number in [a, b): fewest significant digits
for d = -4:12
    u = ceil(a * 10 ^ d) / 10 ^ d;
    if u >= a && u < b, return; end
end
u = a;
end
