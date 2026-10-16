classdef Cutoff
    %CUTOFF Filter cutoffs fine-tuned from the data, with the decisions of a
    %   pipeline held fixed.
    %
    %   [cut, info] = pipecompare.eval.Cutoff.search(D, rest, dec, base, env)
    %
    %   D: the data just before the pipeline's first filter (continuous).
    %   rest: the pipeline's steps from that filter on: one or two filters
    %   (highpass, lowpass), then only IC removal, epoching, baseline
    %   removal, an amplitude rejection that removes epochs and weighting
    %   the epochs.
    %   dec: the decisions the pipeline took on these data: comps (the
    %   components it removed), rejIdx (the epochs it rejected) and
    %   optionally weights (the epochs' weights, trialWeights).
    %   base: the pipeline's measures (m) and signal gains (gain).
    %   env: contract, ref, rank (the ranking's options), S and truth (the
    %   signal-check copy of the starting data), and grid: per filter of
    %   rest, the cutoffs to try (empty: keep the pipeline's).
    %
    %   With the decisions fixed, every step from the first filter to the
    %   scores is linear: the filters act on each channel alone and the
    %   same way on every channel, IC removal is the same matrix at every
    %   sample, epoching selects samples and baseline removal subtracts
    %   their mean. So the scores of a cutoff need only the measured
    %   electrodes: the components are removed first (the operations
    %   commute), the other channels are dropped, and only the measured
    %   ones are filtered, epoched and scored, exactly as the whole data
    %   would be there. The signal gain of a cutoff is the pipeline's gain
    %   times what its filters do to the known signal, relative to the
    %   pipeline's filters (the known signal is a time course times a
    %   scalp field; the filters change only the time course), measured
    %   with the filters-alone check (Injection.filterCheck), which also
    %   excludes cutoffs whose filters alone distort the signal. The
    %   objective of the ranking (gain-corrected SME) is computed for
    %   every pair of cutoffs on the grid. At the pipeline's own cutoffs
    %   it must equal the pipeline's objective; when it does not, the
    %   steps do not commute as assumed and nothing is fine-tuned.
    %
    %   cut: the cutoffs (one per filter of rest) with the lowest
    %   objective, or [] when the pipeline's own are the best or nothing
    %   could be fine-tuned (info.note says why). info.grid and
    %   info.objectives: every pair of cutoffs tried and its objective.

    properties (Constant)
        Tolerance = 1e-3   % relative: the objective at the pipeline's own cutoffs must reproduce it
    end

    methods (Static)
        function [cut, info] = search(D, rest, dec, base, env)
            cut = [];
            contract = env.contract; ref = env.ref; opts = env.rank;
            types = cellfun(@(i) i.type, rest, 'UniformOutput', false);
            fi = find(ismember(types, {'highpass', 'lowpass'}));
            own = cellfun(@(i) i.params.cutoff, rest(fi));
            info = struct('cutoffs', own, 'objectiveBase', NaN, 'objectiveBest', NaN, 'pairs', 0, 'feasible', 0, ...
                'reproduced', NaN, 'grid', [], 'objectives', [], 'note', '');
            gBase = double(base.gain(:)');
            spBase = pipecompare.utils.fieldOr(base, 'spread', []);
            info.objectiveBase = pipecompare.eval.Rank.pointObjective(base.m, ref, opts, gBase, spBase);
            % the measured electrodes, the components already removed
            E = D;
            if ~isempty(dec.comps), [~, E] = evalc('pop_subcomp(E, dec.comps, 0)'); end
            E.icaweights = []; E.icasphere = []; E.icawinv = []; E.icaact = []; E.icachansind = [];
            keep = find(ismember(lower({E.chanlocs.labels}), lower(contract.allRoi())));
            [~, E] = evalc('pop_select(E, ''channel'', keep)');
            E.data = double(E.data);
            % the epochs as samples of the continuous data: the epoching
            % of a dataset holding each sample's index (pop_epoch only
            % selects samples, so any data give the same epochs)
            ep = rest{strcmp(types, 'epoch')};
            ctx = struct('contract', contract);
            [~, T0] = evalc('pipecompare.run.Steps.run(ep, E, ctx)');
            I = E; I.data = repmat(1:E.pnts, E.nbchan, 1);
            [~, I] = evalc('pipecompare.run.Steps.run(ep, I, ctx)');
            idx = reshape(round(double(I.data(1, :, :))), I.pnts, I.trials);
            ok = I.trials == T0.trials && I.pnts == T0.pnts && ...
                max(abs(double(T0.data(:)) - reshape(E.data(:, idx(:)), [], 1))) <= 1e-6 * max(abs(E.data(:)));
            assert(ok, 'PipeCompare:Cutoff', 'The epochs could not be located in the continuous data.');
            bl = rest(strcmp(types, 'baseline'));
            % the grid: every pair of cutoffs, the pipeline's own included
            vals = env.grid;
            for k = 1:numel(fi), if isempty(vals{k}), vals{k} = own(k); end, end
            for k = 1:numel(fi), vals{k} = unique([vals{k}(:)' own(k)]); end
            if numel(fi) == 1, pairs = vals{1}(:); else, [b, a] = ndgrid(vals{2}, vals{1}); pairs = [a(:) b(:)]; end   % first filter slowest
            info.pairs = size(pairs, 1);
            % what the filters alone do to the known signal (gain, check)
            insts = @(r) arrayfun(@(k) pipecompare.plan.Plan.reparam(rest{fi(k)}, 'cutoff', r(k)), 1:numel(fi), 'UniformOutput', false);
            P = cell(1, size(pairs, 1));
            for r = 1:size(pairs, 1), P{r} = insts(pairs(r, :)); end
            leaves = struct('path', P);
            [why, sig] = pipecompare.eval.Injection.filterCheck(env.S, contract, env.truth, leaves, opts);
            sp = pipecompare.eval.Injection.shapeSpreads(env.S, contract, env.truth, leaves);   % the total error's shape part
            ib = find(all(abs(pairs - own) <= 1e-12 * max(1, abs(own)), 2), 1);
            assert(~isempty(sig{ib}), 'PipeCompare:Cutoff', 'The filters alone could not be checked on the known signal.');
            gF0 = double(sig{ib}.gain(:)');
            obj = nan(size(pairs, 1), 1);
            outer = NaN; F1 = [];
            % a low-pass after the first filter is computed only at the
            % epochs' samples, as firfilt computes it there (its kernel from
            % pop_eegfiltnew): the same numbers, without filtering the
            % whole recording for every pair
            atEpochs = numel(fi) == 2 && strcmp(rest{fi(2)}.type, 'lowpass') && exist('findboundaries', 'file') == 2;
            if atEpochs, dc = [findboundaries(E.event) E.pnts + 1]; kernels = containers.Map('KeyType', 'double', 'ValueType', 'any'); end
            for r = 1:size(pairs, 1)
                if ~isempty(why{r}) && r ~= ib, continue; end
                g = gBase .* double(sig{r}.gain(:)') ./ gF0;
                if ~all(g > 0), continue; end
                % the first filter once per value, the second on its output
                if pairs(r, 1) ~= outer
                    in1 = insts(pairs(r, :)); in1 = in1{1};
                    [~, F1] = evalc('pipecompare.run.Steps.run(in1, E, struct())'); outer = pairs(r, 1);
                    Z = double(F1.data);
                end
                if atEpochs
                    v2 = pairs(r, 2);
                    if ~isKey(kernels, v2), kernels(v2) = lowpassKernel(v2, E.srate); end
                    Y = filterAt(Z, kernels(v2), idx, dc);
                elseif numel(fi) == 2
                    in2 = insts(pairs(r, :)); in2 = in2{2};
                    [~, F] = evalc('pipecompare.run.Steps.run(in2, F1, struct())');
                    Y = double(F.data(:, idx(:)));
                else
                    Y = Z(:, idx(:));
                end
                T = T0;
                T.data = cast(reshape(Y, E.nbchan, I.pnts, I.trials), 'like', T0.data);
                if ~isempty(bl), [~, T] = evalc('pipecompare.run.Steps.run(bl{1}, T, ctx)'); end
                if ~isempty(dec.rejIdx), [~, T] = evalc('pop_rejepoch(T, dec.rejIdx, 0)'); end
                if isfield(dec, 'weights') && ~isempty(dec.weights), T.etc.pipecompare.trialWeights = dec.weights; end
                m = pipecompare.eval.Measure.candidate(T, contract, ref, opts);
                spr = sp{r}; if isempty(spBase), spr = []; end   % as the pipeline's own objective
                obj(r) = pipecompare.eval.Rank.pointObjective(m, ref, opts, g, spr);
            end
            info.feasible = sum(isfinite(obj));
            info.reproduced = obj(ib);
            info.grid = pairs; info.objectives = obj;   % per pair of cutoffs (NaN: not a candidate)
            if ~(abs(obj(ib) - info.objectiveBase) <= pipecompare.eval.Cutoff.Tolerance * info.objectiveBase)
                info.note = sprintf(['the objective at the pipeline''s own cutoffs (%.4g) does not reproduce the ', ...
                    'pipeline''s (%.4g), so the steps do not commute as assumed'], obj(ib), info.objectiveBase);
                return;
            end
            obj(ib) = info.objectiveBase;
            [info.objectiveBest, best] = min(obj);
            if best == ib || ~(info.objectiveBest < info.objectiveBase)
                info.objectiveBest = info.objectiveBase;
                info.note = 'the pipeline''s own cutoffs are the best on the grid';
                return;
            end
            cut = pairs(best, :);
            info.cutoffs = cut;
        end

        function v = grid(vals)
            % Cutoffs between the lowest and the highest compared, in steps
            % of a round number (1, 2 or 5 times a power of ten) giving at
            % most about 20 values, the compared ones included.
            vals = unique(vals(:)');
            lo = vals(1); hi = vals(end);
            if hi <= lo, v = vals; return; end
            raw = (hi - lo) / 19;
            p = 10 ^ floor(log10(raw));
            s = p * [1 2 5 10]; s = s(find(s >= raw - 1e-12 * raw, 1));
            v = unique(round([lo:s:hi vals] / (s / 100)) * (s / 100));
            v = v(v >= lo - 1e-12 & v <= hi + 1e-12);
        end
    end
end

function b = lowpassKernel(cutoff, fs)
% the FIR pop_eegfiltnew designs for a low-pass at cutoff (it depends only
% on the cutoff and the sampling rate), as the low-pass step calls it
D = eeg_emptyset(); D.srate = fs; D.nbchan = 1; D.trials = 1;
D.pnts = round(20 * fs); D.data = zeros(1, D.pnts); D.xmin = 0; D.xmax = (D.pnts - 1) / fs;
D = eeg_checkset(D);
[~, ~, ~, b] = evalc('pop_eegfiltnew(D, ''hicutoff'', cutoff, ''plotfreqz'', 0)');
b = double(b(:)');
assert(mod(numel(b), 2) == 1, 'PipeCompare:Cutoff', 'The low-pass FIR has an even length.');
end

function Y = filterAt(Z, b, idx, dc)
% firfilt's output at the samples idx (epoch samples x epochs) of Z
% (channels x samples): y(n) = sum_j b(j) z(n + g - j + 1), g the group
% delay, with each stretch between boundaries (dc: their first samples,
% then one past the end, as firfilt finds them) padded by its first and
% last value
g = (numel(b) - 1) / 2; P = size(idx, 1);
Y = zeros(size(Z, 1), numel(idx));
for e = 1:size(idx, 2)
    n = idx(:, e)';
    s = dc(find(dc <= n(1), 1, 'last')); t = dc(find(dc > n(end), 1)) - 1;
    src = min(max((n(1) - g):(n(end) + g), s), t);
    Y(:, (e - 1) * P + (1:P)) = conv2(Z(:, src), b, 'valid');
end
end
