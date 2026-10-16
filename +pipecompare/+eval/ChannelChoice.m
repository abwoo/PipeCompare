classdef ChannelChoice
    %CHANNELCHOICE Which bad channels to interpolate, chosen by their effect
    %   on the total error (badchannels with repair = 'error').
    %
    %   [chosen, info] = pipecompare.eval.ChannelChoice.choose(EEG, src, bad, ctx)
    %
    %   bad: the channels the bad-channel test flagged (indices); src: the
    %   copy they were tested on (high- and low-passed to the band the
    %   analysis keeps). Interpolation is linear: the data with the set b
    %   interpolated are M_b X, M_b EEGLAB's spherical spline (pop_interp),
    %   read off by interpolating an identity matrix. The steps that follow
    %   in every pipeline below the step (ctx.after) are applied as they
    %   would run where they are linear and the same for all of them: the
    %   re-referencing (a channel map, applied after M_b), epoching and
    %   baseline removal (per channel, so they commute with M_b and are
    %   applied once), and the epoch rejection that follows, at its best
    %   limit when it is chosen from the data. The filters in between are
    %   replaced by the band of src, and ICA is left out (it differs
    %   between the pipelines, and what it does is not linear). The
    %   objective is the ranking's total error (IcThreshold.score), the
    %   gain from the known-signal copy (ctx.signal) with the same channels
    %   interpolated, so interpolating a measured electrode counts with the
    %   signal it loses there.
    %
    %   From every flagged channel interpolated, the channel whose being
    %   left as it is gives the lowest objective is left, one at a time:
    %   nested sets from all to none. Among the sets shown equivalent to the
    %   best one (paired bootstrap, the whole interval of the difference
    %   within equivalenceMargin of the best objective, as IcThreshold),
    %   the one interpolating the most is chosen: a flagged channel is left
    %   as it is only when interpolating it would make the measure
    %   measurably worse. More than MaxChannels flagged channels are all
    %   interpolated. info: flagged, chosen, objectives (per set, all
    %   first), sets, equallyGood, note.

    properties (Constant)
        MaxChannels = 12   % more flagged channels than this: all interpolated
    end

    methods (Static)
        function [chosen, info] = choose(EEG, src, bad, ctx)
            bad = bad(:)'; chosen = bad;
            info = struct('flagged', bad, 'chosen', bad, 'objectives', [], 'sets', {{}}, 'equallyGood', 0, 'note', '');
            if numel(bad) > pipecompare.eval.ChannelChoice.MaxChannels
                info.note = sprintf('%d channels flagged (more than %d): all interpolated', numel(bad), ...
                    pipecompare.eval.ChannelChoice.MaxChannels);
                pipecompare.utils.log('Bad channels: %s.', info.note);
                return;
            end
            contract = ctx.contract;
            ref = pipecompare.utils.fieldOr(ctx, 'ref', []);
            if isempty(ref), ref = pipecompare.eval.Measure.reference(EEG, contract); end
            opts = pipecompare.utils.withDefaults(pipecompare.utils.fieldOr(ctx, 'rank', struct()), ...
                pipecompare.eval.Rank.defaults());
            [timeSteps, chanSteps, rej] = split(pipecompare.utils.fieldOr(ctx, 'after', {}));
            Ep = runSteps(src, timeSteps, ctx);
            sig = pipecompare.utils.fieldOr(ctx, 'signal', []);
            Sp = [];
            if ~isempty(sig) && ~isempty(sig.S), Sp = runSteps(sig.S, timeSteps, ctx); end
            chans = [];
            if ~isempty(rej), chans = pipecompare.run.Steps.testedChannels(Ep, rej.params); end
            score = @(b) scoreOf(EEG, Ep, Sp, b, chanSteps, rej, chans, ctx, ref, opts, sig);
            % the path: every flagged channel interpolated, then one fewer at a time
            cur = bad; sets = {cur};
            [f, u, X1, M1, C1] = score(cur);
            X = {X1}; M = {M1}; C = {C1};
            for step = 1:numel(bad)
                fb = Inf; kb = 0; out = {};
                for c = cur
                    [fc, uc, Xc, Mc, Cc] = score(setdiff(cur, c, 'stable'));
                    if fc < fb, fb = fc; kb = c; out = {uc, Xc, Mc, Cc}; end
                end
                if kb == 0, break; end                     % no further set passes
                cur = setdiff(cur, kb, 'stable');
                sets{end+1} = cur; f(end+1) = fb; u(end+1) = out{1}; %#ok<AGROW>
                X{end+1} = out{2}; M{end+1} = out{3}; C{end+1} = out{4}; %#ok<AGROW>
            end
            info.objectives = f; info.sets = sets;
            F = find(isfinite(f));
            if isempty(F)
                info.note = 'no set of interpolated channels meets the limits: all interpolated';
                pipecompare.utils.log('Bad channels: %s.', info.note);
                return;
            end
            best = F(find(f(F) == min(f(F)), 1));     % ties: more channels interpolated
            Bq = zeros(numel(F), opts.nBoot);
            for i = 1:numel(F)
                Bq(i, :) = pipecompare.eval.Threshold.bootstrap(X{F(i)}, M{F(i)}, u(F(i)), opts, ref, C{F(i)});
            end
            % the most channels interpolated shown equivalent to the best, going towards all from it
            ib = find(F == best);
            [c, nEq] = pipecompare.eval.Rank.fixedSequence(f(F)', Bq, ib, ib-1:-1:1, opts.alpha, opts.equivalenceMargin);
            k = F(c);
            chosen = sort(sets{k}); info.chosen = chosen; info.equallyGood = nEq;
            left = setdiff(bad, chosen);
            if isempty(left)
                pipecompare.utils.log(['Bad channels by their effect on the total error: all %d interpolated (leaving one as ', ...
                    'it is would not lower the error measurably).'], numel(bad));
            else
                pipecompare.utils.log(['Bad channels by their effect on the total error: %d of %d interpolated; left as they ', ...
                    'are: %s (interpolating them would raise the error from %.4g to %.4g).'], numel(chosen), numel(bad), ...
                    strjoin({EEG.chanlocs(left).labels}, ', '), f(k), f(1));
            end
        end

        function M = interpMatrix(EEG, b)
            % The channels x channels matrix of pop_interp(EEG, b,
            % 'spherical'): the data with channels b interpolated are M * X.
            n = EEG.nbchan; M = eye(n);
            if isempty(b), return; end
            E = eeg_emptyset();
            E.chanlocs = EEG.chanlocs; E.chaninfo = EEG.chaninfo; E.nbchan = n; E.srate = EEG.srate;
            E.data = eye(n); E.pnts = n; E.trials = 1; E.xmin = 0; E.xmax = (n - 1) / EEG.srate;
            [~, E] = evalc('pop_interp(E, b, ''spherical'')');
            M = double(E.data);
        end
    end
end

function [f, u, X, M, C] = scoreOf(EEG, Ep, Sp, b, chanSteps, rej, chans, ctx, ref, opts, sig)
% the objective with the channels b interpolated (NaN when that set fails
% the signal check or keeps too few trials)
T = pipecompare.eval.ChannelChoice.interpMatrix(EEG, b);
Ek = runSteps(mapped(Ep, T), chanSteps, ctx);
Sk = [];
if ~isempty(Sp)
    Sk = runSteps(mapped(Sp, T), chanSteps, ctx);
    % the field the check expects is re-referenced as the copy was
    for i = 1:numel(chanSteps), Sk = pipecompare.eval.Injection.noteReference(Sk, chanSteps{i}); end
end
[f, u, ~, X, M, C] = pipecompare.eval.IcThreshold.score(Ek, Sk, rej, chans, ctx.contract, ref, opts, sig);
end

function E = mapped(E, T)
% E with its channels mixed by T (channels x channels)
sz = size(E.data);
E.data = reshape(T * double(reshape(E.data, sz(1), [])), sz);
end

function [timeSteps, chanSteps, rej] = split(after)
% Of the steps that follow, up to the epoch rejection: epoching and
% baseline removal (per channel), the re-referencing (a channel map) and
% the rejection; the others (filters, ICA, ...) are left out.
rej = [];
k = find(cellfun(@(in) strcmp(in.type, 'reject_threshold'), after), 1);
if ~isempty(k), rej = after{k}; after = after(1:k-1); end
types = cellfun(@(in) in.type, after, 'UniformOutput', false);
timeSteps = after(ismember(types, {'epoch', 'baseline'}));
chanSteps = after(strcmp(types, 'reref'));
end

function E = runSteps(E, steps, ctx)
% the steps, then epochs (as the measure takes them) when the data are
% still continuous
for i = 1:numel(steps)
    [~, E] = evalc('pipecompare.run.Steps.run(steps{i}, E, ctx)');
end
if E.trials == 1
    [~, E] = evalc('pop_epoch(E, ctx.contract.allEvents(), ctx.contract.epoch, ''epochinfo'', ''yes'')');
end
end
