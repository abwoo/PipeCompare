classdef Steps
    %STEPS Execute one plan step with native EEGLAB functions.
    %
    %   [EEG, coms, info] = pipecompare.run.Steps.run(inst, EEG, ctx)
    %
    %   coms are the EEGLAB commands that reproduce the step; the executor
    %   appends them to the candidate's EEG.history. info carries the
    %   data-dependent decisions the step took (bad channels, ICA matrices,
    %   removed components, rejected epochs) and counts for the evaluation.
    %
    %   EEG = pipecompare.run.Steps.replayDecision(inst, EEG, info, ctx)
    %
    %   Applies the SAME decisions to another dataset with the same
    %   structure (used by pipecompare.eval.Injection to measure how a known
    %   signal is transferred by exactly the operations applied to the
    %   real data). ASR replays its recorded window-by-window
    %   reconstructions (pipecompare.run.AsrRecord, verified against EEGLAB's
    %   output); captured mark/remove workflows replay the epochs and
    %   components they removed; epoch rejection replays the channels it
    %   interpolated within epochs (same epochs, same channels) and the
    %   epochs it removed. Other native commands are re-run and
    %   reported as not decision-matched. Weighting epochs replays the
    %   weights of the real data.

    properties (Constant)
        MinWindows = 10   % a stretch of a channel bad for part of the recording lasts 10 s or more (badStretches)
        StretchZ = 3      % and its noise is on average 3 robust SDs above the channel's usual level
    end

    methods (Static)
        function [EEG, coms, info] = run(inst, EEG, ctx)
            p = inst.params;
            coms = {};
            info = struct();
            switch inst.type
                case 'resample'
                    [EEG, com] = pop_resample(EEG, p.fs);
                    coms = {com};
                    if EEG.srate ~= p.fs && abs(EEG.srate - p.fs) <= 1e-9 * p.fs
                        % pop_resample leaves e.g. 250.00000000000003 when the
                        % original rate was not exactly representable; only a
                        % floating-point residue is corrected (ICLabel needs it).
                        EEG.srate = p.fs; EEG = eeg_checkset(EEG);
                        coms{end+1} = sprintf('EEG.srate = %g; EEG = eeg_checkset(EEG); %% PipeCompare: exact rate after resampling', p.fs);
                    end
                case 'highpass'
                    [fa, fb] = fftArgs(EEG, p.cutoff);   % fb: the stand-in fftfilt while it runs
                    [EEG, com] = pop_eegfiltnew(EEG, 'locutoff', p.cutoff, 'plotfreqz', 0, fa{:});
                    coms = {plainCommand(com, fb)};
                case 'lowpass'
                    [EEG, com] = pop_eegfiltnew(EEG, 'hicutoff', p.cutoff, 'plotfreqz', 0);
                    coms = {com};
                case 'linenoise'
                    [EEG, com] = pop_eegfiltnew(EEG, 'locutoff', p.freq - p.halfwidth, ...
                        'hicutoff', p.freq + p.halfwidth, 'revfilt', 1, 'plotfreqz', 0);
                    coms = {com};
                case 'asr'
                    assert(exist('pop_clean_rawdata', 'file') == 2, 'PipeCompare:Dependency', 'clean_rawdata plugin not installed');
                    E0 = EEG;
                    [EEG, com] = pop_clean_rawdata(EEG, 'FlatlineCriterion', 'off', 'ChannelCriterion', 'off', ...
                        'LineNoiseCriterion', 'off', 'Highpass', 'off', 'BurstCriterion', p.cutoff, ...
                        'WindowCriterion', 'off', 'BurstRejection', 'off', 'Distance', 'Euclidian', ...
                        'MaxMem', pipecompare.run.AsrRecord.maxMemMB(EEG.nbchan, EEG.srate));
                    coms = {com};
                    info.asr = recordAsr(E0, EEG, p.cutoff);
                case 'badchannels'
                    [EEG, coms, info] = pipecompare.run.Steps.badChannels(EEG, p, ctx);
                case 'channels'
                    [EEG, coms, info] = pipecompare.run.Steps.listedChannels(EEG, p);
                case 'restore'
                    % back to the starting montage, then the channels removed
                    % before PipeCompare (EEG.etc.pipecompare.preRemoved)
                    target = rootChanlocs(EEG);
                    missing = setdiff(lower({target.labels}), lower({EEG.chanlocs.labels}));
                    pre = [];
                    if isfield(EEG.etc.pipecompare, 'preRemoved'), pre = EEG.etc.pipecompare.preRemoved; end
                    if ~isempty(missing) || isempty(pre)
                        requireLocations(target(ismember(lower({target.labels}), missing)), 'restore');
                        EEG = pop_interp(EEG, target, 'spherical');
                        coms = {'EEG = pop_interp(EEG, EEG.etc.pipecompare.rootChanlocs, ''spherical'');'};
                    end
                    if ~isempty(pre)
                        absent = ~ismember(lower({pre.labels}), lower({EEG.chanlocs.labels}));
                        if any(absent)
                            requireLocations(pre(absent), 'restore');
                            EEG = pop_interp(EEG, pre(absent), 'spherical');
                            coms{end+1} = sprintf(['EEG = pop_interp(EEG, EEG.etc.pipecompare.preRemoved(%s), ''spherical''); ', ...
                                '%% PipeCompare: channels removed before the plan'], mat2str(find(absent)));
                            missing = [missing lower({pre(absent).labels})];
                        end
                    end
                    info.interpolated = missing;
                case 'reref'
                    ex = {};
                    if isfield(p, 'exclude') && ~isempty(p.exclude)
                        present = presentChannels(EEG, p.exclude, 'reref exclude');
                        xi = find(ismember(lower({EEG.chanlocs.labels}), lower(present)));
                        if ~isempty(xi), ex = {'exclude', xi}; end
                    end
                    if strcmp(p.mode, 'average')
                        [EEG, com] = pop_reref(EEG, [], ex{:});
                    else
                        ch = cellstr(p.channels);
                        assert(~isempty(ch), 'PipeCompare:Reref', 'reref mode channels needs channel labels');
                        [EEG, com] = pop_reref(EEG, ch, ex{:});
                    end
                    coms = {com};
                    info.reref = p;
                case 'ica'
                    [EEG, coms] = pipecompare.run.Steps.ica(EEG, p, ctx);
                    info.ica = struct('icaweights', EEG.icaweights, 'icasphere', EEG.icasphere, ...
                        'icachansind', EEG.icachansind, 'coms', {coms});
                case 'icremove'
                    [EEG, coms, info] = pipecompare.run.Steps.icRemove(EEG, p, ctx);
                case 'epoch'
                    c = ctx.contract;
                    [EEG, ~, com] = pop_epoch(EEG, c.allEvents(), c.epoch, 'epochinfo', 'yes');
                    coms = {com};
                case 'baseline'
                    c = ctx.contract;
                    assert(~isempty(c.baseline), 'PipeCompare:Baseline', 'The contract defines no baseline window.');
                    % clamp to the epoch limits (they can differ from the
                    % contract by less than one sample after resampling)
                    b = 1000 * c.baseline;
                    assert(b(1) >= 1000 * EEG.xmin - 1000 / EEG.srate && b(2) <= 1000 * EEG.xmax + 1000 / EEG.srate, ...
                        'PipeCompare:Baseline', 'Baseline [%g %g] ms lies outside the epochs [%g %g] ms', b(1), b(2), 1000 * EEG.xmin, 1000 * EEG.xmax);
                    b = [max(b(1), 1000 * EEG.xmin) min(b(2), 1000 * EEG.xmax)];
                    [EEG, com] = pop_rmbase(EEG, b, []);
                    coms = {com};
                case {'reject_threshold','reject_jointprob','reject_kurtosis'}
                    [EEG, coms, info] = pipecompare.run.Steps.rejectEpochs(EEG, inst.type, p, ctx);
                case 'eogregress'
                    c = ctx.contract; groups = {c.conditions.events};
                    [EEG, info.eog] = pipecompare.run.Steps.eogRegress(EEG, p.channels, groups, c.epoch);
                    coms = {sprintf(['EEG = pipecompare.run.Steps.eogRegress(EEG, %s, %s, %s); %% PipeCompare: the eye ', ...
                        'channels %s regressed out of the EEG channels'], strCells(info.eog.eyes), ...
                        ['{' strjoin(cellfun(@strCells, groups, 'UniformOutput', false), ', ') '}'], ...
                        mat2str(c.epoch, 17), strjoin(info.eog.eyes, ', '))};
                case 'weight_trials'
                    c = ctx.contract;
                    assert(strcmp(c.analysis, 'erp'), 'PipeCompare:Weights', ...
                        'Weighting epochs by their noise is for ERPs (a weighted mean of band power would favour the epochs with less power).');
                    roi = c.allRoi(); groups = {c.conditions.events};
                    scores = arrayfun(@(m) {cellstr(m.roi), m.window}, c.components, 'UniformOutput', false);
                    [EEG, info] = pipecompare.run.Steps.weightTrials(EEG, roi, groups, scores);
                    coms = {sprintf(['EEG = pipecompare.run.Steps.weightTrials(EEG, %s, %s, %s); %% PipeCompare: epochs ', ...
                        'weighted by their noise (EEG.etc.pipecompare.trialWeights; weighted average: ', ...
                        'pipecompare.eval.Measure.epochWeights)'], strCells(roi), ...
                        ['{' strjoin(cellfun(@strCells, groups, 'UniformOutput', false), ', ') '}'], ...
                        ['{' strjoin(cellfun(@(q) sprintf('{%s, [%.17g %.17g]}', strCells(q{1}), q{2}), scores, ...
                        'UniformOutput', false), ', ') '}'])};
                case 'native'
                    [EEG, coms, info] = pipecompare.run.Steps.native(EEG, p.command);
                otherwise
                    error('PipeCompare:UnknownStep', 'Unknown step %s', inst.type);
            end
            coms = coms(~cellfun(@isempty, coms));
        end

        function [EEG, matched] = replayDecision(inst, EEG, info, ctx)
            % Same decisions on a structurally identical dataset.
            matched = true;
            EEG = pipecompare.eval.Injection.noteReference(EEG, inst);
            switch inst.type
                case 'badchannels'
                    if strcmp(inst.params.action, 'remove')
                        if ~isempty(info.badIdx), EEG = pop_select(EEG, 'nochannel', info.badIdx); end
                    elseif ~isempty(info.badIdx)
                        EEG = pop_interp(EEG, info.badIdx, 'spherical');
                    end
                    if isfield(info, 'stretches') && ~isempty(info.stretches)
                        EEG = pipecompare.run.Steps.interpStretches(EEG, info.stretches);
                    end
                case 'channels'
                    if ~isempty(info.listedIdx)
                        if strcmp(inst.params.action, 'remove'), EEG = pop_select(EEG, 'nochannel', info.listedIdx);
                        else, EEG = pop_interp(EEG, info.listedIdx, 'spherical'); end
                    end
                case 'ica'
                    EEG.icaweights = info.ica.icaweights; EEG.icasphere = info.ica.icasphere;
                    EEG.icachansind = info.ica.icachansind; EEG.icawinv = []; EEG.icaact = [];
                    EEG = eeg_checkset(EEG);
                case 'icremove'
                    if ~isempty(info.comps), EEG = pop_subcomp(EEG, info.comps, 0); end
                case {'reject_threshold','reject_jointprob','reject_kurtosis'}
                    if isfield(info, 'epochInterp') && ~isempty(info.epochInterp.chans)
                        EEG = pipecompare.run.Steps.interpolateEpochs(EEG, info.epochInterp.chans, info.epochInterp.epochs);
                    end
                    if ~isempty(info.rejIdx), EEG = pop_rejepoch(EEG, info.rejIdx, 0); end
                case 'eogregress'
                    % the coefficients estimated on the real data
                    [~, ti] = ismember(info.eog.targets, {EEG.chanlocs.labels});
                    [~, ei] = ismember(info.eog.eyes, {EEG.chanlocs.labels});
                    EEG = applyEog(EEG, ti, ei, info.eog.B);
                case 'weight_trials'
                    EEG.etc.pipecompare.trialWeights = info.trialWeights;
                case 'native'
                    if isfield(info, 'decisions') && ~isempty(info.decisions)
                        % a workflow of marks and removals: the removals the
                        % real data decided, applied to this copy
                        for d = info.decisions
                            if strcmp(d.kind, 'epochs') && ~isempty(d.idx), EEG = pop_rejepoch(EEG, d.idx, 0); end
                            if strcmp(d.kind, 'comps') && ~isempty(d.idx), EEG = pop_subcomp(EEG, d.idx, 0); end
                        end
                    else
                        % re-run; that is the same operation when every
                        % statement is a fixed transform (no data-driven choice)
                        EEG = pipecompare.run.Steps.run(inst, EEG, ctx);
                        matched = all(cellfun(@isFixedTransform, pipecompare.run.Native.statements(inst.params.command)));
                    end
                case {'highpass', 'lowpass', 'linenoise'}
                    EEG = lowRankRun(inst, EEG, ctx);
                case 'asr'
                    if isfield(info, 'asr') && ~isempty(info.asr)
                        % the windows and reconstructions ASR chose on the real data
                        EEG.data = cast(pipecompare.run.AsrRecord.apply(info.asr, EEG.data), 'like', EEG.data);
                    else
                        matched = false;
                        EEG = pipecompare.run.Steps.run(inst, EEG, ctx);
                    end
                otherwise
                    EEG = pipecompare.run.Steps.run(inst, EEG, ctx);
            end
        end

        function [EEG, coms, info] = badChannels(EEG, p, ctx)
            elec = 1:EEG.nbchan;
            if isfield(p, 'exclude') && ~isempty(p.exclude)   % e.g. EOG / reference channels are not tested
                elec = find(~ismember(lower({EEG.chanlocs.labels}), lower(presentChannels(EEG, p.exclude, 'badchannels exclude'))));
            end
            % measure: one measure, or several joined by '+' (e.g.
            % 'kurt+prob+resid'): a channel is bad when any of them flags it.
            % kurt, prob and spec are pop_rejchan's; resid is
            % PipeCompare's (residualZ: noise the other channels do not
            % explain), which finds what kurtosis and joint probability
            % cannot see: both depend only on the shape of a channel's
            % values, not their size, so a channel that is just noisier
            % (poor contact) passes them.
            measures = strsplit(p.measure, '+');
            % Detect on a high-passed copy (slow drifts distort kurtosis
            % and probability), and low-passed (detectLowpass, the band
            % the analysis keeps: line noise alone that a later low-pass
            % removes does not make a channel bad); interpolate or remove
            % on the data as they are.
            hp = pipecompare.utils.fieldOr(p, 'detectHighpass', 0);
            lp = pipecompare.utils.fieldOr(p, 'detectLowpass', 0);
            applied = 0; if nargin > 2 && isfield(ctx, 'highpass'), applied = ctx.highpass; end
            doHp = hp > 0 && hp > applied;
            doLp = lp > 0 && lp < EEG.srate / 2;
            copy = doHp || doLp;
            lpOnly = bandCopy(EEG, false, hp, doLp, lp);     % (the filters commute: low-pass first)
            src = bandCopy(lpOnly, doHp, hp, false, lp); on = '';
            if doHp, on = sprintf(' on a %g Hz high-passed', hp); end
            if doLp
                on = [pipecompare.utils.ternary(doHp, [on ','], ' on a') sprintf(' %g Hz low-passed', lp)];
            end
            if copy, on = [on ' copy']; end
            % Flat channels (no signal) are bad whatever the measures say:
            % kurtosis and probability are undefined on a constant, and left
            % in, they would also spoil the normalisation over the others.
            flat = elec(pipecompare.run.Steps.flatChannels(src.data(elec, :, :)));
            elec = setdiff(elec, flat, 'stable');
            if ~isempty(flat)
                pipecompare.utils.log('Flat channel(s) (no signal) marked bad: %s', strjoin({EEG.chanlocs(flat).labels}, ', '));
            end
            % kurtosis and joint probability on a copy without the blinks:
            % they make the frontal channels' values heavy-tailed, so
            % those would be marked bad for what ICA is there to remove
            noBlink = src; share = 0;
            if any(~strcmp(measures, 'resid'))
                if copy && EEG.trials == 1
                    % left out of the low-passed data, then high-passed like
                    % the copy: high-passing first would spread each blink's
                    % slow part beyond the stretch left out (eegfiltnew
                    % filters the stretches between boundaries on their
                    % own), while low-passing after would leave what it
                    % removes (line noise) at the edges of every stretch
                    [noBlink, share, eyes] = pipecompare.run.Steps.withoutBlinks(src, lpOnly);
                    if share > 0, noBlink = bandCopy(noBlink, doHp, hp, false, lp); else, noBlink = src; end
                else
                    [noBlink, share, eyes] = pipecompare.run.Steps.withoutBlinks(src);
                end
                if share > 0
                    on = sprintf('%s; %.0f%% of the data with blinks left out for kurtosis and probability', on, 100 * share);
                    pipecompare.utils.log('Bad channels: kurtosis and joint probability measured without the %.0f%% of the data with blinks (in %s).', ...
                        100 * share, strjoin({src.chanlocs(eyes).labels}, ', '));
                end
            end
            bad = [];
            for m = measures
                if strcmp(m{1}, 'resid')
                    z = pipecompare.run.Steps.residualZ(src.data(elec, :, :));
                    if isempty(z)
                        pipecompare.utils.log(['Bad channels: the noise the other channels do not explain is not ', ...
                            'tested (the channels are linearly dependent, e.g. already average-referenced or interpolated).']);
                        continue
                    end
                    b = find(z > p.threshold);
                    if ~isempty(b)
                        pipecompare.utils.log('Bad channel(s) by their noise the other channels do not explain: %s', ...
                            strjoin(arrayfun(@(k) sprintf('%s (z = %.1f)', EEG.chanlocs(elec(k)).labels, z(k)), b, ...
                            'UniformOutput', false), ', '));
                    end
                else
                    args = {'elec', elec, 'threshold', p.threshold, 'norm', 'on', 'measure', m{1}};
                    if strcmp(m{1}, 'spec'), args = [args {'freqrange', [1 min(50, EEG.srate/2 - 1)]}]; end
                    [EEGrem, b, ~, com] = pop_rejchan(noBlink, args{:});
                end
                bad = union(bad, b(:)');
            end
            bad = elec(bad);   % both index the tested channels (pop_rejchan removes opt.elec(indelec) itself)
            on = [on pipecompare.utils.ternary(isempty(flat), '', sprintf('; flat channels [%s] bad', ...
                strjoin({EEG.chanlocs(flat).labels}, ' ')))];
            bad = sort([bad flat]);
            labels = {EEG.chanlocs(bad).labels};
            info.badChannels = labels; info.badIdx = bad;
            if strcmp(p.action, 'remove')
                info.removed = labels;
                if ~copy && share == 0 && isscalar(measures) && ~strcmp(measures{1}, 'resid') && isempty(flat)
                    EEG = EEGrem;
                    coms = {com};
                    return;
                end
                coms = {sprintf('%% PipeCompare: bad-channel test (measure %s, threshold %g, norm on, %d channels tested)%s -> bad channels [%s] removed', ...
                    p.measure, p.threshold, numel(elec), on, strjoin(labels, ' '))};
                if ~isempty(bad)
                    [EEG, com2] = pop_select(EEG, 'nochannel', bad);
                    coms{end+1} = com2;
                end
                return;
            end
            done = sprintf('[%s] interpolated in place', strjoin(labels, ' '));
            if strcmp(pipecompare.utils.fieldOr(p, 'repair', 'all'), 'error') && ~isempty(bad)
                % interpolated only where that does not measurably raise
                % the total error (pipecompare.eval.ChannelChoice)
                if nargin > 2 && isfield(ctx, 'chosenChannels')
                    % re-running a pipeline: the channels the search interpolated
                    keep = bad(ismember(lower(labels), lower(ctx.chosenChannels)));
                elseif nargin > 2 && isfield(ctx, 'after') && isfield(ctx, 'contract') && pipecompare.eval.Threshold.applies(ctx.contract)
                    [keep, info.channelChoice] = pipecompare.eval.ChannelChoice.choose(EEG, src, bad, ctx);
                else
                    keep = bad;
                    pipecompare.utils.log(['Bad channels: every bad channel interpolated (whether interpolating one would ', ...
                        'raise the total error needs mean amplitudes or band power and the steps that follow).']);
                end
                if numel(keep) < numel(bad)
                    done = sprintf(['[%s]: [%s] interpolated in place, [%s] left as they are (interpolating them ', ...
                        'would raise the total error)'], strjoin(labels, ' '), strjoin({EEG.chanlocs(keep).labels}, ' '), ...
                        strjoin({EEG.chanlocs(setdiff(bad, keep)).labels}, ' '));
                end
                info.leftAsIs = {EEG.chanlocs(setdiff(bad, keep)).labels};
                bad = keep; info.badIdx = bad;
            end
            coms = {sprintf('%% PipeCompare: bad-channel test (measure %s, threshold %g, norm on, %d channels tested)%s -> bad channels %s', ...
                p.measure, p.threshold, numel(elec), on, done)};
            info.interpolated = {EEG.chanlocs(bad).labels};
            if ~isempty(bad)
                requireLocations(EEG.chanlocs(bad), 'badchannels');
                EEG = pop_interp(EEG, bad, 'spherical');   % returns no command from the command line
                coms{end+1} = sprintf('EEG = pop_interp(EEG, %s, ''spherical'');', mat2str(bad(:)'));
            end
            if pipecompare.utils.fieldOr(p, 'stretches', false) && EEG.trials == 1
                % channels bad for only part of the recording: found on the
                % same copy among the channels not found bad as a whole, and
                % interpolated in those stretches only
                good = setdiff(elec, find(ismember({EEG.chanlocs.labels}, labels)), 'stable');
                st = pipecompare.run.Steps.badStretches(src.data(good, :), EEG.srate);
                st = struct('channel', arrayfun(@(r) EEG.chanlocs(good(r)).labels, [st.row], 'UniformOutput', false), ...
                    'samples', {st.samples});
                if ~isempty(st)
                    requireLocations(EEG.chanlocs(ismember({EEG.chanlocs.labels}, {st.channel})), 'badchannels');
                    EEG = pipecompare.run.Steps.interpStretches(EEG, st);
                    coms{end+1} = stretchCommand(st);
                    pipecompare.utils.log('Channels bad for part of the recording, interpolated there only: %s', ...
                        strjoin(arrayfun(@(q) sprintf('%s %.0f-%.0f s', q.channel, (q.samples - [1 0]) / EEG.srate), st, ...
                        'UniformOutput', false), ', '));
                end
                info.stretches = st;
            end
        end

        function st = badStretches(X, srate)
            % Channels (rows of X, continuous, srate Hz) bad for only part of
            % the recording. Each channel's noise the others do not explain
            % (its residual from their least-squares prediction, as in
            % residualZ: (P x)_i / P_ii for P the inverse covariance) has its
            % log SD taken in 1 s windows; what all channels share in a
            % window (movement, a change of reference) is taken out (the
            % median over the channels). On that sequence a two-state HMM
            % (pipecompare.eval.Hmm.segment: EM, Viterbi) finds the stretches
            % in its higher state; one counts when it lasts MinWindows (10 s)
            % or more and its mean is StretchZ (3) robust SDs above the
            % channel's usual level. That level: the channel's 20th
            % percentile plus the median channel's gap between its 20th
            % percentile and its median, and the spread the median
            % channel's MAD, so a channel bad for most of the recording is
            % still measured against its good part. st: row (of X) and
            % samples ([first last]) of each stretch; none with fewer than 8
            % channels, 30 windows, or linearly dependent channels.
            st = struct('row', {}, 'samples', {});
            [nc, n] = size(X); w = round(srate); nw = floor(n / w);
            if nc < 8 || nw < 30, return; end
            Xs = double(X(:, 1:max(1, ceil(n / 2e5)):end)); Xs = Xs - median(Xs, 2);
            mad = 1.4826 * median(abs(Xs), 2);
            keep = max(abs(Xs) ./ max(mad, eps), [], 1) < 6;
            if mean(keep) >= 0.5, Xs = Xs(:, keep); end
            Xs = Xs - mean(Xs, 2);
            C = (Xs * Xs') / size(Xs, 2);
            if rcond(C) < 1e-10
                pipecompare.utils.log(['Bad channels: stretches of bad channels not looked for (the channels are ', ...
                    'linearly dependent, e.g. already average-referenced or interpolated).']);
                return;
            end
            P = inv(C); Wr = P ./ diag(P);
            L = zeros(nc, nw);
            for k = 1:nw
                L(:, k) = log(std(Wr * double(X(:, (k - 1) * w + (1:w))), 0, 2) + eps);
            end
            L = L - median(L, 1);
            q20 = zeros(nc, 1); gap = zeros(nc, 1); sc = zeros(nc, 1);
            for c = 1:nc
                v = sort(L(c, :)); q20(c) = v(ceil(0.2 * nw)); m = median(v);
                gap(c) = m - q20(c); sc(c) = 1.4826 * median(abs(v - m));
            end
            scale = max(0.05, median(sc));
            for c = 1:nc
                y = (L(c, :) - q20(c) - median(gap)) / scale;
                on = pipecompare.eval.Hmm.segment(y, -inf, pipecompare.run.Steps.MinWindows);
                d = diff([0 on 0]); s = find(d == 1); e = find(d == -1) - 1;
                for q = 1:numel(s)
                    if mean(y(s(q):e(q))) >= pipecompare.run.Steps.StretchZ
                        last = e(q) * w; if e(q) == nw, last = n; end
                        st(end + 1) = struct('row', c, 'samples', [(s(q) - 1) * w + 1, last]); %#ok<AGROW>
                    end
                end
            end
        end

        function EEG = interpStretches(EEG, st)
            % The channels st(k).channel (labels) interpolated (spherical, as
            % pop_interp) in the samples st(k).samples ([first last]) only,
            % from the channels not bad there; continuous data. A fixed
            % linear operation once the stretches are known, so the signal
            % check replays it exactly.
            if isempty(st), return; end
            [~, ch] = ismember(lower({st.channel}), lower({EEG.chanlocs.labels}));
            S = reshape([st.samples], 2, [])';
            edges = unique([S(:, 1); S(:, 2) + 1; 1; EEG.pnts + 1])';
            done = containers.Map('KeyType', 'char', 'ValueType', 'any');
            for k = 1:numel(edges) - 1
                a = edges(k); b = edges(k + 1) - 1;
                bad = unique(ch(S(:, 1)' <= a & S(:, 2)' >= b));
                if isempty(bad) || a > b, continue; end
                key = mat2str(bad);
                if ~isKey(done, key), done(key) = interpWeights(EEG, bad); end
                W = done(key);
                for q = a:1e5:b
                    cols = q:min(b, q + 1e5 - 1);
                    EEG.data(bad, cols) = cast(W * double(EEG.data(:, cols)), 'like', EEG.data);
                end
            end
        end

        function [T, share, eyes] = withoutBlinks(T, D)
            % T without its blinks: the samples (continuous data) or epochs
            % where an eye channel (Injection.eyeChannels: EOG, else Fp1,
            % Fp2, FPz) is beyond 5 robust SDs of its own (median and MAD),
            % continuous data with 200 ms on either side. With D (the same
            % samples, e.g. T before filtering), those are left out of D
            % instead. share: the share of the data left out; nothing is
            % left out (share 0) without an eye channel, or when that would
            % be over half the data.
            %   A blink too weak to stand out in an eye channel's noise
            %   is clear in the frontal channels, and one blink left in is
            %   enough to make them heavy-tailed. So the blinks found in
            %   the eye channels give the blink's scalp pattern a
            %   (regression of every channel on the eye channel with most
            %   of them), and the data weighted by C \ a, C their
            %   covariance away from the blinks, give its time course: the
            %   matched spatial filter, in which a noisy channel weighs
            %   little, so a bad frontal channel's own noise or spikes are
            %   not taken for blinks. Its samples beyond 5 robust SDs are
            %   blinks too.
            if nargin > 1, S = D; else, S = T; end
            share = 0;
            eyes = pipecompare.eval.Injection.eyeChannels(T.chanlocs);
            if isempty(eyes), return; end
            Y = reshape(T.data, T.nbchan, []);
            H = beyondRobust(double(Y(eyes, :)), 5);
            if T.trials == 1   % weak blinks too: the HMM of the eye channels (Hmm.events)
                H(end + 1, :) = pipecompare.eval.Hmm.events(Y(eyes, :), T.srate);
            end
            hit = any(H, 1);
            w = round(0.2 * T.srate);
            if any(hit)
                [~, j] = max(sum(H(1:numel(eyes), :), 2));
                idx = find(hit); e = double(Y(eyes(j), idx)); a = zeros(T.nbchan, 1);
                for b = 1:1e5:numel(idx)           % a block of samples at a time
                    k = b:min(numel(idx), b + 1e5 - 1);
                    a = a + double(Y(:, idx(k))) * e(k)';
                end
                a = a / (e * e');
                away = ~hit;
                if T.trials == 1, away = ~(conv(double(hit), ones(1, 2 * w + 1), 'same') > 0); end
                C = blockCov(Y, away);
                C = C + 1e-6 * max(trace(C), eps) / size(C, 1) * eye(size(C, 1));
                f = C \ a;
                hit = hit | beyondRobust(double(single(f') * Y) / (f' * a), 5);
            end
            if T.trials > 1
                drop = find(any(reshape(hit, T.pnts, T.trials), 1));
                if isempty(drop) || numel(drop) > T.trials / 2, return; end
                share = numel(drop) / T.trials;
                [~, T] = evalc('pop_select(S, ''notrial'', drop)');
                return;
            end
            hit = conv(double(hit), ones(1, 2 * w + 1), 'same') > 0;
            if ~any(hit) || mean(hit) > 0.5, return; end
            d = diff([0 hit 0]);
            T = pipecompare.utils.selectPoints(S, 'nopoint', [find(d == 1)' find(d == -1)' - 1]);
            share = mean(hit);
        end

        function z = residualZ(X)
            % For each row of X (channels x samples [x epochs]): how far its
            % noise that the other rows do not explain stands out, as a
            % robust z (median and MAD over the rows, MAD at least 0.1) of
            % the log SD of what is left of it after its least-squares
            % prediction from all the others (leave one out:
            % SD^2 / (C^-1)_ii for the covariance C). Brain activity and
            % artifacts reach neighbouring electrodes too (volume
            % conduction), so the others predict most of a good channel;
            % noise of its own (poor contact, a broken lead) they cannot.
            % Size differences the others share (alpha, blinks, a nearby
            % reference) are predicted, so they do not count. Samples where
            % any row is beyond 6 robust SDs of its own (blinks, movement)
            % are left out, unless that is over half of them. [] when the
            % rows are linearly dependent (average-referenced or
            % interpolated): every row is then predicted exactly. Fewer
            % than 8 rows: no outlier can be told, z is all zeros. Epochs:
            % each epoch's mean is removed first.
            n = size(X, 1);
            z = zeros(1, n);
            if n < 8, return; end
            if size(X, 3) > 1, X = X - mean(X, 2); end              % epochs: each one's offset removed
            X = reshape(X, n, []);
            X = double(X(:, 1:max(1, ceil(size(X, 2) / 2e5)):end));   % at most 2e5 samples
            X = X - median(X, 2);
            mad = 1.4826 * median(abs(X), 2);
            keep = max(abs(X) ./ max(mad, eps), [], 1) < 6;
            if mean(keep) >= 0.5, X = X(:, keep); end
            X = X - mean(X, 2);
            C = (X * X') / size(X, 2);
            sd = sqrt(diag(C));
            R = C ./ (sd * sd');
            ev = eig((R + R') / 2);
            if min(ev) < 1e-6 * max(ev), z = []; return; end
            r = log(sd' ./ sqrt(diag(inv(R))'));   % the residual SD: sd * sqrt(1 - R^2)
            z = (r - median(r)) / max(0.1, 1.4826 * median(abs(r - median(r))));
        end

        function flat = flatChannels(X)
            % Which rows of X (channels x samples [x epochs]) are flat: a
            % spread (SD) of zero, or under 1% of the median spread of the
            % rows (a disconnected input records almost nothing).
            X = reshape(double(X), size(X, 1), []);
            sd = std(X, 0, 2)';
            flat = sd == 0 | sd < 0.01 * median(sd);
        end

        function [EEG, coms, info] = listedChannels(EEG, p)
            % Channels you name (e.g. known-bad O1/O2): remove or interpolate.
            labels = cellstr(p.labels);
            [ok, idx] = ismember(lower(labels), lower({EEG.chanlocs.labels}));
            assert(all(ok), 'PipeCompare:Channels', 'Channel(s) not in the data: %s', strjoin(labels(~ok), ', '));
            info.listedIdx = idx(:)';
            switch p.action
                case 'remove'
                    [EEG, com] = pop_select(EEG, 'nochannel', labels);
                    coms = {com}; info.removed = labels;
                case 'interpolate'
                    requireLocations(EEG.chanlocs(idx), 'channels');
                    EEG = pop_interp(EEG, idx(:)', 'spherical');   % returns no command from the command line
                    coms = {sprintf('EEG = pop_interp(EEG, %s, ''spherical'');', mat2str(idx(:)'))}; info.interpolated = labels;
                otherwise
                    error('PipeCompare:Channels', 'channels action must be remove or interpolate');
            end
        end

        function [EEG, coms] = ica(EEG, p, ctx)
            % 'rndreset','no' makes runica start from a fixed random state,
            % so the same data give the same decomposition. Adopting a
            % candidate (replay must reproduce the evaluated scores) relies
            % on this; it holds for runica in EEGLAB 2026.0 (checked by
            % test_engine/testEndToEndRecoversSensibleChoice and the adopt
            % replay check, which refuses a candidate that differs).
            opts = {'icatype', 'runica', 'extended', p.extended, 'rndreset', 'no', 'interrupt', 'off'};
            % pop_runica returns its command only from its dialog, so the
            % history line is written here
            args = vararg2str(opts);
            applied = 0; if isfield(ctx, 'highpass'), applied = ctx.highpass; end
            fitHp = p.fitHighpass > 0 && p.fitHighpass > applied;
            clean = pipecompare.utils.fieldOr(p, 'fitClean', 1) > 0;
            % ctx.icaSeed > 0: runica sees the samples in another order
            % (the check of how much the recommendation depends on where
            % ICA starts; runica resets its own random state, so the order
            % of the samples is what changes its path). Not adopted.
            seed = pipecompare.utils.fieldOr(ctx, 'icaSeed', 0);
            if fitHp || clean || seed > 0
                % Fit on a copy: high-passed (stable decomposition) and/or
                % without the stretches of extreme data (movement, bursts),
                % which would otherwise claim components of their own; the
                % weights are applied to the data as they are (ERP signal kept).
                pre = {}; what = {};
                if fitHp
                    [fa, fb] = fftArgs(EEG, p.fitHighpass); %#ok<ASGLU>
                    [tmp, ~] = pop_eegfiltnew(EEG, 'locutoff', p.fitHighpass, 'plotfreqz', 0, fa{:});
                    fftTxt = ''; if ~isempty(fa) && isempty(fb), fftTxt = ', ''usefftfilt'', 1'; end
                    pre{end+1} = sprintf('EEGica = pop_eegfiltnew(EEG, ''locutoff'', %g, ''plotfreqz'', 0%s);', p.fitHighpass, fftTxt);
                    what{end+1} = sprintf('a %g Hz high-passed copy', p.fitHighpass);
                else
                    tmp = EEG;
                    pre{end+1} = 'EEGica = eeg_checkset(EEG);';
                end
                if clean
                    [tmp, left] = pipecompare.run.Steps.cleanForIca(tmp);
                    pre{end+1} = 'EEGica = pipecompare.run.Steps.cleanForIca(EEGica);';
                    if left > 0
                        what{end+1} = sprintf('%.0f%% of the data left out as extreme', 100 * left);
                        pipecompare.utils.log('ICA fitted without the most extreme %.0f%% of the data.', 100 * left);
                    end
                end
                if isempty(what), what = {'the data as they are (nothing extreme to leave out)'}; end
                if seed > 0
                    sz = size(tmp.data); x = reshape(tmp.data, sz(1), []);
                    x = x(:, randperm(RandStream('mt19937ar', 'Seed', seed), size(x, 2)));
                    tmp.data = reshape(x, sz); clear x
                    what{end+1} = sprintf('its samples in a random order (seed %d)', seed);
                end
                tmp = runIca(tmp, opts);
                EEG.icaweights = tmp.icaweights; EEG.icasphere = tmp.icasphere;
                EEG.icachansind = tmp.icachansind; EEG.icawinv = []; EEG.icaact = [];
                EEG = eeg_checkset(EEG);
                c2 = sprintf('EEGica = pop_runica(EEGica, %s);', args);
                coms = {sprintf(['%s %s ', ...
                    'EEG.icaweights = EEGica.icaweights; EEG.icasphere = EEGica.icasphere; ', ...
                    'EEG.icachansind = EEGica.icachansind; EEG.icawinv = []; EEG.icaact = []; ', ...
                    'EEG = eeg_checkset(EEG); clear EEGica; %% PipeCompare: ICA fitted on %s'], ...
                    strjoin(pre, ' '), c2, strjoin(what, ', '))};
            else
                EEG = runIca(EEG, opts);
                coms = {sprintf('EEG = pop_runica(EEG, %s);', args)};
            end
        end

        function [T, left] = cleanForIca(T)
            % The copy ICA is fitted on, without its most extreme stretches:
            % 1 s windows of continuous data (or epochs) whose spread is far
            % above the others' (robust z of the mean log SD over the EEG
            % channels > 3, median and MAD over the windows), at most 20 %
            % of them, worst first. left: the share of the data left out.
            % EOG and other non-EEG channels do not count (blinks stay in, so
            % that ICA learns them); every channel stays in the copy.
            eeg = ~pipecompare.simple.Presets.isNonEeg(T.chanlocs(:)');
            if T.trials > 1
                n = T.trials; X = double(T.data(eeg, :, :));
                s = squeeze(mean(log(std(X, 0, 2) + eps), 1))';
            else
                w = round(T.srate); n = floor(T.pnts / w);
                s = zeros(1, n);
                for k = 1:n
                    s(k) = mean(log(std(double(T.data(eeg, (k - 1) * w + (1:w))), 0, 2) + eps));
                end
            end
            left = 0;
            if n < 10, return; end
            mad = max(1.4826 * median(abs(s - median(s))), eps);
            z = (s - median(s)) / mad;
            [~, order] = sort(z, 'descend');
            drop = order(z(order) > 3); drop = drop(1:min(numel(drop), floor(0.2 * n)));
            if isempty(drop), return; end
            if T.trials > 1
                [~, T] = evalc('pop_select(T, ''notrial'', drop)');
                left = numel(drop) / n;
            else
                drop = sort(drop);
                T = pipecompare.utils.selectPoints(T, 'nopoint', [(drop(:) - 1) * w + 1, drop(:) * w]);
                left = numel(drop) * w / (n * w);
            end
        end

        function [EEG, coms, info] = icRemove(EEG, p, ctx)
            % ctx.iclabel, when set, is the ICLabel classification of these
            % same data (a sibling step with another threshold made it), so
            % the network does not run again.
            assert(exist('pop_iclabel', 'file') == 2, 'PipeCompare:Dependency', 'ICLabel plugin not installed');
            if nargin > 2 && isfield(ctx, 'iclabel') && ~isempty(ctx.iclabel)
                EEG.etc.ic_classification = ctx.iclabel.classification; c1 = ctx.iclabel.com;
            else
                warnIclabel15();
                [EEG, c1] = pop_iclabel(EEG, 'default');
            end
            info.iclabel = struct('classification', EEG.etc.ic_classification, 'com', c1);
            cats = {'Brain','Muscle','Eye','Heart','Line Noise','Channel Noise','Other'};
            classes = cellstr(p.classes);
            assert(all(ismember(classes, cats)) && ~ismember('Brain', classes), 'PipeCompare:ICLabel', ...
                'classes must be ICLabel artifact classes: %s', strjoin(cats(2:end), ', '));
            if strcmp(p.threshold, 'sme')
                % the components chosen by their effect on the SME
                % (pipecompare.eval.IcThreshold.chooseBySme), never one whose
                % likeliest ICLabel class is Brain
                L = EEG.etc.ic_classification.ICLabel;
                [~, top] = max(L.classifications, [], 2);
                cand = find(~strcmp(cellstr(L.classes(top)), 'Brain'))';
                if isfield(ctx, 'chosenComps')
                    comps = ctx.chosenComps;   % re-running a pipeline: the components the search chose
                else
                    [comps, info.icSme] = pipecompare.eval.IcThreshold.chooseBySme(EEG, cand, ctx);
                    info.icsProtected = numel(info.icSme.protected);   % kept: time-locked activity in a measured window
                    info.curves = info.icSme.curves;   % for the ranking's selection check (Rank.crossfit)
                end
                nIc = size(EEG.icaweights, 1);
                EEG.reject.gcompreject = false(1, nIc); EEG.reject.gcompreject(comps) = true;
                c2 = sprintf(['EEG.reject.gcompreject = false(1, %d); EEG.reject.gcompreject(%s) = true; ', ...
                    '%% PipeCompare: components chosen by their effect on the SME (none whose likeliest ICLabel class is Brain)'], ...
                    nIc, mat2str(comps));
                info.icsConsidered = numel(cand);
            else
                if ischar(p.threshold)
                    % 'auto': the threshold chosen from these data
                    % (pipecompare.eval.IcThreshold), then applied as a fixed one
                    L = EEG.etc.ic_classification.ICLabel;
                    P = L.classifications(:, ismember(cellstr(L.classes), classes));
                    P(P >= 1) = 0;   % pop_icflag flags t < p < 1
                    q = max(P, [], 2);
                    [p.threshold, info.icThreshold] = pipecompare.eval.IcThreshold.choose(EEG, q, ctx);
                    info.thresholdChosen = p.threshold;
                    info.curves = info.icThreshold.curves;
                    info.icsConsidered = info.icThreshold.sets - 1;   % components over IcThreshold.MinProbability
                end
                T = nan(7, 2);
                T(ismember(cats, classes), :) = repmat([p.threshold 1], sum(ismember(cats, classes)), 1);
                [EEG, c2] = pop_icflag(EEG, T);
            end
            comps = find(EEG.reject.gcompreject(:)');
            info.icsTotal = size(EEG.icaweights, 1);
            info.icsRemoved = numel(comps);
            info.icaPoints = EEG.pnts * EEG.trials;   % the data points ICA had (the ICA check)
            % what ICLabel recognised: components it takes for brain
            % activity (Brain >= 0.5) and those whose likeliest class is
            % Other, and the median Other probability (the ICA check)
            L = EEG.etc.ic_classification.ICLabel;
            names = cellstr(L.classes); Pc = L.classifications;
            iB = strcmp(names, 'Brain'); iO = strcmp(names, 'Other');
            if any(iB) && any(iO) && ~isempty(Pc)
                [~, top] = max(Pc, [], 2);
                info.icsBrain = sum(Pc(:, iB) >= 0.5);
                info.icsOther = sum(top == find(iO));
                info.otherMedian = median(Pc(:, iO));
            end
            info.comps = comps;
            coms = {c1, c2};
            if ~isempty(comps)
                [EEG, c3] = pop_subcomp(EEG, comps, 0);
                coms{end+1} = c3;
            end
        end

        function [EEG, coms, info] = rejectEpochs(EEG, type, p, ctx)
            n0 = EEG.trials;
            cInterp = '';
            chans = pipecompare.run.Steps.testedChannels(EEG, p);
            switch type
                case 'reject_threshold'
                    pp = strcmp(pipecompare.utils.fieldOr(p, 'method', 'absolute'), 'peaktopeak');
                    win = pipecompare.utils.fieldOr(p, 'window', 200);
                    if ischar(p.uv)
                        % 'auto': the limit chosen from these data
                        % (pipecompare.eval.Threshold), then applied as a fixed one
                        m = pipecompare.run.Steps.limitStat(EEG, chans, p);
                        [uv, th] = pipecompare.eval.Threshold.choose(EEG, m, ctx);
                        p.uv = uv; info.threshold = th; info.uvChosen = uv;
                        info.curves = th.curves;
                    end
                    if pp
                        EEG = pipecompare.run.Steps.markPeakToPeak(EEG, chans, p.uv, win);
                        c1 = sprintf('EEG = pipecompare.run.Steps.markPeakToPeak(EEG, %s, %g, %g);', mat2str(chans), p.uv, win);
                    else
                        [EEG, ~, c1] = pop_eegthresh(EEG, 1, chans, -p.uv, p.uv, EEG.xmin, EEG.xmax, 0, 0);
                    end
                    marks = EEG.reject.rejthresh; E = pipecompare.utils.fieldOr(EEG.reject, 'rejthreshE');
                case 'reject_jointprob'
                    [EEG, ~, ~, ~, c1] = pop_jointprob(EEG, 1, chans, p.sd, p.sd, 0, 0, 0, [], 0);
                    marks = EEG.reject.rejjp; E = pipecompare.utils.fieldOr(EEG.reject, 'rejjpE');
                case 'reject_kurtosis'
                    [EEG, ~, ~, ~, c1] = pop_rejkurt(EEG, 1, chans, p.sd, p.sd, 0, 0, 0, [], 0);
                    marks = EEG.reject.rejkurt; E = pipecompare.utils.fieldOr(EEG.reject, 'rejkurtE');
            end
            info.epochsInterpolated = 0; info.epochsNotRepaired = 0;
            info.epochInterp = struct('chans', {{}}, 'epochs', {{}});
            nMax = pipecompare.utils.fieldOr(p, 'interpolate', 0);
            if nMax > 0 && ~isempty(E) && size(E, 1) == EEG.nbchan
                % Epochs failed by at most nMax channels keep them, interpolated
                % within the epoch (spherical splines from the other
                % channels), instead of being rejected; the decision is the
                % test's own per-channel marks, taken before interpolating.
                flagged = E ~= 0;
                nf = sum(flagged, 1);
                ok = @(v) isnumeric(v) && isscalar(v) && isfinite(v);
                located = arrayfun(@(c) isfield(c, 'X') && ok(c.X) && ok(c.Y) && ok(c.Z), EEG.chanlocs(:)');
                % never a measured electrode: an interpolated value is a
                % weighted mean of its neighbours, with less noise than a real
                % electrode, so repairing it would make the SME look better
                % than it is; such an epoch is rejected
                measured = false(1, EEG.nbchan);
                if nargin > 3 && isfield(ctx, 'contract') && ~isempty(ctx.contract)
                    measured = ismember(lower({EEG.chanlocs.labels}), lower(ctx.contract.allRoi()));
                end
                kept = find(marks(:)' & nf >= 1 & nf <= nMax & ~any(flagged & ~located(:), 1) & ~any(flagged & measured(:), 1));
                info.epochsNotRepaired = sum(marks(:)' & nf >= 1 & nf <= nMax & any(flagged & measured(:), 1));
                if info.epochsNotRepaired > 0
                    pipecompare.utils.log('%d epoch(s) not repaired: an electrode you measure was over the limit in them.', ...
                        info.epochsNotRepaired);
                end
                if ~isempty(kept)
                    [sigs, ~, g] = unique(arrayfun(@(e) mat2str(find(flagged(:, e))'), kept, 'UniformOutput', false));
                    sets = cellfun(@str2num, sigs(:)', 'UniformOutput', false); %#ok<ST2NM> (mat2str of index vectors)
                    epochs = arrayfun(@(k) kept(g == k), 1:numel(sigs), 'UniformOutput', false);
                    EEG = pipecompare.run.Steps.interpolateEpochs(EEG, sets, epochs);
                    info.epochInterp = struct('chans', {sets}, 'epochs', {epochs});
                    info.epochsInterpolated = numel(kept);
                    field = rejField(type);
                    EEG.reject.(field)(kept) = 0; EEG.reject.([field 'E'])(:, kept) = 0;   % kept: no longer marked
                    marks(kept) = 0; E(:, kept) = 0;
                    cInterp = sprintf(['EEG = pipecompare.run.Steps.interpolateEpochs(EEG, %s, %s); ', ...
                        'EEG.reject.%s(%s) = 0; EEG.reject.%sE(:, %s) = 0; %% PipeCompare: %d epoch(s) kept with ', ...
                        'up to %d flagged channel(s) interpolated in them'], cellText(sets), cellText(epochs), ...
                        field, mat2str(kept), field, mat2str(kept), numel(kept), nMax);
                    pipecompare.utils.log('%d epoch(s) kept by interpolating the 1 to %d channel(s) flagged in them.', numel(kept), nMax);
                end
            end
            idx = find(marks);
            % which channels drive the rejections (a hint for bad channels)
            info.topChannels = '';
            info.overLimit = struct('labels', {{}}, 'counts', []);
            if ~isempty(E) && size(E, 1) == EEG.nbchan && ~isempty(idx)
                n = sum(E(:, idx) ~= 0, 2)';
                info.overLimit = struct('labels', {{EEG.chanlocs(n > 0).labels}}, 'counts', n(n > 0));
                [cnt, order] = sort(n', 'descend');
                keep = order(cnt > 0); keep = keep(1:min(5, numel(keep)));
                info.topChannels = strjoin(arrayfun(@(c) sprintf('%s (%d)', EEG.chanlocs(c).labels, ...
                    sum(E(c, idx) ~= 0)), keep(:)', 'UniformOutput', false), ', ');
                pipecompare.utils.log('%d/%d epochs rejected; channels most often over the limit: %s', ...
                    numel(idx), n0, info.topChannels);
            end
            % whether the epochs kept differ in the measured window from all
            % of them (a low limit can reject epochs for their large real
            % activity there, which makes the average smaller)
            info.keptShift = [];
            if ~isempty(idx) && numel(idx) < n0 && nargin > 3 && isfield(ctx, 'contract') && ~isempty(ctx.contract) && ...
                    strcmp(ctx.contract.analysis, 'erp')
                info.keptShift = keptShift(EEG, idx, ctx.contract);
            end
            coms = {c1, cInterp};
            if ~isempty(idx)
                assert(numel(idx) < n0, 'PipeCompare:AllRejected', 'every epoch would be rejected');
                [EEG, c2] = pop_rejepoch(EEG, idx, 0);
                coms{end+1} = c2;
            end
            info.rejected = numel(idx); info.rejIdx = idx; info.epochsBefore = n0;
        end

        function [EEG, info] = eogRegress(EEG, eyes, groups, epoch)
            % The eye channels regressed out of the EEG channels (blinks and
            % eye movements spread to the scalp in proportion to them):
            % EEG.data(t, :) - B * EEG.data(e, :), B (targets x eyes) the
            % least-squares fit of the EEG channels on the eye channels.
            % The fit uses only the samples where an eye channel is beyond 4
            % robust SD of its median or in the active state of the eye
            % channels' HMM (pipecompare.eval.Hmm.events) (+-100 ms), where
            % the eyes dominate,
            % and the average of the epochs around each event type of
            % groups (cells of event types; epoch [t0 t1] s) is taken out
            % of both first, so that what is time-locked (the ERP, also seen
            % at the eye channels) does not enter B. eyes: labels; {} = the
            % dedicated eye channels (EOG by type or name, not Fp1/Fp2:
            % pipecompare.eval.Injection.dedicatedEyes). Targets: every
            % channel that is EEG (pipecompare.simple.Presets.isNonEeg).
            % info: eyes, targets (labels), B, fitted (samples used), seconds.
            labs = {EEG.chanlocs.labels};
            if isempty(eyes), ei = pipecompare.eval.Injection.dedicatedEyes(EEG.chanlocs);
            else, [~, ei] = ismember(lower(cellstr(eyes)), lower(labs)); ei = ei(ei > 0); end
            assert(~isempty(ei), 'PipeCompare:Eog', ['Regressing the eye channels out needs an eye channel of its ', ...
                'own (EOG, HEOG, VEOG, EYE...; type EOG): Fp1/Fp2 are EEG channels and would be removed with it.']);
            ti = setdiff(find(~pipecompare.simple.Presets.isNonEeg(EEG.chanlocs)), ei);
            n = EEG.pnts * EEG.trials;
            D = reshape(EEG.data, EEG.nbchan, n);
            % where each sample lies in the epochs around the events (its
            % group and place in the epoch; 0 outside), for the average
            if EEG.trials > 1
                grp = ones(1, n); pos = repmat(1:EEG.pnts, 1, EEG.trials); nG = 1;
                epochsOf = {reshape(1:n, EEG.pnts, EEG.trials)'};
            else
                grp = zeros(1, n); pos = zeros(1, n); nG = 0; epochsOf = {};
                if ~isempty(groups) && ~isempty(EEG.event)
                    w = round(epoch(1) * EEG.srate):round(epoch(2) * EEG.srate);
                    types = cellfun(@(x) char(string(x)), {EEG.event.type}, 'UniformOutput', false);
                    for g = 1:numel(groups)
                        lat = round([EEG.event(ismember(types, cellstr(groups{g}))).latency]);
                        lat = lat(lat + w(1) >= 1 & lat + w(end) <= n);
                        if isempty(lat), continue; end
                        nG = nG + 1; epochsOf{nG} = lat(:) + w;   %#ok<AGROW> % epochs x samples
                        grp(epochsOf{nG}) = nG; pos(epochsOf{nG}) = repmat(1:numel(w), numel(lat), 1);
                    end
                end
            end
            % the time-locked average of every channel used, per group
            avg = cell(1, nG); rows = [ti ei];
            for g = 1:nG
                ix = epochsOf{g};
                avg{g} = mean(reshape(double(D(rows, ix')), numel(rows), size(ix, 2), size(ix, 1)), 3);
            end
            noErp = @(cols) double(D(rows, cols)) - erpAt(avg, grp(cols), pos(cols), numel(rows));
            % the samples dominated by the eyes, on the eye channels without the ERP
            E = zeros(numel(ei), n);
            for q = 1:1e5:n
                cols = q:min(n, q + 1e5 - 1); Z = noErp(cols); E(:, cols) = Z(numel(ti)+1:end, :);
            end
            z = abs(E - median(E, 2)) ./ (1.4826 * median(abs(E - median(E, 2)), 2) + eps);
            half = round(0.1 * EEG.srate);
            on = any(z > 4, 1);
            if EEG.trials == 1, on = on | pipecompare.eval.Hmm.events(E, EEG.srate); end   % weak blinks too
            sel = conv(double(on), ones(1, 2 * half + 1), 'same') > 0;
            if sum(sel) < EEG.srate, sel = true(1, n); end   % less than 1 s of eye activity: all samples
            clear E z
            % B from the cross-products over those samples (blocks: memory)
            idx = find(sel); m = numel(rows); C = zeros(m); s1 = zeros(m, 1);
            for q = 1:1e5:numel(idx)
                Z = noErp(idx(q:min(end, q + 1e5 - 1)));
                C = C + Z * Z'; s1 = s1 + sum(Z, 2);
            end
            C = C / numel(idx) - (s1 / numel(idx)) * (s1 / numel(idx))';   % covariance
            k = numel(ti);
            B = C(1:k, k+1:end) * pinv(C(k+1:end, k+1:end));   % (two eye channels can be nearly the same)
            EEG = applyEog(EEG, ti, ei, B);
            info = struct('eyes', {labs(ei)}, 'targets', {labs(ti)}, 'B', B, 'fitted', numel(idx), 'seconds', numel(idx) / EEG.srate);
        end

        function [EEG, info] = weightTrials(EEG, roi, groups, scores)
            % Weights for averaging the epochs, from their noise at the
            % electrodes roi (labels); groups: the event types of each
            % condition. An epoch's noise is the mean square, over those
            % electrodes and all its samples, of the epoch minus the mean of
            % the other epochs of its condition (leave-one-out, so the epoch
            % does not lower its own noise). scores (optional): the measured
            % scores, a cell of {electrodes, window in s}, each the mean
            % over its electrodes and window. The part of that residual
            % that goes with the scores is left out of the noise first (in
            % each condition, the residuals' projection on the epochs'
            % score residuals): EEG noise is slow, so an epoch whose score
            % happens to be far from the others' also has a large residual
            % over the epoch, and weighting by it would count such epochs
            % less for their score, which makes the weighted mean look more
            % precise than it is (its SME a third too low in simulation).
            % What is left is uncorrelated with the scores and still
            % follows each epoch's noise level. Its weight is 1 / noise, scaled
            % to a mean of 1 per condition: with noise that differs between
            % epochs only in size, these weights give the weighted mean with
            % the smallest error (inverse-variance weighting). A noise below
            % a tenth of the condition's median counts as that tenth, so no
            % epoch counts more than 10 times the typical one. Stored in
            % EEG.etc.pipecompare.trialWeights by the urevent of each
            % epoch's time-locking event, so they follow the epochs through
            % later selections; no epoch is removed.
            codes = {};
            for c = 1:numel(groups), codes = [codes cellstr(groups{c})]; end %#ok<AGROW>
            codes = unique(codes, 'stable');
            lock = pipecompare.eval.Measure.lockingEvents(EEG, codes);
            labels = lower({EEG.chanlocs.labels});
            [ok, ch] = ismember(lower(cellstr(roi)), labels);
            assert(all(ok), 'PipeCompare:RoiMissing', 'Electrode(s) missing for the weights: %s', strjoin(cellstr(roi(~ok)), ', '));
            ch = unique(ch);
            u = nan(1, EEG.trials); g = zeros(1, EEG.trials);
            for k = find(lock > 0)
                e = EEG.event(lock(k)); t = strtrim(char(string(e.type)));
                g(k) = find(cellfun(@(ev) any(strcmp(cellstr(ev), t)), groups), 1);
                u(k) = double(e.urevent);
            end
            w = nan(1, EEG.trials); X = double(EEG.data(ch, :, :));
            % each score as a kernel over those electrodes and the samples
            if nargin < 4, scores = {}; end
            times = EEG.xmin + (0:EEG.pnts - 1) / EEG.srate;
            Kn = zeros(numel(ch) * EEG.pnts, 0);
            for j = 1:numel(scores)
                kc = ismember(lower({EEG.chanlocs(ch).labels}), lower(cellstr(scores{j}{1})))';
                kt = times >= scores{j}{2}(1) - 1e-9 & times <= scores{j}{2}(2) + 1e-9;
                if ~any(kc) || ~any(kt), continue; end
                kk = double(kc) / sum(kc) * double(kt) / sum(kt);   % channels x samples
                Kn(:, end+1) = kk(:); %#ok<AGROW>
            end
            ratio = 1;
            for c = 1:numel(groups)
                k = find(g == c); n = numel(k);
                if n < 2, w(k) = 1; continue; end
                Xc = X(:, :, k);
                R = Xc - (sum(Xc, 3) - Xc) / (n - 1);
                if ~isempty(Kn) && n > size(Kn, 2) + 1
                    Rm = reshape(R, [], n)';   % epochs x (channels x samples)
                    U = Rm * Kn;               % the score residuals
                    Rm = Rm - U * (pinv(U) * Rm);
                    R = reshape(Rm', size(R));
                end
                v = reshape(mean(mean(R .^ 2, 1), 2), 1, n);
                v = max(v, median(v) / 10);
                wc = 1 ./ v; wc = wc / mean(wc);
                w(k) = wc;
                ratio = max(ratio, max(wc) / min(wc));
            end
            keep = g > 0;
            tw = struct('urevent', u(keep), 'weight', w(keep), 'events', {codes});
            EEG.etc.pipecompare.trialWeights = tw;
            info = struct('trialWeights', tw, 'weightRatio', ratio);
            pipecompare.utils.log('Epochs weighted by their noise: the cleanest counts %.1f times as much as the noisiest.', ratio);
        end

        function chans = testedChannels(EEG, p)
            % the channels an epoch rejection tests (all but p.exclude)
            chans = 1:EEG.nbchan;
            if isfield(p, 'exclude') && ~isempty(p.exclude)
                chans = find(~ismember(lower({EEG.chanlocs.labels}), lower(cellstr(p.exclude))));
            end
        end

        function m = limitStat(EEG, chans, p)
            % Per epoch, the value reject_threshold compares with its limit
            % (rejected when m > uv): the largest absolute value over chans
            % and the whole epoch, or with method peaktopeak the largest
            % moving-window peak-to-peak value.
            if strcmp(pipecompare.utils.fieldOr(p, 'method', 'absolute'), 'peaktopeak')
                m = max(pipecompare.run.Steps.peakToPeak(EEG, chans, pipecompare.utils.fieldOr(p, 'window', 200)), [], 1);
            else
                m = reshape(max(max(abs(EEG.data(chans, :, :)), [], 1), [], 2), 1, []);
            end
        end

        function EEG = markPeakToPeak(EEG, chans, uv, win)
            % Moving-window peak-to-peak test (ERPLAB's pop_artmwppth, as in
            % ERP CORE): within each window of win ms, moved in 50 ms steps
            % over the epoch, the highest minus the lowest value of a
            % channel; an epoch is marked when this exceeds uv on any of the
            % channels chans. Marks go where pop_eegthresh puts its own
            % (EEG.reject.rejthresh, per channel rejthreshE), so they are
            % removed and repaired the same way.
            pp = pipecompare.run.Steps.peakToPeak(EEG, chans, win);
            E = false(EEG.nbchan, EEG.trials);
            E(chans, :) = pp > uv;
            EEG.reject.rejthreshE = double(E);
            EEG.reject.rejthresh = double(any(E, 1));
            pipecompare.utils.log('%d/%d epochs over %g uV peak-to-peak (%g ms windows).', sum(any(E, 1)), EEG.trials, uv, win);
        end

        function pp = peakToPeak(EEG, chans, win)
            % channels x epochs: the largest peak-to-peak value within win
            % ms windows moved in 50 ms steps (markPeakToPeak's test value)
            n = max(2, round(win / 1000 * EEG.srate));
            step = max(1, round(0.050 * EEG.srate));
            starts = 1:step:max(1, EEG.pnts - n + 1);
            pp = zeros(numel(chans), EEG.trials);
            for s0 = starts
                X = double(EEG.data(chans, s0:min(EEG.pnts, s0 + n - 1), :));
                pp = max(pp, reshape(max(X, [], 2) - min(X, [], 2), numel(chans), EEG.trials));
            end
        end

        function EEG = interpolateEpochs(EEG, chans, epochs)
            % Within the epochs epochs{k}, the channels chans{k} replaced by
            % EEGLAB's spherical spline interpolation (eeg_interp, Perrin et
            % al., 1989) from the other channels of the same epochs.
            % eeg_interp is linear in the data, so its weights are read once
            % per set of channels (its output for unit impulses) and applied
            % to those epochs: the same numbers eeg_interp gives on them.
            for k = 1:numel(chans)
                bad = chans{k}; ep = epochs{k};
                if isempty(bad) || isempty(ep), continue; end
                W = interpWeights(EEG, bad);
                X = reshape(double(EEG.data(:, :, ep)), EEG.nbchan, []);
                EEG.data(bad, :, ep) = cast(reshape(W * X, numel(bad), EEG.pnts, numel(ep)), 'like', EEG.data);
            end
        end

        function [EEG, coms, info] = native(EEG, command)
            % Replay a command captured from an EEGLAB dialog, statement by
            % statement, in an isolated workspace that only contains EEG.
            % When a native step only marks and removes (a captured
            % workflow), the removals it took on these data are recorded
            % (info.decisions) so the signal check can apply the same
            % removals instead of re-deciding on its own copy.
            st = pipecompare.run.Native.statements(command);
            coms = cell(1, numel(st));
            info = struct('decisions', struct('kind', {}, 'idx', {}), 'rejected', 0, 'icsRemoved', 0);
            decided = true;
            for k = 1:numel(st)
                cmd = st{k};
                assert(isempty(regexp(cmd, '\<(ALLEEG|CURRENTSET|STUDY|CURRENTSTUDY)\>', 'once')), ...
                    'PipeCompare:Native', 'The command refers to other datasets (ALLEEG/STUDY); it cannot be replayed on candidates.');
                % Display-only options must not open windows during a search.
                cmd = regexprep(cmd, '''plotfreqz''\s*,\s*1', '''plotfreqz'',0');
                cmd = regexprep(cmd, '''interrupt''\s*,\s*''on''', '''interrupt'',''off''');
                e = pipecompare.live.History.classify(cmd);
                if strcmp(e.step, 'reject_epochs') && EEG.trials > 1
                    n0 = EEG.trials; gone = []; known = false;
                    if strcmp(e.fn, 'pop_rejepoch') && numel(e.args) >= 2
                        % the epochs it removes, read before it runs
                        v = evalExprWithEEG(EEG, e.args{2});
                        if islogical(v) || (numel(v) == n0 && all(ismember(v(:), [0 1]))), gone = find(v(:)');
                        else, gone = unique(double(v(:)')); end
                        known = true;
                        EEG = evalWithEEG(EEG, cmd);
                    elseif all(ismember(1:n0, unique([EEG.event.epoch])))
                        % other removals: every epoch carries an event, so the
                        % surviving events tell which epochs are left
                        for q = 1:numel(EEG.event), EEG.event(q).nqc_epoch = EEG.event(q).epoch; end
                        EEG = evalWithEEG(EEG, cmd);
                        if isfield(EEG.event, 'nqc_epoch')
                            gone = setdiff(1:n0, unique(arrayfun(@(x) double(x.nqc_epoch), EEG.event)));
                            EEG.event = rmfield(EEG.event, 'nqc_epoch'); known = true;
                        end
                    else
                        EEG = evalWithEEG(EEG, cmd);   % an epoch without events: cannot be told apart
                    end
                    assert(EEG.trials > 0, 'PipeCompare:AllRejected', 'every epoch would be rejected');
                    if known && EEG.trials == n0 - numel(gone)
                        info.decisions(end+1) = struct('kind', 'epochs', 'idx', gone);
                    else
                        decided = false;               % re-run on the signal copy, flagged
                    end
                    info.rejected = info.rejected + (n0 - EEG.trials);
                elseif strcmp(e.fn, 'pop_subcomp') && ~isempty(EEG.icaweights)
                    W0 = EEG.icaweights;
                    EEG = evalWithEEG(EEG, cmd);
                    gone = find(~ismember(W0, EEG.icaweights, 'rows'))';
                    info.decisions(end+1) = struct('kind', 'comps', 'idx', gone);
                    info.icsRemoved = info.icsRemoved + numel(gone);
                else
                    if ~strcmp(e.kind, 'mark'), decided = false; end
                    EEG = evalWithEEG(EEG, cmd);
                end
                coms{k} = cmd;
            end
            if ~decided, info.decisions = info.decisions([]); end   % not only marks and removals: re-run instead
        end
    end
end

function f = rejField(type)
% the EEG.reject field a rejection step marks
f = struct('reject_threshold', 'rejthresh', 'reject_jointprob', 'rejjp', 'reject_kurtosis', 'rejkurt');
f = f.(type);
end

function c = stretchCommand(st)
% the command that interpolates the channels st in their stretches
samples = arrayfun(@(q) mat2str(q.samples), st, 'UniformOutput', false);
c = sprintf(['EEG = pipecompare.run.Steps.interpStretches(EEG, struct(''channel'', %s, ''samples'', {%s})); ', ...
    '%% PipeCompare: channels bad for part of the recording interpolated there only'], strCells({st.channel}), strjoin(samples, ' '));
end

function t = strCells(c)
% a cell of strings as MATLAB code
c = cellstr(c);
t = ['{' strjoin(cellfun(@(x) ['''' strrep(x, '''', '''''') ''''], c, 'UniformOutput', false), ' ') '}'];
end

function t = cellText(c)
% a cell of index vectors as MATLAB code
t = ['{' strjoin(cellfun(@mat2str, c, 'UniformOutput', false), ', ') '}'];
end

function W = interpWeights(EEG, bad)
% The rows of eeg_interp's spherical interpolation for channels bad: its
% output on a one-trial dataset holding a unit impulse on each channel.
n = EEG.nbchan;
T = eeg_emptyset();
T.chanlocs = EEG.chanlocs; T.chaninfo = EEG.chaninfo; T.nbchan = n;
T.data = eye(n); T.pnts = n; T.trials = 1; T.srate = 1; T.xmin = 0; T.xmax = n - 1;
[~, T] = evalc('eeg_interp(T, bad, ''spherical'')');
W = double(T.data(bad, :));
end

function EEG = evalWithEEG(EEG, PIPECOMPARE_CMD__)
eval(PIPECOMPARE_CMD__);
end

function v = evalExprWithEEG(EEG, PIPECOMPARE_EXPR__) %#ok<INUSL>
% value of an argument expression of a captured command (e.g.
% EEG.reject.rejthresh), in a workspace that only contains EEG
v = eval(PIPECOMPARE_EXPR__);
end

function tf = isFixedTransform(stmt)
% EEGLAB calls whose effect does not depend on the data values (filters,
% resampling, re-referencing, baseline, epoching, channel selection or
% interpolation by explicit lists): re-running them on another dataset
% applies exactly the same operation.
e = pipecompare.live.History.classify(stmt);
% Filters are linear time-invariant operators fixed by their design
% parameters (Widmann, Schroger & Maess, 2015); selections by events or
% channel lists depend on the events, which the signal copy shares.
tf = any(strcmp(e.fn, {'pop_eegfiltnew','pop_firws','pop_firma','pop_firpm','pop_eegfilt','pop_basicfilter', ...
    'pop_resample','pop_reref','pop_rmbase','pop_epoch','pop_select','pop_selectevent','pop_rmdat', ...
    'pop_interp','pop_chanedit'}));
end

function L = presentChannels(EEG, wanted, what)
% The listed channels that are still in the data. A channel of the
% starting montage that an earlier step removed is skipped (nothing left
% to exclude); a label that never existed is an error (a typo).
wanted = cellstr(wanted);
here = ismember(lower(wanted), lower({EEG.chanlocs.labels}));
known = {EEG.chanlocs.labels};
if isfield(EEG, 'etc') && isstruct(EEG.etc) && isfield(EEG.etc, 'pipecompare') && isfield(EEG.etc.pipecompare, 'rootChanlocs')
    known = [known {EEG.etc.pipecompare.rootChanlocs.labels}];
end
unknown = wanted(~ismember(lower(wanted), lower(known)));
assert(isempty(unknown), 'PipeCompare:Channels', '%s: channel(s) not in the dataset: %s', what, strjoin(unknown, ', '));
if any(~here)
    pipecompare.utils.log('%s: %s already removed by an earlier step; nothing to exclude there.', what, strjoin(wanted(~here), ', '));
end
L = wanted(here);
end

function rec = recordAsr(E0, E1, cutoff)
% ASR's decisions on these data, kept only if replaying them reproduces
% EEGLAB's own output (otherwise the signal check re-runs ASR, flagged).
rec = [];
try
    r = pipecompare.run.AsrRecord.record(E0, cutoff);
    Y = pipecompare.run.AsrRecord.apply(r, E0.data);
    err = max(abs(Y(:) - double(E1.data(:)))) / max(1, max(abs(double(E1.data(:)))));
    if isequal(size(Y), size(E1.data)) && err < 1e-6
        rec = r;
    else
        pipecompare.utils.log('ASR decisions not reproduced exactly (relative difference %.2g); the signal check re-runs ASR.', err);
    end
catch ME
    pipecompare.utils.log('ASR decisions not recorded (%s); the signal check re-runs ASR.', ME.message);
end
end

function [a, fb] = fftArgs(EEG, edge)
% A high-pass with a long FIR (above order 2000, as EEGLAB's dialog
% advises; e.g. 0.1 or 0.3 Hz at 250 Hz, 33,000 taps for 0.1 Hz at 1000
% Hz) runs in the frequency domain (firfilt's fftfilt option): the same
% filter, much faster. fftfilt is the Signal Processing Toolbox's; without
% it, PipeCompare's own (fallback/fftfilt.m, overlap-add, the same result
% to rounding) is put on the path while the filter runs: fb, an object
% that takes it off again when it is cleared ([] when not needed).
a = {}; fb = [];
df = min(max(edge * 0.25, 2), edge);              % pop_eegfiltnew's default transition band (high-pass)
if ~(3.3 * EEG.srate / df > 2000), return; end
if exist('fftfilt', 'file') ~= 2
    d = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'fallback');
    if ~isfolder(d), return; end
    addpath(d); fb = onCleanup(@() rmpath(d));
end
a = {'usefftfilt', 1};
end

function EEG = lowRankRun(inst, EEG, ctx)
% A filter step on the signal-check copy. The copy holds only the known
% signals: one field over the channels per measure, each with its own time
% course, so its channels x samples are U G with U orthonormal (channels x
% k, k the number of fields, often 1 to 3) and G = U' X. The filter is
% linear, time-invariant and the same for every channel, so filtering the
% channels is filtering the k rows of G: k channels filtered instead of
% all of them, the same result to rounding. The factors come from the
% eigenvectors of X X'; the step runs as usual when the copy is not of low
% rank (k above a quarter of the channels) or U G does not give it back.
n = EEG.nbchan;
Y = reshape(EEG.data, n, []); m = size(Y, 2);
blocks = arrayfun(@(b) b:min(m, b + 1e5 - 1), 1:1e5:m, 'UniformOutput', false);   % (memory)
C = zeros(n);
for q = 1:numel(blocks), X = double(Y(:, blocks{q})); C = C + X * X'; end
[V, d] = eig((C + C') / 2, 'vector');
[d, o] = sort(d, 'descend'); V = V(:, o);
k = sum(d > 1e-12 * max(d));
if k == 0 || k > n / 4
    EEG = pipecompare.run.Steps.run(inst, EEG, ctx); return;
end
U = V(:, 1:k); G = zeros(k, m); left = 0;
for q = 1:numel(blocks)
    X = double(Y(:, blocks{q})); G(:, blocks{q}) = U' * X;
    left = left + norm(X - U * G(:, blocks{q}), 'fro') ^ 2;
end
if sqrt(left) > 1e-6 * sqrt(sum(d))
    EEG = pipecompare.run.Steps.run(inst, EEG, ctx); return;
end
F = EEG;
F.data = cast(reshape(G, k, EEG.pnts, EEG.trials), 'like', EEG.data);
F.nbchan = k;
if numel(EEG.chanlocs) >= k, F.chanlocs = EEG.chanlocs(1:k); end
F.icaweights = []; F.icasphere = []; F.icawinv = []; F.icaact = []; F.icachansind = [];
F = pipecompare.run.Steps.run(inst, F, ctx);
EEG.data = cast(reshape(U * double(reshape(F.data, k, [])), n, F.pnts, F.trials), 'like', EEG.data);
end

function A = erpAt(avg, grp, pos, m)
% the averages of eogRegress at these samples (0 outside the epochs)
A = zeros(m, numel(grp));
for g = 1:numel(avg)
    in = grp == g;
    if any(in), A(:, in) = avg{g}(:, pos(in)); end
end
end

function EEG = applyEog(EEG, ti, ei, B)
% EEG.data(ti) minus B times EEG.data(ei), every sample and epoch
n = EEG.pnts * EEG.trials;
D = reshape(EEG.data, EEG.nbchan, n);
D(ti, :) = cast(double(D(ti, :)) - B * double(D(ei, :)), 'like', D);
EEG.data = reshape(D, EEG.nbchan, EEG.pnts, EEG.trials);
end

function com = plainCommand(com, fb)
% the command as a script would run it where fftfilt is missing: without
% the frequency-domain option when PipeCompare's stand-in was used (the
% same filter, only slower)
if ~isempty(fb), com = regexprep(com, ',\s*''usefftfilt''\s*,\s*1', ''); end
end

function requireLocations(locs, what)
% The channels to interpolate must have positions: EEGLAB's eeg_interp
% leaves a channel without one as it was, silently, while the step would
% report it interpolated.
ok = @(v) isnumeric(v) && isscalar(v) && isfinite(v);
has = arrayfun(@(c) isfield(c, 'X') && ok(c.X) && ok(c.Y) && ok(c.Z), locs);
assert(all(has), 'PipeCompare:Chanlocs', ['%s: no channel location for %s, so it cannot be interpolated ', ...
    '(Edit > Channel locations, or exclude/remove it).'], what, strjoin({locs(~has).labels}, ', '));
end

function target = rootChanlocs(EEG)
assert(isfield(EEG, 'etc') && isfield(EEG.etc, 'pipecompare') && isfield(EEG.etc.pipecompare, 'rootChanlocs'), ...
    'PipeCompare:Restore', 'Root channel montage not recorded');
target = EEG.etc.pipecompare.rootChanlocs;
end

function warnIclabel15()
% ICLabel 1.5 (bundled with EEGLAB 2024.0) draws the scalp maps it
% classifies rotated, so eye components are mostly not recognised; 1.6
% fixed it. Said once per MATLAB session.
persistent said
if ~isempty(said), return; end
said = true;
f = which('topoplotFast');
if ~isempty(f) && contains(fileread(f), '[xi,yi] = meshgrid(xi,yi);')
    pipecompare.utils.log(['WARNING: this ICLabel (1.5, as bundled with EEGLAB 2024.0) misclassifies components ', ...
        '(rotated scalp maps; eye components are mostly missed). Update ICLabel to 1.6 or later in ', ...
        'File > Manage EEGLAB extensions.']);
end
end

function hit = beyondRobust(X, z)
% the samples of each row of X beyond z robust SDs (median, MAD) of that row
X = abs(X - median(X, 2));
hit = X > z * max(1.4826 * median(X, 2), eps);
end

function C = blockCov(Y, keep)
% the covariance of the columns of Y (channels x samples) where keep is
% true, a block of samples at a time (no copy of the whole data)
n = size(Y, 1); C = zeros(n); s = zeros(n, 1); m = 0;
for b = 1:1e5:size(Y, 2)
    k = b:min(size(Y, 2), b + 1e5 - 1);
    X = double(Y(:, k(keep(k))));
    C = C + X * X'; s = s + sum(X, 2); m = m + size(X, 2);
end
C = (C - s * s' / max(m, 1)) / max(m - 1, 1);
end

function T = bandCopy(T, doHp, hp, doLp, lp)
% the copy bad channels are tested on: high-passed at hp and/or low-passed at lp
if doHp
    [fa, fb] = fftArgs(T, hp); %#ok<ASGLU>
    T = pop_eegfiltnew(T, 'locutoff', hp, 'plotfreqz', 0, fa{:});
end
if doLp, T = pop_eegfiltnew(T, 'hicutoff', lp, 'plotfreqz', 0); end
end

function EEG = runIca(EEG, opts)
% pop_runica(EEG, opts{:}) for opts = {'icatype', 'runica', runica options}.
% EEGLAB before 2025.1 adds 'interrupt', 'on' to the options even from the
% command line, so runica opens a window with an Interrupt button for
% every fit, and a click stops the fit half way. There the same
% computation is made here, as pop_runica makes it, without that window.
if ~contains(fileread(which('pop_runica')), 'if ismatlab, g.options')
    EEG = pop_runica(EEG, opts{:});
    return
end
o = opts(3:end);
x = reshape(EEG.data, EEG.nbchan, EEG.pnts * EEG.trials);
r = icaRank(double(x(:, 1:min(3000, size(x, 2)))), contains(fileread(which('pop_runica')), 'tmprank2 = min('));
x = x - repmat(mean(x, 2), [1 size(x, 2)]);
x = double(x);
x = x - repmat(mean(x, 2), [1 size(x, 2)]);
if r < size(x, 1), o = [o {'pca', r}]; end
[EEG.icaweights, EEG.icasphere] = runica(x, 'lrate', 0.001, o{:});
clear x
EEG.icachansind = 1:EEG.nbchan; EEG.icaact = [];
EEG.icawinv = pinv(EEG.icaweights * EEG.icasphere);
v = sum(EEG.icawinv .^ 2) .* sum(transpose((EEG.icaweights * EEG.icasphere) * EEG.data(EEG.icachansind, :)) .^ 2) ...
    / ((numel(EEG.icachansind) * EEG.pnts) - 1);
[~, k] = sort(v); k = k(end:-1:1);   % by variance, as pop_runica
EEG.icaweights = EEG.icaweights(k, :);
EEG.icawinv = pinv(EEG.icaweights * EEG.icasphere);
EEG = eeg_checkset(EEG);
end

function r = icaRank(x, takeMin)
% pop_runica's getrank (EEGLAB 2024.0 to 2024.2 take the smaller of the
% two estimates, older releases the larger)
r1 = rank(x);
[~, D] = eig(cov(x', 1));
r = sum(diag(D) > 1e-7);
if r1 ~= r
    if takeMin, r = min(r1, r); else, r = max(r1, r); end
end
end

function s = keptShift(EEG, idx, contract)
% For each mean-amplitude component (not contralateral minus ipsilateral)
% and condition: the mean score of the epochs kept against that of all of
% them, as z = (kept - all) / (SD * sqrt(1/nKept - 1/n)), the spread that
% difference has when the rejected epochs are a random draw, signed so
% that z < 0 means the kept epochs have less of the component (a smaller
% positive or a smaller negative deflection). s: the component and
% condition with the lowest z; [] when none is measured.
s = [];
T = pipecompare.eval.Measure.trials(EEG, contract);
rej = ismember(T.epoch, idx);
for j = 1:numel(contract.components)
    comp = contract.components(j);
    if contract.isLateral(j) || ~strcmp(comp.measure, 'mean'), continue; end
    sg = 1; if strcmp(comp.polarity, 'negative'), sg = -1; end
    x = double(T.data{j}(:, 1));
    for c = unique(T.cond(:))'
        r = T.cond == c; k = r & ~rej;
        n = sum(r); nk = sum(k);
        if nk < 2 || nk == n, continue; end
        d = mean(x(k)) - mean(x(r));
        z = sg * d / max(std(x(r)) * sqrt(1 / nk - 1 / n), eps);
        if isempty(s) || z < s.z
            s = struct('z', z, 'component', comp.name, 'condition', contract.conditions(c).name, ...
                'kept', mean(x(k)), 'all', mean(x(r)), 'nKept', nk, 'n', n);
        end
    end
end
end
