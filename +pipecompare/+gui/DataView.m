classdef DataView < handle
    %DATAVIEW A plot of data inside the PipeCompare window.
    %
    %   v = pipecompare.gui.DataView(parent)
    %   v.showData(EEG)    the dataset: continuous data as one trace per
    %                      channel, 10 s at a time (the slider moves through
    %                      the recording); epoched data as the average of
    %                      the epochs of each event type they are
    %                      time-locked to, at the channel chosen
    %   v.addPipeline(E, contract, name)
    %                      the data of a pipeline (E, rebuilt): the average
    %                      of each condition at the measured electrodes,
    %                      with the measurement window marked; band power:
    %                      the spectrum at them, with the band marked. Added
    %                      to what is shown, to compare pipelines.
    %   v.clearPlot()

    properties
        Root; Axes; Slider; ChanDrop; Note
        EEG = []            % what showData shows
        Names = {}          % the pipelines addPipeline shows
        Seconds = 10        % continuous data: the length shown
    end

    methods
        function obj = DataView(parent)
            g = uigridlayout(parent, [3 1]); g.RowHeight = {34, '1x', 'fit'}; g.Padding = [0 0 0 0];
            obj.Root = g;
            c = uigridlayout(g, [1 4]); c.Padding = [0 0 0 0]; c.ColumnWidth = {60, '1x', 60, 120};
            uilabel(c, 'Text', 'Time (s)');
            obj.Slider = uislider(c, 'Limits', [0 1], 'Enable', 'off', 'ValueChangedFcn', @(~, ~) obj.drawData());
            uilabel(c, 'Text', 'Channel', 'HorizontalAlignment', 'right');
            obj.ChanDrop = uidropdown(c, 'Items', {'-'}, 'Enable', 'off', 'ValueChangedFcn', @(~, ~) obj.drawData());
            obj.Axes = uiaxes(g);
            obj.Note = uilabel(g, 'Text', '', 'WordWrap', 'on', 'FontColor', [0.3 0.3 0.3]);
        end

        function showData(obj, EEG)
            obj.EEG = EEG; obj.Names = {};
            obj.Slider.Enable = 'off'; obj.ChanDrop.Enable = 'off';
            if isempty(EEG), obj.reset(); obj.Note.Text = ''; return; end
            if EEG.trials == 1
                obj.Slider.Limits = [0 max(EEG.pnts / EEG.srate - obj.Seconds, 1e-3)];
                obj.Slider.Value = 0; obj.Slider.Enable = 'on';
            else
                obj.ChanDrop.Items = labelsOf(EEG); obj.ChanDrop.ItemsData = 1:EEG.nbchan;
                obj.ChanDrop.Value = 1; obj.ChanDrop.Enable = 'on';
            end
            obj.drawData();
        end

        function drawData(obj)
            E = obj.EEG; ax = obj.Axes;
            if isempty(E), return; end
            obj.reset();
            labels = labelsOf(E);
            if E.trials == 1
                t0 = obj.Slider.Value;
                i = max(1, round(t0 * E.srate) + 1):min(E.pnts, round((t0 + obj.Seconds) * E.srate));
                i = i(1:max(1, floor(numel(i) / 2000)):end);   % at most about 2000 points per channel
                X = double(E.data(:, i)); X = X - mean(X, 2);
                gap = 4 * median(std(X, 0, 2)); if ~(gap > 0), gap = 1; end
                n = size(X, 1);
                plot(ax, (i - 1) / E.srate, (X + (n:-1:1)' * gap)', 'Color', [0.1 0.25 0.55], 'LineWidth', 0.5);
                ax.YTick = (1:n) * gap; ax.YTickLabel = labels(n:-1:1);   % the first channel at the top
                ax.YLim = [0 (n + 1) * gap]; ax.XLim = t0 + [0 obj.Seconds];
                xlabel(ax, 'Time (s)'); ylabel(ax, ''); title(ax, '');
                obj.Note.Text = sprintf('%d channels, %g s from %.1f s; traces %.3g apart (data units).', n, obj.Seconds, t0, gap);
            else
                ch = obj.ChanDrop.Value;
                ty = lockTypes(E); types = unique(ty, 'stable');
                hold(ax, 'on');
                for k = 1:numel(types)
                    sel = strcmp(ty, types{k});
                    plot(ax, E.times, mean(double(E.data(ch, :, sel)), 3), 'DisplayName', sprintf('%s (%d)', types{k}, sum(sel)));
                end
                hold(ax, 'off');
                legend(ax, 'show', 'Interpreter', 'none', 'Location', 'best');
                xlabel(ax, 'Time (ms)'); ylabel(ax, 'Amplitude'); title(ax, labels{ch}, 'Interpreter', 'none');
                obj.Note.Text = sprintf('Average of the %d epochs of each event type they are time-locked to.', E.trials);
            end
        end

        function addPipeline(obj, E, c, name)
            obj.EEG = []; obj.Slider.Enable = 'off'; obj.ChanDrop.Enable = 'off';
            ax = obj.Axes;
            if isempty(obj.Names), obj.reset(); end
            hold(ax, 'on');
            labels = lower(labelsOf(E));
            note = 'The lines of several pipelines stay, to compare them; Clear plot starts again.';
            if strcmp(c.analysis, 'bandpower')
                b = c.bands(1); roi = find(ismember(labels, lower(b.roi)));
                X = double(E.data(roi, :, :)); n = size(X, 2);
                P = mean(mean(abs(fft(X - mean(X, 2), [], 2)) .^ 2, 3), 1) / (n * E.srate);   % average spectrum
                f = (0:n - 1) * E.srate / n; keep = f >= 0.5 & f <= min(60, E.srate / 2);
                plot(ax, f(keep), 10 * log10(P(keep)), 'DisplayName', name);
                if isempty(obj.Names), markRange(ax, b.freq); end
                xlabel(ax, 'Frequency (Hz)'); ylabel(ax, 'Power (dB)');
                title(ax, sprintf('%s (%g-%g Hz) at %s', b.name, b.freq, strjoin(b.roi, ', ')), 'Interpreter', 'none');
            else
                m = c.components(1); roi = find(ismember(labels, lower(m.roi)));
                if E.trials == 1, [~, E] = evalc('pop_epoch(E, c.allEvents(), c.epoch)'); end
                ew = pipecompare.eval.Measure.epochWeights(E);   % 1 unless the pipeline weights its epochs
                if any(ew ~= 1), note = ['Weighted averages (the pipeline weights the epochs by their noise). ' note]; end
                for k = 1:numel(c.conditions)
                    sel = pipecompare.eval.Measure.lockingEvents(E, c.conditions(k).events) > 0;
                    w = sum(mean(double(E.data(roi, :, sel)), 1) .* reshape(ew(sel), 1, 1, []), 3) / sum(ew(sel));
                    plot(ax, E.times, w, 'DisplayName', sprintf('%s, %s (%d trials)', name, c.conditions(k).name, sum(sel)));
                end
                if isempty(obj.Names), markRange(ax, 1000 * m.window); end
                xlabel(ax, 'Time (ms)'); ylabel(ax, 'Amplitude (uV)');
                title(ax, sprintf('%s at %s', m.name, strjoin(m.roi, ', ')), 'Interpreter', 'none');
                if c.isLateral(1)
                    note = ['Average of the electrodes; the comparison scores contralateral minus ipsilateral. ' note];
                end
            end
            hold(ax, 'off');
            legend(ax, 'show', 'Interpreter', 'none', 'Location', 'best');
            obj.Names{end+1} = name;
            obj.Note.Text = note;
        end

        function clearPlot(obj)
            obj.Names = {}; obj.reset(); obj.Note.Text = '';
        end

        function reset(obj)
            ax = obj.Axes;
            cla(ax); legend(ax, 'off'); title(ax, '');
            ax.YTickMode = 'auto'; ax.YTickLabelMode = 'auto'; ax.XLimMode = 'auto'; ax.YLimMode = 'auto';
        end
    end
end

function L = labelsOf(E)
L = arrayfun(@(k) sprintf('%d', k), 1:E.nbchan, 'UniformOutput', false);
if isfield(E, 'chanlocs') && numel(E.chanlocs) == E.nbchan && isfield(E.chanlocs, 'labels')
    L = arrayfun(@(c) char(string(c.labels)), E.chanlocs, 'UniformOutput', false);
end
L = L(:)';
end

function ty = lockTypes(E)
% the event type at time 0 of each epoch ('' when there is none)
ty = repmat({''}, 1, E.trials);
if ~isfield(E, 'epoch') || numel(E.epoch) ~= E.trials, return; end
for k = 1:E.trials
    lat = E.epoch(k).eventlatency; t = E.epoch(k).eventtype;
    if ~iscell(lat), lat = {lat}; end
    if ischar(t) || isstring(t), t = {char(t)}; elseif ~iscell(t), t = num2cell(t); end
    if isempty(lat), continue; end
    [~, q] = min(abs(cellfun(@(x) double(x(1)), lat)));
    if q <= numel(t), ty{k} = strtrim(char(string(t{q}))); end
end
end

function markRange(ax, r)
% the measurement window (or band), as two dashed lines
for x = r(:)'
    xline(ax, x, '--', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
end
end
