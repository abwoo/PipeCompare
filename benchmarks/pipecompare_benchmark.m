function out = pipecompare_benchmark(nRec, outDir)
%PIPECOMPARE_BENCHMARK Ground-truth benchmark: does the pipeline PipeCompare
%   recommends measure the P3 more accurately than common fixed pipelines?
%
%   out = pipecompare_benchmark(nRec, outDir)
%
%   nRec simulated recordings (bench_synth, seeds 1..nRec; default 8, or
%   the environment variable PIPECOMPARE_BENCH_N), each with a known P3,
%   time-locked blinks, muscle, movement artifacts and a bad channel. On
%   each, the P3 (ERP CORE: Pz, 300-600 ms, mean amplitude) of targets and
%   standards is measured by
%     PipeCompare  the pipeline pop_pipecompare recommends (simple mode,
%                  Standard steps, average reference)
%     minimal      0.1-30 Hz, average reference, epochs, baseline,
%                  rejection at 100 uV
%     common       0.1-30 Hz, bad channels (kurtosis and joint
%                  probability, 5 SD) interpolated, average reference,
%                  ICA with ICLabel (Eye and Muscle above 0.8) removed,
%                  epochs, baseline, rejection at 100 uV
%   and compared with the P3 actually in the data (each trial's P3 at Pz
%   against the average reference, averaged over the window and the
%   trials of its condition). Per
%   pipeline: the bias (mean error) and root mean square error over the
%   recordings, for each condition and for their difference, and the
%   trials kept. Also, among all the pipelines PipeCompare compared, where
%   its recommendation ranks by the actual error of the difference.
%   outDir (default tempdir): benchmark.csv (one row per recording and
%   pipeline) and benchmark_summary.txt.
if nargin < 1 || isempty(nRec)
    nRec = str2double(getenv('PIPECOMPARE_BENCH_N'));
    if isnan(nRec), nRec = 8; end
end
if nargin < 2 || isempty(outDir), outDir = tempdir; end
names = {'PipeCompare', 'minimal', 'common'};
rows = {};
rankOfPick = nan(1, nRec);
for s = 1:nRec
    [EEG, truth] = bench_synth(s);
    setBase(EEG);
    [~, ~, r] = pop_pipecompare(EEG, 'measure', 'P3', 'events', {'11', '31'}, 'reference', 'average', 'show', 'off');
    c = r.contract;
    want = trueP3(truth, c, EEG.srate);
    eog = {'EOG1', 'EOG2'};
    est = nan(3, 2); kept = nan(3, 1); label = cell(3, 1);
    % PipeCompare's recommendation, and the actual error of each pipeline it compared
    k = r.ranking.recommended;
    if ~isempty(k)
        [est(1, :), kept(1)] = estimateOf(r, k(1)); label{1} = r.labels{k(1)};
        feas = find(strcmp(r.ranking.table.status, 'feasible'))';
        err = arrayfun(@(j) abs(diff(estimateOf(r, j)) - diff(want)), feas);
        rankOfPick(s) = sum(err < abs(diff(est(1, :)) - diff(want))) + 1;
        fprintf('[benchmark] recording %d: recommended %s; its difference error ranks %d of %d compared pipelines\n', ...
            s, label{1}, rankOfPick(s), numel(feas));
    end
    % the fixed pipelines
    p = pipecompare.plan.Plan();
    p = p.add('highpass', 'cutoff', 0.1); p = p.add('lowpass', 'cutoff', 30);
    p = p.add('reref', 'mode', 'average', 'exclude', {eog});
    p = p.add('epoch'); p = p.add('baseline');
    p = p.add('reject_threshold', 'uv', 100, 'exclude', {eog});
    [est(2, :), kept(2), label{2}] = fixed(p, c);
    q = pipecompare.plan.Plan();
    q = q.add('highpass', 'cutoff', 0.1); q = q.add('lowpass', 'cutoff', 30);
    q = q.add('badchannels', 'measure', 'kurt+prob', 'threshold', 5, 'exclude', {eog}, 'action', 'interpolate');
    q = q.add('reref', 'mode', 'average', 'exclude', {eog});
    q = q.add('ica'); q = q.add('icremove', 'threshold', 0.8, 'classes', {{'Eye', 'Muscle'}});
    q = q.add('epoch'); q = q.add('baseline');
    q = q.add('reject_threshold', 'uv', 100, 'exclude', {eog});
    [est(3, :), kept(3), label{3}] = fixed(q, c);
    for i = 1:3
        e = est(i, :) - want;
        rows(end+1, :) = {s, names{i}, label{i}, want(1), want(2), est(i, 1), est(i, 2), e(1), e(2), e(1) - e(2), kept(i), ...
            strjoin(truth.badLabels, ' '), truth.lockedBlinkShare}; %#ok<AGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'recording', 'pipeline', 'steps', 'trueTarget', 'trueStandard', ...
    'target', 'standard', 'errTarget', 'errStandard', 'errDifference', 'kept', 'badChannel', 'lockedBlinkShare'});
writetable(T, fullfile(outDir, 'benchmark.csv'));
lines = {sprintf('PipeCompare benchmark: %d simulated recordings (P3 at Pz, 300-600 ms; uV)', nRec), ...
    sprintf('%-12s %22s %22s %22s %8s', 'pipeline', 'target bias / RMSE', 'standard bias / RMSE', 'difference bias / RMSE', 'kept')};
for i = 1:3
    R = T(strcmp(T.pipeline, names{i}), :);
    f = @(x) sprintf('%+8.2f / %6.2f', mean(x, 'omitnan'), sqrt(mean(x .^ 2, 'omitnan')));
    lines{end+1} = sprintf('%-12s %22s %22s %22s %7.0f%%', names{i}, f(R.errTarget), f(R.errStandard), ...
        f(R.errDifference), 100 * mean(R.kept, 'omitnan')); %#ok<AGROW>
end
lines{end+1} = sprintf(['PipeCompare''s recommendation by the actual error of the difference, among the pipelines it ', ...
    'compared: rank %s (1 = the most accurate)'], mat2str(rankOfPick));
txt = strjoin(lines, newline);
fprintf('%s\n', txt);
fid = fopen(fullfile(outDir, 'benchmark_summary.txt'), 'w'); fprintf(fid, '%s\n', txt); fclose(fid);
out = struct('table', T, 'summary', txt, 'rankOfPick', rankOfPick);
end

function want = trueP3(truth, c, fs)
% the P3 actually in the data, against the average of the EEG channels
% (every pipeline here uses the average reference): each trial's P3
% averaged over the component's electrodes and window, then over the
% trials of each condition
comp = c.components(1);
roi = ismember(upper(truth.labels), upper(cellstr(comp.roi)));
eeg = ~startsWith(upper(truth.labels), 'EOG');
topo = truth.topo - mean(truth.topo(eeg));
tw = (ceil(comp.window(1) * fs - 1e-6):floor(comp.window(2) * fs + 1e-6)) / fs;
p3 = mean(topo(roi)) * truth.amp .* arrayfun(@(k) mean(exp(-0.5 * ((tw - truth.latency(k)) / truth.sigma) .^ 2)), 1:numel(truth.amp));
want = arrayfun(@(i) mean(p3(ismember(truth.types, cellstr(c.conditions(i).events)))), 1:numel(c.conditions));
end

function [est, kept, label] = fixed(p, c)
r = pipecompare.PipeCompare.optimize(p, c);
label = r.labels{1};
[est, kept] = estimateOf(r, 1);
end

function [est, kept] = estimateOf(r, k)
% the measured P3 of each condition (as it comes out of the pipeline, not
% corrected by the signal check's gain) and the share of trials kept
est = nan(1, 2); kept = NaN;
if ~isfield(r.cands(k), 'm') || isempty(r.cands(k).m), return; end
est = r.cands(k).m.objectives(1).estimate;
kept = r.ranking.table.minRetention(k);
end

function setBase(EEG)
% EEG becomes the current EEGLAB dataset
assignin('base', 'BENCH_TMP', EEG);
evalin('base', ['global ALLCOM; ALLEEG = []; [ALLEEG, EEG, CURRENTSET] = eeg_store([], BENCH_TMP, 0); ' ...
    'EEG.history = BENCH_TMP.history; ALLEEG(CURRENTSET).history = BENCH_TMP.history; clear BENCH_TMP;']);
end
