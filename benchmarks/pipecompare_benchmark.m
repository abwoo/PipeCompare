function out = pipecompare_benchmark(nRec, outDir, firstSeed, recipe, artifactVariance)
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
%   And whether PipeCompare's own error estimate tells the truth: for every
%   pipeline it compared, the actual error of each condition's P3 (the
%   estimate divided by the signal check's gain, minus the true P3)
%   against its ranked total error (SME, artifact and shape bias); their
%   ratio (root mean squares over recordings and conditions; 1 = right,
%   above 1 = the reported error is understated), for all pipelines and for those
%   with and without each step learned from the data (EOG regression,
%   ICA components chosen by SME, ICA, trial weights, stretches).
%   outDir (default tempdir): benchmark.csv (one row per recording and
%   pipeline), calibration.csv (one row per recording, compared pipeline
%   and condition) and benchmark_summary.txt.
%   Calibration actualError is against the POPULATION P3, matching SME's
%   sampling target; conditionalError is against the realized trial mean.
%   firstSeed (default 1) permits disjoint development/validation seed sets.
if nargin < 1 || isempty(nRec)
    nRec = str2double(getenv('PIPECOMPARE_BENCH_N'));
    if isnan(nRec), nRec = 8; end
end
if nargin < 2 || isempty(outDir), outDir = tempdir; end
if nargin < 3, firstSeed = 1; end
if nargin < 4, recipe='standard'; end
if nargin < 5, artifactVariance='pooled'; end
validateattributes(nRec,{'numeric'},{'scalar','integer','positive'});
validateattributes(firstSeed,{'numeric'},{'scalar','integer','positive'});
names = {'PipeCompare', 'minimal', 'common'};
rows = {}; calib = table();
rankOfPick = nan(1, nRec);
for rec = 1:nRec
    s = firstSeed + rec - 1;
    [EEG, truth] = bench_synth(s);
    setBase(EEG);
    [~, ~, r] = pop_pipecompare(EEG, 'measure', 'P3', 'events', {'11', '31'}, 'recipe', recipe, 'artifactvariance', artifactVariance, 'reference', 'average', 'show', 'off');
    c = r.contract;
    want = trueP3(truth, c, EEG.srate);
    population = populationP3(truth, c, EEG.srate);
    eog = {'EOG1', 'EOG2'};
    est = nan(3, 2); kept = nan(3, 1); label = cell(3, 1);
    % PipeCompare's recommendation, and the actual error of each pipeline it compared
    k = r.ranking.recommended;
    if ~isempty(k)
        [est(1, :), kept(1)] = estimateOf(r, k(1)); label{1} = r.labels{k(1)};
        feas = find(strcmp(r.ranking.table.status, 'feasible'))';
        err = arrayfun(@(j) abs(diff(estimateOf(r, j)) - diff(want)), feas);
        rankOfPick(rec) = sum(err < abs(diff(est(1, :)) - diff(want))) + 1;
        fprintf('[benchmark] recording %d: recommended %s; its difference error ranks %d of %d compared pipelines\n', ...
            s, label{1}, rankOfPick(rec), numel(feas));
    end
    feas = find(strcmp(r.ranking.table.status, 'feasible'))';
    C = pipecompare.eval.Calibration.rows(r,feas,population,want,s);
    if ~isempty(C), if isempty(calib), calib = C; else, calib = [calib; C]; end, end %#ok<AGROW>
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
    [T, txt] = writeResults(rows, rankOfPick, rec, outDir, names, calib);   % after every recording: kept if the run stops early
end
out = struct('table', T, 'summary', txt, 'rankOfPick', rankOfPick);
end

function [T, txt] = writeResults(rows, rankOfPick, nDone, outDir, names, calib)
% benchmark.csv and benchmark_summary.txt for the recordings done so far
T = cell2table(rows, 'VariableNames', {'recording', 'pipeline', 'steps', 'trueTarget', 'trueStandard', ...
    'target', 'standard', 'errTarget', 'errStandard', 'errDifference', 'kept', 'badChannel', 'lockedBlinkShare'});
writetable(T, fullfile(outDir, 'benchmark.csv'));
lines = {sprintf('PipeCompare benchmark: %d simulated recordings (P3 at Pz, 300-600 ms; uV)', nDone), ...
    sprintf('%-12s %22s %22s %22s %8s', 'pipeline', 'target bias / RMSE', 'standard bias / RMSE', 'difference bias / RMSE', 'kept')};
for i = 1:3
    R = T(strcmp(T.pipeline, names{i}), :);
    f = @(x) sprintf('%+8.2f / %6.2f', mean(x, 'omitnan'), sqrt(mean(x .^ 2, 'omitnan')));
    lines{end+1} = sprintf('%-12s %22s %22s %22s %7.0f%%', names{i}, f(R.errTarget), f(R.errStandard), ...
        f(R.errDifference), 100 * mean(R.kept, 'omitnan')); %#ok<AGROW>
end
lines{end+1} = sprintf(['PipeCompare''s recommendation by the actual error of the difference, among the pipelines it ', ...
    'compared: rank %s (1 = the most accurate)'], mat2str(rankOfPick(1:nDone)));
lines{end+1} = sprintf('Recommendations available for %d of %d recordings; unavailable recordings are not successes.', ...
    sum(isfinite(rankOfPick(1:nDone))),nDone);
if ~isempty(calib)
    C = calib;
    writetable(C, fullfile(outDir, 'calibration.csv'));
    lines{end+1} = sprintf(['Does the reported error tell the truth? actual / reported error (root mean squares; ', ...
        'above 1 = reported error understated), against population P3, over %d recordings:'], nDone);
    lines = [lines calibrationLines(C)];
    lines{end+1} = ['Calibration concerns gain-corrected means (estimate/gain). Raw adopted-data means ', ...
        'are not automatically divided by this diagnostic gain; rawError and rawConditionalError are also exported.'];
end
txt = strjoin(lines, newline);
fprintf('%s\n', txt);
fid = fopen(fullfile(outDir, 'benchmark_summary.txt'), 'w'); fprintf(fid, '%s\n', txt); fclose(fid);
end

function lines = calibrationLines(C)
% the ratio for all pipelines, and with and without each learned step
lines = {pipecompare.eval.Calibration.describe(C,'all compared pipelines'), ...
    pipecompare.eval.Calibration.describe(C(C.selected,:),'selected pipelines'), ...
    '  Records have equal weight; with/without groups below are descriptive, not paired ablations.'};
steps = {'EOG regression', 'eogregress'; 'ICA components chosen by SME', 'threshold=sme'; ...
    'ICA', 'icremove'; 'trial weights', 'weight_trials'; 'stretches interpolated', 'stretches=true'};
for k = 1:size(steps, 1)
    has = contains(C.steps, steps{k, 2});
    if any(has), lines{end+1} = pipecompare.eval.Calibration.describe(C(has,:),steps{k,1}); end %#ok<AGROW>
    if any(~has), lines{end+1} = pipecompare.eval.Calibration.describe(C(~has,:),['without ' steps{k,1}]); end %#ok<AGROW>
end
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
tb = [];
if ~isempty(c.baseline), tb = (ceil(c.baseline(1)*fs-1e-6):floor(c.baseline(2)*fs+1e-6))/fs; end
baseline = zeros(size(truth.amp));
if ~isempty(tb), baseline = arrayfun(@(k) mean(exp(-0.5*((tb-truth.latency(k))/truth.sigma).^2)),1:numel(truth.amp)); end
p3 = mean(topo(roi)) * truth.amp .* (arrayfun(@(k) mean(exp(-0.5 * ((tw - truth.latency(k)) / truth.sigma) .^ 2)), 1:numel(truth.amp))-baseline);
want = arrayfun(@(i) mean(p3(ismember(truth.types, cellstr(c.conditions(i).events)))), 1:numel(c.conditions));
end

function want = populationP3(truth, c, fs)
comp = c.components(1); P = truth.population;
roi = ismember(upper(truth.labels),upper(cellstr(comp.roi)));
eeg = ~startsWith(upper(truth.labels),'EOG'); topo = truth.topo-mean(truth.topo(eeg));
tw = (ceil(comp.window(1)*fs-1e-6):floor(comp.window(2)*fs+1e-6))/fs;
sigma = sqrt(truth.sigma^2+P.latencySd^2);
shape = mean(truth.sigma/sigma*exp(-0.5*((tw-P.latencyMean)/sigma).^2));
if ~isempty(c.baseline)
    tb = (ceil(c.baseline(1)*fs-1e-6):floor(c.baseline(2)*fs+1e-6))/fs;
    shape = shape-mean(truth.sigma/sigma*exp(-0.5*((tb-P.latencyMean)/sigma).^2));
end
want = zeros(1,numel(c.conditions));
for k = 1:numel(want)
    assert(numel(c.conditions(k).events)==1,'PipeCompare:Benchmark','Population truth needs one event type per condition.');
    assert(any(strcmp(c.conditions(k).events{1},{'11','31'})),'PipeCompare:Benchmark','Unknown population event type.');
    if strcmp(c.conditions(k).events{1},'11'), a=P.targetAmp; else, a=P.standardAmp; end
    want(k) = mean(topo(roi))*a*shape;
end
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

