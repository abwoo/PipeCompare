function [EEG, truth] = nqc_synth(opts)
%NQC_SYNTH Continuous synthetic EEGLAB dataset (microvolts) with known ERP.
%   Two conditions ('11' target, '31' standard). A parietal P3 (Gaussian,
%   peak at truth.p3Latency) is added to every target trial (smaller to
%   standards). Noise: 1/f background, slow drift, 50 Hz line noise, blinks
%   and, on a known subset of trials, large movement artifacts. Channel
%   locations come from EEGLAB's sample 32-channel montage. p3AmpSd and
%   p3Jitter make the P3's size (relative SD) and latency (SD, s) differ
%   from trial to trial (truth.p3Amp, truth.p3Shift per trial).
%   lockedBlinks: the share of target trials with a blink (lockedBlinkUv
%   at FPz) at the P3's latency, time-locked to the stimulus
%   (truth.lockedBlink per trial).
%   lockedSaccades: the share of target trials with a leftward eye
%   movement (lockedSaccadeUv, left minus right EOG) held from 0.25 to
%   0.55 s; lockedEmg: the share of all trials with a muscle burst
%   (lockedEmgUv, fast noise, near F3) from 0.25 to 0.55 s
%   (truth.lockedSaccade, truth.lockedEmg). Each is drawn from its own
%   random stream, so turning one off leaves the others as they were.
if nargin < 1, opts = struct(); end
d = struct('srate', 250, 'seconds', 400, 'nPerCond', 90, 'seed', 7, 'p3', 6, ...
    'noise', 8, 'drift', 25, 'line', 4, 'blinkRate', 0.25, 'artifactTrials', 0.15, ...
    'p3Jitter', 0, 'p3AmpSd', 0, 'noisyChannels', {{}}, 'noisyUv', 60, 'alphaUv', 0, ...
    'lockedBlinks', 0, 'lockedBlinkUv', 200, 'lockedSaccades', 0, 'lockedSaccadeUv', 100, ...
    'lockedEmg', 0, 'lockedEmgUv', 30);
for f = fieldnames(d)', if ~isfield(opts, f{1}), opts.(f{1}) = d.(f{1}); end, end
rng(opts.seed, 'twister');
locfile = fullfile(fileparts(which('eeglab')), 'sample_data', 'eeglab_chan32.locs');
chanlocs = readlocs(locfile);
nch = numel(chanlocs); fs = opts.srate; n = round(opts.seconds * fs);
xyz = [[chanlocs.X]' [chanlocs.Y]' [chanlocs.Z]']; xyz = xyz ./ vecnorm(xyz, 2, 2);
lab = upper({chanlocs.labels});
topo = @(name, width) exp(-(acos(max(-1, min(1, xyz * xyz(strcmp(lab, upper(name)), :)')))).^2 / (2 * width^2))';

% background: 1/f-ish + drift + line noise
data = zeros(nch, n);
for c = 1:nch
    w = randn(1, n);
    x = filter(1, [1 -0.95], w);
    data(c, :) = opts.noise * x / std(x);
end
t = (0:n-1) / fs;
for c = 1:nch
    data(c, :) = data(c, :) + opts.drift * sin(2*pi*(0.02 + 0.03*rand) * t + 2*pi*rand);
end
data = data + opts.line * sin(2*pi*50*t);

% blinks (frontal), Poisson times
blinkTopo = topo('FPZ', 0.6);
nb = round(opts.blinkRate * opts.seconds);
bt = sort(randi([fs n - fs], 1, nb));
bw = -round(0.2*fs):round(0.2*fs);
blink = 120 * exp(-0.5 * (bw / (0.05*fs)).^2);
for k = 1:nb
    data(:, bt(k) + bw) = data(:, bt(k) + bw) + blinkTopo' * blink;
end

% events and ERP
nTrials = 2 * opts.nPerCond;
types = [repmat({'11'}, 1, opts.nPerCond) repmat({'31'}, 1, opts.nPerCond)];
types = types(randperm(nTrials));
isi = (opts.seconds - 10) / nTrials;
lat = round(fs * (5 + (0:nTrials-1) * isi + 0.2 * rand(1, nTrials)));
p3Topo = topo('PZ', 0.5);
truth.p3Latency = 0.40; truth.p3Sigma = 0.07;
tt = -round(0.2*fs):round(1.0*fs);
wave = exp(-0.5 * ((tt/fs - truth.p3Latency) / truth.p3Sigma).^2);
art = false(1, nTrials); truth.p3Amp = zeros(1, nTrials); truth.p3Shift = zeros(1, nTrials);
for k = 1:nTrials
    amp = opts.p3; if strcmp(types{k}, '31'), amp = opts.p3 / 3; end
    if opts.p3AmpSd > 0, amp = amp * max(0, 1 + opts.p3AmpSd * randn); end   % the P3 differs from trial to trial
    if opts.p3Jitter > 0
        truth.p3Shift(k) = opts.p3Jitter * randn;
        wk = exp(-0.5 * ((tt/fs - truth.p3Latency - truth.p3Shift(k)) / truth.p3Sigma).^2);
    else
        wk = wave;
    end
    truth.p3Amp(k) = amp;
    data(:, lat(k) + tt) = data(:, lat(k) + tt) + p3Topo' * (amp * wk);
    if rand < opts.artifactTrials
        art(k) = true; % movement artifact: large, broad, all channels
        seg = round(fs * (0.1 + 0.5*rand)) + (0:round(0.3*fs));
        bump = 150 * sin(pi * (0:numel(seg)-1) / numel(seg));
        data(:, lat(k) + seg) = data(:, lat(k) + seg) + (0.5 + rand(nch, 1)) * bump;
    end
end
if opts.alphaUv > 0   % occipital 10 Hz rhythm with slow amplitude modulation
    a = opts.alphaUv * (1 + 0.3 * sin(2*pi*0.1*t)) .* sin(2*pi*10*t);
    data = data + topo('OZ', 0.5)' * a;
end
for k = 1:numel(opts.noisyChannels)   % broken electrodes: large broadband noise
    ch = strcmpi(lab, opts.noisyChannels{k});
    data(ch, :) = data(ch, :) + opts.noisyUv * randn(1, n);
end
truth.lockedBlink = false(1, nTrials);
if opts.lockedBlinks > 0   % after the rest: the other draws stay the same
    truth.lockedBlink = rand(1, nTrials) < opts.lockedBlinks & strcmp(types, '11');
    for k = find(truth.lockedBlink)
        at = lat(k) + round(truth.p3Latency * fs);
        data(:, at + bw) = data(:, at + bw) + blinkTopo' * (opts.lockedBlinkUv / 120 * blink);
    end
end
truth.lockedSaccade = false(1, nTrials);
if opts.lockedSaccades > 0
    rs = RandStream('mt19937ar', 'Seed', opts.seed + 101);
    truth.lockedSaccade = rand(rs, 1, nTrials) < opts.lockedSaccades & strcmp(types, '11');
    sTopo = topo('EOG2', 0.5) - topo('EOG1', 0.5);   % EOG2 is on the left, EOG1 on the right
    hold = round(0.25 * fs):round(0.55 * fs);
    step = ones(1, numel(hold)); r = min(numel(hold), round(0.02 * fs));
    step(1:r) = (1:r) / r; step(end-r+1:end) = (r:-1:1) / r;
    for k = find(truth.lockedSaccade)
        data(:, lat(k) + hold) = data(:, lat(k) + hold) + sTopo' * (opts.lockedSaccadeUv * step);
    end
end
truth.lockedEmg = false(1, nTrials);
if opts.lockedEmg > 0
    rs = RandStream('mt19937ar', 'Seed', opts.seed + 202);
    truth.lockedEmg = rand(rs, 1, nTrials) < opts.lockedEmg;
    eTopo = topo('F3', 0.3);
    hold = round(0.25 * fs):round(0.55 * fs);
    for k = find(truth.lockedEmg)
        x = diff(randn(rs, 1, numel(hold) + 1)) / sqrt(2);   % fast: most of its power above 30 Hz
        data(:, lat(k) + hold) = data(:, lat(k) + hold) + eTopo' * (opts.lockedEmgUv * x);
    end
end
truth.alphaUv = opts.alphaUv;
EEG = eeg_emptyset();
EEG.setname = 'nqc_synth'; EEG.srate = fs; EEG.data = single(data);
EEG.nbchan = nch; EEG.pnts = n; EEG.trials = 1; EEG.xmin = 0; EEG.xmax = (n-1)/fs;
EEG.chanlocs = chanlocs;
EEG.event = struct('type', types, 'latency', num2cell(lat), 'duration', 0);
EEG = eeg_checkset(EEG, 'eventconsistency');
EEG = eeg_checkset(EEG, 'makeur');
EEG.history = sprintf('EEG = pop_loadset(''nqc_synth.set''); %% synthetic');
truth.p3Topo = p3Topo; truth.artifactTrials = art; truth.types = types; truth.roi = {'Pz','CPz','POz','P3','P4'};
truth.roi = truth.roi(ismember(upper(truth.roi), lab));
end
