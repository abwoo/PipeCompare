function [EEG, truth] = bench_synth(seed, opts)
%BENCH_SYNTH Simulated oddball recording with a known P3, for the benchmark.
%   [EEG, truth] = bench_synth(seed, opts)
%   Background: 24 cortical sources (smooth scalp fields, 1/f activity) on
%   EEGLAB's sample 32-channel montage, plus sensor noise, an occipital
%   alpha rhythm, slow drifts and 50 Hz line noise. Artifacts: random
%   blinks, blinks time-locked to target stimuli inside the P3 window,
%   temporal muscle bursts, movement artifacts in some trials, and bad
%   channels with broadband noise of their own. P3: parietal (Pz), its
%   size and latency differ from trial to trial. truth holds what is
%   needed to compute the P3 actually in the data for any window and
%   electrodes (pipecompare_benchmark), and the artifacts' positions.
%   opts (all optional): srate, seconds, nPerCond, p3 (target amplitude,
%   uV; standards a third), lockedBlinks (share of target trials),
%   blinkRate (per s), badChannels (count), muscleShare (of the time),
%   artifactTrials (share).
if nargin < 2, opts = struct(); end
d = struct('srate', 250, 'seconds', 420, 'nPerCond', 75, 'p3', 6, 'lockedBlinks', 0.25, ...
    'blinkRate', 0.2, 'badChannels', 1, 'muscleShare', 0.05, 'artifactTrials', 0.1);
for f = fieldnames(d)', if ~isfield(opts, f{1}), opts.(f{1}) = d.(f{1}); end, end
rng(seed, 'twister');
chanlocs = readlocs(fullfile(fileparts(which('eeglab')), 'sample_data', 'eeglab_chan32.locs'));
lab = upper({chanlocs.labels});
eog = find(startsWith(lab, 'EOG'));
[chanlocs.type] = deal('EEG'); [chanlocs(eog).type] = deal('EOG');
nch = numel(chanlocs); fs = opts.srate; n = round(opts.seconds * fs);
xyz = [[chanlocs.X]' [chanlocs.Y]' [chanlocs.Z]']; xyz = xyz ./ vecnorm(xyz, 2, 2);
field = @(v, w) exp(-acos(max(-1, min(1, xyz * (v(:) / norm(v))))) .^ 2 / (2 * w ^ 2))';
at = @(name) xyz(strcmp(lab, upper(name)), :);
t = (0:n-1) / fs;

% background: cortical sources seen through the skull (smooth fields)
K = 24; S = randn(K, 3); S(:, 3) = abs(S(:, 3));
A = zeros(nch, K);
for k = 1:K, A(:, k) = (sign(randn) * field(S(k, :), 0.45))'; end
data = A * band(randn(K, n), fs, 0.3, 45, true);
data = 10 * data / median(std(data, 0, 2));
data = data + 2 * randn(nch, n);                                         % sensor noise
alpha = band(randn(1, n), fs, 8.5, 11.5, false); alpha = 8 * alpha / std(alpha);
data = data + field(at('OZ'), 0.5)' * alpha;
for c = 1:nch                                                            % slow drifts
    data(c, :) = data(c, :) + 20 * sin(2*pi*(0.01 + 0.04*rand) * t + 2*pi*rand);
end
data = data + (3 + 3 * rand(nch, 1)) * sin(2*pi*50*t);

% random blinks (frontal)
blinkTopo = field(at('FPZ'), 0.6);
bw = -round(0.2*fs):round(0.2*fs); blink = exp(-0.5 * (bw / (0.05*fs)) .^ 2);
bt = sort(randi([fs n - fs], 1, round(opts.blinkRate * opts.seconds)));
for k = 1:numel(bt), data(:, bt(k) + bw) = data(:, bt(k) + bw) + blinkTopo' * (150 * blink); end

% events and the P3
nTrials = 2 * opts.nPerCond;
types = [repmat({'11'}, 1, opts.nPerCond) repmat({'31'}, 1, opts.nPerCond)];
types = types(randperm(nTrials));
isi = (opts.seconds - 10) / nTrials;
lat = round(fs * (5 + (0:nTrials-1) * isi + 0.3 * rand(1, nTrials)));
truth.topo = field(at('PZ'), 0.5); truth.sigma = 0.08;
truth.amp = zeros(1, nTrials); truth.latency = zeros(1, nTrials);
tt = -round(0.2*fs):round(1.0*fs);
truth.lockedBlink = false(1, nTrials); truth.artifact = false(1, nTrials);
for k = 1:nTrials
    a = opts.p3; if strcmp(types{k}, '31'), a = opts.p3 / 3; end
    truth.amp(k) = a * max(0, 1 + 0.3 * randn);
    truth.latency(k) = 0.40 + 0.04 * randn;
    wave = exp(-0.5 * ((tt/fs - truth.latency(k)) / truth.sigma) .^ 2);
    data(:, lat(k) + tt) = data(:, lat(k) + tt) + truth.topo' * (truth.amp(k) * wave);
    if strcmp(types{k}, '11') && rand < opts.lockedBlinks              % a blink in the P3 window
        truth.lockedBlink(k) = true;
        b0 = lat(k) + round((0.30 + 0.25 * rand) * fs);
        data(:, b0 + bw) = data(:, b0 + bw) + blinkTopo' * (150 * blink);
    end
    if rand < opts.artifactTrials                                       % movement: large, broad, all channels
        truth.artifact(k) = true;
        seg = round(fs * (0.1 + 0.5 * rand)) + (0:round(0.3*fs));
        data(:, lat(k) + seg) = data(:, lat(k) + seg) + (0.5 + rand(nch, 1)) * (150 * sin(pi * (0:numel(seg)-1) / numel(seg)));
    end
end

% temporal muscle bursts
emg = band(randn(1, n), fs, 20, 100, false); emg = 15 * emg / std(emg);
on = false(1, n); m = 0;
while m < opts.muscleShare * n
    len = round(fs * (0.5 + 1.5 * rand)); s0 = randi(n - len);
    on(s0 + (0:len-1)) = true; m = nnz(on);
end
side = {'T7', 'T8'};
for s = 1:2, data = data + field(at(side{s}), 0.35)' * (emg .* on .* (0.5 + rand)); end

% bad channels: broadband noise of their own (poor contact)
cand = find(~ismember(lab, [{'PZ', 'FPZ'} lab(eog)]));
truth.bad = cand(randperm(numel(cand), opts.badChannels));
for c = truth.bad, data(c, :) = data(c, :) + 30 * randn(1, n); end

EEG = eeg_emptyset();
EEG.setname = sprintf('bench_synth_%d', seed); EEG.srate = fs; EEG.data = single(data);
EEG.nbchan = nch; EEG.pnts = n; EEG.trials = 1; EEG.xmin = 0; EEG.xmax = (n-1)/fs;
EEG.chanlocs = chanlocs;
EEG.event = struct('type', types, 'latency', num2cell(lat), 'duration', 0);
EEG = eeg_checkset(EEG, 'eventconsistency');
EEG = eeg_checkset(EEG, 'makeur');
EEG.history = sprintf('EEG = bench_synth(%d); %% simulated', seed);
truth.types = types; truth.labels = {chanlocs.labels}; truth.badLabels = {chanlocs(truth.bad).labels};
truth.lockedBlinkShare = mean(truth.lockedBlink(strcmp(types, '11')));
end

function x = band(x, fs, lo, hi, pink)
% x (rows) restricted to lo-hi Hz in the frequency domain; pink: 1/f power
n = size(x, 2); F = fft(x, [], 2);
f = (0:n-1) * fs / n; f = min(f, fs - f);
g = double(f >= lo & f <= hi);
if pink, g = g ./ sqrt(max(f, lo)); end
x = real(ifft(F .* g, [], 2));
end
