function [EEG, truth, sources] = bench_synth(seed, opts)
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
%   Optional third output sources separates additive generating components.
%   It does not change the recording, truth, or RNG stream; use for audits,
%   never as input to the operational selector.
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
trace=nargout>2;sources=struct();
if trace
    sources.background=data;
    for name={'blink','p3','movement','muscle','badChannel'},sources.(name{1})=zeros(nch,n);end
end

% random blinks (frontal)
blinkTopo = field(at('FPZ'), 0.6);
bw = -round(0.2*fs):round(0.2*fs); blink = exp(-0.5 * (bw / (0.05*fs)) .^ 2);
bt = sort(randi([fs n - fs], 1, round(opts.blinkRate * opts.seconds)));
for k = 1:numel(bt)
    delta=blinkTopo' * (150 * blink);cols=bt(k)+bw;
    data(:,cols)=data(:,cols)+delta;
    if trace,sources.blink(:,cols)=sources.blink(:,cols)+delta;end
end

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
    delta=truth.topo' * (truth.amp(k) * wave);cols=lat(k)+tt;
    data(:,cols)=data(:,cols)+delta;
    if trace,sources.p3(:,cols)=sources.p3(:,cols)+delta;end
    if strcmp(types{k}, '11') && rand < opts.lockedBlinks              % a blink in the P3 window
        truth.lockedBlink(k) = true;
        b0 = lat(k) + round((0.30 + 0.25 * rand) * fs);
        delta=blinkTopo' * (150 * blink);cols=b0+bw;
        data(:,cols)=data(:,cols)+delta;
        if trace,sources.blink(:,cols)=sources.blink(:,cols)+delta;end
    end
    if rand < opts.artifactTrials                                       % movement: large, broad, all channels
        truth.artifact(k) = true;
        seg = round(fs * (0.1 + 0.5 * rand)) + (0:round(0.3*fs));
        delta=(0.5 + rand(nch, 1)) * (150 * sin(pi * (0:numel(seg)-1) / numel(seg)));cols=lat(k)+seg;
        data(:,cols)=data(:,cols)+delta;
        if trace,sources.movement(:,cols)=sources.movement(:,cols)+delta;end
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
for s = 1:2
    delta=field(at(side{s}),0.35)' * (emg .* on .* (0.5 + rand));
    data=data+delta;
    if trace,sources.muscle=sources.muscle+delta;end
end

% bad channels: broadband noise of their own (poor contact)
cand = find(~ismember(lab, [{'PZ', 'FPZ'} lab(eog)]));
truth.bad = cand(randperm(numel(cand), opts.badChannels));
for c = truth.bad
    delta=30*randn(1,n);data(c,:)=data(c,:)+delta;
    if trace,sources.badChannel(c,:)=sources.badChannel(c,:)+delta;end
end

EEG = eeg_emptyset();
EEG.setname = sprintf('bench_synth_%d', seed); EEG.srate = fs; EEG.data = single(data);
EEG.nbchan = nch; EEG.pnts = n; EEG.trials = 1; EEG.xmin = 0; EEG.xmax = (n-1)/fs;
if trace
    total=sources.background+sources.blink+sources.p3+sources.movement+sources.muscle+sources.badChannel;
    sources.rounding=double(EEG.data)-total;
end
EEG.chanlocs = chanlocs;
EEG.event = struct('type', types, 'latency', num2cell(lat), 'duration', 0);
EEG = eeg_checkset(EEG, 'eventconsistency');
EEG = eeg_checkset(EEG, 'makeur');
EEG.history = sprintf('EEG = bench_synth(%d); %% simulated', seed);
truth.types = types; truth.labels = {chanlocs.labels}; truth.badLabels = {chanlocs(truth.bad).labels};
truth.lockedBlinkShare = mean(truth.lockedBlink(strcmp(types, '11')));
% Population signal, for calibrating SME (sampling error about a population
% mean). The realized trial mean above is a separate denoising target.
z = 1/0.3; ampFactor = 0.5*erfc(-z/sqrt(2)) + 0.3*exp(-z^2/2)/sqrt(2*pi);
truth.population = struct('targetAmp',opts.p3*ampFactor,'standardAmp',opts.p3*ampFactor/3, ...
    'latencyMean',0.40,'latencySd',0.04);
end

function x = band(x, fs, lo, hi, pink)
% x (rows) restricted to lo-hi Hz in the frequency domain; pink: 1/f power
n = size(x, 2); F = fft(x, [], 2);
f = (0:n-1) * fs / n; f = min(f, fs - f);
g = double(f >= lo & f <= hi);
if pink, g = g ./ sqrt(max(f, lo)); end
x = real(ifft(F .* g, [], 2));
end

