function y = fftfilt(b, x, n)
%FFTFILT FIR filtering by overlap-add in the frequency domain.
%   y = fftfilt(b, x) is filter(b, 1, x) (columns of x filtered; a row
%   vector as one signal), computed with FFTs. PipeCompare adds this folder
%   to the path only while one of its own long high-pass filters runs on a
%   computer without the Signal Processing Toolbox, whose fftfilt it stands
%   in for there (EEGLAB's firfilt calls fftfilt when 'usefftfilt' is on):
%   the same result to rounding, in seconds instead of minutes for a filter
%   with thousands of taps. n (optional): the FFT length to start from.
b = double(b(:)); L = numel(b);
row = isrow(x);
if row, x = x(:); end
x = double(x);
N = size(x, 1);
if nargin < 3 || isempty(n), n = 4 * L; end
nfft = 2 ^ nextpow2(max(n, L + 1));
nfft = min(nfft, 2 ^ nextpow2(N + L - 1));   % no longer than the whole signal needs
M = nfft - L + 1;                             % new samples per block
B = fft(b, nfft);
y = zeros(size(x));
for k = 1:M:N
    last = min(k + M - 1, N);
    Y = ifft(fft(x(k:last, :), nfft) .* B, 'symmetric');
    to = min(k + nfft - 1, N);
    y(k:to, :) = y(k:to, :) + Y(1:to - k + 1, :);
end
if row, y = y.'; end
end
