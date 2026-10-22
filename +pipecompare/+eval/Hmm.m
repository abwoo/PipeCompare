classdef Hmm
    % Two-state hidden Markov models (quiet / active) with Gaussian
    % emissions: the parameters fitted to the data by EM (Baum-Welch),
    % the states found by Viterbi (the most probable sequence, exactly,
    % by dynamic programming).

    methods (Static)
        function [active, m] = segment(y, peak, minLen)
            % The samples of the row y in the active state (the state with
            % the higher mean). A stretch of active samples counts only when
            % y reaches peak in it and it is minLen samples or longer, so
            % that the noise of data without any events is not split into
            % two states. m: the fitted model (mu, sd, A, pi, valid, reason).
            y = double(y(:)');
            n = numel(y);
            active = false(1, n);
            m = struct('mu', [], 'sd', [], 'A', [], 'pi', [], 'valid', false, 'reason', 'insufficient finite samples');
            if n < 10 || ~any(isfinite(y)), return; end
            y(~isfinite(y)) = median(y(isfinite(y)));
            q = sort(y);
            mu = [q(ceil(0.5 * n)) q(ceil(0.99 * n))];
            sd = [std(y) / 2 std(y)] + eps;
            A = [0.98 0.02; 0.1 0.9]; p0 = [0.9 0.1];
            if any(~isfinite(sd)), m.reason = 'nonfinite initial scale'; return; end
            ll = -inf;
            for it = 1:30
                [g1, xi, llNew] = forwardBackward(y, A, mu, sd, p0);
                if ~isfinite(llNew) || any(~isfinite(g1)) || any(~isfinite(xi(:)))
                    m.reason = 'nonfinite forward-backward fit'; return;
                end
                g2 = 1 - g1;
                % Positive transition probabilities also protect a state
                % with negligible occupancy from becoming a zero row.
                xi = xi + 1e-12;
                A = xi ./ sum(xi, 2);
                mu = [sum(g1 .* y) / max(sum(g1), eps), sum(g2 .* y) / max(sum(g2), eps)];
                sd = sqrt([sum(g1 .* (y - mu(1)).^2) / max(sum(g1), eps), sum(g2 .* (y - mu(2)).^2) / max(sum(g2), eps)]) + 1e-6;
                p0 = [g1(1) g2(1)] + 1e-6; p0 = p0 / sum(p0);
                if llNew - ll < 1e-6 * abs(llNew), break; end
                ll = llNew;
            end
            if any(~isfinite([mu sd A(:)' p0])) || any(sd <= 0)
                m.reason = 'nonfinite fitted parameters'; return;
            end
            m = struct('mu', mu, 'sd', sd, 'A', A, 'pi', p0, 'valid', true, 'reason', '');
            [~, hi] = max(mu);
            active = viterbi(y, A, mu, sd, p0) == hi;
            d = diff([0 active 0]); s = find(d == 1); e = find(d == -1) - 1;
            for k = 1:numel(s)
                if e(k) - s(k) + 1 < minLen || max(y(s(k):e(k))) < peak, active(s(k):e(k)) = false; end
            end
        end

        function [hit, model] = events(X, srate)
            % Blinks and eye movements in the rows of X (eye channels, or a
            % left minus right pair; continuous, srate Hz), found by a two-
            % state HMM: each row minus its moving median over 1 s, smoothed
            % over 1/16 s, taken at about 25 Hz; the largest over the rows of
            % its robust z (median, MAD) of the absolute value. A blink's
            % samples add their evidence together, so weak blinks that never
            % reach 5 robust SDs in one sample are found as well. A stretch
            % counts when it reaches PeakZ and lasts MinSeconds. hit: the
            % samples (logical row).
            n = size(X, 2);
            hit = false(1, n);
            model = struct('mu',[],'sd',[],'A',[],'pi',[],'valid',false,'reason','insufficient finite samples');
            if isempty(X) || n < 2 * srate, return; end
            d = max(1, floor(srate / 25));
            Z = zeros(size(X, 1), numel(1:d:n));
            for k = 1:size(X, 1)
                r = double(X(k, :)) - movmedian(double(X(k, :)), round(srate), 2);
                s = movmean(r, max(1, round(srate / 16)), 2);
                s = s(1:d:end);
                a = abs(s - median(s));
                Z(k, :) = a / max(1.4826 * median(a), eps);
            end
            [on, model] = pipecompare.eval.Hmm.segment(max(Z, [], 1), pipecompare.eval.Hmm.PeakZ, ...
                max(1, round(pipecompare.eval.Hmm.MinSeconds * srate / d)));
            hit = repelem(on, d); hit = hit(1:n);
        end
    end

    properties (Constant)
        PeakZ = 3.5         % a stretch counts when it reaches 3.5 robust SDs (5 for single samples before)
        MinSeconds = 0.1    % and lasts 100 ms or more (blinks last 200-400 ms)
    end
end

function [g1, xi, ll] = forwardBackward(y, A, mu, sd, p0)
% Scaled forward-backward for two states: g1 the probability of state 1 at
% each sample, xi the expected transitions, ll the log-likelihood.
n = numel(y);
L1 = -0.5 * ((y - mu(1)) / sd(1)).^2 - log(sd(1));
L2 = -0.5 * ((y - mu(2)) / sd(2)).^2 - log(sd(2));
% Scale BEFORE exponentiation: dividing two already-underflowed Gaussian
% densities by realmin leaves both zero and poisons the entire EM fit.
s = max(L1, L2); L1 = exp(L1 - s); L2 = exp(L2 - s);
a1 = zeros(1, n); a2 = zeros(1, n); c = zeros(1, n);
v1 = p0(1) * L1(1); v2 = p0(2) * L2(1); c(1) = v1 + v2; a1(1) = v1 / c(1); a2(1) = v2 / c(1);
A11 = A(1, 1); A12 = A(1, 2); A21 = A(2, 1); A22 = A(2, 2);
for t = 2:n
    v1 = (a1(t-1) * A11 + a2(t-1) * A21) * L1(t);
    v2 = (a1(t-1) * A12 + a2(t-1) * A22) * L2(t);
    c(t) = v1 + v2; a1(t) = v1 / c(t); a2(t) = v2 / c(t);
end
b1 = 1; b2 = 1; g1 = zeros(1, n); g1(n) = a1(n);
x11 = 0; x12 = 0; x21 = 0; x22 = 0;
for t = n-1:-1:1
    e1 = L1(t+1) * b1; e2 = L2(t+1) * b2;
    t11 = a1(t) * A11 * e1; t12 = a1(t) * A12 * e2; t21 = a2(t) * A21 * e1; t22 = a2(t) * A22 * e2;
    z = t11 + t12 + t21 + t22;
    x11 = x11 + t11 / z; x12 = x12 + t12 / z; x21 = x21 + t21 / z; x22 = x22 + t22 / z;
    n1 = A11 * e1 + A12 * e2; n2 = A21 * e1 + A22 * e2; nb = n1 + n2;
    b1 = n1 / nb; b2 = n2 / nb;
    g1(t) = a1(t) * b1 / (a1(t) * b1 + a2(t) * b2);
end
xi = [x11 x12; x21 x22];
ll = sum(log(c)) + sum(s);
end

function path = viterbi(y, A, mu, sd, p0)
% The most probable state sequence (1 or 2), exactly, in log space.
n = numel(y);
L1 = -0.5 * ((y - mu(1)) / sd(1)).^2 - log(sd(1));
L2 = -0.5 * ((y - mu(2)) / sd(2)).^2 - log(sd(2));
lA = log(max(A, realmin));
d1 = log(max(p0(1), realmin)) + L1(1); d2 = log(max(p0(2), realmin)) + L2(1);
from1 = false(1, n); from2 = false(1, n);   % true: came from state 2
for t = 2:n
    u1 = d1 + lA(1, 1); w1 = d2 + lA(2, 1);
    u2 = d1 + lA(1, 2); w2 = d2 + lA(2, 2);
    from1(t) = w1 > u1; from2(t) = w2 > u2;
    d1 = max(u1, w1) + L1(t); d2 = max(u2, w2) + L2(t);
end
path = zeros(1, n); path(n) = 1 + (d2 > d1);
for t = n:-1:2
    if path(t) == 1, path(t-1) = 1 + from1(t); else, path(t-1) = 1 + from2(t); end
end
end

