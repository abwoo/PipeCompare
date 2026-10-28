function [keep, undefined] = bestset_strict_prototype(pts, Bq, ib, alpha)
%BESTSET_STRICT_PROTOTYPE Strict best-set prototype, NOT wired in.
%   STATUS: Round 8 bundle, PROTOTYPE / SEPARATE ARM ONLY.
%   When the simultaneous critical value is infinite (incomplete
%   bootstrap support), the production bestSet keeps every candidate
%   (keep-all) and Rank.run in Round 8 only FLAGS this as
%   bestSetUndefined (report-only, pick unchanged). This prototype
%   shows the strict alternative: keep only the point-best and flag
%   undefined=true. It is NOT called by Rank.run in Round 8, because
%   changing keep would change pick and would need its own held-out
%   arm with paired MSE, calibration and pick-stability gates.
%   Do not wire it in without that arm and without updating METHODS.
pts = pts(:);
K = numel(pts);
keep = false(K, 1);
keep(ib) = true;
undefined = false;
B = size(Bq, 2);
if B == 0
    undefined = true;
    return;
end
complete = all(isfinite(Bq), 1);
maxD = -inf(1, B);
Dhat = pts - pts';
for i = 1:K
    D = Bq(i, :) - Bq;
    D(~isfinite(D)) = NaN;
    maxD = max(maxD, max(D - Dhat(i, :)', [], 1, 'omitnan'));
end
maxD(~complete) = Inf;
% Upper quantile, same percentile rule family as Rank.upperQuantile
% (implemented here standalone so this file has no Rank dependency).
s = sort(maxD);
at = max(1, min(B, ceil((B + 1) * (1 - alpha))));
c = s(at);
if isinf(c)
    undefined = true;
    return;
end
worse = Dhat > c;
keep = keep | ~any(worse, 2);
end
