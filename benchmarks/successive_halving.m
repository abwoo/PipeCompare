function [survivors, audit] = successive_halving(candidateIds, fidelityScores, eta, keepFractions)
%SUCCESSIVE_HALVING Successive-halving prototype, NOT wired in.
%   STATUS: Round 8 bundle, PROTOTYPE / DESIGN, no MATLAB search wiring.
%   candidateIds: vector of candidate indices.
%   fidelityScores: cell array, one score vector per fidelity rung
%     (low fidelity first, full fidelity last). Scores are
%     lower-is-better and must be computed by the caller; this file
%     does NOT compute any EEG score and does NOT call Rank.
%   eta: promotion divisor (e.g. 3); keepFractions: kept fraction per
%     rung, or [] to use 1/eta per rung.
%   RULES (frozen before any held-out use, see MANIFEST.md):
%     - Low-fidelity scores are NEVER the final score. Survivors of
%       the last rung must still be ranked by full-fidelity Rank.run.
%     - The full-fidelity best must be reported as survived/not, per
%       held-out seed (survival rate). A halving rule that silently
%       drops the true best is a failure, not a speedup.
%     - eta, rung definitions and keep fractions are frozen on design
%       seeds before held-out evaluation; tuning them on held-out
%       results is double use.
%   This file is not called by Executor, Rank or the shadow script
%   in Round 8. Wiring it as a new arm requires the gates in
%   MANIFEST.md (paired MSE, survival, wall time, pick stability).
if nargin < 4 || isempty(keepFractions)
    keepFractions = repmat(1 / eta, 1, numel(fidelityScores));
end
survivors = candidateIds(:)';
audit = struct('rung', {}, 'nBefore', {}, 'nAfter', {});
for r = 1:numel(fidelityScores)
    scores = fidelityScores{r};
    if numel(scores) ~= numel(survivors)
        error('PipeCompare:Halving', 'Fidelity rung %d scores do not match survivors.', r);
    end
    [~, order] = sort(scores, 'ascend');
    nKeep = max(1, floor(numel(survivors) * keepFractions(r)));
    survivors = survivors(order(1:nKeep));
    audit(end+1) = struct('rung', r, 'nBefore', numel(scores), 'nAfter', numel(survivors)); %#ok<AGROW>
end
end
