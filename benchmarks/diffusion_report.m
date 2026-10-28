function rep = diffusion_report(cands, pick)
%DIFFUSION_REPORT Report-only reference-diffusion exposure for a pick.
%   STATUS: Round 8 bundle, REPORT-ONLY, not a scoring rule.
%   Returns the selected candidate's interpolatedFraction (0-1) when
%   available, otherwise NaN with a flag. Never 0 for "missing":
%   missing is unknown, not "no interpolation". This column is
%   report-first (step 2 design): it must not enter Rank scoring or
%   any gate until its sign and held-out calibration are frozen.
%   Not wired into the shadow run in Round 8; the shadow run JSON
%   field for this column is a follow-up wiring step, listed in
%   MANIFEST.md as design/prototype, not as shipped evidence.
rep = struct('interpolatedFraction', NaN, 'available', false, 'pick', pick);
if isempty(pick) || ~isscalar(pick) || pick < 1 || pick > numel(cands)
    return;
end
c = cands(pick);
if isfield(c, 'interpolatedFraction') && ~isempty(c.interpolatedFraction) && isfinite(c.interpolatedFraction)
    rep.interpolatedFraction = c.interpolatedFraction;
    rep.available = true;
elseif isfield(c, 'm') && isfield(c.m, 'interpolatedFraction')
    v = c.m.interpolatedFraction;
    if ~isempty(v) && isfinite(v)
        rep.interpolatedFraction = v;
        rep.available = true;
    end
end
end
