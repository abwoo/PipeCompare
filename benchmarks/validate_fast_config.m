function cfg = validate_fast_config(cfg)
%VALIDATE_FAST_CONFIG Frozen fast-mode recipe check (diagnostic only).
%
% This function validates a recipe, it does NOT run a search and it
% does NOT change any pick. A passing validation is not evidence that
% a search is fast, that a pick is good, or that a reduced family
% exists: see the executionReady flag below.
%
% Frozen recipe (Round 10 fast mode, design-seed diagnostic only):
%   seed in 40101-40104 (design seeds, user instruction 2026-10-28
%   23:45: fast tests and micro searches use the existing 4000-series
%   sample data to FIND problems; a design seed is never held-out
%   evidence and a result on it must never be reported as a held-out
%   verdict);
%   nSeeds == 1, nStrata == 1, armMode in {'standard','balanced'}
%   (never 'all': two searches break the 1-hour budget);
%   checkpointDir non-empty, maxSeconds finite and <= 2400,
%   sliceSeconds == 2400 (40 minutes);
%   icaSeeds == 1, refineCutoffs == false, eogCheck == false;
%   reducedFamily and reducedTrials are non-empty text describing
%   the intended reduction rule (the rule itself is recorded, see
%   the execution gap below).
%
% Execution gap, stated plainly (do not remove this comment without
% landing the Presets threading it describes): at head 6035b50,
% pop_pipecompare's documented option list does NOT include
% icaSeeds, refineCutoffs or eogCheck, and
% Presets.searchOptions() hardcodes refineCutoffs=true, icaSeeds=2,
% eogCheck=true. This package therefore does NOT thread those
% overrides through pop_pipecompare, and a caller must NOT claim
% that a search run under this validation used a reduced family.
% cfg.executionReady is therefore always false in this package, and
% cfg.sliceDiagnosticReady is true only for the slice-diagnostic
% run: one design seed, one arm, production family, checkpoint +
% 40-minute slice, labelled slice-diagnostic, never held-out.
% Reduced-family execution is an explicit unfinished item (see
% docs/PR50_STATUS.md); it needs Presets threading plus pre-registration
% of the reduction rule before any reduced-family run is claimed.
if ~isstruct(cfg)
    error('PipeCompare:FastMode', 'fast config must be a struct.');
end
req = {'seed', 'nSeeds', 'nStrata', 'armMode', 'checkpointDir', ...
    'maxSeconds', 'sliceSeconds', 'icaSeeds', 'refineCutoffs', ...
    'eogCheck', 'reducedFamily', 'reducedTrials'};
for k = 1:numel(req)
    if ~isfield(cfg, req{k})
        error('PipeCompare:FastMode', 'fast config is missing field %s.', req{k});
    end
end
if ~isscalar(cfg.seed) || ~ismember(cfg.seed, 40101:40104)
    error('PipeCompare:FastMode', 'micro mode accepts design seeds 40101-40104 only; seed %s is not a design-seed diagnostic.', num2str(cfg.seed));
end
if cfg.nSeeds ~= 1
    error('PipeCompare:FastMode', 'micro mode requires exactly 1 seed, got %d.', cfg.nSeeds);
end
if cfg.nStrata ~= 1
    error('PipeCompare:FastMode', 'micro mode requires exactly 1 stratum, got %d.', cfg.nStrata);
end
if ~ismember(cfg.armMode, {'standard', 'balanced'})
    error('PipeCompare:FastMode', 'micro mode requires one arm (standard or balanced), not ''all''.');
end
if ~ischar(cfg.checkpointDir) || isempty(cfg.checkpointDir)
    error('PipeCompare:FastMode', 'micro mode requires a non-empty checkpointDir; empty checkpoint plus Inf time is forbidden.');
end
if ~isscalar(cfg.maxSeconds) || ~isfinite(cfg.maxSeconds) || cfg.maxSeconds <= 0 || cfg.maxSeconds > 2400
    error('PipeCompare:FastMode', 'micro mode requires finite maxSeconds in (0, 2400], got %s.', num2str(cfg.maxSeconds));
end
if cfg.sliceSeconds ~= 2400
    error('PipeCompare:FastMode', 'micro mode requires sliceSeconds == 2400 (40 minutes); changing it is a new recipe, not micro.');
end
if cfg.icaSeeds ~= 1 || ~isequal(cfg.refineCutoffs, false) || ~isequal(cfg.eogCheck, false)
    error('PipeCompare:FastMode', 'micro recipe requires icaSeeds==1, refineCutoffs==false, eogCheck==false (diagnostic only).');
end
if ~ischar(cfg.reducedFamily) || isempty(cfg.reducedFamily) || ~ischar(cfg.reducedTrials) || isempty(cfg.reducedTrials)
    error('PipeCompare:FastMode', 'micro recipe requires non-empty reducedFamily and reducedTrials rule text.');
end
cfg.mode = 'micro';
cfg.designSeedDiagnostic = true;
cfg.executionReady = false;
cfg.sliceDiagnosticReady = true;
cfg.unwired = {'presetsSearchOptionsOverride', 'trialSubsetOverride'};
cfg.evidence = 'diagnostic only: never a held-out verdict, never merged into a paired mean, never compared by MSE with a full-family pick.';
end
