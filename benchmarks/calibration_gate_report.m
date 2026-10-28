function rep = calibration_gate_report(current, baseline)
%CALIBRATION_GATE_REPORT Calibration-gate report, REPORT-ONLY.
%   STATUS: Round 8 bundle, REPORT-ONLY, not a gate.
%   current/baseline: structs from Calibration.summarize (ratio, ci,
%   actualRms, nRecords, coverage). Without a baseline from the same
%   protocol, no upgrade verdict is possible: rep.verdict is
%   'not-evaluated'. This function never blocks a run, never changes
%   a ranking score, and never voids a run. An upgrade gate (a future
%   step, not this bundle) would additionally require: log-distance
%   |log(r)| shrinks by >= 1/2, latest pooled 95% CI fully inside
%   [0.8, 1.25], actualRms worsens by <= 5%, assessable >= 30.
rep = struct('verdict', 'not-evaluated', 'ratio', NaN, 'ci', [NaN NaN], ...
    'coverage', NaN, 'nRecords', 0, 'logDistanceShrink', NaN, 'actualRmsWorsening', NaN);
if nargin < 1 || isempty(current)
    return;
end
if isfield(current, 'ratio'), rep.ratio = current.ratio; end
if isfield(current, 'ci'), rep.ci = current.ci; end
if isfield(current, 'coverage'), rep.coverage = current.coverage; end
if isfield(current, 'nRecords'), rep.nRecords = current.nRecords; end
if nargin < 2 || isempty(baseline) || ~isfield(baseline, 'ratio') || ~isfield(current, 'ratio')
    return;
end
if baseline.ratio > 0 && current.ratio > 0
    rep.logDistanceShrink = abs(log(baseline.ratio)) - abs(log(current.ratio));
end
if isfield(baseline, 'actualRms') && isfield(current, 'actualRms') && baseline.actualRms > 0
    rep.actualRmsWorsening = current.actualRms / baseline.actualRms - 1;
end
% Report-only checks (no gate action in Round 8).
ciInside = isfinite(rep.ci(1)) && isfinite(rep.ci(2)) && rep.ci(1) >= 0.8 && rep.ci(2) <= 1.25;
if ciInside && rep.logDistanceShrink >= 0.5 * abs(log(baseline.ratio)) && rep.actualRmsWorsening <= 0.05 && rep.nRecords >= 30
    rep.verdict = 'would-meet-upgrade-criteria-report-only';
else
    rep.verdict = 'does-not-meet-upgrade-criteria-report-only';
end
end
