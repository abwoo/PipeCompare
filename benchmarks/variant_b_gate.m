function [allowed, reason] = variant_b_gate(validation)
%VARIANT_B_GATE Gate for a variant B arm: closed until validated.
%
% Report-only gate logic: it never runs a search, never refits ICA,
% and never creates an arm. It only decides whether a recorded
% validation is complete enough to ALLOW a named variant B arm to
% be built (in a later round, with held-out pairing). A passing
% gate is not a variant B result.
%
% Why the gate is closed by default: pop_pipecompare has NO
% protectReferenceCleanup option (an unknown option errors), and
% the first Rank.run defaults the guard to true, with adaptive
% leaves starting from that ranking's recommendation. Running
% pop_pipecompare three times would therefore give three identical
% arms with guarded-minus-unguarded identically zero: false
% evidence. A replay on an already-fitted family is variant A,
% never variant B, and must never be labelled variant B.
if nargin < 1 || isempty(validation) || ~isstruct(validation)
    allowed = false;
    reason = 'Variant B not validated: no validation record. A variant B arm must not be built.';
    return;
end
req = {'seed', 'guardOnPick', 'guardOffPick', 'method', 'planReuseVerified', 'explanation'};
for k = 1:numel(req)
    if ~isfield(validation, req{k}) || isempty(validation.(req{k}))
        allowed = false;
        reason = sprintf('Variant B not validated: validation record is missing field %s.', req{k});
        return;
    end
end
if strcmp(validation.method, 'three-runs') || strcmp(validation.method, 'variant-A-replay')
    allowed = false;
    reason = 'Variant B not validated: method is a known false-evidence path (three identical runs, or a variant A replay relabelled as variant B). Use the bypass-pop_pipecompare-optimize path (M3) or the production pass-through path (M1) with a 1-seed guard on/off validation.';
    return;
end
if ~ismember(validation.method, {'bypass-pop_pipecompare-optimize', 'production-passthrough'})
    allowed = false;
    reason = 'Variant B not validated: method must be bypass-pop_pipecompare-optimize (M3) or production-passthrough (M1).';
    return;
end
if ~isequal(validation.planReuseVerified, true)
    allowed = false;
    reason = 'Variant B not validated: plan construction reuse is not verified in MATLAB yet.';
    return;
end
if isequal(validation.guardOnPick, validation.guardOffPick) && (~ischar(validation.explanation) || isempty(validation.explanation))
    allowed = false;
    reason = 'Variant B not validated: guard on/off picks are identical, so an explanation is required (why the guard did not change this seed).';
    return;
end
allowed = true;
reason = 'Variant B gate passed for a NAMED variant B arm only: the arm still requires held-out pairing, its production Standard/Balanced definitions are unchanged, and its result must never be merged into a variant A mean.';
end
