function design = variant_b_prototype()
%VARIANT_B_PROTOTYPE Variant B design record, NOT a search.
%   STATUS: Round 8 bundle, DESIGN ONLY. No search is run here.
%   DECISION (recorded, not MATLAB-verified): M3 first.
%     M3 = bypass pop_pipecompare and call PipeCompare.optimize with
%       the guard in runOpts, re-fitting each arm separately.
%       Reason: it does not change production API, does not change
%       the production default (guard on), and keeps the guard switch
%       as a run-option contrast. Risk: plan-construction reuse from
%       pop_pipecompare has not been MATLAB-verified, so M3 needs a
%       minimal MATLAB validation (one seed, guard on vs off, picks
%       differ or are explained) before any held-out arm is built on it.
%     M1 = add an optional guard pass-through to pop_pipecompare.
%       Reason to defer: it is a production API change, needs its own
%       tests and its own commit, and would mix a production change
%       with a validation arm in the same bundle.
%   Variant B arms (when built): guarded refit, unguarded refit, each
%   with its own fitted family (unlike variant A, which shares one
%   family). Standard and balancedGuarded definitions are NOT changed
%   by variant B; variant B is a new contrast, not a replacement.
%   Gates: paired MSE vs Standard, 95% CI, calibration, pick
%   stability, and a negative control by selectedKey (not family index).
design = struct('decision', 'M3-first', 'deferred', 'M1-production-passthrough', ...
    'matlabValidationRequired', true, 'wiredIntoShadow', false, ...
    'arms', {{'guardedRefit', 'unguardedRefit'}}, ...
    'mustNotReplace', {{'standard', 'balancedGuarded'}});
end
