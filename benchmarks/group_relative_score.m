function z = group_relative_score(objectives)
%GROUP_RELATIVE_SCORE Group-relative (GRPO-style) diagnostic scores.
%   STATUS: Round 8 bundle, PROTOTYPE / DIAGNOSTIC ONLY.
%   This is a pure arithmetic helper, not a ranking rule and not a
%   production pick rule. It standardises a vector of candidate
%   objectives by the group mean and standard deviation, the same
%   shape as GRPO's group-relative advantage (no critic, no learned
%   policy, no LLM, no external API, no user data sent anywhere).
%   A negative z means "better than the group mean" for a
%   lower-is-better objective. Callers must NOT use z as the final
%   pick: the final pick is still Rank.run on full-fidelity data.
%   Validation gate (see MANIFEST.md): report correlation of z with
%   full-fidelity Rank objectives on held-out seeds, plus survival
%   rate of the full-fidelity best under z-screening, before any
%   use in screening. Not wired into Rank.run in Round 8.
objectives = objectives(:);
finiteVals = objectives(isfinite(objectives));
if numel(finiteVals) < 2
    z = NaN(size(objectives));
    return;
end
mu = mean(finiteVals);
sd = std(finiteVals);
if sd == 0 || ~isfinite(sd)
    z = zeros(size(objectives));
    z(~isfinite(objectives)) = NaN;
    return;
end
z = (objectives - mu) / sd;
end
