function [rawChoice,groupChoice,eligible]=validation_group_choices(y,v,theta)
% Compare the same assessed candidates; no choice is returned for NaN rows.
assert(isvector(y)&&isequal(size(y),size(v),size(theta)), ...
    'PipeCompare:Validation','Supply one recording with aligned scores and variances.');
rawChoice=NaN;groupChoice=NaN;
eligible=find(isfinite(y)&y>0&isfinite(v)&v>=0&isfinite(theta)&theta>0);
if isempty(eligible),return;end
[~,a]=min(y(eligible));rawChoice=eligible(a);
[~,b]=min(theta(eligible));groupChoice=eligible(b);
end
