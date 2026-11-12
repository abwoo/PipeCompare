classdef IcaRobust
%ICAROBUST Exploratory common-feasible ICA-start robustness, not deployment.
% Use original Rank.run statuses for EVERY ICA start, not only finite scores.
% If no candidate is feasible under every start, abstain. Never retroactively
% change the original acceptance constraints, scientific gate, or PR51 ranking.
methods (Static)
    function out = select(ids, obj, status, stratum)
        if nargin < 4, stratum = repmat({'default'},size(status)); end
        if ~isnumeric(ids) || numel(unique(ids))~=numel(ids) ...
                || any(~isfinite(ids)) || any(ids<1) || any(ids~=floor(ids))
            error('PipeCompare:IcaRobust','Candidate IDs must be distinct positive integers');
        end
        ids=ids(:);
        if ~isnumeric(obj) || ndims(obj)>2 || size(obj,1)~=numel(ids) ...
                || size(obj,2)<2 || ~iscell(status) || ~isequal(size(status),size(obj)) ...
                || ~iscell(stratum) || ~isequal(size(stratum),size(obj))
            error('PipeCompare:IcaRobust','Seed objective/status/stratum matrices are misaligned');
        end
        out=struct('status','abstain','candidateId',NaN,'worstObjective',NaN, ...
            'eligibleIds',zeros(1,0),'eligibleWorst',zeros(1,0), ...
            'nEligible',0,'reason','no-common-feasible-candidate', ...
            'allSeedsScientificallyFeasible',false, ...
            'originalRecommendationUnchanged',true);
        valid=false(numel(ids),1);
        for i=1:numel(ids)
            row=status(i,:);
            % A finite objective on a rejected or failed candidate MUST NOT pass.
            isFeasible=all(cellfun(@(s) (ischar(s)||isstring(s)) ...
                && strcmp(char(s),'feasible'),row));
            consistentStratum=all(cellfun(@(s) (ischar(s)||isstring(s)) ...
                && strcmp(char(s),char(stratum{i,1})),stratum(i,:)));
            valid(i)=isFeasible && consistentStratum ...
                && all(isfinite(obj(i,:))) && all(obj(i,:)>=0);
        end
        feas=find(valid);
        if isempty(feas),return;end
        groups=unique(stratum(feas,1),'stable');
        if numel(groups)~=1
            out.reason='multiple-feasible-strata-cannot-be-ranked-together';return;
        end
        w=max(obj(feas,:),[],2);
        out.eligibleIds=ids(feas)';
        out.eligibleWorst=w';
        out.nEligible=numel(feas);
        best=min(w);
        ties=find(abs(w-best) <= 1e-12*max(1,abs(best)));
        if numel(ties)~=1
            out.reason='indistinguishable-minimax-tie-requires-predeclared-tiebreak';return;
        end
        out.status='selected';
        out.reason='all-seed-Rank-run-feasible-and-finite';
        out.candidateId=ids(feas(ties));
        out.worstObjective=best;
        out.allSeedsScientificallyFeasible=true;
    end
end
end
