classdef IcaRobustFamily
%ICAROBUSTFAMILY Exploratory ICA minimax within one P300 estimator family.
% Family is defined solely by whether weight_trials is used. The original
% selected recipe's family is the reference; all other science constraints
% are delegated to IcaRobust.select, preserving abstention on failures.
% NO modifications to the original ERP CORE ranking or validation verdict.
methods (Static)
    function out=select(ids,obj,status,stratum,supportPass,familyBySeed,baselineId)
        if nargin~=7
            error('PipeCompare:IcaRobustFamily','Seven science-frozen inputs required');
        end
        ids=ids(:);
        if ~iscell(familyBySeed) || ~isequal(size(familyBySeed),size(obj)) ...
                || numel(ids)~=size(obj,1)
            error('PipeCompare:IcaRobustFamily','Family rows and seed columns must match science audit');
        end
        if numel(unique(ids))~=numel(ids)
            error('PipeCompare:IcaRobustFamily','Duplicate candidate IDs');
        end
        allowedLabels={'weighted','unweighted','unavailable'};
        for i=1:numel(familyBySeed)
            if ~(ischar(familyBySeed{i}) || isstring(familyBySeed{i})) ...
                    || ~any(strcmp(char(familyBySeed{i}),allowedLabels))
                error('PipeCompare:IcaRobustFamily','Invalid trial-weight estimator provenance');
            end
        end
        base=find(ids==baselineId);
        out=struct('status','abstain','candidateId',NaN, ...
            'baselineId',baselineId,'baselineFamily','unavailable', ...
            'includedCandidateIds',zeros(1,0), ...
            'excludedCandidateIds',ids', ...
            'reason','baseline-candidate-not-present', ...
            'originalRecommendationUnchanged',true, ...
            'estimatorFamilyLocked',true);
        if isempty(base),return;end
        ref=char(familyBySeed{base,1});
        out.baselineFamily=ref;
        if strcmp(ref,'unavailable')
            out.reason='baseline-estimator-unavailable';return;
        end
        allowed=false(numel(ids),1);
        for i=1:numel(ids)
            allowed(i)=all(strcmp(familyBySeed(i,:),ref));
        end
        out.includedCandidateIds=ids(allowed)';
        out.excludedCandidateIds=ids(~allowed)';
        if ~any(allowed)
            out.reason='no-consistent-estimator-family';return;
        end
        robust=pipecompare.eval.IcaRobust.select(ids(allowed),obj(allowed,:), ...
            status(allowed,:),stratum(allowed,:),supportPass(allowed,:));
        out.status=robust.status;
        out.candidateId=robust.candidateId;
        out.reason=robust.reason;
        out.scientific=robust;
    end
end
end
