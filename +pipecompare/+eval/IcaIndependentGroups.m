classdef IcaIndependentGroups
%ICAINDEPENDENTGROUPS Two disjoint triplet ICA recommendation replications.
% Seed groups {0,1,2} and {3,4,5}. All scientific and support gates are
% evaluated independently; no seed appears in both triplets. The original
% Standard ranking remains unchanged. This is research, NOT production.
methods (Static)
    function out=compare(ids,obj,status,stratum,supportPass,family,recommended,measurements)
        if nargin~=8 || size(obj,2)~=6 || numel(recommended)~=6 ...
                || numel(ids)~=size(obj,1) || ~isequal(size(status),size(obj)) ...
                || ~isequal(size(stratum),size(obj)) || ~isequal(size(supportPass),size(obj)) ...
                || ~isequal(size(family),size(obj)) || ~iscell(measurements) ...
                || ~isequal(size(measurements),size(obj))
            error('PipeCompare:IcaIndependentGroups','Exact six-start aligned seed audit required');
        end
        out=struct('status','not-comparable','sameCandidate',false, ...
            'groups',[],'sourceSeedLabels',0:5,'originalSearchRankingUnchanged',true, ...
            'isIndependentHeldoutStudy',false);
        descriptions=struct([]);
        for group=1:2
            cols=(1:3)+3*(group-1);
            anchor=recommended(cols(1));
            z=pipecompare.eval.IcaRobustFamily.select(ids,obj(:,cols), ...
                status(:,cols),stratum(:,cols),supportPass(:,cols), ...
                family(:,cols),anchor);
            z.seedNumbers=cols-1;
            z.independentGroup=group;
            if strcmp(z.status,'selected')
                k=find(ids==z.candidateId,1);
                if isempty(k),error('PipeCompare:IcaIndependentGroups','Selected candidate absent from seed matrix');end
                z.seedMeasurements=measurements(k,cols);
                for sd=1:3
                    item=z.seedMeasurements{sd};
                    if ~isfield(item,'candidateId') || item.candidateId~=z.candidateId ...
                            || ~isfield(item,'seed') || item.seed~=cols(sd)-1
                        error('PipeCompare:IcaIndependentGroups','Candidate/seed measurement provenance mismatch');
                    end
                end
            else
                z.seedMeasurements={};
            end
            if isempty(descriptions),descriptions=z;else,descriptions(end+1)=z;end %#ok<AGROW>
        end
        out.groups=descriptions;
        if all(strcmp({descriptions.status},'selected'))
            out.sameCandidate=descriptions(1).candidateId==descriptions(2).candidateId;
            if out.sameCandidate,out.status='same-across-disjoint-triplets';
            else,out.status='different-across-disjoint-triplets';end
        else
            out.status='at-least-one-group-abstained';
        end
    end
end
end
