classdef IcaRobust
    %ICAROBUST Scientific-gate-aware all-seed worst-objective SHADOW rule.
    % Research-only; never changes Rank.run, Executor's recommendation,
    % the frozen ERP CORE thresholds, or the 20-subject PARTIAL verdict.
    % A recipe must pass all original per-seed Rank statuses and the
    % unflagged effective n>=5 for BOTH conditions for all seeds.
    methods (Static)
        function out = preview(check)
            % CHECK is the unmodified result.icaSeeds / icaSeedCheck output.
            % In-memory MATLAB cell matrices retain their NxS layout.
            % Python JSON readers must independently handle column-major
            % flattening of those cell arrays.
            required = {'ids','objective','rankStatus','strata', ...
                        'candidateMeasurements','recommended'};
            assert(isstruct(check) && all(isfield(check,required)), ...
                'PipeCompare:IcaRobust','Incomplete ICA seed evidence.');
            ids = check.ids(:); obj = check.objective;
            n = numel(ids);
            assert(n>0 && numel(unique(ids))==n && size(obj,1)==n ...
                && size(obj,2)>=2, 'PipeCompare:IcaRobust', ...
                'ICA IDs/objectives not aligned.');
            seeds = size(obj,2);
            statuses = check.rankStatus; strata = check.strata;
            measures = check.candidateMeasurements;
            assert(iscell(statuses) && isequal(size(statuses),[n,seeds]) ...
                && iscell(strata) && isequal(size(strata),[n,seeds]) ...
                && iscell(measures) && isequal(size(measures),[n,seeds]), ...
                'PipeCompare:IcaRobust','Seed/candidate cell matrices not aligned.');
            baseline = check.recommended(1);
            k = find(ids==baseline,1);
            assert(~isempty(k) && isfinite(baseline), ...
                'PipeCompare:IcaRobust','Baseline recommendation is missing.');
            baseStratum = strata{k,1};
            feasible = false(n,1); reasons = cell(n,1);
            for i=1:n
                why={};
                if ~all(cellfun(@(x) isequal(x,baseStratum),strata(i,:)))
                    why{end+1}='not-baseline-stratum'; %#ok<AGROW>
                end
                if ~all(strcmp(statuses(i,:),'feasible'))
                    why{end+1}='not-feasible-every-seed'; %#ok<AGROW>
                end
                if ~all(isfinite(obj(i,:)))
                    why{end+1}='nonfinite-objective'; %#ok<AGROW>
                end
                for seed=1:seeds
                    z=measures{i,seed};
                    if ~isstruct(z) || ~isscalar(z) || ~isfield(z,'candidateId') ...
                            || z.candidateId~=ids(i) || ~isfield(z,'seed') ...
                            || z.seed~=seed-1
                        why{end+1}='measurement-misaligned'; %#ok<AGROW>
                        continue;
                    end
                    if ~isfield(z,'hasMeasurement') || ~z.hasMeasurement ...
                            || ~isfield(z,'correctedDifferenceEstimate') ...
                            || ~isscalar(z.correctedDifferenceEstimate) ...
                            || ~isfinite(z.correctedDifferenceEstimate) ...
                            || ~isfield(z,'gain') || ~isscalar(z.gain) ...
                            || ~isfinite(z.gain) || z.gain<=0
                        why{end+1}='finite-calibrated-P3-missing'; %#ok<AGROW>
                    end
                    if ~isfield(z,'effectiveUnflagged') ...
                            || numel(z.effectiveUnflagged)~=2 ...
                            || ~all(isfinite(z.effectiveUnflagged(:))) ...
                            || ~all(z.effectiveUnflagged(:) >= 5-1e-9)
                        why{end+1}='frozen-effective-n-below-five'; %#ok<AGROW>
                    end
                end
                feasible(i)=isempty(why);
                reasons{i}=unique(why,'stable');
            end
            costs=inf(n,1);
            costs(feasible)=max(obj(feasible,:),[],2);
            eligible=find(feasible);
            chosen=NaN; reason='no-common-feasible-candidate';
            if ~isempty(eligible)
                best=min(costs);
                tied=find(costs==best);
                if numel(tied)==1
                    chosen=ids(tied);
                    reason='unique-all-seed-feasible-minimax';
                else
                    reason='minimax-tie-abstain';
                end
            end
            out=struct('id',chosen,'reason',reason,'eligibleIds',ids(eligible)', ...
                'costs',costs,'reasons',{reasons}, ...
                'frozenTrialFloor',5,'candidateSelectedForProduction',false, ...
                'existingRankRecommendationUnchanged',true);
        end
    end
end
