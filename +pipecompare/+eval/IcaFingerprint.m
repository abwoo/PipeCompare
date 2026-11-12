classdef IcaFingerprint
    %ICAFINGERPRINT Read-only summary of an evaluated ICA pipeline.
    %
    % Experimental diagnostics only. Never changes the fitted EEG,
    % component rejection or Rank selection. The raw ICA scalp maps are
    % normalized by column norm for independent sign/permutation matching
    % in Python. Component numbers across starts are NOT comparable until
    % that matching is performed. Missing or ambiguous evidence abstains.
    methods (Static)
        function z = snapshot(c)
            z = struct('available',false,'reason','missing-ICA', ...
                'channels',[],'maps',[],'removedIndices',[], ...
                'removedCount',NaN,'icSteps',0,'icRemoveSteps',0, ...
                'scientificPromotionAuthorized',false);
            if ~isstruct(c) || ~isscalar(c) || ~isfield(c,'ica') ...
                    || ~iscell(c.ica) || isempty(c.ica)
                return;
            end
            z.icSteps = numel(c.ica);
            if z.icSteps~=1
                z.reason='multiple-ICA-fits-not-comparable';
                return;
            end
            d = c.ica{1};
            fields={'icaweights','icasphere','icachansind'};
            if ~isstruct(d) || ~all(isfield(d,fields))
                z.reason='incomplete-ICA-matrix';
                return;
            end
            W=double(d.icaweights)*double(d.icasphere);
            ch=double(d.icachansind(:)');
            if ~ismatrix(W) || isempty(W) || any(~isfinite(W(:))) ...
                    || numel(ch)~=size(W,2) || numel(unique(ch))~=numel(ch) ...
                    || any(~isfinite(ch)) || any(ch<1)
                z.reason='ICA-matrix-or-channel-invalid';
                return;
            end
            % For rank-deficient ICA, only compare the retained components.
            A=pinv(W);
            normA=sqrt(sum(A.^2,1));
            if any(~isfinite(A(:))) || any(~isfinite(normA)) || any(normA<eps)
                z.reason='ICA-topography-degenerate';
                return;
            end
            steps=c.steps;
            if ~iscell(steps)
                z.reason='step-audit-not-cell';
                return;
            end
            types=cellfun(@(s) s.type,steps,'UniformOutput',false);
            removeIdx=find(strcmp(types,'icremove'));
            z.icRemoveSteps=numel(removeIdx);
            if z.icRemoveSteps~=1
                z.reason='ambiguous-component-rejection-path';
                return;
            end
            s=steps{removeIdx};
            if ~isfield(s,'comps') || ~isfield(s,'icsRemoved') ...
                    || ~isscalar(s.icsRemoved) || ~isfinite(s.icsRemoved)
                z.reason='component-rejection-unrecorded';
                return;
            end
            removed=double(s.comps(:)');
            if any(~isfinite(removed)) || any(removed<1) || ...
                    any(removed~=floor(removed)) || ...
                    any(removed>size(A,2)) || ...
                    numel(unique(removed))~=numel(removed) || ...
                    numel(removed)~=s.icsRemoved
                z.reason='component-rejection-inconsistent';
                return;
            end
            if ~isfield(c,'icsRemoved') || ~isscalar(c.icsRemoved) ...
                    || ~isfinite(c.icsRemoved) || c.icsRemoved~=numel(removed)
                z.reason='component-count-disagrees-with-pipeline';
                return;
            end
            z.available=true;
            z.reason='complete-fingerprint';
            z.channels=ch;
            z.maps=A./normA;
            z.removedIndices=sort(removed);
            z.removedCount=numel(removed);
        end
    end
end
