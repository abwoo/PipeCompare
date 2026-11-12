classdef IcaRankLock
    %ICARANKLOCK Original-pop_runica-compatible unshuffled fitting rank.
    % Development-only ICA-start sensitivity: compute once on the same
    % unpermuted ICA fitting copy, then pass an explicit PCA dimension for
    % each later sample-order start. The seed-zero fit remains unchanged.
    % The native EEGLAB algorithm's first-3000-sample rule is intentionally
    % retained to isolate the permutation confound, not to estimate the
    % universally optimal data rank.
    methods (Static)
        function [r, audit] = estimate(data, takeMin)
            if nargin < 2
                src=which('pop_runica');
                assert(~isempty(src),'PipeCompare:IcaRankLock','EEGLAB pop_runica unavailable');
                takeMin=contains(fileread(src), 'tmprank2 = min(');
            end
            assert(islogical(takeMin)&&isscalar(takeMin), ...
                'PipeCompare:IcaRankLock','Invalid EEGLAB rank rule.');
            assert(isnumeric(data)&&~isempty(data), ...
                'PipeCompare:IcaRankLock','ICA fitting data missing.');
            data=double(reshape(data,size(data,1),[]));
            assert(size(data,1)>=2 && size(data,2)>=2 && all(isfinite(data(:))), ...
                'PipeCompare:IcaRankLock','Invalid nonfinite/too-small ICA data.');
            sampleCount=min(3000,size(data,2));
            sample=data(:,1:sampleCount);
            r1=rank(sample);
            [~,D]=eig(cov(sample',1));
            r2=sum(diag(D)>1e-7);
            if r1~=r2
                if takeMin,r=min(r1,r2);else,r=max(r1,r2);end
            else
                r=r1;
            end
            assert(r>=1 && r<=size(data,1), ...
                'PipeCompare:IcaRankLock','Rank estimator returned invalid dimension.');
            audit=struct('rank',r,'numericRank',r1,'covarianceRank',r2, ...
                'sampleCount',sampleCount,'nChan',size(data,1), ...
                'eeGLABtakeMin',takeMin, ...
                'computedBeforeSeedPermutation',true, ...
                'productionSeedZeroUnchanged',true);
        end
    end
end
