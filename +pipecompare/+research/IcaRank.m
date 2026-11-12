classdef IcaRank
    %ICARANK Order-independent full-fit-data rank for ICA research only.
    % The default EEGLAB path estimates rank on the first 3000 samples.
    % Permuting the fit data therefore changes both ICA initialization
    % AND the estimated PCA dimension. This method freezes rank using
    % all fitted points before permutation, without changing any source
    % EEG samples or original production defaults.
    methods (Static)
        function [r,diagInfo] = wholeData(X)
            assert(isnumeric(X)&&ndims(X)<=3&&size(X,1)>=2, ...
                'PipeCompare:IcaRank','Expected channel x samples fitted EEG.');
            X=reshape(double(X),size(X,1),[]);
            assert(size(X,2)>=size(X,1) && all(isfinite(X(:))), ...
                'PipeCompare:IcaRank','Too few or nonfinite ICA fit samples.');
            % Reproduce pop_runica/getrank's *two* rank checks and exact
            % covariance eigenvalue threshold, but over the same complete
            % fitted observations for every seed (permutation invariant).
            rankSvd=rank(X);
            eigenvalues=eig(cov(X.',1));
            rankCov=sum(eigenvalues>1e-7);
            r=min(rankSvd,rankCov);
            assert(r>=2 && r<=size(X,1), 'PipeCompare:IcaRank', ...
                'Ill-defined full fit-data ICA rank.');
            diagInfo=struct('channelCount',size(X,1), ...
                'fitSampleCount',size(X,2),'svdRank',rankSvd, ...
                'covarianceEigenRank',rankCov,'chosenPcaDimension',r, ...
                'rankTolerance',1e-7, ...
                'policy','development-only-full-fit-data-pre-permutation');
        end
    end
end
