function T=pipecompare_spectral_validation(outDir)
% Independent known-spectrum stress tests, separate from human records.
% Finite circular stationary Gaussian white/1-f backgrounds and on-grid
% deterministic lines. Target: unsmoothed, bin-averaged one-sided PSD.
% Log/Jensen and taper-smoothing bias are reported, not tuned away.
if ~isfolder(outDir),mkdir(outDir);end
rs=RandStream('mt19937ar','Seed',120001);fs=160;B=600;rows={};
for family={'white','pink'}
    for seconds=[2 4 8]
        n=round(seconds*fs);hz=(0:n-1)'*fs/n;folded=min(hz,fs-hz);
        shape=ones(n,1);
        if strcmp(family{1},'pink'),shape=1./sqrt(max(folded,1));end
        for band={[8 12],[10 10.5]}
            f=band{1};sel=hz>=f(1)&hz<=f(2);
            density=2/fs*shape.^2;df=fs/n;
            meanNoise=mean(density(sel));
            positions={'none','nearest-grid center','edge','outside'};
            frequencies=[NaN round(mean(f)/df)*df f(2) f(2)+max(0.5,2*df)];
            for p=1:numel(positions)
                W=fft(randn(rs,n,B),[],1).*shape;
                X=real(ifft(W,[],1));
                target=meanNoise;frequency=frequencies(p);
                if p>1
                    amplitude=sqrt(2*meanNoise*nnz(sel)*df);
                    phase=2*pi*rand(rs,1,B);
                    X=X+amplitude*sin(2*pi*(0:n-1)'/fs*frequency+phase);
                    if frequency>=f(1)&&frequency<=f(2)
                        target=target+amplitude^2/(2*nnz(sel)*df);
                    end
                end
                [scores,info]=pipecompare.eval.Measure.bandScores(reshape(X,1,n,B),fs,f);
                truth=log10(target);err=scores-truth;
                assert(all(isfinite(err)),'PipeCompare:Validation','Nonfinite spectral stress result.');
                rows(end+1,:)={family{1},seconds,f(1),f(2),positions{p},frequency,B, ...
                    info.resolved,info.NW,info.tapers,truth,mean(err),sqrt(mean(err.^2)),std(err)/sqrt(B)}; %#ok<AGROW>
            end
        end
    end
end
T=cell2table(rows,'VariableNames',{'background','seconds','lowHz','highHz','peakPosition','peakHz', ...
    'independentSegments','resolved','NW','tapers','trueLog10BinMeanPSD','biasLog10PSD','rmseLog10PSD','monteCarloSEBias'});
writetable(T,fullfile(outDir,'known_spectrum.csv'));
fprintf('PIPECOMPARE_SPECTRAL_JSON %s\n',jsonencode(table2struct(T)));
end
