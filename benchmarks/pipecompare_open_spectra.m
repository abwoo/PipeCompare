function T=pipecompare_open_spectra(D,subject)
[yes,pz]=ismember('pz',lower(D.labels)); assert(yes);
x=D.data(pz,:); rows={};
for duration=[2 4 8]
    n=round(duration*D.fs); nt=floor(numel(x)/n);
    X=reshape(x(1:n*nt),1,n,nt);
    for band={[8 12],[4 7],[10 10.5]}
        f=band{1}; [actual,info]=pipecompare.eval.Measure.bandScores(X,D.fs,f);
        [H,e]=dpss(n,info.NW,info.tapers);
        if info.tapers==1
            % pmtm requires >=2 supplied tapers when DropLastTaper=false.
            % With one taper the matched independent reference is exactly
            % its windowed periodogram, using MATLAB's PSD normalization.
            [P,hz]=periodogram(reshape(X,n,nt),H,n,D.fs);
        else
            [P,hz]=pmtm(reshape(X,n,nt),H,e,n,D.fs,'unity','DropLastTaper',false);
        end
        reference=log10(mean(P(hz>=f(1)&hz<=f(2),:),1))';
        discrepancy=max(abs(actual-reference));
        assert(discrepancy<1e-6,'PipeCompare:Validation','Multitaper reference mismatch: %.9g',discrepancy);
        for freq=[10 12 13]
            injected=X+reshape(5*sin(2*pi*freq*(0:n-1)/D.fs),1,n,1);
            delta=mean(pipecompare.eval.Measure.bandScores(injected,D.fs,f)-actual);
            rows(end+1,:)={subject,duration,f(1),f(2),freq,info.NW,info.tapers,info.resolved,discrepancy,delta}; %#ok<AGROW>
        end
    end
end
T=cell2table(rows,'VariableNames',{'subject','seconds','lowHz','highHz','injectedHz','NW','tapers','resolved','referenceMaxError','injectedChangeLog10PSD'});
end
