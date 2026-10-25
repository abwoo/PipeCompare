function pipecompare_short_fault_probe(outDir)
% Development stress: amplitude, autocorrelation and window alignment.
mkdir(outDir);fs=100;rows={};
for seed=620001:620005
    rs=RandStream('mt19937ar','Seed',seed);
    for rho=[0 0.6]
        X=filter(1,[1 -rho],randn(rs,12,60*fs),[],2);X=10*X./std(X,0,2);
        for amplitude=[0 10 30 100]
            for duration=[0.25 0.5 1 2 5 10 20]
                for offset=[0 0.25]
                    at=round((20+offset)*fs)+(1:round(duration*fs));D=X;
                    D(3,at)=D(3,at)+amplitude*randn(rs,1,numel(at));
                    truth=false(size(X));if amplitude>0,truth(3,at)=true;end
                    for mode={'long','multiscale'}
                        [st,check]=pipecompare.run.Steps.badStretches(D,fs,mode{1});
                        pred=false(size(X));
                        for k=1:numel(st),pred(st(k).row,st(k).samples(1):st(k).samples(2))=true;end
                        recall=NaN;if any(truth(:)),recall=nnz(pred&truth)/nnz(truth);end
                        rows(end+1,:)={seed,rho,amplitude,duration,offset,mode{1},check.assessed,recall,nnz(pred&~truth)/nnz(~truth)}; %#ok<AGROW>
                    end
                end
            end
        end
    end
end
T=cell2table(rows,'VariableNames',{'seed','rho','amplitude','seconds','offset','mode','assessed','recall','outsideShare'});
writetable(T,fullfile(outDir,'probe.csv'));summary={};
for mode={'long','multiscale'}
    for amplitude=[0 10 30 100]
        for duration=[0.25 0.5 1 2 5 10 20]
            D=T(strcmp(T.mode,mode{1})&T.amplitude==amplitude&T.seconds==duration,:);
            summary(end+1,:)={mode{1},amplitude,duration,mean(D.recall,'omitnan'),mean(D.outsideShare),nnz(D.assessed),height(D)}; %#ok<AGROW>
        end
    end
end
S=cell2table(summary,'VariableNames',{'mode','amplitude','seconds','recall','outsideShare','assessed','cases'});
writetable(S,fullfile(outDir,'probe_summary.csv'));
fprintf('PIPECOMPARE_SHORT_PROBE_JSON %s\n',jsonencode(table2struct(S)));
end
