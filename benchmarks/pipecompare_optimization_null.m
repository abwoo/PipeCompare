function pipecompare_optimization_null(outDir,stage)
% Independent simulated negative controls; never counted as human records.
if ~isfolder(outDir),mkdir(outDir);end
count=20;first=520001;if strcmp(stage,'confirmation'),count=100;first=530001;end
rows={};fs=100;
for seed=first:first+count-1
    rs=RandStream('mt19937ar','Seed',seed);
    X=filter(1,[1 -0.6],randn(rs,12,60*fs),[],2);
    for mode={'long','multiscale'}
        [st,check]=pipecompare.run.Steps.badStretches(X,fs,mode{1});pred=false(size(X));
        for k=1:numel(st),pred(st(k).row,st(k).samples(1):st(k).samples(2))=true;end
        rows(end+1,:)={seed,mode{1},check.assessed,mean(pred(:)),numel(st)}; %#ok<AGROW>
    end
end
T=cell2table(rows,'VariableNames',{'seed','mode','assessed','falsePositiveShare','events'});
writetable(T,fullfile(outDir,'null_controls.csv'));
for mode={'long','multiscale'}
    D=T(strcmp(T.mode,mode{1}),:);
    fprintf('PIPECOMPARE_NULL mode=%s records=%d assessed=%d meanFP=%.8g maxFP=%.8g\n', ...
        mode{1},height(D),nnz(D.assessed),mean(D.falsePositiveShare),max(D.falsePositiveShare));
    assert(all(D.assessed),'PipeCompare:Validation','Null inputs unexpectedly unassessed.');
    assert(mean(D.falsePositiveShare)<0.005&&max(D.falsePositiveShare)<0.05, ...
        'PipeCompare:Validation','Prespecified null marking limit exceeded.');
end
end
