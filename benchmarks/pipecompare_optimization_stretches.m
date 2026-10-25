function T=pipecompare_optimization_stretches(EEG,subject)
% Real background: markings outside injections can include natural artifacts.
% Evaluate actual detector-driven interpolation separately from the oracle.
rs=RandStream('mt19937ar','Seed',180000+subject);rows={};
clean=double(EEG.data);clean=clean/max(median(std(clean,0,2)),eps)*10;
labels={EEG.chanlocs.labels};pz=find(strcmpi(labels,'Pz'));
for duration=[0 0.5 1 2 5 10 20]
    for adjacent=[false true]
        bad=pz;if adjacent,bad=[pz find(strcmpi(labels,'P3'))];end
        interval=20*EEG.srate+(1:round(duration*EEG.srate));
        corrupted=clean;corrupted(bad,interval)=corrupted(bad,interval)+100*randn(rs,numel(bad),numel(interval));
        truth=false(size(clean));truth(bad,interval)=true;
        for mode={'long','multiscale'}
            [st,check]=pipecompare.run.Steps.badStretches(corrupted,EEG.srate,mode{1});
            pred=false(size(clean));
            for k=1:numel(st),pred(st(k).row,st(k).samples(1):st(k).samples(2))=true;end
            repair=struct('channel',{},'samples',{});
            for k=1:numel(st),repair(k)=struct('channel',labels{st(k).row},'samples',st(k).samples);end
            E=EEG;E.data=corrupted;E=pipecompare.run.Steps.interpStretches(E,repair);
            changed=double(E.data)-clean;
            recall=NaN;faultRMSE=NaN;
            if any(truth(:)),recall=nnz(pred&truth)/nnz(truth);faultRMSE=sqrt(mean(changed(truth).^2));end
            outside=nnz(pred&~truth)/nnz(~truth);
            untouched=max(abs(double(E.data(~pred))-corrupted(~pred)));
            if isempty(untouched),untouched=0;end
            assert(untouched==0,'PipeCompare:Validation','Repair changed samples outside detected intervals.');
            rows(end+1,:)={subject,duration,adjacent,mode{1},check.assessed,check.reason, ...
                recall,outside,faultRMSE,sqrt(mean(changed(:).^2)),sqrt(mean(changed(~truth).^2)),untouched}; %#ok<AGROW>
        end
    end
end
T=cell2table(rows,'VariableNames',{'subject','faultSeconds','adjacent','mode','assessed','reason', ...
    'sampleRecall','outsideInjectedShare','detectorRepairFaultRMSE','detectorRepairWholeRMSE','outsideInjectedRepairRMSE','untouchedMaxError'});
end
