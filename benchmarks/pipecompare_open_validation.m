function pipecompare_open_validation(dataDir,outDir,shard)
% Frozen public-data evaluation. Metrics are observations, not success gates.
% Only numerical/reference/replay invariants and execution errors fail CI.
% Real EEG background is NOT certified clean neural ground truth.
if ~isfolder(outDir), mkdir(outDir); end
copyfile(fullfile(dataDir,'provenance.json'),fullfile(outDir,'provenance.json'));
protocol=jsondecode(fileread(fullfile(fileparts(mfilename('fullpath')),'open_data_protocol.json')));
meta=struct('commit',getenv('GITHUB_SHA'),'run',getenv('GITHUB_RUN_ID'), ...
    'matlab',version,'eeglab',eeg_getversion,'fitTimer','explicit-tic-handle','protocol',protocol,'shard',shard);
writeJson(fullfile(outDir,'metadata.json'),meta);
failures={}; S=pipecompare_record_summary(); calibration=table();
allFiles=dir(fullfile(dataDir,'S*.mat'));
files=allFiles(~cellfun('isempty',regexp({allFiles.name},'^S\d\d\.mat$','once')));
expectedBlink=numel(protocol.blink.test_subjects(shard+1:6:end));
assert(numel(files)==expectedBlink,'PipeCompare:Validation','Expected %d blink records, found %d.',expectedBlink,numel(files));
for k=1:numel(files)
    try
        D=load(fullfile(files(k).folder,files(k).name));
        T=pipecompare_open_blinks(D,files(k).name,protocol.blink.matching_margin_seconds);
        appendCsv(T,fullfile(outDir,'blink.csv'));
    catch ex, failures{end+1}=recordFailure(ex,files(k).name,outDir); end %#ok<AGROW>
end
files=allFiles(~cellfun('isempty',regexp({allFiles.name},'^S\d\d\dR01\.mat$','once')));
assert(numel(files)==5,'PipeCompare:Validation','Expected five EEGMMIDB records, found %d.',numel(files));
for k=1:numel(files)
    name=files(k).name; subject=sscanf(name,'S%dR01.mat');
    fprintf('[open validation] start %s, source participant %d\n',name,subject);
    try
        D=load(fullfile(files(k).folder,name));
        EEG=pipecompare_open_background(D);
        try
            appendCsv(pipecompare_open_spectra(D,subject),fullfile(outDir,'spectral.csv'));
        catch ex,failures{end+1}=recordFailure(ex,[name ' spectra'],outDir);end %#ok<AGROW>
        try
            appendCsv(stretchMetrics(EEG,subject),fullfile(outDir,'stretches.csv'));
        catch ex,failures{end+1}=recordFailure(ex,[name ' stretches'],outDir);end %#ok<AGROW>
        [EEG,population,realized,prevalence]=pipecompare_open_injection(EEG,subject);
        setBase(EEG); fitClock=tic;
        transcript=evalc('[~,~,r]=pop_pipecompare(EEG,''measure'',''P3'',''events'',{''11'',''31''},''reference'',''average'',''show'',''off'');');
        seconds=toc(fitClock); saveTranscript(outDir,name,transcript);
        feasible=find(strcmp(r.ranking.table.status,'feasible'))';
        C=pipecompare.eval.Calibration.rows(r,feasible,population,realized,subject);
        if ~isempty(C), if isempty(calibration), calibration=C; else, calibration=[calibration;C]; end, end %#ok<AGROW>
        writetable(calibration,fullfile(outDir,'calibration.csv'));
        item=pipecompare_record_summary(r,name);
        save(fullfile(outDir,sprintf('score_%03d.mat',subject)),'item','-v7');
        S(end+1)=item; %#ok<AGROW>
        save(fullfile(outDir,'derived_scores.mat'),'S','-v7');
        replayError=NaN; pick=r.ranking.recommended;
        if ~isempty(pick)
            transcript=evalc('E=pipecompare.run.Executor.replay(r,pick(1),true);'); %#ok<NASGU>
            m=pipecompare.eval.Measure.candidate(E,r.contract,r.ref,r.options);
            replayError=max(abs(m.objectives(1).estimate-r.cands(pick(1)).m.objectives(1).estimate));
            assert(replayError<1e-6,'PipeCompare:Validation','Native replay changed the scored estimate.');
        end
        T=table(subject,prevalence,numel(r.cands),numel(feasible),~isempty(pick),seconds,replayError, ...
            'VariableNames',{'subject','blinkPrevalence','candidates','feasible','recommended','seconds','replayError'});
        appendCsv(T,fullfile(outDir,'recordings.csv'));
        fprintf('[open validation] finished %s: %d feasible, recommendation %s, %.1f s\n',name,numel(feasible),mat2str(pick),seconds);
    catch ex, failures{end+1}=recordFailure(ex,name,outDir); end %#ok<AGROW>
end
files=dir(fullfile(dataDir,'erp_core','*.set'));
for k=1:numel(files)
    try
        EEG=pipecompare_load_public_set(fullfile(files(k).folder,files(k).name));
        T=pipecompare_open_real_erp(EEG,files(k).name,outDir); appendCsv(T,fullfile(outDir,'real_erp.csv'));
    catch ex, failures{end+1}=recordFailure(ex,files(k).name,outDir); end %#ok<AGROW>
end
writeJson(fullfile(outDir,'execution.json'),struct('failures',{failures},'completed',isempty(failures)));
assert(isempty(failures),'PipeCompare:Validation','%d public-data execution/invariant failures; see artifacts.',numel(failures));
end




function T=stretchMetrics(EEG,subject)
rs=RandStream('mt19937ar','Seed',80000+subject); rows={};
clean=double(EEG.data); scale=median(std(clean,0,2)); clean=clean/max(scale,eps)*10;
labels={EEG.chanlocs.labels}; pz=find(strcmpi(labels,'Pz'));
for duration=[1 5 10 20]
    for adjacent=[false true]
        bad=pz; if adjacent, bad=[pz find(strcmpi(labels,'P3'))]; end
        interval=20*EEG.srate+(1:round(duration*EEG.srate));
        corrupted=clean; corrupted(bad,interval)=corrupted(bad,interval)+100*randn(rs,numel(bad),numel(interval));
        [st,check]=pipecompare.run.Steps.badStretches(corrupted,EEG.srate);
        pred=false(size(clean));
        for j=1:numel(st), pred(st(j).row,st(j).samples(1):st(j).samples(2))=true; end
        truth=false(size(clean)); truth(bad,interval)=true;
        recall=nnz(pred & truth)/nnz(truth); fp=nnz(pred & ~truth)/nnz(~truth);
        E=EEG; E.data=corrupted;
        oracle=struct('channel',labels(bad),'samples',repmat({[interval(1) interval(end)]},1,numel(bad)));
        repaired=pipecompare.run.Steps.interpStretches(E,oracle);
        before=sqrt(mean((corrupted(truth)-clean(truth)).^2));
        after=sqrt(mean((double(repaired.data(truth))-clean(truth)).^2));
        endpoints=double(repaired.data(bad,[interval(1)-1 interval(1) interval(end) interval(end)+1]));
        jumps=diff(endpoints,1,2);boundary=max(abs(jumps(:,[1 3])),[],'all');
        unchanged=max(abs(double(repaired.data(~truth))-corrupted(~truth)));
        assert(unchanged==0,'PipeCompare:Validation','Interpolation changed samples outside its declared intervals.');
        % Test additivity on the actual montage/operator, including the
        % adjacent-channel donor exclusion, rather than on saved scores.
        A=EEG;A.data=clean;B=EEG;B.data=corrupted-clean;
        A=pipecompare.run.Steps.interpStretches(A,oracle);B=pipecompare.run.Steps.interpStretches(B,oracle);
        replay=max(abs(double(repaired.data)-double(A.data)-double(B.data)),[],'all');
        assert(replay<1e-8,'PipeCompare:Validation','Interpolation operator is not replayable.');
        erpBias=mean(double(repaired.data(pz,interval))-clean(pz,interval));
        spectralChange=pipecompare.eval.Measure.bandScores(double(repaired.data(pz,interval)),EEG.srate,[8 12])- ...
            pipecompare.eval.Measure.bandScores(clean(pz,interval),EEG.srate,[8 12]);
        rows(end+1,:)={subject,duration,adjacent,check.assessed,check.reason,recall,fp,before,after,boundary,erpBias,spectralChange,replay}; %#ok<AGROW>
    end
end
T=cell2table(rows,'VariableNames',{'subject','faultSeconds','adjacent','assessed','reason','sampleRecall','falsePositiveShare','corruptedRMSE','oracleRepairRMSE','boundaryJumpUv','meanBiasUv','alphaChangeLog10PSD','replayMaxError'});
end



function setBase(EEG)
assignin('base','EEG',EEG);assignin('base','ALLEEG',EEG);assignin('base','CURRENTSET',1);
end


function appendCsv(T,path)
if isempty(T),return;end
if isfile(path), writetable(T,path,'WriteMode','append','WriteVariableNames',false);
else, writetable(T,path);end
end

function writeJson(path,value)
fid=fopen(path,'w');assert(fid>=0);done=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(value));
end

function message=recordFailure(ex,name,outDir)
message=sprintf('%s: %s',name,getReport(ex,'extended','hyperlinks','off'));
fprintf(2,'[open validation] %s\n',message);
fid=fopen(fullfile(outDir,'failures.txt'),'a');fprintf(fid,'%s\n',message);fclose(fid);
end

function saveTranscript(outDir,name,transcript)
fid=fopen(fullfile(outDir,[name '.log']),'w');fprintf(fid,'%s',transcript);fclose(fid);
end
