function pipecompare_open_validation(dataDir,outDir,shard)
% Frozen public-data evaluation. Metrics are observations, not success gates.
% Only numerical/reference/replay invariants and execution errors fail CI.
% Real EEG background is NOT certified clean neural ground truth.
if ~isfolder(outDir), mkdir(outDir); end
copyfile(fullfile(dataDir,'provenance.json'),fullfile(outDir,'provenance.json'));
protocol=jsondecode(fileread(fullfile(fileparts(mfilename('fullpath')),'open_data_protocol.json')));
meta=struct('commit',getenv('GITHUB_SHA'),'run',getenv('GITHUB_RUN_ID'), ...
    'matlab',version,'eeglab',eeg_getversion,'protocol',protocol,'shard',shard);
writeJson(fullfile(outDir,'metadata.json'),meta);
failures={}; S=pipecompare_record_summary(); calibration=table();
allFiles=dir(fullfile(dataDir,'S*.mat'));
files=allFiles(~cellfun('isempty',regexp({allFiles.name},'^S\d\d\.mat$','once')));
expectedBlink=numel(protocol.blink.test_subjects(shard+1:6:end));
assert(numel(files)==expectedBlink,'PipeCompare:Validation','Expected %d blink records, found %d.',expectedBlink,numel(files));
for k=1:numel(files)
    try
        D=load(fullfile(files(k).folder,files(k).name));
        T=blinkMetrics(D,files(k).name,protocol.blink.matching_margin_seconds);
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
        EEG=backgroundEEG(D);
        try
            appendCsv(spectralMetrics(D,subject),fullfile(outDir,'spectral.csv'));
        catch ex,failures{end+1}=recordFailure(ex,[name ' spectra'],outDir);end %#ok<AGROW>
        try
            appendCsv(stretchMetrics(EEG,subject),fullfile(outDir,'stretches.csv'));
        catch ex,failures{end+1}=recordFailure(ex,[name ' stretches'],outDir);end %#ok<AGROW>
        [EEG,population,realized,prevalence]=injectedERP(EEG,subject);
        setBase(EEG); tic;
        transcript=evalc('[~,~,r]=pop_pipecompare(EEG,''measure'',''P3'',''events'',{''11'',''31''},''reference'',''average'',''show'',''off'');');
        seconds=toc; saveTranscript(outDir,name,transcript);
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
        EEG=loadPublicSet(fullfile(files(k).folder,files(k).name));
        T=realERP(EEG,files(k).name,outDir); appendCsv(T,fullfile(outDir,'real_erp.csv'));
    catch ex, failures{end+1}=recordFailure(ex,files(k).name,outDir); end %#ok<AGROW>
end
writeJson(fullfile(outDir,'execution.json'),struct('failures',{failures},'completed',isempty(failures)));
assert(isempty(failures),'PipeCompare:Validation','%d public-data execution/invariant failures; see artifacts.',numel(failures));
end

function T=blinkMetrics(D,name,margin)
X=D.data; valid=logical(D.valid); X(:,~valid)=NaN;
[hit,model]=pipecompare.eval.Hmm.events(X,D.fs);
% An interval counts once. Matching is one-to-one, so a long detection
% spanning several annotated blinks cannot receive several true positives.
d=diff([false hit false]); a=find(d==1); b=find(d==-1)-1;
gt=D.annotations(:,1)'; matched=false(size(gt)); tp=0; timing=[];
for j=1:numel(a)
    start=D.time(a(j))-margin; stop=D.time(b(j))+margin;
    candidates=find(~matched & gt>=start & gt<=stop);
    if isempty(candidates), continue; end
    centre=(D.time(a(j))+D.time(b(j)))/2;
    [~,at]=min(abs(gt(candidates)-centre)); at=candidates(at);
    matched(at)=true; tp=tp+1; timing(end+1)=abs(gt(at)-centre); %#ok<AGROW>
end
precision=tp/max(numel(a),1); recall=tp/max(numel(gt),1);
f1=2*precision*recall/max(precision+recall,eps);
minutes=nnz(valid)/D.fs/60; falseEvents=numel(a)-tp;
T=table({name},numel(gt),numel(a),tp,falseEvents,precision,recall,f1, ...
    median(timing),falseEvents/max(minutes,eps),mean(~valid),mean(hit(valid)),model.valid, ...
    D.effective_fs,D.timestamp_jitter_seconds, ...
    'VariableNames',{'recording','annotations','detections','truePositive','falsePositive','precision','recall','f1', ...
    'medianTimingErrorSeconds','falseEventsPerMinute','unavailableShare','activeShare','assessed','effectiveFs','timestampJitterSeconds'});
end

function EEG=backgroundEEG(D)
wanted={'Fp1','Fp2','F7','F3','Fz','F4','F8','C3','Cz','C4','P3','Pz','P4','O1','Oz','O2'};
[yes,loc]=ismember(lower(wanted),lower(D.labels)); assert(all(yes),'PipeCompare:Validation','Missing source channels.');
cap=readlocs(fullfile(fileparts(which('eeglab')),'functions','supportfiles','Standard-10-5-Cap385.sfp'));
[yes,at]=ismember(lower(wanted),lower({cap.labels})); assert(all(yes));
EEG=eeg_emptyset(); EEG.data=double(D.data(loc,:)); EEG.srate=D.fs;
EEG.nbchan=numel(wanted); EEG.pnts=size(EEG.data,2); EEG.trials=1;
EEG.xmin=0; EEG.xmax=(EEG.pnts-1)/EEG.srate; EEG.chanlocs=cap(at);
[EEG.chanlocs.type]=deal('EEG'); EEG.setname='EEGMMIDB public background';
EEG=eeg_checkset(EEG);
end

function T=spectralMetrics(D,subject)
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

function [EEG,population,realized,prevalence]=injectedERP(EEG,subject)
rs=RandStream('mt19937ar','Seed',90000+subject); fs=EEG.srate; n=EEG.pnts;
background=double(EEG.data); background=background-mean(background,2);
background=background/max(median(std(background,0,2)),eps)*8;
background=sign(randn(rs))*circshift(background,randi(rs,n),2);
labels={EEG.chanlocs.labels}; xyz=[[EEG.chanlocs.X]' [EEG.chanlocs.Y]' [EEG.chanlocs.Z]'];
xyz=xyz./vecnorm(xyz,2,2); at=find(strcmpi(labels,'Pz')); front=find(strcmpi(labels,'Fp1'));
topo=exp(-acos(max(-1,min(1,xyz*xyz(at,:)'))).^2/(2*0.5^2));
eyeTopo=exp(-acos(max(-1,min(1,xyz*xyz(front,:)'))).^2/(2*0.6^2));
lat=round(fs*(3:1.35:(n/fs-3))); nt=numel(lat); types=repmat({'31'},1,nt); types(1:2:end)={'11'};
types=types(randperm(rs,nt)); amp=8*exp(0.2*randn(rs));
levels=[0 0.3 0.7 0.95 1]; prevalence=levels(mod(subject-21,numel(levels))+1);
data=background; eog=0.5*background(front,:); neural=zeros(size(background));
tt=-round(0.2*fs):round(fs); window=(tt/fs>=0.3 & tt/fs<=0.6); base=(tt/fs>=-0.2 & tt/fs<=0);
trueTrial=zeros(1,nt); correctedTopo=topo-mean(topo); sigma=0.08;
for j=1:nt
    target=strcmp(types{j},'11'); meanAmp=amp*(1/3+2/3*target);
    hasBlink=rand(rs)<prevalence;
    % Genuine neural amplitude can covary with blinking. The known target
    % includes that association; it must not be called causal artifact bias.
    a=meanAmp*max(0,1+0.3*randn(rs))*(1+0.3*hasBlink);
    delay=0.4+0.04*randn(rs); wave=exp(-0.5*((tt/fs-delay)/sigma).^2);
    neural(:,lat(j)+tt)=neural(:,lat(j)+tt)+topo*(a*wave);
    trueTrial(j)=correctedTopo(at)*a*(mean(wave(window))-mean(wave(base)));
    if hasBlink
        blink=150*exp(-0.5*((tt/fs-0.42)/0.06).^2);
        data(:,lat(j)+tt)=data(:,lat(j)+tt)+eyeTopo*blink;
        eog(lat(j)+tt)=eog(lat(j)+tt)+blink;
    end
end
EEG.data=single([data+neural;eog]); EEG.nbchan=size(EEG.data,1);
EEG.chanlocs(end+1)=EEG.chanlocs(front); EEG.chanlocs(end).labels='EOG1';EEG.chanlocs(end).type='EOG';
EEG.event=struct('type',types,'latency',num2cell(lat),'duration',0);
EEG=eeg_checkset(EEG,'eventconsistency');EEG=eeg_checkset(EEG,'makeur');
EEG.history='';EEG.setname=sprintf('Independent_real_background_%d',subject);
realized=[mean(trueTrial(strcmp(types,'11'))) mean(trueTrial(strcmp(types,'31')))];
% E max(0,1+0.3 Z), with Gaussian latency variation integrated analytically.
z=1/0.3; ampFactor=0.5*erfc(-z/sqrt(2))+0.3*exp(-z^2/2)/sqrt(2*pi);
wide=sqrt(sigma^2+0.04^2); shape=sigma/wide*exp(-0.5*((tt/fs-0.4)/wide).^2);
p=correctedTopo(at)*amp*ampFactor*(1+0.3*prevalence)*(mean(shape(window))-mean(shape(base)));
population=[p p/3];
end

function T=realERP(EEG,name,outDir)
% The authors' event mapping, 26-ms display correction and correct response
% criterion (200-1000 ms). This is a six-person external execution/precision
% example, not a known-truth calibration or a population-effect replication.
code=arrayfun(@(e)str2double(string(e.type)),EEG.event);
lat=[EEG.event.latency]; stimulus=ismember(code,[11:15 21:25 31:35 41:45 51:55]);
response=find(ismember(code,[201 202]));
for j=find(stimulus)
    following=response(lat(response)>lat(j)); good=false;
    if ~isempty(following)
        q=following(1); seconds=(lat(q)-lat(j))/EEG.srate;
        good=code(q)==201 && seconds>0.2 && seconds<1;
    end
    if good
        type='pc_frequent';if ismember(code(j),[11 22 33 44 55]),type='pc_rare';end
        EEG.event(j).type=type;EEG.event(j).latency=lat(j)+0.026*EEG.srate;
    end
end
transcript=evalc('EEG=pop_resample(EEG,160);'); %#ok<NASGU>
% Keep the scalp channels with known standard positions, preserving source
% voltage units. No reference waveform from these people tunes the method.
cap=readlocs(fullfile(fileparts(which('eeglab')),'functions','supportfiles','Standard-10-5-Cap385.sfp'));
[yes,at]=ismember(lower({EEG.chanlocs.labels}),lower({cap.labels}));
EEG=pop_select(EEG,'channel',find(yes));EEG.chanlocs=cap(at(yes));[EEG.chanlocs.type]=deal('EEG');
EEG.urevent=[];EEG=eeg_checkset(EEG,'makeur');EEG.history='';setBase(EEG);
transcript=evalc('[~,~,r]=pop_pipecompare(EEG,''measure'',''P3'',''events'',{''pc_rare'',''pc_frequent''},''recipe'',''filters'',''reference'',''average'',''show'',''off'');');
saveTranscript(outDir,name,transcript);
rows={};
for k=find(strcmp(r.ranking.table.status,'feasible'))'
    o=r.cands(k).m.objectives(1); names=r.ranking.table.Properties.VariableNames;
    gains=names(startsWith(names,'gain_'));g=r.ranking.table.(gains{1})(k);
    rows(end+1,:)={name,k,r.labels{k},any(r.ranking.recommended==k),o.estimate(1)-o.estimate(2), ...
        sqrt(sum(o.sme.^2)),g,min(r.cands(k).m.retention)}; %#ok<AGROW>
end
T=cell2table(rows,'VariableNames',{'recording','candidate','steps','selected','rawRareMinusFrequentUv','independentTrialSE','gain','minRetention'});
end

function setBase(EEG)
assignin('base','EEG',EEG);assignin('base','ALLEEG',EEG);assignin('base','CURRENTSET',1);
end

function EEG=loadPublicSet(path)
% This distribution stores integer scalar metadata. EEGLAB's arithmetic
% expects doubles; normalize metadata on the in-memory copy only. Keep the
% original checksum-verified .set/.fdt files untouched and use native loading.
D=load(path,'-mat');EEG=D.EEG;
for field={'setname','filename','filepath','history','ref','comments'}
    if isfield(EEG,field{1}) && (ischar(EEG.(field{1})) || isstring(EEG.(field{1})))
        EEG.(field{1})=reshape(char(EEG.(field{1})),1,[]);
    end
end
for k=1:numel(EEG.chanlocs)
    EEG.chanlocs(k).labels=reshape(char(EEG.chanlocs(k).labels),1,[]);
end
for field={'srate','nbchan','pnts','trials','xmin','xmax','times'}
    if isfield(EEG,field{1}),EEG.(field{1})=double(EEG.(field{1}));end
end
for field={'latency','duration','urevent','epoch'}
    for k=1:numel(EEG.event)
        if isfield(EEG.event,field{1}) && isnumeric(EEG.event(k).(field{1}))
            EEG.event(k).(field{1})=double(EEG.event(k).(field{1}));
        end
    end
end
for k=1:numel(EEG.event)
    if ischar(EEG.event(k).type) || isstring(EEG.event(k).type)
        EEG.event(k).type=reshape(char(EEG.event(k).type),1,[]);
    end
end
[EEG.filepath,base,ext]=fileparts(path);EEG.filename=[base ext];
transcript=evalc('EEG=eeg_checkset(EEG,''loaddata'');'); %#ok<NASGU>
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
