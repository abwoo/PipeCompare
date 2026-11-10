function public_pilot(dataDir,outDir,subject)
% One frozen source record. Keep adverse results and production exceptions.
if ~isfolder(outDir),mkdir(outDir);end
p=jsondecode(fileread(fullfile(fileparts(mfilename('fullpath')),'public_pilot_protocol.json')));
assert(ismember(subject,[p.source.files.subject]),'PipeCompare:Pilot','Unplanned participant.');
assert(strcmp(version('-release'),'2026a')&&strcmp(eeg_getversion(),'2026.0.0'),'PipeCompare:Pilot','Wrong toolchain.');
source=jsondecode(fileread(fullfile(dataDir,'source.json')));
result=struct('subject',subject,'completed',false,'algorithmCommit',p.algorithm_commit, ...
 'executionCommit',getenv('GITHUB_SHA'),'workflowRun',getenv('GITHUB_RUN_ID'), ...
 'matlab',version,'eeglab',eeg_getversion(),'eeglabSha',getenv('EEGLAB_SHA'), ...
 'protocolSha256',source.protocol_sha256,'source',source,'errors',{{}}, ...
 'erp',emptyReadout(),'minimal',emptyReadout(),'common',emptyReadout(), ...
 'alpha',struct('completed',false,'recommended',false),'replayMaxError',NaN,'exportMaxError',NaN);
algorithmOverride=getenv('PILOT_ALGORITHM_COMMIT');
if ~isempty(algorithmOverride),result.algorithmCommit=algorithmOverride;end
writeJson(fullfile(outDir,'result.json'),result);
clock=tic;
try
 D=load(fullfile(dataDir,'input.mat')); EEG=background(D,p);
 result.sourceChannels={EEG.chanlocs.labels};
 try
  setBase(EEG);alphaClock=tic;
  transcript=evalc('[~,~,a]=pop_pipecompare(EEG,''measure'',''alpha'',''recipe'',''filters'',''reference'',''asis'',''show'',''off'');');
  saveText(fullfile(outDir,'alpha.log'),transcript);
  writetable(a.ranking.table,fullfile(outDir,'alpha_candidates.csv'));
  result.alpha=struct('completed',true,'recommended',~isempty(a.ranking.recommended), ...
   'candidates',numel(a.cands),'seconds',toc(alphaClock));
 catch ex,result.errors{end+1}=failure('unaltered-alpha',ex);end
 try
  T=faults(EEG,subject);writetable(T,fullfile(outDir,'faults.csv'));
 catch ex,result.errors{end+1}=failure('channel-faults',ex);end
 [E,population,realized,prevalence]=inject(EEG,subject,p);
 result.populationTruth=population;result.realizedTruth=realized;result.blinkPrevalence=prevalence;
 try
  setBase(E);fitClock=tic;
  transcript=evalc('[~,~,r]=pop_pipecompare(E,''measure'',''P3'',''events'',{''11'',''31''},''recipe'',''standard'',''reference'',''average'',''show'',''off'');');
  saveText(fullfile(outDir,'standard.log'),transcript);
  result.fitSeconds=toc(fitClock);result.candidates=numel(r.cands);result.notRun=r.notRun;
  assert(r.notRun==0,'PipeCompare:Pilot','Incomplete search.');
  result.feasible=sum(strcmp(r.ranking.table.status,'feasible'));
  result.searchOptions=r.options;result.searchOptions.progress=[];
  result.icaSeedCheck=r.icaSeeds;
  writetable(r.ranking.table,fullfile(outDir,'candidates.csv'));
  feasible=find(strcmp(r.ranking.table.status,'feasible'));
  if ~isempty(feasible)
   C=pipecompare.eval.Calibration.rows(r,feasible,population,realized,subject);
   writetable(C,fullfile(outDir,'calibration.csv'));
  end
  pick=r.ranking.recommended;
  if ~isempty(pick)
   result.erp=readout(r,pick,population,realized);
   losses=nan(1,numel(feasible));
   for q=1:numel(feasible),v=readout(r,feasible(q),population,realized);losses(q)=v.mse;end
   result.oracleFeasibleMSE=min(losses);result.selectionRegret=result.erp.mse-result.oracleFeasibleMSE;
   transcript=evalc('R=pipecompare.run.Executor.replay(r,pick,true);');
   saveText(fullfile(outDir,'replay.log'),transcript);
   m=pipecompare.eval.Measure.candidate(R,r.contract,r.ref,r.options);
   result.replayMaxError=max(abs(m.objectives(1).estimate-r.cands(pick).m.objectives(1).estimate));
   assert(result.replayMaxError<p.evaluation.replay_tolerance_uv,'PipeCompare:Pilot','Native replay estimate changed.');
   exportName=sprintf('pilot_export_%03d',subject);exportFile=fullfile(outDir,[exportName '.m']);
   pipecompare.PipeCompare.writeScript(r,pick,exportFile);
   priorPath=path;cleanPath=onCleanup(@()path(priorPath));addpath(outDir); %#ok<NASGU>
   f=str2func(exportName);transcript=evalc('X=f(E);');
   saveText(fullfile(outDir,'export.log'),transcript);
   mx=pipecompare.eval.Measure.candidate(X,r.contract,r.ref,r.options);
   result.exportMaxError=max(abs(mx.objectives(1).estimate-r.cands(pick).m.objectives(1).estimate));
   assert(result.exportMaxError<p.evaluation.replay_tolerance_uv,'PipeCompare:Pilot','Executable export changed the estimate.');
   clear cleanPath
  else
   result.noRecommendationReasons=r.ranking.table.reason;
  end
 catch ex,result.errors{end+1}=failure('standard-erp-and-replay',ex);end
 % Fixed comparators are declared before outcomes, never chosen by truth.
 for arm={'minimal','common'}
  try
   setBase(E);con=pipecompare.eval.Contract('conditions',{'target',{'11'};'standard',{'31'}}, ...
    'epoch',[-0.2 0.8],'baseline',[-0.2 0],'components',{'P3',[0.3 0.6],{'Pz'}});
   plan=pipecompare.plan.Plan();
   plan=plan.add('highpass','cutoff',0.1);plan=plan.add('lowpass','cutoff',30);
   if strcmp(arm{1},'common')
    plan=plan.add('badchannels','measure','kurt+prob','threshold',5,'action','interpolate','exclude',{{'EOG1'}});
   end
   plan=plan.add('reref','mode','average','exclude',{{'EOG1'}});
   if strcmp(arm{1},'common')
    plan=plan.add('ica');plan=plan.add('icremove','threshold',0.8,'classes',{{'Eye','Muscle'}});
   end
   plan=plan.add('epoch');plan=plan.add('baseline');
   plan=plan.add('reject_threshold','uv',100,'exclude',{{'EOG1'}});
   transcript=evalc('fixed=pipecompare.PipeCompare.optimize(plan,con,struct(''nSplits'',0));');
   saveText(fullfile(outDir,[arm{1} '.log']),transcript);
   if strcmp(fixed.cands(1).status,'ok')
    result.(arm{1})=readout(fixed,1,population,realized);
    result.(arm{1}).passedConstraints=strcmp(fixed.ranking.table.status{1},'feasible');
   end
  catch ex,result.errors{end+1}=failure([arm{1} '-fixed'],ex);end
 end
catch ex,result.errors{end+1}=failure('source-or-injection',ex);
end
result.seconds=toc(clock);result.completed=isempty(result.errors);
writeJson(fullfile(outDir,'result.json'),result);
fprintf('PIPECOMPARE_PILOT_RESULT %s\n',jsonencode(result));
% The collector checks all planned records and distinguishes execution from science.
end

function EEG=background(D,p)
wanted=cellstr(p.excerpt.channels);labels=cellstr(D.labels(:));
[yes,loc]=ismember(lower(wanted),lower(labels));assert(all(yes),'Missing source channels.');
cap=readlocs(fullfile(fileparts(which('eeglab')),'functions','supportfiles','Standard-10-5-Cap385.sfp'));
[yes,at]=ismember(lower(wanted),lower({cap.labels}));assert(all(yes),'Missing montage sites.');
EEG=eeg_emptyset();EEG.data=double(D.data(loc,:));EEG.srate=double(D.fs);
EEG.nbchan=numel(wanted);EEG.pnts=size(EEG.data,2);EEG.trials=1;
assert(EEG.pnts==p.excerpt.duration_seconds*EEG.srate,'Wrong excerpt length.');
EEG.xmin=0;EEG.xmax=(EEG.pnts-1)/EEG.srate;EEG.chanlocs=cap(at);
[EEG.chanlocs.type]=deal('EEG');EEG.setname='Public resting EEG excerpt';EEG=eeg_checkset(EEG);
end

function [EEG,population,realized,prevalence]=inject(EEG,subject,p)
s=p.simulation;rs=RandStream('mt19937ar','Seed',s.seed_base+subject);fs=EEG.srate;n=EEG.pnts;
bg=double(EEG.data);bg=bg-mean(bg,2);bg=bg/max(median(std(bg,0,2)),eps)*s.background_scale_uv;
bg=sign(randn(rs))*circshift(bg,randi(rs,n),2);
labels={EEG.chanlocs.labels};xyz=[[EEG.chanlocs.X]' [EEG.chanlocs.Y]' [EEG.chanlocs.Z]'];
xyz=xyz./vecnorm(xyz,2,2);at=find(strcmpi(labels,'Pz'));front=find(strcmpi(labels,'Fp1'));
topo=exp(-acos(max(-1,min(1,xyz*xyz(at,:)'))).^2/(2*0.5^2));
eyeTopo=exp(-acos(max(-1,min(1,xyz*xyz(front,:)'))).^2/(2*0.6^2));
lat=round(fs*(3:s.trial_spacing_seconds:(n/fs-3)));nt=numel(lat);
types=repmat({'31'},1,nt);types(1:2:end)={'11'};types=types(randperm(rs,nt));
amp=s.target_amplitude_uv*exp(0.2*randn(rs));
prevalence=s.blink_prevalence(mod(subject-1,numel(s.blink_prevalence))+1);
data=bg;neural=zeros(size(bg));eog=0.5*bg(front,:);
tt=-round(0.2*fs):round(fs);times=tt/fs;
win=times>=0.3 & times<=0.6;base=times>=-0.2 & times<=0;
canonicalTopo=topo-mean(topo);trueTrial=zeros(1,nt);
for j=1:nt
 target=strcmp(types{j},'11');mu=amp*(s.standard_amplitude_ratio+(1-s.standard_amplitude_ratio)*target);
 blink=rand(rs)<prevalence;
 a=mu*max(0,1+s.amplitude_cv*randn(rs))*(1+(s.neural_blink_multiplier-1)*blink);
 wave=exp(-0.5*((times-s.latency_mean_seconds-s.latency_sd_seconds*randn(rs))/s.width_seconds).^2);
 neural(:,lat(j)+tt)=neural(:,lat(j)+tt)+topo*(a*wave);
 trueTrial(j)=canonicalTopo(at)*a*(mean(wave(win))-mean(wave(base)));
 if blink
  b=s.blink_amplitude_uv*exp(-0.5*((times-0.42)/0.06).^2);
  data(:,lat(j)+tt)=data(:,lat(j)+tt)+eyeTopo*b;eog(lat(j)+tt)=eog(lat(j)+tt)+b;
 end
end
EEG.data=single([data+neural;eog]);EEG.nbchan=size(EEG.data,1);
EEG.chanlocs(end+1)=EEG.chanlocs(front);EEG.chanlocs(end).labels='EOG1';EEG.chanlocs(end).type='EOG';
EEG.event=struct('type',types,'latency',num2cell(lat),'duration',0);
EEG=eeg_checkset(EEG,'eventconsistency');EEG=eeg_checkset(EEG,'makeur');
EEG.history='';EEG.setname=sprintf('Public_background_injection_%03d',subject);
realized=[mean(trueTrial(strcmp(types,'11'))) mean(trueTrial(strcmp(types,'31')))];
z=1/s.amplitude_cv;amplitudeFactor=0.5*erfc(-z/sqrt(2))+s.amplitude_cv*exp(-z^2/2)/sqrt(2*pi);
wide=sqrt(s.width_seconds^2+s.latency_sd_seconds^2);
shape=s.width_seconds/wide*exp(-0.5*((times-s.latency_mean_seconds)/wide).^2);
a=canonicalTopo(at)*amp*amplitudeFactor*(1+(s.neural_blink_multiplier-1)*prevalence)*(mean(shape(win))-mean(shape(base)));
population=[a a*s.standard_amplitude_ratio];
end

function a=emptyReadout()
a=struct('available',false,'estimate',[],'correctedEstimate',[],'actualError',[], ...
 'rawError',[],'conditionalError',[],'reportedError',[],'mse',NaN,'rawMse',NaN, ...
 'conditionalMse',NaN,'reportedMse',NaN,'gain',NaN,'pipeline','');
end
function a=readout(r,k,population,realized)
a=emptyReadout();o=r.cands(k).m.objectives(1);
gcol=matlab.lang.makeValidName(['gain_' r.ref.objectives{1}]);
scol=matlab.lang.makeValidName(['shapeSpread_' r.ref.objectives{1}]);
g=r.ranking.table.(gcol)(k);sp=r.ranking.table.(scol)(k);
if ~(isfinite(g)&&g>0&&all(isfinite(o.estimate))),return;end
[reported,~,~]=pipecompare.eval.Rank.conditionErrors(o,g,sp);
a.available=true;a.estimate=o.estimate;a.correctedEstimate=o.estimate/g;a.gain=g;a.pipeline=r.labels{k};
a.actualError=a.correctedEstimate-population;a.rawError=o.estimate-population;
a.conditionalError=a.correctedEstimate-realized;a.reportedError=reported;
a.mse=mean(a.actualError.^2);a.rawMse=mean(a.rawError.^2);
a.conditionalMse=mean(a.conditionalError.^2);a.reportedMse=mean(reported.^2);
end

function T=faults(EEG,subject)
rs=RandStream('mt19937ar','Seed',961000+subject);clean=double(EEG.data);
clean=clean/max(median(std(clean,0,2)),eps)*10;labels={EEG.chanlocs.labels};pz=find(strcmpi(labels,'Pz'));rows={};
for seconds=[1 5 10 20]
 idx=20*EEG.srate+(1:round(seconds*EEG.srate));x=clean;x(pz,idx)=x(pz,idx)+100*randn(rs,1,numel(idx));
 [st,check]=pipecompare.run.Steps.badStretches(x,EEG.srate);pred=false(size(x));truth=false(size(x));truth(pz,idx)=true;
 for k=1:numel(st),pred(st(k).row,st(k).samples(1):st(k).samples(2))=true;end
 rows(end+1,:)={subject,seconds,check.assessed,check.reason,nnz(pred&truth)/nnz(truth),nnz(pred&~truth)/nnz(~truth)}; %#ok<AGROW>
end
T=cell2table(rows,'VariableNames',{'subject','faultSeconds','assessed','reason','sampleRecall','outsideInjectedShare'});
end
function setBase(EEG)
assignin('base','EEG',EEG);assignin('base','ALLEEG',EEG);assignin('base','CURRENTSET',1);
end
function s=failure(stage,ex)
s=struct('stage',stage,'identifier',ex.identifier,'message',ex.message,'trace',getReport(ex,'extended','hyperlinks','off'));
fprintf(2,'PILOT_FAILURE %s\n',jsonencode(s));
end
function writeJson(file,v)
fid=fopen([file '.tmp'],'w');assert(fid>=0,'Cannot write result.');c=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(v));clear c;movefile([file '.tmp'],file,'f');
end
function saveText(file,s)
fid=fopen(file,'w');assert(fid>=0);c=onCleanup(@()fclose(fid));fprintf(fid,'%s',s); %#ok<NASGU>
end
