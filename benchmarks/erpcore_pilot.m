function erpcore_pilot(dataDir,outDir,subject)
% One frozen ERP CORE P3 record. Keep adverse results and production exceptions.
% Real P3 data: no injected truth, no simulation. The primary endpoint is the
% bootstrap SE of the P3 difference amplitude (target minus standard), computed
% by the collector from the saved per-trial scores; this script records the
% point estimates, reported errors and the trial scores that make that
% recomputation possible. Run with no arguments to exercise the support
% recorder without EEG, fitting or performance outcomes.
if nargin==0,checkSupportRecorder();return;end
if ~isfolder(outDir),mkdir(outDir);end
p=jsondecode(fileread(fullfile(fileparts(mfilename('fullpath')),'erpcore_protocol.json')));
assert(ismember(subject,[p.source.files.subject]),'PipeCompare:Erpcore','Unplanned participant.');
assert(strcmp(version('-release'),'2026a')&&strcmp(eeg_getversion(),'2026.0.0'),'PipeCompare:Erpcore','Wrong toolchain.');
source=jsondecode(fileread(fullfile(dataDir,'source.json')));
assert(source.subject==subject,'PipeCompare:Erpcore','Source subject mismatch.');
% The protocol JSON cannot pin its own hash (a file cannot contain its own
% SHA256), so hash the frozen file at runtime and compare against the digest
% recorded by fetch at download time.
protocolPath=fullfile(fileparts(mfilename('fullpath')),'erpcore_protocol.json');
[hashStatus,hashOut]=system(sprintf('sha256sum "%s"',protocolPath));
assert(hashStatus==0,'PipeCompare:Erpcore','Protocol checksum unavailable.');
protocolSha256=strtok(hashOut);
assert(strcmp(source.protocol_sha256,protocolSha256),'PipeCompare:Erpcore','Protocol mismatch.');
result=struct('subject',subject,'completed',false,'algorithmCommit',p.algorithm_commit, ...
 'executionCommit',getenv('GITHUB_SHA'),'workflowRun',getenv('GITHUB_RUN_ID'), ...
 'runAttempt',getenv('GITHUB_RUN_ATTEMPT'),'matlab',version,'eeglab',eeg_getversion(), ...
 'eeglabSha',getenv('EEGLAB_SHA'),'measureSha256',getenv('PILOT_MEASURE_SHA256'), ...
 'protocolSha256',source.protocol_sha256,'source',source, ...
 'errors',{{}},'erp',emptyReadout(),'minimal',emptyReadout(),'common',emptyReadout(), ...
 'reference',struct('available',false),'hmmVeog',struct('available',false), ...
 'alpha',struct('completed',false,'recommended',false),'replayMaxError',NaN,'exportMaxError',NaN);
writeJson(fullfile(outDir,'result.json'),result);
clock=tic;
try
 D=load(fullfile(dataDir,'input.mat')); EEG=background(D,p);
 result.sourceChannels={EEG.chanlocs.labels};
 assert(EEG.pnts==p.excerpt.duration_seconds*EEG.srate,'Wrong excerpt length.');
 try
  T=faults(EEG,subject);writetable(T,fullfile(outDir,'faults.csv'));
 catch ex,result.errors{end+1}=failure('channel-faults',ex);end
 try
  result.reference=referenceAmplitude(D,p);
  writetable(struct2table(result.reference.counts),fullfile(outDir,'reference_counts.csv'));
 catch ex,result.errors{end+1}=failure('reference-amplitude',ex);end
 try
  setBase(EEG);fitClock=tic;
  transcript=evalc('[~,~,r]=pop_pipecompare(EEG,''measure'',''P3'',''events'',{''11'',''31''},''recipe'',''standard'',''reference'',''average'',''show'',''off'');');
  saveText(fullfile(outDir,'standard.log'),transcript);
  result.fitSeconds=toc(fitClock);result.candidates=numel(r.cands);result.notRun=r.notRun;
  assert(r.notRun==0,'PipeCompare:Erpcore','Incomplete search.');
  result.feasible=sum(strcmp(r.ranking.table.status,'feasible'));
  result.searchOptions=r.options;result.searchOptions.progress=[];
  result.icaSeedCheck=r.icaSeeds;
  writetable(r.ranking.table,fullfile(outDir,'candidates.csv'));
  feasible=find(strcmp(r.ranking.table.status,'feasible'));
  if ~isempty(feasible)
   audit=struct('cands',r.cands(feasible),'ref',r.ref,'options',r.ranking.options, ...
    'objective',r.ranking.objective,'ids',feasible,'ranking',r.ranking); %#ok<NASGU>
   save(fullfile(outDir,'scoring_audit.mat'),'audit','-v7');
  end
  pick=r.ranking.recommended;
  if ~isempty(pick)
   result.erp=readout(r,pick);
   result.selectedSupport=selectedSupport(r,pick);
   try
    result.hmmVeog=hmmVeogDiagnostic(D,result.selectedSupport,p);
   catch ex,result.errors{end+1}=failure('hmm-veog-diagnostic',ex);end
   transcript=evalc('R=pipecompare.run.Executor.replay(r,pick,true);');
   saveText(fullfile(outDir,'replay.log'),transcript);
   m=pipecompare.eval.Measure.candidate(R,r.contract,r.ref,r.options);
   result.replayMaxError=max(abs(m.objectives(1).estimate-r.cands(pick).m.objectives(1).estimate));
   assert(result.replayMaxError<p.evaluation.replay_tolerance_uv,'PipeCompare:Erpcore','Native replay estimate changed.');
   exportName=sprintf('erpcore_export_%03d',subject);exportFile=fullfile(outDir,[exportName '.m']);
   pipecompare.PipeCompare.writeScript(r,pick,exportFile);
   priorPath=path;cleanPath=onCleanup(@()path(priorPath));addpath(outDir); %#ok<NASGU>
   f=str2func(exportName);transcript=evalc('X=f(E);');
   saveText(fullfile(outDir,'export.log'),transcript);
   mx=pipecompare.eval.Measure.candidate(X,r.contract,r.ref,r.options);
   result.exportMaxError=max(abs(mx.objectives(1).estimate-r.cands(pick).m.objectives(1).estimate));
   assert(result.exportMaxError<p.evaluation.replay_tolerance_uv,'PipeCompare:Erpcore','Executable export changed the estimate.');
   clear cleanPath
  else
   result.noRecommendationReasons=r.ranking.table.reason;
  end
 catch ex,result.errors{end+1}=failure('standard-erp-and-replay',ex);end
 % Fixed comparators are declared before outcomes, never chosen by truth.
 for arm={'minimal','common'}
  try
   setBase(EEG);con=pipecompare.eval.Contract('conditions',{'target',{'11'};'standard',{'31'}}, ...
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
    result.(arm{1})=readout(fixed,1);
    result.(arm{1}).passedConstraints=strcmp(fixed.ranking.table.status{1},'feasible');
    result.(arm{1}).support=selectedSupport(fixed,1);
   end
  catch ex,result.errors{end+1}=failure([arm{1} '-fixed'],ex);end
 end
catch ex,result.errors{end+1}=failure('source',ex);
end
result.seconds=toc(clock);result.completed=isempty(result.errors);
writeJson(fullfile(outDir,'result.json'),result);
fprintf('PIPECOMPARE_ERPCORE_RESULT %s\n',jsonencode(result));
% The collector checks all planned records and distinguishes execution from science.
end

function EEG=background(D,p)
% 16 PR51 channels from input.mat, plus real VEOG as EOG1 and the relabeled events.
wanted=cellstr(p.excerpt.channels);labels=cellstr(D.labels(:));
[yes,loc]=ismember(lower(wanted),lower(labels));assert(all(yes),'Missing source channels.');
cap=readlocs(fullfile(fileparts(which('eeglab')),'functions','supportfiles','Standard-10-5-Cap385.sfp'));
[yes,at]=ismember(lower(wanted),lower({cap.labels}));assert(all(yes),'Missing montage sites.');
EEG=eeg_emptyset();EEG.data=double(D.data(loc,:));EEG.srate=double(D.fs);
EEG.nbchan=numel(wanted);EEG.pnts=size(EEG.data,2);EEG.trials=1;
EEG.xmin=0;EEG.xmax=(EEG.pnts-1)/EEG.srate;EEG.chanlocs=cap(at);
[EEG.chanlocs.type]=deal('EEG');
% Append real VEOG as EOG1 (the pilot appended a synthetic EOG1 here).
EEG.data(end+1,:)=double(D.veog);EEG.nbchan=EEG.nbchan+1;
EEG.chanlocs(end+1)=EEG.chanlocs(find(strcmpi(wanted,'Fp1'),1));
EEG.chanlocs(end).labels='EOG1';EEG.chanlocs(end).type='EOG';
EEG.event=struct('type',cellstr(D.eventTypes(:)),'latency',num2cell(double(D.eventLatency(:))),'duration',0);
EEG=eeg_checkset(EEG);EEG=eeg_checkset(EEG,'eventconsistency');EEG=eeg_checkset(EEG,'makeur');
EEG.history='';EEG.setname='ERP_CORE_P3_excerpt';
end

function T=faults(EEG,subject)
% Same known-duration channel-fault check as the PR51 pilot: a +100 uV noise
% burst on Pz, scored by the pipeline's own bad-stretch detector.
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

function ref=referenceAmplitude(D,p)
% Secondary reference only: the published ERP CORE pipeline's P3 difference
% amplitude (Pz, 300-600 ms, P9/P10 average reference, correct trials only),
% computed directly from the excerpt without the tool. The excerpt carries
% the recording's large DC offsets and slow drift (sub-001 Pz std ~696 uV
% before filtering), so the paper's 0.1 Hz high-pass is applied first via
% EEGLAB's pop_eegfiltnew; without it the fixed +-100 uV rejection below
% would reject nearly every trial. Simplifications vs the paper are
% documented in the frozen protocol: fixed +-100 uV amplitude rejection
% instead of the paper's individualized thresholds and ICA.
fs=double(D.fs);data=double(D.data);veog=double(D.veog);heog=double(D.heog);
tmp=eeg_emptyset();tmp.data=[data;veog;heog;double(D.refP9);double(D.refP10)];
tmp.srate=fs;tmp.nbchan=size(tmp.data,1);tmp.pnts=size(tmp.data,2);tmp.trials=1;
tmp=eeg_checkset(tmp);tmp=pop_eegfiltnew(tmp,0.1,[]);
data=tmp.data(1:16,:);veog=tmp.data(17,:);heog=tmp.data(18,:);
reref=data-(tmp.data(19,:)+tmp.data(20,:))/2;
pz=find(strcmpi(cellstr(D.labels(:)),'Pz'));
epochSamp=-round(0.2*fs):round(0.8*fs);baseSamp=epochSamp(epochSamp<0);
winSamp=epochSamp(epochSamp>=round(0.3*fs)&epochSamp<=round(0.6*fs));
types=cellstr(D.eventTypes(:));lat=double(D.eventLatency(:));
counts=struct('targetKept',0,'standardKept',0,'targetRejectedAmp',0,'standardRejectedAmp',0);
ampT=[];ampS=[];
for j=1:numel(lat)
 cols=lat(j)+epochSamp;if cols(1)<1||cols(end)>size(reref,2),continue;end
 ep=reref(:,cols);ep=ep-mean(ep(:,ismember(epochSamp,baseSamp)),2);
 vep=veog(cols);hep=heog(cols);
 rejected=any(abs(ep(:))>100)||any(abs(vep)>100)||any(abs(hep)>100);
 isTarget=strcmp(types{j},'11');
 if rejected
  if isTarget,counts.targetRejectedAmp=counts.targetRejectedAmp+1;else,counts.standardRejectedAmp=counts.standardRejectedAmp+1;end
  continue;
 end
 a=mean(ep(pz,ismember(epochSamp,winSamp)));
 if isTarget,ampT(end+1)=a;counts.targetKept=counts.targetKept+1; %#ok<AGROW>
 else,ampS(end+1)=a;counts.standardKept=counts.standardKept+1;end %#ok<AGROW>
end
ref=struct('available',true,'targetAmplitude',mean(ampT),'standardAmplitude',mean(ampS), ...
 'differenceAmplitude',mean(ampT)-mean(ampS),'counts',counts);
end

function d=hmmVeogDiagnostic(D,support,p)
% Secondary diagnostic only: do the HMM blink flags from the selected
% pipeline line up with large VEOG deflections in the same epochs?
fs=double(D.fs);veog=double(D.veog);lat=double(D.eventLatency(:));
epochSamp=-round(0.2*fs):round(0.8*fs);
d=struct('available',true,'veogThresholdUv',50,'nFlagged',0,'nUnflagged',0, ...
 'flaggedWithVeog',0,'unflaggedWithVeog',0,'flaggedMedianMaxVeog',NaN,'unflaggedMedianMaxVeog',NaN);
flaggedMax=[];unflaggedMax=[];
for c=1:numel(support)
 if strcmp(support(c).condition,'target'),rows=find(strcmp(cellstr(D.eventTypes(:)),'11'));
 else,rows=find(strcmp(cellstr(D.eventTypes(:)),'31'));end
 for k=1:numel(rows)
  j=rows(k);cols=lat(j)+epochSamp;if cols(1)<1||cols(end)>numel(veog),continue;end
  % Baseline-subtract within the epoch first: the raw VEOG carries the
  % recording's DC offset and slow drift, which would swamp the 50 uV
  % threshold. This mirrors how the tool's blink reference treats EOG.
  seg=veog(cols);seg=seg-mean(seg(epochSamp<0));
  mv=max(abs(seg));
  pos=find(support(c).referenceIds==j,1);
  if isempty(pos),continue;end
  if support(c).flagged(pos),flaggedMax(end+1)=mv; %#ok<AGROW>
  else,unflaggedMax(end+1)=mv;end %#ok<AGROW>
 end
end
d.nFlagged=numel(flaggedMax);d.nUnflagged=numel(unflaggedMax);
d.flaggedWithVeog=sum(flaggedMax>d.veogThresholdUv);
d.unflaggedWithVeog=sum(unflaggedMax>d.veogThresholdUv);
d.flaggedMedianMaxVeog=median(flaggedMax);d.unflaggedMedianMaxVeog=median(unflaggedMax);
end

function a=emptyReadout()
a=struct('available',false,'estimate',[],'correctedEstimate',[],'reportedError',[], ...
 'sme',[],'gain',NaN,'pipeline','','differenceAmplitude',NaN, ...
 'differenceSEAnalytic',NaN,'reportedDifferenceSE',NaN,'kept',[],'passedConstraints',false);
end

function a=readout(r,k)
% No known truth in real data: record the per-condition estimates, their
% SMEs and reported errors, and the difference amplitude with its analytic
% and reported SEs. The primary bootstrap SE is recomputed by the collector
% from the saved per-trial scores (selectedSupport), never from these SMEs.
a=emptyReadout();o=r.cands(k).m.objectives(1);
gcol=matlab.lang.makeValidName(['gain_' r.ref.objectives{1}]);
scol=matlab.lang.makeValidName(['shapeSpread_' r.ref.objectives{1}]);
g=r.ranking.table.(gcol)(k);sp=r.ranking.table.(scol)(k);
if ~(isfinite(g)&&g>0&&all(isfinite(o.estimate))),return;end
[reported,~,~]=pipecompare.eval.Rank.conditionErrors(o,g,sp);
a.available=true;a.estimate=o.estimate;a.correctedEstimate=o.estimate/g;a.gain=g;a.pipeline=r.labels{k};
a.sme=o.sme;a.reportedError=reported;a.kept=r.cands(k).m.kept;
a.differenceAmplitude=a.correctedEstimate(1)-a.correctedEstimate(2);
a.differenceSEAnalytic=sqrt(sum(o.sme.^2))/g;
a.reportedDifferenceSE=sqrt(sum(reported.^2));
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

function rows=selectedSupport(r,k)
% Observation only: do not alter reference labels, scoring, ranking or weights.
% Copied from the held-out pilot (same semantics, same checks).
m=r.cands(k).m;o=m.objectives(1);rows=struct([]);
assert(~isempty(r.ref.blink)&&~isempty(r.ref.blink{1}),'Missing artifact reference labels.');
for c=1:numel(r.ref.ids)
 x=o.X{c}(:,1);flags=logical(r.ref.blink{1}{c}(:));w=ones(size(x));
 if ~isempty(o.w),w=o.w{c}(:);end
 good=isfinite(x)&isfinite(w)&w>0;w(~good)=0;
 wn=w(~flags);wb=w(flags);effective=@(v) sum(v)^2/max(sum(v.^2),eps);
 a=struct('condition',r.ref.names{c},'referenceN',r.ref.n(c), ...
  'kept',m.kept(c),'retention',m.retention(c),'positiveN',sum(good), ...
  'flaggedN',sum(good&flags),'unflaggedN',sum(good&~flags), ...
  'effectiveN',effective(w),'effectiveFlaggedN',effective(wb),'effectiveUnflaggedN',effective(wn), ...
  'eyeBiasAssessed',logical(o.eyeBiasAssessed(c)), ...
  'eyeBiasFinite',isfinite(o.eyeBias(c))&&isfinite(o.eyeBiasSq(c)), ...
  'referenceIds',r.ref.ids{c},'scores',x','weights',w','flagged',flags');
 if c==1,rows=repmat(a,1,numel(r.ref.ids));else,rows(c)=a;end
end
end

function checkSupportRecorder()
% Exercise the actual recorder without EEG, fitting or performance outcomes.
r.ref.ids={[1 2 3],[4 5 6]};r.ref.names={'target','standard'};
r.ref.n=[3 3];r.ref.blink={{logical([0;1;0]),logical([1;0;0])}};
o.X={[2;4;NaN],[3;5;7]};o.w={};o.eyeBiasAssessed=[true true];
o.eyeBias=[0 0];o.eyeBiasSq=[0 0];
r.cands.m=struct('objectives',o,'kept',[2 3],'retention',[2/3 1]);
a=selectedSupport(r,1);
assert(numel(a)==2&&strcmp(a(1).condition,'target')&&strcmp(a(2).condition,'standard'));
assert(isequal([a.positiveN],[2 3])&&isequal([a.unflaggedN],[1 2]));
assert(isequal(a(1).weights,[1 1 0])&&isequal(a(2).referenceIds,[4 5 6]));
r.cands.m.objectives.w={[.5;1;NaN],[0;1;2]};a=selectedSupport(r,1);
assert(isequal([a.positiveN],[2 2])&&isequal([a.flaggedN],[1 0]));
assert(abs(a(2).effectiveUnflaggedN-1.8)<1e-12);
decoded=jsondecode(jsonencode(a));assert(numel(decoded)==2&&decoded(2).kept==3);
fprintf('HELDOUT_SUPPORT_RECORDER_CHECK_PASSED\n');
end
