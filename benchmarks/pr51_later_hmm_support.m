function pr51_later_hmm_support(baseDir,fullDir,outDir,subject)
% Development-only original-vs-full-task HMM source reference comparison.
% No ICA, no threshold optimization, no ranking, no new subjects, no merge.
if nargin==0
    assert(isequal(expectedDuration(4),544) && isequal(expectedDuration(8),460));
    fprintf('PR51_LATER_HMM_RECORDER_SMOKE_OK\n');
    return;
end
assert(nargin==4 && any(subject==[4,8]), 'PR51:Support','Unknown support subject.');
if ~isfolder(outDir),mkdir(outDir);end
p=jsondecode(fileread(fullfile(fileparts(mfilename('fullpath')),'erpcore_protocol.json')));
assert(p.excerpt.duration_seconds==300,'PR51:Support','Original protocol modified.');
assert(strcmp(version('-release'),'2026a') && strcmp(eeg_getversion(),'2026.0.0'), ...
    'PR51:Support','Pinned EEG toolchain mismatch.');
a=jsondecode(fileread(fullfile(baseDir,'source.json')));
b=jsondecode(fileread(fullfile(fullDir,'source.json')));
assert(a.subject==subject && b.subject==subject && ...
       a.duration_seconds==300 && b.duration_seconds==expectedDuration(subject) && ...
       ~a.exploratory_duration_override && b.exploratory_duration_override, ...
       'PR51:Support','Input length/provenance mismatch.');
assert(strcmp(a.source_set_sha256,b.source_set_sha256) && ...
       strcmp(a.source_fdt_sha256,b.source_fdt_sha256) && ...
       strcmp(a.source_events_sha256,b.source_events_sha256) && ...
       strcmp(a.protocol_sha256,b.protocol_sha256), ...
       'PR51:Support','Different raw source used between 300s and full data.');

con=pipecompare.eval.Contract('conditions',{'target',{'11'};'standard',{'31'}}, ...
    'epoch',[-0.2 0.8],'baseline',[-0.2 0], ...
    'components',{'P3',[0.3 0.6],{'Pz'}});
d0=load(fullfile(baseDir,'input.mat'));
d1=load(fullfile(fullDir,'input.mat'));
x=referenceTrialSupport(d0,p,con);
y=referenceTrialSupport(d1,p,con);
% Source trial counts are independently pinned by the task-event inventory.
assert(isequal([x.referenceN],[a.n_valid_targets,a.n_valid_standards]), ...
    'PR51:Support','Original source/event alignment mismatch.');
assert(isequal([y.referenceN],[b.n_valid_targets,b.n_valid_standards]), ...
    'PR51:Support','Extended source/event alignment mismatch.');
known=struct('subject',subject,'target',0,'standard',0);
if subject==4,known.target=4;known.standard=3;else,known.target=5;known.standard=20;end
if ~isequal([x.unflaggedN],[known.target known.standard])
    error('PR51:SupportBaseline', ...
        '300s HMM reference mismatch subject %d; observed %s expected [%d %d].', ...
        subject,mat2str([x.unflaggedN]),known.target,known.standard);
end

out=struct('subject',subject,'sourceSetSha256',a.source_set_sha256, ...
    'sourceFdtSha256',a.source_fdt_sha256,'sourceEventsSha256',a.source_events_sha256, ...
    'originalProtocolSha256',a.protocol_sha256, ...
    'originalSeconds',300,'fullSeconds',b.duration_seconds, ...
    'baseline',x,'extended',y, ...
    'addedValidTargetEvents',b.n_valid_targets-a.n_valid_targets, ...
    'addedValidStandardEvents',b.n_valid_standards-a.n_valid_standards, ...
    'fullWindowUnflaggedFloorMet',all([y.unflaggedN]>=5), ...
    'originalUnflaggedFloorMet',all([x.unflaggedN]>=5), ...
    'sameHMMScienceAndTrialRules',true,'notAnAlgorithmAccuracyTrial',true, ...
    'notAnIndependentParticipant',true,'productionPromotionAuthorized',false, ...
    'original20SubjectVerdict','PARTIAL');
f=fullfile(outDir,'later_hmm_support.json');
fid=fopen([f '.tmp'],'w');
assert(fid>=0,'PR51:Support','Could not save result');
fprintf(fid,'%s\n',jsonencode(out));fclose(fid);
movefile([f '.tmp'],f,'f');
fprintf('PR51_LATER_HMM_SCIENCE subject=%d original=[%d %d] extended=[%d %d] full=%d floor=%d\n', ...
    subject,x(1).unflaggedN,x(2).unflaggedN,y(1).unflaggedN,y(2).unflaggedN, ...
    b.duration_seconds,out.fullWindowUnflaggedFloorMet);
end

function duration=expectedDuration(subject)
if subject==4,duration=544;elseif subject==8,duration=460;
else,error('PR51:Support','Unplanned participant');end
end

function rows=referenceTrialSupport(D,p,con)
% Use exactly the source dataset reconstruction in erpcore_pilot.m; no
% preprocessing, rereference, ICA, stimulus reinterpretation or thresholds.
EEG=background(D,p);
ref=pipecompare.eval.Measure.reference(EEG,con);
assert(numel(ref.ids)==2 && numel(ref.n)==2 && numel(ref.blink)==1 && ...
    numel(ref.blink{1})==2, 'PR51:Support','Reference HMM evidence unavailable.');
rows=repmat(struct('condition','','referenceN',0,'flaggedN',0, ...
    'unflaggedN',0,'unflaggedIds',[],'flaggedIds',[]),1,2);
for k=1:2
    flag=ref.blink{1}{k};
    flag=logical(flag(:)');
    ids=ref.ids{k}(:)';
    assert(numel(flag)==numel(ids),'PR51:Support','HMM flags misaligned.');
    rows(k).condition=ref.names{k};
    rows(k).referenceN=ref.n(k);
    rows(k).flaggedN=nnz(flag);
    rows(k).unflaggedN=nnz(~flag);
    rows(k).unflaggedIds=ids(~flag);
    rows(k).flaggedIds=ids(flag);
end
end

function EEG=background(D,p)
wanted=cellstr(p.excerpt.channels);
labels=cellstr(D.labels(:));
[yes,loc]=ismember(lower(wanted),lower(labels));
assert(all(yes),'PR51:Support','Source missing EEG channels.');
cap=readlocs(fullfile(fileparts(which('eeglab')),'functions','supportfiles','Standard-10-5-Cap385.sfp'));
[yes,at]=ismember(lower(wanted),lower({cap.labels}));
assert(all(yes),'PR51:Support','Montage missing.');
EEG=eeg_emptyset();
EEG.data=double(D.data(loc,:));
EEG.srate=double(D.fs);
EEG.nbchan=numel(wanted);
EEG.pnts=size(EEG.data,2);
EEG.trials=1;
EEG.xmin=0;
EEG.xmax=(EEG.pnts-1)/EEG.srate;
EEG.chanlocs=cap(at);
[EEG.chanlocs.type]=deal('EEG');
EEG.data(end+1,:)=double(D.veog);
EEG.nbchan=EEG.nbchan+1;
EEG.chanlocs(end+1)=EEG.chanlocs(find(strcmpi(wanted,'Fp1'),1));
EEG.chanlocs(end).labels='EOG1';
EEG.chanlocs(end).type='EOG';
EEG.event=struct('type',cellstr(D.eventTypes(:)), ...
    'latency',num2cell(double(D.eventLatency(:))),'duration',0);
EEG=eeg_checkset(EEG);
EEG=eeg_checkset(EEG,'eventconsistency');
EEG=eeg_checkset(EEG,'makeur');
EEG.history='';
EEG.setname='ERP_CORE_P3_support_extension';
end
