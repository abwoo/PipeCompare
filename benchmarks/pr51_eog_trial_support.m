function pr51_eog_trial_support(dataDir, outDir, subject)
% PR51 exploration only: fixed minimal pipeline, with and without existing
% EOG regression. Never rewrites the frozen protocol, flags or eligibility.
% The fit removes event-locked ERP averages before estimating the EOG-to-EEG
% transfer coefficients (Steps.eogRegress); assess the P3 measurements as
% well as the remaining HMM-unflagged trial support.
if nargin==0
    assert(strcmp(conditionCode('11'),'11'));
    assert(strcmp(conditionCode('31'),'31'));
    fprintf('PR51_EOG_SUPPORT_SELF_CHECK_PASSED\n');
    return;
end
if ~isfolder(outDir), mkdir(outDir); end
base=fileparts(mfilename('fullpath'));
p=jsondecode(fileread(fullfile(base,'erpcore_protocol.json')));
D=load(fullfile(dataDir,'input.mat'));
EEG=buildEeg(D,p);
assert(subject==round(subject),'Subject must be an integer.');
result=struct('subject',subject,'purpose','exploratory EOG-regression feasibility, NOT confirmatory validation', ...
    'implementationCommit',getenv('GITHUB_SHA'),'protocolUnchanged',true, ...
    'thresholdsUnchanged',true,'sourceEventIdentityPreserved',true, ...
    'arms',struct(),'errors',{{}},'completed',false);
saveJson(fullfile(outDir,'eog_experiment.json'),result);
for mode={'minimal','eog_regressed'}
    arm=mode{1};
    try
        setBase(EEG);
        con=pipecompare.eval.Contract('conditions',{'target',{'11'};'standard',{'31'}}, ...
          'epoch',[-0.2 0.8],'baseline',[-0.2 0], ...
          'components',{'P3',[0.3 0.6],{'Pz'}});
        plan=pipecompare.plan.Plan();
        plan=plan.add('highpass','cutoff',0.1);
        plan=plan.add('lowpass','cutoff',30);
        if strcmp(arm,'eog_regressed')
            plan=plan.add('eogregress','channels',{{'EOG1'}});
        end
        plan=plan.add('reref','mode','average','exclude',{{'EOG1'}});
        plan=plan.add('epoch'); plan=plan.add('baseline');
        plan=plan.add('reject_threshold','uv',100,'exclude',{{'EOG1'}});
        transcript=evalc('r=pipecompare.PipeCompare.optimize(plan,con,struct(''nSplits'',0));');
        saveText(fullfile(outDir,[arm '.log']),transcript);
        assert(numel(r.cands)==1,'Unexpected exploratory pipeline multiplicity.');
        c=r.cands(1);
        out=struct('candidateStatus',c.status, 'scientificStatus',char(r.ranking.table.status{1}), ...
          'eyeRegressionEnabled',strcmp(arm,'eog_regressed'), 'estimate',[], ...
          'differenceAmplitude',NaN, 'sme',[], 'support',[]);
        if strcmp(c.status,'ok')
            m=c.m; o=m.objectives(1);
            out.estimate=o.estimate; out.sme=o.sme;
            if numel(o.estimate)==2 && all(isfinite(o.estimate))
                out.differenceAmplitude=o.estimate(1)-o.estimate(2);
            end
            out.support=reportSupport(r,m,o);
        end
        result.arms.(arm)=out;
    catch ex
        err=struct('arm',arm,'identifier',ex.identifier,'message',ex.message, ...
                   'trace',getReport(ex,'extended','hyperlinks','off'));
        result.errors{end+1}=err;
        fprintf(2,'PR51_EOG_EXPLORATORY_ERROR %s\n',jsonencode(err));
    end
    saveJson(fullfile(outDir,'eog_experiment.json'),result);
end
result.completed=isempty(result.errors) && isfield(result.arms,'minimal') && ...
    isfield(result.arms,'eog_regressed');
saveJson(fullfile(outDir,'eog_experiment.json'),result);
fprintf('PR51_EOG_EXPLORATORY_RESULT %s\n',jsonencode(result));
end

function EEG=buildEeg(D,p)
% Same source channel and event assembly as the frozen real-P3 recorder.
wanted=cellstr(p.excerpt.channels);labels=cellstr(D.labels(:));
[yes,loc]=ismember(lower(wanted),lower(labels));assert(all(yes),'Missing source channels.');
cap=readlocs(fullfile(fileparts(which('eeglab')),'functions','supportfiles','Standard-10-5-Cap385.sfp'));
[yes,at]=ismember(lower(wanted),lower({cap.labels}));assert(all(yes),'Missing montage sites.');
EEG=eeg_emptyset();EEG.data=double(D.data(loc,:));EEG.srate=double(D.fs);
EEG.nbchan=numel(wanted);EEG.pnts=size(EEG.data,2);EEG.trials=1;
EEG.xmin=0;EEG.xmax=(EEG.pnts-1)/EEG.srate;EEG.chanlocs=cap(at);
[EEG.chanlocs.type]=deal('EEG');
EEG.data(end+1,:)=double(D.veog);EEG.nbchan=EEG.nbchan+1;
EEG.chanlocs(end+1)=EEG.chanlocs(find(strcmpi(wanted,'Fp1'),1));
EEG.chanlocs(end).labels='EOG1';EEG.chanlocs(end).type='EOG';
EEG.event=struct('type',cellstr(D.eventTypes(:)), ...
  'latency',num2cell(double(D.eventLatency(:))),'duration',0);
EEG=eeg_checkset(EEG);EEG=eeg_checkset(EEG,'eventconsistency');
EEG=eeg_checkset(EEG,'makeur');
EEG.history='';EEG.setname='ERP_CORE_P3_EOG_exploration';
end

function support=reportSupport(r,m,o)
support=struct([]);
for k=1:numel(r.ref.ids)
    names=r.ref.names; if isstring(names),names=cellstr(names);end
    rawName=names{k};
    code=conditionCode(rawName);
    xx=o.X{k}(:,1);
    flagged=logical(r.ref.blink{1}{k}(:));
    w=ones(size(xx));
    if ~isempty(o.w), w=o.w{k}(:); end
    if numel(xx)~=numel(w) || numel(xx)~=numel(flagged)
        error('PR51:EogAlignment','Mismatched frozen source trial vectors.');
    end
    good=isfinite(xx)&isfinite(w)&w>0;w(~good)=0;
    clean=w(~flagged);
    en=sum(clean)^2/max(sum(clean.^2),eps);
    z=struct('condition',code,'referenceN',numel(r.ref.ids{k}), ...
        'sourceUnflaggedN',nnz(~flagged), ...
        'kept',m.kept(k),'effectiveUnflaggedN',en, ...
        'hmmFlaggedN',nnz(flagged),'positiveN',nnz(good), ...
        'referenceIds',r.ref.ids{k},'scores',xx','weights',w','flagged',flagged');
    if isempty(support),support=z;else,support(end+1)=z;end %#ok<AGROW>
end
end

function code=conditionCode(name)
if any(strcmp(char(string(name)),{'11','target'})),code='11';
elseif any(strcmp(char(string(name)),{'31','standard'})),code='31';
else,error('PR51:EogCondition','Unknown P3 condition');
end
end

function setBase(EEG)
assignin('base','EEG',EEG);assignin('base','ALLEEG',EEG);assignin('base','CURRENTSET',1);
end

function saveJson(file,v)
fid=fopen([file '.tmp'],'w');assert(fid>=0,'Cannot save EOG diagnostic JSON.');
c=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(v));clear c;
movefile([file '.tmp'],file,'f');
end

function saveText(file,s)
fid=fopen(file,'w');assert(fid>=0);c=onCleanup(@()fclose(fid));fprintf(fid,'%s',s);
end
