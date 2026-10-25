function pipecompare_optimization_public(dataDir,outDir,subject,stage,legacyDir)
% Paired actual-program comparison. Truth never enters either search.
if nargin<5,legacyDir='';end
if ~isfolder(outDir),mkdir(outDir);end
protocol=jsondecode(fileread(fullfile(fileparts(mfilename('fullpath')),'optimization_protocol.json')));
assert(ismember(subject,protocol.([stage '_subjects'])),'PipeCompare:Validation','Participant not assigned to this stage.');
copyfile(fullfile(dataDir,'frozen_inputs.json'),fullfile(outDir,'frozen_inputs.json'));
D=load(fullfile(dataDir,sprintf('S%03dR01.mat',subject)));
background=pipecompare_open_background(D);
[EEG,population,realized,prevalence]=pipecompare_open_injection(background,subject);
meta=struct('commit',getenv('GITHUB_SHA'),'stage',stage,'subject',subject, ...
    'prevalence',prevalence,'fitTimer','explicit-tic-handle','protocol',protocol, ...
    'matlab',version,'eeglab',eeg_getversion);
writeJson(fullfile(outDir,'metadata.json'),meta);
T=pipecompare_open_spectra(D,subject);writetable(T,fullfile(outDir,'spectral.csv'));
T=pipecompare_optimization_stretches(background,subject);writetable(T,fullfile(outDir,'stretches.csv'));
for recipe={'standard','balanced'}
    folder=fullfile(outDir,recipe{1});if ~isfolder(folder),mkdir(folder);end
    varianceMode='pooled';
    if strcmp(recipe{1},'standard')
        varianceMode='legacy-total';
        if strcmp(stage,'development')&&validation_reuse_legacy(legacyDir,folder,'public',subject,population,realized),continue;end
    end
    setBase(EEG);fitTimer=tic;
    transcript=evalc('[~,~,r]=pop_pipecompare(EEG,''measure'',''P3'',''events'',{''11'',''31''},''recipe'',recipe{1},''artifactvariance'',varianceMode,''reference'',''average'',''show'',''off'');');
    seconds=toc(fitTimer);fid=fopen(fullfile(folder,'native.log'),'w');fprintf(fid,'%s',transcript);fclose(fid);
    feasible=find(strcmp(r.ranking.table.status,'feasible'))';
    C=pipecompare.eval.Calibration.rows(r,feasible,population,realized,subject);
    writetable(C,fullfile(folder,'calibration.csv'));
    item=pipecompare_record_summary(r,sprintf('S%03dR01.mat',subject));
    save(fullfile(folder,'score.mat'),'item','-v7');
    pick=r.ranking.recommended;replay=NaN;
    if ~isempty(pick)
        transcript=evalc('E=pipecompare.run.Executor.replay(r,pick(1),true);'); %#ok<NASGU>
        m=pipecompare.eval.Measure.candidate(E,r.contract,r.ref,r.options);
        replay=max(abs(m.objectives(1).estimate-r.cands(pick(1)).m.objectives(1).estimate));
        assert(replay<1e-6,'PipeCompare:Validation','Native replay changed the selected readout.');
    end
    writeJson(fullfile(folder,'execution.json'),struct('completed',true,'recommended',~isempty(pick), ...
        'seconds',seconds,'replayMaxError',replay,'candidates',numel(r.cands), ...
        'artifactVariance',varianceMode,'legacyReused',false));
end
writeJson(fullfile(outDir,'execution.json'),struct('completed',true));
end

function setBase(EEG)
assignin('base','EEG',EEG);assignin('base','ALLEEG',EEG);assignin('base','CURRENTSET',1);
assignin('base','ALLCOM','');assignin('base','STUDY',[]);assignin('base','CURRENTSTUDY',0);
end

function writeJson(path,value)
fid=fopen(path,'w');assert(fid>=0);finish=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(value));
end
