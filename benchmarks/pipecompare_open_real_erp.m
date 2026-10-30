function T=pipecompare_open_real_erp(EEG,name,outDir)
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

function saveTranscript(outDir,name,transcript)
fid=fopen(fullfile(outDir,[name '.log']),'w');fprintf(fid,'%s',transcript);fclose(fid);
end
