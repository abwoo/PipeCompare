function EEG=pipecompare_load_public_set(path)
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

