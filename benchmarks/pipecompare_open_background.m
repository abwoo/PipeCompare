function EEG=pipecompare_open_background(D)
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
