function [E,population,realized,truth]=audit_input(dataDir,subject)
p=jsondecode(fileread(fullfile(fileparts(fileparts(mfilename('fullpath'))),'public_pilot_protocol.json')));
D=load(fullfile(dataDir,'input.mat'));EEG=background(D,p);
[E,population,realized,~,truth]=inject(EEG,subject,p);
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

function [EEG,population,realized,prevalence,truth]=inject(EEG,subject,p)
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
canonicalTopo=topo-mean(topo);trueTrial=zeros(1,nt);truth=false(1,nt);
for j=1:nt
 target=strcmp(types{j},'11');mu=amp*(s.standard_amplitude_ratio+(1-s.standard_amplitude_ratio)*target);
 blink=rand(rs)<prevalence;truth(j)=blink;
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

