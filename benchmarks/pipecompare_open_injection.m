function [EEG,population,realized,prevalence]=pipecompare_open_injection(EEG,subject)
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
