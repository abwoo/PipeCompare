function out=audit_hmm(E,control,truth,outDir)
fs=E.srate;n=E.pnts;e=find(strcmp({E.chanlocs.labels},'EOG1'));x=double(E.data(e,:));xc=double(control.data(e,:));
[on,model]=pipecompare.eval.Hmm.events(x,fs);[controlOn,controlModel]=pipecompare.eval.Hmm.events(xc,fs);
d=max(1,floor(fs/25));r=x-movmedian(x,round(fs),2);smooth=movmean(r,max(1,round(fs/16)),2);a=abs(smooth(1:d:n)-median(smooth(1:d:n)));z=a/max(1.4826*median(a),eps);
[check,~]=pipecompare.eval.Hmm.segment(z,pipecompare.eval.Hmm.PeakZ,max(1,round(pipecompare.eval.Hmm.MinSeconds*fs/d)));
check=repelem(check,d);assert(isequal(check(1:n),on),'HMM decomposition differs');
zFull=repelem(z,d);zFull=zFull(1:n);
w=round(.05*fs);padded=conv(double(on),ones(1,2*w+1),'same')>0;
controlPadded=conv(double(controlOn),ones(1,2*w+1),'same')>0;
y=abs(r-median(r));robust=y>5*max(1.4826*median(y),eps);
robustPadded=conv(double(robust),ones(1,2*w+1),'same')>0;
lat=round([E.event.latency]);spans=intervals(on);padspans=intervals(padded);q={};
for id=[3 31 38]
 c=lat(id);bl=c+round(-.2*fs):c;win=c+round(.3*fs):c+round(.6*fs);sel=unique([bl win]);
 hits=find(any(spans(:,1)<=sel & spans(:,2)>=sel,2));phits=find(any(padspans(:,1)<=sel & padspans(:,2)>=sel,2));
 h=struct('trial',id,'eventTimeSeconds',(c-1)/fs,'injectedBlink',truth(id),'sampleRate',fs, ...
 'baselineRawHmmSamples',nnz(on(bl)),'windowRawHmmSamples',nnz(on(win)), ...
 'baselinePaddedHmmSamples',nnz(padded(bl)),'windowPaddedHmmSamples',nnz(padded(win)), ...
 'rawHmmOverlap',any(on(sel)),'paddedHmmOverlap',any(padded(sel)), ...
 'dilationOnly',~any(on(sel))&&any(padded(sel)),'robustOverlap',any(robustPadded(sel)), ...
 'controlRawOverlap',any(controlOn(sel)),'controlPaddedOverlap',any(controlPadded(sel)), ...
 'eogMatchesNoInjectionControlInWindows',max(abs(x(sel)-xc(sel)))==0, ...
 'maxEogDifferenceInWindows',max(abs(x(sel)-xc(sel))),'rawIntervals',{{}},'paddedIntervalsRelative',((padspans(phits,:)-c)/fs));
 for j=hits(:)'
  ix=spans(j,1):spans(j,2);[nearest,at]=min(abs(lat(truth)+.42*fs-mean(ix)));
  injectedIds=find(truth);
  h.rawIntervals{end+1}=struct('relativeStart',(ix(1)-c)/fs,'relativeEnd',(ix(end)-c)/fs, ...
   'durationSeconds',numel(ix)/fs,'maxRobustZ',max(zFull(ix)), ...
   'nearestInjectedTrial',injectedIds(at),'distanceToNearestInjectedCenterSeconds',nearest/fs);
 end
 assert(~truth(id)&&h.paddedHmmOverlap&&~h.robustOverlap);
 q{end+1}=h;fprintf('HMM_WINDOW_AUDIT %s\n',jsonencode(h));
end
T=table((0:n-1)'/fs,x',xc',r',smooth',zFull',on',padded',robustPadded',controlOn',controlPadded', ...
 'VariableNames',{'timeSeconds','eogUV','noInjectedEyeEogUV','detrendedUV','smoothedUV','hmmRobustZ','hmmRaw','hmmPadded','robustPadded','controlHmmRaw','controlHmmPadded'});
writetable(T,fullfile(outDir,'hmm_trace.csv'));
write(fullfile(outDir,'events.json'),struct('latencySamples',lat,'injectedBlink',truth));
out=struct('cases',{q},'model',model,'noInjectedEyeModel',controlModel,'hmmFeatureDownsample',d,'paddingSamples',w);
write(fullfile(outDir,'hmm_audit.json'),out);
end
function s=intervals(mask)
z=diff([0 mask 0]);s=[find(z==1)' find(z==-1)'-1];
end
function write(f,x)
h=fopen(f,'w');assert(h>=0);fprintf(h,'%s\n',jsonencode(x));fclose(h);
end
