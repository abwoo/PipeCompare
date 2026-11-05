function run_validation(out)
% Calls unchanged production HMM. External blocks are NOT event-level truth.
rows={};
for sub=1:10
 for name={'open','blink'}
  tag=sprintf('S%03d_%s',sub,name{1});D=load(fullfile(out,[tag '.mat']));
  for channels=[1 2]
   X=D.X(1:channels,:);[h,m]=pipecompare.eval.Hmm.events(X,D.fs);
   r=X-movmedian(X,round(D.fs),2);a=abs(r-median(r,2));rob=any(a>5*max(1.4826*median(a,2),eps),1);
   pad=round(.05*D.fs);hp=conv(double(h),ones(1,2*pad+1),'same')>0;rp=conv(double(rob),ones(1,2*pad+1),'same')>0;
   row=struct('subject',sub,'block',name{1},'channels',channels,'hmmValid',m.valid, ...
    'hmmFraction',mean(h),'hmmPaddedFraction',mean(hp),'robustPaddedFraction',mean(rp), ...
    'unionFraction',mean(hp|rp),'hmmOnlyFraction',mean(hp&~rp),'hmmSegments',nnz(diff([false h])==1));
   rows{end+1}=row;fprintf('EXTERNAL_BLOCK %s\n',jsonencode(row));
   writetable(table((0:numel(h)-1)'/D.fs,h',hp',rp','VariableNames',{'seconds','hmm','hmmPadded','robustPadded'}),fullfile(out,sprintf('%s_ch%d.csv',tag,channels)));
  end
 end
end
write(fullfile(out,'external_summary.json'),rows);
controls={};fs=250;t=(0:fs*60-1)/fs;centers=5:5:50;truth=false(size(t));for c=centers,truth=truth|abs(t-c)<=.3;end
for seed=0:9
 rng(1000+seed,'twister');noise=filter(1,[1 -.8],randn(size(t)));noise=noise/std(noise);base=noise+sin(2*pi*10*t);
 for amp=[0 4 8 16]
  x=base;
  for c=centers
   ix=abs(t-c)<=.3; x(ix)=x(ix)+amp*.5*(1+cos(pi*(t(ix)-c)/.3));
  end
  [h,m]=pipecompare.eval.Hmm.events(x,fs);pad=round(.05*fs);h=conv(double(h),ones(1,2*pad+1),'same')>0;
  detected=false(size(centers));for j=1:numel(centers),detected(j)=any(h(abs(t-centers(j))<=.1));end
  if amp==0,neg=true(size(t));hits=0;total=0;else,neg=~truth;hits=nnz(detected);total=numel(centers);end
  row=struct('seed',seed,'amplitude',amp,'hmmValid',m.valid,'hitEvents',hits,'totalEvents',total,'outsideSupportFraction',mean(h(neg)),'markedFraction',mean(h));
  controls{end+1}=row;fprintf('KNOWN_TRUTH_CONTROL %s\n',jsonencode(row));
 end
end
write(fullfile(out,'known_truth_summary.json'),controls);
assert(numel(rows)==40&&numel(controls)==40);fprintf('DETECTOR_VALIDATION_EXECUTION_PASSED\n');
end
function write(f,x)
h=fopen(f,'w');assert(h>=0);fprintf(h,'%s\n',jsonencode(x));fclose(h);
end
