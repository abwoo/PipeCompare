function T=pipecompare_open_blinks(D,name,margin)
X=D.data; valid=logical(D.valid); X(:,~valid)=NaN;
[hit,model]=pipecompare.eval.Hmm.events(X,D.fs);
% An interval counts once. Matching is one-to-one, so a long detection
% spanning several annotated blinks cannot receive several true positives.
d=diff([false hit false]); a=find(d==1); b=find(d==-1)-1;
gt=D.annotations(:,1)'; matched=false(size(gt)); tp=0; timing=[];
for j=1:numel(a)
    start=D.time(a(j))-margin; stop=D.time(b(j))+margin;
    candidates=find(~matched & gt>=start & gt<=stop);
    if isempty(candidates), continue; end
    centre=(D.time(a(j))+D.time(b(j)))/2;
    [~,at]=min(abs(gt(candidates)-centre)); at=candidates(at);
    matched(at)=true; tp=tp+1; timing(end+1)=abs(gt(at)-centre); %#ok<AGROW>
end
precision=tp/max(numel(a),1); recall=tp/max(numel(gt),1);
f1=2*precision*recall/max(precision+recall,eps);
minutes=nnz(valid)/D.fs/60; falseEvents=numel(a)-tp;
T=table({name},numel(gt),numel(a),tp,falseEvents,precision,recall,f1, ...
    median(timing),falseEvents/max(minutes,eps),mean(~valid),mean(hit(valid)),model.valid, ...
    D.effective_fs,D.timestamp_jitter_seconds, ...
    'VariableNames',{'recording','annotations','detections','truePositive','falsePositive','precision','recall','f1', ...
    'medianTimingErrorSeconds','falseEventsPerMinute','unavailableShare','activeShare','assessed','effectiveFs','timestampJitterSeconds'});
end

