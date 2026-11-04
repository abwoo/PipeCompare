function records=audit_rejection(E,folder,outDir)
contract=pipecompare.simple.Presets.contract(E,'P3',{'11','31'});
root=pipecompare.run.Executor.prepareRoot(E,contract,struct('dataUnit','uV'));
records={};prefix=[];prefixSteps={};
for arm={'baseline','variant'}
 d=load(fullfile(folder,arm{1},'after_eog.mat'));k=d.ranking.recommended;c=d.cands(k);
 assert(d.cfg.seed==0&&d.cfg.subject==5&&strcmp(d.cfg.arm,arm{1}));
 ctx=struct('contract',contract,'highpass',0,'icaSeed',0,'ref',d.ref,'rank',d.options);
 G=root;C=[];stages={};stageNames={};
 for j=1:numel(c.steps)
  st=c.steps{j};in=struct('type',st.type,'params',st.params);
  if strcmp(in.type,'badchannels')
   assert(strcmp(in.params.action,'interpolate'));
   [yes,idx]=ismember(lower(st.interpolated),lower({G.chanlocs.labels}));assert(all(yes));
   if ~isempty(idx),G=pop_interp(G,idx,'spherical');end
   if ~isempty(st.stretches),G=pipecompare.run.Steps.interpStretches(G,st.stretches);end
  elseif strcmp(in.type,'ica')
   if isempty(prefix)
    G=pipecompare.run.Steps.run(in,G,ctx);prefix=G;prefixSteps=c.steps(1:j);
   else
    for t=1:j,assert(isequaln(c.steps{t}.params,prefixSteps{t}.params));end
    assert(isequal(G.data,prefix.data),'ICA prefix input differs');G=prefix;
   end
  elseif strcmp(in.type,'eogregress')
   before=epochBaseline(G,contract);stages{end+1}=before;stageNames{end+1}='before_eog';C=G;
   G=pipecompare.run.Steps.replayDecision(in,G,st,ctx);
   stages{end+1}=epochBaseline(G,contract);stageNames{end+1}='after_eog';
  elseif strcmp(in.type,'icremove')
   G=pipecompare.run.Steps.replayDecision(in,G,st,ctx);
   if ~isempty(C),C=pipecompare.run.Steps.replayDecision(in,C,st,ctx);end
  elseif strcmp(in.type,'reject_threshold')
   assert(strcmp(in.params.method,'absolute')&&in.params.interpolate==0);
   assert(isfinite(st.uvChosen)&&st.uvChosen>0);uv=st.uvChosen;in.params.uv=uv;
   stages{end+1}=G;stageNames{end+1}='before_rejection';
   stages{end+1}=C;stageNames{end+1}='without_eog_same_other_decisions';
   T=pipecompare.eval.Measure.trials(G,contract,true);ep=find(T.id==23);assert(isscalar(ep));ep=T.epoch(ep);
   chans=pipecompare.run.Steps.testedChannels(G,in.params);
   peak=max(abs(double(G.data(chans,:,ep))),[],2);over=chans(peak>uv);
   perChannel={};times=G.xmin+(0:G.pnts-1)/G.srate;
   for ci=chans(:)'
    [mx,at]=max(abs(double(G.data(ci,:,ep))));
    perChannel{end+1}=struct('channel',G.chanlocs(ci).labels,'maxAbsUV',mx,'signedUV',double(G.data(ci,at,ep)), ...
     'timeSeconds',times(at),'thresholdUV',uv,'excessUV',mx-uv,'crosses',mx>uv);
   end
   out=struct('arm',arm{1},'seed',0,'candidate',k,'thresholdUV',uv,'removedComponents',[], ...
    'crossingChannels',{{G.chanlocs(over).labels}},'channelPeaks',{perChannel},'stagePeaks',{{}}, ...
    'savedRejectedIds',T.id(st.rejIdx)');
   ic=find(cellfun(@(x)strcmp(x.type,'icremove'),c.steps));out.removedComponents=c.steps{ic}.comps;
   for s=1:numel(stages)
    S=stages{s};tt=pipecompare.eval.Measure.trials(S,contract,true);e=tt.epoch(tt.id==23);assert(isscalar(e));
    cc=pipecompare.run.Steps.testedChannels(S,in.params);v=double(S.data(cc,:,e));
    [val,idx]=max(abs(v(:)));[ch,tm]=ind2sub(size(v),idx);ts=S.xmin+(0:S.pnts-1)/S.srate;
    out.stagePeaks{s}=struct('stage',stageNames{s},'maxAbsUV',val,'channel',S.chanlocs(cc(ch)).labels, ...
     'timeSeconds',ts(tm),'crossesSameFrozenThreshold',val>uv,'perChannelMaxAbsUV',max(abs(v),[],2)', ...
     'channels',{{S.chanlocs(cc).labels}});
    [ti,chi]=meshgrid(ts,1:numel(cc));labels={S.chanlocs(cc).labels};
    tab=table(repmat(stageNames(s),numel(v),1),reshape(labels(chi(:)),[],1),ti(:),v(:), ...
     'VariableNames',{'stage','channel','timeSeconds','valueUV'});
    fn=fullfile(outDir,['trial23_' arm{1} '_' stageNames{s} '.csv']);writetable(tab,fn);
   end
   [G,~,info]=pipecompare.run.Steps.run(in,G,ctx);
   assert(isequal(sort(info.rejIdx(:)),sort(st.rejIdx(:))),'Frozen threshold did not reproduce rejection');
   m=pipecompare.eval.Measure.candidate(G,contract,d.ref,d.options);original=c.m.objectives(1);obs=m.objectives(1);
   err=max(abs(original.estimate(:)-obs.estimate(:)));
   for cond=1:numel(original.X)
    assert(isequal(isnan(original.X{cond}),isnan(obs.X{cond})));
    diff=abs(original.X{cond}-obs.X{cond});diff=diff(isfinite(diff));err=max([err;diff(:)]);
   end
   assert(err<1e-6,'Saved selected output not reproduced');out.replayMaxError=err;out.rejectionExactlyReproduced=true;
   records{end+1}=out;fprintf('REJECTION_WAVEFORM_AUDIT %s\n',jsonencode(out));
  else
   G=pipecompare.run.Steps.run(in,G,ctx);
   if ~isempty(C),C=pipecompare.run.Steps.run(in,C,ctx);end
  end
  if strcmp(in.type,'highpass'),ctx.highpass=max(ctx.highpass,in.params.cutoff);end
 end
end
f=fopen(fullfile(outDir,'rejection_audit.json'),'w');assert(f>=0);fprintf(f,'%s\n',jsonencode(records));fclose(f);
end
function G=epochBaseline(G,contract)
ctx=struct('contract',contract);G=pipecompare.run.Steps.run(struct('type','epoch','params',struct()),G,ctx);
G=pipecompare.run.Steps.run(struct('type','baseline','params',struct()),G,ctx);
end
