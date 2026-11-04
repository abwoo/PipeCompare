function audit_support(dataRoot,fragments,outDir)
% Observation only: regenerate frozen inputs, not fitted candidate searches.
if ~isfolder(outDir),mkdir(outDir);end
reports={};retention={};
for subject=[4 5]
 [E,population,realized,truth]=audit_input(fullfile(dataRoot,num2str(subject)),subject);
 contract=pipecompare.simple.Presets.contract(E,'P3',{'11','31'});
 E=pipecompare.run.Executor.prepareRoot(E,contract,struct('dataUnit','uV'));
 ref=pipecompare.eval.Measure.reference(E,contract);
 [flags,info]=audit_blinks(E,contract,ref);
 assert(isequal(flags,ref.blink),'Detector instrumentation changed flags');
 T=pipecompare.eval.Measure.trials(E,contract,true);
 assert(numel(unique(T.id))==numel(truth)&&isequal(sort(T.id(:)),(1:numel(truth))'));
 a=info.audit{1};target=truth(a.ids(:))';target=target(:);
 assert(isequal(a.ids(:),T.id(:))&&isequal(a.conditions(:),T.cond(:)));
 rows=[a.ids(:),a.conditions(:),target,a.epochKinds,a.sampleRobust,a.sampleHmm,a.finalKinds,a.combined];
 names={'trialId','condition','injectedBlink','epochBlink','epochHorizontal','epochMuscle','sampleRobustBlink','sampleRobustHorizontal','sampleHmmBlink','sampleHmmHorizontal','finalBlink','finalHorizontal','finalMuscle','combinedFlag'};
 writetable(array2table(rows,'VariableNames',names),fullfile(outDir,sprintf('S%03d_trial_flags.csv',subject)));
 q=struct('subject',subject,'trials',numel(truth),'injectedBlinks',nnz(truth),'flagged',nnz(a.combined), ...
  'channels',{ref.blinkInfo.channels},'horizontal',ref.blinkInfo.horizontal, ...
  'kindShare',ref.blinkInfo.kindShare,'noInjectedBlinkButFlaggedIds',a.ids(~target&a.combined)', ...
  'injectedBlinkButUnflaggedIds',a.ids(target&~a.combined)','conditions',{{}});
 for c=1:numel(ref.ids)
  use=a.conditions==c;b=a.combined(use);tr=target(use);
  q.conditions{c}=struct('condition',c,'referenceTrials',nnz(use),'unflagged',nnz(~b),'injectedBlinks',nnz(tr), ...
   'injectedNonblinks',nnz(~tr),'trueInjectedFlagged',nnz(tr&b),'injectedNotFlagged',nnz(tr&~b), ...
   'noInjectionFlagged',nnz(~tr&b),'noInjectionUnflagged',nnz(~tr&~b));
 end
 reports{end+1}=q;fprintf('REFERENCE_SUPPORT %s\n',jsonencode(q));
 for seed=0:4
  for arm={'baseline','variant'}
   for stage={'initial','after_cutoffs','after_eog'}
    d=load(fullfile(fragments,sprintf('eog-pair-%d-%d',subject,seed),arm{1},[stage{1} '.mat']));
    assert(d.cfg.subject==subject&&d.cfg.seed==seed&&strcmp(d.cfg.arm,arm{1}));
    assert(max(abs(d.cfg.population(:)-population(:)))<1e-12&&max(abs(d.cfg.realized(:)-realized(:)))<1e-12,'Input generator differs');
    assert(isequal(d.ref.ids,ref.ids)&&isequal(d.ref.blink,ref.blink),'Saved reference does not match regeneration');
    for k=1:numel(d.cands)
     cand=d.cands(k);if ~strcmp(cand.status,'ok'),continue;end
     ob=cand.m.objectives(1);keptIds=[];
     entry=struct('subject',subject,'seed',seed,'arm',arm{1},'stage',stage{1},'candidate',k, ...
      'selected',any(d.ranking.recommended==k),'pipeline',d.labels{k},'conditions',{{}},'rejectionSteps',{{}});
     for c=1:numel(ref.ids)
      ids=ref.ids{c}(:);b=ref.blink{1}{c}(:);tr=truth(ids);tr=tr(:);present=~isnan(ob.X{c}(:,1));
      wt=ones(size(present));if ~isempty(ob.w),wt=ob.w{c};end
      usable=present&isfinite(wt)&wt>0;keptIds=[keptIds;ids(present)]; %#ok<AGROW>
      entry.conditions{c}=struct('condition',c,'referenceTrials',numel(ids),'referenceUnflagged',nnz(~b), ...
       'retained',nnz(present),'retainedUnflagged',nnz(present&~b),'positiveUnflagged',nnz(usable&~b), ...
       'removedUnflaggedIds',ids(~present&~b)','retainedUnflaggedIds',ids(present&~b)', ...
       'removedInjectedBlink',nnz(~present&tr),'removedNoInjectedBlink',nnz(~present&~tr), ...
       'retainedInjectedBlink',nnz(present&tr),'retainedNoInjectedBlink',nnz(present&~tr));
     end
     attributed=[];mapped=false;
     for j=1:numel(cand.steps)
      st=cand.steps{j};if ~isstruct(st)||~isfield(st,'rejIdx'),continue;end
      z=struct('step',j,'epochsBefore',st.epochsBefore,'rejected',numel(st.rejIdx),'indices',st.rejIdx,'fields',{fieldnames(st)});
      if isfield(st,'topChannels'),z.topChannels=st.topChannels;end
      if isfield(st,'keptShift'),z.keptShift=st.keptShift;end
      if st.epochsBefore==numel(T.id)
       z.ureventIds=T.id(st.rejIdx)';attributed=[attributed;T.id(st.rejIdx)];mapped=true; %#ok<AGROW>
      end
      entry.rejectionSteps{end+1}=z;
     end
     entry.absentIds=setdiff(T.id,keptIds)';
     entry.rejectionIdsExplainAbsence=mapped&&isequal(sort(unique(attributed(:))),sort(entry.absentIds(:)));
     if entry.selected,assert(entry.rejectionIdsExplainAbsence,'Selected absence not explained by rejection records');end
     retention{end+1}=entry;
     if entry.selected
      brief=rmfield(entry,{'pipeline','rejectionSteps'});fprintf('SELECTED_SUPPORT %s\n',jsonencode(brief));
     end
    end
   end
  end
 end
end
write(fullfile(outDir,'reference.json'),reports);write(fullfile(outDir,'retention.json'),retention);
fprintf('SUPPORT_AUDIT_CHECKS_PASSED subjects=2 pairs=10 stages=3 arms=2 candidateSnapshots=%d\n',numel(retention));
end
function write(path,x)
f=fopen(path,'w');assert(f>=0);fprintf(f,'%s\n',jsonencode(x));fclose(f);
end
