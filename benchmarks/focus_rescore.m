function focus_rescore(inputDir,outDir,expectedArm)
% Re-score exactly the same fitted candidate outputs under the active formula.
if ~isfolder(outDir),mkdir(outDir);end
for stage={'initial','after_cutoffs','after_eog'}
 d=load(fullfile(inputDir,[stage{1} '.mat']));assert(strcmp(d.cfg.arm,expectedArm));
 cands=d.cands;
 for k=1:numel(cands)
  if ~strcmp(cands(k).status,'ok'),continue;end
  objs=cands(k).m.objectives;
  for o=1:numel(objs)
   if ~strcmp(objs(o).kind,'scalar')||numel(d.ref.blink)<o||isempty(d.ref.blink{o}),continue;end
   for c=1:numel(objs(o).estimate)
    K=ones(d.ref.n(c),1);
    if ~isempty(objs(o).w),K=objs(o).w{c};K(isnan(K))=0;end
    [objs(o).eyeBias(c),objs(o).eyeBiasSq(c),objs(o).eyeBiasAssessed(c)]= ...
     pipecompare.eval.Measure.eyeBias(objs(o).X{c}(:,1),d.ref.blink{o}{c},K);
   end
  end
  cands(k).m.objectives=objs;
 end
 R=pipecompare.eval.Rank.run(cands,d.ref,d.options,d.levels);
 r=struct('cands',cands,'ranking',R,'ref',d.ref,'labels',{d.labels});
 rows=cell(1,numel(cands));
 for k=1:numel(cands)
  q=struct('id',k,'pipeline',d.labels{k},'ranking',table2struct(R.table(k,:)), ...
    'selected',any(R.recommended==k),'available',false);
  if strcmp(cands(k).status,'ok'),q.readout=focus_readout(r,k,d.cfg.population,d.cfg.realized);q.available=q.readout.available;end
  rows{k}=q;
 end
 out=struct('stage',stage{1},'sourceArm',expectedArm,'subject',d.cfg.subject,'seed',d.cfg.seed, ...
  'recommended',R.recommended,'populationTruth',d.cfg.population,'realizedTruth',d.cfg.realized,'rows',{rows});
 fid=fopen(fullfile(outDir,[stage{1} '.json']),'w');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(out));fclose(fid);
end
fprintf('FOCUS_RESCORE_PASSED sourceArm=%s\n',expectedArm);
end
