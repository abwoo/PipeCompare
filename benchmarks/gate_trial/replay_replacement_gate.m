function replay_replacement_gate(folder,outDir)
if ~isfolder(outDir),mkdir(outDir);end
check_replacement_gate();rows={};
for subject=[4 5]
 for seed=0:4
  p=fullfile(folder,sprintf('eog-pair-%d-%d',subject,seed));
  base=jsondecode(fileread(fullfile(p,'baseline','result.json')));
  expanded=jsondecode(fileread(fullfile(p,'variant','result.json')));
  assert(base.completed&&expanded.completed&&base.subject==subject&&base.seed==seed);
  assert(expanded.subject==subject&&expanded.seed==seed);
  assert(strcmp(base.executionCommit,'b05b1833597d76276337a9e93621293e6e18f510'));
  assert(strcmp(base.executionCommit,expanded.executionCommit));
  assert(isequal(base.source,expanded.source)&&isequal(base.populationTruth,expanded.populationTruth));
  d=load(fullfile(p,'variant','after_eog.mat'));
  assert(d.cfg.subject==subject&&d.cfg.seed==seed&&strcmp(d.cfg.arm,'variant'));
  assert(isequal(d.cfg.population(:),base.populationTruth(:)));
  ids=find(strcmp(d.ranking.table.status,'feasible'));
  old=find(strcmp(d.labels,base.erp.pipeline));assert(isscalar(old));
  ib=find(ids==old);ic=find(ids==d.ranking.recommended);assert(isscalar(ib)&&isscalar(ic));
  r=struct('cands',d.cands,'ranking',d.ranking,'ref',d.ref,'labels',{d.labels});
  oldout=focus_readout(r,old,d.cfg.population,d.cfg.realized);
  assert(max(abs(oldout.estimate(:)-base.erp.estimate(:)))<1e-7,'Incumbent output changed');
  B=pipecompare.eval.Rank.bootstrap(d.cands(ids),d.ref,d.options,{d.ranking.objective});
  % Frozen alpha .02 and 1999 replicates match existing scalar defaults.
  assert(size(B{1},2)==1999);alpha=.02;
  g=replacement_gate(d.ranking.table.objective(ids),B{1},ib,ic,alpha);
  selected=ids(g.selected);a=focus_readout(r,selected,d.cfg.population,d.cfg.realized);
  q=struct('subject',subject,'seed',seed,'baseline',base.erp,'expanded',expanded.erp, ...
   'gated',a,'gate',g,'selectedId',selected,'candidateIds',ids(:)', ...
   'populationTruth',d.cfg.population,'source',base.source,'sourceExecution',base.executionCommit);
  rows{end+1}=q;
  fprintf('REPLACEMENT_GATE_RECORD %s\n',jsonencode(q));
 end
end
fid=fopen(fullfile(outDir,'records.json'),'w');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(rows));fclose(fid);
end
