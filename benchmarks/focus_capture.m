function focus_capture(stage,result,cands,ranking,ref)
% Observe only; never modify ranking or use truth to choose a recipe.
cfg=getappdata(0,'PipeCompareFocus');assert(~isempty(cfg));
r=struct('cands',cands,'ranking',ranking,'ref',ref,'labels',{result.labels});
rows=cell(1,numel(cands));
for k=1:numel(cands)
 q=struct('id',k,'pipeline',result.labels{k},'ranking',table2struct(ranking.table(k,:)), ...
  'selected',any(ranking.recommended==k),'available',false);
 if strcmp(cands(k).status,'ok')
  q.readout=focus_readout(r,k,cfg.population,cfg.realized);q.available=q.readout.available;
  o=cands(k).m.objectives(1);
  q.estimate=o.estimate;q.sme=o.sme;q.eyeBias=o.eyeBias;q.eyeBiasSq=o.eyeBiasSq;
  q.retention=cands(k).m.retention;q.kept=cands(k).m.kept;
  q.decisions=cands(k).steps;
  for z=1:numel(q.decisions)
   if isfield(q.decisions{z},'curves'),q.decisions{z}=rmfield(q.decisions{z},'curves');end
  end
 end
 rows{k}=q;
end
summary=struct('stage',stage,'subject',cfg.subject,'seed',cfg.seed,'arm',cfg.arm, ...
 'recommended',ranking.recommended,'populationTruth',cfg.population,'realizedTruth',cfg.realized,'rows',{rows});
fid=fopen(fullfile(cfg.outDir,[stage '.json']),'w');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(summary));fclose(fid);
% Keep score arrays for opposite-formula re-ranking, no raw EEG or ICA matrices.
for k=1:numel(cands)
 cands(k).ica={};cands(k).coms={};
end
levels=cell(numel(result.leaves),1);
for k=1:numel(result.leaves)
 p=result.leaves(k).path;
 levels{k}=[cellfun(@(i)i.slot,p(:),'UniformOutput',false),cellfun(@(i)i.key,p(:),'UniformOutput',false)];
end
options=result.options;options.progress=[];labels=result.labels;
save(fullfile(cfg.outDir,[stage '.mat']),'cands','ref','ranking','levels','options','labels','cfg','-v7');
end
