function S=pipecompare_record_summary(r,name)
% Persist only candidate scores/metadata, never raw EEG. The empty schema
% must have the same fields as populated summaries (MATLAB struct arrays).
S=struct('name',{},'keys',{},'labels',{},'table',{},'cands',{},'ref',{},'ranking',{},'levels',{});
if nargin==0,return;end
keep=intersect(fieldnames(r.cands),{'id','key','stratum','status','message','m','signal','interpolatedFraction'});
cands=rmfield(r.cands,setdiff(fieldnames(r.cands),keep));
[~,order]=sort([cands.id]);cands=cands(order);
levels=cell(numel(cands),1);
for k=1:numel(cands)
    i=order(k);
    if isfield(r.leaves,'path')
        p=r.leaves(i).path;
        levels{k}=[cellfun(@(v)v.slot,p(:),'UniformOutput',false),cellfun(@(v)v.key,p(:),'UniformOutput',false)];
    end
    st=pipecompare.utils.fieldOr(r.cands(i),'steps',{});
    if any(cellfun(@(v)pipecompare.utils.fieldOr(v,'cutoffRefined',false),st)),levels{k}={};end
end
S=struct('name',name,'keys',{{r.leaves.key}},'labels',{r.labels},'table',r.ranking.table, ...
    'cands',{cands},'ref',r.ref,'ranking',rmfield(r.ranking,setdiff(fieldnames(r.ranking),{'objective','options','byStratum'})),'levels',{levels});
end
