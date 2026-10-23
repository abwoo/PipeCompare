function S=pipecompare_record_summary(r,name)
% Persist only candidate scores/metadata, never raw EEG. The empty schema
% must have the same fields as populated summaries (MATLAB struct arrays).
S=struct('name',{},'keys',{},'labels',{},'table',{},'cands',{},'ref',{},'ranking',{});
if nargin==0,return;end
keep=intersect(fieldnames(r.cands),{'id','key','stratum','status','m','signal','interpolatedFraction'});
cands=rmfield(r.cands,setdiff(fieldnames(r.cands),keep));
[~,order]=sort([cands.id]);cands=cands(order);
S=struct('name',name,'keys',{{r.leaves.key}},'labels',{r.labels},'table',r.ranking.table, ...
    'cands',{cands},'ref',r.ref,'ranking',rmfield(r.ranking,setdiff(fieldnames(r.ranking),{'objective','options'})));
end
