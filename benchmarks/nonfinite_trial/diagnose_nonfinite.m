function diagnose_nonfinite(folder,previousDir,outDir)
% Read-only reconstruction using the same paired counts and active formula.
if ~isfolder(outDir),mkdir(outDir);end
prior=jsondecode(fileread(fullfile(previousDir,'records.json')));summaries={};events={};
for subject=[4 5]
 for seed=0:4
  p=fullfile(folder,sprintf('eog-pair-%d-%d',subject,seed));
  d=load(fullfile(p,'variant','after_eog.mat'));
  base=jsondecode(fileread(fullfile(p,'baseline','result.json')));
  assert(d.cfg.subject==subject&&d.cfg.seed==seed&&strcmp(d.cfg.arm,'variant'));
  assert(strcmp(base.executionCommit,'b05b1833597d76276337a9e93621293e6e18f510'));
  opts=pipecompare.utils.withDefaults(d.options,pipecompare.eval.Rank.defaults());
  assert(opts.nBoot==1999&&numel(d.ref.objectives)==1);
  ids=find(strcmp(d.ranking.table.status,'feasible'));cs=d.cands(ids);
  B=pipecompare.eval.Rank.bootstrap(cs,d.ref,opts,{d.ranking.objective});actual=B{1};
  nC=numel(d.ref.ids);W=cell(1,nC);indices=W;rs=RandStream('mt19937ar','Seed',opts.seed);
  segmented=pipecompare.utils.fieldOr(d.ref,'segmented',false);assert(~segmented);
  for c=1:nC,[W{c},indices{c}]=pipecompare.eval.Rank.resampleWeights(rs,d.ref.n(c),1999,segmented);end
  rebuilt=zeros(size(actual));candidateRows={};
  for k=1:numel(cs)
   ob=cs(k).m.objectives(1);assert(strcmp(ob.kind,'scalar'));ob.segmented=false;
   gain=1;if isfield(cs(k).signal,'gain'),gain=cs(k).signal.gain;end
   spread=0;if isfield(cs(k).signal,'shapeSpread'),spread=cs(k).signal.shapeSpread;if ~isfinite(spread),spread=0;end;end
   assert(isscalar(gain)&&isfinite(gain)&&gain>0);
   sm=pipecompare.eval.Measure.smeBoot(ob,W,opts,7919,ones(1,nC),indices)/gain;
   b2=zeros(1,1999);conditionRows={};
   for c=1:nC
    x=ob.X{c}(:,1);present=~isnan(x);z=x;z(~present)=0;
    wt=ones(size(x));if ~isempty(ob.w),wt=ob.w{c};wt(isnan(wt))=0;end
    K=W{c}.*wt;K=K.*double(present);
    blink=logical(d.ref.blink{1}{c}(:))&present;clean=~logical(d.ref.blink{1}{c}(:))&present;
    nb=double(blink)'*K;nn=double(clean)'*K;n=nb+nn;drawn=double(present)'*W{c};
    [~,eq,assessed]=pipecompare.eval.Measure.eyeBias(x,d.ref.blink{1}{c},W{c}.*wt);
    one=ob;one.X=ob.X(c);one.w={};if ~isempty(ob.w),one.w=ob.w(c);end
    smc=pipecompare.eval.Measure.smeBoot(one,W(c),opts,7919,1,indices(c));
    sb=(z.*double(blink))'*K;sn=(z.*double(clean))'*K;s2=(z.^2)'*K;
    residual=s2-sb.^2./nb-sn.^2./nn;
    % Original formula, diagnostic comparison only (not used for selection).
    mu=(sb+sn)./n;vo=max(s2./n-mu.^2,0).*n./max(n-1,eps);
    bias=nb./n.*(sb./nb-sn./nn);oldEq=max(abs(bias)-2*(nb./n).*sqrt(vo.*(1./nb+1./nn)),0).^2;
    oldEq(nb==0&nn>0)=0;oldEq(~(nn>0))=NaN;
    est=ob.estimate(c);if ~isfinite(est),est=0;end
    b2=b2+((sqrt(eq)+abs(est)*spread)/gain).^2;
    bad=find(~isfinite(eq)|~isfinite(smc));causes=cell(1,numel(bad));
    for j=1:numel(bad)
     t=bad(j);
     if any(~isfinite(K(:,t))),cause='nonfinite_weights';
     elseif drawn(t)<2,cause='fewer_than_two_retained_draws';
     elseif nn(t)==0,cause='no_nonblink_comparison';
     elseif nb(t)>0&&n(t)<=2,cause='pooled_weight_mass_at_most_two';
     elseif nb(t)==0,cause='no_blink_group_unexpected_nonfinite';
     else,cause='unexplained_numerical_or_other';end
     causes{j}=cause;
     positive=K(:,t)>0;
     events(end+1,:)={subject,seed,ids(k),c,t,cause,drawn(t),nb(t),nn(t),n(t), ...
      nnz(positive&blink),nnz(positive&clean),sum(K(:,t))^2/sum(K(:,t).^2), ...
      eq(t),smc(t),oldEq(t),assessed(t),residual(t),mat2str(W{c}(:,t)'),mat2str(wt')}; %#ok<AGROW>
    end
    validClean=nnz(clean&isfinite(wt)&wt>0);
    conditionRows{c}=struct('condition',c,'referenceTrials',d.ref.n(c),'retainedTrials',nnz(present), ...
     'positiveNonblinkTrials',validClean,'positiveBlinkTrials',nnz(blink&isfinite(wt)&wt>0), ...
     'probabilityNoNonblinkInIIDDraw',(1-validClean/d.ref.n(c))^d.ref.n(c), ...
     'badEyeCount',nnz(~isfinite(eq)),'badSMECount',nnz(~isfinite(smc)), ...
     'oldFormulaBadEyeCount',nnz(~isfinite(oldEq)),'badReplicates',bad,'causes',{causes});
   end
   rebuilt(k,:)=sqrt(sm.^2+b2/nC);
   candidateRows{k}=struct('id',ids(k),'pipeline',d.labels{ids(k)},'badScoreCount',nnz(~isfinite(actual(k,:))), ...
    'badScoreReplicates',find(~isfinite(actual(k,:))),'conditions',{conditionRows});
  end
  assert(isequal(isfinite(actual),isfinite(rebuilt)),'Nonfinite mask mismatch');
  ok=isfinite(actual);assert(max([0;abs(actual(ok)-rebuilt(ok))])<1e-10,'Score reconstruction mismatch');
  old=find(strcmp(d.labels,base.erp.pipeline));ib=find(ids==old);ic=find(ids==d.ranking.recommended);
  pi=find([prior.subject]==subject&[prior.seed]==seed);assert(isscalar(pi));
  invalid=any(~isfinite(actual(:)))||any(~isfinite(d.ranking.table.objective(ids)));
  assert(invalid==strcmp(prior(pi).gate.reason,'nonfinite or insufficient bootstrap evidence'));
  q=struct('subject',subject,'seed',seed,'candidates',numel(ids),'bootstrapColumns',1999, ...
   'badScoreCells',nnz(~isfinite(actual)),'badBootstrapColumns',nnz(any(~isfinite(actual),1)), ...
   'incumbentId',old,'challengerId',d.ranking.recommended,'incumbentBadCount',nnz(~isfinite(actual(ib,:))), ...
   'challengerBadCount',nnz(~isfinite(actual(ic,:))),'matchesPreviousGate',true,'candidateDiagnostics',{candidateRows});
  summaries{end+1}=q;brief=rmfield(q,'candidateDiagnostics');fprintf('NONFINITE_DIAGNOSIS %s\n',jsonencode(brief));
 end
end
T=cell2table(events,'VariableNames',{'subject','seed','candidate','condition','replicate','cause','retainedDraws','blinkMass','nonblinkMass','totalMass','uniqueBlink','uniqueNonblink','kishEffectiveN','eyeSq','sme','oldFormulaEyeSq','assessed','pooledResidual','drawCounts','analyticWeights'});
writetable(T,fullfile(outDir,'invalid_draws.csv'));
fid=fopen(fullfile(outDir,'diagnosis.json'),'w');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(summaries));fclose(fid);
fprintf('NONFINITE_RECONSTRUCTION_CHECKS_PASSED records=%d invalidConditionRows=%d\n',numel(summaries),height(T));
end
