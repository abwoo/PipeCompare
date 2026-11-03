function diagnose_error_estimation(outDir)
% Fixed-label known-truth probes of unchanged production eyeBias.
% Alternative SE is a diagnostic, not a production replacement.
if ~isfolder(outDir),mkdir(outDir);end
rs=RandStream('mt19937ar','Seed',520113);B=4000;n=40;
b=[false(20,1);true(20,1)];cases=struct([]);
for weighted=[false true]
 w=ones(n,1);if weighted,w([1 21])=9;w=w/mean(w);end
 p=sum(w(b))/sum(w);nb=sum(w(b));nn=sum(w(~b));
 exactSE=p*sqrt(sum(w(b).^2)/nb^2+sum(w(~b).^2)/nn^2);
 for delta=[0 3 8]
  err=zeros(B,1);old=zeros(B,1);alt=zeros(B,1);qo=zeros(B,1);qa=zeros(B,1);
  for k=1:B
   x=randn(rs,n,1)+double(b)*delta;
   [bias,qo(k)]=pipecompare.eval.Measure.eyeBias(x,b,w);
   mu=sum(w.*x)/sum(w);v=sum(w.*(x-mu).^2)/(sum(w)-1);
   old(k)=p*sqrt(v*(1/nb+1/nn));
   assert(abs(qo(k)-max(abs(bias)-2*old(k),0)^2)<1e-10);
   mb=sum(w(b).*x(b))/nb;mn=sum(w(~b).*x(~b))/nn;
   df=sum(w)-sum(w(b).^2)/nb-sum(w(~b).^2)/nn;
   vw=(sum(w(b).*(x(b)-mb).^2)+sum(w(~b).*(x(~b)-mn).^2))/df;
   alt(k)=p*sqrt(vw*(sum(w(b).^2)/nb^2+sum(w(~b).^2)/nn^2));
   qa(k)=max(abs(bias)-2*alt(k),0)^2;err(k)=bias-p*delta;
  end
  z=struct('weighted',weighted,'delta',delta,'replications',B,'trueBias',p*delta,'exactSE',exactSE,'empiricalSE',std(err),'legacyRmsSE',sqrt(mean(old.^2)),'withinGroupWeightedRmsSE',sqrt(mean(alt.^2)),'legacyMeanPenaltySquared',mean(qo),'alternativeMeanPenaltySquared',mean(qa),'legacy95Coverage',mean(abs(err)<=1.96*old),'alternative95Coverage',mean(abs(err)<=1.96*alt));
  cases=[cases z]; %#ok<AGROW>
 end
end
x=[zeros(20,1);3*ones(20,1)];w=ones(40,1);scales=[.1 1 10];scaled=struct([]);
for a=scales
 [v,q]=pipecompare.eval.Measure.eyeBias(x,b,a*w);
 scaled=[scaled struct('scale',a,'bias',v,'penaltySquared',q)]; %#ok<AGROW>
end
counts=ones(40,1);counts([2 23])=3;
[bc,qc]=pipecompare.eval.Measure.eyeBias(x,b,counts);
xx=repelem(x,counts);bb=repelem(b,counts);[be,qe]=pipecompare.eval.Measure.eyeBias(xx,bb,ones(size(xx)));
assert(abs(bc-be)<1e-12&&abs(qc-qe)<1e-12,'Frequency representation mismatch');
r=struct('productionCommit','7151821edf8a5a54c3167d7f18b66c530911b568','executionCommit',getenv('GITHUB_SHA'),'run',getenv('GITHUB_RUN_ID'),'matlab',version,'cases',cases,'weightScaleProbe',scaled,'frequencyExpansionEquivalent',true,'productionModified',false,'scope','Fixed known blink labels, independent homoscedastic Gaussian scores, deterministic analytic weights. No causal separation of neural and artifact blink effects; no end-to-end improvement claim.');
fid=fopen(fullfile(outDir,'diagnosis.json'),'w');fprintf(fid,'%s\n',jsonencode(r));fclose(fid);
fprintf('PIPECOMPARE_ERROR_DIAGNOSIS %s\n',jsonencode(r));
end
