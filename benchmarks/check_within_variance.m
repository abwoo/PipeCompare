function check_within_variance()
% Known-truth and representation checks for the isolated variance change.
x=[0;1;-1;2;3;4;6;8];b=[false(4,1);true(4,1)];
[a,q]=pipecompare.eval.Measure.eyeBias(x,b,ones(size(x)));
p=mean(b);sv=((sum(b)-1)*var(x(b))+(sum(~b)-1)*var(x(~b)))/(numel(x)-2);
want=max(abs(a)-2*p*sqrt(sv*(1/sum(b)+1/sum(~b))),0)^2;
assert(abs(q-want)<1e-12);
[aa,qq]=pipecompare.eval.Measure.eyeBias(x+10,b,ones(size(x)));
assert(abs(a-aa)<1e-12&&abs(q-qq)<1e-12);
x=[zeros(20,1);3*ones(20,1)];b=[false(20,1);true(20,1)];
[a,q]=pipecompare.eval.Measure.eyeBias(x,b,ones(size(x)));assert(abs(a-1.5)<1e-12&&abs(q-2.25)<1e-12);
c=ones(40,1);c([2 23])=3;[a,q]=pipecompare.eval.Measure.eyeBias(x,b,c);
xx=repelem(x,c);bb=repelem(b,c);[aa,qq]=pipecompare.eval.Measure.eyeBias(xx,bb,ones(size(xx)));
assert(abs(a-aa)<1e-12&&abs(q-qq)<1e-12);
[~,q]=pipecompare.eval.Measure.eyeBias(x,false(size(b)),ones(size(x)));assert(q==0);
[~,q]=pipecompare.eval.Measure.eyeBias(x,true(size(b)),ones(size(x)));assert(isnan(q));
rs=RandStream('mt19937ar','Seed',3);q=zeros(1,2000);
for k=1:2000
 [~,q(k)]=pipecompare.eval.Measure.eyeBias(randn(rs,40,1),rand(rs,40,1)<.7,ones(40,1));
end
assert(mean(q)<.05/40,'Null penalty exceeds existing statistical tolerance');
s=testsuite('PipeCompare/tests/test_statistics.m');s=s(contains({s.Name},'testFactorialShrinkageChoosesBetter'));assert(numel(s)==1);res=run(s);assertSuccess(res);
fprintf('WITHIN_VARIANCE_CHECKS_PASSED meanNullPenalty=%.12g\n',mean(q));
end
