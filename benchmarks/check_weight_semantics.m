function check_weight_semantics()
x=[0;1;-1;2;3;4;6;8];b=[false(4,1);true(4,1)];w=[1;3;1;2;2;1;4;1];c=[1;2;0;3;1;0;2;1];
[a,q]=pipecompare.eval.Measure.eyeBias(x,b,w,ones(size(w)));
for scale=[.1 10]
 [aa,qq]=pipecompare.eval.Measure.eyeBias(x,b,scale*w,ones(size(w)));
 assert(abs(a-aa)<1e-12&&abs(q-qq)<1e-12,'Analytic scaling changed result');
end
[a,q]=pipecompare.eval.Measure.eyeBias(x,b,c.*w,c);
xx=repelem(x,c);bb=repelem(b,c);ww=repelem(w,c);
[aa,qq]=pipecompare.eval.Measure.eyeBias(xx,bb,ww,ones(size(ww)));
assert(abs(a-aa)<1e-12&&abs(q-qq)<1e-12,'Mixed frequency/weight expansion differs');
[a,q]=pipecompare.eval.Measure.eyeBias(x,b,c);
[aa,qq]=pipecompare.eval.Measure.eyeBias(xx,bb,ones(size(xx)));
assert(abs(a-aa)<1e-12&&abs(q-qq)<1e-12,'Frequency expansion differs');
s=testsuite('PipeCompare/tests/test_statistics.m');
s=s(contains({s.Name},'testWhatBlinksAddToTheMean')|contains({s.Name},'testFactorialShrinkageChoosesBetter'));
assert(numel(s)==2);res=run(s);assertSuccess(res);
fprintf('WEIGHT_SEMANTICS_CHECKS_PASSED\n');
end
