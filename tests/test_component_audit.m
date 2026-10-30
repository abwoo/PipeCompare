function tests=test_component_audit
tests=functiontests(localfunctions);
end
function setupOnce(~)
here=fileparts(mfilename('fullpath'));addpath(here);addpath(fullfile(here,'..'));addpath(fullfile(here,'..','benchmarks'));
end

function testTracingPreservesDefaultSimulationAndRandomStream(tc)
opts=struct('seconds',24,'nPerCond',4);
[plain,t0]=bench_synth(871301,opts);rng0=rng;
[traced,t1,sources]=bench_synth(871301,opts);rng1=rng;
verifyEqual(tc,traced.data,plain.data);verifyEqual(tc,traced.event,plain.event);
verifyEqual(tc,t0,t1);verifyEqual(tc,rng0,rng1);
names=fieldnames(sources);total=zeros(size(traced.data));
for k=1:numel(names),total=total+sources.(names{k});end
verifyEqual(tc,total,double(traced.data),'AbsTol',1e-10);
verifyGreaterThan(tc,norm(sources.p3,'fro'),0);verifyGreaterThan(tc,norm(sources.badChannel,'fro'),0);
end

function testFixedReplayNeverFitsIcaOrReselectsComponents(tc)
E=eeg_emptyset();E.data=single([1:100;101:200;201:300]);E.nbchan=3;E.pnts=100;E.trials=1;E.srate=100;
E.chanlocs=struct('labels',{'A','B','C'});E=eeg_checkset(E);
ica=struct('icaweights',eye(3),'icasphere',eye(3),'icachansind',1:3);
path={struct('type','ica','params',struct()),struct('type','icremove','params',struct())};
facts={struct('type','ica','params',struct()),struct('type','icremove','params',struct(),'comps',2)};
c=struct('steps',{facts},'ica',{{ica}});
before=rng;[got,audit]=pipecompare_fixed_replay(E,path,c,pipecompare.eval.Contract(),struct(),struct());after=rng;
verifyEqual(tc,double(got.data([1 3],:)),double(E.data([1 3],:)),'AbsTol',1e-10);
verifyEqual(tc,double(got.data(2,:)),zeros(1,100),'AbsTol',1e-10);
verifyEqual(tc,before,after);verifyEqual(tc,audit.icaFits,0);verifyFalse(tc,audit.preprocessingRefitted);
end

function testFixedReplayRejectsIncompleteOrUnsupportedDecisions(tc)
path={struct('type','ica','params',struct())};facts={struct('type','ica','params',struct())};
c=struct('steps',{facts},'ica',{{}});
verifyError(tc,@()pipecompare_fixed_replay(struct(),path,c,pipecompare.eval.Contract(),struct(),struct()),'PipeCompare:FixedReplay');
path={struct('type','native','params',struct())};c.steps={struct('type','native','params',struct())};
verifyError(tc,@()pipecompare_fixed_replay(struct(),path,c,pipecompare.eval.Contract(),struct(),struct()),'PipeCompare:FixedReplay');
end
