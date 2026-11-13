function tests = test_ica_contribution
% Pre-removal ICA sensor projection observer: no EEG/algorithm modification.
tests=functiontests(localfunctions);
end

function testFullRankKnownPzProjection(t)
[E,contract]=fixture();
prior=E;
out=pipecompare.research.IcaContribution.capturePreRemove(E,contract);
verifyEqual(t,out.status,'complete');
verifyEqual(t,out.nEvents,6);
verifyEqual(t,out.nComponents,3);
verifyEqual(t,out.eventIds,1:6);
verifyEqual(t,out.conditionIndex,[1 1 1 2 2 2]);
verifyEqual(t,out.componentUv(3,:),1:6,'AbsTol',1e-10);
verifyEqual(t,out.componentUv(1:2,:),zeros(2,6),'AbsTol',1e-10);
verifyEqual(t,out.maxAbsResidualUv,0,'AbsTol',1e-10);
verifyEqual(t,E,prior);
end

function testRankDeficientResidualCannotBeHidden(t)
[E,contract]=fixture();
E.icaweights=[1 0 0;0 1 1];
E.icawinv=pinv(E.icaweights*E.icasphere);
out=pipecompare.research.IcaContribution.capturePreRemove(E,contract);
verifyEqual(t,out.status,'complete');
verifyEqual(t,out.rank,2);
verifyGreaterThan(t,out.maxAbsResidualUv,1);
verifyEqual(t,out.observedPzUv,sum(out.componentUv,1)+out.residualPzUv,'AbsTol',1e-10);
end

function testNativeMixingAndPinvAreNotInterchangeable(t)
[E,contract]=fixture();
E.icawinv=diag([1,1,0.9]);
original=E;
out=pipecompare.research.IcaContribution.capturePreRemove(E,contract);
verifyEqual(t,out.status,'complete');
verifyGreaterThan(t,out.nativeVsPinvMixingMaxAbs,0.09);
verifyGreaterThan(t,out.nativeMaxAbsResidualUv,0.1);
verifyEqual(t,out.nativeICComponentUv(3,:),0.9*(1:6),'AbsTol',1e-10);
verifyEqual(t,E,original);
end

function testMissingNativeEEGLABMixingMustFailClosed(t)
[E,contract]=fixture();
E.icawinv=[];
verifyError(t,@()pipecompare.research.IcaContribution.capturePreRemove(E,contract), ...
 'PipeCompare:IcaContribution');
end

function testNonfiniteICAMustFailClosed(t)
[E,contract]=fixture();
E.icaweights(1,1)=NaN;
verifyError(t,@()pipecompare.research.IcaContribution.capturePreRemove(E,contract), ...
 'PipeCompare:IcaContribution');
end

function testMissingPzICMontageMustFailClosed(t)
[E,contract]=fixture();
E.icachansind=[1 2];
E.icaweights=eye(2);
E.icasphere=eye(2);
verifyError(t,@()pipecompare.research.IcaContribution.capturePreRemove(E,contract), ...
 'PipeCompare:IcaContribution');
end

function testDuplicateUreventMustFailClosed(t)
[E,contract]=fixture();
E.event(2).urevent=1;
verifyError(t,@()pipecompare.research.IcaContribution.capturePreRemove(E,contract), ...
 'PipeCompare:IcaContribution');
end

function [E,contract]=fixture()
contract=pipecompare.eval.Contract('conditions',{'target',{'11'};'standard',{'31'}}, ...
 'epoch',[-0.2 0.8],'baseline',[-0.2 0],'components',{'P3',[0.3 0.6],{'Pz'}});
E=struct();E.srate=100;E.nbchan=3;E.pnts=101;E.trials=6;
E.xmin=-0.2;E.xmax=0.8;
E.chanlocs=struct('labels',{'Fp1','Cz','Pz'});
E.data=zeros(3,101,6);
tm=E.xmin+(0:E.pnts-1)/E.srate;
win=tm>=0.3-1e-8 & tm<=0.6+1e-8;
for k=1:6,E.data(3,win,k)=k;end
E.icaweights=eye(3);E.icasphere=eye(3);E.icawinv=eye(3);E.icachansind=1:3;
E.event=repmat(struct('type','','urevent',0),1,6);
E.epoch=repmat(struct('event',[],'eventlatency',[],'eventtype',{{}}),1,6);
for k=1:6
 if k<=3,code='11';else,code='31';end
 E.event(k).type=code;E.event(k).urevent=k;
 E.epoch(k).event={k};E.epoch(k).eventlatency={0};E.epoch(k).eventtype={code};
end
end
