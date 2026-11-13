function tests = test_ica_postremoval
% Actual after-ICA-deletion Pz event-level recorder must not change EEG.
tests=functiontests(localfunctions);
end

function testRecordedPzAfterSyntheticDeletion(t)
[E,con]=fixture();
original=E;
z=pipecompare.research.IcaRemovalStage.capturePost(E,con);
verifyEqual(t,z.status,'complete');
verifyEqual(t,z.nEvents,6);
verifyEqual(t,z.eventIds,1:6);
verifyEqual(t,z.conditionIndex,[1 1 1 2 2 2]);
verifyEqual(t,z.observedPzUv,1:6,'AbsTol',1e-10);
verifyEqual(t,z.channel,'Pz');
verifyEqual(t,E,original);
% Artificially subtract a known IC's event-locked Pz signature;
% record REAL returned EEG, never infer it from the original source.
E.data(3,:,:)=E.data(3,:,:)*0.5;
z2=pipecompare.research.IcaRemovalStage.capturePost(E,con);
verifyEqual(t,z2.observedPzUv,0.5*(1:6),'AbsTol',1e-10);
end

function testUniqueEventProvenance(t)
[E,con]=fixture();
E.event(2).urevent=1;
verifyError(t,@()pipecompare.research.IcaRemovalStage.capturePost(E,con), ...
    'PipeCompare:IcaRemovalStage');
end

function testMissingPzFails(t)
[E,con]=fixture();
E.chanlocs(3).labels='CPz';
verifyError(t,@()pipecompare.research.IcaRemovalStage.capturePost(E,con), ...
    'PipeCompare:IcaRemovalStage');
end

function testNonfinitePostRemovalEEGFails(t)
[E,con]=fixture();
E.data(3,55,1)=NaN;
verifyError(t,@()pipecompare.research.IcaRemovalStage.capturePost(E,con), ...
    'PipeCompare:IcaRemovalStage');
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
E.event=repmat(struct('type','','urevent',0),1,6);
E.epoch=repmat(struct('event',[],'eventlatency',[],'eventtype',{{}}),1,6);
for k=1:6
 if k<=3,code='11';else,code='31';end
 E.event(k).type=code;E.event(k).urevent=k;
 E.epoch(k).event={k};E.epoch(k).eventlatency={0};E.epoch(k).eventtype={code};
end
end
