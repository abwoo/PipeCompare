function tests = test_ica_science
% Isolated read-only scientific recorder regression, no EEGLAB fitting.
tests = functiontests(localfunctions);
end

function testCapturesTruthfulICAAndWaveformWithoutMutation(t)
[E,c,ref,contract]=fixture();
beforeE=E;beforeC=c;beforeRef=ref;
z=pipecompare.research.IcaScience.capture(E,c,ref,contract);
verifyEqual(t,z.status,'complete');
verifyEqual(t,z.decompositions{1}.nComponents,4);
verifyEqual(t,z.decompositions{1}.removed,1);
verifyEqual(t,z.decompositions{1}.classNames,{'Brain','Eye'});
verifyEqual(t,z.decompositions{1}.mixing,eye(4));
verifyEqual(t,z.conditions{1}.retainedIds,[1 2 3]);
verifyEqual(t,z.conditions{2}.retainedIds,[4 5 6]);
verifyEqual(t,z.conditions{1}.windowMean,c.m.objectives(1).estimate(1),'AbsTol',1e-10);
verifyEqual(t,z.conditions{2}.windowMean,c.m.objectives(1).estimate(2),'AbsTol',1e-10);
verifyEqual(t,E,beforeE);
verifyEqual(t,c,beforeC);
verifyEqual(t,ref,beforeRef);
end

function testIncorrectScoreIsNotSilentlyRecorded(t)
[E,c,ref,contract]=fixture();
c.m.objectives(1).estimate(1)=999;
verifyError(t,@()pipecompare.research.IcaScience.capture(E,c,ref,contract), ...
    'PipeCompare:IcaScience');
end

function testMissingPositiveWeightedTrialFailsClosed(t)
[E,c,ref,contract]=fixture();
% Source reference ID exists but no surviving epoch.
ref.ids{1}=[1 2 3 101];
c.m.objectives(1).X{1}=[1;2;3;5];
c.m.objectives(1).estimate(1)=mean([1 2 3 5]);
verifyError(t,@()pipecompare.research.IcaScience.capture(E,c,ref,contract), ...
    'PipeCompare:IcaScience');
end

function testICLabelComponentDimensionMustMatch(t)
[E,c,ref,contract]=fixture();
c.steps{1}.iclabelProbabilities=ones(3,2);
verifyError(t,@()pipecompare.research.IcaScience.capture(E,c,ref,contract), ...
    'PipeCompare:IcaScience');
end

function [E,c,ref,contract] = fixture()
contract=pipecompare.eval.Contract('conditions',{'target',{'11'};'standard',{'31'}}, ...
 'epoch',[-0.2 0.8],'baseline',[-0.2 0],'components',{'P3',[0.3 0.6],{'Pz'}});
E=struct();
E.srate=100;E.trials=6;E.pnts=101;E.nbchan=4;
E.xmin=-0.2;E.xmax=0.8;
E.chanlocs=struct('labels',{'Fz','Cz','Pz','Oz'});
E.data=zeros(4,101,6);
tm=E.xmin+(0:E.pnts-1)/E.srate;
idx=tm>=0.3-1e-8 & tm<=0.6+1e-8;
for k=1:6
 E.data(3,idx,k)=k;
end
E.event=repmat(struct('type','','urevent',0),1,6);
E.epoch=repmat(struct('event',[],'eventlatency',[],'eventtype',{{}}),1,6);
for k=1:6
 if k<=3,code='11';else,code='31';end
 E.event(k).type=code;E.event(k).urevent=k;
 E.epoch(k).event={k};E.epoch(k).eventlatency={0};E.epoch(k).eventtype={code};
end
ref=struct('names',{{'target','standard'}}, ...
 'ids',{{[1;2;3],[4;5;6]}}, ...
 'blink',{{{false(3,1),false(3,1)}}});
o=struct('kind','scalar','X',{{(1:3)',(4:6)'}},'w',{{}}, ...
 'estimate',[2 5]);
c=struct('status','ok','m',struct('objectives',o), ...
 'ica',{{struct('icaweights',eye(4),'icasphere',eye(4), ...
 'icachansind',1:4)}}, ...
 'steps',{{struct('type','icremove','comps',1,'iclabelClasses', ...
 {{'Brain','Eye'}},'iclabelProbabilities',repmat([.8 .2],4,1))}});
end
