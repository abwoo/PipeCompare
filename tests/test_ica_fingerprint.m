function tests = test_ica_fingerprint
% Read-only, fail-closed ICA maps and component rejection summaries.
tests=functiontests(localfunctions);
end

function testRecordKnownTopologyAndRemovedIndex(t)
c=fixture();
z=pipecompare.eval.IcaFingerprint.snapshot(c);
verifyTrue(t,z.available);
verifyEqual(t,z.channels,[1 2 3]);
verifyEqual(t,z.maps,eye(3),'AbsTol',1e-12);
verifyEqual(t,z.removedIndices,2);
verifyEqual(t,z.removedCount,1);
verifyFalse(t,z.scientificPromotionAuthorized);
end

function testNothingRemovedIsExplicit(t)
c=fixture();
c.icsRemoved=0;
c.steps{2}.icsRemoved=0;
c.steps{2}.comps=[];
z=pipecompare.eval.IcaFingerprint.snapshot(c);
verifyTrue(t,z.available);
verifyEqual(t,z.removedCount,0);
verifyEmpty(t,z.removedIndices);
end

function testMissingScientificDecisionAbstains(t)
c=fixture();c.steps{2}.icsRemoved=NaN;
z=pipecompare.eval.IcaFingerprint.snapshot(c);
verifyFalse(t,z.available);
verifyEqual(t,z.reason,'component-rejection-unrecorded');
end

function testComponentCountMismatchAbstains(t)
c=fixture();c.icsRemoved=3;
z=pipecompare.eval.IcaFingerprint.snapshot(c);
verifyFalse(t,z.available);
verifyEqual(t,z.reason,'component-count-disagrees-with-pipeline');
end

function testComponentIndexOutOfRangeAbstains(t)
c=fixture();c.steps{2}.comps=4;
z=pipecompare.eval.IcaFingerprint.snapshot(c);
verifyFalse(t,z.available);
verifyEqual(t,z.reason,'component-rejection-inconsistent');
end

function testDuplicateChannelIndexAbstains(t)
c=fixture();c.ica{1}.icachansind=[1 1 3];
z=pipecompare.eval.IcaFingerprint.snapshot(c);
verifyFalse(t,z.available);
verifyEqual(t,z.reason,'ICA-matrix-or-channel-invalid');
end

function testMoreThanOneIcaFitDoesNotPretendComparable(t)
c=fixture();c.ica{2}=c.ica{1};
z=pipecompare.eval.IcaFingerprint.snapshot(c);
verifyFalse(t,z.available);
verifyEqual(t,z.reason,'multiple-ICA-fits-not-comparable');
end

function testOriginalInputsUntouched(t)
c=fixture();before=c;
pipecompare.eval.IcaFingerprint.snapshot(c);
verifyEqual(t,c,before);
end

function c=fixture()
d=struct('icaweights',eye(3),'icasphere',eye(3),'icachansind',[1 2 3]);
c=struct();
c.ica={d};
c.steps={struct('type','ica','icsRemoved',NaN,'comps',[]), ...
         struct('type','icremove','icsRemoved',1,'comps',2)};
c.icsRemoved=1;
end
