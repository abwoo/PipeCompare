function tests=test_ica_independent_groups
tests=functiontests(localfunctions);
end
function [ids,obj,status,stratum,support,family,recommended,measures]=fixture()
ids=[3;5];
obj=[1.11 1.12 1.12 1.20 1.21 1.22;
     1.10 1.11 1.10 1.25 1.26 1.25];
status=repmat({'feasible'},2,6);
stratum=repmat({'ERP_CORE_P3'},2,6);
support=true(2,6);
family=repmat({'weighted'},2,6);
recommended=[5 5 5 5 5 5];
measures=cell(2,6);
for j=1:2
 for sd=1:6
  measures{j,sd}=struct('seed',sd-1,'candidateId',ids(j), ...
      'correctedDifferenceEstimate',j+.01*sd,'effectiveUnflagged',[8 20]);
 end
end
end
function testDisjointTripletsMayDisagree(t)
[ids,obj,status,stratum,support,family,recommended,measures]=fixture();
z=pipecompare.eval.IcaIndependentGroups.compare(ids,obj,status,stratum, ...
    support,family,recommended,measures);
verifyEqual(t,z.status,'different-across-disjoint-triplets');
verifyFalse(t,z.sameCandidate);
verifyEqual(t,[z.groups{1}.candidateId,z.groups{2}.candidateId],[5 3]);
verifyEqual(t,z.groups{1}.seedNumbers,[0 1 2]);
verifyEqual(t,z.groups{2}.seedNumbers,[3 4 5]);
verifyFalse(t,z.isIndependentHeldoutStudy);
end
function testSameCandidateDoesNotImplyHoldout(t)
[ids,obj,status,stratum,support,family,recommended,measures]=fixture();
obj(1,4:6)=[1.3 1.31 1.32];
z=pipecompare.eval.IcaIndependentGroups.compare(ids,obj,status,stratum, ...
    support,family,recommended,measures);
verifyTrue(t,z.sameCandidate);
verifyEqual(t,z.status,'same-across-disjoint-triplets');
verifyEqual(t,z.groups{1}.candidateId,5);
verifyEqual(t,z.groups{2}.candidateId,5);
verifyFalse(t,z.isIndependentHeldoutStudy);
end
function testSecondGroupAbstainsOnSupport(t)
[ids,obj,status,stratum,support,family,recommended,measures]=fixture();
support(:,4:6)=false;
z=pipecompare.eval.IcaIndependentGroups.compare(ids,obj,status,stratum, ...
    support,family,recommended,measures);
verifyEqual(t,z.status,'at-least-one-group-abstained');
verifyEqual(t,z.groups{2}.status,'abstain');
verifyFalse(t,z.sameCandidate);
end
function testIncompleteSixSeedsRejected(t)
[ids,obj,status,stratum,support,family,recommended,measures]=fixture();
verifyError(t,@()pipecompare.eval.IcaIndependentGroups.compare(ids,obj(:,1:5), ...
    status(:,1:5),stratum(:,1:5),support(:,1:5),family(:,1:5), ...
    recommended(1:5),measures(:,1:5)),'PipeCompare:IcaIndependentGroups');
end
function testSeedMeasurementProvenanceMismatch(t)
[ids,obj,status,stratum,support,family,recommended,measures]=fixture();
measures{1,4}.seed=1;
verifyError(t,@()pipecompare.eval.IcaIndependentGroups.compare(ids,obj, ...
    status,stratum,support,family,recommended,measures),'PipeCompare:IcaIndependentGroups');
end
function testBaselineFamilySwitchInSecondGroup(t)
[ids,obj,status,stratum,support,family,recommended,measures]=fixture();
family(1,4:6)={'unweighted','unweighted','unweighted'};
z=pipecompare.eval.IcaIndependentGroups.compare(ids,obj,status,stratum, ...
    support,family,recommended,measures);
verifyEqual(t,z.groups{2}.candidateId,5);
verifyEqual(t,z.groups{2}.baselineFamily,'weighted');
end
function testDefaultSimpleModeUnchangedUnlessExplicitDevelopmentFlag(t)
prior=getenv('PR51_DEV_SIX_SEED_AUDIT');
guard=onCleanup(@()setenv('PR51_DEV_SIX_SEED_AUDIT',prior));
setenv('PR51_DEV_SIX_SEED_AUDIT','');
p=pipecompare.simple.Presets.searchOptions();
verifyEqual(t,p.icaSeeds,2);
setenv('PR51_DEV_SIX_SEED_AUDIT','1');
p=pipecompare.simple.Presets.searchOptions();
verifyEqual(t,p.icaSeeds,5);
clear guard;
end
