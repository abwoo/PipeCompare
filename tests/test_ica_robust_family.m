function tests=test_ica_robust_family
tests=functiontests(localfunctions);
end
function [ids,obj,status,strata,support,family]=fixture()
ids=[2,3,5];
obj=[0.916 0.911 0.912;0.950 0.948 0.947;0.952 0.949 0.947];
status=repmat({'feasible'},3,3);
strata=repmat({'ERP_CORE_P3'},3,3);
support=true(3,3);
family={'unweighted','unweighted','unweighted'; ...
        'weighted','weighted','weighted'; ...
        'weighted','weighted','weighted'};
end
function testPreservesBaselineWeightedFamily(t)
[ids,obj,status,strata,support,family]=fixture();
z=pipecompare.eval.IcaRobustFamily.select(ids,obj,status,strata,support,family,5);
verifyEqual(t,z.baselineFamily,'weighted');
verifyEqual(t,z.candidateId,3);
verifyEqual(t,z.includedCandidateIds,[3 5]);
verifyEqual(t,z.excludedCandidateIds,2);
verifyEqual(t,z.status,'selected');
verifyTrue(t,z.originalRecommendationUnchanged);
end
function testBaselineUnweightedFamily(t)
[ids,obj,status,strata,support,family]=fixture();
z=pipecompare.eval.IcaRobustFamily.select(ids,obj,status,strata,support,family,2);
verifyEqual(t,z.candidateId,2);
verifyEqual(t,z.includedCandidateIds,2);
end
function testSeedWeightedInconsistencyExcluded(t)
[ids,obj,status,strata,support,family]=fixture();
family{2,3}='unweighted';
z=pipecompare.eval.IcaRobustFamily.select(ids,obj,status,strata,support,family,5);
verifyEqual(t,z.candidateId,5);
verifyEqual(t,z.includedCandidateIds,5);
end
function testFailedScienceConstraintsCannotBeRescued(t)
[ids,obj,status,strata,support,family]=fixture();
status(2:3,2)={'rejected';'failed'};
z=pipecompare.eval.IcaRobustFamily.select(ids,obj,status,strata,support,family,5);
verifyEqual(t,z.status,'abstain');
verifyTrue(t,isnan(z.candidateId));
end
function testInsufficientEffectiveHMMUnflaggedTrialsCannotBeRescued(t)
[ids,obj,status,strata,support,family]=fixture();
support(2:3,3)=false;
z=pipecompare.eval.IcaRobustFamily.select(ids,obj,status,strata,support,family,5);
verifyEqual(t,z.status,'abstain');
end
function testMissingOriginalRecipeAbstains(t)
[ids,obj,status,strata,support,family]=fixture();
z=pipecompare.eval.IcaRobustFamily.select(ids,obj,status,strata,support,family,42);
verifyEqual(t,z.reason,'baseline-candidate-not-present');
end
function testUnknownFamilyFails(t)
[ids,obj,status,strata,support,family]=fixture();
family{1,2}='other';
verifyError(t,@()pipecompare.eval.IcaRobustFamily.select(ids,obj,status,strata,support,family,5), ...
    'PipeCompare:IcaRobustFamily');
end
function testShapeMisalignmentFails(t)
[ids,obj,status,strata,support,family]=fixture();
verifyError(t,@()pipecompare.eval.IcaRobustFamily.select(ids,obj,status,strata,support,family(:,1:2),5), ...
    'PipeCompare:IcaRobustFamily');
end
