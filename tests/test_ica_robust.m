function tests=test_ica_robust
% Scientific gating is absolute: finite scores never rescue failed constraints.
tests=functiontests(localfunctions);
end
function [ids,obj,status,strata]=fixture()
ids=[4,25];
obj=[1.12 1.127 1.133;1.121 1.350 1.354];
status=repmat({'feasible'},2,3);
strata=repmat({'ERP_CORE_P3'},2,3);
end
function testCommonFeasibleMinimax(t)
[ids,obj,status,strata]=fixture();
z=pipecompare.eval.IcaRobust.select(ids,obj,status,strata);
verifyEqual(t,z.candidateId,4);verifyEqual(t,z.status,'selected');
verifyEqual(t,z.worstObjective,1.133,'AbsTol',1e-12);
verifyTrue(t,z.allSeedsScientificallyFeasible);
verifyTrue(t,z.originalRecommendationUnchanged);
end
function testFailedSeedExcludesSmallestObjective(t)
[ids,obj,status,strata]=fixture();
status{1,3}='rejected';
z=pipecompare.eval.IcaRobust.select(ids,obj,status,strata);
verifyEqual(t,z.candidateId,25);
verifyEqual(t,z.eligibleIds,25);
end
function testAllInfeasibleAbstains(t)
[ids,obj,status,strata]=fixture();
status(:,2)={'rejected';'failed'};
z=pipecompare.eval.IcaRobust.select(ids,obj,status,strata);
verifyTrue(t,isnan(z.candidateId));
verifyEqual(t,z.reason,'no-common-feasible-candidate');
end
function testNonfiniteCannotWin(t)
[ids,obj,status,strata]=fixture();
obj(1,2)=NaN;obj(2,3)=Inf;
z=pipecompare.eval.IcaRobust.select(ids,obj,status,strata);
verifyEqual(t,z.status,'abstain');
end
function testScientificStratumChangesAbstain(t)
[ids,obj,status,strata]=fixture();
strata{1,3}='incomparable';
z=pipecompare.eval.IcaRobust.select(ids,obj,status,strata);
verifyEqual(t,z.candidateId,25);
end
function testMultipleFeasibleStrataAbstain(t)
[ids,obj,status,strata]=fixture();
strata(2,:)={'different','different','different'};
z=pipecompare.eval.IcaRobust.select(ids,obj,status,strata);
verifyTrue(t,isnan(z.candidateId));
verifyEqual(t,z.reason,'multiple-feasible-strata-cannot-be-ranked-together');
end
function testNumericallyTiedAbstains(t)
[ids,obj,status,strata]=fixture();
obj(2,:)=[1.13,1.133,1.133];
z=pipecompare.eval.IcaRobust.select(ids,obj,status,strata);
verifyEqual(t,z.reason,'indistinguishable-minimax-tie-requires-predeclared-tiebreak');
verifyTrue(t,isnan(z.candidateId));
end
function testMismatchedDimensionsFail(t)
[ids,obj,status,strata]=fixture();
verifyError(t,@()pipecompare.eval.IcaRobust.select(ids,obj,status(:,1:2),strata),'PipeCompare:IcaRobust');
end
function testDuplicateIdsFail(t)
[~,obj,status,strata]=fixture();
verifyError(t,@()pipecompare.eval.IcaRobust.select([4,4],obj,status,strata),'PipeCompare:IcaRobust');
end
function testNegativeObjectiveCannotPass(t)
[ids,obj,status,strata]=fixture();
obj(1,1)=-0.1;
z=pipecompare.eval.IcaRobust.select(ids,obj,status,strata);
verifyEqual(t,z.candidateId,25);
end
