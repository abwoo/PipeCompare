function tests = test_ica_robust
% Scientific-gate-aware ICA minimax shadow method, isolated unit tests.
tests = functiontests(localfunctions);
end

function testAllSeedFeasibleChoice(testCase)
x=fixture();
y=pipecompare.eval.IcaRobust.preview(x);
verifyEqual(testCase,y.id,4);
verifyEqual(testCase,y.eligibleIds,[4 25]);
verifyFalse(testCase,y.candidateSelectedForProduction);
verifyTrue(testCase,y.existingRankRecommendationUnchanged);
end

function testOneSeedRejectedExcludesLowScore(testCase)
x=fixture();x.rankStatus{1,3}='rejected';
y=pipecompare.eval.IcaRobust.preview(x);
verifyEqual(testCase,y.id,25);
verifyEqual(testCase,y.eligibleIds,25);
end

function testMissingEffectiveTrialCannotBePromoted(testCase)
x=fixture();
x.candidateMeasurements{1,2}.effectiveUnflagged=[4 20];
x.candidateMeasurements{2,2}.effectiveUnflagged=[3 20];
y=pipecompare.eval.IcaRobust.preview(x);
verifyTrue(testCase,isnan(y.id));
verifyEmpty(testCase,y.eligibleIds);
verifyEqual(testCase,y.reason,'no-common-feasible-candidate');
end

function testCrossStratumCandidateNotComparable(testCase)
x=fixture();x.strata{1,2}='different';
y=pipecompare.eval.IcaRobust.preview(x);
verifyEqual(testCase,y.id,25);
end

function testTiesMustAbstain(testCase)
x=fixture();x.objective=ones(2,3);
y=pipecompare.eval.IcaRobust.preview(x);
verifyTrue(testCase,isnan(y.id));
verifyEqual(testCase,y.reason,'minimax-tie-abstain');
end

function testBrokenCandidateAlignmentFailsClosed(testCase)
x=fixture();x.candidateMeasurements{1,2}.candidateId=17;
y=pipecompare.eval.IcaRobust.preview(x);
verifyEqual(testCase,y.id,25);
end

function testMissingCalibrationDoesNotPass(testCase)
x=fixture();x.candidateMeasurements{1,1}.gain=NaN;
y=pipecompare.eval.IcaRobust.preview(x);
verifyEqual(testCase,y.id,25);
end

function testNoUnintendedPublicRankChange(testCase)
x=fixture();before=x.recommended;
pipecompare.eval.IcaRobust.preview(x);
verifyEqual(testCase,x.recommended,before);
end

function x=fixture()
ids=[4 25];scores=[1.12 1.13 1.14;1.10 1.35 1.36];
status=repmat({'feasible'},2,3);
strata=repmat({'reref:average'},2,3);
ms=cell(2,3);
for i=1:2
    for seed=0:2
        ms{i,seed+1}=struct('seed',seed,'candidateId',ids(i), ...
          'hasMeasurement',true,'gain',1, ...
          'correctedDifferenceEstimate',2+seed/100, ...
          'effectiveUnflagged',[10 25]);
    end
end
x=struct('ids',ids,'objective',scores,'rankStatus',{status}, ...
         'strata',{strata},'candidateMeasurements',{ms}, ...
         'recommended',[25 4 4]);
end
