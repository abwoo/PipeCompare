function tests = test_ica_rank_lock
% Regression for rank changes caused by shuffling before 3000-sample EEG
% getrank, not a test of source reconstruction or ERP science.
tests=functiontests(localfunctions);
end

function testRankOfSameUnshuffledICAInputIsDeterministic(t)
n=6000;tt=1:n;
x=[sin(.04*tt);cos(.07*tt);zeros(1,n)];
x(3,3001:end)=sin(.14*tt(3001:end));
[r,z]=pipecompare.run.IcaRankLock.estimate(x,true);
verifyEqual(t,r,2);
verifyEqual(t,z.sampleCount,3000);
verifyTrue(t,z.computedBeforeSeedPermutation);
verifyTrue(t,z.productionSeedZeroUnchanged);
% The naive first-3000 approach sees a different, full-rank sample subset.
stream=RandStream('mt19937ar','Seed',11);
xs=x(:,randperm(stream,n));
[rShuffled,~]=pipecompare.run.IcaRankLock.estimate(xs,true);
verifyEqual(t,rShuffled,3);
% New Steps.ica always selects r=2 BEFORE applying any seed permutation.
end

function testFullRankWithoutArtificialDimensionalityReduction(t)
tt=1:3600;
x=[sin(.012*tt);cos(.028*tt);sin(.044*tt)];
[r,report]=pipecompare.run.IcaRankLock.estimate(x,false);
verifyEqual(t,r,3);
verifyEqual(t,report.nChan,3);
end

function testDuplicatedChannelCannotIncreaseRank(t)
tt=1:3500;
x=[sin(.04*tt);cos(.05*tt);sin(.04*tt)];
[r,~]=pipecompare.run.IcaRankLock.estimate(x,true);
verifyEqual(t,r,2);
end

function testNonfiniteInputsAbstain(t)
x=[1 2 3;1 2 NaN];
verifyError(t,@()pipecompare.run.IcaRankLock.estimate(x,true),'PipeCompare:IcaRankLock');
end

function testInvalidRankRuleAbstains(t)
x=[1 2 3;3 2 1];
verifyError(t,@()pipecompare.run.IcaRankLock.estimate(x,1),'PipeCompare:IcaRankLock');
end
