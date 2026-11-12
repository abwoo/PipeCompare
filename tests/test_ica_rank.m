function tests = test_ica_rank
% Regression for order-dependent first-3000-sample rank confound.
tests = functiontests(localfunctions);
end

function testFullRankDoesNotDependOnSampleOrdering(testCase)
X=fixture();
[r0,d0]=pipecompare.research.IcaRank.wholeData(X);
order=[3001:6000 1:3000];
[r1,d1]=pipecompare.research.IcaRank.wholeData(X(:,order));
verifyEqual(testCase,r0,16);
verifyEqual(testCase,r1,r0);
verifyEqual(testCase,d0.fitSampleCount,6000);
verifyEqual(testCase,d1.chosenPcaDimension,16);
end

function testOriginalWindowCanHaveLowerRank(testCase)
X=fixture();
% A factual example of the same first-3000 sampling bug that the
% unmodified runica/pop_runica path can encounter after permutation.
verifyEqual(testCase,rank(X(:,1:3000)),15);
verifyEqual(testCase,rank(X(:,3001:6000)),16);
[r,~]=pipecompare.research.IcaRank.wholeData(X);
verifyEqual(testCase,r,16);
end

function testNonfiniteDataFailsClosed(testCase)
X=fixture();X(2,5)=NaN;
verifyError(testCase,@()pipecompare.research.IcaRank.wholeData(X),'PipeCompare:IcaRank');
end

function testNoArtificialRankWhenDataDegenerate(testCase)
X=ones(17,4000);
verifyError(testCase,@()pipecompare.research.IcaRank.wholeData(X),'PipeCompare:IcaRank');
end

function X=fixture()
s=rng;cleanup=onCleanup(@()rng(s)); %#ok<NASGU>
rng(81501,'twister');
Y=randn(16,6000);
Y(16,1:3000)=Y(15,1:3000);
X=[Y;Y(1,:)+Y(2,:)];
end
