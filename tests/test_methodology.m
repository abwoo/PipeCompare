function tests = test_methodology
tests = functiontests(localfunctions);
end

function setupOnce(~)
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
end

function testArtifactBiasIsInvariantToWeightScale(tc)
x = [zeros(20,1); 10*ones(20,1)]; b = [false(20,1); true(20,1)];
K = [ones(40,1), linspace(0.2,1,40)'];
[bias, sq] = pipecompare.eval.Measure.eyeBias(x,b,K);
[scaledBias, scaledSq] = pipecompare.eval.Measure.eyeBias(x,b,100*K);
verifyEqual(tc,scaledBias,bias,'AbsTol',1e-12);
verifyEqual(tc,scaledSq,sq,'AbsTol',1e-12);
end

function testArtifactBiasBootstrapCountsMatchExpandedSample(tc)
x = [0; 1; 2; 7; 9; 10]; b = logical([0;0;0;1;1;1]);
w = [0.4;1;0.7;0.5;0.8;1]; counts = [1;3;2;2;1;3];
idx = repelem((1:numel(x))',counts);
[bias,sq] = pipecompare.eval.Measure.eyeBias(x,b,counts.*w,counts);
[expandedBias,expandedSq] = pipecompare.eval.Measure.eyeBias(x(idx),b(idx),w(idx));
verifyEqual(tc,bias,expandedBias,'AbsTol',1e-12);
verifyEqual(tc,sq,expandedSq,'AbsTol',1e-12);
end

function testGroupVarianceRespectsDisconnectedDesign(tc)
A = [1 2; 2 1]; Y = nan(4); Y(1:2,1:2)=A; Y(3:4,3:4)=A;
[~,~,one] = pipecompare.Group.twoLevel(A,0.01*ones(2));
[~,~,two] = pipecompare.Group.twoLevel(Y,0.01*ones(4));
verifyEqual(tc,two,one,'AbsTol',1e-12);
end

function testGroupPairedCovariancePreservesKnownContrasts(tc)
rs = RandStream('mt19937ar','Seed',41);
Y = exp(0.3*randn(rs,20,1)+0.03*randn(rs,1,12)+0.08*randn(rs,20,12));
V = (0.3*Y).^2;
C = arrayfun(@(i) 0.3^2*(Y(i,:)'*Y(i,:)),(1:size(Y,1))','UniformOutput',false);
[theta,~,tau,model] = pipecompare.Group.twoLevel(Y,V,C);
% The only estimation error is a shared recording scale: pipeline
% contrasts are exact and must not be shrunk towards a common choice.
verifyGreaterThan(tc,tau,0);
verifyEqual(tc,log(theta)-mean(log(theta),2),log(Y)-mean(log(Y),2),'AbsTol',1e-9);
verifyTrue(tc,model.pairedCovariance);
end

function testGroupDiagonalCovarianceMatchesIndependentModel(tc)
rs = RandStream('mt19937ar','Seed',19); Y=exp(0.2*randn(rs,8,5)); V=(0.05*Y).^2;
C=arrayfun(@(i) diag(V(i,:)),(1:size(Y,1))','UniformOutput',false);
[a,sa,ta]=pipecompare.Group.twoLevel(Y,V); [b,sb,tb]=pipecompare.Group.twoLevel(Y,V,C);
verifyEqual(tc,b,a,'AbsTol',1e-12); verifyEqual(tc,[sb tb],[sa ta],'AbsTol',1e-12);
end

function testMultitaperEndpointsAreNotDoubled(tc)
n=256; fs=256;
for band = {[0 0.1],[127.9 128]}
    f=band{1}; X=ones(1,n); if f(1)>0, X=(-1).^(0:n-1); end
    H=pipecompare.eval.Measure.tapers(n,fs,f); F=fft(X'.*H(:,1));
    at=1; if f(1)>0, at=n/2+1; end
    expected=log10(abs(F(at))^2/fs);
    verifyEqual(tc,pipecompare.eval.Measure.bandScores(X,fs,f),expected,'AbsTol',1e-10);
end
end

function testHmmExtremeImpulseDoesNotMarkWholeRecording(tc)
rs = RandStream('mt19937ar', 'Seed', 93);
fs = 250; X = randn(rs, 1, 150000); X(75000) = 1e4;
hit = pipecompare.eval.Hmm.events(X, fs);
verifyLessThan(tc, mean(hit), 0.1);
[~, m] = pipecompare.eval.Hmm.segment([zeros(1,4999) 100 zeros(1,5000)], 3.5, 3);
verifyTrue(tc, m.valid);
verifyTrue(tc, all(isfinite([m.mu m.sd m.A(:)' m.pi])));
verifyEqual(tc, sum(m.A, 2), ones(2, 1), 'AbsTol', 1e-12);
end

function testHmmInvalidFitIsReported(tc)
[hit, m] = pipecompare.eval.Hmm.segment(nan(1,100), 3.5, 3);
verifyFalse(tc, any(hit)); verifyFalse(tc, m.valid); verifyNotEmpty(tc, m.reason);
[hit, m] = pipecompare.eval.Hmm.segment(zeros(1,100), 3.5, 3);
verifyFalse(tc, any(hit)); verifyTrue(tc, m.valid);
end

function testAllBlinkTrialsHaveUnassessedBias(tc)
[b, sq, ok] = pipecompare.eval.Measure.eyeBias((1:40)', true(40,1), ones(40,1));
verifyTrue(tc, isnan(b) && isnan(sq)); verifyFalse(tc, ok);
[b, sq, ok] = pipecompare.eval.Measure.eyeBias((1:40)', false(40,1), ones(40,1));
verifyEqual(tc, [b sq], [0 0]); verifyTrue(tc, ok);
end

function testGroupBootstrapStreamsAreIndependentAndStable(tc)
names = {'recording_a', 'recording_b', 'recording_c'};
seeds = pipecompare.Group.recordingSeeds(names, 1);
verifyEqual(tc, numel(unique(seeds)), 3);
p = [3 1 2]; reordered = pipecompare.Group.recordingSeeds(names(p), 1);
verifyEqual(tc, reordered, seeds(p));
verifyEqual(tc, numel(unique(pipecompare.Group.recordingSeeds({'same','same',''}, 1))), 3);
N = 60; c = nqc_cand(exp(linspace(-2,2,N))'); ref = nqc_ref(N);
o = pipecompare.eval.Rank.defaults(); o.nBoot = 999; o.seed = seeds(1);
b1 = pipecompare.eval.Rank.bootstrap(c, ref, o);
o.seed = seeds(2); b2 = pipecompare.eval.Rank.bootstrap(c, ref, o);
a=b1{1}-mean(b1{1}); b=b2{1}-mean(b2{1});
verifyLessThan(tc, abs((a*b')/(norm(a)*norm(b))), 0.15);
end

function testSegmentSmeMatchesCorrelatedReplications(tc)
rs = RandStream('mt19937ar', 'Seed', 785); N = 500; reps = 600;
for rho = [0 0.6]
    av = zeros(1,reps); se = av;
    for r = 1:reps
        x = filter(1,[1 -rho],randn(rs,N+500,1)); x = x(501:end);
        [av(r), se(r)] = pipecompare.eval.Measure.meanError(x, [], true);
    end
    verifyEqual(tc, sqrt(mean(se.^2)), std(av), 'RelTol', 0.15);
end
end

function testSegmentSmeKeepsMissingPositionsAndBootstrapOrder(tc)
x = [1; 2; NaN; 4; 5; 2; 3; 6; 8; 9; 3; 4];
[~, se] = pipecompare.eval.Measure.meanError(x, [], true);
o = struct('kind','scalar','X',{{x}},'segmented',true);
opts = pipecompare.eval.Rank.defaults(); opts.nBoot = 20;
s = RandStream('mt19937ar','Seed',12);
[W, I] = pipecompare.eval.Rank.resampleWeights(s,numel(x),20,true);
v = pipecompare.eval.Measure.smeBoot(o,{W},opts,0,1,{I});
for k = 1:20
    [~, expected] = pipecompare.eval.Measure.meanError(x(I(:,k)), [], true);
    verifyEqual(tc, v(k), expected, 'AbsTol', 1e-12);
end
verifyGreaterThan(tc, se, 0);
end

function testUndefinedArtifactBiasCannotPassRanking(tc)
c = nqc_cand((1:40)'); o = c.m.objectives;
[o.eyeBias,o.eyeBiasSq,o.eyeBiasAssessed] = pipecompare.eval.Measure.eyeBias(o.X{1},true(40,1),ones(40,1));
c.m.objectives = o;
opts = pipecompare.eval.Rank.defaults(); opts.nBoot=99; opts.nSplits=0;
R = pipecompare.eval.Rank.run(c,nqc_ref(40),opts);
verifyEmpty(tc,R.recommended);
verifyTrue(tc,contains(R.table.reason{1},'artifact bias unassessed'));
end

function testCalibrationUsesTheRankedErrorAndSeparateTruths(tc)
c = nqc_cand({(1:20)',(21:40)'}); o=c.m.objectives;
o.eyeBias=[2 3]; o.eyeBiasSq=[0.25 1]; c.m.objectives=o;
g=0.9; sp=0.1;
T=table(g,sp,'VariableNames',{'gain_P3','shapeSpread_P3'});
r=struct('cands',c,'labels',{{'example'}},'ranking',struct('table',T,'recommended',1));
C=pipecompare.eval.Calibration.rows(r,1,[8 28],[9 29],10);
[reported,~,~]=pipecompare.eval.Rank.conditionErrors(o,g,sp);
verifyEqual(tc,C.reportedError,reported(:),'AbsTol',1e-12);
verifyGreaterThan(tc,C.reportedError,sqrt(o.sme(:).^2+o.eyeBiasSq(:))/g);
verifyEqual(tc,C.actualError-C.conditionalError,ones(2,1),'AbsTol',1e-12);
verifyTrue(tc,all(C.selected));
E=pipecompare.eval.Calibration.rows(r,[],[8 28],[9 29],10); S=pipecompare.eval.Calibration.summarize(E);
verifyEqual(tc,height(E),0); verifyTrue(tc,isnan(S.ratio));
end

function testAncillaryUnknownBiasDoesNotBlockTheChosenObjective(tc)
x=(1:40)'; c=nqc_cand({{x},{x}}); ref=nqc_ref(40,{'obj1','obj2'},{'uV','uV'});
c.m.objectives(2).eyeBias=NaN; c.m.objectives(2).eyeBiasSq=NaN;
c.m.objectives(2).eyeBiasAssessed=false;
opts=pipecompare.eval.Rank.defaults(); opts.objective='obj1'; opts.nBoot=99; opts.nSplits=0;
R=pipecompare.eval.Rank.run(c,ref,opts);
verifyEqual(tc,R.recommended,1); verifyTrue(tc,isnan(R.table.eyeBias_obj2));
end

function testCalibrationRecordsAreTheSamplingUnits(tc)
C=table([1;2],[1;3],ones(2,1),'VariableNames',{'recording','actualError','reportedError'});
S=pipecompare.eval.Calibration.summarize(C,199,4);
% More correlated rows for the first recording must not change its weight.
D=[repmat(C(1,:),20,1);C(2,:)]; Q=pipecompare.eval.Calibration.summarize(D,199,4);
verifyEqual(tc,S.ratio,sqrt(5),'AbsTol',1e-12); verifyEqual(tc,Q.ratio,S.ratio,'AbsTol',1e-12);
verifyEqual(tc,Q.ci,S.ci,'AbsTol',1e-12); verifyEqual(tc,Q.nRecords,2);
verifyEqual(tc,S.interpretation,'reported error understated'); verifyTrue(tc,S.exploratory);
C.actualError=[0.1;0.3]; S=pipecompare.eval.Calibration.summarize(C);
verifyEqual(tc,S.interpretation,'reported error overstated');
C.reportedError(1)=NaN; S=pipecompare.eval.Calibration.summarize(C);
verifyEqual(tc,S.unassessedRecords,1); verifyEqual(tc,S.nRecords,1); verifyTrue(tc,all(isnan(S.ci)));
end

function testCrossfitStatesItsConditionalScope(tc)
c = nqc_cand((1:40)'); opts=pipecompare.eval.Rank.defaults(); opts.nSplits=2;
cf=pipecompare.eval.Rank.crossfit(c,nqc_ref(40),opts,'composite',c.m.composite);
verifyEqual(tc,cf.scope,'conditional-selection'); verifyFalse(tc,cf.fullyHeldOut);
verifyTrue(tc,any(strcmp(cf.fixedDecisions,'factorial shrinkage')));
end

function testBadStretchesReportUnsupportedRank(tc)
rs=RandStream('mt19937ar','Seed',8); X=randn(rs,8,40*50); X=X-mean(X,1);
[st, check]=pipecompare.run.Steps.badStretches(X,50);
verifyEmpty(tc,st); verifyFalse(tc,check.assessed); verifyTrue(tc,contains(check.reason,'rank-deficient'));
[st, check]=pipecompare.run.Steps.badStretches(X(1:4,:),50);
verifyEmpty(tc,st); verifyFalse(tc,check.assessed); verifyTrue(tc,contains(check.reason,'8'));
[st, check]=pipecompare.run.Steps.badStretches(X(:,1:500),50);
verifyEmpty(tc,st); verifyFalse(tc,check.assessed); verifyTrue(tc,contains(check.reason,'30'));
end

function testMultitaperReportsMeanPsdAndEffectiveResolution(tc)
[~, info]=pipecompare.eval.Measure.tapers(500,250,[9.9 10.1]);
verifyFalse(tc,info.resolved); verifyGreaterThan(tc,info.smoothingWidth,info.bandWidth);
verifyEqual(tc,info.unit,'log10(uV^2/Hz)');
if exist('dpss','file')
    for n=[250 500 2000]
        for NW=[1 1.5 2]
            [H, info]=pipecompare.eval.Measure.tapers(n,250,[8 8+2*NW*250/n]);
            D=dpss(n,NW,size(H,2)); verifyGreaterThan(tc,min(svd(H'*D)),1-1e-8);
            verifyTrue(tc,info.resolved);
        end
    end
end
end
