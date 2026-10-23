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

function testWeightedRankingAndBootstrapIgnoreWeightNormalization(tc)
x=[linspace(-1,1,20)'; 10+linspace(-1,1,20)']; b=[false(20,1);true(20,1)];
c=nqc_cand(x); ref=nqc_ref(40); ref.blink={{b}};
opts=pipecompare.eval.Rank.defaults(); opts.nBoot=99; opts.nSplits=0;
w=linspace(0.2,1,40)';
for scale=[1 100]
    o=c.m.objectives; o.w={scale*w};
    [o.estimate,o.sme]=pipecompare.eval.Measure.estimate(o,x,opts,0,o.w{1});
    [o.eyeBias,o.eyeBiasSq]=pipecompare.eval.Measure.eyeBias(x,b,o.w{1});
    c.m.objectives=o;
    R=pipecompare.eval.Rank.run(c,ref,opts); B=pipecompare.eval.Rank.bootstrap(c,ref,opts);
    if scale==1, point=R.table.objective; boot=B{1};
    else
        verifyEqual(tc,R.table.objective,point,'AbsTol',1e-12);
        verifyEqual(tc,B{1},boot,'AbsTol',1e-12);
    end
end
end

function testGroupRankPassesPairedCovarianceToOwnChoices(tc)
rs=RandStream('mt19937ar','Seed',3007); x=randn(rs,60,1);
Y=[1 1.05; 1.08 1; 1.01 1.02]; opts=pipecompare.eval.Rank.defaults();
opts.nBoot=199; opts.nSplits=0;
S=struct('name',{},'keys',{},'labels',{},'table',{},'cands',{},'ref',{},'ranking',{});
for i=1:size(Y,1)
    c=[nqc_cand(x*Y(i,1)),nqc_cand(x*Y(i,2))]; ref=nqc_ref(60);
    R=pipecompare.eval.Rank.run(c,ref,opts);
    S(i)=struct('name',sprintf('record%d',i),'keys',{{'a','b'}},'labels',{{'a','b'}}, ...
        'table',R.table,'cands',c,'ref',ref,'ranking',struct('options',opts,'objective',R.objective));
end
G=pipecompare.Group.rank(S);
verifyEqual(tc,G.own.id,[1;2;1]);
verifyTrue(tc,G.twoLevelModel.pairedCovariance);
verifyEqual(tc,G.ownShrunk./G.ownShrunk(:,1),Y./Y(:,1),'AbsTol',1e-9);
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

function testFrozenGroupPriorDoesNotPoolNewParticipants(tc)
rs=RandStream('mt19937ar','Seed',6701);
Y=exp(0.2*randn(rs,12,3)+[0 0.1 -0.1]); V=(0.1*Y).^2;
[~,~,~,prior]=pipecompare.Group.twoLevel(Y,V);
x=[1.2 0.9 1.1]; v=(0.08*x).^2;
[a,info]=pipecompare.Group.applyPrior(prior,x,v);
b=pipecompare.Group.applyPrior(prior,[x;1000 1 1e-3],[v;ones(1,3)]);
verifyEqual(tc,b(1,:),a,'AbsTol',1e-12);
verifyFalse(tc,info.priorRefitted); verifyTrue(tc,info.assessed);
verifyFalse(tc,info.fullyHeldOut);
scaled=pipecompare.Group.applyPrior(prior,10*x,100*v);
verifyEqual(tc,scaled,10*a,'AbsTol',1e-10);
end

function testFrozenGroupPriorPreservesExactPairedContrasts(tc)
Y=[1 1.1;1.3 1.2;0.9 1.15;1.2 1.1]; V=(0.2*Y).^2;
C=arrayfun(@(s) 0.2^2*(Y(s,:)'*Y(s,:)),(1:4)','UniformOutput',false);
[~,~,~,prior]=pipecompare.Group.twoLevel(Y,V,C);
x=[1.05 1.02]; v=(0.2*x).^2;
[a,info]=pipecompare.Group.applyPrior(prior,x,v,{0.2^2*(x'*x)});
verifyTrue(tc,info.assessed);
verifyEqual(tc,a(2)/a(1),x(2)/x(1),'AbsTol',1e-9);
end

function testFrozenGroupPriorRejectsDisconnectedScales(tc)
Y=nan(4);Y(1:2,1:2)=[1 2;2 1];Y(3:4,3:4)=[1 2;2 1];
[~,~,~,prior]=pipecompare.Group.twoLevel(Y,0.01*ones(4));
x=[1 2 3 4];[a,info]=pipecompare.Group.applyPrior(prior,x,0.01*ones(1,4));
verifyEqual(tc,a,x);verifyFalse(tc,info.assessed);verifyNotEmpty(tc,info.reason{1});
end

function testMultitaperRejectsAliasedBand(tc)
verifyError(tc,@()pipecompare.eval.Measure.bandScores(ones(1,256),256,[120 140]),'PipeCompare:Measure');
verifyError(tc,@()pipecompare.eval.Measure.bandScores(ones(1,256),256,[-1 10]),'PipeCompare:Measure');
end

function testBadStretchMissingSamplesAreUnassessed(tc)
rs=RandStream('mt19937ar','Seed',6702);X=randn(rs,8,3000);X(1,100)=NaN;
[st,check]=pipecompare.run.Steps.badStretches(X,50);
verifyEmpty(tc,st);verifyFalse(tc,check.assessed);verifyEqual(tc,check.reason,'nonfinite input samples');
end

function testHmmKeepsFiniteEventsWhenASampleIsMissing(tc)
rs=RandStream('mt19937ar','Seed',317); fs=250; t=(0:30*fs-1)/fs;
X=0.2*randn(rs,size(t));
for at=5:5:25, X=X+12*exp(-0.5*((t-at)/0.08).^2); end
X(100)=NaN;
[hit,m]=pipecompare.eval.Hmm.events(X,fs);
verifyTrue(tc,m.valid); verifyGreaterThan(tc,nnz(hit),100);
verifyFalse(tc,hit(100)); verifyEqual(tc,m.inputMissingSamples,1);
[hit,m]=pipecompare.eval.Hmm.segment([0:8 nan(1,100)],3.5,3);
verifyFalse(tc,m.valid); verifyFalse(tc,any(hit));
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
verifyEqual(tc,C.rawEstimate,o.estimate(:),'AbsTol',1e-12);
verifyEqual(tc,C.correctedEstimate,o.estimate(:)/g,'AbsTol',1e-12);
verifyEqual(tc,C.rawError,o.estimate(:)-[8;28],'AbsTol',1e-12);
verifyEqual(tc,C.actualError,C.correctedEstimate-C.populationTruth,'AbsTol',1e-12);
verifyNotEqual(tc,C.rawError,C.actualError);
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

function testCalibrationRetainsInfiniteRatioDraws(tc)
% Two equal-error recordings, one reporting zero: a quarter of
% two-record resamples have an infinite ratio. The upper CI must be Inf.
C=table([1;2],ones(2,1),[0;1], ...
    'VariableNames',{'recording','actualError','reportedError'});
S=pipecompare.eval.Calibration.summarize(C,1999,4);
verifyEqual(tc,S.ratio,sqrt(2),'AbsTol',1e-12);
verifyEqual(tc,S.ci(1),1,'AbsTol',1e-12);verifyEqual(tc,S.ci(2),Inf);
verifyEqual(tc,S.bootstrapDraws,1999);
verifyGreaterThan(tc,S.bootstrapInfiniteShare,0.15);
verifyLessThan(tc,S.bootstrapInfiniteShare,0.35);
verifyEqual(tc,S.bootstrapUndefinedShare,0);
end

function testCalibrationRetainsUndefinedRatioMass(tc)
% Zero/zero resamples convey no ratio. Returning [1,1] after deleting
% them would certify precision that these two recordings cannot supply.
C=table([1;2],[0;1],[0;1], ...
    'VariableNames',{'recording','actualError','reportedError'});
S=pipecompare.eval.Calibration.summarize(C,1999,4);
verifyEqual(tc,S.ratio,1);verifyEqual(tc,S.ci,[0 Inf]);
verifyGreaterThan(tc,S.bootstrapUndefinedShare,0.15);
verifyEqual(tc,S.bootstrapInfiniteShare,0);
end

function testCalibrationZeroReportingHasExplicitPointResult(tc)
C=table([1;2],[1;2],zeros(2,1), ...
    'VariableNames',{'recording','actualError','reportedError'});
S=pipecompare.eval.Calibration.summarize(C);
verifyEqual(tc,S.ratio,Inf);verifyEqual(tc,S.actualRms,sqrt(2.5),'AbsTol',1e-12);
verifyEqual(tc,S.reportedRms,0);
verifyTrue(tc,contains(S.interpretation,'zero error'));
C.actualError(:)=0;S=pipecompare.eval.Calibration.summarize(C);
verifyTrue(tc,isnan(S.ratio));verifyEqual(tc,S.actualRms,0);
verifyEqual(tc,S.reportedRms,0);
end

function testCrossfitStatesItsConditionalScope(tc)
c = nqc_cand((1:40)'); opts=pipecompare.eval.Rank.defaults(); opts.nSplits=2;
cf=pipecompare.eval.Rank.crossfit(c,nqc_ref(40),opts,'composite',c.m.composite);
verifyEqual(tc,cf.scope,'conditional-selection'); verifyFalse(tc,cf.fullyHeldOut);
verifyTrue(tc,any(strcmp(cf.fixedDecisions,'factorial shrinkage')));
opts.nBoot=99; R=pipecompare.eval.Rank.run(c,nqc_ref(40),opts);
verifyEqual(tc,R.validationScope,'conditional-on-fitted-candidates'); verifyFalse(tc,R.fullyHeldOut);
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

function testFrozenGroupPriorHandlesNoAssessedTrainingRecordings(tc)
[~,~,~,prior]=pipecompare.Group.twoLevel(nan(0,3),nan(0,3));
x=[1 NaN -1];v=[0.01 0.01 0.01];
[a,info]=pipecompare.Group.applyPrior(prior,x,v);
verifyEqual(tc,a(1),1);verifyTrue(tc,all(isnan(a(2:3))));
verifyFalse(tc,info.assessed);verifyEqual(tc,info.reason{1},'training interaction variance unassessed');
end

function testValidationScoresSurviveCheckpointAndGroupRanking(tc)
old=path;restore=onCleanup(@()path(old)); %#ok<NASGU>
addpath(fullfile(fileparts(mfilename('fullpath')),'..','benchmarks'));
opts=pipecompare.eval.Rank.defaults();opts.nBoot=99;opts.nSplits=0;
x=sin((1:40)');ref=nqc_ref(40);S=pipecompare_record_summary();
for j=1:2
    c=[nqc_cand(j*x),nqc_cand(1.1*j*x)];
    c(1).id=1;c(1).key='a';c(2).id=2;c(2).key='b';
    rank=pipecompare.eval.Rank.run(c,ref,opts);
    r=struct('cands',c,'leaves',struct('key',{'a','b'}),'labels',{{'a','b'}},'ranking',rank,'ref',ref);
    S(end+1)=pipecompare_record_summary(r,sprintf('checkpoint%d',j)); %#ok<AGROW>
end
file=[tempname '.mat'];clean=onCleanup(@()delete(file)); %#ok<NASGU>
save(file,'S');D=load(file);G=pipecompare.Group.rank(D.S);
verifyEqual(tc,{D.S.name},{'checkpoint1','checkpoint2'});
verifySize(tc,G.perRecording,[2 2]);verifyTrue(tc,all(isfinite(G.perRecording(:))));
end

function testValidationScriptsHaveNoAnalyzerErrors(tc)
assumeTrue(tc,exist('codeIssues','file')>0||exist('codeIssues','class')==8,'Code Analyzer requires R2022b or newer.');
folder=fullfile(fileparts(mfilename('fullpath')),'..','benchmarks');
issues=codeIssues(folder,'CodeAnalyzerConfiguration','factory');
errors=issues.Issues(string(issues.Issues.Severity)=="error",:);
verifyEmpty(tc,errors,evalc('disp(errors)'));
end

function testMultitaperFractionalSampleRateKeepsNyquist(tc)
% These grids round the Nyquist frequency in opposite directions. The
% endpoint must be included once, using its FFT index rather than equality.
cases=[100.1 328; 160.1 416; 200.1 328; 250.003 264; 250.1 264];
for q=1:size(cases,1)
    fs=cases(q,1);n=cases(q,2);band=[fs/2-fs/n/10 fs/2];
    X=(-1).^(0:n-1);H=pipecompare.eval.Measure.tapers(n,fs,band);
    F=fft(X'.*H(:,1));
    expected=log10(abs(F(n/2+1))^2/(fs*sum(H(:,1).^2)));
    verifyEqual(tc,pipecompare.eval.Measure.bandScores(X,fs,band),expected,'AbsTol',1e-10);
end
end

function testMissingEyeChannelsCannotBecomeQuietEvidence(tc)
EEG=nqc_synth(struct('seconds',90,'nPerCond',20,'artifactTrials',0,'blinkRate',0));
eyes=pipecompare.eval.Injection.eyeChannels(EEG.chanlocs);
assert(~isempty(eyes));
c=nqc_contract();
quiet=EEG;quiet.data(eyes,:,:)=0;
refQuiet=pipecompare.eval.Measure.reference(quiet,c);
mQuiet=pipecompare.eval.Measure.candidate(quiet,c,refQuiet);
verifyTrue(tc,all(isfinite(mQuiet.objectives.eyeBiasSq)));
EEG.data(eyes,:,:)=NaN;
ref=pipecompare.eval.Measure.reference(EEG,c);
m=pipecompare.eval.Measure.candidate(EEG,c,ref);
verifyTrue(tc,all(isfinite(m.objectives.estimate)));
verifyFalse(tc,any(m.objectives.eyeBiasAssessed));
verifyTrue(tc,all(isnan(m.objectives.eyeBiasSq)));
candidate=nqc_cand(m.objectives.X{1});candidate.m=m;
opts=pipecompare.eval.Rank.defaults();opts.nBoot=99;opts.nSplits=0;
R=pipecompare.eval.Rank.run(candidate,ref,opts);
verifyEmpty(tc,R.recommended);
verifyTrue(tc,contains(R.table.reason{1},'artifact bias unassessed'));
end

function testOcularMissingWindowIsAssessedOnlyAfterItsTrialIsRemoved(tc)
EEG=nqc_synth(struct('seconds',90,'nPerCond',20,'artifactTrials',0,'blinkRate',0));
eyes=pipecompare.eval.Injection.eyeChannels(EEG.chanlocs);EEG.data(eyes,:,:)=0;
c=nqc_contract();
[~,EEG]=evalc('pop_epoch(EEG,c.allEvents(),c.epoch,''epochinfo'',''yes'')');
times=EEG.xmin+(0:EEG.pnts-1)/EEG.srate;
EEG.data(eyes,times>=0.3&times<=0.6,1)=NaN;
EEG.data(eyes,times>=0.3&times<=0.6,2:4)=100;
ref=pipecompare.eval.Measure.reference(EEG,c);
verifyTrue(tc,all(ref.blink{1}{1}(2:4)));
m=pipecompare.eval.Measure.candidate(EEG,c,ref);
verifyTrue(tc,all(isfinite(m.objectives.estimate)));
verifyFalse(tc,any(m.objectives.eyeBiasAssessed));
verifyTrue(tc,all(isnan(m.objectives.eyeBiasSq)));
[~,remaining]=evalc('pop_select(EEG,''notrial'',1)');
m=pipecompare.eval.Measure.candidate(remaining,c,ref);
verifyTrue(tc,all(m.objectives.eyeBiasAssessed));
verifyTrue(tc,all(isfinite(m.objectives.eyeBiasSq)));
end

function testGroupValidationUsesActualRecommendationInsteadOfRawMinimum(tc)
oldPath=path;cleanup=onCleanup(@()path(oldPath)); %#ok<NASGU>
addpath(fullfile(fileparts(mfilename('fullpath')),'..','benchmarks'));
C=table(21*ones(4,1),{'raw-minimum';'raw-minimum';'selected';'selected'},[1;2;1;2], ...
    [0.1;0.1;2;3],logical([0;0;1;1]),[0.1;0.1;2;3],ones(4,1), ...
    'VariableNames',{'recording','steps','condition','actualError','selected','rawEstimate','gain'});
S=struct('name','S021R01.mat','keys',{{'a','b'}},'labels',{{'raw-minimum','selected'}});
S.table=table(ones(2,1),'VariableNames',{'gain_P3'});
S.cands=struct('m',{struct('objectives',struct('estimate',[0.1 0.1])), ...
    struct('objectives',struct('estimate',[2 3]))});
[loss,contrast,inFamily]=validation_recommended_loss(C,S,{'a'});
verifyEqual(tc,loss,6.5);verifyEqual(tc,contrast,-1);verifyFalse(tc,inFamily);
[again,~,inFamily]=validation_recommended_loss(C,S,{'b'});
verifyEqual(tc,again,loss);verifyTrue(tc,inFamily);
C.selected(:)=false;
[missing,contrast,inFamily]=validation_recommended_loss(C,S,{'a','b'});
verifyTrue(tc,isnan(missing));verifyTrue(tc,isnan(contrast));verifyFalse(tc,inFamily);
end

function testUnassessedCovarianceCannotCreateAnApparentPriorChange(tc)
oldPath=path;cleanup=onCleanup(@()path(oldPath)); %#ok<NASGU>
addpath(fullfile(fileparts(mfilename('fullpath')),'..','benchmarks'));
Y=[1 1.2;1.1 1.1;0.9 1.3];V=0.01*ones(size(Y));
[~,~,~,prior]=pipecompare.Group.twoLevel(Y,V);
y=[1 2];v=[NaN 0.01];
[theta,~]=pipecompare.Group.applyPrior(prior,y,v);
[a,b,eligible]=validation_group_choices(y,v,theta);
verifyEqual(tc,a,2);verifyEqual(tc,b,2);verifyEqual(tc,eligible,2);
v(:)=NaN;
[theta,info]=pipecompare.Group.applyPrior(prior,y,v);
[a,b,eligible]=validation_group_choices(y,v,theta);
verifyEmpty(tc,eligible);verifyTrue(tc,isnan(a));verifyTrue(tc,isnan(b));
verifyFalse(tc,info.assessed);
end

function testMissingPairedReplicatesCannotCertifyEquivalence(tc)
[keep,eqv]=pipecompare.eval.Rank.confidenceSet([1;1],nan(2,99),1,0.02,0.05);
verifyTrue(tc,all(keep));verifyTrue(tc,eqv(1));verifyFalse(tc,eqv(2));
B=repmat([1;2],1,99);B(2,1:20)=NaN;
[keep,eqv]=pipecompare.eval.Rank.confidenceSet([1;2],B,1,0.02,0.05);
verifyTrue(tc,all(keep));verifyFalse(tc,eqv(2));
end

function testNestedEquivalencePreservesMissingTailMass(tc)
pts=[1;1.01];B=repmat(pts,1,99);
[chosen,nEq]=pipecompare.eval.Rank.fixedSequence(pts,B,1,2,0.02,0.05);
verifyEqual(tc,[chosen nEq],[2 2]);
B(2,1:20)=NaN;
[chosen,nEq]=pipecompare.eval.Rank.fixedSequence(pts,B,1,2,0.02,0.05);
verifyEqual(tc,[chosen nEq],[1 1]);
end

function testIncompleteBootstrapCannotFitFactorialVariance(tc)
rs=RandStream('mt19937ar','Seed',20261023);
pts=(1:5)';B=pts+0.1*randn(rs,5,99);levels=repmat({{'step','same'}},1,5);
F=pipecompare.eval.Rank.factorialModel(levels,B,pts);verifyNotEmpty(tc,F);
B(1,1)=NaN;
F=pipecompare.eval.Rank.factorialModel(levels,B,pts);verifyEmpty(tc,F);
verifyEqual(tc,pipecompare.eval.Rank.shrink(F,pts),pts);
end

function testBootstrapStabilityCountsUnavailableRecommendations(tc)
b=[false(2,1);true(18,1)];c=nqc_cand([zeros(2,1);10*ones(18,1)]);
o=c.m.objectives;[o.eyeBias,o.eyeBiasSq,o.eyeBiasAssessed]=pipecompare.eval.Measure.eyeBias(o.X{1},b,ones(20,1));
c.m.objectives=o;ref=nqc_ref(20);ref.blink={{b}};
opts=pipecompare.eval.Rank.defaults();opts.nBoot=99;opts.nSplits=0;
R=pipecompare.eval.Rank.run([c c],ref,opts);
missing=R.byStratum.bootstrapUnavailableShare;
verifyGreaterThan(tc,missing,0);verifyLessThan(tc,R.byStratum.bootstrapCompleteShare,1);
verifyEqual(tc,sum(R.table.chosenShare)+missing,1,'AbsTol',1e-12);
end
