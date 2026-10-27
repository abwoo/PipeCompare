function tests=test_algorithm_regressions
tests=functiontests(localfunctions);
end

function setupOnce(~)
addpath(fileparts(mfilename('fullpath')));
addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
addpath(fullfile(fileparts(mfilename('fullpath')),'..','benchmarks'));
end

function testUnequalVarianceMatchesIndependentWelchCalculation(tc)
x=[-2;0;2;11;15;19;23];b=logical([0;0;0;1;1;1;1]);n=numel(x);p=mean(b);
[effect,~,ok,se]=pipecompare.eval.Measure.eyeBias(x,b,ones(n,1),ones(n,1),'heteroscedastic');
delta=mean(x(b))-mean(x(~b));
expected=p^2*(var(x(b))/nnz(b)+var(x(~b))/nnz(~b))+delta^2*p*(1-p)/(n-1);
verifyTrue(tc,ok);verifyEqual(tc,effect,p*delta,'AbsTol',1e-12);
verifyEqual(tc,se^2,expected,'AbsTol',1e-12);
end

function testUnequalVarianceTracksIndependentSampling(tc)
% Independent records, not candidate rows: common-variance pooling fails
% with rare, much noisier artifact trials. No benchmark seeds are used.
rs=RandStream('mt19937ar','Seed',871226);n=240;reps=1800;
for prevalence=[0.1 0.9]
    effect=zeros(reps,1);robust=effect;pooled=effect;
    for j=1:reps
        b=rand(rs,n,1)<prevalence;
        x=randn(rs,n,1).*(1+19*b)+5*b;
        [effect(j),~,ok,se]=pipecompare.eval.Measure.eyeBias(x,b,ones(n,1),ones(n,1),'heteroscedastic');
        verifyTrue(tc,ok);robust(j)=se^2;
        [~,~,~,se]=pipecompare.eval.Measure.eyeBias(x,b,ones(n,1),ones(n,1),'pooled');pooled(j)=se^2;
    end
    ratio=mean(robust)/var(effect);
    fprintf('UNEQUAL_VARIANCE prevalence=%.1f robust_variance_ratio=%.4f pooled_variance_ratio=%.4f\n',prevalence,ratio,mean(pooled)/var(effect));
    verifyGreaterThan(tc,ratio,0.85);verifyLessThan(tc,ratio,1.15);
    verifyGreaterThan(tc,abs(mean(pooled)/var(effect)-1),0.4);
end
end

function testUnequalVarianceHonorsWeightsAndPhysicalCounts(tc)
x=[-2;0;3;7;15;26];b=logical([0;0;0;1;1;1]);
w=[1;2;4;2;1;3];counts=[2;1;3;1;2;2];at=repelem((1:numel(x))',counts);
[a,sq,ok,se]=pipecompare.eval.Measure.eyeBias(x,b,w.*counts,counts,'heteroscedastic');
[c,csq,cok,cse]=pipecompare.eval.Measure.eyeBias(x(at),b(at),w(at),ones(numel(at),1),'heteroscedastic');
[d,dsq,dok,dse]=pipecompare.eval.Measure.eyeBias(x,b,100*w.*counts,counts,'heteroscedastic');
verifyEqual(tc,[a sq se],[c csq cse],'AbsTol',1e-10);verifyEqual(tc,ok,cok);
verifyEqual(tc,[a sq se],[d dsq dse],'AbsTol',1e-10);verifyEqual(tc,ok,dok);
end

function testUnequalVarianceDoesNotInventMissingGroupVariance(tc)
for b={logical([0;0;0;1]),logical([0;1;1;1]),true(4,1)}
    [a,q,ok,se]=pipecompare.eval.Measure.eyeBias((1:4)',b{1},ones(4,1),ones(4,1),'heteroscedastic');
    verifyFalse(tc,ok);verifyTrue(tc,all(isnan([a q se])));
end
[a,q,ok,se]=pipecompare.eval.Measure.eyeBias((1:4)',false(4,1),ones(4,1),ones(4,1),'heteroscedastic');
verifyTrue(tc,ok);verifyEqual(tc,[a q se],[0 0 0]);
end

function testPooledAndLegacyControlsKeepTheirDefinitions(tc)
x=[-2;0;2;11;15;19;23];b=logical([0;0;0;1;1;1;1]);p=mean(b);n=numel(x);
delta=mean(x(b))-mean(x(~b));
pooled=((nnz(b)-1)*var(x(b))+(nnz(~b)-1)*var(x(~b)))/(n-2);
expected=p^2*pooled*(1/nnz(b)+1/nnz(~b))+delta^2*p*(1-p)/(n-1);
[~,~,~,se]=pipecompare.eval.Measure.eyeBias(x,b,ones(n,1),ones(n,1),'pooled');
verifyEqual(tc,se^2,expected,'AbsTol',1e-12);
[~,~,~,se]=pipecompare.eval.Measure.eyeBias(x,b,ones(n,1),ones(n,1),'legacy-total');
verifyEqual(tc,se^2,p^2*var(x)*(1/nnz(b)+1/nnz(~b)),'AbsTol',1e-12);
end

function testUpgradeGuardRequiresEvidenceBeyondThePairedBound(tc)
% Round 1 flipped the guard default: the proposal is maintained unless the baseline is shown significantly better (reverse bound < 0); 'promote' is now 'maintain', and missing evidence is 'unassessed'.
[p,g]=pipecompare.eval.Rank.guardUpgrade([1;0.99],2,1,0.02);
verifyEqual(tc,p,2);verifyEqual(tc,g.decision,'maintain');
[p,g]=pipecompare.eval.Rank.guardUpgrade([1;0.5],2,1,0.1);
verifyEqual(tc,p,2);verifyEqual(tc,g.decision,'maintain');verifyLessThan(tc,g.upperDifference,0);
[p,g]=pipecompare.eval.Rank.guardUpgrade([1;0.5],2,1,Inf);
verifyEqual(tc,p,2);verifyEqual(tc,g.upperDifference,Inf);
[p,g]=pipecompare.eval.Rank.guardUpgrade([1;NaN],2,1,0.1);
verifyEqual(tc,p,2);verifyEqual(tc,g.upperDifference,Inf);
end

function testUpgradeGuardIsUsedByRankingAndStability(tc)
% Round 1 flipped the guard default (maintain unless the baseline is significantly better); the stability share can only move toward the baseline, so equality is the expected outcome when no resample shows the baseline significantly better.
n=80;x=sin((1:n)'*1.7);z=cos((1:n)'*1.3);z=z*std(x)/std(z);
c=[nqc_cand(x);nqc_cand(0.99*z)];ref=nqc_ref(n);
opts=pipecompare.eval.Rank.defaults();opts.nBoot=199;opts.nSplits=0;
R=pipecompare.eval.Rank.run(c,ref,opts);verifyEqual(tc,R.recommended,2);
opts.upgradeBaseline=1;G=pipecompare.eval.Rank.run(c,ref,opts);
verifyEqual(tc,G.recommended,2);verifyEqual(tc,G.byStratum.upgradeGuard.decision,'maintain');
verifyGreaterThanOrEqual(tc,G.table.chosenShare(1),R.table.chosenShare(1));
c(2)=nqc_cand(0.2*x);G=pipecompare.eval.Rank.run(c,ref,opts);
verifyEqual(tc,G.recommended,2);verifyEqual(tc,G.byStratum.upgradeGuard.decision,'maintain');
end

function testUpgradeGuardRefusesAnIneligibleBaseline(tc)
c=[nqc_cand((1:40)');nqc_cand(sin((1:40)'))];c(1).interpolatedFraction=0.9;
opts=pipecompare.eval.Rank.defaults();opts.upgradeBaseline=1;opts.nBoot=99;opts.nSplits=0;
verifyError(tc,@()pipecompare.eval.Rank.run(c,nqc_ref(40),opts),'PipeCompare:Baseline');
end

function testScoreArchiveRetainsReRankingInputsWithoutRawEeg(tc)
c=nqc_cand((1:40)');c.id=1;c.key='highpass(cutoff=0.3)';
c.root=rand(3,100);c.steps={struct('cutoffRefined',false)};
ref=nqc_ref(40);opts=pipecompare.eval.Rank.defaults();opts.nBoot=99;opts.nSplits=0;
R=pipecompare.eval.Rank.run(c,ref,opts);
leaf=struct('key',c.key,'path',{{struct('slot','hp','key',c.key)}});
r=struct('cands',c,'leaves',leaf,'labels',{{c.key}},'ranking',R,'ref',ref);
S=pipecompare_record_summary(r,'unit');
verifyEqual(tc,S.levels{1},{'hp',c.key});verifyFalse(tc,isfield(S.cands,'root'));
verifyTrue(tc,isfield(S.cands,'message'));verifyEqual(tc,S.cands.m.objectives.X,c.m.objectives.X);
file=[tempname '.mat'];cleanup=onCleanup(@()delete(file)); %#ok<NASGU>
save(file,'S','-v7');saved=load(file);
G=pipecompare.eval.Rank.run(saved.S.cands,saved.S.ref,opts,saved.S.levels);
verifyEqual(tc,G.recommended,R.recommended);verifyEqual(tc,G.table.objective,R.table.objective);
end


function testReferenceGuardRecognizesOnlyRepairBeforeAverageReference(tc)
keys={'reref(mode=average)','badchannels(action=interpolate) > reref(mode=average)', ...
    'reref(mode=average) > badchannels(action=interpolate)', ...
    'badchannels(action=mark) > reref(mode=average)','badchannels(action=interpolate) > reref(mode=channels)'};
c=struct('key',keys);
[repaired,needed]=pipecompare.eval.Rank.referenceProtection(c);
verifyEqual(tc,repaired,logical([0;1;0;0;0]));verifyEqual(tc,needed,logical([1;1;1;1;0]));
end

function testReferenceGuardProtectsCleanupOnlyWithoutDemonstratedImprovement(tc)
% Round 1 flipped the guard default: the cleanup guard now maintains the skip unless the repaired candidate is shown significantly better.
repaired=logical([1;0]);needed=true(2,1);
[k,g]=pipecompare.eval.Rank.guardReference([1;0.99],[1;0.99],[0;0],repaired,needed,2,0.1);
verifyEqual(tc,k,2);verifyEqual(tc,g.origin,'average-reference-cleanup');
[k,g]=pipecompare.eval.Rank.guardReference([1;0.5],[1;0.5],[0;0],repaired,needed,2,0.1);
verifyEqual(tc,k,2);verifyEqual(tc,g.decision,'maintain');
[k,g]=pipecompare.eval.Rank.guardReference([NaN;0.5],[NaN;0.5],[0;0],repaired,needed,2,0.1);
verifyEqual(tc,k,2);verifyEmpty(tc,g);
end

function testReferenceGuardIsAppliedToTheDefaultRankAndResamples(tc)
% Round 1 flipped the guard default (maintain unless significantly better); the guard still runs in the default rank (origin assertion), and the share comparison is one-directional, so equality is expected when no resample flips.
n=80;x=sin((1:n)'*1.7);z=cos((1:n)'*1.3);z=z*std(x)/std(z);
c=[nqc_cand(x,'key','badchannels(action=interpolate) > reref(mode=average)'); ...
    nqc_cand(0.99*z,'key','reref(mode=average)')];
opts=pipecompare.eval.Rank.defaults();opts.nBoot=199;opts.nSplits=0;
R=pipecompare.eval.Rank.run(c,nqc_ref(n),opts);
verifyEqual(tc,R.recommended,2);verifyEqual(tc,R.byStratum.upgradeGuard.origin,'average-reference-cleanup');
opts.protectReferenceCleanup=false;old=pipecompare.eval.Rank.run(c,nqc_ref(n),opts);
verifyEqual(tc,old.recommended,2);verifyGreaterThanOrEqual(tc,R.table.chosenShare(1),old.table.chosenShare(1));
opts.protectReferenceCleanup=true;c(1).interpolatedFraction=0.8;
R=pipecompare.eval.Rank.run(c,nqc_ref(n),opts);
verifyEqual(tc,R.recommended,2);verifyEqual(tc,R.table.status{1},'rejected');
end
