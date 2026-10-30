function tests = test_checkpoint_delivery
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(here); addpath(fullfile(here,'..')); addpath(fullfile(here,'..','benchmarks'));
end

function testSliceResumeMatchesNativeSearchAndReusesSavedWork(tc)
EEG = nqc_synth(struct('seconds',60,'nPerCond',10,'artifactTrials',0));
nqc_setBase(EEG);
p = pipecompare.plan.Plan(); p = p.add('reref','mode','average');
p = p.add('ica','fitHighpass',0,'fitClean',0);
p = p.add('highpass','cutoff',{0.1,0.5}); p = p.add('epoch'); p = p.add('baseline');
opts = struct('icaSeeds',0,'filterCheck',false,'nBoot',99,'nSplits',0, ...
    'artifactVariance','pooled','maxSeconds',600);
folder = tempname; mkdir(folder); cleanup = onCleanup(@()rmdir(folder,'s')); %#ok<NASGU>
fullOpts = opts; fullOpts.maxSeconds = Inf;
full = pipecompare.PipeCompare.optimize(p,nqc_c(),fullOpts);
opts.checkpoint = folder; opts.stopAfter = 1;
passport = struct('commit','checkpoint-proof','seed',40101,'arm','balanced', ...
    'artifactVariance','pooled','eeglabSha','fixed-test-eeglab');
[partial, state] = pipecompare_checkpoint_search(folder,passport, ...
    @()pipecompare.PipeCompare.optimize(p,nqc_c(),opts),600,'start');
verifyEmpty(tc,partial); verifyEqual(tc,state.status,'partial');
leaves = dir(fullfile(folder,'leaf_*.mat')); verifyEqual(tc,numel(leaves),1);
savedLeaf = readBytes(fullfile(leaves(1).folder,leaves(1).name));
cache = dir(fullfile(folder,'ica_cache','*.mat')); verifyEqual(tc,numel(cache),1);
savedIca = readBytes(fullfile(cache(1).folder,cache(1).name));
[txt, result, resumed] = evalc('pipecompare_checkpoint_search(folder,passport,@()mustNotStart(),600,''resume'')');
verifyTrue(tc,resumed.complete); verifyTrue(tc,contains(txt,'Reused exact ICA checkpoint'));
verifyEqual(tc,result.ranking.recommended,full.ranking.recommended);
verifyEqual(tc,result.ranking.table.objective,full.ranking.table.objective,'AbsTol',1e-10);
verifyEqual(tc,readBytes(fullfile(leaves(1).folder,leaves(1).name)),savedLeaf);
verifyEqual(tc,readBytes(fullfile(cache(1).folder,cache(1).name)),savedIca);
verifyEqual(tc,numel(dir(fullfile(folder,'ica_cache','*.mat'))),1);
[again, final] = pipecompare_checkpoint_search(folder,passport,@()mustNotStart(),600,'resume');
verifyTrue(tc,final.complete); verifyEqual(tc,again.ranking.recommended,result.ranking.recommended);
bad = passport; bad.commit = 'other-code';
verifyError(tc,@()pipecompare_checkpoint_search(folder,bad,@()mustNotStart(),600,'resume'),'PipeCompare:Slice');
bad = passport; bad.seed = 40102;
verifyError(tc,@()pipecompare_checkpoint_search(folder,bad,@()mustNotStart(),600,'resume'),'PipeCompare:Slice');
bad = passport; bad.artifactVariance = 'legacy-total';
verifyError(tc,@()pipecompare_checkpoint_search(folder,bad,@()mustNotStart(),600,'resume'),'PipeCompare:Slice');
end

function testResumeRejectsMissingPassportBeforeComputing(tc)
folder = tempname; mkdir(folder); cleanup = onCleanup(@()rmdir(folder,'s')); %#ok<NASGU>
passport = struct('commit','test','seed',40101,'arm','balanced','artifactVariance','pooled','eeglabSha','test');
verifyError(tc,@()pipecompare_checkpoint_search(folder,passport,@()mustNotStart(),10,'resume'),'PipeCompare:Slice');
end

function testShadowSlicesRejectInvalidInputsBeforeSearch(tc)
out = tempname; mkdir(out); cleanup = onCleanup(@()rmdir(out,'s')); %#ok<NASGU>
verifyError(tc,@()pipecompare_shadow_selection(out,20031,2,'mode','slice','arm','balanced', ...
    'checkpointDir',out,'maxSeconds',60,'headSha','testcode'),'PipeCompare:Shadow');
verifyError(tc,@()pipecompare_shadow_selection(out,20031,1,'mode','resume','arm','balanced', ...
    'checkpointDir',out,'maxSeconds',60,'headSha','testcode'),'PipeCompare:Shadow');
end

function bytes = readBytes(file)
fid = fopen(file,'rb'); clean = onCleanup(@()fclose(fid)); %#ok<NASGU>
bytes = fread(fid,Inf,'*uint8');
end

function r = mustNotStart() %#ok<STOUT>
error('test:UnexpectedCompute','Resume must not start a new search.');
end
