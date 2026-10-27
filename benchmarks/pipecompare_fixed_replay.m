function [EEG, audit]=pipecompare_fixed_replay(root,path,cand,contract,ref,opts)
% Replay a cached candidate's decisions on a structurally identical source.
% This is a conditional operator audit. It must never learn ICA, regressions,
% bad channels, rejection limits, component thresholds or trial weights.
if nargin<6,opts=struct();end
assert(numel(path)==numel(cand.steps),'PipeCompare:FixedReplay','Incomplete cached step facts.');
EEG=root;ica=0;ctx=struct('contract',contract,'ref',ref,'rank',opts,'highpass',0);
allowed={'badchannels','reref','ica','highpass','lowpass','linenoise','icremove','epoch','baseline', ...
    'reject_threshold','reject_jointprob','reject_kurtosis','eogregress','weight_trials'};
for k=1:numel(path)
    in=path{k};fact=cand.steps{k};
    assert(strcmp(in.type,fact.type)&&isequal(in.params,fact.params),'PipeCompare:FixedReplay','Path and cached decisions differ.');
    assert(any(strcmp(in.type,allowed)),'PipeCompare:FixedReplay','Unsupported step must not be silently refitted.');
    info=struct();
    switch in.type
        case 'ica'
            ica=ica+1;
            assert(isfield(cand,'ica')&&numel(cand.ica)>=ica,'PipeCompare:FixedReplay','Cached ICA is missing.');
            info.ica=cand.ica{ica};
        case 'badchannels'
            action=in.params.action;
            if strcmp(action,'remove'),labels=fact.removed;else,labels=fact.interpolated;end
            assert(any(strcmp(action,{'remove','interpolate'})),'PipeCompare:FixedReplay','Unsupported bad-channel action.');
            [ok,idx]=ismember(lower(labels),lower({EEG.chanlocs.labels}));
            assert(all(ok),'PipeCompare:FixedReplay','Cached repair channel is missing.');
            info.badIdx=idx;info.stretches=fact.stretches;
        case 'icremove'
            assert(isfield(fact,'comps'),'PipeCompare:FixedReplay','Cached component choice is missing.');
            info.comps=fact.comps;
        case {'reject_threshold','reject_jointprob','reject_kurtosis'}
            assert(~pipecompare.utils.fieldOr(in.params,'interpolate',0),'PipeCompare:FixedReplay', ...
                'This audit requires stored per-epoch interpolation decisions, unavailable in historical facts.');
            assert(isfield(fact,'rejIdx'),'PipeCompare:FixedReplay','Cached rejection decision is missing.');
            info.rejIdx=fact.rejIdx;
        case 'eogregress'
            assert(isstruct(fact.eog)&&isfield(fact.eog,'B'),'PipeCompare:FixedReplay','Cached EOG fit is missing.');
            info.eog=fact.eog;
        case 'weight_trials'
            assert(~isempty(fact.trialWeights),'PipeCompare:FixedReplay','Cached trial weights are missing.');
            info.trialWeights=fact.trialWeights;
    end
    [EEG,matched]=pipecompare.run.Steps.replayDecision(in,EEG,info,ctx);
    assert(matched,'PipeCompare:FixedReplay','A step was not decision-matched.');
    if strcmp(in.type,'highpass'),ctx.highpass=max(ctx.highpass,in.params.cutoff);end
end
EEG=pipecompare.run.Executor.selectEligible(EEG,contract);
audit=struct('scope','fixed-fitted-operator','fullyHeldOut',false,'preprocessingRefitted',false, ...
    'icaFits',0,'steps',numel(path));
end
