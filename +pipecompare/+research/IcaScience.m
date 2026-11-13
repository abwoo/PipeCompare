classdef IcaScience
    %ICASCIENCE Read-only diagnostics for developmental ICA-start research.
    % Never chooses a candidate, fits ICA, changes EEG, or changes Rank.
    % Captures the actually evaluated candidate, after its decisions.
    methods (Static)
        function z = capture(EEG, candidate, ref, contract)
            z = struct('status','unavailable','reason','', ...
                'decompositions',{{}},'timeMs',[], 'conditions',{{}}, ...
                'measurementKind','','srate',NaN,'channel','', ...
                'preRemovalContributions',{{}},'postRemovalSignals',{{}});
            assert(strcmp(contract.analysis,'erp') && numel(ref.ids)==2, ...
                'PipeCompare:IcaScience','This readout is only for two-condition ERP.');
            assert(isfield(candidate,'ica') && ~isempty(candidate.ica), ...
                'PipeCompare:IcaScience','No fitted ICA in this candidate.');
            assert(strcmp(candidate.status,'ok') && isfield(candidate,'m') && ...
                ~isempty(candidate.m.objectives), 'PipeCompare:IcaScience', ...
                'Candidate has no actual evaluated measurement.');
            assert(EEG.trials>1 && EEG.pnts>1 && EEG.nbchan==numel(EEG.chanlocs), ...
                'PipeCompare:IcaScience','Expected epoched evaluated EEG.');
            z.decompositions = cell(1,numel(candidate.ica));
            for q=1:numel(candidate.ica)
                d=candidate.ica{q};
                W=double(d.icaweights)*double(d.icasphere);
                ix=double(d.icachansind(:)');
                assert(ismatrix(W) && all(isfinite(W(:))) && ...
                    size(W,2)==numel(ix) && numel(unique(ix))==numel(ix) && ...
                    all(ix>=1) && all(ix<=numel(EEG.chanlocs)), ...
                    'PipeCompare:IcaScience','ICA decomposition channels misaligned.');
                A=pinv(W);
                assert(all(isfinite(A(:))), 'PipeCompare:IcaScience','ICA maps nonfinite.');
                labels={EEG.chanlocs(ix).labels};
                assert(numel(unique(lower(string(labels))))==numel(labels), ...
                    'PipeCompare:IcaScience','Duplicate ICA montage channel labels.');
                removed=[];
                classes={}; probs=[];
                for k=1:numel(candidate.steps)
                    f=candidate.steps{k};
                    if ~strcmp(f.type,'icremove'),continue;end
                    removed=union(removed,double(f.comps(:)'));
                    if isfield(f,'iclabelClasses') && ~isempty(f.iclabelClasses)
                        classes=cellstr(f.iclabelClasses);
                        probs=double(f.iclabelProbabilities);
                    end
                end
                assert(all(removed>=1) && all(removed<=size(W,1)) && ...
                    all(removed==fix(removed)), 'PipeCompare:IcaScience', ...
                    'Removed component ID is out of decomposition range.');
                if ~isempty(probs)
                    assert(size(probs,1)==size(W,1) && ...
                        size(probs,2)==numel(classes) && all(isfinite(probs(:))), ...
                        'PipeCompare:IcaScience','ICLabel class/component alignment failed.');
                end
                z.decompositions{q}=struct('nComponents',size(W,1), ...
                    'nChannels',numel(ix),'channelLabels',{labels}, ...
                    'mixing',A,'removed',removed,'classNames',{classes}, ...
                    'classProbabilities',probs, ...
                    'componentUnits','columns of ICA mixing map; normalized before matching');
            end
            o=candidate.m.objectives(1);
            assert(strcmp(o.kind,'scalar') && numel(o.X)==2 && ...
                numel(contract.components)==1 && ...
                strcmpi(contract.components(1).name,'P3'), ...
                'PipeCompare:IcaScience','Requires scalar ERP CORE P3 objective.');
            z.measurementKind='P3 mean amplitude, uncorrected waveform';
            labels={EEG.chanlocs.labels};
            at=find(strcmpi(labels,'Pz'));
            assert(numel(at)==1,'PipeCompare:IcaScience','Pz channel unavailable.');
            z.channel='Pz';
            z.srate=double(EEG.srate);
            z.timeMs=1000*(double(EEG.xmin)+(0:EEG.pnts-1)/double(EEG.srate));
            assert(all(isfinite(z.timeMs)), 'PipeCompare:IcaScience', ...
                'Missing epoch time coordinates.');
            T=pipecompare.eval.Measure.trials(EEG, contract, true);
            conditions=cell(1,2);
            for c=1:2
                ids=double(ref.ids{c}(:));
                rows=find(T.cond==c);
                matched=double(T.id(rows));
                [found,loc]=ismember(ids,matched);
                assert(numel(unique(matched))==numel(matched), ...
                    'PipeCompare:IcaScience','Duplicate time-locking urevents.');
                scores=double(o.X{c}(:,1));
                assert(numel(scores)==numel(ids), ...
                    'PipeCompare:IcaScience','Scores do not align to reference identities.');
                weights=ones(size(scores));
                if isfield(o,'w') && ~isempty(o.w)
                    weights=double(o.w{c}(:));
                end
                assert(numel(weights)==numel(ids), ...
                    'PipeCompare:IcaScience','Weights do not align to reference identities.');
                valid=isfinite(scores)&isfinite(weights)&weights>0;
                assert(all(found(valid)),'PipeCompare:IcaScience', ...
                    'A positive-weight scored trial has no surviving epoch.');
                assert(all(weights(valid)>=0),'PipeCompare:IcaScience', ...
                    'Invalid weighted epoch.');
                flags=false(size(ids));
                if ~isempty(ref.blink) && ~isempty(ref.blink{1}) && ...
                        numel(ref.blink{1})>=c
                    flags=logical(ref.blink{1}{c}(:));
                end
                assert(numel(flags)==numel(ids),'PipeCompare:IcaScience', ...
                    'Blink flags not aligned to urevent IDs.');
                use=find(valid);
                X=zeros(numel(use),EEG.pnts);
                for j=1:numel(use)
                    epoch=T.epoch(rows(loc(use(j))));
                    X(j,:)=reshape(double(EEG.data(at,:,epoch)),1,[]);
                end
                w=weights(use);
                wave=(w'*X)/sum(w);
                p3=z.timeMs>=300-1e-8 & z.timeMs<=600+1e-8;
                assert(any(p3),'PipeCompare:IcaScience','Missing P300 window.');
                windowMean=mean(wave(p3));
                assert(isfinite(windowMean) && ...
                    abs(windowMean-o.estimate(c))<=1e-5*max(1,abs(o.estimate(c))), ...
                    'PipeCompare:IcaScience', ...
                    'Independent full-waveform P300 integral differs from original score.');
                conditions{c}=struct('name',ref.names{c},'referenceIds',ids', ...
                    'retainedIds',ids(use)', 'referenceFlagged',flags', ...
                    'trialScores',scores','trialWeights',weights', ...
                    'usedCount',numel(use),'effectiveN',sum(w)^2/sum(w.^2), ...
                    'meanWaveform',wave,'windowMean',windowMean, ...
                    'originalEstimate',double(o.estimate(c)));
            end
            z.conditions=conditions;
            % Copy only observer evidence already recorded by the search;
            % never replay, infer, or alter any IC cleaning decision here.
            for step=1:numel(candidate.steps)
                info=candidate.steps{step};
                if strcmp(info.type,'icremove') && ...
                        isfield(info,'icaContribution')
                    z.preRemovalContributions{end+1}=info.icaContribution;
                end
                if strcmp(info.type,'icremove') && ...
                        isfield(info,'icaPostRemoval')
                    z.postRemovalSignals{end+1}=info.icaPostRemoval;
                end
            end
            z.status='complete';z.reason='';
        end
    end
end
