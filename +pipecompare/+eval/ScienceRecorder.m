classdef ScienceRecorder
    %SCIENCERECORDER Read-only observational ICA/ERP diagnostics.
    % Cannot alter optimization, eligibility, source EEG or ICA fit.
    % Called on the exact already-evaluated candidate, never on a refit.
    methods (Static)
        function z = capture(EEG, cand, contract, ref, seed)
            z = struct('schema',1,'seed',seed,'candidateId',cand.id, ...
                'available',false,'reason','','ica',struct([]), ...
                'epochs',struct([]),'waveform',struct([]), ...
                'bootstrap',struct([]));
            if ~strcmp(cand.status,'ok')
                z.reason='candidate-not-ok'; return;
            end
            if ~strcmp(contract.analysis,'erp') || numel(ref.ids)~=2 ...
                    || isempty(contract.components)
                z.reason='not-two-condition-erp'; return;
            end
            % Preserve each decomposition separately. Topographies are
            % columns of pinv(W*S); column signs and scale are arbitrary.
            ica=struct('weights',{},'sphere',{},'channelIndices',{}, ...
                       'mixing',{},'componentCount',{},'removedIndices',{});
            for j=1:numel(cand.ica)
                d=cand.ica{j};
                assert(isfield(d,'icaweights') && isfield(d,'icasphere') ...
                    && isfield(d,'icachansind'), 'PipeCompare:ScienceRecorder', ...
                    'Missing evaluated ICA decomposition.');
                W=double(d.icaweights);S=double(d.icasphere);
                mix=pinv(W*S);
                assert(size(mix,1)==numel(d.icachansind) && ...
                    size(mix,2)==size(W,1),'PipeCompare:ScienceRecorder', ...
                    'ICA channel/component topology mismatch.');
                removed=[];
                for k=1:numel(cand.steps)
                    step=cand.steps{k};
                    if strcmp(step.type,'icremove') && ~isempty(step.comps)
                        removed=union(removed,double(step.comps(:)'));
                    end
                end
                assert(all(removed>=1 & removed<=size(mix,2)), ...
                    'PipeCompare:ScienceRecorder','Removed IC outside decomposition.');
                ica(j)=struct('weights',W,'sphere',S, ...
                    'channelIndices',double(d.icachansind(:)'), ...
                    'mixing',mix,'componentCount',size(W,1), ...
                    'removedIndices',removed); %#ok<AGROW>
            end
            z.ica=ica;
            % Full Pz epoch waveform from the already processed E. Unlike
            % Measure.trials' scalar ROI mean, preserve every time sample.
            T=pipecompare.eval.Measure.trials(EEG,contract,true);
            labels={EEG.chanlocs.labels};
            pz=find(strcmpi(labels,'Pz'));
            assert(isscalar(pz) && EEG.trials>1, ...
                'PipeCompare:ScienceRecorder','Pz or epoched data unavailable.');
            t=double(EEG.xmin)+(0:EEG.pnts-1)/double(EEG.srate);
            b=t>=contract.baseline(1)-1e-9 & t<=contract.baseline(2)+1e-9;
            assert(any(b),'PipeCompare:ScienceRecorder','Baseline unavailable.');
            o=cand.m.objectives(1);
            gainCol=matlab.lang.makeValidName(['gain_' ref.objectives{1}]);
            % Gain comes from Rank, not from this recording stage. Reuse the
            % raw waveform and mark calibration absent; caller may supply
            % independently pinned Rank gain during analysis.
            waves=struct('condition',{},'sourceIds',{},'keptIds',{}, ...
                'weights',{},'waveform',{}, ...
                'meanAmplitude',{},'analyticSME',{});
            for c=1:2
                ids=double(ref.ids{c}(:));
                sel=find(T.cond==c);
                [known,loc]=ismember(double(T.id(sel)),ids);
                sel=sel(known);loc=loc(known);
                [loc,first]=unique(loc,'stable');sel=sel(first);
                w=zeros(numel(ids),1);
                if ~isempty(o.w) && numel(o.w)>=c && ~isempty(o.w{c})
                    raw=double(o.w{c}(:));raw(~isfinite(raw)|raw<=0)=0;
                    w=raw;
                else
                    w(loc)=1;
                end
                assert(numel(w)==numel(ids),'PipeCompare:ScienceRecorder', ...
                    'Trial weight/reference length mismatch.');
                x=zeros(numel(ids),numel(t));x(:)=NaN;
                E=double(EEG.data(pz,:,T.epoch(sel)));
                E=reshape(E,numel(t),[])';
                E=E-mean(E(:,b),2);
                x(loc,:)=E;
                usable=isfinite(x(:,1)) & w>0;
                assert(any(usable),'PipeCompare:ScienceRecorder', ...
                    'No weighted valid waveforms.');
                avg=(w(usable)'*x(usable,:))/sum(w(usable));
                waves(c)=struct('condition',char(string(ref.names{c})), ...
                    'sourceIds',ids','keptIds',ids(usable)', ...
                    'weights',w','waveform',avg, ...
                    'meanAmplitude',o.estimate(c), ...
                    'analyticSME',o.sme(c)); %#ok<AGROW>
            end
            z.waveform=struct('timesMs',1000*t,'site','Pz', ...
                'condition',{waves},'difference',waves(1).waveform-waves(2).waveform, ...
                'gainSource',gainCol,'gainApplied',false);
            z.epochs=struct('referenceCount',[numel(ref.ids{1}) numel(ref.ids{2})], ...
                'trialIdentitySource','EEGLAB urevent via Measure.trials', ...
                'sourceIDs',{{double(ref.ids{1}(:)') double(ref.ids{2}(:)')}});
            z.bootstrap=struct('available',false, ...
                'reason','bootstrap SME must be recomputed from per-trial scores by separately pinned collector');
            z.available=true;
        end
    end
end
