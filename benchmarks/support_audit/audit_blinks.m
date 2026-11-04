        function [B, info] = audit_blinks(EEG, contract, ref)
            % Which reference trials have an artifact that can move the mean
            % of a mean-amplitude component's score: B{j}{c} (ref.n(c) x 1
            % logical) for component j and condition c, {} when it cannot
            % be told (band power, or no channel to tell it from). Three
            % kinds, in the data as they are (before any processing), each
            % found when its robust z (median and MAD over all reference
            % trials) exceeds BlinkZ; a trial with any of them counts:
            %   blinks: the largest, over the eye channels
            %   (Injection.eyeChannels: EOG by type or name, else Fp1, Fp2,
            %   FPz; never a measured electrode), of the absolute window
            %   mean minus the baseline mean (what an eye artifact adds to
            %   the score is carried by that difference);
            %   eye movements: the same of the difference of a left and a
            %   right channel (horizontalPair), in which a blink, the same
            %   on both sides, cancels and a horizontal movement doubles;
            %   muscle: the log of the largest, over the measured
            %   electrodes, SD of the first difference in the window (fast
            %   activity, which a slow ERP hardly has).
            % In continuous data, a trial also counts when a blink or eye
            % movement found sample by sample (sampleHits) falls in its
            % window or baseline: comparing the trials with each other
            % finds few when most trials have a blink (n-back tasks with a
            % response on every trial).
            % info: channels (eye channels' labels), horizontal ('L-R' or
            % ''), share (of the trials, per component), kinds and
            % kindShare (components x kinds).
            kinds = {'blinks', 'eye movements', 'muscle'};
            B = {}; info = struct('channels', {{}}, 'horizontal', '', 'share', [], 'kinds', {kinds}, 'kindShare', []);
            if ~strcmp(contract.analysis, 'erp') || isempty(contract.components), return; end
            labs = lower({EEG.chanlocs.labels});
            roi = find(ismember(labs, lower(contract.allRoi())));
            eyes = setdiff(pipecompare.eval.Injection.eyeChannels(EEG.chanlocs), roi);
            [hl, hr] = horizontalPair(EEG.chanlocs, eyes, roi);
            use = unique([eyes hl hr roi]);
            if isempty(use), return; end
            info.channels = {EEG.chanlocs(eyes).labels};
            if ~isempty(hl), info.horizontal = sprintf('%s-%s', EEG.chanlocs(hl).labels, EEG.chanlocs(hr).labels); end
            E = EEG; E.icaact = []; E.icaweights = []; E.icasphere = []; E.icawinv = []; E.icachansind = [];
            [~, E] = evalc('pop_select(E, ''channel'', use)');
            H = [];
            if E.trials == 1
                H = sampleHits(E, ismember(use, eyes), ismember(use, hl), ismember(use, hr), contract);
                [~, E] = evalc('pop_epoch(E, contract.allEvents(), contract.epoch, ''epochinfo'', ''yes'')');
            end
            T = pipecompare.eval.Measure.trials(E, contract, true);
            times = E.xmin + (0:E.pnts-1) / E.srate;
            D = double(E.data);
            De = D(ismember(use, eyes), :, :); Dh = [];
            if ~isempty(hl), Dh = D(use == hl, :, :) - D(use == hr, :, :); end
            bsel = false(size(times));
            if ~isempty(contract.baseline)
                bsel = times >= contract.baseline(1) - 1e-9 & times <= contract.baseline(2) + 1e-9;
                if any(bsel)
                    De = De - mean(De(:, bsel, :), 2);
                    if ~isempty(Dh), Dh = Dh - mean(Dh(:, bsel, :), 2); end
                end
            end
            z = @(v) (v - median(v)) / max(1.4826 * median(abs(v - median(v))), eps);
            perTrial = @(v) reshape(v, [], 1);
            nJ = numel(contract.components); nC = numel(ref.ids);
            B = cell(1, nJ); info.share = nan(1, nJ); info.kindShare = nan(nJ, numel(kinds));
            for j = 1:nJ
                B{j} = arrayfun(@(c) false(ref.n(c), 1), 1:nC, 'UniformOutput', false);
                comp = contract.components(j);
                w = times >= comp.window(1) - 1e-9 & times <= comp.window(2) + 1e-9;
                if ~strcmp(comp.measure, 'mean') || ~any(w), continue; end
                K = false(numel(T.epoch), numel(kinds));
                if ~isempty(De)
                    v = perTrial(max(abs(mean(De(:, w, :), 2)), [], 1)); K(:, 1) = z(v(T.epoch)) > pipecompare.eval.Measure.BlinkZ;
                end
                if ~isempty(Dh)
                    v = perTrial(abs(mean(Dh(:, w, :), 2))); K(:, 2) = z(v(T.epoch)) > pipecompare.eval.Measure.BlinkZ;
                end
                epochKinds=K;
                robustHit=false(size(K,1),2);hmmHit=robustHit;
                if ~isempty(H)
                    % also the trials with a blink or eye movement found in
                    % the continuous data, in the window or the baseline
                    rh=reshape(any(H.auditRobust(:,w | bsel,:),2),size(H.auditRobust,1),[]);
                    hh=reshape(any(H.auditHmm(:,w | bsel,:),2),size(H.auditHmm,1),[]);
                    robustHit(:,1:size(rh,1))=rh(:,T.epoch)';hmmHit(:,1:size(hh,1))=hh(:,T.epoch)';
                    hit = reshape(any(H.data(:, w | bsel, :), 2), size(H.data, 1), []);
                    K(:, 1) = K(:, 1) | hit(1, T.epoch)';
                    if size(hit, 1) > 1, K(:, 2) = K(:, 2) | hit(2, T.epoch)'; end
                end
                r = ismember(use, roi(ismember(labs(roi), lower(comp.roi))));
                if any(r) && sum(w) > 2
                    v = perTrial(log(max(std(diff(D(r, w, :), 1, 2), 0, 2), [], 1) + eps));
                    K(:, 3) = z(v(T.epoch)) > pipecompare.eval.Measure.BlinkZ;
                end
                epochKinds(:,3)=K(:,3);
                isB = any(K, 2);
                info.audit{j}=struct('ids',T.id,'conditions',T.cond,'epochKinds',epochKinds,'finalKinds',K,'sampleRobust',robustHit,'sampleHmm',hmmHit,'combined',isB);
                for c = 1:nC
                    rr = find(T.cond == c);
                    [tf, loc] = ismember(T.id(rr), ref.ids{c});
                    B{j}{c}(loc(tf)) = isB(rr(tf));
                end
                info.share(j) = mean(isB); info.kindShare(j, :) = mean(K, 1);
            end
        end


function H = sampleHits(E, isEye, isL, isR, contract)
% The blinks and eye movements of the continuous data E, found sample by
% sample (the trial rule compares the trials with each other, so when
% most trials have a blink the median trial has one, and few are found):
% where an eye channel, or the left minus the right channel, minus its
% moving median over 1 s (slow drifts out) is beyond 5 robust SDs (median
% and MAD over the samples: blinks fill far fewer samples than trials),
% or in the active state of a two-state HMM fitted to it
% (pipecompare.eval.Hmm.events: weak blinks too), with 50 ms on either
% side. Epoched like the data: one row per kind (blinks; eye movements
% when there is a pair). [] without eye channels.
H = [];
if ~any(isEye), return; end
rows = {double(E.data(isEye, :))};
if any(isL) && any(isR), rows{2} = double(E.data(isL, :)) - double(E.data(isR, :)); end
w = round(0.05 * E.srate); M = zeros(3*numel(rows), E.pnts);
for k = 1:numel(rows)
    Y = rows{k} - movmedian(rows{k}, round(E.srate), 2);
    Y = abs(Y - median(Y, 2));
    robust=any(Y > 5 * max(1.4826 * median(Y,2),eps),1);
    hmm=pipecompare.eval.Hmm.events(rows{k},E.srate);on=robust|hmm;
    M(k,:)=conv(double(on),ones(1,2*w+1),'same')>0;
    M(numel(rows)+k,:)=conv(double(robust),ones(1,2*w+1),'same')>0;
    M(2*numel(rows)+k,:)=conv(double(hmm),ones(1,2*w+1),'same')>0;
end
H = E; H.data = single(M); H.nbchan = size(M, 1);
H.chanlocs = struct('labels', arrayfun(@(k) sprintf('hit%d', k), 1:H.nbchan, 'UniformOutput', false));
[~, H] = evalc('pop_epoch(H, contract.allEvents(), contract.epoch, ''epochinfo'', ''yes'')');
H.data = H.data > 0.5;
H.auditRobust=H.data(numel(rows)+(1:numel(rows)),:,:);H.auditHmm=H.data(2*numel(rows)+(1:numel(rows)),:,:);H.data=H.data(1:numel(rows),:,:);H.nbchan=numel(rows);
end

function [l, r] = horizontalPair(chanlocs, eyes, roi)
% The left and the right channel whose difference shows horizontal eye
% movements: of the EOG channels, the leftmost and the rightmost by
% position (EEGLAB's Y axis points to the left ear) when they lie on
% either side, else those whose names end in L and R (e.g. EYEL and EYER);
% without EOG channels, or when those give no pair, F7 and F8. Empty when
% there is none or one of them is measured.
l = []; r = [];
eog = [];
if numel(eyes) >= 2, eog = eyes(pipecompare.simple.Presets.isNonEeg(chanlocs(eyes))); end
if numel(eog) >= 2
    y = nan(1, numel(eog));
    for k = 1:numel(eog)
        c = chanlocs(eog(k));
        if isfield(c, 'Y') && ~isempty(c.X) && ~isempty(c.Y) && ~isempty(c.Z)
            y(k) = c.Y / max(norm([c.X c.Y c.Z]), eps);
        end
    end
    [ym, a] = max(y); [yn, b] = min(y);
    if ym > 0.2 && yn < -0.2, l = eog(a); r = eog(b); return; end
    labs = upper(strtrim({chanlocs(eog).labels}));
    a = find(~cellfun(@isempty, regexp(labs, 'L$', 'once'))); b = find(~cellfun(@isempty, regexp(labs, 'R$', 'once')));
    if isscalar(a) && isscalar(b), l = eog(a); r = eog(b); return; end
end
labs = upper({chanlocs.labels});
a = find(strcmp(labs, 'F7')); b = find(strcmp(labs, 'F8'));
if isscalar(a) && isscalar(b) && ~any(ismember([a b], roi)), l = a; r = b; end
end

