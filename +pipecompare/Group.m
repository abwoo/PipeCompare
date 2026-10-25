classdef Group
    %GROUP One pipeline for a group of recordings.
    %
    %   G = pipecompare.Group.run(data, 'measure', 'P3', 'events', {...}, ...)
    %
    %   data: a cell of .set file names, or EEGLAB datasets (ALLEEG or some
    %   of it). The other arguments are those of pop_pipecompare (the
    %   simple mode), without 'show' and 'progress'. Every recording is
    %   compared on its own first, exactly as pop_pipecompare does with
    %   'show', 'off': each pipeline makes its decisions (bad channels,
    %   components, rejected epochs, limits) on that recording's own data.
    %   Then the pipelines run on every recording are ranked together:
    %
    %   - A pipeline passes when it passes on every recording. Recordings on
    %     which no pipeline passes are left out and named (G.leftOut).
    %   - Its group objective is the root mean square over the recordings
    %     of its total error (Rank step 3b: the gain-corrected SME with the
    %     bias of what blinks add and of the signal's unknown shape).
    %     This summarizes individual measurement errors in the observed
    %     group. It is not the MSE of a group mean: same-direction bias
    %     across recordings does not average away like independent noise.
    %   - Its uncertainty is that of each recording's paired bootstrap of
    %     its trials (Rank.bootstrap, as one recording's ranking), combined
    %     the same way; the pipelines not distinguished from the best and
    %     those equivalent to it follow as in the ranking
    %     (Rank.confidenceSet).
    %   - The recommendation, among those not distinguished from the best:
    %     the most trials kept in its worst recording, then the smallest
    %     largest signal change, then the lowest objective.
    %
    %   G.table: one row per pipeline run on every recording (id, pipeline,
    %   status, objective, notDistinguished, equivalent, minRetention: in
    %   the worst recording, ampError: the largest, passedIn: on how many
    %   recordings it passed); G.perRecording: the objective of each
    %   pipeline (columns) on each recording (rows); G.recommended (a row of
    %   G.table, [] when no pipeline passes on every recording);
    %   G.recordings, G.leftOut. Also stored in the base variable
    %   pipecompare_group.
    %
    %   Each recording's own pipeline, helped by the group (G.own; two-level
    %   empirical Bayes, Group.twoLevel): the log of each pipeline's
    %   objective on each recording is a recording effect plus a pipeline
    %   effect (fitted over all recordings) plus what this pipeline does on
    %   this recording beyond them. Its posterior gain retains the full
    %   within-recording paired bootstrap covariance (the scalar factor
    %   tau^2/(tau^2+v) when that covariance is diagonal). tau^2 is fitted
    %   after projecting out recording and pipeline effects. Each recording
    %   gets the pipeline lowest by this estimate among those passing on it.
    %   With tau^2 = 0 recordings with the same eligible options get the
    %   same pipeline; with
    %   precise recordings each gets its own best. G.own: a table (recording,
    %   id, pipeline, sameAsGroup); G.ownShrunk (recordings x pipelines);
    %   G.share: the mean diagonal of the posterior gain matrices.
    %
    %   pipecompare.Group.writeScript(G, file) writes the recommended
    %   pipeline as a function to run on each recording (its decisions are
    %   made on that recording's own data); writeScript(G, file, i) the
    %   pipeline of recording i in G.own.

    methods (Static)
        function G = run(data, varargin)
            if iscell(data) || isstring(data)
                files = cellstr(data); n = numel(files); names = files;
            else
                files = {}; n = numel(data);
                names = arrayfun(@(e) char(e.setname), data, 'UniformOutput', false);
            end
            assert(n >= 2, 'PipeCompare:Group', 'A group needs at least two recordings.');
            if ~isempty(varargin), [varargin{:}] = convertStringsToChars(varargin{:}); end
            drop = false(1, numel(varargin));
            for k = 1:2:numel(varargin)
                if any(strcmpi(varargin{k}, {'show', 'progress'})), drop(k:k+1) = true; end
            end
            args = [varargin(~drop) {'show', 'off'}];
            saved = saveBase(); restore = onCleanup(@() restoreBase(saved)); %#ok<NASGU>
            S = struct('name', {}, 'keys', {}, 'labels', {}, 'table', {}, 'cands', {}, 'ref', {}, 'ranking', {});
            first = []; firstIdx = 0;
            for i = 1:n
                pipecompare.utils.log('Recording %d of %d: %s', i, n, names{i});
                if isempty(files)
                    EEG = data(i);
                else
                    [~, EEG] = evalc('pop_loadset(files{i})');
                end
                assignin('base', 'EEG', EEG);
                [~, ~, res] = pop_pipecompare(EEG, args{:});
                keep = intersect(fieldnames(res.cands), {'id', 'key', 'stratum', 'status', 'm', 'signal', 'interpolatedFraction'});
                cands = rmfield(res.cands, setdiff(fieldnames(res.cands), keep));
                [~, order] = sort([cands.id]); cands = cands(order);   % cands(k): pipeline k
                S(i) = struct('name', names{i}, 'keys', {{res.leaves.key}}, 'labels', {res.labels}, ...
                    'table', res.ranking.table, 'cands', {cands}, 'ref', res.ref, 'ranking', ...
                    rmfield(res.ranking, setdiff(fieldnames(res.ranking), {'objective', 'options'})));
                if isempty(first) && any(strcmp(res.ranking.table.status, 'feasible')), first = res; firstIdx = i; end
            end
            G = pipecompare.Group.rank(S);
            G.first = first; G.firstRecording = firstIdx;
            assignin('base', 'pipecompare_group', G);
            pipecompare.Group.print(G);
        end

        function G = rank(S)
            % The group ranking of per-recording summaries S (see run).
            n = numel(S);
            keys = S(1).keys;
            for i = 2:n, keys = intersect(keys, S(i).keys, 'stable'); end
            used = arrayfun(@(s) any(strcmp(s.table.status, 'feasible')), S);
            nK = numel(keys); nU = sum(used);
            G = struct('recordings', {{S.name}}, 'leftOut', {{S(~used).name}}, 'keys', {keys}, 'labels', {{}}, ...
                'perRecording', nan(nU, nK), 'table', table(), 'recommended', [], 'note', '', ...
                'validationScope','observed-recordings','fullyHeldOut',false);
            assert(nK > 0, 'PipeCompare:Group', 'No pipeline was run on every recording.');
            [~, loc1] = ismember(keys, S(1).keys);
            G.labels = S(1).labels(loc1);
            feas = true(1, nK); ret = inf(1, nK); amp = -inf(1, nK); passed = zeros(1, nK);
            strata = {}; u = find(used);
            for j = 1:nU
                s = S(u(j)); [~, loc] = ismember(keys, s.keys); T = s.table;
                G.perRecording(j, :) = T.objective(loc)';
                f = strcmp(T.status(loc), 'feasible')';
                feas = feas & f; passed = passed + f;
                ret = min(ret, T.minRetention(loc)'); amp = max(amp, T.ampError(loc)');
                strata = union(strata, T.stratum(loc(f)));
            end
            assert(numel(strata) <= 1, 'PipeCompare:Group', ['The pipelines differ in what is measured (e.g. the ', ...
                'reference); compare one such choice at a time.']);
            J = sqrt(mean(G.perRecording .^ 2, 1));
            status = repmat({'excluded'}, nK, 1); status(feas) = {'passed'};
            nd = false(nK, 1); eq = false(nK, 1);
            F = find(feas);
            % each recording's bootstrap of the pipelines passing on it
            Fm = false(nU, nK); V = nan(nU, nK); boots = cell(1, nU); covariance = cell(nU, 1);
            baseOpts = pipecompare.utils.withDefaults(S(1).ranking.options, pipecompare.eval.Rank.defaults());
            seeds = pipecompare.Group.recordingSeeds({S.name}, baseOpts.seed);
            G.bootstrapSeeds = seeds(u);
            for j = 1:nU
                s = S(u(j)); [~, loc] = ismember(keys, s.keys);
                Fm(j, :) = strcmp(s.table.status(loc), 'feasible')';
                if ~any(Fm(j, :)), continue; end
                bootOpts = s.ranking.options; bootOpts.seed = seeds(u(j));
                b = pipecompare.eval.Rank.bootstrap(s.cands(loc(Fm(j, :))), s.ref, bootOpts, {s.ranking.objective});
                boots{j} = nan(nK, size(b{1}, 2)); boots{j}(Fm(j, :), :) = b{1};
                V(j, Fm(j, :)) = var(b{1}, 0, 2)';
                covariance{j} = nan(nK);
                covariance{j}(Fm(j,:), Fm(j,:)) = cov(b{1}');
            end
            G.own = table(); G.ownShrunk = nan(nU, nK); G.share = NaN;
            if nU >= 2
                Y = G.perRecording; Y(~Fm) = NaN;
                [G.ownShrunk, G.share, ~, G.twoLevelModel] = pipecompare.Group.twoLevel(Y, V, covariance);
                th = G.ownShrunk; th(~Fm) = Inf;
                [best, own] = min(th, [], 2); own(~isfinite(best)) = NaN;
                G.ownChoice = own;
            end
            if nU == 0
                G.note = 'No pipeline passed on any recording.';
            elseif isempty(F)
                G.note = 'No pipeline passed on every recording.';
            else
                acc = 0;
                for j = 1:nU, acc = acc + boots{j}(F, :) .^ 2; end
                Bq = sqrt(acc / nU);
                [~, ib] = min(J(F));
                opts = pipecompare.utils.withDefaults(S(u(1)).ranking.options, pipecompare.eval.Rank.defaults());
                [keep, equiv] = pipecompare.eval.Rank.confidenceSet(J(F)', Bq, ib, opts.alpha, opts.equivalenceMargin);
                nd(F(keep)) = true; eq(F(equiv)) = true;
                t = F(keep);
                [~, o] = sortrows(round([-ret(t)', amp(t)', J(t)'], 10));
                G.recommended = t(o(1));
            end
            G.table = table((1:nK)', G.labels(:), status, J', nd, eq, ret', amp', passed', ...
                'VariableNames', {'id', 'pipeline', 'status', 'objective', 'notDistinguished', 'equivalent', ...
                'minRetention', 'ampError', 'passedIn'});
            if isfield(G, 'ownChoice')
                own = G.ownChoice; ok = isfinite(own);
                lab = repmat({''}, nU, 1); lab(ok) = G.labels(own(ok));
                same = false(nU, 1); if ~isempty(G.recommended), same = own == G.recommended; end
                G.own = table(G.recordings(u)', own, lab, same, ...
                    'VariableNames', {'recording', 'id', 'pipeline', 'sameAsGroup'});
                G = rmfield(G, 'ownChoice');
            end
        end

        function seeds = recordingSeeds(names, seed)
            % Independent streams across recordings, paired across pipelines
            % within one recording. Unique names retain their stream when
            % input order changes; repeated names get separate occurrences.
            keys = cell(size(names));
            for i = 1:numel(names)
                occurrence = sum(strcmp(names(1:i), names{i}));
                keys{i} = sprintf('%s|%d', names{i}, occurrence);
            end
            [~, order] = sort(keys); seeds = zeros(size(names)); used = [];
            % Leave room for the deterministic salts used by inner checks.
            modulus = 2147483647;
            for i = order(:)'
                h = mod(double(seed), modulus);
                for c = double(keys{i}), h = mod(h * 65599 + c, modulus); end
                while any(used == h), h = mod(h + 1, modulus); end
                seeds(i) = h; used(end+1) = h; %#ok<AGROW>
            end
        end

        function [theta, share, tau2, model] = twoLevel(Y, V, covariance)
            % Two-level empirical Bayes estimate of each recording's
            % objective for each pipeline (see the class help). Y, V:
            % recordings x pipelines objectives and their bootstrap
            % variances, NaN where the pipeline does not pass. theta: the
            % estimates (NaN where Y is); share: the mean shrinkage factor
            % over the cells; tau2: the variance between recordings of what
            % a pipeline does beyond the two effects (log scale).
            % covariance (optional): one pipelines x pipelines raw-objective
            % covariance per recording, from its PAIRED bootstrap. With no
            % covariance supplied, V defines a diagonal observation model.
            % Effects use GLS; residual uncertainty and posterior gains keep
            % the covariance. This is a plug-in model for observed records,
            % not independent validation on unseen participants.
            assert(isequal(size(Y),size(V)), 'PipeCompare:Group', 'Objectives and variances must have the same size.');
            paired = nargin >= 3;
            if paired
                assert(iscell(covariance) && numel(covariance)==size(Y,1), ...
                    'PipeCompare:Group', 'Supply one covariance matrix per recording.');
            end
            ok = isfinite(Y) & Y > 0 & isfinite(V) & V >= 0;
            [nR, nP] = size(Y);
            [record, pipeline] = find(ok); m = numel(record);
            theta = nan(size(Y)); share = NaN; tau2 = 0;
            model = struct('scope','observed-recordings','fullyHeldOut',false, ...
                'pairedCovariance',paired,'cells',m,'components',0,'df',0,'varianceAssessed',false, ...
                'pipelineEffects',nan(1,nP),'pipelineComponents',nan(1,nP),'tau2',NaN, ...
                'trainingRecords',nnz(any(ok,2)),'recordsPerPipeline',sum(ok,1),'minimumTransferRecords',3);
            if m == 0, return; end
            [D, components, columns, component] = groupDesign(record,pipeline,nR,nP);
            model.pipelineEffects = nan(1,nP);
            model.pipelineComponents = component(nR+(1:nP));
            model.tau2 = NaN;
            model.components = components; model.df = m-size(D,2);
            model.varianceAssessed = model.df > 0;
            y = log(Y(ok)); residual = y-D*(D\y);
            C = cell(nR,1); ix = cell(nR,1);
            for s = unique(record)'
                ix{s} = find(record==s); p = pipeline(ix{s});
                raw = diag(V(s,p));
                if paired
                    assert(isequal(size(covariance{s}),[nP nP]), 'PipeCompare:Group', 'Covariance matrix has the wrong size.');
                    raw = covariance{s}(p,p);
                    assert(all(isfinite(raw(:))), 'PipeCompare:Group', 'Covariance contains nonfinite observed entries.');
                    tolerance = 1e-10*max([V(s,p),realmin]);
                    assert(max(abs(diag(raw)-V(s,p)')) <= tolerance, 'PipeCompare:Group', 'Covariance diagonal must match V.');
                end
                C{s} = groupCovariance(raw./(Y(s,p)'*Y(s,p)));
            end
            if model.df == 0
                theta(ok) = Y(ok); share = 0; return;
            end
            % E||R y||^2 = tau^2 tr(R)+tr(R Sigma). Diagonal variances
            % alone over-subtract noise shared by a recording's pipelines.
            normal = decomposition(D'*D,'chol'); noise = 0;
            for s = unique(record)'
                Ds = D(ix{s},:); Hss = Ds*(normal\Ds');
                noise = noise+trace(C{s})-trace(Hss*C{s});
            end
            tau2 = max(0,(sum(residual.^2)-max(noise,0))/model.df);
            Z = sparse(size(D,1),size(D,2)); z = zeros(m,1); cholC = cell(nR,1);
            for s = unique(record)'
                rows = ix{s}; np = numel(rows); S = C{s}+tau2*eye(np);
                [L,failed] = chol(S,'lower');
                if failed
                    jitter = max(eps,1e-12*max(diag(S)));
                    L = chol(S+jitter*eye(np),'lower');
                end
                cholC{s} = L;
                cols = find(any(D(rows,:),1));
                Z(rows,cols) = L\full(D(rows,cols)); z(rows) = L\y(rows);
            end
            beta = Z\z;
            effects = zeros(nR+nP,1); effects(columns) = beta;
            model.pipelineEffects(unique(pipeline)) = effects(nR+unique(pipeline));
            model.tau2 = tau2;
            fit = D*beta; estimate = fit; kept = 0;
            for s = unique(record)'
                rows = ix{s}; L = cholC{s};
                gain = tau2*(L'\(L\eye(numel(rows))));
                estimate(rows) = fit(rows)+gain*(y(rows)-fit(rows));
                kept = kept+trace(gain);
            end
            theta(ok) = exp(estimate); share = kept/m;
        end

        function [theta, info] = applyPrior(prior, Y, V, covariance)
            % Apply a prior fitted by twoLevel on OTHER participants. The
            % pipeline effects and interaction variance are frozen; only a
            % new recording's nuisance scale and posterior are estimated.
            % This permits independent-recording adaptation, not inference
            % on unseen trials. Unseen/disconnected candidates are not
            % compared using arbitrary gauges and retain their raw scores.
            assert(isequal(size(Y),size(V)), 'PipeCompare:Group', 'Objectives and variances must have the same size.');
            [nR,nP] = size(Y); theta = nan(size(Y));
            observed=isfinite(Y)&Y>0&isfinite(V)&V>=0;theta(observed)=Y(observed);
            assert(isfield(prior,'pipelineEffects') && numel(prior.pipelineEffects)==nP, ...
                'PipeCompare:Group', 'Prior pipeline order must match the new observations.');
            paired = nargin >= 4;
            if paired, assert(iscell(covariance) && numel(covariance)==nR, 'PipeCompare:Group', 'Supply one covariance matrix per recording.'); end
            info = struct('scope','new-recording-with-frozen-prior','fullyHeldOut',false, ...
                'priorRefitted',false,'assessed',false(nR,1),'reason',{repmat({''},nR,1)});
            for s = 1:nR
                p = find(isfinite(Y(s,:)) & Y(s,:)>0 & isfinite(V(s,:)) & V(s,:)>=0);
                if ~prior.varianceAssessed || ~isfinite(prior.tau2)
                    info.reason{s}='training interaction variance unassessed'; continue;
                end
                if isfield(prior,'recordsPerPipeline') && ...
                        (numel(p)<2 || any(prior.recordsPerPipeline(p)<prior.minimumTransferRecords))
                    info.reason{s}='insufficient independent training support for compared candidates';continue;
                end
                if isempty(p) || any(~isfinite(prior.pipelineEffects(p)))
                    info.reason{s}='candidate absent from training'; continue;
                end
                if numel(unique(prior.pipelineComponents(p)))~=1
                    info.reason{s}='training candidates have disconnected scales'; continue;
                end
                raw = diag(V(s,p));
                if paired
                    assert(isequal(size(covariance{s}),[nP nP]), 'PipeCompare:Group', 'Covariance matrix has the wrong size.');
                    raw = covariance{s}(p,p);
                    assert(all(isfinite(raw(:))) && max(abs(diag(raw)-V(s,p)'))<=1e-10*max([V(s,p),realmin]), ...
                        'PipeCompare:Group', 'Covariance diagonal must match V and be finite.');
                end
                C = groupCovariance(raw./(Y(s,p)'*Y(s,p)))+prior.tau2*eye(numel(p));
                [L,failed] = chol(C,'lower');
                if failed, L=chol(C+max(eps,1e-12*max(diag(C)))*eye(numel(p)),'lower'); end
                solve = @(x) L'\(L\x);
                y = log(Y(s,p))'; b = prior.pipelineEffects(p)'; one = ones(numel(p),1);
                offset = (one'*solve(y-b))/(one'*solve(one));
                fit = b+offset;
                theta(s,p) = exp(fit+prior.tau2*solve(y-fit))';
                info.assessed(s)=true;
            end
        end

        function print(G)
            T = G.table;
            pipecompare.utils.log('Group: %d recordings, %d pipelines run on all of them, %d passed on every recording.', ...
                numel(G.recordings), height(T), sum(strcmp(T.status, 'passed')));
            if ~isempty(G.leftOut)
                pipecompare.utils.log('Left out (no pipeline passed on them): %s', strjoin(G.leftOut, ', '));
            end
            if isempty(G.recommended), pipecompare.utils.log('%s', G.note); return; end
            k = G.recommended;
            pipecompare.utils.log(['Recommended for the group: pipeline %d, %s (noise %.3g, root mean square over the ', ...
                'recordings; at least %.0f%% of the trials kept in every recording; %d pipeline(s) not distinguished ', ...
                'from the best).'], k, G.labels{k}, T.objective(k), 100 * T.minRetention(k), sum(T.notDistinguished));
            pipecompare.utils.log('%s', pipecompare.Group.ownText(G));
        end

        function t = ownText(G)
            % Each recording's own pipeline, in words ('' when not made)
            t = '';
            if ~isfield(G, 'own') || isempty(G.own) || height(G.own) == 0, return; end
            O = G.own; ok = isfinite(O.id);
            if ~any(ok), return; end
            t = sprintf(['Each recording''s own pipeline, its estimate helped by the others (two-level empirical ', ...
                'Bayes; mean posterior gain %.0f%%): '], 100 * G.share);
            if ~isempty(G.recommended) && all(O.sameAsGroup(ok))
                t = [t sprintf('every recording gets the group''s pipeline %d, so one pipeline fits them all.', G.recommended)];
            else
                parts = arrayfun(@(i) sprintf('%s: %d', O.recording{i}, O.id(i)), find(ok)', 'UniformOutput', false);
                t = [t strjoin(parts, '; ') '.'];
                if ~isempty(G.recommended)
                    t = [t sprintf(' %d of %d use the group''s pipeline %d.', sum(O.sameAsGroup(ok)), sum(ok), G.recommended)];
                end
            end
        end

        function writeScript(G, file, i)
            % The recommended pipeline as a function file (PipeCompare.writeScript);
            % i: recording i's own pipeline (a row of G.own) instead.
            if nargin < 3
                assert(~isempty(G.recommended) && ~isempty(G.first), 'PipeCompare:Group', ...
                    'There is no recommended pipeline for this group.');
                k = G.recommended;
            else
                assert(isfield(G, 'own') && i >= 1 && i <= height(G.own) && isfinite(G.own.id(i)) && ~isempty(G.first), ...
                    'PipeCompare:Group', 'Recording %d has no pipeline of its own.', i);
                k = G.own.id(i);
            end
            idx = find(strcmp({G.first.leaves.key}, G.keys{k}), 1);
            pipecompare.PipeCompare.writeScript(G.first, idx, file);
        end
    end
end

function [D, components, active, component] = groupDesign(record,pipeline,nR,nP)
% One gauge per connected component; missing candidate patterns can split
% the recording/pipeline graph, so rank is not always nR+nP-1.
m = numel(record); rows = (1:m)';
D = sparse([rows;rows],[record;nR+pipeline],ones(2*m,1),m,nR+nP);
component = conncomp(graph(record,nR+pipeline,[],nR+nP));
active = false(1,nR+nP); active(unique([record;nR+pipeline])) = true;
labels = unique(component(nR+pipeline)); components = numel(labels);
for c = labels(:)'
    at = find(component(nR+pipeline)==c,1);
    active(nR+pipeline(at)) = false;
end
D = D(:,active);
end

function C = groupCovariance(C)
% Covariances from paired bootstraps may have tiny negative eigenvalues
% from rounding; reject larger violations instead of inventing precision.
scale = max(norm(C,'fro'),realmin);
assert(norm(C-C','fro')<=1e-10*scale, 'PipeCompare:Group', 'Covariance must be symmetric.');
C = (C+C')/2; [U,E] = eig(C); v = diag(E);
assert(min(v)>=-1e-10*scale, 'PipeCompare:Group', 'Covariance must be positive semidefinite.');
if any(v<0), C = U*diag(max(v,0))*U'; C = (C+C')/2; end
end

function s = saveBase()
% The base variables pop_pipecompare reads or sets, to put back afterwards.
s = struct();
for v = {'EEG', 'pipecompare_result'}
    if evalin('base', sprintf('exist(''%s'', ''var'')', v{1})), s.(v{1}) = evalin('base', v{1}); end
end
end

function restoreBase(s)
for v = {'EEG', 'pipecompare_result'}
    if isfield(s, v{1}), assignin('base', v{1}, s.(v{1})); else, evalin('base', ['clear ' v{1}]); end
end
end

