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
    %     Measurement error adds the mean of the squared SMEs, divided by
    %     the number of recordings, to the variance of a group mean, so
    %     this is what the group-level test loses to noise (and bias).
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
    %   pipecompare.Group.writeScript(G, file) writes the recommended
    %   pipeline as a function to run on each recording (its decisions are
    %   made on that recording's own data).

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
                'perRecording', nan(nU, nK), 'table', table(), 'recommended', [], 'note', '');
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
            if nU == 0
                G.note = 'No pipeline passed on any recording.';
            elseif isempty(F)
                G.note = 'No pipeline passed on every recording.';
            else
                acc = [];
                for j = 1:nU
                    s = S(u(j)); [~, loc] = ismember(keys(F), s.keys);
                    b = pipecompare.eval.Rank.bootstrap(s.cands(loc), s.ref, s.ranking.options, {s.ranking.objective});
                    if isempty(acc), acc = zeros(size(b{1})); end
                    acc = acc + b{1} .^ 2;
                end
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
        end

        function writeScript(G, file)
            % The recommended pipeline as a function file (PipeCompare.writeScript).
            assert(~isempty(G.recommended) && ~isempty(G.first), 'PipeCompare:Group', ...
                'There is no recommended pipeline for this group.');
            idx = find(strcmp({G.first.leaves.key}, G.keys{G.recommended}), 1);
            pipecompare.PipeCompare.writeScript(G.first, idx, file);
        end
    end
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
