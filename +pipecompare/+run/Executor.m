classdef Executor
    %EXECUTOR Run every legal pipeline of a plan on the current EEGLAB dataset.
    %
    %   result = pipecompare.run.Executor.run(plan, contract, opts)
    %   result = pipecompare.run.Executor.resume(checkpointDir)
    %
    %   - Reads the dataset that is current in EEGLAB (no loading step),
    %     prints its state and parsed EEG.history. The live dataset is never
    %     modified: everything runs on PipeCompare's own copy ("root").
    %   - Enumerates the legal pipelines (exhaustive, or explicitly
    %     never truncated or sampled) as a prefix tree and runs it depth first, so pipelines
    %     that share their first steps share that work.
    %   - Each EEGLAB command is printed and appended to the candidate's
    %     EEG.history.
    %   - Scores each finished pipeline (pipecompare.eval.Measure), checks that
    %     the signal survives (pipecompare.eval.Injection) and
    %     ranks (pipecompare.eval.Rank).
    %
    %   opts (all optional; see also pipecompare.eval.Rank.defaults):
    %     maxLeaves     refuse above this many pipelines (default 500)
    %     injectUv      injected amplitude (default 5 uV)
    %     dataUnit      'uV' (default) | 'V' (scaled to uV on the copy)
    %     checkpoint    folder: every finished candidate is saved there and
    %                   the search can be resumed with Executor.resume
    %     parallel      true: independent subtrees run on a parallel pool
    %     verbose       'normal' (commands, warnings, results) | 'full'
    %                   (everything EEGLAB prints)
    %     dryRun        enumerate and report, run nothing
    %     filterCheck   true (default): pipelines whose filters alone break
    %                   a signal-check limit are excluded before anything
    %                   runs (Injection.filterCheck)
    %     refineCutoffs false (default; true in the simple mode): when filter cutoffs are compared,
    %                   the recommended pipeline's are fine-tuned on a finer
    %                   grid with its decisions held fixed
    %                   (pipecompare.eval.Cutoff), and the pipeline with the
    %                   best ones is run and ranked with the others
    %     eogCheck      false (default; true in the simple mode): the
    %                   recommended pipeline also runs with the eye channels
    %                   regressed out (eogregress, before its components are
    %                   removed or its epochs made) and joins the ranking:
    %                   blinks ICA left in (e.g. no eye component found) are
    %                   taken out when that lowers the total error. Only with
    %                   an eye channel of its own (EOG) and when the pipeline
    %                   does not already regress them out (eogCheck)
    %     icaSeeds      0 (default; 2 in the simple mode): how many times
    %                   the pipelines that decide the recommendation are
    %                   run again with ICA started differently (icaSeedCheck);
    %                   result.icaSeeds says whether the recommendation stays
    %     signalTemplate 'bump' (default) | 'data' (the simple mode's
    %                   default): the time course of the known signal of the
    %                   signal check, a bump centred in each window or this
    %                   recording's own average waveform at the ROI where it
    %                   is reliable (pipecompare.eval.Injection)
    %     progress      function stop = f(n), called with the number of
    %                   candidates just evaluated (0 = only asking); when it
    %                   returns true, no further step runs and the
    %                   candidates not run are 'failed' ("not run"), so the
    %                   ranking covers the finished ones. Serial runs only.
    %                   A callback that takes two arguments is also told
    %                   each step as it starts: f(0, inst), inst being the
    %                   plan's step (type, params, label).
    %     artifactVariance 'pooled' (default) | 'legacy-total' |
    %                   'heteroscedastic': how the artifact contribution
    %                   to the paired bootstrap variance is pooled
    %                   (pipecompare.eval.Measure.defaults); run() and
    %                   execute() share this one default.

    properties (Constant)
        % The one default for opts.artifactVariance, shared by run() and
        % execute(): a resumed search must not silently change its
        % variance model, and it matches pipecompare.eval.Measure.defaults.
        ArtifactVarianceDefault = 'pooled';
    end

    methods (Static)
        function result = run(plan, contract, opts)
            if nargin < 3, opts = struct(); end
            opts = pipecompare.utils.withDefaults(opts, struct('dryRun', false, 'injectUv', 5, 'filterCheck', true, ...
                'refineCutoffs', false, 'eogCheck', false, 'icaSeeds', 0, 'signalTemplate', 'bump', 'dataUnit', 'uV', 'checkpoint', '', 'maxSeconds', Inf, 'parallel', false, 'verbose', 'normal', 'progress', []));
            validateattributes(opts.maxSeconds,{'numeric'},{'scalar','real','positive','nonnan'});
            assert(isinf(opts.maxSeconds)||~isempty(opts.checkpoint),'PipeCompare:Options', ...
                'A time-limited search requires a checkpoint folder.');
            assert(isinf(opts.maxSeconds)||~opts.parallel,'PipeCompare:Options', ...
                'Time slices use serial prefix sharing; run independent recordings in parallel.');
            opts.artifactVariance=pipecompare.utils.fieldOr(opts,'artifactVariance',pipecompare.run.Executor.ArtifactVarianceDefault);
            assert(any(strcmp(opts.artifactVariance,{'pooled','legacy-total','heteroscedastic'})),'PipeCompare:Options','Unknown artifact variance model.');
            assert(isnumeric(opts.icaSeeds) && isscalar(opts.icaSeeds) && opts.icaSeeds >= 0 && opts.icaSeeds == round(opts.icaSeeds), ...
                'PipeCompare:Options', 'icaSeeds must be a whole number >= 0.');
            assert(any(strcmp(opts.signalTemplate, {'bump', 'data'})), 'PipeCompare:Options', ...
                'signalTemplate must be ''bump'' or ''data''.');
            [EEG, live] = pipecompare.live.Session.current();
            assert(~isempty(EEG), 'PipeCompare:NoDataset', 'No dataset is loaded in EEGLAB.');
            if ~live.stored
                pipecompare.utils.log('WARNING: base EEG differs from ALLEEG(CURRENTSET) (changed but not stored); using base EEG.');
            end
            state = pipecompare.live.DataState.fromEEG(EEG);
            pipecompare.live.DataState.print(state);
            contract.validate(state);
            plan.print();
            [leaves, tree, rep] = plan.enumerate(state, contract, opts);
            pipecompare.run.Executor.printReport(rep);
            labels = arrayfun(@(l) l.key, leaves, 'UniformOutput', false);
            result = struct('plan', plan, 'contract', contract, 'options', opts, 'state', state, ...
                'rootFingerprint', live.fingerprint, 'leaves', leaves, 'tree', tree, 'report', rep, ...
                'labels', {labels}, 'cands', [], 'ranking', [], 'marginal', [], 'robustness', [], 'root', [], 'ref', [], ...
                'rootComs', {{}}, 'versions', versions(), 'identity', '');
            if opts.dryRun, return; end
            [root, rootComs] = pipecompare.run.Executor.prepareRoot(EEG, contract, opts);
            result.root = root; result.rootComs = rootComs;
            result = pipecompare.run.Executor.execute(result, false(numel(leaves), 1));
        end

        function result = resume(dirName, runtime)
            % Continue an interrupted search from its checkpoint folder.
            if nargin<2,runtime=struct();end
            % Round 8 strict checkpoint: the default resume path is
            % unchanged (version warning + pooled default). Only when
            % the caller explicitly passes strictCheckpoint=true does
            % resume refuse a version mismatch or a missing
            % artifactVariance, so old checkpoints cannot be silently
            % mixed under a new variance model.
            assert(all(ismember(fieldnames(runtime),{'maxSeconds','strictCheckpoint'})),'PipeCompare:Options', ...
                'Resume overrides may only change maxSeconds and strictCheckpoint.');
            strictCheckpoint = pipecompare.utils.fieldOr(runtime, 'strictCheckpoint', false);
            M = load(fullfile(dirName, 'manifest.mat'));
            R = load(fullfile(dirName, 'root.mat'));
            assert(isfield(M, 'identity') && isfield(R, 'identity') && strcmp(M.identity, R.identity), ...
                'PipeCompare:Checkpoint', ['%s: the starting dataset (root.mat) and the search description ', ...
                '(manifest.mat) do not belong to the same search (or were written by an older PipeCompare).'], dirName);
            result = M.result; result.root = R.root;
            result.identity = M.identity;
            w = what(dirName); result.options.checkpoint = w(1).path;   % absolute, as writeManifest stores it
            assert(strcmp(searchIdentity(result), M.identity), 'PipeCompare:Checkpoint', ...
                '%s: the stored starting dataset no longer matches the search identity.', dirName);
            v = versions();
            if strictCheckpoint
                assert(isequal(v.eeglab, result.versions.eeglab) && isequal(v.matlab, result.versions.matlab), ...
                    'PipeCompare:Checkpoint', 'Strict resume refuses a version mismatch (%s / %s now, %s / %s then).', ...
                    v.matlab, v.eeglab, result.versions.matlab, result.versions.eeglab);
                assert(isfield(result.options, 'artifactVariance') && ~isempty(result.options.artifactVariance), ...
                    'PipeCompare:Checkpoint', 'Strict resume refuses a checkpoint without an explicit artifactVariance (it must not silently default to pooled).');
            elseif ~isequal(v.eeglab, result.versions.eeglab) || ~isequal(v.matlab, result.versions.matlab)
                pipecompare.utils.log('WARNING: software versions differ from the interrupted run (%s / %s now, %s / %s then).', ...
                    v.matlab, v.eeglab, result.versions.matlab, result.versions.eeglab);
            end
            done = false(numel(result.leaves), 1);
            files = dir(fullfile(dirName, 'leaf_*.mat'));
            files=files(~cellfun(@isempty,regexp({files.name},'^leaf_[0-9]{6}\.mat$','once')));
            prev = repmat(emptyCand(), numel(result.leaves), 1);
            for f = 1:numel(files)
                d = load(fullfile(dirName, files(f).name));
                if isfield(d, 'id') && isscalar(d.id) && d.id > numel(result.leaves)
                    % An adaptive candidate (refineCutoffs/eogCheck): saved
                    % for audit and replay, but the manifest predates it, so
                    % resume cannot restore it by id. Skipping it here keeps
                    % resume working; the resumed search re-derives the
                    % adaptive candidates when it reaches those phases.
                    % The identity check below still applies, unchanged, to
                    % every file that is actually loaded.
                    pipecompare.utils.log('Skipping adaptive leaf file %s (id %d beyond the manifest leaves; saved for audit, re-derived on resume).', files(f).name, d.id);
                    continue;
                end
                assert(isfield(d, 'identity') && strcmp(d.identity, M.identity) && d.id <= numel(result.leaves) && ...
                    strcmp(d.key, result.leaves(d.id).key), 'PipeCompare:Checkpoint', ...
                    'Checkpoint file %s does not belong to this search.', files(f).name);
                c = emptyCand();   % a file from an earlier version may lack newer fields
                for n = fieldnames(d.cand)', c.(n{1}) = d.cand.(n{1}); end
                prev(d.id) = c; done(d.id) = true;
            end
            pipecompare.utils.log('Resuming: %d of %d candidates already evaluated.', sum(done), numel(done));
            result.prevCands = prev;
            if isfield(result.options, 'stopAfter'), result.options = rmfield(result.options, 'stopAfter'); end
            if isfield(runtime,'maxSeconds')
                validateattributes(runtime.maxSeconds,{'numeric'},{'scalar','real','positive','nonnan'});
                assert(isinf(runtime.maxSeconds)||~result.options.parallel,'PipeCompare:Options','Time slices require serial prefix sharing.');
                result.options.maxSeconds=runtime.maxSeconds;
            end
            result = pipecompare.run.Executor.execute(result, done);
        end

        function result = execute(result, done)
            opts = result.options;
            % Single default (was 'legacy-total' here): a checkpoint made
            % before this option existed now resumes under 'pooled', the
            % same model run() uses, rather than silently changing models
            % between the two entry points.
            opts.artifactVariance=pipecompare.utils.fieldOr(opts,'artifactVariance',pipecompare.run.Executor.ArtifactVarianceDefault);
            result.options=opts;
            contract = result.contract; leaves = result.leaves; tree = result.tree;
            root = result.root;
            ref = pipecompare.eval.Measure.reference(root, contract,opts);
            result.ref = ref;
            pipecompare.utils.log('Reference trials per condition: %s', ...
                strjoin(arrayfun(@(c) sprintf('%s=%d', ref.names{c}, ref.n(c)), 1:numel(ref.n), 'UniformOutput', false), ', '));
            env = struct('tree', tree, 'leaves', leaves, 'contract', contract, 'ref', ref, 'opts', opts, ...
                'nbchan', root.nbchan + numel(pipecompare.utils.fieldOr(root.etc.pipecompare, 'preRemoved', [])), 'done', done, 'checkpoint', opts.checkpoint, 'tStart', tic, ...
                'truth', [], 'identity', '');
            % signal check: a copy holding only a known signal goes through every
            % candidate with the same operations and decisions as the real data
            [S, env.truth] = pipecompare.eval.Injection.prepare(root, contract, ref, opts);
            pipecompare.utils.log(['Signal check: a known %g uV signal is carried through every candidate with the ', ...
                'same decisions (matched-decision injection).'], opts.injectUv);
            % how much each pipeline's filters make the gain depend on the
            % signal's shape (part of the bias of the total error)
            try
                env.spread = pipecompare.eval.Injection.shapeSpreads(S, contract, env.truth, leaves);
            catch ME
                    if startsWith(ME.identifier,'PipeCompare:Checkpoint') || strcmp(ME.identifier,'PipeCompare:Interrupted'),rethrow(ME);end
                env.spread = {};
                pipecompare.utils.log('Signal shapes: not checked (%s); the bias leaves out what the filters do to other shapes.', ME.message);
            end
            bi = pipecompare.utils.fieldOr(ref, 'blinkInfo', struct('channels', {{}}, 'share', []));
            if ~isempty(pipecompare.utils.fieldOr(bi, 'kindShare', []))
                for j = find(isfinite(bi.share))
                    pipecompare.utils.log(['Artifacts in the %s window (what the trials with them add to each mean counts ', ...
                        'in the total error): blinks %s, eye movements %s, muscle %s; %.0f%% of the trials have one.'], ...
                        contract.components(j).name, ...
                        pipecompare.utils.ternary(isempty(bi.channels), 'not checked (no eye channel)', ...
                            sprintf('%.0f%% (%s)', 100 * bi.kindShare(j, 1), strjoin(bi.channels, ', '))), ...
                        pipecompare.utils.ternary(isempty(bi.horizontal), 'not checked (no left-right pair)', ...
                            sprintf('%.0f%% (%s)', 100 * bi.kindShare(j, 2), bi.horizontal)), ...
                        sprintf('%.0f%%', 100 * bi.kindShare(j, 3)), 100 * bi.share(j));
                end
            end
            if strcmp(contract.analysis, 'erp') && isempty(bi.channels)
                pipecompare.utils.log('Blinks: no eye channel (EOG, or Fp1/Fp2/FPz) that is not measured, so what blinks add is not estimated.');
            end
            result.signalTemplate = [];
            if isfield(env.truth, 'shape') && ~isempty(env.truth.shape)
                result.signalTemplate = env.truth.shape.info;
                for j = 1:numel(contract.components)
                    if result.signalTemplate.used(j)
                        pipecompare.utils.log('Signal check, %s: the shape of this recording''s own average (split-half reliability %.2f).', ...
                            contract.components(j).name, result.signalTemplate.reliability(j));
                    else
                        pipecompare.utils.log('Signal check, %s: the standard bump (this recording''s average: %s).', ...
                            contract.components(j).name, result.signalTemplate.why{j});
                    end
                end
            end
            % depth-first execution holds the data and the signal copy once per
            % level of the current branch (the price of sharing prefixes)
            % (a level after a resample step holds data smaller by fs_new/fs)
            b = whos('root'); gb = b.bytes / 2^30;
            levels = max([1 arrayfun(@(l) dataLevels(l.path, root.srate), leaves)]);
            if gb * 2 * levels > 2
                pipecompare.utils.log(['Memory: up to about %.1f GB (%.2f GB of starting data, 2 copies, %.1f data-sized ', ...
                    'levels on the deepest branch). Resampling early in the plan reduces it.'], gb * 2 * levels, gb, levels);
            end
            if ~isempty(opts.checkpoint) && ~any(done)
                result = pipecompare.run.Executor.writeManifest(result);
                env.identity = result.identity; env.checkpoint = result.options.checkpoint;
            elseif isfield(result, 'identity')
                env.identity = result.identity;
            end
            ctx0 = struct('contract', contract, 'highpass', rootHighpass(result.state), 'ref', ref, 'rank', opts);
            acc0 = struct('interpolated', {{}}, 'icsRemoved', 0, 'rejected', 0, 'coms', {{}}, 'seconds', 0, 'unmatched', {{}}, ...
                'ica', {{}}, 'overLimit', noOverLimit(), 'steps', {{}});
            % pipelines whose filters alone distort the known signal are
            % excluded before anything runs (filters are linear and
            % time-invariant: that does not depend on the data)
            pre = repmat(emptyCand(), 0, 1);
            if pipecompare.utils.fieldOr(opts, 'filterCheck', true) && ~isempty(S)
                why = pipecompare.eval.Injection.filterCheck(S, contract, env.truth, leaves, opts);
                out = find(~cellfun(@isempty, why(:)) & ~done(:))';
                for li = out
                    c = emptyCand(); c.id = li; c.key = leaves(li).key; c.stratum = leaves(li).stratum;
                    c.status = 'excluded'; c.message = why{li};
                    pre(end+1, 1) = c; env.done(li) = true; tick(env, 1); %#ok<AGROW>
                    saveLeaf(env, c);   % an excluded candidate is still part of the ranking input
                    pipecompare.utils.log('candidate %d excluded before running: %s', li, c.message);
                end
                if ~isempty(out)
                    pipecompare.utils.log(['Filter check: %d pipeline(s) excluded before running, because their filters ', ...
                        'alone distort the known signal (e.g. %s).'], numel(out), why{out(1)});
                end
            end
            % pipelines with no step at all (every slot chose 'none') are the starting copy itself
            for li = tree(1).leaves
                if env.done(li), continue; end
                c = evaluateLeaf(env, li, root, S, acc0); pre(end+1, 1) = c; saveLeaf(env, c); %#ok<AGROW>
                tick(env, 1);
            end
            if opts.parallel && canParallel()
                cands = pipecompare.run.Executor.runParallel(env, root, S, ctx0, acc0);
            else
                if opts.parallel, pipecompare.utils.log('Parallel evaluation unavailable (Parallel Computing Toolbox/pool); running serially.'); end
                cands = pipecompare.run.Executor.runSubtree(env, 1, root, S, ctx0, acc0);
            end
            cands = [pre; cands];
            allc = repmat(emptyCand(), numel(leaves), 1);
            if isfield(result, 'prevCands'), allc = result.prevCands; result = rmfield(result, 'prevCands'); end
            for k = 1:numel(cands), allc(cands(k).id) = cands(k); end
            notRun = find([allc.id] == 0);   % only after a stop (opts.progress)
            for li = notRun
                allc(li).id = li; allc(li).key = leaves(li).key; allc(li).stratum = leaves(li).stratum;
                allc(li).status = 'failed'; allc(li).message = 'not run: the search was stopped';
            end
            if ~isempty(notRun)
                pipecompare.utils.log('Stopped: %d of %d candidates were not run; the ranking covers the others.', ...
                    numel(notRun), numel(allc));
            end
            result.notRun = numel(notRun);
            R = pipecompare.eval.Rank.run(allc, ref, opts, levelsOf(result.leaves));
            if pipecompare.utils.fieldOr(opts, 'refineCutoffs', false) && isempty(notRun)
                [result, allc, added] = refineCutoffs(result, env, root, S, ctx0, acc0, allc, R);
                if added, R = pipecompare.eval.Rank.run(allc, ref, opts, levelsOf(result.leaves)); end
            end
            if pipecompare.utils.fieldOr(opts, 'eogCheck', false) && isempty(notRun)
                env.leaves = result.leaves;   % with a pipeline added by fine-tuning the cutoffs
                [result, allc, added] = eogCheck(result, env, root, S, ctx0, acc0, allc, R);
                if added, R = pipecompare.eval.Rank.run(allc, ref, opts, levelsOf(result.leaves)); end
            end
            result.icaSeeds = [];
            if pipecompare.utils.fieldOr(opts, 'icaSeeds', 0) > 0 && isempty(notRun)
                env.leaves = result.leaves;   % with a pipeline added by fine-tuning the cutoffs
                result.icaSeeds = icaSeedCheck(env, root, S, ctx0, acc0, allc, R, ref, opts);
            end
            result.cands = allc;
            result.ranking = R;
            if isempty(notRun)
                saveRanking(env, R, result.labels);
            end
            result.marginal = pipecompare.eval.Rank.marginal(R, result.leaves, result.report.searched);
            result.robustness = pipecompare.eval.Rank.robustness(allc, R, result.leaves, result.report.searched, ref);
            pipecompare.eval.Rank.print(R, result.labels);
            notes = unique([allc.unmatched]);
            if ~isempty(notes)
                pipecompare.utils.log(['Signal check note: %s step(s) were re-run rather than replayed on the injected copy; ', ...
                    'their signal effect is measured but not decision-matched.'], strjoin(notes, ', '));
            end
            if ~isempty(result.marginal) && height(result.marginal) > 0
                pipecompare.utils.log('Effect of each searched choice (medians over evaluated candidates; best among feasible):');
                disp(result.marginal);
            end
            if ~isempty(result.robustness) && height(result.robustness) > 0
                pipecompare.utils.log(['Multiverse summary: each measure over the feasible pipelines, and the searched ', ...
                    'choice that accounts for most of its spread (share = eta^2). Sensitivity only; not used for the ranking.']);
                disp(result.robustness);
            end
            pipecompare.utils.log('Finished in %.1f s.', toc(env.tStart));
            if isfield(result.options, 'progress'), result.options.progress = []; end   % do not keep the caller's window
        end

        function out = runSubtree(env, node, E, S, ctx, acc)
            % Depth-first execution of the subtree below node. Returns one
            % record per evaluated leaf. No shared state: safe for parfor.
            out = repmat(emptyCand(), 0, 1);
            tree = env.tree;
            labels = [];   % ICLabel classification of E, shared by sibling icremove steps
            for child = tree(node).children
                under = leavesUnder(tree, child);
                if all(env.done(under)), continue; end
                in = tree(child).inst;
                if tick(env, 0, in), continue; end   % stopped: nothing more runs
                t0 = tic;
                try
                    pipecompare.utils.log('step %s', in.label);
                    cx = ctx;
                    if strcmp(in.type, 'icremove'), cx.iclabel = labels; end
                    cx = choiceContext(env, child, cx, S);
                    if strcmp(in.type,'ica') && ~isempty(env.opts.checkpoint)
                        [E2,coms,info]=pipecompare.run.IcaCheckpoint.run(E,in,cx, ...
                            env.opts.checkpoint,env.identity,@()runStep(in,E,cx,env.opts.verbose));
                    else
                        [E2, coms, info] = runStep(in, E, cx, env.opts.verbose);
                    end
                    if isfield(info, 'iclabel'), labels = info.iclabel; end
                    for q = 1:numel(coms)
                        fprintf('    %s\n', coms{q});
                        E2 = eeg_hist(E2, coms{q});
                    end
                catch ME
                    if startsWith(ME.identifier,'PipeCompare:Checkpoint') || strcmp(ME.identifier,'PipeCompare:Interrupted'),rethrow(ME);end
                    out = [out; failUnder(env, under, in, ME)]; %#ok<AGROW>
                    continue;
                end
                [S2, unmatched] = replayCopy(in, S, info, ctx);
                [ctx2, acc2] = advance(in, ctx, acc, info, coms, unmatched, toc(t0));
                for li = tree(child).leaves
                    if env.done(li), continue; end
                    c = evaluateLeaf(env, li, E2, S2, acc2);
                    out(end+1, 1) = c; %#ok<AGROW>
                    saveLeaf(env, c);
                    tick(env, 1);
                end
                out = [out; pipecompare.run.Executor.runSubtree(env, child, E2, S2, ctx2, acc2)]; %#ok<AGROW>
            end
        end

        function out = runParallel(env, root, S, ctx, acc)
            % Run the shared single-path trunk serially, then the independent
            % subtrees below it on pool workers. Each worker receives its
            % own copies of the datasets (MATLAB value semantics), so no
            % candidate can see another candidate's data.
            tree = env.tree; node = 1; E = root;
            out = repmat(emptyCand(), 0, 1);
            while isscalar(tree(node).children)
                child = tree(node).children;
                in = tree(child).inst;
                if tick(env, 0, in), return; end   % stopped: nothing more runs
                pipecompare.utils.log('step %s (shared trunk)', in.label);
                t0 = tic;
                try
                    [E, coms, info] = runStep(in, E, choiceContext(env, child, ctx, S), env.opts.verbose);
                    for q = 1:numel(coms), fprintf('    %s\n', coms{q}); E = eeg_hist(E, coms{q}); end
                catch ME   % same handling as runSubtree: the candidates below fail, the search goes on
                    out = [out; failUnder(env, leavesUnder(tree, child), in, ME)];
                    return;
                end
                [S, unmatched] = replayCopy(in, S, info, ctx);
                [ctx, acc] = advance(in, ctx, acc, info, coms, unmatched, toc(t0));
                for li = tree(child).leaves
                    if env.done(li), continue; end
                    c = evaluateLeaf(env, li, E, S, acc); out(end+1, 1) = c; saveLeaf(env, c); %#ok<AGROW>
                    tick(env, 1);
                end
                node = child;
            end
            kids = tree(node).children;
            if isempty(kids), return; end
            pipecompare.utils.log('Parallel: %d independent subtrees on the pool.', numel(kids));
            % The caller's progress window cannot be reached from the workers:
            % they report to it through a DataQueue (each finished pipeline,
            % and each check for Stop), and its Stop reaches them as a file
            % they look for (the pool's workers run on this computer).
            prog = env.opts.progress; env.opts.progress = [];
            stopFile = [tempname '_pipecompare_stop'];
            cleanStop = onCleanup(@() deleteIfThere(stopFile)); %#ok<NASGU>
            if ~isempty(prog)
                try
                    q = parallel.pool.DataQueue;
                    afterEach(q, @(args) relayProgress(prog, args, stopFile));
                    env.opts.progress = @(varargin) workerProgress(q, stopFile, varargin);
                catch ME
                    if startsWith(ME.identifier,'PipeCompare:Checkpoint') || strcmp(ME.identifier,'PipeCompare:Interrupted'),rethrow(ME);end
                    pipecompare.utils.log('Progress and Stop are not available while the pool runs (%s).', ME.message);
                end
            end
            parts = cell(1, numel(kids));
            parfor k = 1:numel(kids)
                sub = env;
                sub.tree(node).children = kids(k);
                parts{k} = pipecompare.run.Executor.runSubtree(sub, node, E, S, ctx, acc);
            end
            out = [out; vertcat(parts{:})];
        end

        function [root, coms] = prepareRoot(EEG, contract, opts)
            % PipeCompare's own copy of the live dataset; the original is untouched.
            if nargin < 3, opts = struct('dataUnit', 'uV'); end
            root = EEG; coms = {};
            if ~isempty(root.chanlocs), root.chanlocs = root.chanlocs(:)'; end   % a column after some importers (pop_biosig)
            function add(c), root = eeg_hist(root, c); coms{end+1} = c; end
            if isfield(opts, 'dataUnit') && strcmp(opts.dataUnit, 'V')
                % the script line does exactly what is done here (ICA weights
                % too, so existing activations keep their scale)
                cmd = ['EEG.data = EEG.data * 1e6; if ~isempty(EEG.icaweights), EEG.icaweights = EEG.icaweights / 1e6; ', ...
                    'EEG.icawinv = []; EEG.icaact = []; end; EEG = eeg_checkset(EEG); % PipeCompare: volts -> microvolts'];
                root = evalWithEEG(root, cmd);
                add(cmd);
            end
            r = round(root.srate);
            if root.srate ~= r && abs(root.srate - r) <= 1e-9 * r
                pipecompare.utils.log('EEG.srate is %.17g (floating-point residue); using %d on PipeCompare''s copy.', root.srate, r);
                root.srate = r; root = eeg_checkset(root);
                add(sprintf('EEG.srate = %d; EEG = eeg_checkset(EEG); %% PipeCompare: exact sampling rate', r));
            end
            if ~isfield(root, 'urevent') || isempty(root.urevent)
                root = eeg_checkset(root, 'makeur');
                add('EEG = eeg_checkset(EEG, ''makeur''); % PipeCompare: trial identities');
            end
            if contract.isSegmented() && root.trials == 1
                % consecutive analysis segments, marked with EEGLAB's own function
                % (events and urevents, so each segment is paired across candidates)
                cmd = sprintf(['EEG = eeg_regepochs(EEG, ''recurrence'', %g, ''limits'', [0 %g], ''eventtype'', ''%s'', ', ...
                    '''extractepochs'', ''off''); %% PipeCompare: %g s analysis segments'], contract.segment, contract.segment, ...
                    pipecompare.eval.Contract.SegmentEvent, contract.segment);
                [~, root] = evalc('evalWithEEG(root, cmd)');
                add(cmd);
            end
            % a trial list left in the dataset by an earlier search (e.g. an
            % adopted candidate) must never restrict this one; datasets adopted
            % with NeuroQC (<= 0.7) keep it under EEG.etc.neuroqc
            for f = {'pipecompare', 'neuroqc'}
                if isfield(root, 'etc') && isstruct(root.etc) && isfield(root.etc, f{1}) && isfield(root.etc.(f{1}), 'eligibleUrevents')
                    root.etc.(f{1}) = rmfield(root.etc.(f{1}), 'eligibleUrevents');
                    add(sprintf('EEG.etc.%s = rmfield(EEG.etc.%s, ''eligibleUrevents''); %% PipeCompare: drop an earlier trial rule', ...
                        f{1}, f{1}));
                end
            end
            elig = eligibleTrials(root, contract);
            if ~isempty(elig)
                root.etc.pipecompare.eligibleUrevents = elig;
                add(sprintf('EEG.etc.pipecompare.eligibleUrevents = %s; %% PipeCompare: trial rule %s (urevent ids)', ...
                    mat2str(elig(:)'), contract.trials.mode));
                pipecompare.utils.log('Trial rule %s: %d eligible time-locking events.', contract.trials.mode, numel(elig));
            end
            root.etc.pipecompare.rootChanlocs = root.chanlocs;
            add('EEG.etc.pipecompare.rootChanlocs = EEG.chanlocs; % PipeCompare: montage before the plan');
            pre = pipecompare.live.DataState.fromEEG(root).restorableChannels;
            if ~isempty(pre)
                root.etc.pipecompare.preRemoved = root.chaninfo.removedchans(pre);
                add(sprintf(['EEG.etc.pipecompare.preRemoved = EEG.chaninfo.removedchans(%s); ', ...
                    '%% PipeCompare: channels removed before the plan, restorable by interpolation'], mat2str(pre)));
                pipecompare.utils.log('Channels removed before PipeCompare that a restore step can interpolate: %s', ...
                    strjoin({root.etc.pipecompare.preRemoved.labels}, ', '));
            end
        end

        function root = rootOf(result)
            % The starting dataset of a search. With a checkpoint folder it
            % is not kept in the result (long searches start from large
            % data); it is read back from root.mat, checked against the
            % search identity.
            root = result.root;
            if ~isempty(root), return; end
            d = '';
            if isfield(result.options, 'checkpoint'), d = result.options.checkpoint; end
            f = fullfile(d, 'root.mat');
            assert(~isempty(d) && isfile(f), 'PipeCompare:Checkpoint', ['The starting dataset of this result is not in ', ...
                'memory and its checkpoint file %s is missing.'], f);
            R = load(f);
            assert(isfield(R, 'identity') && isfield(result, 'identity') && strcmp(R.identity, result.identity), ...
                'PipeCompare:Checkpoint', '%s belongs to another search; the starting dataset cannot be restored.', f);
            root = R.root;
        end

        function EEG = replay(result, idx, recordGlobal)
            % Re-run candidate idx from the stored starting dataset.
            if nargin < 3, recordGlobal = false; end
            EEG = pipecompare.run.Executor.rootOf(result);
            path = result.leaves(idx).path;
            ctx = struct('contract', result.contract, 'highpass', rootHighpass(result.state), 'ref', pipecompare.utils.fieldOr(result, 'ref', []), ...
                'rank', result.options);
            if recordGlobal
                for q = 1:numel(result.rootComs), eegh(result.rootComs{q}); end
            end
            % ICA decompositions the search computed on this path (runica is
            % deterministic here, so fitting again would give the same ones)
            ica = {}; nIca = 0;
            if isfield(result, 'cands') && numel(result.cands) >= idx && isfield(result.cands, 'ica')
                ica = result.cands(idx).ica;
            end
            for k = 1:numel(path)
                in = path{k};
                if strcmp(in.type, 'ica'), nIca = nIca + 1; end
                if strcmp(in.type, 'ica') && nIca <= numel(ica)
                    pipecompare.utils.log('step %s: the decomposition from the search is reused', in.label);
                    d = ica{nIca};
                    EEG.icaweights = d.icaweights; EEG.icasphere = d.icasphere;
                    EEG.icachansind = d.icachansind; EEG.icawinv = []; EEG.icaact = [];
                    EEG = eeg_checkset(EEG);
                    coms = d.coms;
                else
                    pipecompare.utils.log('step %s', in.label);
                    if strcmp(in.type, 'icremove') && ischar(in.params.threshold) && isfield(result, 'cands') && ...
                            numel(result.cands) >= idx && numel(result.cands(idx).steps) >= k && ...
                            isfinite(pipecompare.utils.fieldOr(result.cands(idx).steps{k}, 'thresholdChosen', NaN))
                        % the threshold the search chose on these same data
                        in.params.threshold = result.cands(idx).steps{k}.thresholdChosen;
                    end
                    cx = ctx;
                    if isRepairChoice(in) && isfield(result, 'cands') && numel(result.cands) >= idx && numel(result.cands(idx).steps) >= k
                        % the bad channels the search interpolated on these same data
                        cx.chosenChannels = result.cands(idx).steps{k}.interpolated;
                    end
                    if strcmp(in.type, 'icremove') && isequal(in.params.threshold, 'sme') && isfield(result, 'cands') && ...
                            numel(result.cands) >= idx && numel(result.cands(idx).steps) >= k
                        % the components the search chose on these same data
                        cx.chosenComps = result.cands(idx).steps{k}.comps;
                    end
                    [EEG, coms] = runStep(in, EEG, cx, result.options.verbose);
                end
                for q = 1:numel(coms)
                    fprintf('    %s\n', coms{q});
                    % eeg_hist always appends; eegh(com, EEG) would skip a command equal to
                    % the previous session command, losing it from the dataset history
                    EEG = eeg_hist(EEG, coms{q});
                    if recordGlobal, eegh(coms{q}); end
                end
                ctx = advance(in, ctx, [], struct(), coms, {}, 0);
            end
            [EEG, com] = pipecompare.run.Executor.selectEligible(EEG, result.contract);
            if ~isempty(com)
                fprintf('    %s\n', com);
                EEG = eeg_hist(EEG, com);
                if recordGlobal, eegh(com); end
            end
        end

        function EEG = apply(EEG, steps, contract, opts)
            % Run steps (rows of type, params) on any dataset as the search
            % runs them, so each step decides from these data (bad
            % channels, components, rejected epochs); used by saved
            % scripts. The preparation is the search's own (unit judged
            % from the amplitude scale); every command goes to EEG.history.
            % opts.signalTemplate: the known signal's shape for the
            % choices that use the signal check, as in the search.
            if nargin < 3 || isempty(contract), contract = pipecompare.eval.Contract(); end
            if nargin < 4, opts = struct(); end
            state = pipecompare.live.DataState.fromEEG(EEG);
            EEG = pipecompare.run.Executor.prepareRoot(EEG, contract, struct('dataUnit', state.unitGuess));
            ctx = struct('contract', contract, 'highpass', rootHighpass(state));
            chosen = @(f) cellfun(@(p) isfield(p, f) && ischar(p.(f)), steps(:, 2));
            repair = cellfun(@(t, p) strcmp(t, 'badchannels') && strcmp(pipecompare.utils.fieldOr(p, 'repair', 'all'), 'error'), ...
                steps(:, 1), steps(:, 2));
            if any(chosen('uv') | chosen('threshold') | repair)
                % a limit or threshold chosen from the data is judged against these data's own trials
                ctx.ref = pipecompare.eval.Measure.reference(EEG, contract,opts);
            end
            S = []; truth = [];
            if any(chosen('threshold') | repair)
                % an ICLabel threshold or bad channels chosen from the data use the signal check
                [S, truth] = pipecompare.eval.Injection.prepare(EEG, contract, ctx.ref, opts);
            end
            insts = arrayfun(@(k) struct('type', steps{k, 1}, 'params', steps{k, 2}), 1:size(steps, 1), 'UniformOutput', false);
            for k = 1:numel(insts)
                in = insts{k};
                cx = ctx;
                if strcmp(in.type, 'icremove') && ischar(in.params.threshold)
                    cx.after = insts(k+1:end); cx.signal = struct('S', S, 'truth', truth);
                elseif isRepairChoice(in)
                    cx.after = linearAfter(insts(k+1:end)); cx.signal = struct('S', S, 'truth', truth);
                end
                [EEG, coms, info] = pipecompare.run.Steps.run(in, EEG, cx);
                for q = 1:numel(coms), EEG = eeg_hist(EEG, coms{q}); end
                if ~isempty(S), S = replayCopy(in, S, info, ctx); end
                ctx = advance(in, ctx, [], struct(), coms, {}, 0);
            end
            [EEG, com] = pipecompare.run.Executor.selectEligible(EEG, contract);
            if ~isempty(com), EEG = eeg_hist(EEG, com); end
            EEG = storeWeightedAverages(EEG, contract);
        end

        function [elig, EEG] = eligibleUrevents(EEG, contract)
            % urevent ids the contract's trial rule keeps ([] = all trials),
            % for display before a search (EEG is a copy; urevents are
            % built the same way the search builds them).
            if ~isfield(EEG, 'urevent') || isempty(EEG.urevent)
                [~, EEG] = evalc('eeg_checkset(EEG, ''makeur'')');
            end
            elig = eligibleTrials(EEG, contract);
        end

        function [EEG, com] = selectEligible(EEG, contract)
            % Output only the trials the trial rule keeps. The rule decides
            % which trials are scored; the data a candidate hands on must
            % hold the same trials. Epochs time-locked to a condition event
            % outside the rule are removed (pop_select, in EEG.history).
            com = '';
            if ~isfield(EEG, 'etc') || ~isstruct(EEG.etc) || ~isfield(EEG.etc, 'pipecompare') || ...
                    ~isfield(EEG.etc.pipecompare, 'eligibleUrevents')
                return;
            end
            if EEG.trials == 1 && (~isfield(EEG, 'epoch') || isempty(EEG.epoch))
                pipecompare.utils.log(['Trial rule: the output is continuous, so it still holds every trial; the ', ...
                    'eligible urevent ids are kept in EEG.etc.pipecompare.eligibleUrevents for later epoching.']);
                return;
            end
            elig = EEG.etc.pipecompare.eligibleUrevents;
            lock = pipecompare.eval.Measure.lockingEvents(EEG, contract.allEvents());
            drop = false(1, EEG.trials);
            for k = find(lock > 0)
                drop(k) = ~ismember(double(EEG.event(lock(k)).urevent), elig);
            end
            if ~any(drop), return; end
            assert(~all(drop), 'PipeCompare:TrialRule', 'The trial rule keeps none of the output epochs.');
            [EEG, com] = pop_select(EEG, 'notrial', find(drop));
            com = sprintf('%s %% PipeCompare: keep only trials of the trial rule (%d removed)', com, sum(drop));
        end

        function idx = pickCandidate(result, idx)
            % idx, or the recommended candidate when idx is empty. With
            % several strata there is no single recommendation, and with no
            % feasible candidate there is none at all: both are errors that
            % name what to do (adopt, script and writeScript share this).
            if ~isempty(idx), return; end
            idx = result.ranking.recommended;
            if isempty(idx) && numel(result.ranking.byStratum) > 1
                error('PipeCompare:Adopt', ['The candidates fall into %d strata (different references) that are not ', ...
                    'comparable; choose one: pass the id of a stratum''s recommendation %s.'], ...
                    numel(result.ranking.byStratum), mat2str([result.ranking.byStratum.recommended]));
            end
            assert(~isempty(idx), 'PipeCompare:Adopt', ['No recommended candidate: none satisfies the constraints. ', ...
                'Pass a candidate id explicitly.']);
        end

        function id = identity(result)
            % The checkpoint identity of a search (see searchIdentity).
            id = searchIdentity(result);
        end

        function adopt(result, idx, force)
            % Store candidate idx as a NEW EEGLAB dataset with its full history.
            if nargin < 2, idx = []; end
            if nargin < 3, force = false; end
            idx = pipecompare.run.Executor.pickCandidate(result, idx);
            st = result.ranking.table.status{idx};
            if ~strcmp(st, 'feasible')
                why = result.ranking.table.reason{idx};
                assert(force, 'PipeCompare:Adopt', ['Candidate %d is %s, not feasible: %s. Call adopt(result, %d, true) ', ...
                    'to adopt it anyway.'], idx, st, why, idx);
                pipecompare.utils.log('WARNING: adopting candidate %d although it is %s: %s', idx, st, why);
            end
            [cur, ~] = pipecompare.live.Session.current();
            here = strcmp(pipecompare.live.Session.fingerprint(cur), result.rootFingerprint);
            if ~here && ~force
                inAll = false;
                if evalin('base', 'exist(''ALLEEG'',''var'')')
                    A = evalin('base', 'ALLEEG');
                    for k = 1:numel(A)
                        if strcmp(pipecompare.live.Session.fingerprint(A(k)), result.rootFingerprint), inAll = true; break; end
                    end
                end
                assert(inAll, 'PipeCompare:StaleState', ['The dataset this search started from is no longer in EEGLAB ', ...
                    '(it was changed or removed). Re-run the search, or call adopt(result, id, true) to build the ', ...
                    'candidate from the stored starting copy anyway.']);
                pipecompare.utils.log('The current dataset is not the search''s starting dataset; the candidate is built from that starting dataset.');
            end
            pipecompare.utils.log('Rebuilding candidate %d: %s', idx, result.labels{idx});
            EEG = pipecompare.run.Executor.replay(result, idx, true);
            if strcmp(result.cands(idx).status, 'ok')
                m = pipecompare.eval.Measure.candidate(EEG, result.contract, result.ref, result.options);
                a = [m.objectives.agg]; b = [result.cands(idx).m.objectives.agg];
                same = all(abs(a - b) <= 1e-6 * max(1, abs(b))) && isequal(m.kept, result.cands(idx).m.kept);
                pipecompare.utils.log('Replay check: objectives %s, trials %s (search: %s, %s) -> %s', mat2str(a, 5), ...
                    mat2str(m.kept), mat2str(b, 5), mat2str(result.cands(idx).m.kept), ...
                    pipecompare.utils.ternary(same, 'identical', 'DIFFERENT (non-deterministic step?)'));
                % the evaluation (and the constraints it passed) belongs to the
                % searched data; a rebuilt dataset that differs is not that candidate
                assert(same || force, 'PipeCompare:ReplayMismatch', ['The rebuilt candidate %d differs from the one ', ...
                    'evaluated in the search, so its evaluation does not apply to it. Nothing was stored. Call ', ...
                    'adopt(result, %d, true) to store it anyway.'], idx, idx);
                if ~same, pipecompare.utils.log('WARNING: storing a rebuilt candidate that differs from the evaluated one.'); end
            end
            EEG = storeWeightedAverages(EEG, result.contract);
            EEG.setname = sprintf('%s PipeCompare#%d', result.state.setname, idx);
            assignin('base', 'PIPECOMPARE_ADOPT__', EEG);
            evalin('base', ['[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, PIPECOMPARE_ADOPT__, CURRENTSET, ', ...
                '''gui'', ''off''); clear PIPECOMPARE_ADOPT__;']);
            % (without EEGLAB's main window, eeglab redraw would open it)
            if ~isempty(findobj(groot, 'Tag', 'EEGLAB')), evalin('base', 'eeglab redraw;'); end
            pipecompare.utils.log('Candidate %d stored as a new EEGLAB dataset; its EEG.history lists every step.', idx);
        end

        function result = writeManifest(result)
            % Bind the checkpoint folder to this search. A folder holding
            % files of another search (other data, contract, plan or
            % options) is refused rather than mixed.
            d = result.options.checkpoint;
            if ~isfolder(d), mkdir(d); end
            w = what(d); d = w(1).path;        % absolute (from MATLAB's current folder): found from any folder later
            result.options.checkpoint = d;
            identity = searchIdentity(result);
            result.identity = identity;
            used = ~isempty(dir(fullfile(d, '*.mat')));
            if used
                prev = '';
                if isfile(fullfile(d, 'manifest.mat'))
                    M = load(fullfile(d, 'manifest.mat'), 'identity');
                    if isfield(M, 'identity'), prev = M.identity; end
                end
                assert(strcmp(prev, identity), 'PipeCompare:Checkpoint', ['Checkpoint folder %s already holds files ', ...
                    'of a different search (other data, contract, plan or options). Use an empty folder, or ', ...
                    'pipecompare.PipeCompare.resume(folder) to continue that search.'], d);
                % same search started again: earlier candidates are discarded, not mixed
                delete(fullfile(d, 'leaf_*.mat'));
            end
            root = result.root; %#ok<NASGU>
            save(fullfile(d, 'root.mat'), 'root', 'identity', '-v7.3');
            result.root = [];          % read back on demand (Executor.rootOf)
            M = struct('result', result, 'identity', identity);   % without the caller's progress window
            if isfield(M.result.options, 'progress'), M.result.options.progress = []; end
            save(fullfile(d, 'manifest.mat'), '-struct', 'M', '-v7.3');
            pipecompare.utils.log('Checkpoint folder: %s (resume with pipecompare.PipeCompare.resume).', d);
        end

        function printReport(rep)
            pipecompare.utils.log('Exhaustive search: %d legal pipeline(s) from %d order(s).', rep.nLeaves, rep.nOrders);
            pipecompare.utils.log('%d step executions with prefix sharing (%d without).', rep.nNodes, rep.nStepsUnshared);
            for k = 1:numel(rep.rejected)
                pipecompare.utils.log('  excluded combinations: %s (x%d)', rep.rejected(k).reason, rep.rejected(k).count);
            end
        end
    end
end

% ---------------------------------------------------------------------
function EEG = storeWeightedAverages(EEG, contract)
% Epochs weighted by weight_trials: their weighted average per condition in
% EEG.etc.pipecompare.weightedERP (Measure.weightedAverages), since EEGLAB's
% and ERPLAB's own averages would count every epoch equally.
if isempty(pipecompare.eval.Measure.storedWeights(EEG)) || ~strcmp(contract.analysis, 'erp') || ...
        isempty(contract.conditions) || EEG.trials < 2
    return;
end
EEG.etc.pipecompare.weightedERP = pipecompare.eval.Measure.weightedAverages(EEG, contract);
pipecompare.utils.log(['The epochs are weighted by their noise: the weighted average of each condition is in ', ...
    'EEG.etc.pipecompare.weightedERP.data (channels x times x conditions). EEGLAB''s and ERPLAB''s own ', ...
    'averages count every epoch equally and do not give the ERP this pipeline was scored on.']);
end

function c = emptyCand()
% One evaluated pipeline. Every field is set here and only here; the
% producers and consumers are:
%   id, key, stratum        leaf index, pipeline key and measure-defining
%                           choices (Plan.enumerate leaves)
%   status, message         'ok' | 'rejected' (e.g. all epochs removed) |
%                           'excluded' (not run: its filters alone distort
%                           the known signal) | 'failed', with the reason
%   m                       Measure.candidate: kept, retention, objectives
%                           (per-measure SME), composite
%   signal                  Injection.compare: source, amplitudeError,
%                           latencyShiftMs, artifactPct, waveformCorr,
%                           topoCorr, notApplicable
%   interpolatedFraction    share of the montage interpolated
%   icsRemoved, rejectedEpochs   counts from the steps' decisions
%   coms                    EEGLAB commands that rebuild it (script/adopt)
%   seconds                 run time of its steps
%   unmatched               step types re-run (not decision-matched) on
%                           the signal copy
%   notes                   facts the user must know that are not
%                           failures (e.g. epochs marked but not removed)
%   ica                     the ICA decompositions of its path (weights,
%                           sphere, channels, commands), reused by replay
%   overLimit               per channel (labels, counts), the rejected
%                           epochs in which it was over the limit: names
%                           bad channels the detection missed
%   steps                   per step of its path, in order: the step's
%                           type and settings and what it did on these
%                           data (channels interpolated, components and
%                           epochs removed), for the list of steps
% Rank.run reads status, m, signal, interpolatedFraction and stratum.
c = struct('id', 0, 'key', '', 'stratum', '', 'status', 'pending', 'message', '', 'm', [], 'signal', [], ...
    'interpolatedFraction', NaN, 'icsRemoved', 0, 'rejectedEpochs', 0, 'coms', {{}}, 'seconds', 0, 'unmatched', {{}}, ...
    'notes', {{}}, 'ica', {{}}, 'overLimit', noOverLimit(), 'steps', {{}});
end

function o = noOverLimit()
% per channel, the rejected epochs in which it was over the limit
o = struct('labels', {{}}, 'counts', []);
end

function out = failUnder(env, under, in, ME)
% A step failed: every not yet evaluated candidate below it fails with the
% reason (the search goes on with the other branches).
pipecompare.utils.log('FAILED %s: %s (%d candidate(s) affected)', in.label, ME.message, numel(under));
out = repmat(emptyCand(), 0, 1);
for li = under
    if env.done(li), continue; end
    c = emptyCand(); c.id = li; c.key = env.leaves(li).key; c.stratum = env.leaves(li).stratum;
    c.status = failStatus(ME); c.message = sprintf('%s: %s', in.label, ME.message);
    out(end+1, 1) = c; %#ok<AGROW>
    saveLeaf(env, c);
end
end

function s = failStatus(ME)
% A step that removes every epoch is an outcome (retention 0), not an error.
if strcmp(ME.identifier, 'PipeCompare:AllRejected'), s = 'rejected'; else, s = 'failed'; end
end

function [ctx, acc] = advance(in, ctx, acc, info, coms, unmatched, secs)
% State carried to the next step. One definition for the serial search,
% the parallel trunk and replay, so all three run identical steps.
if strcmp(in.type, 'highpass'), ctx.highpass = max(ctx.highpass, in.params.cutoff); end
if strcmp(in.type, 'native')
    for stmt = pipecompare.run.Native.statements(in.params.command)
        e = pipecompare.live.History.classify(stmt{1});
        if isfield(e.params, 'locutoff') && isfinite(e.params.locutoff) && ~(isfield(e.params, 'revfilt') && e.params.revfilt)
            ctx.highpass = max(ctx.highpass, e.params.locutoff);
        end
    end
end
if isempty(acc), return; end
acc.seconds = acc.seconds + secs;
acc.coms = [acc.coms coms];
acc.unmatched = unique([acc.unmatched unmatched]);
if isfield(info, 'interpolated'), acc.interpolated = union(acc.interpolated, lower(info.interpolated)); end
if isfield(info, 'icsRemoved'), acc.icsRemoved = acc.icsRemoved + info.icsRemoved; end
if isfield(info, 'rejected'), acc.rejected = acc.rejected + info.rejected; end
if isfield(info, 'ica'), acc.ica{end+1} = info.ica; end
if isfield(info, 'overLimit') && ~isempty(info.overLimit.counts)
    [labels, ~, j] = unique([acc.overLimit.labels(:); info.overLimit.labels(:)]);
    acc.overLimit = struct('labels', {labels(:)'}, 'counts', accumarray(j(:), [acc.overLimit.counts(:); info.overLimit.counts(:)])');
end
acc.steps{end+1} = stepFacts(in, info, secs);
end

function f = stepFacts(in, info, secs)
% Round 8 timing: per-step wall seconds. Rank, replayDecision and
% searchIdentity do not read this field, so ranking, replay and
% checkpoint identity are unchanged by adding it.
if nargin < 3, secs = NaN; end
% What a step did, for the results' list of steps: its settings and the
% decisions it made on these data (channels, components, epochs, epochs
% kept by interpolating channels in them), and for ICLabel how many
% components it took for brain activity and for 'Other' and how many data
% points ICA had (the ICA check); for a rejection limit chosen from the
% data, the limit (uvChosen); for an ICLabel threshold chosen from the
% data, the threshold (thresholdChosen) and how many components it
% considered (icsConsidered); the components and epochs removed (comps,
% rejIdx), which fine-tuning the cutoffs holds fixed; for a filter,
% whether its cutoff was fine-tuned from the data (cutoffRefined); for
% weighted epochs, how much more the cleanest counts than the noisiest
% (weightRatio) and the weights (trialWeights), which fine-tuning the
% cutoffs also holds fixed; for a choice made from the data (rejection
% limit, ICLabel threshold, components by SME), the scores of every option
% it chose among (curves), with which the ranking's selection check makes
% the choice again on half of the trials (Rank.crossfit); for a rejection,
% how far the mean score of the epochs kept is from that of all of them
% (keptShift, Presets.keptShiftText); for bad channels chosen by the total
% error, those left as they are (leftAsIs), and those bad for part of the
% recording (stretches: channel, samples); for the eye channels regressed
% out, the channels and coefficients (eog).
f = struct('type', in.type, 'params', in.params, 'seconds', secs, 'interpolated', {{}}, 'removed', {{}}, ...
    'icsRemoved', NaN, 'icsTotal', NaN, 'rejected', NaN, 'epochsBefore', NaN, 'epochsInterpolated', NaN, ...
    'icsBrain', NaN, 'icsOther', NaN, 'otherMedian', NaN, 'icaPoints', NaN, 'uvChosen', NaN, 'thresholdChosen', NaN, 'icsConsidered', NaN, 'icsProtected', NaN, ...
    'comps', [], 'rejIdx', [], 'cutoffRefined', false, 'weightRatio', NaN, 'trialWeights', [], 'curves', [], 'keptShift', [], 'leftAsIs', {{}}, 'eog', [], 'stretches', [], 'stretchCheck', []);
for n = {'interpolated', 'removed', 'leftAsIs'}
    if isfield(info, n{1}), f.(n{1}) = cellstr(info.(n{1})); end
end
for n = {'icsRemoved', 'icsTotal', 'rejected', 'epochsBefore', 'epochsInterpolated', 'icsBrain', 'icsOther', 'otherMedian', 'icaPoints', 'uvChosen', 'thresholdChosen', 'icsConsidered', 'icsProtected', ...
        'comps', 'rejIdx', 'cutoffRefined', 'weightRatio', 'trialWeights', 'curves', 'keptShift', 'eog', 'stretches', 'stretchCheck'}
    if isfield(info, n{1}), f.(n{1}) = info.(n{1}); end
end
end

function [S, unmatched] = replayCopy(in, S, info, ctx)
% The step's decisions applied to the signal-check copy. An error here
% costs only the signal check of the candidates below (they are excluded
% with 'signal check missing'), not the candidates: their real data ran.
unmatched = {};
if isempty(S), return; end
try
    [~, S, matched] = evalc('pipecompare.run.Steps.replayDecision(in, S, info, ctx)');
    if ~matched, unmatched = {in.type}; end
catch ME
                    if startsWith(ME.identifier,'PipeCompare:Checkpoint') || strcmp(ME.identifier,'PipeCompare:Interrupted'),rethrow(ME);end
    pipecompare.utils.log('Signal check: the copy failed at step %s (%s); the pipelines below have no signal check.', ...
        in.label, ME.message);
    S = [];
end
end

function cx = choiceContext(env, node, ctx, S)
% An ICLabel threshold chosen from the data (pipecompare.eval.IcThreshold)
% also needs the steps that follow it and the signal-check copy; the
% steps that follow must be the same in every pipeline below it.
cx = ctx;
in = env.tree(node).inst;
d = env.tree(node).depth;
under = leavesUnder(env.tree, node);
if isRepairChoice(in)
    % bad channels interpolated by their effect on the total error
    % (ChannelChoice): the linear steps up to the rejection, when they are
    % the same in every pipeline below (otherwise all are interpolated)
    after = linearAfter(env.leaves(under(1)).path(d+1:end));
    if all(arrayfun(@(li) isequal(plainSteps(linearAfter(env.leaves(li).path(d+1:end))), plainSteps(after)), under))
        cx.after = after; cx.signal = struct('S', S, 'truth', env.truth);
    else
        pipecompare.utils.log(['Bad channels: the pipelines below differ in the reference, epochs or rejection, so the ', ...
            'channels to interpolate are not chosen by the total error.']);
    end
    return;
end
if ~(strcmp(in.type, 'icremove') && ischar(in.params.threshold)), return; end
% the steps it is chosen for: up to the epoch rejection that follows (the
% steps after that do not change the choice)
upto = @(path) path(1:min([numel(path) find(cellfun(@(x) strcmp(x.type, 'reject_threshold'), path), 1)]));
after = upto(env.leaves(under(1)).path(d+1:end));
plain = @(path) cellfun(@(x) {x.type, x.params}, path, 'UniformOutput', false);
same = all(arrayfun(@(li) isequal(plain(upto(env.leaves(li).path(d+1:end))), plain(after)), under));
assert(same, 'PipeCompare:IcThreshold', ['An ICLabel threshold chosen from the data (''auto'') needs the same steps ', ...
    'after it in every pipeline that shares it (it is chosen for the rejection that follows); compare the steps ', ...
    'after it in separate runs, or set the thresholds (e.g. 0.7 | 0.8 | 0.9).']);
cx.after = after;
cx.signal = struct('S', S, 'truth', env.truth);
end

function tf = isRepairChoice(in)
% bad channels interpolated by their effect on the total error
tf = strcmp(in.type, 'badchannels') && strcmp(pipecompare.utils.fieldOr(in.params, 'repair', 'all'), 'error');
end

function after = linearAfter(path)
% the steps ChannelChoice applies: re-referencing, epoching and baseline
% removal up to the epoch rejection, and that rejection
k = find(cellfun(@(x) strcmp(x.type, 'reject_threshold'), path), 1);
if ~isempty(k), path = path(1:k); end
after = path(cellfun(@(x) any(strcmp(x.type, {'reref', 'epoch', 'baseline', 'reject_threshold'})), path));
end

function p = plainSteps(path)
p = cellfun(@(x) {x.type, x.params}, path, 'UniformOutput', false);
end

function cx = pathContext(path, q, cx, S, truth)
% What a choice made from the data in step q of one pipeline needs: the
% steps after it and the signal-check copy (as choiceContext for the tree)
in = path{q};
if strcmp(in.type, 'icremove') && ischar(in.params.threshold)
    after = path(q+1:end);
    after = after(1:min([numel(after) find(cellfun(@(x) strcmp(x.type, 'reject_threshold'), after), 1)]));
    cx.after = after; cx.signal = struct('S', S, 'truth', truth);
elseif isRepairChoice(in)
    cx.after = linearAfter(path(q+1:end)); cx.signal = struct('S', S, 'truth', truth);
end
end

function [E2, coms, info] = runStep(in, E, ctx, verbose)
if strcmp(verbose, 'full')
    [E2, coms, info] = pipecompare.run.Steps.run(in, E, ctx);
    return;
end
[txt, E2, coms, info] = evalc('pipecompare.run.Steps.run(in, E, ctx)');
lines = regexp(txt, '[^\n]+', 'match');
keep = lines(~cellfun(@isempty, regexpi(lines, 'warning|error|\[PipeCompare\]', 'once')));
for k = 1:numel(keep), fprintf('    %s\n', strtrim(keep{k})); end
end

function L = levelsOf(leaves)
% per pipeline, its settings: one {slot, key} row per step it runs (for
% the ranking's factorial model, Rank.factorialModel)
L = cell(numel(leaves), 1);
for i = 1:numel(leaves)
    p = leaves(i).path;
    L{i} = [cellfun(@(in) in.slot, p(:), 'UniformOutput', false), cellfun(@(in) in.key, p(:), 'UniformOutput', false)];
end
end

function c = evaluateLeaf(env, li, E, S, acc)
c = emptyCand();
c.id = li; c.key = env.leaves(li).key; c.stratum = env.leaves(li).stratum;
c.coms = acc.coms; c.seconds = acc.seconds;
c.icsRemoved = acc.icsRemoved; c.rejectedEpochs = acc.rejected;
c.interpolatedFraction = numel(acc.interpolated) / env.nbchan;
c.unmatched = acc.unmatched; c.ica = acc.ica; c.overLimit = acc.overLimit; c.steps = acc.steps;
try
    c.m = pipecompare.eval.Measure.candidate(E, env.contract, env.ref, env.opts);
    nm = markedNotRemoved(E);
    if nm > 0
        c.notes{end+1} = sprintf(['%d epoch(s) are marked for rejection (EEG.reject) but still in the data; ', ...
            'marks do not remove epochs, so they count in the scores. Remove them in the plan ', ...
            '(e.g. pop_rejepoch) if they should not count.'], nm);
    end
    [~, ~, tcom] = evalc('pipecompare.run.Executor.selectEligible(E, env.contract)');
    if ~isempty(tcom), c.coms{end+1} = tcom; end
    if ~isempty(S)
        c.signal = pipecompare.eval.Injection.compare(S, env.contract, env.truth, env.leaves(li).path);
        if isfield(env, 'spread') && numel(env.spread) >= li, c.signal.shapeSpread = env.spread{li}; end
    end
    c.status = 'ok';
catch ME
                    if startsWith(ME.identifier,'PipeCompare:Checkpoint') || strcmp(ME.identifier,'PipeCompare:Interrupted'),rethrow(ME);end
    c.status = 'failed'; c.message = ['evaluation: ' ME.message];
end
el = toc(env.tStart);
if strcmp(c.status, 'ok')
    pipecompare.utils.log('candidate %d: objectives %s, retention %s, %.0f s elapsed', li, ...
        mat2str([c.m.objectives.agg], 4), mat2str(round(100 * c.m.retention)), el);
else
    pipecompare.utils.log('candidate %d failed: %s', li, c.message);
end
end

function [result, allc, added] = refineCutoffs(result, env, root, S, ctx, acc, allc, R)
% The filter cutoffs of the recommended pipeline fine-tuned from the data
% (pipecompare.eval.Cutoff): with its decisions held fixed, the objective
% is computed on a finer grid of cutoffs between those compared. When
% other cutoffs are better, the pipeline with them runs like every other
% one (its steps take their own decisions on these data) and joins the
% ranking. Its steps before the first filter run again from the starting
% data (the search does not keep intermediate data), its ICA reused.
added = false;
ib = R.recommended;
if isempty(ib) || isempty(S) || root.trials > 1, return; end
c = allc(ib); path = env.leaves(ib).path;
types = cellfun(@(i) i.type, path, 'UniformOutput', false);
fi = find(ismember(types, {'highpass', 'lowpass'}));
if isempty(fi) || numel(fi) > 2 || ~all(cellfun(@(i) isnumeric(i.params.cutoff), path(fi))), return; end
% the cutoffs compared in each of its filters' slots
cuts = cell(1, numel(fi));
for k = 1:numel(fi)
    in = path{fi(k)}; v = [];
    for L = env.leaves(:)'
        for j = find(cellfun(@(x) strcmp(x.slot, in.slot) && strcmp(x.type, in.type), L.path))
            if isnumeric(L.path{j}.params.cutoff), v(end+1) = L.path{j}.params.cutoff; end %#ok<AGROW>
        end
    end
    if numel(unique(v)) > 1, cuts{k} = pipecompare.eval.Cutoff.grid(v); end
end
if all(cellfun(@isempty, cuts)), return; end
first = fi(1); rest = path(first:end); rt = types(first:end);
why = '';
if ~strcmp(env.contract.analysis, 'erp') || ~pipecompare.eval.Threshold.applies(env.contract)
    why = 'it needs mean amplitudes';
elseif ~all(ismember(rt, {'highpass', 'lowpass', 'icremove', 'epoch', 'baseline', 'reject_threshold', 'weight_trials'})) || ...
        sum(strcmp(rt, 'epoch')) ~= 1
    why = 'only IC removal, epoching, baseline removal, an amplitude rejection and weighting the epochs may follow the filters';
elseif any(cellfun(@(i) strcmp(i.type, 'reject_threshold') && pipecompare.utils.fieldOr(i.params, 'interpolate', 0) > 0, rest))
    why = 'epochs are repaired (the rejection then also interpolates channels)';
elseif any(ismember(types(1:first-1), {'icremove', 'resample'}))
    why = 'components are removed or the data resampled before the filters';
elseif numel(c.steps) ~= numel(path)
    why = 'the steps of the pipeline were not recorded';
end
if ~isempty(why)
    pipecompare.utils.log('Filter cutoffs not fine-tuned: %s.', why);
    return;
end
dec = struct('comps', [], 'rejIdx', []);
k = find(strcmp(types, 'icremove'), 1); if ~isempty(k), dec.comps = c.steps{k}.comps; end
k = find(strcmp(types, 'reject_threshold'), 1); if ~isempty(k), dec.rejIdx = c.steps{k}.rejIdx; end
k = find(strcmp(types, 'weight_trials'), 1); if ~isempty(k), dec.weights = c.steps{k}.trialWeights; end
pipecompare.utils.log(['Fine-tuning the filter cutoffs of pipeline %d: its decisions held fixed, only the measured ', ...
    'electrodes are filtered for each pair of cutoffs on a finer grid.'], ib);
S0 = S; E = root; nIca = 0;
try
    for k = 1:first-1
        in = path{k};
        if tick(env, 0, in), return; end
        if strcmp(in.type, 'ica') && nIca < numel(c.ica)
            nIca = nIca + 1; d = c.ica{nIca};
            E.icaweights = d.icaweights; E.icasphere = d.icasphere;
            E.icachansind = d.icachansind; E.icawinv = []; E.icaact = [];
            E = eeg_checkset(E);
            coms = d.coms; info = struct('ica', d);
        else
            cx = ctx;
            if isRepairChoice(in), cx.chosenChannels = c.steps{k}.interpolated; end   % the channels the search interpolated
            [E, coms, info] = runStep(in, E, cx, env.opts.verbose);
        end
        for q = 1:numel(coms), E = eeg_hist(E, coms{q}); end
        [S, unmatched] = replayCopy(in, S, info, ctx);
        [ctx, acc] = advance(in, ctx, acc, info, coms, unmatched, 0);
    end
    cenv = struct('contract', env.contract, 'ref', env.ref, 'rank', env.opts, 'S', S0, 'truth', env.truth, 'grid', {cuts});
    [cut, cinfo] = pipecompare.eval.Cutoff.search(E, rest, dec, struct('m', c.m, 'gain', c.signal.gain, ...
        'spread', pipecompare.utils.fieldOr(c.signal, 'shapeSpread', [])), cenv);
catch ME
                    if startsWith(ME.identifier,'PipeCompare:Checkpoint') || strcmp(ME.identifier,'PipeCompare:Interrupted'),rethrow(ME);end
    pipecompare.utils.log('Filter cutoffs not fine-tuned: %s', ME.message);
    return;
end
own = cellfun(@(i) i.params.cutoff, path(fi));
if isempty(cut)
    pipecompare.utils.log('Filter cutoffs: %s (%d pairs of cutoffs computed).', cinfo.note, cinfo.pairs);
    return;
end
path2 = path;
for k = 1:numel(fi), path2{fi(k)} = pipecompare.plan.Plan.reparam(path{fi(k)}, 'cutoff', cut(k)); end
key2 = strjoin(cellfun(@(i) i.key, path2, 'UniformOutput', false), ' > ');
pipecompare.utils.log(['Filter cutoffs: %s Hz instead of %s Hz would lower the objective from %.4g to %.4g with the ', ...
    'same decisions (%d pairs computed); that pipeline now runs with its own decisions.'], mat2str(cut), mat2str(own), ...
    cinfo.objectiveBase, cinfo.objectiveBest, cinfo.pairs);
if any(strcmp(result.labels, key2))
    pipecompare.utils.log('That pipeline was already compared.');
    return;
end
li = numel(env.leaves) + 1;
leaf = struct('path', {path2}, 'key', key2, 'order', env.leaves(ib).order, 'stratum', env.leaves(ib).stratum);
env.leaves(li) = leaf;
try
    sp = pipecompare.eval.Injection.shapeSpreads(S0, env.contract, env.truth, leaf); env.spread{li} = sp{1};
catch
    env.spread{li} = [];   % the shape part of its bias left out
end
try
    for k = first:numel(path2)
        in = path2{k};
        if tick(env, 0, in), return; end
        pipecompare.utils.log('step %s', in.label);
        cx = ctx;
        if strcmp(in.type, 'icremove') && ischar(in.params.threshold)
            cx.after = path2(k+1:end); cx.signal = struct('S', S, 'truth', env.truth);
        end
        t0 = tic;
        [E, coms, info] = runStep(in, E, cx, env.opts.verbose);
        if any(k == fi) && in.params.cutoff ~= path{k}.params.cutoff, info.cutoffRefined = true; end
        for q = 1:numel(coms), fprintf('    %s\n', coms{q}); E = eeg_hist(E, coms{q}); end
        [S, unmatched] = replayCopy(in, S, info, ctx);
        [ctx, acc] = advance(in, ctx, acc, info, coms, unmatched, toc(t0));
    end
catch ME
                    if startsWith(ME.identifier,'PipeCompare:Checkpoint') || strcmp(ME.identifier,'PipeCompare:Interrupted'),rethrow(ME);end
    pipecompare.utils.log('The pipeline with the fine-tuned cutoffs failed at step %s: %s', in.label, ME.message);
    return;
end
cand = evaluateLeaf(env, li, E, S, acc);
cand.notes{end+1} = sprintf(['filter cutoffs fine-tuned from these data, starting from pipeline %d (%s Hz; ', ...
    'objective %.4g with its decisions, %.4g with the fine-tuned cutoffs)'], ib, mat2str(own), ...
    cinfo.objectiveBase, cinfo.objectiveBest);
result.leaves(li) = leaf; result.labels{li} = key2;
allc(li, 1) = cand;
saveLeaf(env, cand);   % an adaptive candidate is ranking input too: without
                       % this file the recommendation cannot be replayed
added = true;
end

function [result, allc, added] = eogCheck(result, env, root, S, ctx, acc, allc, R)
% The recommended pipeline run again with the eye channels regressed out
% (Steps.eogRegress) just before its first step that removes components,
% makes epochs or rejects or weights them; its steps before that make the
% same decisions as in the search (ICA reused, the channels the search
% interpolated), the steps from there take their own. It joins the
% ranking, so it is recommended when it lowers the total error (noise and
% what blinks add).
added = false;
ib = R.recommended;
if isempty(ib), return; end
if isempty(pipecompare.eval.Injection.dedicatedEyes(root.chanlocs))
    pipecompare.utils.log('Eye channels: none of their own (EOG), so regressing them out is not tried.');
    return;
end
c = allc(ib); path = env.leaves(ib).path;
types = cellfun(@(i) i.type, path, 'UniformOutput', false);
if any(strcmp(types, 'eogregress')) || numel(c.steps) ~= numel(path), return; end
q = find(ismember(types, {'icremove', 'epoch', 'reject_threshold', 'reject_jointprob', 'reject_kurtosis', 'weight_trials'}), 1);
if isempty(q), q = numel(path) + 1; end
d = pipecompare.plan.Catalog.get('eogregress');
in0 = struct('slot', 'eog', 'type', 'eogregress', 'params', struct('channels', {{}}), 'key', 'eogregress', ...
    'label', d.label, 'searched', {{}}, 'defining', '');
path2 = [path(1:q-1) {in0} path(q:end)];
key2 = strjoin(cellfun(@(i) i.key, path2, 'UniformOutput', false), ' > ');
if any(strcmp(result.labels, key2)), return; end
pipecompare.utils.log('Eye channels: pipeline %d runs again with them regressed out, to see whether that lowers the total error.', ib);
li = numel(env.leaves) + 1;
leaf = struct('path', {path2}, 'key', key2, 'order', env.leaves(ib).order, 'stratum', env.leaves(ib).stratum);
env.leaves(li) = leaf;
if isfield(env, 'spread') && numel(env.spread) >= ib, env.spread{li} = env.spread{ib}; end   % the same filters
E = root; nIca = 0;
try
    for k = 1:numel(path2)
        in = path2{k};
        if tick(env, 0, in), return; end
        t0 = tic;
        if k < q && strcmp(in.type, 'ica') && nIca < numel(c.ica)
            nIca = nIca + 1; dd = c.ica{nIca};
            E.icaweights = dd.icaweights; E.icasphere = dd.icasphere;
            E.icachansind = dd.icachansind; E.icawinv = []; E.icaact = [];
            E = eeg_checkset(E);
            coms = dd.coms; info = struct('ica', dd);
        else
            cx = ctx;
            if k < q && isRepairChoice(in), cx.chosenChannels = c.steps{k}.interpolated; end   % as in the search
            if k >= q, cx = pathContext(path2, k, cx, S, env.truth); end
            pipecompare.utils.log('step %s', in.label);
            [E, coms, info] = runStep(in, E, cx, env.opts.verbose);
        end
        if k < q && pipecompare.utils.fieldOr(c.steps{k}, 'cutoffRefined', false), info.cutoffRefined = true; end
        for j = 1:numel(coms), E = eeg_hist(E, coms{j}); end
        [S, unmatched] = replayCopy(in, S, info, ctx);
        [ctx, acc] = advance(in, ctx, acc, info, coms, unmatched, toc(t0));
    end
catch ME
                    if startsWith(ME.identifier,'PipeCompare:Checkpoint') || strcmp(ME.identifier,'PipeCompare:Interrupted'),rethrow(ME);end
    pipecompare.utils.log('Eye channels: pipeline %d with them regressed out failed (%s).', ib, ME.message);
    return;
end
cand = evaluateLeaf(env, li, E, S, acc);
cand.notes{end+1} = sprintf('pipeline %d with the eye channels regressed out (tried for the recommended pipeline)', ib);
if strcmp(cand.status, 'ok') && strcmp(c.status, 'ok')
    pipecompare.utils.log('Eye channels: pipeline %d with them regressed out is pipeline %d.', ib, li);
end
result.leaves(li) = leaf; result.labels{li} = key2;
allc(li, 1) = cand;
saveLeaf(env, cand);   % an adaptive candidate is ranking input too: without
                       % this file the recommendation cannot be replayed
added = true;
end

function out = icaSeedCheck(env, root, S, ctx, acc, allc, R, ref, opts)
% How much the recommendation depends on where ICA starts. runica starts
% from a fixed random state, so another start is another order of the
% samples it sees (Steps.ica, ctx.icaSeed). With each seed 1..opts.icaSeeds
% every pipeline with an ICA step runs again, every step taking its own
% decisions, all of them under the same seed: the tree is walked as in the
% search, so pipelines that share the steps up to ICA share one ICA per
% seed (Standard: one ICA per seed for all its pipelines). The ranking is
% then made again (no cross-fitting) and the recommendation compared. A
% pipeline added by fine-tuning the cutoffs runs again too, reusing the
% ICA of a pipeline with the same steps before it. out: [] when not
% checked; else seeds, ids (re-run), objective (ids x 1+seeds, the first
% column the search's), recommended (1+seeds), same (true when every seed
% recommends the same pipeline).
out = [];
if ~isscalar(R.byStratum), return; end
k = R.byStratum.recommended; T = R.table;
hasIca = @(li) any(cellfun(@(i) strcmp(i.type, 'ica'), env.leaves(li).path));
nTree = numel(env.done);
run = arrayfun(hasIca, 1:numel(env.leaves)) & ~strcmp({allc.status}, 'excluded');
ids = find(run);
if isempty(ids), return; end
n = opts.icaSeeds;
pipecompare.utils.log(['ICA check: the %d pipeline(s) with ICA run again with ICA started in %d other way(s) ', ...
    '(one ICA per seed for pipelines sharing the steps before it), to see whether the recommendation depends on it.'], ...
    numel(ids), n);
e2 = env; e2.checkpoint = ''; e2.done = ~run(1:nTree);
if isfield(env.opts, 'progress') && ~isempty(env.opts.progress)
    f0 = env.opts.progress;   % steps shown, candidates not counted again
    e2.opts.progress = @(varargin) callProgress(f0, [{0} varargin(2:end)]);
end
obj = [T.objective(ids) nan(numel(ids), n)]; rec = [k nan(1, n)];
ro = opts; ro.nSplits = 0;
for seed = 1:n
    c2 = allc; cx = ctx; cx.icaSeed = seed;
    cands = pipecompare.run.Executor.runSubtree(e2, 1, root, S, cx, acc);
    if tick(e2, 0), return; end   % stopped
    for j = 1:numel(cands), c2(cands(j).id) = cands(j); end
    for li = ids(ids > nTree)
        c = rerunWithIca(e2, li, root, S, cx, acc, c2(1:nTree), allc(li));
        if isempty(c), out = []; return; end
        c2(li) = c;
    end
    R2 = pipecompare.eval.Rank.run(c2, ref, ro, levelsOf(env.leaves));
    obj(:, seed + 1) = R2.table.objective(ids);
    if ~isempty(R2.recommended), rec(seed + 1) = R2.recommended; end
end
out = struct('seeds', 1:n, 'ids', ids(:)', 'objective', obj, 'recommended', rec, 'same', all(rec == k));
ch = max(abs(obj(:, 2:end) ./ obj(:, 1) - 1), [], 2);
if out.same
    pipecompare.utils.log(['ICA check: with ICA started in %d other way(s), the same pipeline (%d) is recommended; ', ...
        'the total error of the pipelines run again changes by at most %.0f%%.'], n, k, 100 * max([0; ch(isfinite(ch))]));
else
    alt = unique(rec(rec ~= k & isfinite(rec)));
    alt = strjoin(arrayfun(@(x) sprintf('pipeline %d', x), alt, 'UniformOutput', false), ', ');
    if isempty(alt), alt = 'no pipeline'; end
    pipecompare.utils.log(['ICA check: with ICA started in another way, %s would be recommended instead of %d: ', ...
        'this recommendation depends on where ICA starts.'], alt, k);
end
end

function c = rerunWithIca(env, li, root, S, ctx, acc, done, orig)
% A pipeline that is not in the tree (fine-tuned cutoffs) run again under
% ctx.icaSeed. Its ICA is taken from a re-run pipeline of the tree with the
% same steps before it (so not fitted again) when there is one; what the
% search noted about it (cutoffs fine-tuned, shape spread) is kept from
% orig. [] when it fails.
c = [];
path = env.leaves(li).path; keys = cellfun(@(i) i.key, path, 'UniformOutput', false);
E = root; nIca = 0;
try
    for q = 1:numel(path)
        in = path{q};
        if tick(env, 0, in), return; end
        d = [];
        if strcmp(in.type, 'ica')
            nIca = nIca + 1;
            for j = 1:numel(done)
                pj = env.leaves(j).path;
                if numel(pj) >= q && numel(done(j).ica) >= nIca && ...
                        isequal(cellfun(@(i) i.key, pj(1:q), 'UniformOutput', false), keys(1:q))
                    d = done(j).ica{nIca}; break;
                end
            end
        end
        if isempty(d)
            cx = pathContext(path, q, ctx, S, env.truth);
            [E, coms, info] = runStep(in, E, cx, env.opts.verbose);
        else
            E.icaweights = d.icaweights; E.icasphere = d.icasphere;
            E.icachansind = d.icachansind; E.icawinv = []; E.icaact = [];
            E = eeg_checkset(E);
            coms = d.coms; info = struct('ica', d);
        end
        [S, unmatched] = replayCopy(in, S, info, ctx);
        [ctx, acc] = advance(in, ctx, acc, info, coms, unmatched, 0);
    end
    c = evaluateLeaf(env, li, E, S, acc);
    for q = 1:min(numel(c.steps), numel(orig.steps))
        if pipecompare.utils.fieldOr(orig.steps{q}, 'cutoffRefined', false), c.steps{q}.cutoffRefined = true; end
    end
    if isstruct(orig.signal) && isfield(orig.signal, 'shapeSpread') && isstruct(c.signal)
        c.signal.shapeSpread = orig.signal.shapeSpread;
    end
catch ME
                    if startsWith(ME.identifier,'PipeCompare:Checkpoint') || strcmp(ME.identifier,'PipeCompare:Interrupted'),rethrow(ME);end
    pipecompare.utils.log('ICA check: pipeline %d failed (%s); not checked.', li, ME.message);
    c = [];
end
end

function stop = tick(env, n, in)
% Tell the caller's opts.progress that n more candidates are evaluated (and,
% before a step, which step starts) and return its request to stop.
stop = false;
if toc(env.tStart)>=pipecompare.utils.fieldOr(env.opts,'maxSeconds',Inf)
    error('PipeCompare:Interrupted','Time slice ended; saved candidates and ICA fits can be resumed.');
end
if ~isfield(env.opts, 'progress') || isempty(env.opts.progress), return; end
if nargin < 3, stop = callProgress(env.opts.progress, {n}); else, stop = callProgress(env.opts.progress, {n, in}); end
end

function stop = callProgress(f, args)
% a callback of one argument is only told the count
if nargin(f) == 1, args = args(1); end
stop = f(args{:});
end

function relayProgress(prog, args, stopFile)
% client side of a parallel run: show the workers' progress; Stop pressed
% there leaves the file the workers look for
if callProgress(prog, args) && ~isfile(stopFile), fclose(fopen(stopFile, 'w')); end
end

function stop = workerProgress(q, stopFile, args)
% worker side of a parallel run (the progress callback tick() calls)
send(q, args);
stop = isfile(stopFile);
end

function deleteIfThere(file)
if isfile(file), delete(file); end
end

function saveLeaf(env, cand)
if isempty(env.checkpoint), return; end
id = cand.id; key = cand.key; identity = env.identity; %#ok<NASGU>
tmp = fullfile(env.checkpoint, sprintf('leaf_%06d.tmp.mat', id));
save(tmp, 'id', 'key', 'cand', 'identity', '-v7.3');
movefile(tmp, fullfile(env.checkpoint, sprintf('leaf_%06d.mat', id)), 'f');
if isfield(env.opts, 'stopAfter') && numel(dir(fullfile(env.checkpoint, 'leaf_*.mat'))) >= env.opts.stopAfter
    % simulated interruption (tests): everything saved so far is kept
    error('PipeCompare:Interrupted', 'Search interrupted after %d candidates (opts.stopAfter).', env.opts.stopAfter);
end
end

function saveRanking(env, R, labels)
% The final ranking, saved with the checkpoint for audit and replay: it
% covers the adaptive candidates too, which the manifest predates (the
% manifest is written before the search runs, so it cannot list them).
% resume() never reads this file — it rebuilds the ranking from the leaf
% files — so saving it cannot disturb a resumed search. Only saved for a
% completed search: a stopped search has no final ranking.
if isempty(env.checkpoint), return; end
identity = env.identity; %#ok<NASGU>
tmp = fullfile(env.checkpoint, 'ranking.tmp.mat');
save(tmp, 'R', 'labels', 'identity', '-v7.3');
movefile(tmp, fullfile(env.checkpoint, 'ranking.mat'), 'f');
end

function n = dataLevels(path, fs0)
% Copies of the data a depth-first branch holds, in units of the starting
% data (the start, then one per step; resampling shrinks what follows).
f = 1; n = 1;
for q = 1:numel(path)
    if strcmp(path{q}.type, 'resample'), f = f * path{q}.params.fs / fs0; fs0 = path{q}.params.fs; end
    n = n + f;
end
end

function n = markedNotRemoved(EEG)
% Epochs that a marking step (EEGLAB's pop_eegthresh/pop_jointprob/...,
% ERPLAB's artifact detection, which marks EEG.reject.rejmanual) flagged
% but no step removed. EEGLAB keeps marks and data apart: only
% pop_rejepoch removes epochs (Delorme & Makeig, 2004); ERPLAB's averager
% honours its own flags (Lopez-Calderon & Luck, 2014), PipeCompare's scores
% do not.
n = 0;
if EEG.trials <= 1 || ~isfield(EEG, 'reject') || ~isstruct(EEG.reject), return; end
f = {'rejmanual','rejthresh','rejconst','rejjp','rejkurt','rejfreq'};
any_ = false(1, EEG.trials);
for k = 1:numel(f)
    if isfield(EEG.reject, f{k}) && numel(EEG.reject.(f{k})) == EEG.trials
        any_ = any_ | logical(EEG.reject.(f{k})(:)');
    end
end
n = sum(any_);
end

function hp = rootHighpass(state)
hp = 0;
if ~isempty(state.filters.highpass), hp = max(state.filters.highpass); end
end

function out = leavesUnder(tree, node)
out = tree(node).leaves;
for c = tree(node).children
    out = [out leavesUnder(tree, c)]; %#ok<AGROW>
end
end

function tf = canParallel()
tf = license('test', 'Distrib_Computing_Toolbox') && ~isempty(ver('parallel'));
if tf
    try
        % 'Processes' is the profile's name from R2022b on, 'local' before
        p = gcp('nocreate'); if isempty(p), try, parpool('Processes'); catch, parpool('local'); end, end
    catch
        tf = false;
    end
end
end

function v = versions()
v = struct('matlab', version, 'eeglab', 'unknown', 'pipecompare', pipecompare.PipeCompare.Version, ...
    'iclabel', which('pop_iclabel'), 'firfilt', which('pop_eegfiltnew'));
try, v.eeglab = eeg_getversion; catch, end
end

function elig = eligibleTrials(root, contract)
% urevent ids of the time-locking events the trial rule keeps ([] = all).
elig = [];
r = contract.trials;
if strcmp(r.mode, 'all'), return; end
codes = contract.allEvents();
ev = root.event;
types = arrayfun(@(e) strtrim(char(string(e.type))), ev, 'UniformOutput', false);
isCode = ismember(types, codes);
ure = arrayfun(@(e) double(e.urevent), ev);
switch r.mode
    case 'urevents'
        keep = isCode & ismember(ure, r.ids);
    case 'time_ranges'
        t = ([ev.latency] - 1) / root.srate; keep = false(size(t));
        for k = 1:size(r.ranges, 1), keep = keep | (t >= r.ranges(k, 1) & t < r.ranges(k, 2)); end
        keep = keep & isCode;
    case 'marker_ranges'
        active = false; keep = false(size(types));
        for k = 1:numel(types)
            if strcmp(types{k}, char(string(r.startCode))), active = true;
            elseif strcmp(types{k}, char(string(r.endCode))), active = false;
            elseif isCode(k), keep(k) = active; end
        end
end
elig = unique(ure(keep));
assert(~isempty(elig), 'PipeCompare:TrialRule', 'The trial rule keeps no trials.');
end


function EEG = evalWithEEG(EEG, PIPECOMPARE_CMD__)
% Run a script line on EEG exactly as the exported script will.
eval(PIPECOMPARE_CMD__);
end

function id = searchIdentity(result)
% What makes two searches the same: the starting data (full content), the
% analysis contract, the enumerated pipelines and the options that change
% results. A checkpoint folder belongs to exactly one identity.
o = result.options;
o = rmfield(o, intersect(fieldnames(o), {'checkpoint','stopAfter','maxSeconds','parallel','verbose','dryRun','progress'}));
md = java.security.MessageDigest.getInstance('SHA-256');
root = result.root;
% the data in blocks: the same bytes as double(root.data(:)), without a
% double copy of the whole recording in memory
N = numel(root.data); step = 2^22;
for i0 = 1:step:N
    md.update(typecast(double(root.data(i0:min(N, i0 + step - 1))), 'uint8'));
end
parts = {getByteStreamFromArray(result.rootComs), getByteStreamFromArray(result.contract.toStruct()), ...
    getByteStreamFromArray(result.labels), getByteStreamFromArray(o), ...
    getByteStreamFromArray({root.chanlocs.labels}), getByteStreamFromArray(root.srate), ...
    getByteStreamFromArray(arrayfun(@(e) {char(string(e.type)), double(e.latency)}, root.event, 'UniformOutput', false))};
for k = 1:numel(parts), md.update(parts{k}); end
id = lower(reshape(dec2hex(typecast(md.digest(), 'uint8'))', 1, []));
end


