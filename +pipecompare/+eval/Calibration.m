classdef Calibration
    % Calibration of the gain-corrected ranking error against independent
    % simulated records. Records, not candidate/condition rows, are the
    % sampling units. This measures RMS error, not interval coverage.
    methods (Static)
        function C = rows(r, feasible, population, realized, recording)
            data = cell(0,18); T = r.ranking.table;
            vars = T.Properties.VariableNames;
            gv = vars(startsWith(vars, 'gain_')); sv = vars(startsWith(vars, 'shapeSpread_'));
            for j = feasible(:)'
                o = r.cands(j).m.objectives(1);
                g = 1; if ~isempty(gv), g = T.(gv{1})(j); end
                sp = 0; if ~isempty(sv), sp = T.(sv{1})(j); end
                if ~(isfinite(g) && g > 0), continue; end
                [reported, ~, bias] = pipecompare.eval.Rank.conditionErrors(o, g, sp);
                selected = any(r.ranking.recommended == j);
                for c = 1:numel(population)
                    data(end+1,:) = {recording, r.labels{j}, c, o.estimate(c)/g-population(c), ...
                        reported(c), o.sme(c), o.eyeBias(c), g, selected, o.estimate(c)/g-realized(c), sp, bias(c), ...
                        o.estimate(c)-population(c), o.estimate(c)-realized(c), o.estimate(c), o.estimate(c)/g, ...
                        population(c), realized(c)}; %#ok<AGROW>
                end
            end
            C = cell2table(data, 'VariableNames', {'recording','steps','condition','actualError','reportedError', ...
                'sme','eyeBias','gain','selected','conditionalError','shapeSpread','totalBias', ...
                'rawError','rawConditionalError','rawEstimate','correctedEstimate','populationTruth','realizedTruth'});
        end

        function S = summarize(C, nBoot, seed)
            if nargin < 2, nBoot = 1999; end
            if nargin < 3, seed = 7407; end
            records = unique(C.recording);
            a = []; p = [];
            for id = records(:)'
                R = C(C.recording == id,:);
                ok = isfinite(R.actualError) & isfinite(R.reportedError) & R.reportedError >= 0;
                % Do not silently drop a condition with undefined bias.
                if ~all(ok), continue; end
                a(end+1,1) = mean(R.actualError.^2); p(end+1,1) = mean(R.reportedError.^2); %#ok<AGROW>
            end
            S = struct('ratio',NaN,'ci',[NaN NaN],'actualRms',NaN,'reportedRms',NaN, ...
                'nRecords',numel(a),'nRows',height(C),'unassessedRecords',numel(records)-numel(a), ...
                'exploratory',numel(a)<30,'interpretation','undefined');
            if isempty(a) || mean(p) <= 0, return; end
            S.actualRms = sqrt(mean(a)); S.reportedRms = sqrt(mean(p));
            S.ratio = S.actualRms/S.reportedRms;
            if S.ratio > 1, S.interpretation = 'reported error understated';
            elseif S.ratio < 1, S.interpretation = 'reported error overstated';
            else, S.interpretation = 'RMS errors agree'; end
            if numel(a) < 2 || nBoot < 20, return; end
            stream = RandStream('mt19937ar','Seed',seed);
            I = randi(stream,numel(a),numel(a),nBoot);
            q = sqrt(mean(a(I),1)./mean(p(I),1)); q = sort(q(isfinite(q)));
            if numel(q) < 20, return; end
            at = max(1,min(numel(q),ceil((numel(q)+1)*[0.025 0.975])));
            S.ci = q(at);
        end

        function t = describe(C, label)
            S = pipecompare.eval.Calibration.summarize(C);
            t = sprintf('  %-34s %.2f, record-bootstrap 95%% CI [%.2f, %.2f]; %d records, %d rows; %s', ...
                label,S.ratio,S.ci(1),S.ci(2),S.nRecords,S.nRows,S.interpretation);
            if S.exploratory, t = [t '; exploratory sample size']; end
            if S.unassessedRecords > 0, t = sprintf('%s; %d records unassessed',t,S.unassessedRecords); end
        end
    end
end
