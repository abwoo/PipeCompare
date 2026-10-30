function pipecompare_shadow_summary_fixtures(outDir)
% MATLAB/Python wire-format fixtures: zero, one and multiple complete pairs.
cases = {'zero','single','multiple'};
for k = 1:numel(cases)
    n = 1; if k == 3, n = 3; end
    arms = struct([]);
    for j = 1:n
        cand = struct('m',struct('objectives',struct('estimate',[1 1])), ...
            'signal',struct('gain',[1 1]));
        pick = 1; if k == 1, pick = []; end
        a = struct('seed',20030+j,'truth',[0 0], ...
            'standardCands',cand,'standardPick',1,'balancedCands',cand, ...
            'balancedGuardedPick',pick,'balancedUnguardedPick',pick, ...
            'guardDecision','maintain');
        if isempty(arms), arms = a; else, arms(end+1) = a; end %#ok<AGROW>
    end
    dest = fullfile(outDir,cases{k}); mkdir(dest);
    out = pipecompare_shadow_selection(dest,20031,n,'mode','test','testArms',arms, ...
        'headSha','fixture1234567890','workflowRunId','fixture');
    run = jsondecode(fileread(out.runJson));
    fragments = fullfile(dest,'fragments'); mkdir(fragments);
    for j = 1:numel(out.perSeed)
        e = out.perSeed(j);
        f = struct('seed',e.seed,'arms',e.arms,'guardChanged',e.guardChanged, ...
            'guardDecision',e.guardDecision,'labels',run.labels,'provenance',run.provenance);
        file = fullfile(fragments,sprintf('shadow_selection_seed_%d.json',e.seed));
        fid = fopen(file,'w'); assert(fid >= 0);
        fprintf(fid,'%s\n',jsonencode(f)); fclose(fid);
    end
end
end
