function pipecompare_optimization_synthetic(outDir,seed,stage)
% Same seed and full native search for both legacy and expanded recipes.
protocol=jsondecode(fileread(fullfile(fileparts(mfilename('fullpath')),'optimization_protocol.json')));
assert(ismember(seed,protocol.([stage '_synthetic_seeds'])),'PipeCompare:Validation','Seed outside frozen stage.');
for recipe={'standard','balanced'}
    folder=fullfile(outDir,recipe{1});if ~isfolder(folder),mkdir(folder);end
    pipecompare_benchmark(1,folder,seed,recipe{1});
    if ~isfile(fullfile(folder,'calibration.csv'))
        writetable(pipecompare.eval.Calibration.empty(),fullfile(folder,'calibration.csv'));
    end
    fid=fopen(fullfile(folder,'execution.json'),'w');
    fprintf(fid,'%s\n',jsonencode(struct('completed',true,'seed',seed,'recipe',recipe{1})));
    fclose(fid);
end
end
