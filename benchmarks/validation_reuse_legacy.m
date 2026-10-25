function used=validation_reuse_legacy(legacyDir,outDir,kind,id,population,realized)
% Reuse completed, frozen legacy controls; never reuse a partial fit.
used=false;source=fullfile(legacyDir,'standard');
if ~isfolder(source)||~isfile(fullfile(source,'execution.json')),return;end
execution=jsondecode(fileread(fullfile(source,'execution.json')));
if ~execution.completed,return;end
if strcmp(kind,'public')
    meta=jsondecode(fileread(fullfile(legacyDir,'metadata.json')));
    assert(meta.subject==id&&strcmp(meta.stage,'development'),'PipeCompare:Validation','Wrong legacy participant.');
    oldInput=jsondecode(fileread(fullfile(legacyDir,'frozen_inputs.json')));
    currentInput=jsondecode(fileread(fullfile(fileparts(outDir),'frozen_inputs.json')));
    a=oldInput.files([oldInput.files.subject]==id);
    b=currentInput.files([currentInput.files.subject]==id);
    assert(isscalar(a)&&isscalar(b)&&strcmp(a.path,b.path)&&strcmp(a.sha256,b.sha256), ...
        'PipeCompare:Validation','Legacy input hash does not match the current participant.');
    C=readtable(fullfile(source,'calibration.csv'),'TextType','char');
    for c=1:numel(population)
        rows=C.condition==c;
        assert(all(abs(C.populationTruth(rows)-population(c))<1e-9)&& ...
            all(abs(C.realizedTruth(rows)-realized(c))<1e-9),'PipeCompare:Validation','Legacy truth differs from the current fixture.');
    end
    execution.sourceCommit=meta.commit;
else
    assert(execution.seed==id&&strcmp(execution.recipe,'standard'),'PipeCompare:Validation','Wrong legacy seed.');
end
if ~isfolder(outDir),mkdir(outDir);end
parts=dir(source);
for k=1:numel(parts)
    if ~parts(k).isdir,copyfile(fullfile(source,parts(k).name),fullfile(outDir,parts(k).name));end
end
execution.legacyReused=true;execution.sourceRun=40153112163;
execution.artifactVariance='legacy-total';execution.currentFitSeconds=0;
fid=fopen(fullfile(outDir,'execution.json'),'w');fprintf(fid,'%s\n',jsonencode(execution));fclose(fid);
fprintf('PIPECOMPARE_REUSED_LEGACY kind=%s id=%d source=40153112163\n',kind,id);
used=true;
end
