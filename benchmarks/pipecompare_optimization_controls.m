function pipecompare_optimization_controls(dataDir,outDir)
% Reused annotated/real-task component controls, not additional participants.
mkdir(outDir);protocol=jsondecode(fileread(fullfile(fileparts(mfilename('fullpath')),'open_data_protocol.json')));rows=table();files=dir(fullfile(dataDir,'S*.mat'));
files=files(~cellfun('isempty',regexp({files.name},'^S\d\d\.mat$','once')));
assert(numel(files)==15,'PipeCompare:Validation','Missing blink controls.');
for k=1:numel(files)
    D=load(fullfile(files(k).folder,files(k).name));
    T=pipecompare_open_blinks(D,files(k).name,protocol.blink.matching_margin_seconds);
    if isempty(rows),rows=T;else,rows=[rows;T];end %#ok<AGROW>
end
writetable(rows,fullfile(outDir,'blink.csv'));
files=dir(fullfile(dataDir,'erp_core','*.set'));assert(numel(files)==6,'PipeCompare:Validation','Missing ERP controls.');
rows=table();
for k=1:numel(files)
    E=pipecompare_load_public_set(fullfile(files(k).folder,files(k).name));
    T=pipecompare_open_real_erp(E,files(k).name,outDir);
    if isempty(rows),rows=T;else,rows=[rows;T];end %#ok<AGROW>
end
writetable(rows,fullfile(outDir,'real_erp.csv'));
fprintf('PIPECOMPARE_COMPONENT_CONTROLS_COMPLETE blinkParticipants=15 ERPParticipants=6\n');
end
