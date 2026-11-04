function run_waveform_audit(input,folder,outDir)
if ~isfolder(outDir),mkdir(outDir);end
[E,pop,real,truth,control]=waveform_input(input,5);
d=load(fullfile(folder,'baseline','after_eog.mat'));
assert(max(abs(d.cfg.population(:)-pop(:)))<1e-12&&max(abs(d.cfg.realized(:)-real(:)))<1e-12);
c=pipecompare.simple.Presets.contract(E,'P3',{'11','31'});ref=pipecompare.eval.Measure.reference(E,c);
assert(isequal(ref.ids,d.ref.ids)&&isequal(ref.blink,d.ref.blink),'Input/reference changed');
audit_hmm(E,control,truth,outDir);
audit_rejection(E,folder,outDir);
fprintf('WAVEFORM_AUDIT_CHECKS_PASSED hmmCases=3 selectedRecipes=2 seed=0\n');
end
