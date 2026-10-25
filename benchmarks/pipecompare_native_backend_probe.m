function pipecompare_native_backend_probe(outDir)
% Verify CPU MEX inference against ICLabel's original MATLAB fallback.
if ~isfolder(outDir),mkdir(outDir);end
plugin=fileparts(which('pop_iclabel'));assert(~isempty(plugin),'ICLabel is required.');
addpath(fullfile(plugin,'matconvnet','matlab'));vl_setupnn;
mex('-setup','C++');
vl_compilenn('enableGpu',false,'enableImreadJpeg',false,'verbose',1);
vl_setupnn;rehash;
assert(exist('vl_nnconv','file')==3,'PipeCompare:Backend','CPU convolution MEX was not activated.');
rng(734821,'twister');n=8;
images=single(.1*randn(32,32,1,n));psds=single(.1*randn(1,100,1,n));autocorrs=single(.1*randn(1,100,1,n));
tic;compiled=run_ICL('default',images,psds,autocorrs);compiledSeconds=toc;
netStruct=load(fullfile(plugin,'netICL'));
try,net=dagnn.DagNN.loadobj(netStruct);catch,net=dagnn_bc.DagNN.loadobj(netStruct);end
aug=cat(4,images,-images,images(:,end:-1:1,:,:),-images(:,end:-1:1,:,:));
input={'in_image',aug,'in_psdmed',repmat(psds,[1 1 1 4]),'in_autocorr',repmat(autocorrs,[1 1 1 4])};
tic;net=uncompiled_network_evaluation(net,input);fallbackSeconds=toc;
labels=squeeze(net.getVar(net.getOutputs()).value)';
fallback=reshape(mean(reshape(labels',[],4),2),7,[])';
delta=max(abs(double(compiled(:))-double(fallback(:))));
[~,a]=max(compiled,[],2);[~,b]=max(fallback,[],2);
assert(delta<1e-5&&isequal(a,b),'PipeCompare:Backend','CPU MEX and MATLAB fallback classifications differ.');
report=struct('components',n,'maxProbabilityDifference',delta,'sameWinningClasses',isequal(a,b), ...
 'compiledSeconds',compiledSeconds,'fallbackSeconds',fallbackSeconds,'mex',which('vl_nnconv'), ...
 'matlab',version,'matconvnetCommit','092ecb8c04031ddcb1d9a7ba6ebe76a250c4672c');
fid=fopen(fullfile(outDir,'backend.json'),'w');fprintf(fid,'%s\n',jsonencode(report));fclose(fid);
fprintf('PIPECOMPARE_NATIVE_BACKEND_JSON %s\n',jsonencode(report));
end
