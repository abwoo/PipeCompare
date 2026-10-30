function pipecompare_component_audit(rootDir,outDir)
% Source-resolved fixed-operator audit of the two existing development cases.
% Truth is used exclusively for diagnosis, never for pipeline selection.
if ~isfolder(outDir),mkdir(outDir);end
d=load(fullfile(rootDir,'guard','reference_guard_replay.mat'),'reports');
C=readtable(fullfile(rootDir,'summary','balanced_calibration.csv'),'TextType','string','VariableNamingRule','preserve');
allReports={};
for seed=[40102 40103]
    rr=d.reports(cellfun(@(r)strcmp(r.population,'synthetic')&&r.record==seed,d.reports));
    assert(numel(rr)==1);rr=rr{1};
    manifests=dir(fullfile(rootDir,sprintf('seed%d',seed),'**','manifest.mat'));
    manifests=manifests(endsWith({manifests.folder},[filesep 'search']));assert(numel(manifests)==1);
    folder=manifests(1).folder;M=load(fullfile(folder,'manifest.mat'),'result','identity');
    R=load(fullfile(folder,'root.mat'),'root','identity');assert(strcmp(R.identity,M.identity));r=M.result;
    [generated,truth,sources]=bench_synth(seed); %#ok<ASGLU>
    assert(isequal(generated.data,R.root.data),'PipeCompare:ComponentAudit','Simulation bytes changed.');
    assert(isequal([generated.event.latency],[R.root.event.latency]),'PipeCompare:ComponentAudit','Event identities changed.');
    files=dir(fullfile(folder,'leaf_*.mat'));cands=[];globalIds=[];
    for q=1:numel(files)
        x=load(fullfile(files(q).folder,files(q).name),'cand','id','key','identity');
        assert(strcmp(x.identity,M.identity)&&strcmp(x.key,r.leaves(x.id).key));
        if isempty(cands),cands=x.cand;else,cands(end+1)=x.cand;end %#ok<AGROW>
        globalIds(end+1)=x.id; %#ok<AGROW>
    end
    [globalIds,order]=sort(globalIds);cands=cands(order);
    assert(strcmp(cands(rr.proposed).key,rr.originalKey)&&strcmp(cands(rr.selected).key,rr.selectedKey));
    ids=[rr.proposed rr.selected];roles={'original','guarded'};
    a=counterpart(cands,cands(rr.proposed).key,'badchannels',true);
    if ~isempty(a),ids(end+1)=a;roles{end+1}='repair_only_counterpart';end
    a=counterpart(cands,cands(rr.selected).key,'icremove',false);
    if ~isempty(a),ids(end+1)=a;roles{end+1}='without_ic_removal_of_guarded';end
    frozen=C(C.recording==seed,:);population=zeros(1,numel(r.ref.ids));realized=population;
    for q=1:numel(population)
        v=unique(frozen.populationTruth(frozen.condition==q));assert(isscalar(v));population(q)=v;
        v=unique(frozen.realizedTruth(frozen.condition==q));assert(isscalar(v));realized(q)=v;
    end
    clean=R.root;clean.data=sources.p3;clean.icaact=[];
    path=r.leaves(globalIds(rr.proposed)).path;types=cellfun(@(v)v.type,path,'UniformOutput',false);
    minimal=path(ismember(types,{'reref','epoch','baseline'}));
    assert(nnz(strcmp(types,'reref'))==1);
    ctx=struct('contract',r.contract,'ref',r.ref,'rank',r.options,'highpass',0);
    for q=1:numel(minimal),clean=pipecompare.run.Steps.run(minimal{q},clean,ctx);end
    cleanM=pipecompare.eval.Measure.candidate(clean,r.contract,r.ref,r.options);
    assert(max(abs(cleanM.objectives(1).estimate-realized))<1e-5,'PipeCompare:ComponentAudit','Realized clean target mismatch.');
    for id=unique(ids,'stable')
        cand=cands(id);path=r.leaves(globalIds(id)).path;
        role=roles(ids==id);
        [observed,replayScope]=pipecompare_fixed_replay(R.root,path,cand,r.contract,r.ref,r.options);
        measured=pipecompare.eval.Measure.candidate(observed,r.contract,r.ref,r.options);clear observed
        o=cand.m.objectives(1);
        assert(max(abs(measured.objectives(1).estimate-o.estimate))<2e-5,'PipeCompare:ComponentAudit','Cached estimate was not reproduced.');
        g=cand.signal.gain(1);assert(isfinite(g)&&g>0);
        sourceNames=fieldnames(sources);sourceEst=zeros(numel(sourceNames),numel(population));
        for q=1:numel(sourceNames)
            part=R.root;part.data=sources.(sourceNames{q});part.icaact=[];
            part=pipecompare_fixed_replay(part,path,cand,r.contract,r.ref,r.options);
            pm=pipecompare.eval.Measure.candidate(part,r.contract,r.ref,r.options);
            sourceEst(q,:)=pm.objectives(1).estimate/g;
            clear part pm
        end
        closure=sum(sourceEst,1)-o.estimate/g;
        assert(max(abs(closure))<5e-5,'PipeCompare:ComponentAudit','Fixed operator does not close additively.');
        retained=zeros(size(population));
        for q=1:numel(retained)
            x=cleanM.objectives(1).X{q}(:,1);keep=isfinite(o.X{q}(:,1));
            w=ones(size(x));if isfield(o,'w')&&~isempty(o.w),w=o.w{q};end
            keep=keep&isfinite(w)&w>0;retained(q)=sum(w(keep).*x(keep))/sum(w(keep));
        end
        at=find(strcmp(sourceNames,'p3'));other=setdiff(1:numel(sourceNames),at,'stable');
        partNames=[{'finite_sample','retention_weighting','signal_transform'},sourceNames(other)'];
        errorParts=[realized-population;retained-realized;sourceEst(at,:)-retained;sourceEst(other,:)];
        actual=o.estimate/g-population;
        assert(max(abs(sum(errorParts,1)-actual))<5e-5);
        report=struct('record',seed,'roles',{role},'candidateId',globalIds(id),'key',cand.key, ...
            'population',population,'realized',realized,'estimate',o.estimate/g,'error',actual, ...
            'mse',mean(actual.^2),'nominalGain',g,'trueSignalGain',sourceEst(at,:)*g./retained, ...
            'parts',{partNames},'signedErrorParts',errorParts,'squaredErrorAllocation',mean(errorParts.*actual,2), ...
            'eyeAssociation',o.eyeBias/g,'eyeAssociationSE',o.eyeBiasSE/g,'softBiasRisk',sqrt(o.eyeBiasSq)/g, ...
            'kept',cand.m.kept,'retention',cand.m.retention,'closureError',max(abs(closure)), ...
            'sourceRun',40282756637,'sourceRevision','4116ea11f3d5d39e4d3d282894dcf487172bd50d', ...
            'scope',replayScope,'interpretation','Fixed fitted-operator development diagnosis; signed allocations include cross terms and can be negative.');
        allReports{end+1}=report; %#ok<AGROW>
        fprintf('COMPONENT_AUDIT %s\n',jsonencode(report));
    end
    clear sources cands R generated clean cleanM
end
save(fullfile(outDir,'component_audit.mat'),'allReports','-v7');
fid=fopen(fullfile(outDir,'component_audit.json'),'w');cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(allReports,PrettyPrint=true));
fprintf('COMPONENT_AUDIT_COMPLETE=%d; ICA_FITS=0\n',numel(allReports));
end

function id=counterpart(cands,key,type,present)
target=without(key,type);id=[];
for q=1:numel(cands)
    if ~strcmp(cands(q).status,'ok'),continue;end
    parts=strsplit(cands(q).key,' > ');has=any(startsWith(parts,[type '('])|strcmp(parts,type));
    if has~=present||~strcmp(without(cands(q).key,type),target),continue;end
    id=q;return;
end
end
function value=without(key,type)
parts=strsplit(key,' > ');parts=parts(~startsWith(parts,[type '('])&~strcmp(parts,type));value=strjoin(parts,' > ');
end
