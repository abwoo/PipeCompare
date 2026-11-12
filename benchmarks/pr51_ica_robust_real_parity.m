function pr51_ica_robust_real_parity(evidenceDir)
% Demonstrate that Matlab IcaRobust and Python shadow interpret the same
% REAL development EEG source recordings and original scientific gates.
% No data fitting, tuning, eligibility change or frozen study launched.
assert(nargin==1 && isfolder(evidenceDir));
files=dir(fullfile(evidenceDir,'**','result.json'));
assert(numel(files)==3,'PR51:Parity','Expected three preserved real ERP CORE P3 jobs.');
expected=containers.Map([11 15 21],[25 2 4]);
seen=zeros(1,3);
for k=1:numel(files)
    p=fullfile(files(k).folder,files(k).name);
    r=jsondecode(fileread(p));
    assert(r.completed && isempty(r.errors),'PR51:Parity','Source P3 record incomplete.');
    subj=double(r.subject);
    assert(isKey(expected,subj),'PR51:Parity','Unexpected source cohort subject.');
    z=r.icaSeedCheck;
    n=numel(z.ids); seeds=size(z.objective,2);
    assert(seeds==3 && numel(z.rankStatus)==n*seeds && ...
           numel(z.strata)==n*seeds && ...
           numel(z.candidateMeasurements)==n*seeds, ...
           'PR51:Parity','Frozen cell serialization shape not verified.');
    if isstring(z.rankStatus), z.rankStatus=cellstr(z.rankStatus); end
    if isstring(z.strata), z.strata=cellstr(z.strata); end
    z.rankStatus=reshape(z.rankStatus,n,seeds);
    z.strata=reshape(z.strata,n,seeds);
    z.candidateMeasurements=reshape(num2cell(z.candidateMeasurements),n,seeds);
    q=pipecompare.eval.IcaRobust.preview(z);
    assert(q.id==expected(subj),'PR51:Parity', ...
           'MATLAB robust shadow differs from Python real-EEG preview subject %d.',subj);
    assert(~q.candidateSelectedForProduction && q.existingRankRecommendationUnchanged);
    fprintf('PR51_ICA_MATLAB_REAL_PARITY subject=%d original=%s shadow=%d eligible=%d\n', ...
        subj,mat2str(z.recommended),q.id,numel(q.eligibleIds));
    seen(k)=subj;
end
assert(isequal(sort(seen),[11 15 21]),'PR51:Parity','Missing source participant.');
end
