function [loss,contrast,inFamily]=validation_recommended_loss(C,S,keys)
% Evaluate the frozen per-record recommendation, rather than a raw minimum.
% Truth rows are read only after selection. The actual recommendation may
% use a wider candidate family than the training-supported group comparison.
loss=NaN;contrast=NaN;inFamily=false;
id=sscanf(S.name,'S%dR01.mat');
D=C(C.recording==id&C.selected==1,:);
if isempty(D),return;end
conditions=unique(D.condition);errors=nan(size(conditions));raw=errors;
for j=1:numel(conditions)
    at=D.condition==conditions(j);values=D.actualError(at);
    assert(all(isfinite(values))&&max(values)-min(values)<1e-8, ...
        'PipeCompare:Validation','Selected truth rows are missing or inconsistent.');
    errors(j)=mean(values);raw(j)=mean(D.rawEstimate(at));
end
loss=mean(errors.^2);
if isequal(conditions(:)',[1 2]),contrast=errors(1)-errors(2);end
names=S.table.Properties.VariableNames;gains=names(startsWith(names,'gain_'));
found=false;
for k=1:numel(S.cands)
    if ~strcmp(S.labels{k},D.steps{1})||~isstruct(S.cands(k).m)|| ...
            ~isfield(S.cands(k).m,'objectives')||isempty(S.cands(k).m.objectives)
        continue;
    end
    estimate=S.cands(k).m.objectives(1).estimate(:);
    gain=1;if ~isempty(gains),gain=S.table.(gains{1})(k);end
    if numel(estimate)==numel(raw)&&all(isfinite(estimate))&& ...
            max(abs(estimate-raw(:)))<=1e-9*max(1,max(abs(raw)))&& ...
            abs(gain-D.gain(1))<=1e-9*max(1,abs(gain))
        found=true;inFamily=inFamily||ismember(S.keys(k),keys);
    end
end
assert(found,'PipeCompare:Validation','Cannot match the actual recommendation to its candidate metadata.');
end
