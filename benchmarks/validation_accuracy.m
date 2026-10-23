function S=validation_accuracy(C,field)
% Accuracy for a declared truth/readout axis; never invent a reported-error
% calibration for an axis whose uncertainty model has a different target.
records=unique(C.recording); mse=[]; bias=[];
for id=records(:)'
    x=C.(field)(C.recording==id);
    if isempty(x) || any(~isfinite(x)),continue;end
    mse(end+1,1)=mean(x.^2);bias(end+1,1)=mean(x); %#ok<AGROW>
end
S=struct('field',field,'nRecords',numel(mse),'unassessedRecords',numel(records)-numel(mse), ...
    'rmse',NaN,'bias',NaN,'interpretation','accuracy only; no matched reported-error claim');
if ~isempty(mse),S.rmse=sqrt(mean(mse));S.bias=mean(bias);end
end
