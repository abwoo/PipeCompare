function P=validation_precision(S)
% Prespecified pilot precision, distinct from calibration agreement or
% distributional stability. A precise biased ratio still fails calibration.
P=struct('relativeHalfWidth',NaN,'targetRelativeHalfWidth',0.10, ...
    'minimumRecords',30,'meetsPilotPrecision',false, ...
    'scope','Pooled RMS-ratio precision only; not proof of calibration or stability across settings');
if isfinite(S.ratio)&&S.ratio>0&&all(isfinite(S.ci))
    P.relativeHalfWidth=diff(S.ci)/(2*S.ratio);
    P.meetsPilotPrecision=S.nRecords>=P.minimumRecords&&P.relativeHalfWidth<=P.targetRelativeHalfWidth;
end
end
