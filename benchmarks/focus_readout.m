function a=focus_readout(r,k,population,realized)
a=emptyReadout();o=r.cands(k).m.objectives(1);
gcol=matlab.lang.makeValidName(['gain_' r.ref.objectives{1}]);
scol=matlab.lang.makeValidName(['shapeSpread_' r.ref.objectives{1}]);
g=r.ranking.table.(gcol)(k);sp=r.ranking.table.(scol)(k);
if ~(isfinite(g)&&g>0&&all(isfinite(o.estimate))),return;end
[reported,~,~]=pipecompare.eval.Rank.conditionErrors(o,g,sp);
a.available=true;a.estimate=o.estimate;a.correctedEstimate=o.estimate/g;a.gain=g;a.pipeline=r.labels{k};
a.actualError=a.correctedEstimate-population;a.rawError=o.estimate-population;
a.conditionalError=a.correctedEstimate-realized;a.reportedError=reported;
a.mse=mean(a.actualError.^2);a.rawMse=mean(a.rawError.^2);
a.conditionalMse=mean(a.conditionalError.^2);a.reportedMse=mean(reported.^2);
end


function a=emptyReadout()
a=struct('available',false,'estimate',[],'correctedEstimate',[],'actualError',[], ...
 'rawError',[],'conditionalError',[],'reportedError',[],'mse',NaN,'rawMse',NaN, ...
 'conditionalMse',NaN,'reportedMse',NaN,'gain',NaN,'pipeline','');
end
