function d=p3_hmm_veog_diagnostic(D,support,varargin)
%P3_HMM_VEOG_DIAGNOSTIC Exploratory v2: full-condition VEOG/HMM alignment.
% Diagnose ERP CORE P3 eye proxy only; does not reject trials, refit ICA,
% change P3 results, or change original experiment's inclusion criteria.
% Uses baseline-corrected max(abs(VEOG)) > 50 uV as a PROXY, not blink truth.
% Source event indices are 1-based rows of D.eventTypes, not original urevents.
if ~isfield(D,'fs')||double(D.fs)~=256
 error('PR51:VeogDiagnostic','Expected frozen 256 Hz excerpt.');
end
fs=double(D.fs);veog=double(D.veog(:)');
lat=double(D.eventLatency(:));types=cellstr(D.eventTypes(:));
assert(numel(lat)==numel(types),'PR51:VeogDiagnostic','Event rows mismatch.');
offsets=-round(0.2*fs):round(0.8*fs);
baseline=offsets<0;
assert(any(baseline),'PR51:VeogDiagnostic','Baseline interval empty.');
d=struct('available',true,'diagnosticVersion',2,'veogThresholdUv',50, ...
 'nFlagged',0,'nUnflagged',0,'flaggedWithVeog',0,'unflaggedWithVeog',0, ...
 'flaggedMedianMaxVeog',NaN,'unflaggedMedianMaxVeog',NaN, ...
 'conditionScope','all supplied target and standard support rows', ...
 'proxyIsNotBlinkGroundTruth',true,'byCondition',[]);
flaggedMax=[];unflaggedMax=[];
by=struct('condition',{},'referenceN',{},'matchedN',{},'notMatchedN',{}, ...
 'flaggedN',{},'unflaggedN',{},'flaggedWithVeog',{}, ...
 'unflaggedWithVeog',{},'flaggedMedianMaxVeog',{},'unflaggedMedianMaxVeog',{});
for c=1:numel(support)
 label=char(string(support(c).condition));
 switch lower(label)
  case {'11','target'}, code='11';
  case {'31','standard'}, code='31';
  otherwise
   error('PR51:VeogDiagnostic','Unexpected condition label: %s',label);
 end
 ids=double(support(c).referenceIds(:));
 flags=logical(support(c).flagged(:));
 assert(numel(ids)==numel(flags),'PR51:VeogDiagnostic','Flag/reference index mismatch.');
 assert(numel(unique(ids))==numel(ids),'PR51:VeogDiagnostic','Duplicate reference indices.');
 rows=find(strcmp(types,code));
 fvals=[];uvals=[];nMatch=0;
 for z=1:numel(rows)
  j=rows(z);
  pos=find(ids==j,1);
  if isempty(pos),continue;end
  nMatch=nMatch+1;
  cols=lat(j)+offsets;
  if any(cols<1)||any(cols>numel(veog))
   error('PR51:VeogDiagnostic','Incomplete VEOG window at event %d',j);
  end
  seg=veog(cols);seg=seg-mean(seg(baseline));
  mv=max(abs(seg));
  if ~isfinite(mv)
   error('PR51:VeogDiagnostic','Nonfinite VEOG at event %d',j);
  end
  if flags(pos),fvals(end+1)=mv; %#ok<AGROW>
  else,uvals(end+1)=mv;end %#ok<AGROW>
 end
 if nMatch~=numel(ids)
  error('PR51:VeogDiagnostic','Source support event IDs do not align with condition %s',code);
 end
 q=struct('condition',code,'referenceN',numel(ids),'matchedN',nMatch, ...
   'notMatchedN',numel(ids)-nMatch,'flaggedN',numel(fvals), ...
   'unflaggedN',numel(uvals),'flaggedWithVeog',sum(fvals>50), ...
   'unflaggedWithVeog',sum(uvals>50), ...
   'flaggedMedianMaxVeog',median(fvals),'unflaggedMedianMaxVeog',median(uvals));
 by(end+1)=q; %#ok<AGROW>
 flaggedMax=[flaggedMax fvals];unflaggedMax=[unflaggedMax uvals]; %#ok<AGROW>
end
d.byCondition=by;
d.nFlagged=numel(flaggedMax);d.nUnflagged=numel(unflaggedMax);
d.flaggedWithVeog=sum(flaggedMax>d.veogThresholdUv);
d.unflaggedWithVeog=sum(unflaggedMax>d.veogThresholdUv);
d.flaggedMedianMaxVeog=median(flaggedMax);
d.unflaggedMedianMaxVeog=median(unflaggedMax);
end
