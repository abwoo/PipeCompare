function tests = test_p3_hmm_veog_diagnostic
% The old ERP CORE summary silently omitted target ('11') events because
% it only checked for the literal condition name 'target'. This test
% guards the new versioned diagnostic without changing the frozen pilot.
tests = functiontests(localfunctions);
end

function setupOnce(~)
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'benchmarks'));
end

function testBothNumericConditionCodes(tc)
[D, support] = fixture();
d = p3_hmm_veog_diagnostic(D, support);
verifyTrue(tc, d.available);
verifyEqual(tc, d.diagnosticVersion, 2);
verifyEqual(tc, d.nFlagged, 1);  % target '11', which legacy missed
verifyEqual(tc, d.nUnflagged, 1);
verifyEqual(tc, d.flaggedWithVeog, 1);
verifyEqual(tc, d.unflaggedWithVeog, 0);
verifyEqual(tc, numel(d.byCondition), 2);
verifyEqual(tc, d.byCondition(1).condition, '11');
verifyEqual(tc, d.byCondition(1).flaggedN, 1);
verifyEqual(tc, d.byCondition(2).condition, '31');
verifyEqual(tc, d.byCondition(2).unflaggedN, 1);
verifyTrue(tc, d.proxyIsNotBlinkGroundTruth);
end

function testAliasConditionNamesPreserveValues(tc)
[D, support] = fixture();
support(1).condition = 'target';
support(2).condition = 'standard';
a = p3_hmm_veog_diagnostic(D,support);
verifyEqual(tc,a.nFlagged,1);
verifyEqual(tc,a.nUnflagged,1);
verifyEqual(tc,a.flaggedWithVeog,1);
end

function testUnrecognizedConditionFails(tc)
[D,support] = fixture();
support(1).condition = 'undefined';
verifyError(tc,@()p3_hmm_veog_diagnostic(D,support),'PR51:VeogDiagnostic');
end

function testMissingReferenceIsNotSilentlyDiscarded(tc)
[D,support] = fixture();
support(1).referenceIds = 99;
verifyError(tc,@()p3_hmm_veog_diagnostic(D,support),'PR51:VeogDiagnostic');
end

function testDuplicateReferenceIndexFails(tc)
[D,support] = fixture();
support(1).referenceIds = [1 1];
support(1).flagged = [true true];
verifyError(tc,@()p3_hmm_veog_diagnostic(D,support),'PR51:VeogDiagnostic');
end

function testProxyThresholdAndBaselineCorrection(tc)
[D,support] = fixture();
% Add a large DC offset without changing the baseline-subtracted result.
base = p3_hmm_veog_diagnostic(D,support);
D.veog = D.veog + 300;
shifted = p3_hmm_veog_diagnostic(D,support);
verifyEqual(tc,base.flaggedWithVeog,shifted.flaggedWithVeog);
verifyEqual(tc,base.nFlagged,shifted.nFlagged);
verifyEqual(tc,base.flaggedMedianMaxVeog,shifted.flaggedMedianMaxVeog,'AbsTol',1e-9);
end

function [D,support] = fixture()
D = struct('fs',256,'veog',zeros(1,2200), ...
    'eventLatency',[500,1500],'eventTypes',{{'11','31'}});
D.veog(520) = 80;
support(1) = struct('condition','11','referenceIds',1,'flagged',true);
support(2) = struct('condition','31','referenceIds',2,'flagged',false);
end
