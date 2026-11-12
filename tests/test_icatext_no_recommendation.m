function tests = test_icatext_no_recommendation
% No viable pipeline is a legitimate research outcome. The ICA summary
% must not crash when the ranking's recommended candidate ID is [].
tests = functiontests(localfunctions);
end

function setupOnce(~)
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
end

function testEmptyRecommendationDoesNotCrash(tc)
r = fixture();
r.ranking.recommended = [];
report = pipecompare.simple.Presets.icaText(r);
verifyClass(tc, report, 'char');
verifyTrue(tc, isrow(report));
end

function testScalarRecommendationStillWorks(tc)
r = fixture();
r.ranking.recommended = 1;
report = pipecompare.simple.Presets.icaText(r);
verifyClass(tc, report, 'char');
verifyTrue(tc, isrow(report));
end

function testNoIcaEvidenceReturnsEmpty(tc)
r = struct('cands', struct('id', 1, 'steps', {{}}), ...
    'ranking', struct('recommended', []));
verifyEqual(tc, pipecompare.simple.Presets.icaText(r), '');
end

function r = fixture()
f = struct('type', 'icremove', 'icsTotal', 4, 'icsRemoved', 0, ...
    'icsBrain', 1, 'icsOther', 3, 'icaPoints', 100, ...
    'icsConsidered', 2, 'otherMedian', 0.4);
r = struct('cands', struct('id', 1, 'steps', {{f}}), ...
    'ranking', struct('recommended', []));
end
