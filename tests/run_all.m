function results = run_all(names)
%RUN_ALL Run the PipeCompare test suite and print a summary.
%   results = run_all()              all suites
%   results = run_all({'test_plan'}) selected suites
%   Needs EEGLAB on the path. Set PIPECOMPARE_REAL_SET to a real .set file
%   (never committed) to include the real-data test.
here = fileparts(mfilename('fullpath'));
addpath(here); addpath(fullfile(here, '..'));
if nargin < 1
    names = {'test_history', 'test_plan', 'test_panel_text', 'test_statistics', 'test_methodology', 'test_signal', 'test_engine', 'test_eeglab', 'test_simple', 'test_main', 'test_invariants'};
end
results = [];
for k = 1:numel(names)
    % Name each running test, so a stalled CI job identifies its last case.
    results = [results, runtests(fullfile(here, [names{k} '.m']), 'OutputDetail', 3)]; %#ok<AGROW>
end
disp(table(results));
skipped = sum([results.Incomplete] & ~[results.Failed]);
fprintf('\n%d passed, %d failed, %d incomplete (skipped) of %d\n', sum([results.Passed]), ...
    sum([results.Failed]), skipped, numel(results));
% Skipped tests are reported separately from failed tests.
fprintf('TOTAL=%d PASSED=%d FAILED=%d SKIPPED=%d\n', numel(results), sum([results.Passed]), ...
    sum([results.Failed]), skipped);
end
