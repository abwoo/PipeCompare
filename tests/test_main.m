function tests = test_main
%TEST_MAIN The PipeCompare window (pipecompare.start): opening data,
%   channel locations, the Simple and Advanced tabs, a comparison with its
%   progress and results in the window, and export.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
tc.TestData.EEG = nqc_synth(struct('seconds', 120, 'nPerCond', 30, 'alphaUv', 10));
end

function w = openWindow(tc)
w = pipecompare.start();
tc.addTeardown(@() closeIf(w));
end

function testWindowShowsTheCurrentDataset(tc)
EEG = tc.TestData.EEG;
nqc_setBase(EEG);
w = openWindow(tc);
verifyEqual(tc, w.Nav.Items, {'1. Import', '2. Settings', '3. Run', '4. Results', '5. Export'});
verifyTrue(tc, any(contains(w.DataText.Value, 'nqc_synth')));
verifyTrue(tc, contains(w.LocText.Text, 'All 32 channels have positions'));
verifyEqual(tc, size(w.ChanTable.Data), [32 4]);               % the channels, editable
verifyNumElements(tc, findobj(w.DataView.Axes, 'Type', 'line'), 32);   % the data plotted, one trace per channel
verifyEqual(tc, w.DataList.Value, evalin('base', 'CURRENTSET'));   % the datasets, the current one selected
verifyClass(tc, w.Simple, 'pipecompare.gui.SimpleDialog');     % the simple mode, in the Simple tab
verifyEqual(tc, ancestor(w.Simple.Grid, 'figure'), w.Fig);
verifyEmpty(tc, findall(w.Fig, 'Type', 'uibutton', 'Text', 'Cancel'));   % in the window: no Cancel
gh = findall(w.Fig, 'Type', 'uimenu', 'Text', 'PipeCompare on GitHub');
verifyNumElements(tc, gh, 1);                                   % Help > PipeCompare on GitHub
verifyTrue(tc, isa(gh.MenuSelectedFcn, 'function_handle'));
verifyTrue(tc, pipecompare.start() == w);                      % one window: opening again shows it
w.showPage(5);
verifyTrue(tc, logical(w.Pages(5).Visible) && ~logical(w.Pages(1).Visible));
verifyTrue(tc, contains(w.ExportData.Text, 'nqc_synth'));
end

function testOpenFileAndStandardLocations(tc)
EEG = tc.TestData.EEG;
nqc_setBase(EEG);
L = EEG.chanlocs;                                               % a recording without positions
for f = {'X', 'Y', 'Z', 'theta', 'radius', 'sph_theta', 'sph_phi', 'sph_radius'}
    if isfield(L, f{1}), [L.(f{1})] = deal([]); end
end
E = EEG; E.chanlocs = L; E.setname = 'nolocs';
folder = tempname; mkdir(folder); tc.addTeardown(@() rmdir(folder, 's'));
evalc('pop_saveset(E, ''filename'', ''nolocs.set'', ''filepath'', folder);');
w = openWindow(tc);
n = evalin('base', 'numel(ALLEEG)');
evalc('w.openFile(fullfile(folder, ''nolocs.set''))');
verifyEqual(tc, evalin('base', 'numel(ALLEEG)'), n + 1);       % a new EEGLAB dataset, the current one
verifyEqual(tc, evalin('base', 'EEG.setname'), 'nolocs');
verifyTrue(tc, contains(evalin('base', 'EEG.history'), 'nolocs.set'));   % how it was opened
verifyTrue(tc, startsWith(w.LocText.Text, 'None of the 32 channels has a position'));
evalc('w.setLocations(''standard'')');
s = pipecompare.live.DataState.fromEEG(evalin('base', 'EEG'));
verifyGreaterThanOrEqual(tc, s.nLocated, 30);                   % EOG1 and EOG2 are not 10-05 labels
verifyTrue(tc, contains(evalin('base', 'EEG.history'), 'pop_chanedit'));
verifyFalse(tc, startsWith(w.LocText.Text, 'None'));
verifyEqual(tc, evalin('base', 'numel(ALLEEG)'), n + 1);       % stored in place of the current dataset
d = w.ChanTable.Data; d{1, 3} = 'EOG'; d{2, 2} = 'TEST2'; w.ChanTable.Data = d;
evalc('w.applyChannelEdits()');                                 % labels and types corrected in the table
L = evalin('base', 'EEG.chanlocs');
verifyEqual(tc, L(1).type, 'EOG'); verifyEqual(tc, L(2).labels, 'TEST2');
verifyTrue(tc, contains(evalin('base', 'EEG.history'), 'changefield'));
verifyNumElements(tc, w.DataList.Items, n + 1);                 % the datasets, under the steps
evalc('w.selectDataset(n)');                                    % clicking one makes it current
verifyEqual(tc, evalin('base', 'CURRENTSET'), n);
verifyEqual(tc, evalin('base', 'EEG.setname'), 'nqc_synth');
verifyTrue(tc, any(contains(w.DataText.Value, 'nqc_synth')));
end

function testSimpleRunShowsProgressAndResultsInTheWindow(tc)
EEG = tc.TestData.EEG;
nqc_setBase(EEG);
w = openWindow(tc);
nFig = numel(findall(groot, 'Type', 'figure'));
d = w.Simple;
d.EventList.Value = {'11', '31'}; d.MeasureDrop.Value = 'P3'; d.measureChanged(); d.usePreset('filters');
verifyEqual(tc, char(d.RunButton.Enable), 'on');
evalc('d.run()');
verifyEqual(tc, numel(findall(groot, 'Type', 'figure')), nFig);   % no progress or results window of their own
verifyClass(tc, w.Results, 'pipecompare.gui.SimpleResults');
verifyEqual(tc, w.Nav.Value, 4);                                % the Results page
verifyEqual(tc, char(w.Nav.Enable), 'on');
verifyEqual(tc, numel(w.Result.cands), 12);
verifyTrue(tc, startsWith(strjoin(w.Results.Headline.Value, ' '), sprintf('Use pipeline %d', w.Result.ranking.recommended)));
verifyEmpty(tc, w.RunBox.Children);                             % the progress is gone
evalc('w.plotPipeline()');                                      % the recommended pipeline's ERP
verifyNumElements(tc, findobj(w.ResultView.Axes, 'Type', 'line'), 2);   % one line per condition
verifyEqual(tc, w.ResultView.Names, {sprintf('pipeline %d', w.Result.ranking.recommended)});
f = [tempname '.m']; tc.addTeardown(@() deleteIf(f));
evalc('w.saveScript(f)');
verifyTrue(tc, isfile(f));
nMain = numel(findobj(groot, 'Tag', 'EEGLAB'));
evalc('w.Results.adopt()');                                     % Use this pipeline
verifyEqual(tc, w.Nav.Value, 5);                                % what next, on the Export page
verifyTrue(tc, startsWith(w.ExportNote.Text, sprintf('Pipeline %d is now the current dataset', w.Result.ranking.recommended)));
verifyTrue(tc, contains(evalin('base', 'EEG.setname'), 'PipeCompare#'));
verifyTrue(tc, contains(w.ExportData.Text, 'PipeCompare#'));
verifyNumElements(tc, findobj(groot, 'Tag', 'EEGLAB'), nMain);  % EEGLAB's main window was not opened
folder = tempname; mkdir(folder); tc.addTeardown(@() rmdir(folder, 's'));
evalc('w.saveDataset(fullfile(folder, ''clean.set''))');
verifyTrue(tc, isfile(fullfile(folder, 'clean.set')));
verifyEqual(tc, evalin('base', 'EEG.filename'), 'clean.set');
nqc_setBase(EEG);
end

function testAdvancedTabTakesTheSimpleChoices(tc)
EEG = tc.TestData.EEG;
nqc_setBase(EEG);
w = openWindow(tc);
nFig = numel(findall(groot, 'Type', 'figure'));
w.Tabs.SelectedTab = w.AdvancedTab; w.tabChanged();
verifyClass(tc, w.Advanced, 'pipecompare.gui.Panel');           % the panel, in the Advanced tab
verifyEqual(tc, w.Advanced.Root.Parent, w.AdvancedTab);
first = w.Advanced;
open = ~isempty(findobj(groot, 'Tag', 'EEGLAB'));
verifyEqual(tc, char(first.MenuStepButton.Enable), pipecompare.utils.ternary(open, 'on', 'off'));   % EEGLAB's menus: only with its window
d = w.Simple;
d.EventList.Value = {'11', '31'}; d.MeasureDrop.Value = 'P3'; d.measureChanged();
evalc('w.openAdvanced(d)');
verifyFalse(tc, isvalid(first));                                % a new panel with the choices
verifyTrue(tc, w.Tabs.SelectedTab == w.AdvancedTab);
verifyTrue(tc, ismember('ica', {w.Advanced.Plan.Slots.id}));    % the Standard recipe
verifyEqual(tc, {w.Advanced.contract().conditions.name}, {'11', '31'});
verifyTrue(tc, isvalid(d) && isvalid(d.Grid));                  % the Simple tab stays
verifyEqual(tc, numel(findall(groot, 'Type', 'figure')), nFig);
p = w.Advanced; t = p.Timer;
delete(w);
verifyFalse(tc, isvalid(p) || isvalid(t));                      % closing the window stops the panel
end

function closeIf(w)
if isvalid(w), delete(w); end
end

function deleteIf(f)
if isfile(f), delete(f); end
end
