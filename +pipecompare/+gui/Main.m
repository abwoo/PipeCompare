classdef Main < handle
    %MAIN The PipeCompare window: every step in one window, with or without
    %   EEGLAB's main window.
    %
    %   pipecompare.start                (or EEGLAB > Tools > PipeCompare > PipeCompare window)
    %   w = pipecompare.gui.Main.open()  the window (the one already open, if any)
    %
    %   On the left, the steps in order, each showing its page on the right,
    %   and under them the datasets open (clicking one makes it current):
    %   1. Import    open a data file (.set; EDF, BDF and GDF through
    %                EEGLAB's BIOSIG import; other formats through its
    %                FILE-IO import), or use the dataset already current in
    %                EEGLAB; add channel locations (the standard 10-05
    %                positions looked up by channel label, or a locations
    %                file); correct channel labels and types in a table; a
    %                plot of the data (pipecompare.gui.DataView).
    %   2. Settings  two tabs: Simple (the simple mode's choices,
    %                pipecompare.gui.SimpleDialog) and Advanced (the full
    %                panel, pipecompare.gui.Panel); Advanced... in Simple
    %                carries the choices over to Advanced.
    %   3. Run       the progress and Stop of the comparison running now.
    %   4. Results   the recommendation in plain words, the pipelines, Use
    %                this pipeline and Save script
    %                (pipecompare.gui.SimpleResults), also after a search of
    %                the Advanced tab (whose own table has the details); a
    %                plot of the selected pipeline's ERP (or spectrum), to
    %                compare pipelines.
    %   5. Export    save the current dataset (after Use this pipeline: the
    %                cleaned data) as a .set file, and the pipeline as a
    %                script for other recordings.
    %   The data are EEGLAB's: the window works on the dataset current in
    %   EEGLAB (base EEG, ALLEEG, CURRENTSET) and stores what it opens or
    %   makes there, so EEGLAB's menus and the command line see the same
    %   data; it re-reads the dataset when a page is shown. EEGLAB must be
    %   on the MATLAB path; it is started without its main window if it has
    %   not been started yet.

    properties
        Fig
        Nav                     % the list of steps
        DataList                % the datasets open in EEGLAB, under the steps
        Pages                   % one panel per step
        DataText; LocText; LocButtons
        ChanTable               % channel labels and types, editable
        DataView                % plot of the current dataset (Import)
        ResultView              % plot of pipelines (Results)
        Tabs; SimpleTab; AdvancedTab; SimpleNote; SimpleBox
        Simple = []             % pipecompare.gui.SimpleDialog in the Simple tab
        Advanced = []           % pipecompare.gui.Panel in the Advanced tab
        RunBox; RunNote
        ResultsBox; ResultsNote
        Results = []            % pipecompare.gui.SimpleResults on the Results page
        Result = []             % the last comparison
        ExportData; ExportNote
        Quick = '-'             % Session.quickPrint of the dataset shown
        Busy = false            % a comparison is running
    end

    properties (Constant)
        PageNames = {'1. Import', '2. Settings', '3. Run', '4. Results', '5. Export'}
    end

    methods (Static)
        function w = open()
            f = findall(groot, 'Type', 'figure', 'Tag', 'PipeCompareMain');
            if ~isempty(f) && isappdata(f(1), 'PipeCompareMain')
                w = getappdata(f(1), 'PipeCompareMain');
                figure(f(1));
                return;
            end
            w = pipecompare.gui.Main();
        end

        function startEeglab()
            % EEGLAB's functions on the path, without its main window
            if exist('pop_loadset', 'file') == 2 && exist('eeg_store', 'file') == 2, return; end
            assert(exist('eeglab', 'file') == 2, 'PipeCompare:Dependency', ['PipeCompare needs EEGLAB: add the EEGLAB ', ...
                'folder to the MATLAB path (addpath(''/path/to/eeglab'')), then start PipeCompare again.']);
            evalin('base', 'eeglab nogui;');
        end

        function [EEG, com] = readFile(file)
            % An EEG recording read with EEGLAB's import for its format. A
            % format plugin that is missing is offered by EEGLAB itself; if
            % it is not installed, EEG is empty.
            [p, n, e] = fileparts(file);
            switch lower(e)
                case '.set'
                    [EEG, com] = pop_loadset('filename', [n e], 'filepath', p);
                case {'.edf', '.bdf', '.gdf'}
                    [EEG, com] = pop_biosig(file);
                otherwise
                    [EEG, com] = pop_fileio(file);
            end
        end

        function f = standardLocations()
            % The standard 10-05 positions shipped with EEGLAB (DIPFIT)
            f = '';
            if exist('dipfitdefs', 'file'), f = fullfile(fileparts(which('dipfitdefs')), 'standard_BEM', 'elec', 'standard_1005.elc'); end
            if ~isfile(f), f = which('Standard-10-5-Cap385.sfp'); end
            assert(~isempty(f), 'PipeCompare:Locations', ['EEGLAB''s standard channel locations were not found (DIPFIT ', ...
                'plugin); load the locations file of your cap instead.']);
        end

        function storeNew(EEG, com)
            % EEG as a new EEGLAB dataset, the current one (as EEGLAB's
            % import menus store it)
            if ~isempty(com), EEG = eegh(com, EEG); end
            assignin('base', 'PIPECOMPARE_NEW__', EEG);
            evalin('base', ['if ~exist(''ALLEEG'', ''var''), ALLEEG = []; end; if ~exist(''CURRENTSET'', ''var''), CURRENTSET = 0; end; ', ...
                '[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, PIPECOMPARE_NEW__, CURRENTSET, ''gui'', ''off''); clear PIPECOMPARE_NEW__;']);
            pipecompare.gui.Main.told(com);
        end

        function storeCurrent(EEG, com, inHistory)
            % EEG in place of the current EEGLAB dataset; com in its history
            % (not for saving, which changes nothing in the data)
            if nargin < 3, inHistory = true; end
            if inHistory && ~isempty(com), EEG = eegh(com, EEG); end
            assignin('base', 'PIPECOMPARE_NEW__', EEG);
            evalin('base', ['EEG = PIPECOMPARE_NEW__; clear PIPECOMPARE_NEW__; if exist(''ALLEEG'', ''var'') && ', ...
                'exist(''CURRENTSET'', ''var'') && isscalar(CURRENTSET) && CURRENTSET > 0, ', ...
                '[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET); end']);
            pipecompare.gui.Main.told(com);
        end

        function told(com)
            % EEGLAB's session history (ALLCOM) and its main window, if open
            if ~isempty(com)
                try, eegh(com); catch, end
            end
            if ~isempty(findobj(groot, 'Tag', 'EEGLAB')), evalin('base', 'eeglab redraw;'); end
        end

        function ended(w)
            % a search of the window has ended (the window may be closed)
            if ~isempty(w) && isvalid(w), w.runEnded(); end
        end
    end

    methods
        function obj = Main()
            pipecompare.gui.Main.startEeglab();
            obj.Fig = uifigure('Name', sprintf('PipeCompare %s', pipecompare.PipeCompare.version()), 'Tag', 'PipeCompareMain', ...
                'Position', pipecompare.utils.onScreen([40 40 1500 900]), 'CloseRequestFcn', @(~, ~) delete(obj));
            setappdata(obj.Fig, 'PipeCompareMain', obj);
            pipecompare.gui.helpMenu(obj.Fig);   % GitHub, user guide, latest version
            g = uigridlayout(obj.Fig, [1 2]); g.ColumnWidth = {170, '1x'};
            left = uigridlayout(g, [4 1]); left.Padding = [0 0 0 0]; left.RowHeight = {30, 130, 22, '1x'};
            uilabel(left, 'Text', 'PipeCompare', 'FontWeight', 'bold', 'FontSize', 16);
            obj.Nav = uilistbox(left, 'Items', obj.PageNames, 'ItemsData', 1:numel(obj.PageNames), 'Value', 1, ...
                'ValueChangedFcn', @(h, ~) obj.showPage(h.Value));
            uilabel(left, 'Text', 'Datasets', 'FontWeight', 'bold');
            obj.DataList = uilistbox(left, 'Items', {}, 'Tooltip', 'The datasets open (in EEGLAB too); click one to use it', ...
                'ValueChangedFcn', @(h, ~) obj.selectDataset(h.Value));
            right = uigridlayout(g, [1 1]); right.Padding = [0 0 0 0];
            obj.Pages = gobjects(1, numel(obj.PageNames));
            for k = 1:numel(obj.PageNames)   % all in one place; the one shown is visible
                obj.Pages(k) = uipanel(right, 'Title', obj.PageNames{k}, 'FontWeight', 'bold', 'Visible', k == 1);
                obj.Pages(k).Layout.Row = 1; obj.Pages(k).Layout.Column = 1;
            end
            obj.buildImport(); obj.buildSettings(); obj.buildRun(); obj.buildResults(); obj.buildExport();
            obj.Fig.AutoResizeChildren = 'off';   % the grids resize themselves; this lets SizeChangedFcn run
            obj.Fig.SizeChangedFcn = @(~, ~) obj.resized();
            obj.refreshData();
        end

        function delete(obj)
            if ~isempty(obj.Advanced) && isvalid(obj.Advanced), delete(obj.Advanced); end
            if ~isempty(obj.Simple) && isvalid(obj.Simple), delete(obj.Simple); end
            if ~isempty(obj.Fig) && isvalid(obj.Fig), delete(obj.Fig); end
        end

        % ------------------------------------------------------------ pages
        function buildImport(obj)
            grey = [0.3 0.3 0.3];
            two = uigridlayout(obj.Pages(1), [1 2]); two.ColumnWidth = {'1x', '1.3x'};
            g = uigridlayout(two, [9 1]); g.Padding = [0 0 0 0]; g.RowHeight = {'fit', 30, 120, 22, 'fit', 30, '1x', 30, 30};
            obj.DataView = pipecompare.gui.DataView(two);
            uilabel(g, 'WordWrap', 'on', 'Text', ['PipeCompare works on one EEG recording: open its file here, or load ', ...
                'it in EEGLAB (this window uses the dataset that is current in EEGLAB). Start from the raw recording: ', ...
                'PipeCompare compares the filters, reference, bad channels, ICA and epoch rejection itself.']);
            b = uigridlayout(g, [1 2]); b.Padding = [0 0 0 0]; b.ColumnWidth = {170, '1x'};
            uibutton(b, 'Text', 'Open data file...', 'FontWeight', 'bold', 'ButtonPushedFcn', @(~, ~) obj.openFile());
            uilabel(b, 'FontColor', grey, 'Text', ['.set (EEGLAB); EDF, BDF, GDF (EEGLAB''s BIOSIG import); other ', ...
                'formats through EEGLAB''s FILE-IO import. EEGLAB offers to install an import it lacks.']);
            obj.DataText = uitextarea(g, 'Editable', 'off', 'Value', {''});
            uilabel(g, 'Text', 'Channel locations', 'FontWeight', 'bold');
            obj.LocText = uilabel(g, 'Text', '', 'WordWrap', 'on');
            l = uigridlayout(g, [1 3]); l.Padding = [0 0 0 0]; l.ColumnWidth = {200, 120, '1x'};
            obj.LocButtons = [ ...
                uibutton(l, 'Text', 'Look up standard positions', 'Tooltip', ['The positions of the standard 10-05 system, ', ...
                    'found by channel label (EEGLAB''s channel lookup)'], 'ButtonPushedFcn', @(~, ~) obj.setLocations('standard')), ...
                uibutton(l, 'Text', 'From a file...', 'Tooltip', 'The locations file of your cap (.ced, .elc, .elp, .sfp, .loc, ...)', ...
                    'ButtonPushedFcn', @(~, ~) obj.setLocations('file'))];
            uilabel(l, 'Text', '');
            obj.ChanTable = uitable(g, 'ColumnName', {'#', 'label', 'type', 'position'}, 'RowName', {}, ...
                'ColumnEditable', [false true true false], 'ColumnWidth', {35, 'auto', 'auto', 70}, 'Tooltip', ...
                ['Correct a label (then look up the positions again) or a type (EOG, ECG, EMG, REF: not used as scalp ', ...
                 'channels); Apply stores the changes']);
            n = uigridlayout(g, [1 3]); n.Padding = [0 0 0 0]; n.ColumnWidth = {200, '1x', 170};
            obj.LocButtons(end+1) = uibutton(n, 'Text', 'Apply label and type changes', 'ButtonPushedFcn', @(~, ~) obj.applyChannelEdits());
            uilabel(n, 'Text', '');
            uilabel(n, 'Text', '');
            n = uigridlayout(g, [1 2]); n.Padding = [0 0 0 0]; n.ColumnWidth = {'1x', 170};
            uilabel(n, 'Text', '');
            uibutton(n, 'Text', 'Next: 2. Settings', 'ButtonPushedFcn', @(~, ~) obj.showPage(2));
        end

        function buildSettings(obj)
            g = uigridlayout(obj.Pages(2), [1 1]); g.Padding = [0 0 0 0];
            obj.Tabs = uitabgroup(g, 'SelectionChangedFcn', @(~, ~) obj.tabChanged());
            obj.SimpleTab = uitab(obj.Tabs, 'Title', 'Simple');
            obj.AdvancedTab = uitab(obj.Tabs, 'Title', 'Advanced');
            s = uigridlayout(obj.SimpleTab, [1 1]); s.Padding = [0 0 0 0];
            obj.SimpleNote = uilabel(s, 'Text', 'Open a data file first (1. Import).', 'HorizontalAlignment', 'center');
            obj.SimpleBox = uipanel(s, 'BorderType', 'none');
            obj.SimpleNote.Layout.Row = 1; obj.SimpleNote.Layout.Column = 1;
            obj.SimpleBox.Layout.Row = 1; obj.SimpleBox.Layout.Column = 1;
        end

        function buildRun(obj)
            g = uigridlayout(obj.Pages(3), [2 1]); g.RowHeight = {190, '1x'};
            obj.RunBox = uipanel(g, 'BorderType', 'none');
            obj.RunNote = uilabel(g, 'WordWrap', 'on', 'VerticalAlignment', 'top', 'Text', ['Nothing is running. Choose ', ...
                'what to compare in 2. Settings, then press Run (Run search in the Advanced tab).']);
            obj.RunBox.Layout.Row = 1; obj.RunNote.Layout.Row = 1;
        end

        function buildResults(obj)
            g = uigridlayout(obj.Pages(4), [1 2]); g.Padding = [0 0 0 0]; g.ColumnWidth = {'1.15x', '1x'};
            obj.ResultsNote = uilabel(g, 'Text', 'No comparison yet (2. Settings, then Run).', 'HorizontalAlignment', 'center');
            obj.ResultsBox = uipanel(g, 'BorderType', 'none');
            obj.ResultsNote.Layout.Row = 1; obj.ResultsNote.Layout.Column = 1;
            obj.ResultsBox.Layout.Row = 1; obj.ResultsBox.Layout.Column = 1;
            p = uigridlayout(g, [2 1]); p.RowHeight = {30, '1x'}; p.Layout.Row = 1; p.Layout.Column = 2;
            b = uigridlayout(p, [1 3]); b.Padding = [0 0 0 0]; b.ColumnWidth = {170, 100, '1x'};
            uibutton(b, 'Text', 'Plot selected pipeline', 'Tooltip', ['Runs the selected pipeline (else the ', ...
                'recommended one) again and plots its ERP or spectrum next to those plotted before'], ...
                'ButtonPushedFcn', @(~, ~) obj.plotPipeline());
            uibutton(b, 'Text', 'Clear plot', 'ButtonPushedFcn', @(~, ~) obj.ResultView.clearPlot());
            uilabel(b, 'Text', '');
            obj.ResultView = pipecompare.gui.DataView(p);
        end

        function buildExport(obj)
            g = uigridlayout(obj.Pages(5), [8 1]); g.RowHeight = {22, 'fit', 30, 'fit', 22, 'fit', 30, '1x'};
            uilabel(g, 'Text', 'The data', 'FontWeight', 'bold');
            obj.ExportData = uilabel(g, 'Text', '', 'WordWrap', 'on');
            b = uigridlayout(g, [1 2]); b.Padding = [0 0 0 0]; b.ColumnWidth = {200, '1x'};
            uibutton(b, 'Text', 'Save dataset (.set)...', 'ButtonPushedFcn', @(~, ~) obj.saveDataset());
            uilabel(b, 'Text', '');
            obj.ExportNote = uilabel(g, 'Text', '', 'WordWrap', 'on', 'FontColor', [0 0.4 0.1]);
            uilabel(g, 'Text', 'The pipeline as a script', 'FontWeight', 'bold');
            uilabel(g, 'WordWrap', 'on', 'Text', ['A MATLAB function that runs the same steps on other recordings ', ...
                '(the recommended pipeline, or the one selected on 4. Results); bad channels, ICA components and ', ...
                'rejected epochs are decided from each recording''s own data.']);
            s = uigridlayout(g, [1 2]); s.Padding = [0 0 0 0]; s.ColumnWidth = {200, '1x'};
            uibutton(s, 'Text', 'Save script...', 'ButtonPushedFcn', @(~, ~) obj.saveScript());
            uilabel(s, 'Text', '');
        end

        function showPage(obj, k)
            obj.refreshIfChanged();
            for i = 1:numel(obj.Pages), obj.Pages(i).Visible = i == k; end
            obj.Nav.Value = k;
            if k == 2, obj.tabChanged(); end
        end

        function tabChanged(obj)
            % the panel is made when the Advanced tab is first shown
            if obj.Tabs.SelectedTab == obj.AdvancedTab && (isempty(obj.Advanced) || ~isvalid(obj.Advanced))
                obj.newPanel();
            end
        end

        function openAdvanced(obj, d)
            % Advanced... in the Simple tab: a new panel with its choices
            app = d.advanced(@() obj.newPanel());
            if ~isempty(app), obj.Tabs.SelectedTab = obj.AdvancedTab; end
        end

        function app = newPanel(obj)
            % a new panel in the Advanced tab (in place of the one there)
            if ~isempty(obj.Advanced) && isvalid(obj.Advanced), delete(obj.Advanced); end
            obj.Advanced = pipecompare.gui.Panel(obj.AdvancedTab, obj);
            app = obj.Advanced;
        end

        function resized(obj)
            if ~isempty(obj.Advanced) && isvalid(obj.Advanced), obj.Advanced.fit(); end
        end

        % ------------------------------------------------------------- data
        function refreshIfChanged(obj)
            try, [~, q] = pipecompare.live.Session.peek(); catch, q = ''; end
            if ~strcmp(q, obj.Quick), obj.refreshData(); end
        end

        function refreshData(obj)
            % What the window shows of the current EEGLAB dataset; the
            % Simple tab is made again for it.
            [EEG, q] = pipecompare.live.Session.peek(); obj.Quick = q;
            if numel(EEG) ~= 1 || isempty(obj.Quick), EEG = []; end
            if ~isempty(obj.Simple) && isvalid(obj.Simple), delete(obj.Simple); end
            obj.Simple = [];
            [on, off] = deal('on', 'off');
            obj.showDatasets();
            obj.ChanTable.Data = obj.channelRows(EEG);
            obj.DataView.showData(EEG);
            if isempty(EEG)
                obj.DataText.Value = {'No dataset yet: open a data file.'};
                obj.LocText.Text = '';
                obj.ExportData.Text = 'No dataset yet.';
                set(obj.LocButtons, 'Enable', off);
                obj.SimpleNote.Text = 'Open a data file first (1. Import).'; obj.SimpleNote.Visible = on;
                return;
            end
            s = pipecompare.live.DataState.fromEEG(EEG);
            obj.DataText.Value = obj.dataText(EEG, s);
            obj.LocText.Text = obj.locationsText(s);
            set(obj.LocButtons, 'Enable', on);
            obj.ExportData.Text = sprintf('Current dataset: %s (%s). Save it to keep it: it is in memory only until saved.', ...
                s.setname, obj.shapeText(s));
            try
                obj.Simple = pipecompare.gui.SimpleDialog(EEG, obj.SimpleBox, obj);
                obj.SimpleNote.Visible = off;
            catch ME
                obj.SimpleNote.Text = ME.message; obj.SimpleNote.Visible = on;
            end
        end

        function showDatasets(obj)
            % the datasets in EEGLAB's ALLEEG, the current one selected
            names = {}; idx = [];
            try
                if evalin('base', 'exist(''ALLEEG'', ''var'') && exist(''CURRENTSET'', ''var'')')
                    A = evalin('base', 'ALLEEG'); cur = evalin('base', 'CURRENTSET');
                    idx = find(arrayfun(@(e) isfield(e, 'data') && ~isempty(e.data), A));
                    names = arrayfun(@(k) sprintf('%d: %s', k, char(string(A(k).setname))), idx, 'UniformOutput', false);
                end
            catch
            end
            obj.DataList.Items = names; obj.DataList.ItemsData = idx;
            if ~isempty(idx) && isscalar(cur) && ismember(cur, idx), obj.DataList.Value = cur; end
        end

        function selectDataset(obj, k)
            % dataset k of ALLEEG becomes the current one (as in EEGLAB's
            % Datasets menu)
            if obj.Busy, obj.showDatasets(); return; end
            try
                evalin('base', sprintf('[EEG, ALLEEG, CURRENTSET] = eeg_retrieve(ALLEEG, %d);', k));
                pipecompare.gui.Main.told('');
            catch ME
                uialert(obj.Fig, ME.message, 'PipeCompare');
            end
            obj.refreshData();
        end

        function rows = channelRows(~, EEG)
            % number, label, type and whether it has a position, per channel
            rows = cell(0, 4);
            if isempty(EEG) || ~isfield(EEG, 'chanlocs') || numel(EEG.chanlocs) ~= EEG.nbchan, return; end
            L = EEG.chanlocs;
            rows = cell(numel(L), 4);
            for k = 1:numel(L)
                ty = ''; if isfield(L, 'type') && ~isempty(L(k).type), ty = char(string(L(k).type)); end
                has = isfield(L, 'X') && ~isempty(L(k).X) && ~isnan(double(L(k).X(1)));
                rows(k, :) = {k, char(string(L(k).labels)), ty, pipecompare.utils.ternary(has, 'yes', 'no')};
            end
        end

        function applyChannelEdits(obj)
            % The labels and types changed in the table, as one EEGLAB
            % channel edit of the current dataset (in its history)
            try
                EEG = pipecompare.live.Session.current();
                assert(~isempty(EEG), 'PipeCompare:NoDataset', 'Open a data file first.');
                was = obj.channelRows(EEG); edited = obj.ChanTable.Data;
                assert(isequal(size(was), size(edited)), 'PipeCompare:Channels', 'The dataset changed; the table is shown again.');
                args = {};
                for k = 1:size(edited, 1)
                    for f = {2, 3; 'labels', 'type'}
                        v = strtrim(char(string(edited{k, f{1}})));
                        if ~strcmp(v, was{k, f{1}}), args = [args {'changefield', {k, f{2}, v}}]; end %#ok<AGROW>
                    end
                end
                if isempty(args), return; end
                [EEG, ~, ~, com] = pop_chanedit(EEG, args{:});
                pipecompare.gui.Main.storeCurrent(EEG, com);
            catch ME
                uialert(obj.Fig, ME.message, 'PipeCompare');
            end
            obj.refreshData();
        end

        function v = dataText(~, EEG, s)
            v = {sprintf('Dataset: %s', s.setname)};
            if isfield(EEG, 'filename') && ~isempty(EEG.filename)
                v{end+1} = sprintf('File: %s', fullfile(char(EEG.filepath), char(EEG.filename)));
            end
            v{end+1} = sprintf('%d channels, %g Hz, %s', s.nbchan, s.srate, pipecompare.gui.Main.shapeText(s));
            types = s.eventTypes(~strcmpi(s.eventTypes, 'boundary'));
            v{end+1} = sprintf('%d events of %d types%s', s.nEvents, numel(types), ...
                pipecompare.utils.ternary(isempty(types), '', [': ' strjoin(types(1:min(12, end)), ', ') ...
                pipecompare.utils.ternary(numel(types) > 12, ', ...', '')]));
            v = v(:);
        end

        function t = locationsText(~, s)
            if s.hasLocations
                t = sprintf('All %d channels have positions.', s.nbchan);
            elseif s.nLocated == 0
                t = sprintf(['None of the %d channels has a position. Interpolating bad channels and ICLabel need them: ', ...
                    'look up the standard positions by channel label, or load the locations file of your cap.'], s.nbchan);
            else
                t = sprintf('%d of %d channels have positions (none for %s).', s.nLocated, s.nbchan, strjoin(s.unlocated, ', '));
            end
        end

        function openFile(obj, file)
            % A recording as a new EEGLAB dataset (file: without the dialog)
            if nargin < 2
                [f, p] = uigetfile({'*.set;*.edf;*.bdf;*.gdf;*.vhdr;*.cnt;*.mff;*.eeg', 'EEG data'; '*.*', 'All files'}, ...
                    'Open an EEG recording');
                if isequal(f, 0), return; end
                file = fullfile(p, f);
            end
            try
                [EEG, com] = pipecompare.gui.Main.readFile(file);
                if ~isempty(EEG), pipecompare.gui.Main.storeNew(EEG, com); end   % (empty: an import not installed)
            catch ME
                uialert(obj.Fig, ME.message, 'PipeCompare');
            end
            obj.refreshData();
        end

        function setLocations(obj, how, file)
            % Channel locations of the current dataset: 'standard' (looked
            % up by label) or 'file' (file: without the dialog)
            try
                EEG = pipecompare.live.Session.current();
                assert(~isempty(EEG), 'PipeCompare:NoDataset', 'Open a data file first.');
                switch how
                    case 'standard'
                        [EEG, ~, ~, com] = pop_chanedit(EEG, 'lookup', pipecompare.gui.Main.standardLocations());
                    case 'file'
                        if nargin < 3
                            [f, p] = uigetfile({'*.ced;*.elc;*.elp;*.sfp;*.loc;*.locs;*.xyz;*.txt', 'Channel locations'; ...
                                '*.*', 'All files'}, 'Channel locations of the cap');
                            if isequal(f, 0), return; end
                            file = fullfile(p, f);
                        end
                        [EEG, ~, ~, com] = pop_chanedit(EEG, 'load', {file, 'filetype', 'autodetect'});
                end
                if isempty(com), return; end   % cancelled
                pipecompare.gui.Main.storeCurrent(EEG, com);
            catch ME
                uialert(obj.Fig, ME.message, 'PipeCompare');
            end
            obj.refreshData();
        end

        % -------------------------------------------------------------- run
        function runSimple(obj, o)
            % The Simple tab's Run: the comparison, with its progress on the
            % Run page and its result on the Results page.
            if isempty(o) || obj.Busy, return; end
            args = reshape([fieldnames(o)'; struct2cell(o)'], 1, []);
            where = obj.runStarted();
            ended = onCleanup(@() pipecompare.gui.Main.ended(obj)); %#ok<NASGU>   % also on an error
            try
                EEG = pipecompare.live.Session.current();
                [~, com, r] = pop_pipecompare(EEG, args{:}, 'progress', where);
                if ~isvalid(obj), return; end   % the window was closed meanwhile
                if ~isempty(com), try, eegh(com); catch, end, end   % EEGLAB's session history, as from its menu
                obj.showResult(r);
            catch ME
                if isvalid(obj) && isvalid(obj.Fig), uialert(obj.Fig, ME.message, 'PipeCompare'); end
            end
        end

        function where = runStarted(obj)
            % the Run page, where the progress goes; the steps wait
            obj.Busy = true;
            obj.showPage(3);
            obj.Nav.Enable = 'off'; obj.RunNote.Visible = 'off';
            where = obj.RunBox;
        end

        function runEnded(obj)
            obj.Busy = false;
            obj.Nav.Enable = 'on'; obj.RunNote.Visible = 'on';
            if obj.Nav.Value == 3, obj.showPage(2); end   % no result (an error): back to the settings
        end

        function showResult(obj, r)
            obj.Result = r;
            if ~isempty(obj.Results) && isvalid(obj.Results.Root), delete(obj.Results.Root); end
            obj.Results = pipecompare.gui.SimpleResults(r, obj.ResultsBox, obj);
            obj.ResultsNote.Visible = 'off';
            obj.ResultView.clearPlot();   % the lines of another comparison
            obj.RunNote.Text = 'The last comparison is on 4. Results.';
            obj.showPage(4);
        end

        function plotPipeline(obj)
            % The selected pipeline (else the recommended one), run again
            % from the start of the comparison (ICA reused), on the plot
            if isempty(obj.Results) || ~isvalid(obj.Results.Root)
                uialert(obj.Fig, 'Run a comparison first (2. Settings, then Run).', 'PipeCompare'); return;
            end
            dlg = [];
            try
                r = obj.Result;
                k = pipecompare.run.Executor.pickCandidate(r, obj.Results.chosen());
                dlg = uiprogressdlg(obj.Fig, 'Title', 'Plot', 'Indeterminate', 'on', 'Message', ...
                    sprintf('Running pipeline %d again to plot it (ICA comes from the comparison)...', k));
                [~, E] = evalc('pipecompare.run.Executor.replay(r, k)');
                obj.ResultView.addPipeline(E, r.contract, sprintf('pipeline %d', k));
                close(dlg);
            catch ME
                if ~isempty(dlg) && isvalid(dlg), close(dlg); end
                uialert(obj.Fig, ME.message, 'PipeCompare');
            end
        end

        function adopted(obj, k, result)
            % Use this pipeline (here or in the Advanced tab) made it the
            % current dataset: what next, on the Export page
            obj.refreshData();
            obj.ExportNote.Text = sprintf(['Pipeline %d is now the current dataset. Its preprocessing is complete: ', ...
                '%s. Filtering, re-referencing or rejecting epochs again would change what was compared.'], k, ...
                pipecompare.simple.Presets.nextText(result, k));
            obj.showPage(5);
        end

        % ----------------------------------------------------------- export
        function saveDataset(obj, file)
            % The current dataset as a .set file (file: without the dialog)
            try
                EEG = pipecompare.live.Session.current();
                assert(~isempty(EEG), 'PipeCompare:NoDataset', 'There is no dataset to save.');
                if nargin < 2
                    [f, p] = uiputfile('*.set', 'Save the dataset', [regexprep(char(EEG.setname), '[^\w-]+', '_') '.set']);
                    if isequal(f, 0), return; end
                    file = fullfile(p, f);
                end
                [p, n, e] = fileparts(file);
                [EEG, com] = pop_saveset(EEG, 'filename', [n e], 'filepath', p);
                pipecompare.gui.Main.storeCurrent(EEG, com, false);
                obj.ExportNote.Text = sprintf('Saved %s.', file);
            catch ME
                uialert(obj.Fig, ME.message, 'PipeCompare');
            end
            obj.refreshIfChanged();
        end

        function saveScript(obj, varargin)
            % As Save script... on the Results page (file: without the dialog)
            if isempty(obj.Results) || ~isvalid(obj.Results.Root)
                uialert(obj.Fig, 'Run a comparison first (2. Settings, then Run).', 'PipeCompare'); return;
            end
            obj.Results.saveScript(varargin{:});
        end
    end

    methods (Static)
        function t = shapeText(s)
            if s.isEpoched, t = sprintf('%d epochs of %g to %g s', s.trials, s.xmin, s.xmax);
            else, t = sprintf('continuous, %.1f min', s.pnts / s.srate / 60); end
        end
    end
end
