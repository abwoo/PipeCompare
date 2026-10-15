classdef ChannelPicker < handle
    %CHANNELPICKER Choose electrodes from a list (the PipeCompare window's
    %   own list, instead of EEGLAB's).
    %
    %   [ok, chosen] = pipecompare.gui.ChannelPicker.ask(labels, selected)
    %   ok is false when cancelled; chosen are labels, in the order listed.

    properties
        Fig; List
        OK = false
    end

    methods (Static)
        function [ok, chosen] = ask(labels, selected)
            p = pipecompare.gui.ChannelPicker(labels, selected);
            c = onCleanup(@() delete(p));
            uiwait(p.Fig);
            ok = p.OK; chosen = {};
            if ok, chosen = labels(ismember(1:numel(labels), p.List.Value)); end
        end
    end

    methods
        function obj = ChannelPicker(labels, selected)
            obj.Fig = uifigure('Name', 'Electrodes', 'Position', pipecompare.utils.onScreen([300 200 300 460]), ...
                'WindowStyle', 'modal', 'CloseRequestFcn', @(~, ~) uiresume(obj.Fig));
            g = uigridlayout(obj.Fig, [3 1]); g.RowHeight = {'fit', '1x', 30};
            uilabel(g, 'Text', 'Choose the electrodes (Ctrl or Cmd + click for several); they are averaged.', 'WordWrap', 'on');
            obj.List = uilistbox(g, 'Items', labels, 'ItemsData', 1:numel(labels), 'Multiselect', 'on', ...
                'Value', find(ismember(lower(labels), lower(selected))));
            b = uigridlayout(g, [1 2]); b.Padding = [0 0 0 0];
            uibutton(b, 'Text', 'OK', 'FontWeight', 'bold', 'ButtonPushedFcn', @(~, ~) obj.done(true));
            uibutton(b, 'Text', 'Cancel', 'ButtonPushedFcn', @(~, ~) obj.done(false));
        end

        function done(obj, ok)
            obj.OK = ok;
            uiresume(obj.Fig);
        end

        function delete(obj)
            if ~isempty(obj.Fig) && isvalid(obj.Fig), delete(obj.Fig); end
        end
    end
end
