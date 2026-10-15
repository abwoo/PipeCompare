function m = helpMenu(fig)
%HELPMENU The Help menu of PipeCompare's windows: the GitHub repository,
%   the user guide, the latest version, reporting a problem, the version.
%   m = pipecompare.gui.helpMenu(fig)    fig: a uifigure
url = pipecompare.PipeCompare.Url;
m = uimenu(fig, 'Text', 'Help');
uimenu(m, 'Text', 'PipeCompare on GitHub', 'MenuSelectedFcn', @(~, ~) web(url, '-browser'));
uimenu(m, 'Text', 'User guide (README)', 'MenuSelectedFcn', @(~, ~) web([url '#readme'], '-browser'));
uimenu(m, 'Text', 'Download the latest version', 'MenuSelectedFcn', @(~, ~) web([url '/releases/latest'], '-browser'));
uimenu(m, 'Text', 'Report a problem', 'MenuSelectedFcn', @(~, ~) web([url '/issues'], '-browser'));
uimenu(m, 'Text', 'About PipeCompare', 'Separator', 'on', 'MenuSelectedFcn', @(~, ~) uialert(fig, ...
    sprintf('PipeCompare %s\n%s', pipecompare.PipeCompare.version(), url), 'About PipeCompare', 'Icon', 'info'));
end
