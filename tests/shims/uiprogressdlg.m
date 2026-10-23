function dlg = uiprogressdlg(parent, varargin)
%UIPROGRESSDLG CI stand-in for a native dialog that blocks on legacy Xvfb.
%   Only the R2022b/R2023a CI commands add tests/shims to the path.
%   Callers use the returned handle only with isvalid/close. Scientific
%   computation, candidate selection, replay and export still run in full.
assert(isvalid(parent), 'PipeCompare:TestProgress', 'Invalid parent figure.');
assert(mod(numel(varargin), 2) == 0, 'PipeCompare:TestProgress', 'Expected name/value pairs.');
title = 'Working'; message = '';
for k = 1:2:numel(varargin)
    key = lower(char(varargin{k}));
    switch key
        case 'title', title = char(string(varargin{k+1}));
        case 'message', message = char(string(varargin{k+1}));
        case 'indeterminate'
            assert(any(strcmpi(string(varargin{k+1}), ["on", "off"])), ...
                'PipeCompare:TestProgress', 'Expected Indeterminate on/off.');
        otherwise
            error('PipeCompare:TestProgress', 'Unsupported dialog option: %s', key);
    end
end
fprintf('[uiprogressdlg] %s: %s\n', title, message);
dlg = figure('Visible', 'off', 'HandleVisibility', 'off', 'NumberTitle', 'off', ...
    'Name', title, 'Tag', 'PipeCompareTestProgress');
end
