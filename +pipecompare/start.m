function w = start()
%START Open the PipeCompare window.
%   pipecompare.start        the window: open the data, add channel
%                            locations, choose what to compare (Simple or
%                            Advanced), run, see the results and export, all
%                            in one window (pipecompare.gui.Main)
%   w = pipecompare.start;   ... and return it
%   EEGLAB must be on the MATLAB path; its main window need not be open.
%   EEGLAB's menu Tools > PipeCompare stays as it is.
win = pipecompare.gui.Main.open();
if nargout > 0, w = win; end
end
