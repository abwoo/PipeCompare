function w = pcomp()
%PCOMP Open the PipeCompare window (the same as pipecompare.start).
%   pcomp                the window: open the data, add channel locations,
%                        choose what to compare (Simple or Advanced), run,
%                        see the results and export, all in one window
%   w = pcomp;           ... and return it
%   PipeCompare and EEGLAB must be on the MATLAB path (once: start EEGLAB,
%   then savepath); EEGLAB's main window need not be open.
win = pipecompare.start();
if nargout > 0, w = win; end
end
