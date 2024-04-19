function hf = func_createFigure(plotProps,plotParams)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 19/2/2018
%
% Creates a matlab figure and sets size properties, returns the handle to
% the figure.

% If not specified assume we only have a single subplot
if nargin < 2
    plotParams.Rows = 1;
    plotParams.Cols = 1;
end

hf = figure('Visible','On');
set(hf,'Resize','on','Units',plotProps.units,'PaperPositionMode','auto')
set(hf,'Position', [plotProps.locXcm,plotProps.locYcm,...
    plotParams.Cols*plotProps.sizePerFigXcm,plotParams.Rows*plotProps.sizePerFigYcm])
set(hf,'Color',[1 1 1])

end

