function func_plotFullFieldVideos(savePath)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 13/3/2017
%
% Plots 2D frames of data as an image sequence


plotProps = func_initPlotPropsStruct();
figIndex = 1;

labelStrs.t = 'Displacement X, \delta_{x}, (mm)';
plotVar = disp.x*10^3;
plotParams.cutPxX = 10;
plotParams.cutPxY = 0;
plotParams.cRange = [-1.8,0];



if plotParams.saveFigs
    print(hf,[savePath,'\',num2str(figIndex),'_DispXMaps',],plotProps.format,'-r0')
    figIndex = figIndex+1;
end


end