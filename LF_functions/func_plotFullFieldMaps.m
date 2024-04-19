function func_plotFullFieldMaps(plotParams,savePath,time,disp,strain)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 9/3/2017
%
% Plots frames of the displacement and strain data

plotProps = func_initPlotPropsStruct();
figIndex = 1;

labelStrs.t = 'Displacement X, \delta_{x}, (mm)';
plotVar = disp.x*10^3;
plotParams.cutPxX = 10;
plotParams.cutPxY = 0;
plotParams.cRange = [-1.8,0];
hf = func_plotFrameMaps(plotParams,labelStrs,time,plotVar);
if plotParams.saveFigs
    print(hf,[savePath,'\',num2str(figIndex),'_DispXMaps',],plotProps.format,'-r0')
    figIndex = figIndex+1;
end

%{
labelStrs.t = 'Displacement Y, \delta_{y}, (mm)';
plotVar = disp.y*10^3;
plotParams.cutPxX = 0;
plotParams.cutPxY = 10;
hf = func_plotFrameMaps(plotParams,labelStrs,time,plotVar);
if plotParams.saveFigs
    print(hf,[savePath,'\',num2str(figIndex),'_DispYMaps',],plotProps.format,'-r0')
    figIndex = figIndex+1;
end
%}

labelStrs.t = 'Strain X, \epsilon_{x}, (mm/m)';
plotVar = strain.x*10^3;
plotParams.cutPxX = 20;
plotParams.cutPxY = 10;
plotParams.cRange = [-20,20];
hf = func_plotFrameMaps(plotParams,labelStrs,time,plotVar);
if plotParams.saveFigs
    print(hf,[savePath,'\',num2str(figIndex),'_StrainXMaps',],plotProps.format,'-r0')
    figIndex = figIndex+1;
end

%{
labelStrs.t = 'Strain Y, \epsilon_{y}, (mm/m)';
plotVar = strain.y*10^3;
plotParams.cutPxX = 0;
plotParams.cutPxY = 20;
hf = func_plotFrameMaps(plotParams,labelStrs,time,plotVar);
if plotParams.saveFigs
    print(hf,[savePath,'\',num2str(figIndex),'_StrainYMaps',],plotProps.format,'-r0')
    figIndex = figIndex+1;
end

labelStrs.t = 'Strain XY, \epsilon_{xy}, (mm/m)';
plotVar = strain.s*10^3;
plotParams.cutPxX = 0;
plotParams.cutPxY = 0;
hf = func_plotFrameMaps(plotParams,labelStrs,time,plotVar);
if plotParams.saveFigs
    print(hf,[savePath,'\',num2str(figIndex),'_StrainXYMaps',],plotProps.format,'-r0')
    figIndex = figIndex+1;
end
%}

end