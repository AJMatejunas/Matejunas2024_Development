function func_plotAllFullFieldMaps(plotParams,savePath,pos,time,disp,accel,strain,strainRate)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 9/3/2017
%
% Plots frames of the displacement and strain data

% Create the plot properties struct for formatting
plotProps = func_initPlotPropsStruct(plotParams.formatType);
figIndex = 1;

% Specify the plot parameters for all maps
plotParams.cAxisType = 'MaxQuantile';
plotParams.cutEdgePx = true;
plotParams.imageAxis = true;
% Label strings for the image axis
labelStrs.x = 'X (mm)';
labelStrs.y = 'Y (mm)';

%--------------------------------------------------------------------------
% Plot Displacement Maps
labelStrs.t = 'Displacement X, \delta_{x}, (mm)';
[deformDisp,rigidDisp] = func_calcDispRemoveRigidBody(disp.x);
plotVar = deformDisp*10^3;
hf = func_plotFrameMaps(plotParams,labelStrs,pos,time,plotVar);
print(hf,[savePath,'\',num2str(figIndex),'_DispXMaps',],plotProps.format,'-r0')
figIndex = figIndex+1;

labelStrs.t = 'Displacement Y, \delta_{y}, (mm)';
plotVar = disp.y*10^3;
plotParams.cutPxX = 0;
plotParams.cutPxY = 10;
hf = func_plotFrameMaps(plotParams,labelStrs,pos,time,plotVar);
print(hf,[savePath,'\',num2str(figIndex),'_DispYMaps',],plotProps.format,'-r0')
figIndex = figIndex+1;


%--------------------------------------------------------------------------
% Plot Acceleration Maps
labelStrs.t = 'Accel X, a_{x}, (m.s^{-2})';
plotVar = accel.x;
hf = func_plotFrameMaps(plotParams,labelStrs,pos,time,plotVar);
print(hf,[savePath,'\',num2str(figIndex),'_AccelXMaps',],plotProps.format,'-r0')
figIndex = figIndex+1;

labelStrs.t = 'Accel Y, a_{y}, (m.s^{-2})';
plotVar = accel.y;
hf = func_plotFrameMaps(plotParams,labelStrs,pos,time,plotVar);
print(hf,[savePath,'\',num2str(figIndex),'_AccelYMaps',],plotProps.format,'-r0')
figIndex = figIndex+1;

%--------------------------------------------------------------------------
% Plot Strain Maps
labelStrs.t = 'Strain X, \epsilon_{x}, (mm.m^{-1})';
plotVar = strain.x*10^3;
plotParams.cutPxX = plotParams.cutPxX+5;
plotParams.cutPxY = 0;
hf = func_plotFrameMaps(plotParams,labelStrs,pos,time,plotVar);
print(hf,[savePath,'\',num2str(figIndex),'_StrainXMaps',],plotProps.format,'-r0')
figIndex = figIndex+1;

labelStrs.t = 'Strain Y, \epsilon_{y}, (mm/m)';
plotVar = strain.y*10^3;
plotParams.cutPxX = 0;
plotParams.cutPxY = plotParams.cutPxX+5;
hf = func_plotFrameMaps(plotParams,labelStrs,pos,time,plotVar);
print(hf,[savePath,'\',num2str(figIndex),'_StrainYMaps',],plotProps.format,'-r0')
figIndex = figIndex+1;

labelStrs.t = 'Strain XY, \epsilon_{xy}, (mm/m)';
plotVar = strain.s*10^3;
plotParams.cutPxX = 0;
plotParams.cutPxY = 0;
hf = func_plotFrameMaps(plotParams,labelStrs,pos,time,plotVar);
print(hf,[savePath,'\',num2str(figIndex),'_StrainXYMaps',],plotProps.format,'-r0')
figIndex = figIndex+1;

%--------------------------------------------------------------------------
% Plot Strain Rate Maps
labelStrs.t = 'Strain Rate X, d\epsilon_{x}/dt, (s^{-1})';
plotVar = strainRate.x;
plotParams.cutPxX = 0;
hf = func_plotFrameMaps(plotParams,labelStrs,pos,time,plotVar);
print(hf,[savePath,'\',num2str(figIndex),'_StrainRateXMaps',],plotProps.format,'-r0')
figIndex = figIndex+1;
 
end