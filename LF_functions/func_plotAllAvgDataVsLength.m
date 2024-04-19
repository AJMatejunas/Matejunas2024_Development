function func_plotAllAvgDataVsLength(plotParams,savePath,time,pos,...
    disp,accel,strain,stress)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 9/3/2017
% Plots average over width data (disp,accel,strain,stress) vs length for a
% a selected number of points along the specimen length

figIndex = 1;
plotProps = func_initPlotPropsStruct();
labelStrs.x = 'X Position (mm)';
plotParams.markStrs = {'-xb'};
pos.x = pos.x*10^3;
pos.y = pos.y*10^3;
posVar = {pos};

%----------------------------------------------------------------------
% Plot: Avg Displacement Vs Length
labelStrs.y = 'Disp Avg, d_{x} (mmm)';
plotVar = {disp.xAvg*10^3};
hf = func_plotAvgDataVsLength(labelStrs,time,plotParams,posVar,plotVar);
print(hf,[savePath,'\',num2str(figIndex),'_dispXAvg_vsLeng',],plotProps.format,'-r0')
figIndex = figIndex+1;

%----------------------------------------------------------------------
% Plot: Avg Accel Vs Length
labelStrs.y = 'Avg (Line) Accel, a_{x} (m/s^2)';
plotVar = {accel.xAvg};
hf = func_plotAvgDataVsLength(labelStrs,time,plotParams,posVar,plotVar);
print(hf,[savePath,'\',num2str(figIndex),'_accelXAvg_vsLeng',],plotProps.format,'-r0')
figIndex = figIndex+1;

%--------------------------------------------------------------------------
% Plot: Avg Stress Vs Length
labelStrs.y = 'Avg Stress, \sigma_{x} (MPa)';
plotVar = {stress.xAvg*10^-6};
hf = func_plotAvgDataVsLength(labelStrs,time,plotParams,posVar,plotVar);
print(hf,[savePath,'\',num2str(figIndex),'_stressXAvg_vsLeng',],plotProps.format,'-r0')
figIndex = figIndex+1;

%----------------------------------------------------------------------
% Plot: Avg Strain X Vs Length
labelStrs.y = 'Avg Strain, \epsilon_{x}, (mm/m)';
plotVar = {strain.xAvg*10^3};
hf = func_plotAvgDataVsLength(labelStrs,time,plotParams,posVar,plotVar);
print(hf,[savePath,'\',num2str(figIndex),'_strainXAvg_vsLeng',],plotProps.format,'-r0')
figIndex = figIndex+1;


end

