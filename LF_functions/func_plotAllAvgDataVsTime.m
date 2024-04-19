function func_plotAllAvgDataVsTime(plotParams,savePath,time,pos,disp,accel,strain,stress)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 13/3/2017
%
% Plots average over width data (disp,accel,strain,stress) vs time for a
% a selected number of points along the specimen length

figIndex = 1;
plotProps = func_initPlotPropsStruct();
labelStrs.x = 'Time, (\mus)';
plotParams.markStrs = {'-xb'};
posVar = {pos};

%----------------------------------------------------------------------
% Plot: Avg Displacement Vs Time
labelStrs.y = 'Disp Avg, d_{x} (mm)';    
plotVar = {disp.xAvg*10^3};
hf = func_plotAvgDataVsTime(labelStrs,time,plotParams,posVar,plotVar);
print(hf,[savePath,'\',num2str(figIndex),'_dispXAvg_vsTime',],plotProps.format,'-r0')
figIndex = figIndex+1;

%----------------------------------------------------------------------
% Plot: Avg Accel Vs Time
labelStrs.y = 'Avg (Line) Accel, a_{x} (m/s^2)';
plotVar = {accel.xAvg};
hf = func_plotAvgDataVsTime(labelStrs,time,plotParams,posVar,plotVar);
print(hf,[savePath,'\',num2str(figIndex),'_accelXAvg_vsTime',],plotProps.format,'-r0')
figIndex = figIndex+1;

%--------------------------------------------------------------------------
% Plot: Avg Stress Vs Time
labelStrs.y = 'Avg Stress, \sigma_{x} (MPa)';
plotVar = {stress.xAvg*10^-6};
hf = func_plotAvgDataVsTime(labelStrs,time,plotParams,posVar,plotVar);
print(hf,[savePath,'\',num2str(figIndex),'_stressXAvg_vsTime',],plotProps.format,'-r0')
figIndex = figIndex+1;

%----------------------------------------------------------------------
% Plot: Avg Strain X Vs Time
labelStrs.y = 'Avg Strain, \epsilon_{x}, (mm/m)';
plotVar = {strain.xAvg*10^3};
hf = func_plotAvgDataVsTime(labelStrs,time,plotParams,posVar,plotVar);
print(hf,[savePath,'\',num2str(figIndex),'_strainXAvg_vsTime',],plotProps.format,'-r0')
figIndex = figIndex+1;


end

