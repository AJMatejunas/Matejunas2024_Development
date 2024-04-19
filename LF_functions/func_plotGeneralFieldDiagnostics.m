function func_plotGeneralFieldDiagnostics(savePath,plotParams,...
    specimen,time,pos,disp,accel,strain,strainRate,stress)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 28/1/2019
% Date Edited: 28/1/2019
%
% Plots field averages and loading pulse
    
    % Set plot properties for figure formatting
    plotProps = func_initPlotPropsStruct(plotParams.formatType);
    
    % Plot field averages over time on a single plot
    fprintf('\tPlotting field averages over time.\n')
    plotParams.Rows = 2;
    plotParams.Cols = 2;
    plotParams.tRange = 1:time.numFrames;
    hf = func_plotFieldAvgVsTime(plotParams,time,pos,disp,accel,strain,strainRate);
    saveFile = [savePath,'\','GenDiag_FieldAveragesVsTime'];
    func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)
    
    % Plot the loading pulse
    plotParams.Rows = 1;
    plotParams.Cols = 1;
    labelStrs.y1Str = 'Force $\overline{F_{x}}^{y}$ ($kN$)';
    labelStrs.y2Str = 'Stress $\overline{\sigma_{xx}}^{y}$ ($MPa$)';
    labelStrs.xStr = 'Time ($\mu s$)';
    labelStrs.legStrs = {'Exp.'};
    hf = func_plotInputLoadingPulse(plotParams,labelStrs,time,specimen,stress);
    saveFile = [savePath,'\','GenDiag_LoadingPulse'];
    func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)

end

