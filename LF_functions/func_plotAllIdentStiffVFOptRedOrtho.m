function func_plotAllIdentStiffVFOptRedOrtho(savePath,plotParams,...
    time,material,VFPlotOpts,identStiffVFOpt,VFOpts)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 28/1/2019
% Date Edited: 7/6/2019
%
% Plots all identified stiffness components against time for the reduced
% orthotropic case using optimised virtual fields.

    % Set plot properties for figure formatting
    plotProps = func_initPlotPropsStruct(VFPlotOpts.formatType);
    labelStrs.x = 'Time $t$, ($\mu s$)';
    
    % Plot identified Qxx vs time
    labelStrs.t = ['$Q_{xx,avg} =$ ',...
        sprintf('%3.1f',identStiffVFOpt.QxxAvgOverT*VFPlotOpts.unitConv),'$~GPa$'];
    labelStrs.y = '$Q_{xx}$ (GPa)';
    VFPlotOpts.QSRef = material.Qxx;
    hf = func_plotStiffnessVsTime(labelStrs,VFPlotOpts,time,...
        identStiffVFOpt.QxxVsT,identStiffVFOpt.QxxAvgOverT,VFOpts);
    saveFile = [savePath,'\','VFOpt_StiffVsT_Qxx'];
    func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)

    % Plot identified Exx vs time
    labelStrs.t = ['$E_{xx,avg} =$ ',...
        sprintf('%3.1f',identStiffVFOpt.ExxAvgOverT*VFPlotOpts.unitConv),'$~GPa$'];
    labelStrs.y = '$E_{xx}$ (GPa)';
    VFPlotOpts.QSRef = material.Exx;
    hf = func_plotStiffnessVsTime(labelStrs,VFPlotOpts,time,...
        identStiffVFOpt.ExxVsT,identStiffVFOpt.ExxAvgOverT,VFOpts);
    saveFile = [savePath,'\','VFOpt_StiffVsT_Exx'];
    func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)
end

