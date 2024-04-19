function func_plotAllIdentStiffVFOptIso(savePath,plotParams,...
    time,material,VFPlotOpts,identStiffVFOpt,VFOpts)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 28/1/2019
% Date Edited: 28/1/2019
%
% Plots all identified stiffness components against time for the isotropic
% case using optimised virtual fields.

    % Set plot properties for figure formatting
    plotProps = func_initPlotPropsStruct(VFPlotOpts.formatType);
    labelStrs.x = 'Time, $t$ ($\mu s$)';
    
    % Plot identified Qxx vs time
    labelStrs.t = ['$Q_{xx,avg} =$ ',...
        sprintf('%3.1f',identStiffVFOpt.QxxAvgOverT*VFPlotOpts.unitConv),'$~GPa$'];
    labelStrs.y = '$Q_{xx}$ (GPa)';
    VFPlotOpts.QSRef = material.Qxx*10^-9;
    hf = func_plotStiffnessVsTime(labelStrs,VFPlotOpts,time,...
        identStiffVFOpt.QxxVsT,identStiffVFOpt.QxxAvgOverT,VFOpts);
    saveFile = [savePath,'\','VFOpt_StiffVsT_Qxx'];
    func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)

    % Plot identified Qxy vs time
    labelStrs.t = ['$Q_{xy,avg} =$ ',...
        sprintf('%3.1f',identStiffVFOpt.QxyAvgOverT*VFPlotOpts.unitConv),'$~GPa$'];
    labelStrs.y = '$Q_{xy}$ (GPa)';
    VFPlotOpts.QSRef = material.Qxy*10^-9;
    hf = func_plotStiffnessVsTime(labelStrs,VFPlotOpts,time,...
        identStiffVFOpt.QxyVsT,identStiffVFOpt.QxyAvgOverT,VFOpts);
    saveFile = [savePath,'\','VFOpt_StiffVsT_Qxy'];
    func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)

     % Plot identified Exx vs time
     labelStrs.t = ['$E_{avg} =$ ',...
        sprintf('%3.1f',identStiffVFOpt.ExxAvgOverT*VFPlotOpts.unitConv),'$~GPa$'];
    labelStrs.y = '$E$ (GPa)';
    VFPlotOpts.QSRef = material.Exx*10^-9;
    hf = func_plotStiffnessVsTime(labelStrs,VFPlotOpts,time,...
        identStiffVFOpt.ExxVsT,identStiffVFOpt.ExxAvgOverT,VFOpts);
    saveFile = [savePath,'\','VFOpt_StiffVsT_E'];
    func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)

    % Plot identified Nuxy vs length
    labelStrs.t = ['$\nu_{avg} =$ ',...
        sprintf('%0.3f',identStiffVFOpt.NuxyAvgOverT)];
    labelStrs.y = '$\nu_{xy}$ $(m.m^{-1})$';
    VFPlotOpts.unitConv = 1;
    VFPlotOpts.specifyAxisLims = true;
    VFPlotOpts.QSRef = material.nuxy;
    hf = func_plotStiffnessVsTime(labelStrs,VFPlotOpts,time,...
        identStiffVFOpt.NuxyVsT,identStiffVFOpt.NuxyAvgOverT,VFOpts);
    saveFile = [savePath,'\','VFOpt_StiffVsT_Nu'];
    func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)

    VFPlotOpts.unitConv = 10^-9;
    VFPlotOpts.specifyAxisLims = false;
end

