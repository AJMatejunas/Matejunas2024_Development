function func_plotAllIdentStiffVFManIso(savePath,plotParams,...
    time,material,VFPlotOpts,identStiffVFMan,VFOpts)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 28/1/2019
% Date Edited: 28/1/2019
%
% Plots all identified stiffness components against time for the isotropic
% case using manually defined virtual fields.
        
    % Set plot properties for figure formatting
    plotProps = func_initPlotPropsStruct(VFPlotOpts.formatType);
    labelStrs.x = 'Time $t$, ($\mu s$)';

    % Plot identified Qxx vs time
    labelStrs.t = ['$Q_{xx,avg} =$ ',...
        sprintf('%3.1f',identStiffVFMan.QxxAvgOverT*VFPlotOpts.unitConv),'$~GPa$'];
    labelStrs.y = '$Q_{xx}$ (GPa)';
    VFPlotOpts.QSRef = material.Qxx;
    hf = func_plotStiffnessVsTime(labelStrs,VFPlotOpts,time,...
        identStiffVFMan.QxxVsT,identStiffVFMan.QxxAvgOverT,VFOpts);
    saveFile = [savePath,'\','VFMan_StiffVsT_Qxx'];
    func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)

    % Plot identified Qxx vs time
    labelStrs.t = ['$Q_{xy,avg} =$ ',...
        sprintf('%3.1f',identStiffVFMan.QxyAvgOverT*VFPlotOpts.unitConv),'$~GPa$'];
    labelStrs.y = '$Q_{xy}$ (GPa)';
    VFPlotOpts.QSRef = material.Qxy;
    hf = func_plotStiffnessVsTime(labelStrs,VFPlotOpts,time,...
        identStiffVFMan.QxyVsT,identStiffVFMan.QxyAvgOverT,VFOpts);
    saveFile = [savePath,'\','VFMan_StiffVsT_Qxy'];
    func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)

     % Plot identified Exx vs time
     labelStrs.t = ['$E_{avg} =$ ',...
        sprintf('%3.1f',identStiffVFMan.ExxAvgOverT*VFPlotOpts.unitConv),'$~GPa$'];
    labelStrs.y = '$E$ (GPa)';
    VFPlotOpts.QSRef = material.Exx;
    hf = func_plotStiffnessVsTime(labelStrs,VFPlotOpts,time,...
        identStiffVFMan.ExxVsT,identStiffVFMan.ExxAvgOverT,VFOpts);
    saveFile = [savePath,'\','VFMan_StiffVsT_E'];
    func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)

    % Plot identified Nuxy vs length
    labelStrs.t = ['$\nu_{avg} =$ ',...
        sprintf('%.3f',identStiffVFMan.NuxyAvgOverT)];
    labelStrs.y = '$\nu$ $(m.m^{-1})$';
    VFPlotOpts.unitConv = 1;
    VFPlotOpts.specifyAxisLims = true;
    VFPlotOpts.QSRef = material.nuxy;
    hf = func_plotStiffnessVsTime(labelStrs,VFPlotOpts,time,...
        identStiffVFMan.NuxyVsT,identStiffVFMan.NuxyAvgOverT,VFOpts);
    saveFile = [savePath,'\','VFMan_StiffVsT_Nu'];
    func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)

    VFPlotOpts.unitConv = 10^-9;
    VFPlotOpts.specifyAxisLims = false;
end

