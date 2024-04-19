function func_plotSGIdentStiffVsLIso(savePath,plotParams,...
    ssCurveIdentPlotOpts,pos,material,identStiffSG,identOpts)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 28/1/2019
% Date Edited: 27/6/2019
%
% Plots identified stiffness against length using the standard stress-gauge

    % Set plot properties for figure formatting
    plotProps = func_initPlotPropsStruct(ssCurveIdentPlotOpts.formatType);
    labelStrs.x = 'Position, $x$ ($mm$)';
    unitConv = ssCurveIdentPlotOpts.unitConv; 
    
    % Plot identified stiffness vs length from SS curves as Qxx
    labelStrs.t = ['$Q_{xx,avg} =$ ',...
        sprintf('%3.1f',identStiffSG.QxxAvgOverL*unitConv),'$~GPa$'];
    labelStrs.y = 'Stiffness, $Q_{xx}$ (GPa)';
    ssCurveIdentPlotOpts.QSRef = material.Qxx;
    hf = func_plotStiffnessVsLength(labelStrs,ssCurveIdentPlotOpts,pos,...
        identStiffSG.QxxVsL,identStiffSG.QxxAvgOverL,identOpts);
    saveFile = [savePath,'\','SG_StiffVsL_Qxx'];
    func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)
    
    % Plot identified stiffness vs length from SS curves as E
    labelStrs.t = ['$E_{avg} =$ ',...
        sprintf('%3.1f',identStiffSG.ExxAvgOverL*unitConv),'$~GPa$'];
    labelStrs.y = '$E$ (GPa)';
    ssCurveIdentPlotOpts.QSRef = material.Exx;
    hf = func_plotStiffnessVsLength(labelStrs,ssCurveIdentPlotOpts,pos,...
        identStiffSG.ExxVsL,identStiffSG.ExxAvgOverL,identOpts);
    saveFile = [savePath,'\','SG_StiffVsL_Exx'];
    func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)
end

