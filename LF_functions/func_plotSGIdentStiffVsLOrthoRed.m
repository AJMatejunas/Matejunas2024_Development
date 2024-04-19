function func_plotSGIdentStiffVsLOrthoRed(savePath,plotParams,...
    ssCurveIdentPlotOpts,pos,material,identStiffSG,identOpts)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 28/1/2019
% Date Edited: 9/6/2019
%
% Plots identified stiffness against length using the standard stress-gauge

    % Set plot properties for figure formatting
    plotProps = func_initPlotPropsStruct(ssCurveIdentPlotOpts.formatType);
    labelStrs.x = 'Position, $x$ ($mm$)';
    unitConv = ssCurveIdentPlotOpts.unitConv;
    
    % Plot the identified stiffness with the SG as a function of length
    labelStrs.t = ['$E_{xx,avg} =$ ',...
        sprintf('%3.1f',identStiffSG.ExxAvgOverL*unitConv),'$~GPa$'];
    labelStrs.y = '$E_{xx}$ (GPa)';
    ssCurveIdentPlotOpts.QSRef = material.Exx;
    hf = func_plotStiffnessVsLength(labelStrs,ssCurveIdentPlotOpts,pos,...
        identStiffSG.ExxVsL,identStiffSG.ExxAvgOverL,identOpts);
    saveFile = [savePath,'\','SG_StiffVsL_Exx'];
    func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)
end

