function func_plotGenSSIdentStiffVsLIso(savePath,plotParams,material,pos,...
    ssCurveIdentPlotOpts,identStiffGenSS,genSSCurveOpts)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 29/1/2019
% Date Edited: 27/6/2019
%
% Plots the identified stiffness components for each equation set using the
% generalised stress-strain curve method.

    % Create the plot properties struct to format figures
    plotProps = func_initPlotPropsStruct(ssCurveIdentPlotOpts.formatType);

    %----------------------------------------------------------------------
    % PLOT: Generalised SS Curves, Plot Identfied Stiffness Over Length,
    % Isotropic Case
    labelStrs.x = 'Position, $x$ (mm)';
    unitConv = ssCurveIdentPlotOpts.unitConv;
    
    % Qxx
    ssCurveIdentPlotOpts.QSRef = material.Qxx;
    
    labelStrs.t = ['$Q_{xx,avg} =$ ',...
        sprintf('%3.1f',identStiffGenSS.QxxAvgOverL_Eq12*unitConv),'$~GPa$'];
    labelStrs.y = '$Q_{xx}$, Eq.1+2 (GPa)';
    hf = func_plotStiffnessVsLength(labelStrs,ssCurveIdentPlotOpts,pos,...
        identStiffGenSS.QxxVsL_Eq12,identStiffGenSS.QxxAvgOverL_Eq12,genSSCurveOpts);
    saveFile = [savePath,'\','GenSS_QVsL_Qxx_Eq12'];
    func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)
    
    labelStrs.t = ['$Q_{xx,avg} =$ ',...
        sprintf('%3.1f',identStiffGenSS.QxxAvgOverL_Eq13*unitConv),'$~GPa$'];
    labelStrs.y = '$Q_{xx}$, Eq.1+3 (GPa)';
    hf = func_plotStiffnessVsLength(labelStrs,ssCurveIdentPlotOpts,pos,...
        identStiffGenSS.QxxVsL_Eq13,identStiffGenSS.QxxAvgOverL_Eq13,genSSCurveOpts);
    saveFile = [savePath,'\','GenSS_QVsL_Qxx_Eq13'];
    func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)
    
    labelStrs.t = ['$Q_{xx,avg} =$ ',...
        sprintf('%3.1f',identStiffGenSS.QxxAvgOverL_Eq23*unitConv),'$~GPa$'];
    labelStrs.y = '$Q_{xx}$, Eq.2+3 (GPa)';
    hf = func_plotStiffnessVsLength(labelStrs,ssCurveIdentPlotOpts,pos,...
        identStiffGenSS.QxxVsL_Eq23,identStiffGenSS.QxxAvgOverL_Eq23,genSSCurveOpts);
    saveFile = [savePath,'\','GenSS_QVsL_Qxx_Eq23'];
    func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)
    
    % Qxy
    ssCurveIdentPlotOpts.QSRef = material.Qxy;
    
    labelStrs.t = ['$Q_{xy,avg} =$ ',...
        sprintf('%3.1f',identStiffGenSS.QxyAvgOverL_Eq12*unitConv),'$~GPa$'];
    labelStrs.y = '$Q_{xy}$, Eq.1+2 (GPa)';
    hf = func_plotStiffnessVsLength(labelStrs,ssCurveIdentPlotOpts,pos,...
        identStiffGenSS.QxyVsL_Eq12,identStiffGenSS.QxyAvgOverL_Eq12,genSSCurveOpts);
    saveFile = [savePath,'\','GenSS_QVsL_Qxy_Eq12'];
    func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)
    
    labelStrs.t = ['$Q_{xy,avg} =$ ',...
        sprintf('%3.1f',identStiffGenSS.QxyAvgOverL_Eq13*unitConv),'$~GPa$'];
    labelStrs.y = '$Q_{xy}$, Eq.1+3 (GPa)';
    hf = func_plotStiffnessVsLength(labelStrs,ssCurveIdentPlotOpts,pos,...
        identStiffGenSS.QxyVsL_Eq13,identStiffGenSS.QxyAvgOverL_Eq13,genSSCurveOpts);
    saveFile = [savePath,'\','GenSS_QVsL_Qxy_Eq13'];
    func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)
    
    labelStrs.t = ['$Q_{xy,avg} =$ ',...
        sprintf('%3.1f',identStiffGenSS.QxyAvgOverL_Eq23*unitConv),'$~GPa$'];
    labelStrs.y = '$Q_{xy}$, Eq.2+3 (GPa)';
    hf = func_plotStiffnessVsLength(labelStrs,ssCurveIdentPlotOpts,pos,...
        identStiffGenSS.QxyVsL_Eq23,identStiffGenSS.QxyAvgOverL_Eq23,genSSCurveOpts);
    saveFile = [savePath,'\','GenSS_QVsL_Qxy_Eq23'];
    func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)

end

