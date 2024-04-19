function func_plotGenSSIdentStiffVsLOrtho(savePath,plotParams,material,pos,...
    ssCurveIdentPlotOpts,angledSlice1,angledSlice2,identStiffGenSS,genSSCurveOpts)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 29/1/2019
% Date Edited: 27/6/2019
%
% Plots the identified stiffness components for each equation set using the
% generalised stress-strain curve method.

    % Create the plot properties struct to format figures
    plotProps = func_initPlotPropsStruct(ssCurveIdentPlotOpts.formatType);
    unitConv = ssCurveIdentPlotOpts.unitConv; 
    ssCurveIdentPlotOpts.plotAvgPosRange = false; % Turn this off for now, implement later
    
    %//////////////////////////////////////////////////////////////////
    % Slice 1: plot identified stiffness over the length
    if angledSlice1.xMaxInd > 0
        ssCurveIdentPlotOpts.specifyLRange = false;
        labelStrs.x = 'Position, $x$ (mm)';
        
        % Equations 14.1 to 14.3
        % Eq 14.1 Q22
        labelStrs.t = ['$Q_{22,avg} =$ ',...
            sprintf('%3.1f',identStiffGenSS.Q22AvgOverL_Eq141*unitConv),'$~GPa$'];
        ssCurveIdentPlotOpts.QSRef = material.Q22;
        labelStrs.y = '$Q_{22}$, Eq.14.1 (GPa)';
        hf = func_plotStiffnessVsLength(labelStrs,ssCurveIdentPlotOpts,pos,...
            identStiffGenSS.Q22VsL_Eq141,identStiffGenSS.Q22AvgOverL_Eq141,genSSCurveOpts);
        saveFile = [savePath,'\','GenSS_QVsL_Q22_Eq141'];
        func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)
        
        % Eq 14.2 Q12
        labelStrs.t = ['$Q_{12,avg} =$ ',...
            sprintf('%3.1f',identStiffGenSS.Q12AvgOverL_Eq142*unitConv),'$~GPa$'];
        ssCurveIdentPlotOpts.QSRef = material.Q12;
        labelStrs.y = '$Q_{12}$, Eq.14.2 (GPa)';
        hf = func_plotStiffnessVsLength(labelStrs,ssCurveIdentPlotOpts,pos,...
            identStiffGenSS.Q12VsL_Eq142,identStiffGenSS.Q12AvgOverL_Eq142,genSSCurveOpts);
        saveFile = [savePath,'\','GenSS_QVsL_Q12_Eq142'];
        func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)

        % Eq 14.3 Q66
        labelStrs.t = ['$Q_{66,avg} =$ ',...
            sprintf('%3.1f',identStiffGenSS.Q66AvgOverL_Eq143*unitConv),'$~GPa$'];
        ssCurveIdentPlotOpts.QSRef = material.Q66;
        labelStrs.y = '$Q_{66}$, Eq.14.3 (GPa)';
        hf = func_plotStiffnessVsLength(labelStrs,ssCurveIdentPlotOpts,pos,...
            identStiffGenSS.Q66VsL_Eq143,identStiffGenSS.Q66AvgOverL_Eq143,genSSCurveOpts);
        saveFile = [savePath,'\','GenSS_QVsL_Q66_Eq143'];
        func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)
    end

    %//////////////////////////////////////////////////////////////////
    % Slice 2: plot identified stiffness over the length
    if angledSlice2.xMaxInd > 0
        ssCurveIdentPlotOpts.specifyLRange = false;
        labelStrs.x = 'Position, $x$ (mm)';
        % Equations 14.1 to 14.3
        % Eq 18.1 Q11
        labelStrs.t = ['$Q_{11,avg} =$ ',...
            sprintf('%3.1f',identStiffGenSS.Q11AvgOverL_Eq181*unitConv),'$~GPa$'];
        ssCurveIdentPlotOpts.QSRef = material.Q11;
        labelStrs.y = '$Q_{11}$, Eq.18.1 (GPa)';
        hf = func_plotStiffnessVsLength(labelStrs,ssCurveIdentPlotOpts,pos,...
            identStiffGenSS.Q11VsL_Eq181,identStiffGenSS.Q11AvgOverL_Eq181,genSSCurveOpts);
        saveFile = [savePath,'\','GenSS_QVsL_Q11_Eq181'];
        func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)

        % Eq 18.2 Q12
        labelStrs.t = ['$Q_{12,avg} =$ ',...
            sprintf('%3.1f',identStiffGenSS.Q12AvgOverL_Eq182*unitConv),'$~GPa$'];
        ssCurveIdentPlotOpts.QSRef = material.Q12;
        labelStrs.y = '$Q_{12}$, Eq.18.2 (GPa)';
        hf = func_plotStiffnessVsLength(labelStrs,ssCurveIdentPlotOpts,pos,...
            identStiffGenSS.Q12VsL_Eq182,identStiffGenSS.Q12AvgOverL_Eq182,genSSCurveOpts);
        saveFile = [savePath,'\','GenSS_QVsL_Q12_Eq182'];
        func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)

        % Eq 18.3 Q66
        labelStrs.t = ['$Q_{66,avg} =$ ',...
            sprintf('%3.1f',identStiffGenSS.Q66AvgOverL_Eq183*unitConv),'$~GPa$'];
        ssCurveIdentPlotOpts.QSRef = material.Q66;
        labelStrs.y = '$Q_{66}$, Eq.18.3 (GPa)';
        hf = func_plotStiffnessVsLength(labelStrs,ssCurveIdentPlotOpts,pos,...
            identStiffGenSS.Q66VsL_Eq183,identStiffGenSS.Q66AvgOverL_Eq183,genSSCurveOpts);
        saveFile = [savePath,'\','GenSS_QVsL_Q66_Eq183'];
        func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)
    end

end

