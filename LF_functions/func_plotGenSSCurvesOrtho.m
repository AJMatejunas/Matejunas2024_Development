function func_plotGenSSCurvesOrtho(savePath,plotParams,ssCurvePlotParams,ssCurveIdentPlotOpts...
    ,pos,angledSlice1,angledSlice2,genSSCurves)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 31/1/2019
% Date Edited: 31/1/2019
%
% Plots stress- strain curves using the generalised formulation at
% different angled slices along the sample length.

    % Create the plot properties struct to format figures
    plotProps = func_initPlotPropsStruct(ssCurvePlotParams.formatType);
    

    %//////////////////////////////////////////////////////////////////
    % Slice 1: plot generalised stress-strain curves
    if angledSlice1.xMaxInd > 0
        % Calculate indices to average over the middle slices
        ssCurvePlotParams.plotOnlyFittedRegion = false;
        ssCurvePlotParams.plotFittedCurve = false;
        ssCurveIdentPlotOpts.indLRange = 1:angledSlice1.xMaxInd;
        ssCurvePlotParams.locXInd = round(angledSlice1.xMaxInd*ssCurvePlotParams.locXPcVec); 
        labelStrs.y = 'Accel. Avg.';
        labelStrs.x = 'Strain. Avg.';

        % Equations 14.1 to 14.3
        % Eq 14.1 Q22
        labelStrs.t = '$Q_{22}$ Eq.14.1';
        hf = func_plotStressStrainCurvesGeneric(labelStrs,ssCurvePlotParams,pos,...
            genSSCurves.Q22_accelAvg_Eq141,genSSCurves.Q22_strainAvg_Eq141);
        saveFile = [savePath,'\','GenSS_Curve_Q22_Eq141'];
        func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)

        % Eq 14.2, Q12
        labelStrs.t = '$Q_{12}$ Eq.14.2';
        hf = func_plotStressStrainCurvesGeneric(labelStrs,ssCurvePlotParams,pos,...
            genSSCurves.Q12_accelAvg_Eq142,genSSCurves.Q12_strainAvg_Eq142);
        saveFile = [savePath,'\','GenSS_Curve_Q12_Eq142'];
        func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)

        % Eq 14.3, Q66
        labelStrs.t = '$Q_{66}$ Eq.14.3';
        hf = func_plotStressStrainCurvesGeneric(labelStrs,ssCurvePlotParams,pos,...
            genSSCurves.Q66_accelAvg_Eq143,genSSCurves.Q66_strainAvg_Eq143);
        saveFile = [savePath,'\','GenSS_Curve_Q66_Eq143'];
        func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)
    end

    %//////////////////////////////////////////////////////////////////
    % Slice 2: plot generalised stress-strain curves
    if angledSlice2.xMaxInd > 0
        % Calculate indices to average over the middle slices
        ssCurvePlotParams.plotOnlyFittedRegion = false;
        ssCurvePlotParams.plotFittedCurve = false;
        ssCurveIdentPlotOpts.indLRange = 1:angledSlice2.xMaxInd;
        ssCurvePlotParams.locXInd = round(angledSlice2.xMaxInd*ssCurvePlotParams.locXPcVec); 
        labelStrs.y = 'Accel. Avg.';
        labelStrs.x = 'Strain. Avg.';

        % Equations 18.1 to 18.3
        % Eq 18.1 Q11
        labelStrs.t = '$Q_{11}$ Eq.18.1';
        hf = func_plotStressStrainCurvesGeneric(labelStrs,ssCurvePlotParams,pos,...
            genSSCurves.Q11_accelAvg_Eq181,genSSCurves.Q11_strainAvg_Eq181);
        saveFile = [savePath,'\','GenSS_Curve_Q11_Eq181'];
        func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)

        % Eq 18.2, Q12
        labelStrs.t = '$Q_{12}$ Eq.18.2';
        hf = func_plotStressStrainCurvesGeneric(labelStrs,ssCurvePlotParams,pos,...
            genSSCurves.Q12_accelAvg_Eq182,genSSCurves.Q12_strainAvg_Eq182);
        saveFile = [savePath,'\','GenSS_Curve_Q12_Eq182'];
        func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)

        % Eq 18.3, Q66
        labelStrs.t = '$Q_{66}$ Eq.18.3';
        hf = func_plotStressStrainCurvesGeneric(labelStrs,ssCurvePlotParams,pos,...
            genSSCurves.Q66_accelAvg_Eq183,genSSCurves.Q66_strainAvg_Eq183);
        saveFile = [savePath,'\','GenSS_Curve_Q66_Eq183'];
        func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)
    end
end

