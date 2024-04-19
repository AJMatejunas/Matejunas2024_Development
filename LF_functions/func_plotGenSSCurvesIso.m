function func_plotGenSSCurvesIso(savePath,plotParams,ssCurvePlotParams,...
    ssCurveIdentPlotOpts,pos,genSSCurves)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 29/1/2019
% Date Edited: 27/6/2019
%
% Plots stress- strain curves using the generalised formulation at
% different angled slices along the sample length.
    
    % Create the plot properties struct to format figures
    plotProps = func_initPlotPropsStruct(ssCurvePlotParams.formatType);

    %----------------------------------------------------------------------
    % PLOT: Generalised SS Curves, Isotropic Case - FULL CURVES
    fprintf('\tPlotting stress-strain curves.\n')
    ssCurveIdentPlotOpts.indLRange = 1:length(pos.x);
    ssCurvePlotParams.locXInd = round(length(pos.x)*ssCurvePlotParams.locXPcVec);
    
    % Plot the full stress-strain curves
    ssCurvePlotParams.plotOnlyFittedRegion = false;
    ssCurvePlotParams.plotFittedCurve = false;
    
    labelStrs.y = 'Accel. Avg.';
    labelStrs.x = 'Strain. Avg.';
    
    % Qxx
    labelStrs.t = '$Q_{xx}$ Eq.1+2';
    hf = func_plotStressStrainCurvesGeneric(labelStrs,ssCurvePlotParams,pos,...
        genSSCurves.Q11_accelAvg_Eq12,genSSCurves.Q11_strainAvg_Eq12);
    saveFile = [savePath,'\','GenSS_Curve_Qxx_Eq12'];
    func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)

    labelStrs.t = '$Q_{xx}$ Eq.1+3';
    hf = func_plotStressStrainCurvesGeneric(labelStrs,ssCurvePlotParams,pos,...
        genSSCurves.Q11_accelAvg_Eq13,genSSCurves.Q11_strainAvg_Eq13);
    saveFile = [savePath,'\','GenSS_Curve_Qxx_Eq13'];
    func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)
    
    labelStrs.t = '$Q_{xx}$ Eq.2+3';
    hf = func_plotStressStrainCurvesGeneric(labelStrs,ssCurvePlotParams,pos,...
        genSSCurves.Q11_accelAvg_Eq23,genSSCurves.Q11_strainAvg_Eq23);
    saveFile = [savePath,'\','GenSS_Curve_Qxx_Eq23'];
    func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)
    
    % Qxy
    labelStrs.t = '$Q_{xy}$ Eq.1+2';
    hf = func_plotStressStrainCurvesGeneric(labelStrs,ssCurvePlotParams,pos,...
        genSSCurves.Q12_accelAvg_Eq12,genSSCurves.Q12_strainAvg_Eq12);
    saveFile = [savePath,'\','GenSS_Curve_Qxy_Eq12'];
    func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)
    
    labelStrs.t = '$Q_{xy}$ Eq.1+3';
    hf = func_plotStressStrainCurvesGeneric(labelStrs,ssCurvePlotParams,pos,...
        genSSCurves.Q12_accelAvg_Eq13,genSSCurves.Q12_strainAvg_Eq13);
    saveFile = [savePath,'\','GenSS_Curve_Qxy_Eq13'];
    func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)
    
    labelStrs.t = '$Q_{xy}$ Eq.2+3';
    hf = func_plotStressStrainCurvesGeneric(labelStrs,ssCurvePlotParams,pos,...
        genSSCurves.Q12_accelAvg_Eq23,genSSCurves.Q12_strainAvg_Eq23);
    saveFile = [savePath,'\','GenSS_Curve_Qxy_Eq23'];
    func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)
    
    %----------------------------------------------------------------------
    % PLOT: Generalised SS Curves, Isotropic Case - FITTED CURVES
    %{
    fprintf('\tPlotting stress-strain curves.\n')
    ssCurveIdentPlotOpts.indLRange = 1:length(pos.x);
    ssCurvePlotParams.locXInd = round(length(pos.x)*ssCurvePlotParams.locXPcVec);
    
    labelStrs.y = 'Accel. Avg.';
    labelStrs.x = 'Strain. Avg.';
    
    % Qxx
    labelStrs.t = '$Q_{xx}$ Eq.1+2';
    hf = func_plotStressStrainCurvesGeneric(labelStrs,ssCurvePlotParams,pos,...
        genSSCurves.Q11_accelAvg_Eq12,genSSCurves.Q11_strainAvg_Eq12);
    saveFile = [savePath,'\','GenSS_Curve_Qxx_Eq12'];
    func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)

    labelStrs.t = '$Q_{xx}$ Eq.1+3';
    hf = func_plotStressStrainCurvesGeneric(labelStrs,ssCurvePlotParams,pos,...
        genSSCurves.Q11_accelAvg_Eq13,genSSCurves.Q11_strainAvg_Eq13);
    saveFile = [savePath,'\','GenSS_Curve_Qxx_Eq13'];
    func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)
    
    labelStrs.t = '$Q_{xx}$ Eq.2+3';
    hf = func_plotStressStrainCurvesGeneric(labelStrs,ssCurvePlotParams,pos,...
        genSSCurves.Q11_accelAvg_Eq23,genSSCurves.Q11_strainAvg_Eq23);
    saveFile = [savePath,'\','GenSS_Curve_Qxx_Eq23'];
    func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)
    
    % Qxy
    labelStrs.t = '$Q_{xy}$ Eq.1+2';
    hf = func_plotStressStrainCurvesGeneric(labelStrs,ssCurvePlotParams,pos,...
        genSSCurves.Q12_accelAvg_Eq12,genSSCurves.Q12_strainAvg_Eq12);
    saveFile = [savePath,'\','GenSS_Curve_Qxy_Eq12'];
    func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)
    
    labelStrs.t = '$Q_{xy}$ Eq.1+3';
    hf = func_plotStressStrainCurvesGeneric(labelStrs,ssCurvePlotParams,pos,...
        genSSCurves.Q12_accelAvg_Eq13,genSSCurves.Q12_strainAvg_Eq13);
    saveFile = [savePath,'\','GenSS_Curve_Qxy_Eq13'];
    func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)
    
    labelStrs.t = '$Q_{xy}$ Eq.2+3';
    hf = func_plotStressStrainCurvesGeneric(labelStrs,ssCurvePlotParams,pos,...
        genSSCurves.Q12_accelAvg_Eq23,genSSCurves.Q12_strainAvg_Eq23);
    saveFile = [savePath,'\','GenSS_Curve_Qxy_Eq23'];
    func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)
    %}

end

