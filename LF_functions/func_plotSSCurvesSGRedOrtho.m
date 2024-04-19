function func_plotSSCurvesSGRedOrtho(savePath,plotParams,ssCurvePlotParams,...
    pos,stress,strain,identStiffSG)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 28/1/2019
% Date Edited: 5/6/2019
%
% Plots stress-strain curves for the reduced orthotropic case.

    % Set plot properties for figure formatting
    plotProps = func_initPlotPropsStruct(ssCurvePlotParams.formatType);

    % Work out the indices to plot based on the specimen length
    ssCurvePlotParams.locXInd = round(length(pos.x)*ssCurvePlotParams.locXPcVec);
    
    % Full stress-strain curves
    labelStrs.y = 'Stress, $\overline{\sigma_{xx}}$ (MPa)';
    labelStrs.x = 'Strain, $\overline{\epsilon_{xx}}$ (mm/m)';
    ssCurvePlotParams.plotOnlyFittedRegion = false;
    ssCurvePlotParams.plotFittedCurve = false;
    hf = func_plotStressStrainCurvesGeneric(labelStrs,ssCurvePlotParams,...
        pos,stress.xAvg,strain.xAvg);
    saveFile = [savePath,'\','SG_SSCurvesFull'];
    func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)
    
    % Fitted stress-strain curves
    ssCurvePlotParams.plotOnlyFittedRegion = true;
    ssCurvePlotParams.plotFittedCurve = true;
    hf = func_plotStressStrainCurvesGeneric(labelStrs,ssCurvePlotParams,pos,...
        stress.xAvg,strain.xAvg,identStiffSG);
    saveFile = [savePath,'\','SG_SSCurvesFitted'];
    func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)
end

