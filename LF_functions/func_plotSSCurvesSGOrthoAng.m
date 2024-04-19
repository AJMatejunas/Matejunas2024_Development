function func_plotSSCurvesSGOrthoAng(savePath,plotParams,ssCurvePlotParams,material,...
    angledSlice1,angledSlice2,pos,stressSlice1,strainSlice1,stressSlice2,strainSlice2)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 28/1/2019
% Date Edited: 27/6/2019
%
% Plots stress-strain curves for the orthotropic angled case.

    % Set plot properties for figure formatting
    plotProps = func_initPlotPropsStruct(ssCurvePlotParams.formatType);
    
    if angledSlice1.xMaxInd > 0 
        % Work out which slice to plot based on the number of slices
        ssCurvePlotParams.locXInd = round(angledSlice1.xMaxInd*ssCurvePlotParams.locXPcVec);
        ssCurvePlotParams.plotOnlyFittedRegion = false;
        ssCurvePlotParams.plotFittedCurve = false;

        labelStrs.y = 'Stress, $\overline{\sigma_{22}}$ (MPa)';
        labelStrs.x = 'Strain, $\overline{\epsilon_{22}}$ (mm/m)';
        hf = func_plotStressStrainCurvesGeneric(labelStrs,ssCurvePlotParams,pos,...
            stressSlice1.avg22,strainSlice1.avg22);
        saveFile = [savePath,'\','SG_SSCurves_Q22'];
        func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)

        labelStrs.y = 'Stress, $\overline{\sigma_{12}}$ (MPa)';
        labelStrs.x = 'Strain, $\overline{\gamma_{12}}$ (mm/m)';
        hf = func_plotStressStrainCurvesGeneric(labelStrs,ssCurvePlotParams,pos,...
            stressSlice1.avg12,strainSlice1.avg12);
        saveFile = [savePath,'\','SG_SSCurves_Q66_Slice1'];
        func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)
    end
    
    if angledSlice2.xMaxInd > 0
        % Work out which slice to plot based on the number of slices
        ssCurvePlotParams.locXInd = round(angledSlice2.xMaxInd*ssCurvePlotParams.locXPcVec);
        ssCurvePlotParams.plotOnlyFittedRegion = false;
        ssCurvePlotParams.plotFittedCurve = false;

        labelStrs.y = 'Stress, $\overline{\sigma_{11}}$ (MPa)';
        labelStrs.x = 'Strain, $\overline{\epsilon_{11}}$ (mm/m)';
        hf = func_plotStressStrainCurvesGeneric(labelStrs,ssCurvePlotParams,pos,...
            stressSlice2.avg11,strainSlice2.avg11);
        saveFile = [savePath,'\','SG_SSCurves_Q11'];
        func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)

        labelStrs.y = 'Stress, $\overline{\sigma_{12}}$ (MPa)';
        labelStrs.x = 'Strain, $\overline{\gamma_{12}}$ (mm/m)';
        hf = func_plotStressStrainCurvesGeneric(labelStrs,ssCurvePlotParams,pos,...
            stressSlice2.avg12,strainSlice2.avg12);
        saveFile = [savePath,'\','SG_SSCurves_Q66_Slice2'];
        func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)
    end

end

