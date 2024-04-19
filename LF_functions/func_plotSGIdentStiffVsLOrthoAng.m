function func_plotSGIdentStiffVsLOrthoAng(savePath,plotParams,...
    ssCurveIdentPlotOpts,pos,material,identStiffSG,...
    angledSlice1,angledSlice2,stressGaugeOpts)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 28/1/2019
% Date Edited: 27/6/2019
%
% Plots identified stiffness against length using the standard stress-gauge

    % Set plot properties for figure formatting
    plotProps = func_initPlotPropsStruct(ssCurveIdentPlotOpts.formatType);
    ssCurveIdentPlotOpts.specifyLRange = false;
    ssCurveIdentPlotOpts.plotAvgPosRange = false; % turn this off for now, fix later.
    unitConv = ssCurveIdentPlotOpts.unitConv; 
    
    if angledSlice1.xMaxInd > 0 
        labelStrs.x = 'Position, $x$ (mm)';
        
        labelStrs.t = ['$E_{22,avg} =$ ',...
            sprintf('%3.1f',identStiffSG.Q22AvgOverL*unitConv),'$~GPa$'];
        ssCurveIdentPlotOpts.QSRef = material.Eyy;
        labelStrs.y = '$E_{22}$ (GPa)';
        hf = func_plotStiffnessVsLength(labelStrs,ssCurveIdentPlotOpts,pos,...
            identStiffSG.Q22VsL,identStiffSG.Q22AvgOverL,stressGaugeOpts);
        saveFile = [savePath,'\','SG_StiffVsL_E22'];
        func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)
        
        labelStrs.t = ['$G_{12,avg} =$ ',...
            sprintf('%3.1f',identStiffSG.Q66AvgOverL_Slice1*unitConv),'$~GPa$'];
        ssCurveIdentPlotOpts.QSRef = material.Gxy;
        labelStrs.y = '$G_{12}$ (GPa)';
        hf = func_plotStiffnessVsLength(labelStrs,ssCurveIdentPlotOpts,pos,...
            identStiffSG.Q66VsL_Slice1,identStiffSG.Q66AvgOverL_Slice1,stressGaugeOpts);
        saveFile = [savePath,'\','SG_StiffVsL_G12_Slice1'];
        func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)
    end
    
    if angledSlice2.xMaxInd > 0
        labelStrs.x = 'Position, $x$ (mm)';
        
        labelStrs.t = ['$E_{11,avg} =$ ',...
            sprintf('%3.1f',identStiffSG.Q11AvgOverL*unitConv),'$~GPa$'];
        ssCurveIdentPlotOpts.QSRef = material.Exx;
        labelStrs.y = '$E_{11}$ (GPa)';
        hf = func_plotStiffnessVsLength(labelStrs,ssCurveIdentPlotOpts,pos,...
            identStiffSG.Q11VsL,identStiffSG.Q11AvgOverL,stressGaugeOpts);
        saveFile = [savePath,'\','SG_StiffVsL_E11'];
        func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)
        
        labelStrs.t = ['$G_{12,avg} =$ ',...
            sprintf('%3.1f',identStiffSG.Q66AvgOverL_Slice1*unitConv),'$~GPa$'];
        ssCurveIdentPlotOpts.QSRef  = material.Gxy;
        labelStrs.y = '$G_{12}$ (GPa)';
        hf = func_plotStiffnessVsLength(labelStrs,ssCurveIdentPlotOpts,pos,...
            identStiffSG.Q66VsL_Slice2,identStiffSG.Q66AvgOverL_Slice2,stressGaugeOpts);
        saveFile = [savePath,'\','SG_StiffVsL_G12_Slice2'];
        func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)
    end
end

