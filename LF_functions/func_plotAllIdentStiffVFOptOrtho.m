function func_plotAllIdentStiffVFOptOrtho(savePath,plotParams,time,material,...
    VFOptPlotOpts,identStiffVFOpt,VFOpts)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 28/1/2019
% Date Edited: 28/1/2019
%
% Plots all identified stiffness components against time for the
% orthotropic case using optimised virtual fields.

    % Set plot properties for figure formatting
    plotProps = func_initPlotPropsStruct(VFOptPlotOpts.formatType);
    labelStrs.x = 'Time $t$, ($\mu s$)';
    
    % Plot identified Q11 vs Time
    labelStrs.y = 'Stiffness $Q_{11}$, (GPa)';
    VFOptPlotOpts.targVal = material.Q11*10^-9;
    hf = func_plotStiffnessVsTime(labelStrs,VFOptPlotOpts,time,...
        identStiffVFOpt.Q11VsT,median(identStiffVFOpt.Q11VsT));
    saveFile = [savePath,'\','VFOpt_StiffVsT_Q11'];
    print(hf,saveFile,plotProps.format,plotProps.saveRes)
    saveas(hf,saveFile,'fig')

    % Plot identified Q22 vs Time
    labelStrs.y = 'Stiffness $Q_{22}$, (GPa)';
    VFOptPlotOpts.targVal = material.Q22*10^-9;
    hf = func_plotStiffnessVsTime(labelStrs,VFOptPlotOpts,time,...
        identStiffVFOpt.Q22VsT,median(identStiffVFOpt.Q22VsT));
    saveFile = [savePath,'\','VFOpt_StiffVsT_Q22'];
    print(hf,saveFile,'-dpng','-r0')
    saveas(hf,saveFile,'fig')

    % Plot identified Q12 vs Time
    labelStrs.y = 'Stiffness $Q_{12}$, (GPa)';
    VFOptPlotOpts.targVal = material.Q12*10^-9;
    hf = func_plotStiffnessVsTime(labelStrs,VFOptPlotOpts,time,...
        identStiffVFOpt.Q12VsT,median(identStiffVFOpt.Q12VsT));
    saveFile = [savePath,'\','VFOpt_StiffVsT_Q12'];
    print(hf,saveFile,plotProps.format,plotProps.saveRes)
    saveas(hf,saveFile,'fig')

    % Plot identified Q66 vs Time
    labelStrs.y = 'Stiffness $Q_{66}$, (GPa)';
    VFOptPlotOpts.targVal = material.Q66*10^-9;
    hf = func_plotStiffnessVsTime(labelStrs,VFOptPlotOpts,time,...
        identStiffVFOpt.Q66VsT,median(identStiffVFOpt.Q66VsT));
    saveFile = [savePath,'\','VFOpt_StiffVsT_Q66'];
    print(hf,saveFile,plotProps.format,plotProps.saveRes)
    saveas(hf,saveFile,'fig')  

end

