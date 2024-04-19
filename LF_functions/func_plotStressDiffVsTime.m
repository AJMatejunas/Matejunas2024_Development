function hf = func_plotStressDiffVsTime(plotParams,savePath,time,...
    stress,fracture)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 19/2/2018
% Date Edited: 1/5/2019
%
% Plots the difference between stress measures at the fracture location vs
% time
    
    % Create path for saving the image 
    if exist(savePath,'file') ~= 7
        mkdir(savePath);
    end

    % Create a struct of formatting properties for the figure
    plotProps = func_initPlotPropsStruct(plotParams.formatType);
    plotParams.Rows = 2;
    plotParams.Cols = 2;
    
    % Calculate the width averaged stress for the strain recon stress
    stress.QRecon.xAvg = func_avgFFVarOverWidth(stress.QRecon.x);
    
    % Don't plot the 3 last frames due to numerical diff edge effects
    fRange = 1:size(stress.xAvg,2)-3;
    
    %----------------------------------------------------------------------
    % Create and size the figure
    plotProps.sizePerFigXcm = 1.2*plotProps.sizePerFigXcm;
    hf = func_createFigure(plotProps,plotParams);
    axisScaleFact = 1.2;
    
    %----------------------------------------------------------------------
    % Width Averaged Stress Difference Subplots
    titleStr = '(a)';
    yStr = '$\overline{\sigma_{xx}}^{y}$ [$MPa$]';
    subplot(plotParams.Rows,plotParams.Cols,1)
    plot(stress.xAvg(fracture.locX,fRange)*10^-6,...
        '-ok','linewidth',plotProps.lw,'markersize',plotProps.ms/2)
    hold on
    plot(stress.QRecon.xAvg(fracture.locX,fRange)*10^-6,...
        '-xb','linewidth',plotProps.lw,'markersize',plotProps.ms/2)
    hold off
    th = xlabel('Frame');
    set(th,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
    th = ylabel(yStr);
    set(th,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
    title(titleStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
    xlim([0,time.numFrames])
    ylim([round(axisScaleFact*min(stress.xAvg(fracture.locX,fRange)),2,'significant')...
        ,round(axisScaleFact*max(stress.xAvg(fracture.locX,fRange)),2,'significant')]*10^-6)
    lh = legend({'$\overline{\sigma_{xx}}^{y}(SG)$','$\overline{\sigma_{xx}}^{y}(\epsilon)$'},'location','northwest');
    set(lh,'fontsize',plotProps.lfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
    set(gca, 'fontsize', plotProps.fs,'fontname',plotProps.ft)
    set(gca,'XMinorTick','on','YMinorTick','on')
    box on
    grid on
    
    titleStr = '(c)';
    widthAvgDiff = stress.xAvg(fracture.locX,fRange) - stress.QRecon.xAvg(fracture.locX,fRange);
    yStr = '$\overline{\sigma_{xx}}^{y}(SG)- \overline{\sigma_{xx}}^{y}(\epsilon)$ [$MPa$]';
    subplot(plotParams.Rows,plotParams.Cols,3)
    plot(widthAvgDiff*10^-6,...
        '-or','linewidth',plotProps.lw,'markersize',plotProps.ms/2)
    th = xlabel('Frame');
    set(th,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
    th = ylabel(yStr);
    set(th,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
    title(titleStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
    xlim([0,time.numFrames])
    set(gca, 'fontsize', plotProps.fs,'fontname',plotProps.ft)
    set(gca,'XMinorTick','on','YMinorTick','on')
    box on
    grid on
    
    %----------------------------------------------------------------------
    % Area Averaged Stress Difference Subplots
    titleStr = '(b)';
    yStr = '$\overline{\sigma_{xx}}^{VG}$ [$MPa$]';
    subplot(plotParams.Rows,plotParams.Cols,2)
    plot(stress.virtGaugeAvg.x(fRange)*10^-6,...
        '-ok','linewidth',plotProps.lw,'markersize',plotProps.ms/2)
    hold on
    plot(stress.QRecon.virtGaugeAvg.x(fRange)*10^-6,...
        '-xb','linewidth',plotProps.lw,'markersize',plotProps.ms/2)
    hold off
    th = xlabel('Frame');
    set(th,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
    th = ylabel(yStr);
    set(th,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
    title(titleStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
    xlim([0,time.numFrames])
    ylim([round(axisScaleFact*min(stress.virtGaugeAvg.x(fRange)),2,'significant')...
        ,round(axisScaleFact*max(stress.virtGaugeAvg.x(fRange)),2,'significant')]*10^-6)
    lh = legend({'$\overline{\sigma_{xx}}^{VG}(LSG)$','$\overline{\sigma_{xx}}^{VG}(\epsilon)$'},'location','northwest');
    set(lh,'fontsize',plotProps.lfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
    set(gca, 'fontsize', plotProps.fs,'fontname',plotProps.ft)
    set(gca,'XMinorTick','on','YMinorTick','on')
    box on
    grid on
    
    titleStr = '(d)';
    stressAvgOverGaugeAreaDiff = stress.virtGaugeAvg.x(fRange) - stress.QRecon.virtGaugeAvg.x(fRange);
    yStr = '$\overline{\sigma_{xx}}^{VG}(LSG)- \overline{\sigma_{xx}}^{VG}(\epsilon)$ [$MPa$]';
    subplot(plotParams.Rows,plotParams.Cols,4)
    plot(stressAvgOverGaugeAreaDiff*10^-6,...
        '-or','linewidth',plotProps.lw,'markersize',plotProps.ms/2)
    th = xlabel('Frame','fontsize',plotProps.hfs,'fontname',plotProps.ft);
    set(th,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
    th = ylabel(yStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft);
    set(th,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
    title(titleStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
    xlim([0,time.numFrames])
    set(gca, 'fontsize', plotProps.fs,'fontname',plotProps.ft)
    set(gca,'XMinorTick','on','YMinorTick','on')
    box on
    grid on
     
end

