function hf = func_plotFieldAvgVsTime(plotParams,time,pos,disp,accel,strain,strainRate)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 19/6/2017
% Plots field averaged variables vs time

    % Load the plot properties structure
    plotProps = func_initPlotPropsStruct(plotParams.formatType);
    
    % Setup and plot the figure
    hf = func_createFigure(plotProps,plotParams);
    
    dispPlot = squeeze(nanmean(nanmean(disp.x)))*10^3;
    accelPlot = squeeze(nanmean(nanmean(accel.x)));
    strainPlot = squeeze(nanmean(nanmean(strain.x)))*10^3;
    strainRatePlot = squeeze(nanmean(nanmean(strainRate.x)));
    
    xStr = 'Time $t$ [$\mu s$]';
    %----------------------------------------------------------------------
    % Subplot 1: Displacement
    yStr = 'Disp. Avg. $\overline{\delta_{x}}^{S}$ [$mm$]';
    subplot(plotParams.Rows,plotParams.Cols,1);
    plot(time.vec*10^6,dispPlot,'-ob',...
        'linewidth',plotProps.lw,'markersize',plotProps.ms/2)
    xlabel(xStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText)
    ylabel(yStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText)
    set(gca, 'fontsize', plotProps.fs,'fontname',plotProps.ft)
    set(gca,'XMinorTick','on','YMinorTick','on')
    box on
    
    %----------------------------------------------------------------------
    % Subplot 2: Acceleration
    yStr = 'Accel. Avg. $\overline{a_{x}}^{S}$ [$m.s^{-2}$]';
    subplot(plotParams.Rows,plotParams.Cols,2);
    plot(time.vec*10^6,accelPlot,'-ob',...
        'linewidth',plotProps.lw,'markersize',plotProps.ms/2)
    xlabel(xStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText)
    ylabel(yStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText)
    set(gca, 'fontsize', plotProps.fs,'fontname',plotProps.ft)
    set(gca,'XMinorTick','on','YMinorTick','on')
    box on
    
    %----------------------------------------------------------------------
    % Subplot 3: Strain 
    yStr = 'Strain Avg. $\overline{\epsilon_{xx}}^{S}$, [$mm.m^{-1}$]';
    subplot(plotParams.Rows,plotParams.Cols,3);
    plot(time.vec*10^6,strainPlot,'-ob',...
        'linewidth',plotProps.lw,'markersize',plotProps.ms/2)
    xlabel(xStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText)
    ylabel(yStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText)
    set(gca, 'fontsize', plotProps.fs,'fontname',plotProps.ft)
    set(gca,'XMinorTick','on','YMinorTick','on')
    box on
    
    %----------------------------------------------------------------------
    % Subplot 4: Strain Rate
    yStr = 'Strain Rate Avg. $\overline{\epsilon_{xx}/dt}^{S}$, [$s^{-1}$]';
    subplot(plotParams.Rows,plotParams.Cols,4);
    plot(time.vec*10^6,strainRatePlot,'-ob',...
        'linewidth',plotProps.lw,'markersize',plotProps.ms/2)
    xlabel(xStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText)
    ylabel(yStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText)
    set(gca, 'fontsize', plotProps.fs,'fontname',plotProps.ft)
    set(gca,'XMinorTick','on','YMinorTick','on')
    box on

end

