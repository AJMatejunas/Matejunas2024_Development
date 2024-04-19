function [strStrainRate,hf] = func_plotStrainRateAtFracLoc(plotParams,savePath,time,virtualGauge,...
    fracture,strainRate)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 22/6/2017
%
% Plots the strain rate at the identified fracture location vs time

    % Create path for saving the image 
    if exist(savePath,'file') ~= 7
        mkdir(savePath);
    end

    % Calculate the average strain rates over each region
    strStrainRate.atFracLocYFracLocX = squeeze(strainRate.x(fracture.locY,fracture.locX,:));
    strStrainRate.avgOverYFracLocX = squeeze(mean(strainRate.x(:,fracture.locX,:)));
    strStrainRate.avgGaugeYGaugeX = squeeze(mean(mean(strainRate.x(virtualGauge.yRange,virtualGauge.xRange,:))));
    
    strStrainRate.avgGaugeXYFracTime = strStrainRate.avgGaugeYGaugeX(fracture.strengthFrame);
    
    SRR = 1:fracture.strengthFrame;
    maxSR = max([max(strStrainRate.atFracLocYFracLocX(SRR)),max(strStrainRate.avgOverYFracLocX(SRR)),max(strStrainRate.avgGaugeYGaugeX(SRR))]);
    minSR = min([min(strStrainRate.atFracLocYFracLocX(SRR)),min(strStrainRate.avgOverYFracLocX(SRR)),min(strStrainRate.avgGaugeYGaugeX(SRR))]);
    SRAxisMax = round(maxSR,2,'significant');
    SRAxisMin = round(minSR,2,'significant');
    
    % Create a line for the fracture location
    fracLineY = linspace(SRAxisMin,SRAxisMax,10);
    fracLineT(1:length(fracLineY)) = time.vec(fracture.strengthFrame)*10^6;
    fracLineF(1:length(fracLineY)) = fracture.strengthFrame;
    
    % Create a struct of formatting properties for the figure
    plotProps = func_initPlotPropsStruct(plotParams.formatType);
    plotParams.Rows = 1;
    plotParams.Cols = 2;
    
    %----------------------------------------------------------------------
    % Create and size the figure
    plotProps.sizePerFigXcm = 1.2*plotProps.sizePerFigXcm;
    hf = func_createFigure(plotProps,plotParams);
    
    %----------------------------------------------------------------------
    % Strain Rate @ Frac Loc vs Time
    tStr = ['$\dot{\epsilon_{xx}}$ $(t=',sprintf('%.1f',time.vec(fracture.strengthFrame)*10^6)...
        ,') = ',sprintf('%.0f',strStrainRate.avgGaugeYGaugeX(fracture.strengthFrame)),'$ $s^{-1}$'];
    yStr = '$\dot{\epsilon_{xx}}$ $s^{-1}$';
    xStr = 'Time $t$ ($\mu s$)';
    subplot(plotParams.Rows,plotParams.Cols,1)
    plot(time.vec*10^6,strStrainRate.atFracLocYFracLocX ,...
        '-ok','linewidth',plotProps.lw,'markersize',plotProps.ms/2)
    hold on
    plot(time.vec*10^6,strStrainRate.avgOverYFracLocX,...
        '-xb','linewidth',plotProps.lw,'markersize',plotProps.ms/2)
    plot(time.vec*10^6,strStrainRate.avgGaugeYGaugeX,...
        '-+r','linewidth',plotProps.lw,'markersize',plotProps.ms/2)
    plot(fracLineT,fracLineY,'--k','linewidth',plotProps.lw)
    hold off
    title(tStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText)
    xlabel(xStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText)
    ylabel(yStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText)
    axis([0,time.vec(end)*10^6,SRAxisMin,SRAxisMax])
    legend({'$\dot{\epsilon}(x,y)$','$\overline{\dot{\epsilon}}^{y}$','$\overline{\dot{\epsilon}}^{A}$','Frac.'}...
        ,'location','best','Interpreter',plotProps.interpText)
    set(gca, 'fontsize', plotProps.fs,'fontname',plotProps.ft)
    set(gca,'XMinorTick','on','YMinorTick','on')
    box on
    grid on
    
    %----------------------------------------------------------------------
    % Strain Rate @ Frac Loc vs Frame
    tStr = ['$\dot{\epsilon_{xx}}$ $(f=',sprintf('%i',fracture.strengthFrame)...
        ,') = ',sprintf('%.0f',strStrainRate.avgGaugeYGaugeX(fracture.strengthFrame)),'$ $s^{-1}$'];
   yStr = '$\dot{\epsilon_{xx}}$ $s^{-1}$';
    xStr = 'Frame';
    subplot(plotParams.Rows,plotParams.Cols,2)
    plot(strStrainRate.atFracLocYFracLocX ,...
        '-ok','linewidth',plotProps.lw,'markersize',plotProps.ms/2)
    hold on
    plot(strStrainRate.avgOverYFracLocX,...
        '-xb','linewidth',plotProps.lw,'markersize',plotProps.ms/2)
    plot(strStrainRate.avgGaugeYGaugeX,...
        '-+r','linewidth',plotProps.lw,'markersize',plotProps.ms/2)
    plot(fracLineF,fracLineY,'--k','linewidth',plotProps.lw)
    hold off
    title(tStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText)
    xlabel(xStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText)
    ylabel(yStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText)
    axis([0,time.numFrames,SRAxisMin,SRAxisMax])
    legend({'$\dot{\epsilon}(x,y)$','$\overline{\dot{\epsilon}}^{y}$','$\overline{\dot{\epsilon}}^{A}$','Frac.'}...
        ,'location','best','Interpreter',plotProps.interpText)
    set(gca, 'fontsize', plotProps.fs,'fontname',plotProps.ft)
    set(gca,'XMinorTick','on','YMinorTick','on')
    box on
    grid on

    saveFile = [savePath,'\','StrIdent_StrainRate'];
    func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)
end

