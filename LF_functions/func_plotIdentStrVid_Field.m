function hf = func_plotIdentStrVid_Field(plotParams,savePath,pos,time,...
    stress,ssCurveStress,ssCurveStrain,ssCurveLabels,fracture,strOpts,...
    plotFracFrameOnly,plotSSG)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 19/2/2018
% Date Edited: 1/4/2019
%
% Plots stress fields calculated from strain and the stress fields from the
% linear stress gauage. Shows the virtual gauge area and plots the stress 
% divergence and local stress-strain curve.

if nargin < 12
    % Flag for plotting the standard stress gauge as well as the linear gauge
    plotSSG = false; 
end

%----------------------------------------------------------------------
% Initialise Variables for Plotting
virtualGauge = strOpts.virtualGauge;

% Set the frame range for plotting to only have a few frames after fracture
if fracture.strengthFrame < (time.numFrames-7)
    fRange = 1:fracture.strengthFrame+7;
else
    fRange = 1:time.numFrames;
end

% Create cropped position vectors for plotting
plotRangeX = plotParams.cutPxX+1:length(pos.x)-plotParams.cutPxX;
plotRangeY = plotParams.cutPxY+1:length(pos.y)-plotParams.cutPxY;
plotXPos = pos.x(plotParams.cutPxX+1:end-plotParams.cutPxX)*10^3;
plotYPos = pos.y(plotParams.cutPxY+1:end-plotParams.cutPxY)*10^3;

% Create the variables for the virtual gauge rectangle
XPos = pos.x(virtualGauge.xRange(1))*10^3;
YPos = pos.y(virtualGauge.yRange(1))*10^3;
Width = (pos.x(virtualGauge.xRange(end))-pos.x(virtualGauge.xRange(1)))*10^3;
Height = (pos.y(virtualGauge.yRange(end))-pos.y(virtualGauge.yRange(1)))*10^3;
gaugeRectangle = [XPos,YPos,Width,Height];

% Create a struct of formatting properties for the figure
plotProps = func_initPlotPropsStruct(plotParams.formatType);
plotParams.Cols = 2;
plotParams.Rows = 2;

% Calculate the stress difference based on local gauge over the virtual
% gauge area
stressAvgFromStrainOverGaugeArea = squeeze(nanmean(nanmean(...
    stress.QRecon.x(virtualGauge.yRange,virtualGauge.xRange,:))));   
stressAvgFromLinSGOverGaugeArea = squeeze(nanmean(nanmean(...
    stress.xLinearGauge(virtualGauge.yRange,virtualGauge.xRange,:))));

% Create variables for setting colour bar axis
if strcmp(plotParams.cAxisType,'Specified')
    plotParams.cRange{1} = plotParams.cAxisStress;
    plotParams.cRange{2} = plotParams.cAxisStress;
else
    plotParams.cRange{1} = func_calcColourBarRange(plotParams.cAxisType,...
        stress.QRecon.x(plotRangeY,plotRangeX,:))*10^-6;
    plotParams.cRange{2} = func_calcColourBarRange(plotParams.cAxisType,...
        stress.QRecon.x(plotRangeY,plotRangeX,:))*10^-6;
end

%----------------------------------------------------------------------
% Create and size the figure
plotProps.sizePerFigXcm = 1.2*plotProps.sizePerFigXcm;
hf = func_createFigure(plotProps,plotParams);
axisScaleFact = 1.2;

% Check if only the fracture frame will be plotted
if plotFracFrameOnly
    fStart = fracture.strengthFrame;
    fEnd = fracture.strengthFrame; 
else
    fStart = 1;
    fEnd = fRange(end);
end

for ff = fStart:fEnd
    %--------------------------------------------------------------------------
    % 1) Stress Field from Strains
    subplot(plotParams.Rows,plotParams.Cols,1)
    titleStr = {['(a) ','$\sigma_{xx} (\epsilon)$',...
        ', $\overline{\sigma_{xx}}^{VG}(\epsilon)$ = ',sprintf('%.0f',stress.QRecon.virtGaugeAvg.x(ff)*10^-6),' $MPa$'],...
        [sprintf('t = %.2f',time.vec(ff)*10^6),'$\mu s$']};
    imagesc(plotXPos,plotYPos,stress.QRecon.x(plotRangeY,plotRangeX,ff)*10^-6)
    hold on
        rectangle('Position',gaugeRectangle,'linewidth',plotProps.lw,'linestyle','-')
    hold off
    xlabel('X [mm]','fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
    ylabel('Y [mm]','fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
    title(titleStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText,'horizontalAlignment', plotProps.titleAlign);
    set(gca,'fontsize',plotProps.fs,'fontname',plotProps.ft,'linewidth',plotProps.lw)
    colorbar
    colormap(jet)
    caxis(plotParams.cRange{1})
    axis image
    
    %--------------------------------------------------------------------------
    % 2) Stress Field from Linear Stress Gauge
    subplot(plotParams.Rows,plotParams.Cols,2)
    titleStr = {['(b) ','$\sigma_{xx} (LSG)$',...
        ', $\overline{\sigma_{xx}}^{VG}$ = ',sprintf('%.0f',stress.virtGaugeAvg.x(ff)*10^-6),' $MPa$'],...
        [sprintf('t = %.2f',time.vec(ff)*10^6),'$\mu s$']};
    imagesc(plotXPos,plotYPos,stress.xLinearGauge(plotRangeY,plotRangeX,ff)*10^-6)
    hold on
        rectangle('Position',gaugeRectangle,'linewidth',plotProps.lw,'linestyle','-')
    hold off
    xlabel('X [mm]','fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
    ylabel('Y [mm]','fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
    title(titleStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText,'horizontalAlignment', plotProps.titleAlign);
    set(gca,'fontsize',plotProps.fs,'fontname',plotProps.ft,'linewidth',plotProps.lw)
    colorbar
    colormap(jet)
    caxis(plotParams.cRange{2})
    axis image

    %--------------------------------------------------------------------------
    % 3) Plot the average stress in the guage area as a function of time
    subplot(plotParams.Rows,plotParams.Cols,3)
    titleStr = ['(c) ',sprintf('t = %.2f',time.vec(ff)*10^6),'$\mu s$'];
    yStr = '$\overline{\sigma_{xx}}^{VG}$ [MPa]';
    xStr = 'Time $[\mu s]$';
    
    hold on
    if plotSSG
        plot(time.vec*10^6,stress.xAvg(fracture.locX,:)*10^-6,...
            '-k','linewidth',plotProps.lw,'markersize',plotProps.ms/2)
    end
    plot(time.vec*10^6,stressAvgFromLinSGOverGaugeArea*10^-6,...
        '--ok','linewidth',plotProps.lw,'markersize',plotProps.ms/2)
    plot(time.vec*10^6,stressAvgFromStrainOverGaugeArea*10^-6,...
        '--xb','linewidth',plotProps.lw,'markersize',plotProps.ms/2)
    plot(time.vec(ff)*10^6,stressAvgFromLinSGOverGaugeArea(ff)*10^-6,'or',...
        'linewidth',plotProps.lw*2,'markersize',plotProps.ms*1.5)
    hold off
    
    th = title(titleStr);
    set(th,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'interpreter',plotProps.interpText,'horizontalAlignment', plotProps.titleAlign);
    th = xlabel(xStr);
    set(th,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'interpreter',plotProps.interpText);
    th = ylabel(yStr);
    set(th,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'interpreter',plotProps.interpText);
    if plotSSG
        lh = legend({'$\overline{\sigma_{xx}}^{y}$','$\overline{\sigma_{xx}}^{VG} (LSG)$',...
            '$\overline{\sigma_{xx}}^{VG} (\epsilon)$'},'location','northwest');
    else
        lh = legend({'$\overline{\sigma_{xx}}^{VG} (LSG)$',...
            '$\overline{\sigma_{xx}}^{VG} (\epsilon)$'},'location','northwest');
    end
    set(lh,'fontsize',plotProps.lfs,'fontname',plotProps.ft,'interpreter',plotProps.interpText);
    set(gca, 'fontsize', plotProps.fs,'fontname',plotProps.ft)
    set(gca,'XMinorTick','on','YMinorTick','on')
    
    xlim([0,time.vec(end)*10^6])
    yLimStress = max(abs(round(axisScaleFact*min(stress.virtGaugeAvg.x),2,'significant')),...
        abs(round(axisScaleFact*max(stress.virtGaugeAvg.x),2,'significant')));
    ylim([-yLimStress*10^-6,yLimStress*10^-6]);
    box on
    grid on
    
    %--------------------------------------------------------------------------
    % 4) Plot the local stress-strain curve over the gauge area
    subplot(plotParams.Rows,plotParams.Cols,4)
    titleStr = ['(d) ',sprintf('t = %.2f',time.vec(ff)*10^6),'$\mu s$'];
    
    plot(ssCurveStrain(fRange)*10^3,ssCurveStress(fRange)*10^-6,...
        '-+b','linewidth',plotProps.lw,'markersize',plotProps.ms/2)
    hold on
    if plotSSG
        plot(ssCurveStrain(fRange)*10^3,stress.xAvg(fracture.locX,fRange)*10^-6,...
            '-ok','linewidth',plotProps.lw,'markersize',plotProps.ms/2)
    end
    plot(ssCurveStrain(ff)*10^3,ssCurveStress(ff)*10^-6,'or',...
        'linewidth',plotProps.lw*2,'markersize',plotProps.ms*1.5)
    hold off
    
    th = title(titleStr);
    set(th,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'interpreter',plotProps.interpText,'horizontalAlignment', plotProps.titleAlign);
    th = xlabel(ssCurveLabels.xStr);
    set(th,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'interpreter',plotProps.interpText);
    th = ylabel(ssCurveLabels.yStr);
    set(th,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'interpreter',plotProps.interpText);
    set(gca, 'fontsize', plotProps.fs,'fontname',plotProps.ft)
    set(gca,'XMinorTick','on','YMinorTick','on')
    box on
    grid on

    %--------------------------------------------------------------------------
    % Save this frame to file
    if ~plotFracFrameOnly
        saveFile = [savePath,'\Frame_',num2str(ff)];
        func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)
        clf(hf);
    end
end
end


