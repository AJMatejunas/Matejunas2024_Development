function hf = func_plotStiffnessVsTime(labelStrs,stiffPlotOpts,time,...
    identStiffVsT,identStiffAvg,identOpts)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 9/2/2018
% Date Edited: 7/6/2019
%
% Plots stiffness vs time from VFM identification

% Get some plot parameters for nice formatting
plotProps = func_initPlotPropsStruct(stiffPlotOpts.formatType);

% Setup the plotting vectors for the identified average and error bounds
if stiffPlotOpts.plotAvg
    identAvgPlot(1:length(time.vec)) = identStiffAvg;
    identAvgUpperBound(1:length(time.vec)) = identStiffAvg*(1+stiffPlotOpts.plotAvgErrBoundPc);
    identAvgLowerBound(1:length(time.vec)) = identStiffAvg*(1-stiffPlotOpts.plotAvgErrBoundPc);
end

% Setup the plotting vectors for quasi-static reference stiffness and error
% bounds
if stiffPlotOpts.plotQS
    QSRefPlot(1:length(time.vec)) = stiffPlotOpts.QSRef;
    QSRefUpperBound(1:length(time.vec)) = stiffPlotOpts.QSRef*(1+stiffPlotOpts.plotQSErrBoundPc);
    QSRefLowerBound(1:length(time.vec)) = stiffPlotOpts.QSRef*(1-stiffPlotOpts.plotQSErrBoundPc);
end

% Create the axis limits
if strcmp(stiffPlotOpts.axisLimVar,'avg')
     axisLims = [0,time.vec(end)*10^6,...
        identStiffAvg*stiffPlotOpts.unitConv*(1-stiffPlotOpts.axisRangeFactor*stiffPlotOpts.plotAvgErrBoundPc),...
        identStiffAvg*stiffPlotOpts.unitConv*(1+stiffPlotOpts.axisRangeFactor*stiffPlotOpts.plotAvgErrBoundPc)];
else
    axisLims = [0,time.vec(end)*10^6,...
        stiffPlotOpts.QSRef*stiffPlotOpts.unitConv*(1-stiffPlotOpts.axisRangeFactor*stiffPlotOpts.plotQSErrBoundPc),...
        stiffPlotOpts.QSRef*stiffPlotOpts.unitConv*(1+stiffPlotOpts.axisRangeFactor*stiffPlotOpts.plotQSErrBoundPc)];
end

% Plot time bounds over which the average is taken
if stiffPlotOpts.plotAvgTimeRange
    if stiffPlotOpts.specifyAxisLims
        avgWindowY = linspace(stiffPlotOpts.axisLims(3),stiffPlotOpts.axisLims(4),10);
    else
        avgWindowY = linspace(axisLims(3),axisLims(4),10);
    end
    avgWindowStart(1:length(avgWindowY)) = identOpts.avgStartTime;
    avgWindowEnd(1:length(avgWindowY)) = identOpts.avgEndTime;
end

% Setup the legend strings in the correct order
ll = 1;
legendStrs{ll} = 'Ident.';
ll = ll+1;
if stiffPlotOpts.plotAvg
    legendStrs{ll} = 'Ident. Avg.';
    ll = ll+1;
end
if stiffPlotOpts.plotAvgErrBound
    legendStrs{ll} = ['Ident. Avg. $\pm',num2str(stiffPlotOpts.plotAvgErrBoundPc*100),'\%$'];
    ll = ll+1;
end
if stiffPlotOpts.plotQS
    legendStrs{ll} = 'QS Ref.';
    ll = ll+1;
end
if stiffPlotOpts.plotQSErrBound
    legendStrs{ll} = ['QS Ref. $\pm',num2str(stiffPlotOpts.plotQSErrBoundPc*100),'\%$'];
    ll = ll+1;
end
if stiffPlotOpts.plotAvgTimeRange
    legendStrs{ll} = ['Avg. Time'];
    ll = ll+1;
end

%--------------------------------------------------------------------------
% Create and size the figure
plotProps.sizePerFigXcm = plotProps.singleColFigFactor*plotProps.sizePerFigXcm;
hf = func_createFigure(plotProps);

%--------------------------------------------------------------------------
% Plot all required lines for identification and error bounds
hold on
% Plot the identified stiffness as a function of time
plot(time.vec(stiffPlotOpts.tRange)*10^6, identStiffVsT(stiffPlotOpts.tRange)*stiffPlotOpts.unitConv,'-xb','linewidth',plotProps.lw,'markersize',plotProps.ms/2)

if stiffPlotOpts.plotAvg
% Plot the average identified value as a horizontal line
plot(time.vec*10^6,identAvgPlot*stiffPlotOpts.unitConv,'-.b','linewidth',plotProps.lw,'markersize',plotProps.ms)
end

if stiffPlotOpts.plotAvgErrBound
% Plot the first average error bound
plot(time.vec*10^6,identAvgUpperBound*stiffPlotOpts.unitConv,'-b','linewidth',plotProps.lw,'markersize',plotProps.ms)
end

if stiffPlotOpts.plotQS
% Plot the QS ref. horizontal line
plot(time.vec*10^6,QSRefPlot*stiffPlotOpts.unitConv,'--k','linewidth',plotProps.lw,'markersize',plotProps.ms)
end

if stiffPlotOpts.plotQSErrBound
% Plot the first QS ref. bound
plot(time.vec*10^6,QSRefUpperBound*stiffPlotOpts.unitConv,'-k','linewidth',plotProps.lw,'markersize',plotProps.ms)
end

if stiffPlotOpts.plotAvgTimeRange
% Plot the start of the time window for averaging    
plot(avgWindowStart*10^6,avgWindowY,'--r','linewidth',plotProps.lw,'markersize',plotProps.ms)
end

%--------------------------------------------------------------------------
% Final bounding bars plotted last to remove from legend
if stiffPlotOpts.plotQSErrBound
% Plot final bound for QS ref
plot(time.vec*10^6,QSRefLowerBound*stiffPlotOpts.unitConv,'-k','linewidth',plotProps.lw,'markersize',plotProps.ms)
end

if stiffPlotOpts.plotAvgErrBound
% Plot final bound for the average error bounds
plot(time.vec*10^6,identAvgLowerBound*stiffPlotOpts.unitConv,'-b','linewidth',plotProps.lw,'markersize',plotProps.ms)
end

if stiffPlotOpts.plotAvgTimeRange
% Plot the end of the time window for averaging    
plot(avgWindowEnd*10^6,avgWindowY,'--r','linewidth',plotProps.lw,'markersize',plotProps.ms)
end

hold off

%--------------------------------------------------------------------------
% Format axes and labels

% Set the title to the identified average stiffness value
if stiffPlotOpts.printAvgInTitle
    th = title(labelStrs.t);
    set(th,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
end

% Set labels and plot formatting parameters
th = xlabel(labelStrs.x);
set(th,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
th = ylabel(labelStrs.y);
set(th,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
set(gca, 'fontsize', plotProps.fs,'fontname',plotProps.ft)
set(gca,'XMinorTick','on','YMinorTick','on')
lh = legend(legendStrs,'location',plotProps.lLoc);
set(lh,'fontsize',plotProps.lfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
if stiffPlotOpts.specifyAxisLims
    axis(stiffPlotOpts.axisLims)
else
    axis(axisLims)
end
box on
grid on

end

