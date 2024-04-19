function hf = func_plotStiffnessVsLength(labelStrs,stiffPlotOpts,pos,...
    identStiffVsL,identStiffAvg,identOpts)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 9/2/2018
%
% Plots stiffness vs the length of the specimen for the stress-gauge
% stress-strain curve fitting or the generalised stress-strain curves


% Get some plot parameters for nice formatting
plotProps = func_initPlotPropsStruct(stiffPlotOpts.formatType);

% Check how many slices the identification was performed on, useful for
% plotting the identified stiffness for off-axis samples
xR = 1:length(identStiffVsL);
xMax = max(pos.x(xR))*10^3;
if isfield(stiffPlotOpts,'specifyLRange')
    if stiffPlotOpts.specifyLRange
        xR = stiffPlotOpts.indLRange;
        xMax = max(pos.x)*10^3;
    end
end

% Setup the plotting vectors for the identified average and error bounds
if stiffPlotOpts.plotAvg
    identAvgPlot(1:length(pos.x)) = identStiffAvg;
    identAvgUpperBound(1:length(pos.x)) = identStiffAvg*(1+stiffPlotOpts.plotAvgErrBoundPc);
    identAvgLowerBound(1:length(pos.x)) = identStiffAvg*(1-stiffPlotOpts.plotAvgErrBoundPc);
end

% Setup the plotting vectors for quasi-static reference stiffness and error
% bounds
if stiffPlotOpts.plotQS
    QSRefPlot(1:length(pos.x)) = stiffPlotOpts.QSRef;
    QSRefUpperBound(1:length(pos.x)) = stiffPlotOpts.QSRef*(1+stiffPlotOpts.plotQSErrBoundPc);
    QSRefLowerBound(1:length(pos.x)) = stiffPlotOpts.QSRef*(1-stiffPlotOpts.plotQSErrBoundPc);
end

% Create the axis limits
if strcmp(stiffPlotOpts.axisLimVar,'avg')
    yLim1 = identStiffAvg*stiffPlotOpts.unitConv*...
        (1-stiffPlotOpts.axisRangeFactor*stiffPlotOpts.plotAvgErrBoundPc);
    yLim2 = identStiffAvg*stiffPlotOpts.unitConv*...
        (1+stiffPlotOpts.axisRangeFactor*stiffPlotOpts.plotAvgErrBoundPc);
else
    yLim1 = stiffPlotOpts.QSRef*stiffPlotOpts.unitConv*...
        (1-stiffPlotOpts.axisRangeFactor*stiffPlotOpts.plotQSErrBoundPc);
    yLim2 = stiffPlotOpts.QSRef*stiffPlotOpts.unitConv*...
        (1+stiffPlotOpts.axisRangeFactor*stiffPlotOpts.plotQSErrBoundPc);
end
if yLim2 > yLim1
    axisLims = [0,xMax,yLim1,yLim2];
else
    axisLims = [0,xMax,yLim2,yLim1];
end
if isfield(stiffPlotOpts,'yAxisStartZero')
    if stiffPlotOpts.yAxisStartZero
        axisLims(3) = 0;
    end
end

% Create vectors for vertical lines showing the averaging window
if stiffPlotOpts.plotAvgPosRange
    avgWindowY = linspace(axisLims(3),axisLims(4),10);
    avgWindowStart(1:length(avgWindowY)) = pos.x(identOpts.avgQVsLRange(1));
    avgWindowEnd(1:length(avgWindowY)) = pos.x(identOpts.avgQVsLRange(end));
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
if stiffPlotOpts.plotAvgPosRange
    legendStrs{ll} = ['Avg. Range'];
    ll = ll+1;
end

%--------------------------------------------------------------------------
% Create and size the figure
plotProps.sizePerFigXcm = plotProps.singleColFigFactor*plotProps.sizePerFigXcm;
hf = func_createFigure(plotProps);

%--------------------------------------------------------------------------
% Plot the identified values
hold on
% Plot the identified stiffness as a function of length
plot(pos.x(xR)*10^3,identStiffVsL(xR)*stiffPlotOpts.unitConv,'-xb','linewidth',plotProps.lw,'markersize',plotProps.ms/2)

if stiffPlotOpts.plotAvg
% Plot the average identified value as a horizontal line
plot(pos.x*10^3,identAvgPlot*stiffPlotOpts.unitConv,'--b','linewidth',plotProps.lw,'markersize',plotProps.ms)
end

if stiffPlotOpts.plotAvgErrBound
% Plot the first average error bound
plot(pos.x*10^3,identAvgUpperBound*stiffPlotOpts.unitConv,'-b','linewidth',plotProps.lw,'markersize',plotProps.ms)
end

if stiffPlotOpts.plotQS
% Plot the QS ref. horizontal line
plot(pos.x*10^3,QSRefPlot*stiffPlotOpts.unitConv,'--k','linewidth',plotProps.lw,'markersize',plotProps.ms)
end

if stiffPlotOpts.plotQSErrBound
% Plot the first QS ref. bound
plot(pos.x*10^3,QSRefUpperBound*stiffPlotOpts.unitConv,'-k','linewidth',plotProps.lw,'markersize',plotProps.ms)
end

if stiffPlotOpts.plotAvgPosRange
% Plot the start of the time window for averaging    
plot(avgWindowStart*10^3,avgWindowY,'--r','linewidth',plotProps.lw,'markersize',plotProps.ms)
end

%--------------------------------------------------------------------------
% Final bounding bars plotted last to remove from legend

if stiffPlotOpts.plotQSErrBound
% Plot final bound for QS ref
plot(pos.x*10^3,QSRefLowerBound*stiffPlotOpts.unitConv,'-k','linewidth',plotProps.lw,'markersize',plotProps.ms)
end

if stiffPlotOpts.plotAvgErrBound
% Plot final bound for the average error bounds
plot(pos.x*10^3,identAvgLowerBound*stiffPlotOpts.unitConv,'-b','linewidth',plotProps.lw,'markersize',plotProps.ms)
end

if stiffPlotOpts.plotAvgPosRange
% Plot the start of the time window for averaging    
plot(avgWindowEnd*10^3,avgWindowY,'--r','linewidth',plotProps.lw,'markersize',plotProps.ms)
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
axis(axisLims)
box on
grid on

end

