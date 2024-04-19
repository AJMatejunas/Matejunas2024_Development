function hf = func_plotAllSpecQxxVsT(QxxPlotOpts,time,material,identStiffnesses)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 7/9/2017

% Get some plot parameters for nice formatting
if isfield(QxxPlotOpts,'format')
    plotProps = func_initPlotPropsStruct(QxxPlotOpts.format);
else
    plotProps = func_initPlotPropsStruct('presentation');
end
plotProps.ms = plotProps.ms/2;

% How many specimens are we plotting?
numSpecs = length(identStiffnesses);
plotTime = time.vec*10^6;   % Conv to micro s 

% Create the line style vectors based on the number of specimens
s = 1;
for lm = 1:min([length(plotProps.marker),length(plotProps.line)])    
    for c = 1:length(plotProps.colour)
        if s > numSpecs
            break;
        end
        LS{s} = [plotProps.line{lm},plotProps.marker{lm},plotProps.colour{c}];
        avgLS{s} = [plotProps.line{lm},plotProps.colour{c}];
        s = s+1;
    end
end

% Select the target value based on the method chosen
if QxxPlotOpts.method == 2
    labelStrs.y = 'Stiffness Q_{xx} (GPa)';
    targVal = material.Qxx*10^-9;
    for s = 1:numSpecs
        plotStiffness{s} = identStiffnesses{s}.QxxVsT*10^-9;
        identAvg{s}(1:length(time.vec)) = identStiffnesses{s}.QxxAvgOverT*10^-9;
    end
else
    labelStrs.y = 'Elastic Modulus E_{xx} (GPa)';
    targVal = material.Exx*10^-9;
    for s = 1:numSpecs
        plotStiffness{s} = identStiffnesses{s}.ExxVsT*10^-9;
        identAvg{s}(1:length(time.vec)) = identStiffnesses{s}.ExxAvgOverT*10^-9;
    end
end

% Setup the plotting vectors for the target and the error bounds
target(1:length(time.vec)) = targVal;
upperLim(1:length(time.vec)) = targVal*(1+QxxPlotOpts.targPc);
lowerLim(1:length(time.vec)) = targVal*(1-QxxPlotOpts.targPc);
if isfield(QxxPlotOpts,'axisType')
    if strcmp(QxxPlotOpts.axisType,'Specified')
        axisLims = [0,max(time.vec(max(QxxPlotOpts.fRange)))*10^6,...
            min(QxxPlotOpts.QxxRange)*10^-9,max(QxxPlotOpts.QxxRange)*10^-9];
    elseif strcmp(QxxPlotOpts.axisType,'StartZero')
        axisLims = [0,max(time.vec(max(QxxPlotOpts.fRange)))*10^6,...
            0,targVal*(1+QxxPlotOpts.rangePc)];    
    else
        axisLims = [0,max(time.vec(max(QxxPlotOpts.fRange)))*10^6,...
            targVal*(1-QxxPlotOpts.rangePc),targVal*(1+QxxPlotOpts.rangePc)];
    end
else
    axisLims = [0,max(time.vec(max(QxxPlotOpts.fRange)))*10^6,...
        targVal*(1-QxxPlotOpts.rangePc),targVal*(1+QxxPlotOpts.rangePc)];
end

% Setup the legend strings
for s = 1:numSpecs
    legendStrs{s} = ['S',num2str(s)];
end
legendStrs{numSpecs+1} = 'QS Ref';
legendStrs{numSpecs+2} = ['QS Ref \pm',num2str(QxxPlotOpts.targPc*100),'%'];

%--------------------------------------------------------------------------
% Create the figure
hf = figure;
set(hf,'Position', [plotProps.locX,plotProps.locY,...
        2*plotProps.sizePerFigX,1.6*plotProps.sizePerFigY])
set(hf,'PaperPositionMode','auto')
labelStrs.x = 'Time (\mus)';

%--------------------------------------------------------------------------
% PLOT - stiffness vs time for each spec
hold on
for s = 1:numSpecs
    plot(plotTime(QxxPlotOpts.fRange), plotStiffness{s}(QxxPlotOpts.fRange),...
        LS{s},'linewidth',plotProps.lw,'markersize',plotProps.ms)
end

% PLOT - QS ref
plot(plotTime,target,'--k','linewidth',plotProps.lw,'markersize',plotProps.ms)

% PLOT - QS ref bounds
plot(plotTime,upperLim,'-k',time.vec*10^6,lowerLim,'-k','linewidth',plotProps.lw,'markersize',plotProps.ms)

% PLOT - identified average for each spec
if QxxPlotOpts.plotAverage
for s = 1:numSpecs
    plot(plotTime,identAvg{s},avgLS{s},'linewidth',plotProps.lw,'markersize',plotProps.ms)
end
end
hold off
%--------------------------------------------------------------------------
xlabel(labelStrs.x,'fontsize',plotProps.hfs,'fontname',plotProps.ft)
ylabel(labelStrs.y,'fontsize',plotProps.hfs,'fontname',plotProps.ft)
set(gca, 'fontsize', plotProps.fs,'fontname',plotProps.ft)
set(gca,'XMinorTick','on','YMinorTick','on')
legend(legendStrs,'location','eastoutside')
axis(axisLims)
box on
grid on

end

