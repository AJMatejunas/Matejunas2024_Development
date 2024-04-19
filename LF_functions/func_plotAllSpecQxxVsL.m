function hf = func_plotAllSpecQxxVsL(QxxPlotOpts,pos,material,identStiffnesses)
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

% Create the line style vectors based on the number of specimens
plotPosMax = 0;
s = 1;
for lm = 1:min([length(plotProps.marker),length(plotProps.line)])    
    for c = 1:length(plotProps.colour)
        if s > numSpecs
            break;
        end
        LS{s} = [plotProps.line{lm},plotProps.marker{lm},plotProps.colour{c}];
        avgLS{s} = [plotProps.line{lm},plotProps.colour{c}];
        
        % Get the position vector for each specimen for plotting
        plotPos{s} = pos{s}.x*10^3;   % Convert to mm
        if max(plotPos{s}) > plotPosMax
            plotPosMax = max(plotPos{s});
        end
        
        % Calculate the
        [~,minIndex] = min(abs(plotPos{s}-(QxxPlotOpts.pRangePc(1)*max(plotPos{s}))));
        [~,maxIndex] = min(abs(plotPos{s}-(QxxPlotOpts.pRangePc(2)*max(plotPos{s}))));
        posRange{s} = minIndex:maxIndex;
        
        % Increment
        s = s+1;
    end
end 
plotPosMaxRng = linspace(0,plotPosMax,10);

% Select the target value based on the method chosen
if QxxPlotOpts.method == 2
    labelStrs.y = 'Stiffness Q_{xx} (GPa)';
    targVal = material.Qxx*10^-9;
    for s = 1:numSpecs
        plotIdentAvg{s}(1:length(plotPosMaxRng)) = identStiffnesses{s}.QxxAvgOverL*10^-9;
        plotIdentStiff{s} = identStiffnesses{s}.QxxVsL*10^-9;
    end
else
    labelStrs.y = 'Elastic Modulus E_{xx} (GPa)';
    targVal = material.Exx*10^-9;
    for s = 1:numSpecs
        plotIdentAvg{s}(1:length(plotPosMaxRng)) = identStiffnesses{s}.ExxAvgOverL*10^-9;
        plotIdentStiff{s} = identStiffnesses{s}.ExxVsL*10^-9;
    end
end

% Setup the plotting vectors for the target and the error bounds
target(1:length(plotPosMaxRng)) = targVal;
upperLim(1:length(plotPosMaxRng)) = targVal*(1+QxxPlotOpts.targPc);
lowerLim(1:length(plotPosMaxRng)) = targVal*(1-QxxPlotOpts.targPc);
if isfield(QxxPlotOpts,'axisType')
    if strcmp(QxxPlotOpts.axisType,'Specified')
        axisLims = [0,max(plotPosMaxRng),min(QxxPlotOpts.QxxRange)*10^-9,max(QxxPlotOpts.QxxRange)*10^-9];
    elseif strcmp(QxxPlotOpts.axisType,'StartZero')
        axisLims = [0,max(plotPosMaxRng),0,targVal*(1+QxxPlotOpts.rangePc)];    
    else
        axisLims = [0,max(plotPosMaxRng),targVal*(1-QxxPlotOpts.rangePc),targVal*(1+QxxPlotOpts.rangePc)];
    end
else
    axisLims = [0,max(plotPosMaxRng),targVal*(1-QxxPlotOpts.rangePc),targVal*(1+QxxPlotOpts.rangePc)];
end

% Setup the legend strings
for s = 1:numSpecs
    legendStrs{s} = ['S',num2str(s)];
end
legendStrs{numSpecs+1} = 'QS Ref';
legendStrs{numSpecs+2} = ['QS Ref \pm',num2str(QxxPlotOpts.targPc*100),'%'];

%--------------------------------------------------------------------------
% Create the figure and size it
hf = figure;
set(hf,'Position', [plotProps.locX,plotProps.locY,...
        2*plotProps.sizePerFigX,1.6*plotProps.sizePerFigY])
set(hf,'PaperPositionMode','auto')
labelStrs.x = 'X position (mm)';

%--------------------------------------------------------------------------
hold on
% PLOT - indentification over the length
for s = 1:numSpecs
    plot(plotPos{s}(posRange{s}),plotIdentStiff{s}(posRange{s}),...
        LS{s},'linewidth',plotProps.lw,'markersize',plotProps.ms)
end

% PLOT QS ref and bounds
plot(plotPosMaxRng,target,'--k','linewidth',plotProps.lw,'markersize',plotProps.ms)
plot(plotPosMaxRng,upperLim,'-k',plotPosMaxRng,lowerLim,'-k','linewidth',plotProps.lw,'markersize',plotProps.ms)

% PLOT - identified average
if QxxPlotOpts.plotAverage
for s = 1:numSpecs
    plot(plotPosMaxRng,plotIdentAvg{s},avgLS{s},'linewidth',plotProps.lw,'markersize',plotProps.ms)
end
end

%--------------------------------------------------------------------------
xlabel(labelStrs.x,'fontsize',plotProps.hfs,'fontname',plotProps.ft)
ylabel(labelStrs.y,'fontsize',plotProps.hfs,'fontname',plotProps.ft)
set(gca, 'fontsize', plotProps.fs,'fontname',plotProps.ft)
set(gca,'XMinorTick','on','YMinorTick','on')
legend(legendStrs,'location','eastoutside')
axis(axisLims)
box on
grid on
hold off

end

