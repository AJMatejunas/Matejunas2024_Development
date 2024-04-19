function hf = func_plotQxxVsLength(QxxPlotOpts,pos,material,identStiffness)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 10/3/2017
%
% Plots the identified Qxx vs length against the target value value giving
% +/-pc bounds on the target value 
% Method=1: linearly fit Stress,x(avg) vs Strain,x(avg) giving E (small
% poisson effect assumption)
% Method=2: linearly fit Stress,x(avg) vs Strain,x+nu*Strain,y(avg) giving
% Qxx, accounts for the poisson effect

% Get some plot parameters for nice formatting
plotProps = func_initPlotPropsStruct('presentation');

% Select the target value based on the method chosen
if QxxPlotOpts.method == 2
    labelStrs.y = 'Stiffness Qxx (GPa)';
    targVal = material.Qxx*10^-9;
    identAvg(1:length(pos.x)) = identStiffness.QxxAvgOverL;
    identStiff = identStiffness.QxxVsL;
else
    labelStrs.y = 'Elastic Modulus E_{x} (GPa)';
    targVal = material.Exx*10^-9;
    identAvg(1:length(pos.x)) = identStiffness.ExxAvgOverL;
    identStiff = identStiffness.ExxVsL;
end

% Setup the plotting vectors for the target and the error bounds
target(1:length(pos.x)) = targVal;
upperLim(1:length(pos.x)) = targVal*(1+QxxPlotOpts.targPc);
lowerLim(1:length(pos.x)) = targVal*(1-QxxPlotOpts.targPc);
if isfield(QxxPlotOpts,'yAxisStartZero')
    if QxxPlotOpts.yAxisStartZero
        axisLims = [0,max(pos.x)*10^3,0,targVal*(1+QxxPlotOpts.rangePc)];
    else
        axisLims = [0,max(pos.x)*10^3,targVal*(1-QxxPlotOpts.rangePc),targVal*(1+QxxPlotOpts.rangePc)];
    end
else
    axisLims = [0,max(pos.x)*10^3,targVal*(1-QxxPlotOpts.rangePc),targVal*(1+QxxPlotOpts.rangePc)];
end

% Setup the legend strings
if QxxPlotOpts.plotAverage
    legendStrs{1} = 'Identified';
    legendStrs{2} = 'Identified Avg';
    legendStrs{3} = 'QS Ref';
    legendStrs{4} = ['QS Ref \pm',num2str(QxxPlotOpts.targPc*100),'%'];
else
    legendStrs{1} = 'Identified';
    legendStrs{2} = 'QS Ref';
    legendStrs{3} = ['QS Ref \pm',num2str(QxxPlotOpts.targPc*100),'%'];
end

hf = figure;
set(hf,'Position', [plotProps.locX,plotProps.locY,...
        2*plotProps.sizePerFigX,1.6*plotProps.sizePerFigY])
set(hf,'PaperPositionMode','auto')
labelStrs.x = 'Specimen Slice X (mm)';
hold on
plot(pos.x*10^3,identStiff*10^-9,'-xb','linewidth',plotProps.lw,'markersize',plotProps.ms)
if QxxPlotOpts.plotAverage
    plot(pos.x*10^3,identAvg*10^-9,'-.b','linewidth',plotProps.lw,'markersize',plotProps.ms)
end
plot(pos.x*10^3,target,'--k','linewidth',plotProps.lw,'markersize',plotProps.ms)
plot(pos.x*10^3,upperLim,'-k',pos.x*10^3,lowerLim,'-k','linewidth',plotProps.lw,'markersize',plotProps.ms)
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

