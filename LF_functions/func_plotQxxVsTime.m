function hf = func_plotQxxVsTime(QxxPlotOpts,time,material,identStiffness)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 29/3/2017
%
% Plots the identified Qxx vs time against the target value value giving
% +/-pc bounds on the target value 

% Get some plot parameters for nice formatting
plotProps = func_initPlotPropsStruct('presentation');

% Select the target value based on the method chosen
if QxxPlotOpts.method == 2
    labelStrs.y = 'Stiffness, Qxx, (GPa)';
    targVal = material.Qxx*10^-9;
    plotStiffness = identStiffness.QxxVsT*10^-9;
    identAvg(1:length(time.vec)) = identStiffness.QxxAvgOverT;
else
    labelStrs.y = 'Elastic Modulus, E_{x}, (GPa)';
    targVal = material.Exx*10^-9;
    plotStiffness = identStiffness.ExxVsT*10^-9;
    identAvg(1:length(time.vec)) = identStiffness.ExxAvgOverT;
end

% Setup the plotting vectors for the target and the error bounds
target(1:length(time.vec)) = targVal;
upperLim(1:length(time.vec)) = targVal*(1+QxxPlotOpts.targPc);
lowerLim(1:length(time.vec)) = targVal*(1-QxxPlotOpts.targPc);
axisLims = [0,time.vec(end)*10^6,targVal*(1-QxxPlotOpts.rangePc),targVal*(1+QxxPlotOpts.rangePc)];

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
labelStrs.x = 'Time, (\mus)';
hold on
plot(time.vec*10^6, plotStiffness,'-xb','linewidth',plotProps.lw,'markersize',plotProps.ms)
if QxxPlotOpts.plotAverage
    plot(time.vec*10^6,identAvg*10^-9,'-.b','linewidth',plotProps.lw,'markersize',plotProps.ms)
end
plot(time.vec*10^6,target,'--k','linewidth',plotProps.lw,'markersize',plotProps.ms)
plot(time.vec*10^6,upperLim,'-k',time.vec*10^6,lowerLim,'-k','linewidth',plotProps.lw,'markersize',plotProps.ms)
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

