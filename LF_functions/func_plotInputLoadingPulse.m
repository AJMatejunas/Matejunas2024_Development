function hf = func_plotInputLoadingPulse(plotParams,labelStrs,time,specimen,stress,ID_specimen,ID_stress)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 10/5/2018

% Check if the image deformation data needs to be plotted.
plotID = true;
if nargin < 7
    plotID = false;
end

% Create a struct of formatting properties for the figure
plotProps = func_initPlotPropsStruct('article_v2');
plotProps.ms = plotProps.ms/2;
%lColour = linspecer(2,'qualitative');
lColour = 'b';

% Calculate the pulses
if plotID
    IDStressPulse = ID_stress.xAvg(end,:)*10^-6;
    IDForcePulse = ID_stress.xAvg(end,:)*ID_specimen.height*ID_specimen.thickness*10^-3;
end


% Get the stress pulse and calculate the force given the specimen geometry
specStressPulse = stress.xAvg(end,:)*10^-6;
specForcePulse = stress.xAvg(end,:)*specimen.height*specimen.thickness*10^-3;

% If the axis limits have not been set then set them based on the max force
if ~isfield(plotParams,'forceLims')
    plotParams.forceLims = 1.1*[min(specForcePulse),max(specForcePulse)];
end

plotParams.stressLims = (plotParams.forceLims*10^3/(specimen.height*specimen.thickness))*10^-6;

%//////////////////////////////////////////////////////////////////////////
% Create the figure
hf = func_createFigure(plotProps,plotParams);
left_color = [0 0 0];
right_color = [0 0 0];
set(hf,'defaultAxesColorOrder',[left_color; right_color])
hold on 

yyaxis left
plot(time.vec*10^6,specForcePulse,plotProps.singleLS,'color',lColour(1,:)...
    ,'linewidth',plotProps.lw,'markersize',plotProps.ms)
if plotID 
    plot(time.vec*10^6,IDForcePulse,plotProps.singleLS,'color',lColour(2,:)...
        ,'linewidth',plotProps.lw,'markersize',plotProps.ms)
end
ylabel(labelStrs.y1Str,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'interpreter','LaTeX')
ylim([plotParams.forceLims])

yyaxis right
plot(time.vec*10^6,specStressPulse,plotProps.singleLS,'color',lColour(1,:)...
    ,'linewidth',plotProps.lw,'markersize',plotProps.ms)
if plotID
    plot(time.vec*10^6,IDStressPulse,plotProps.singleLS,'color',lColour(2,:)...
        ,'linewidth',plotProps.lw,'markersize',plotProps.ms)
end
ylabel(labelStrs.y2Str,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'interpreter','LaTeX')
ylim([plotParams.stressLims])

hold off
%title(tStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'interpreter','LaTeX')
xlabel(labelStrs.xStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'interpreter','LaTeX')
if plotID
    legend(labelStrs.legStrs,'location','best','interpreter','LaTeX')
end
set(gca,'fontsize', plotProps.fs,'fontname',plotProps.ft)
set(gca,'XMinorTick','on','YMinorTick','on')
box on
grid on

end

