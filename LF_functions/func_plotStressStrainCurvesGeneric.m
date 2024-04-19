function hf = func_plotStressStrainCurvesGeneric(labelStrs,plotParams,pos,...
    stressAvg,strainAvg,identStiffSG)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 9/2/2018
% Date Edited: 23/1/2019
%
% Plots stress strain curves at different locations on the specimen, this
% is a generic plotting function, it is up to the user to ensure the input
% stress and strain averages are consistent.

if nargin < 6
    identStiffSG = nan;
end

% Get the plot properties structure for  consistent formatting
plotProps = func_initPlotPropsStruct(plotParams.formatType);

% Setup the figure
hf = func_createFigure(plotProps,plotParams);

for xx = 1:length(plotParams.locXInd)
    currLocX = plotParams.locXInd(xx);
    
    subplot(plotParams.Rows,plotParams.Cols,xx)
    hold on
    if plotParams.plotOnlyFittedRegion && isstruct(identStiffSG)
       plot(strainAvg(currLocX,identStiffSG.QxxFitFrameRange{currLocX})*10^3,...
             stressAvg(currLocX,identStiffSG.QxxFitFrameRange{currLocX})*10^-6,...
                            '-+b','linewidth',plotProps.lw,'markersize',plotProps.ms) 
    else
        plot(strainAvg(plotParams.locXInd(xx),plotParams.tRange)*10^3,...
             stressAvg(plotParams.locXInd(xx),plotParams.tRange)*10^-6,...
                            '-+b','linewidth',plotProps.lw,'markersize',plotProps.ms)
    end
    
    if plotParams.plotFittedCurve && isstruct(identStiffSG)
        xMin = min(strainAvg(currLocX,identStiffSG.QxxFitFrameRange{currLocX}));
        xMax = max(strainAvg(currLocX,identStiffSG.QxxFitFrameRange{currLocX}));
        xVec = linspace(xMin,xMax,10);
        yVec = identStiffSG.QxxLinFitCoeffs{currLocX}(1).*xVec + identStiffSG.QxxLinFitCoeffs{currLocX}(2);
        plot(xVec*10^3,yVec*10^-6,'-k','linewidth',plotProps.lw,'markersize',plotProps.ms)
    end
    
    
    if isfield(labelStrs,'t')
        title([labelStrs.t,', x = ',sprintf('%0.2f',pos.x(plotParams.locXInd(xx))*10^3),'mm'],'Interpreter',plotProps.interpText)
    else
        title(['x = ',sprintf('%0.2f',pos.x(plotParams.locXInd(xx))*10^3),'mm'],'Interpreter',plotProps.interpText)
    end
    xlabel(labelStrs.x,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText)
    ylabel(labelStrs.y,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText)
    set(gca, 'fontsize', plotProps.fs,'fontname',plotProps.ft)
    set(gca,'XMinorTick','on','YMinorTick','on')
    box on
    grid on
    hold off 
end

end

