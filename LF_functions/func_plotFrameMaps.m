function fh = func_plotFrameMaps(plotParams,labelStrs,pos,time,plotVar)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 13/3/2017
%
% Plots 2D frames of data 
% Inputs: 
% labelStrs.t - string giving the variable a name
% plotVar - 3D matrix, 3rd dim represents the frame number i.e. time         
% plotParams.frames - frame numbers to plot
    
% Create a struct of formatting properties for the figure
plotProps = func_initPlotPropsStruct();

% Specify the range of data points to plot on the image
if isfield(plotParams,'cutEdgePx')
    if plotParams.cutEdgePx
        rangeX = (plotParams.cutPxX+1):size(plotVar,2)-plotParams.cutPxX;
        rangeY = (plotParams.cutPxY+1):size(plotVar,1)-plotParams.cutPxY;
    else
        rangeX = 1:size(plotVar,2);
        rangeY = 1:size(plotVar,1);
    end
else
    rangeX = 1:size(plotVar,2);
    rangeY = 1:size(plotVar,1);
end

% Check that the specificed frames to plot actually exist
[~,~,totalFrames] = size(plotVar);
numPlots = length(plotParams.frames);   
checkFrames = sum(plotParams.frames>totalFrames);
if checkFrames > 0
    error('Frame specified in plotParams.frames does not exist in plotVar')
    return
end

% Calculate the colorbar range for each variable to be plotted
if isfield(plotParams,'cAxisType')  
    plotParams.cRange{1} = func_calcColourBarRange(plotParams.cAxisType,...
        plotVar(rangeY,rangeX,:));
end

% Create and size the figure
fh = figure;
set(fh,'Position', [plotProps.locX,plotProps.locY,...
    plotParams.Cols*plotProps.sizePerFigX,plotParams.Rows*plotProps.sizePerFigY])
set(fh,'PaperPositionMode','auto')
for i = 1:numPlots
    hold on
    subplot(plotParams.Rows,plotParams.Cols,i)
    imagesc(pos.x(rangeX),pos.y(rangeY),squeeze(plotVar(rangeY,rangeX,plotParams.frames(i))))
    titleStr = {labelStrs.t,['Frame: ',num2str(plotParams.frames(i)),...
        ', Time: ',num2str(time.vec(plotParams.frames(i))*10^6,3)]};
    title(titleStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft)
    set(gca,'fontsize',plotProps.fs,'fontname',plotProps.ft)
    colorbar
    colormap(jet)
    % Determine the colour bar range
    if isfield(plotParams,'cAxisType')
        if ~strcmp(plotParams.cAxisType,'Auto')
            if isfield(plotParams,'cRange')
                caxis(plotParams.cRange{1})
            end
        end
    end
    % Remove the image axes if required
    axis image
    if ~plotParams.imageAxis
        axis off
    else
        xlabel(labelStrs.x,'fontsize',plotProps.hfs,'fontname',plotProps.ft);
        ylabel(labelStrs.y,'fontsize',plotProps.hfs,'fontname',plotProps.ft);
    end
    hold off
end
end

