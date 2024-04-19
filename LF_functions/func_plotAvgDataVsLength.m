function hf = func_plotAvgDataVsLength(labelStrs,time,plotParams,posVars,plotVars)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 9/3/2017
% Plots a width averaged full-field variable over the specimen length at a 
% particular time (frame)  
    plotProps = func_initPlotPropsStruct();
    
    % Check if there are legend labels, if not turn off the legend
    legendFlag = false;
    if isfield(labelStrs,'legStrs')  
        legendFlag = true;
    end
    
    % Check that the number of variables to plot is the same in all
    % relevant vars
    if length(posVars) ~= length(plotVars)
        error('The position variable and plot variable must be the same length')
    end
    
    hf = figure;
    set(hf,'Position',[plotProps.locX,plotProps.locY,...
        plotParams.Cols*plotProps.sizePerFigX,plotParams.Rows*plotProps.sizePerFigY])
    set(hf,'PaperPositionMode','auto')
    for i = 1:length(plotParams.frames)
        subplot(plotParams.Rows,plotParams.Cols,i)
        hold on
        
        for j = 1:length(plotVars)
            plot(posVars{j}.x,plotVars{j}(:,plotParams.frames(i)),plotParams.markStrs{j},...
                'linewidth',plotProps.lw,'markersize',plotProps.ms)
        end
   
        title({['Frame = ',num2str(plotParams.frames(i)),', '],...
            ['Time t = ',num2str((plotParams.frames(i)-1)*time.step*10^6),'\mus']})
        xlabel(labelStrs.x,'fontsize',plotProps.hfs,'fontname',plotProps.ft)
        ylabel(labelStrs.y,'fontsize',plotProps.hfs,'fontname',plotProps.ft)
        if legendFlag
            legend(labelStrs.legStrs)
        end
        set(gca, 'fontsize', plotProps.fs,'fontname',plotProps.ft)
        set(gca,'XMinorTick','on','YMinorTick','on')
        box on
        grid on
        hold off 
    end

end

