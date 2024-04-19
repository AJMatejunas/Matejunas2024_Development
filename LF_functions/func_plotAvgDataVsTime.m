function hf = func_plotAvgDataVsTime(labelStrs,time,plotParams,posVars,plotVars)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 9/3/2017
% Plots a width averaged full-field variable over time at a particular
% specimen location
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
    
    % Find positions in the position vecs and check they are the same
    varLocX = zeros(length(plotParams.locXmm),length(posVars));
    varLocXmm = zeros(length(plotParams.locXmm),length(posVars));
    for n = 1:length(plotParams.locXmm)
        if length(posVars) > 1
            for m = 1:length(posVars)
                [c,index] = min(abs(posVars{m}.x-plotParams.locXmm(n)));
                varLocX(n,m) = index;
                varLocXmm(n,m) = c;     
            end
            if range(varLocXmm(n,:)) ~= 0
                fprintf('WARNING: position is not the same in given variables')
            end
        else
            [c,index] = min(abs(posVars{1}.x-plotParams.locXmm(n)));
            varLocX(n) = index;
            varLocXmm(n) = c;     
        end
    end
    
    if length(posVars) > 1
        varLocX = varLocX(:,1);
        % Check that the time and plotVar are the same length
        if size(plotVars{2},2) ~= length(time.vec)
            cutFrames = time.cutFrames;
        else
            cutFrames = 0;
        end
    else
        cutFrames = 3;
    end
    
    % Setup and plot the figure
    hf = figure;
    set(hf,'Position',[plotProps.locX,plotProps.locY,...
        plotParams.Cols*plotProps.sizePerFigX,plotParams.Rows*plotProps.sizePerFigY])
    set(hf,'PaperPositionMode','auto')
    for i = 1:length(varLocX)
        subplot(plotParams.Rows,plotParams.Cols,i)
        hold on

        for j = 1:length(plotVars)
            plot(time.vec(1:end-cutFrames)*10^6,plotVars{j}(varLocX(i),1:end-cutFrames),plotParams.markStrs{j},...
                'linewidth',plotProps.lw,'markersize',plotProps.ms)
        end
        
        title(['x = ',num2str(posVars{1}.x(varLocX(i))*10^3),'mm'])
        xlabel(labelStrs.x,'fontsize',plotProps.hfs,'fontname',plotProps.ft)
        ylabel(labelStrs.y,'fontsize',plotProps.hfs,'fontname',plotProps.ft)
        if legendFlag
            legend(labelStrs.legStrs)
        end
        set(gca, 'fontsize', plotProps.fs,'fontname',plotProps.ft)
        set(gca,'XMinorTick','on','YMinorTick','on')
        box on
        hold off
    end

end

