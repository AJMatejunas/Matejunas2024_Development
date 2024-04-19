function func_plotImageSeqAvgDataVsLength(plotParams,labelStrs,savePath,...
    pos,time,plotVars)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 9/3/2017
%
% Plots and image sequence of average data over the width

    % Create a plot props structure for nice formatting
    plotProps = func_initPlotPropsStruct();
    
    % Check that the number of subplots and number of vars in the plot vars 
    % cell is consistent
    numPlots = length(plotVars);
    checkNumPlots = plotParams.Rows*plotParams.Cols;
    if numPlots ~= checkNumPlots
        disp('WARNING: number of plot vars is not equal to specified number of subplots.')
    end
    
    % Create and size the figure 
    hf = figure('Visible','On');
    set(hf,'Position', [plotProps.locX,plotProps.locY,...
        1.5*plotParams.Cols*plotProps.sizePerFigX,1.5*plotParams.Rows*plotProps.sizePerFigY])
    set(hf,'PaperPositionMode','auto')
    % Loop over the frames to create the video sequence
    for f = 1:time.numFrames
        for p = 1:length(plotVars)
            subplot(plotParams.Rows,plotParams.Cols,p)
            hold on
            plot(pos.x,plotVars{p}(:,f),'-xb','linewidth',plotProps.lw,'markersize',plotProps.ms)
            title([sprintf('Time t = %.1f',time.vec(f)*10^6),'\mus'])
            xlabel(labelStrs.x{1},'fontsize',plotProps.hfs,'fontname',plotProps.ft)
            ylabel(labelStrs.y{p},'fontsize',plotProps.hfs,'fontname',plotProps.ft)
            axis(plotParams.axisRange{p})
            set(gca, 'fontsize', plotProps.fs,'fontname',plotProps.ft)
            set(gca,'XMinorTick','on','YMinorTick','on')
            box on
            grid on
            hold off  
        end
        % Save this frame of the video to file
        print(hf,[savePath,'\Frame',num2str(f),'_AvgData',],plotProps.format,'-r0')
        clf(hf)
    end
end

