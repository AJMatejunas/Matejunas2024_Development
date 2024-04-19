function func_plotFullFieldImageSeq(savePath,plotParams,labelStrs,pos,time,...
    plotVars)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 13/3/2017
% Date Edited: 25/4/2019
%
% Plots 2D frames of kinematic data as an image sequence
    
    % Create a struct of formatting properties for the figure
    plotProps = func_initPlotPropsStruct(plotParams.formatType);
    
    % Specify the range of data points to plot on the image  
    if  plotParams.cutEdgePx
        rangeX = (plotParams.cutPxX+1):size(plotVars{1},2)-plotParams.cutPxX;
        rangeY = (plotParams.cutPxY+1):size(plotVars{1},1)-plotParams.cutPxY;
    else
        rangeX = 1:size(plotVars{1},2);
        rangeY = 1:size(plotVars{1},1);
    end
    
    % Check that the number of subplots and number of vars in the plot vars 
    % cell is consistent
    numPlots = length(plotVars);
    checkNumPlots = plotParams.Rows*plotParams.Cols;
    if numPlots ~= checkNumPlots
        disp('WARNING: number of plot vars is not equal to specified number of subplots.')
    end
    
    % Calculate the colorbar range for each variable to be plotted
    if isfield(plotParams,'cAxisType')
        if ~strcmp(plotParams.cAxisType,'Specified')                
            for p = 1:numPlots 
                plotParams.cRange{p} = func_calcColourBarRange(plotParams.cAxisType,...
                    plotVars{p}(rangeY,rangeX,:));
            end
        end
    end
        
    % Create and size the figure
    hf = func_createFigure(plotProps,plotParams);
    
    % Loop over each frame and plot it
    for f = 1:time.numFrames
        for p = 1:numPlots
            hold on
            % Create the subplot
            subplot(plotParams.Rows,plotParams.Cols,p)
            % Plot the image variable frame by fram
            imagesc(pos.x(rangeX)*10^3,pos.y(rangeY)*10^3,plotVars{p}(rangeY,rangeX,f))
            if isfield(plotParams,'titleFrameNum')
                if plotParams.titleFrameNum
                    titleStr = {labelStrs.t{p},[sprintf('Time: %.1f',time.vec(f)*10^6),'$\mu s$, ',sprintf('Frame: %i',f)]};
                else
                    titleStr =[labelStrs.t{p},sprintf(', Time: %.1f',time.vec(f)*10^6),'$\mu s$'];
                end
            else
                titleStr = [labelStrs.t{p},sprintf(', Time: %.1f',time.vec(f)*10^6),'$\mu s$'];
            end
            title(titleStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText) 
            set(gca, 'fontsize', plotProps.fs,'fontname',plotProps.ft)
            axis image
            % Create a colour bar and set the colour map
            if isfield(plotParams,'cAxisLoc')
                colorbar(plotParams.cAxisLoc)
            else
                colorbar
            end
            colormap(jet)
            % Determine the colour bar range 
            if isfield(plotParams,'cAxisType')
                if ~strcmp(plotParams.cAxisType,'Auto')
                    if isfield(plotParams,'cRange')
                        caxis(plotParams.cRange{p})
                    end
                end
            end
            % Remove the image axes if required
            if ~plotParams.imageAxis
                axis off
            else
                xlabel(labelStrs.x,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
                ylabel(labelStrs.y,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
            end
            hold off  
        end
        
        % Save this frame to file
        saveFile = [savePath,'\Frame_',num2str(f)];
        func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)
        clf(hf);
    end

end

