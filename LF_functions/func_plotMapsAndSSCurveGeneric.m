function func_plotMapsAndSSCurveGeneric(savePath,plotParams,labelStrs...
    ,pos,time,ssCurve,mapVars,locLine)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 29/3/2017
% Date Updated: 30/1/2019
%
% Plots 2D frames of data as an image sequence along with a stress strain
% curve at a specified location on the specimen

    if nargin < 8
        locLine.specify = false;
    end

    % Create a struct of formatting properties for the figure
    plotProps = func_initPlotPropsStruct(plotParams.formatType);
    
    % Specify the range of data points to plot on the image  
    if  plotParams.cutEdgePx
        rangeX = (plotParams.cutPxX+1):length(pos.x)-plotParams.cutPxX;
        rangeY = (plotParams.cutPxY+1):length(pos.y)-plotParams.cutPxY;
    else
        rangeX = 1:length(pos.x);
        rangeY = 1:length(pos.y);
    end
    plotXPos = pos.x(rangeX)*10^3;
    plotYPos = pos.y(rangeY)*10^3;

    % Calculate the colorbar range for each variable to be plotted
    if isfield(plotParams,'cAxisType')
        if ~strcmp(plotParams.cAxisType,'Specified')                
            for p = 2:length(mapVars) 
                plotParams.cRange{p} = func_calcColourBarRange(plotParams.cAxisType,...
                    mapVars{p}(rangeY,rangeX,:));
            end
        end
    end
  
    % Create variables for plotting the trendline on the SS curve
    if ssCurve.plotTrendLine
        trendLineStrain = linspace(min(ssCurve.strain(plotParams.tRange)),...
            max(ssCurve.strain(plotParams.tRange)),10);
        trendLineStress = ssCurve.linFitCoeffs(1).*trendLineStrain +...
            ssCurve.linFitCoeffs(2);
    end
    
    % Setup vars for location lines on the image
    if ~locLine.specify
        locLine.x(1:length(rangeY)) = pos.x(ssCurve.locX)*10^3;
        locLine.y = pos.y(rangeY)*10^3;
    end
    
    % Calculate the colorbar range for each variable to be plotted
    if isfield(plotParams,'cAxisType')
        if ~strcmp(plotParams.cAxisType,'Specified')                
            for p = 2:length(mapVars)
                plotParams.cRange{p} = func_calcColourBarRange(plotParams.cAxisType,...
                    mapVars{p}(rangeY,rangeX,:));
            end
        end
    end
    
    % Create and size the figure
    hf = func_createFigure(plotProps,plotParams);
    
    % Loop over each frame and plot it
    for ff = plotParams.tRange
        %------------------------------------------------------------------
        % SUBPLOT 1: Plot the stress strain curve 
        plotNum = 1;
        subplot(plotParams.Rows,plotParams.Cols,plotNum)
        hold on
        plot(ssCurve.strain(plotParams.tRange)*10^3,ssCurve.stress(plotParams.tRange)*10^-6,...
                '-+b','linewidth',plotProps.lw,'markersize',plotProps.ms)
        if ssCurve.plotTrendLine
            plot(trendLineStrain*10^3,trendLineStress*10^-6,...
                '-k','linewidth',plotProps.lw)    
        end
        hp = plot(ssCurve.strain(ff)*10^3,ssCurve.stress(ff)*10^-6,...
                    'or','markersize',plotProps.ms);
        set(hp, 'MarkerFaceColor', get(hp, 'Color'));
        hold off
        title(labelStrs.t{plotNum},'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
        xlabel(labelStrs.x{plotNum},'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
        ylabel(labelStrs.y{plotNum},'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
        xlim([round(1.25*min(ssCurve.strain(plotParams.tRange)),2,'significant'),...
            round(1.25*max(ssCurve.strain(plotParams.tRange)),2,'significant')]*10^3);
        ylim([round(1.1*min(ssCurve.stress(plotParams.tRange)),2,'significant'),...
            round(1.1*max(ssCurve.stress(plotParams.tRange)),2,'significant')]*10^-6);
        set(gca, 'fontsize', plotProps.fs,'fontname',plotProps.ft)
        set(gca,'XMinorTick','on','YMinorTick','on')
        box on
        grid on
        
        for mm = 2:length(mapVars) 
            %--------------------------------------------------------------------------
            % SUBPLOT X: Map Variables
            plotNum = mm;
            subplot(plotParams.Rows,plotParams.Cols,plotNum)
            imagesc(plotXPos,plotYPos,mapVars{plotNum}(rangeY,rangeX,ff))
            hold on
            plot(locLine.x,locLine.y,'--k','linewidth',plotProps.lw)
            hold off
            titleStr = [labelStrs.t{plotNum},', ',sprintf('t = %.2f',time.vec(ff)*10^6),'$\mu s$'];
            title(titleStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
            xlabel(labelStrs.x{plotNum},'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
            ylabel(labelStrs.y{plotNum},'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
            set(gca,'fontsize',plotProps.fs,'fontname',plotProps.ft,'linewidth',plotProps.lw)
            colorbar
            colormap(jet)
            caxis(plotParams.cRange{plotNum})
            axis image
        end
        
        % Save this frame to file
        saveFile = [savePath,'\Frame_',num2str(ff)];
        print(hf,saveFile,plotProps.imageSeqFormat,plotProps.imageSeqSaveRes);
        if plotParams.saveImageMatFig
            saveas(hf,saveFile,'fig');
        end
        if plotParams.saveImageVecFig
            export_fig(saveFile,'-eps','-pdf','-c[Inf,Inf,Inf,Inf]',hf)
        end
        clf(hf);
    end
end
