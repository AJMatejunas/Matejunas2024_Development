function func_plotMapAndSSCurve(imageSeqSavePath,plotParams,labelStrs,pos,time,...
    mapVar,ssCurveInd,ssCurveMethod,stress,strain)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 13/3/2017
%
% Plots 2D frames of data as an image sequence along with a stress strain
% curve at a specified location on the specimen
    
    % Create a struct of formatting properties for the figure
    plotProps = func_initPlotPropsStruct(plotParams.formatType);

    % Specify the strain variable to use
    if ssCurveMethod == 2
        strainVar = strain.xnyAvg;
    else
        strainVar = strain.xAvg;
    end
    
    % Specify ranges for cropping pixels on image maps
    if  plotParams.cutEdgePx
        rangeX = (plotParams.cutPxX+1):size(mapVar,2)-plotParams.cutPxX;
        rangeY = (plotParams.cutPxY+1):size(mapVar,1)-plotParams.cutPxY;
    else
        rangeX = 1:size(mapVar,2);
        rangeY = 1:size(mapVar,1);
    end
    
    % Setup vars for location lines on the image
    locLineX(1:size(mapVar,1)) = pos.x(ssCurveInd)*10^3;
    locLineY = pos.y(rangeY)*10^3;
    
    % Calculate the colorbar range for each variable to be plotted
    if isfield(plotParams,'cAxisType')
        plotParams.cRange{1} = func_calcColourBarRange(plotParams.cAxisType,...
            mapVar(rangeY,rangeX,:));
    end
    
    % Create and size the figure
    hf = figure('Visible','On');
    set(hf,'Position', [plotProps.locX,plotProps.locY,...
        1.5*plotParams.Cols*plotProps.sizePerFigX,1.5*plotParams.Rows*plotProps.sizePerFigY])
    set(hf,'PaperPositionMode','auto')
    % Loop over each frame and plot it
    for f = plotParams.tRange
        %------------------------------------------------------------------
        % Plot the heat map variable in the left hand subplot
        subplot(plotParams.Rows,plotParams.Cols,1)
        % Create the heat map
        imagesc(pos.x(rangeX)*10^3,pos.y(rangeY)*10^3,mapVar(rangeY,rangeX,f))
        % Plot the location of the stress strain curve on the heat map
        hold on
        plot(locLineX,locLineY,'--k','linewidth',plotProps.lw)
        hold off
        titleStr = [labelStrs.t{1},sprintf(', Time: %.1f',time.vec(f)*10^6),'\mus'];
        title(titleStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft)
        set(gca, 'fontsize', plotProps.fs,'fontname',plotProps.ft)
        % Create the colour bar and set the colour map
        colorbar
        colormap(jet)
        % Determine the colour bar range 
        if ~strcmp(plotParams.cAxisType,'Auto')
            if isfield(plotParams,'cRange')
                caxis(plotParams.cRange{1})
            end
        end
        % Remove the image axes if required
        if ~plotParams.imageAxis
            axis off
        else
            xlabel(labelStrs.x{1},'fontsize',plotProps.hfs,'fontname',plotProps.ft);
            ylabel(labelStrs.y{1},'fontsize',plotProps.hfs,'fontname',plotProps.ft);
        end

        %------------------------------------------------------------------
        % Plot the stress strain curve in the right hand subplot
        subplot(plotParams.Rows,plotParams.Cols,2)
        hold on
        plot(strainVar(ssCurveInd,plotParams.tRange)*10^3,stress.xAvg(ssCurveInd,plotParams.tRange)*10^-6,...
                '-xb','linewidth',plotProps.lw,'markersize',plotProps.ms)
        hp = plot(strainVar(ssCurveInd,f)*10^3,stress.xAvg(ssCurveInd,f)*10^-6,...
                    'or','markersize',plotProps.ms*1.5);
        set(hp, 'MarkerFaceColor', get(hp, 'Color'));
        hold off
        title(['x = ',num2str(pos.x(ssCurveInd)*10^3),'mm'],'fontsize',plotProps.hfs,'fontname',plotProps.ft)
        xlabel(labelStrs.x{2},'fontsize',plotProps.hfs,'fontname',plotProps.ft)
        ylabel(labelStrs.y{2},'fontsize',plotProps.hfs,'fontname',plotProps.ft)
        xlim([round(1.25*min(strainVar(ssCurveInd,plotParams.tRange)),2,'significant'),...
            round(1.25*max(strainVar(ssCurveInd,plotParams.tRange)),2,'significant')]*10^3);
        ylim([round(1.1*min(stress.xAvg(ssCurveInd,plotParams.tRange)),2,'significant'),...
            round(1.1*max(stress.xAvg(ssCurveInd,plotParams.tRange)),2,'significant')]*10^-6);
        set(gca, 'fontsize', plotProps.fs,'fontname',plotProps.ft)
        set(gca,'XMinorTick','on','YMinorTick','on')
        box on
        grid on 
        
        if plotParams.saveFigs
            print(hf,[imageSeqSavePath,'\SSMapVideo_Frame',num2str(f)],plotProps.format,'-r0')
        end
        % Clear the figure
        clf(hf)
    end
end
