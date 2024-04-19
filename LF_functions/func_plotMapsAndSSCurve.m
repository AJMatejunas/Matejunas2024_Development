function func_plotMapsAndSSCurve(imageSeqSavePath,plotParams,labelStrs,pos,time,...
   ssCurveLoc,ssCurveMethod,stress,strain,accel,strainRate,plotTrendLine,identStiffness)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 29/3/2017
%
% Plots 2D frames of data as an image sequence along with a stress strain
% curve at a specified location on the specimen

    
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
    %//////////////////////////////////////////////////////////////////////


    if nargin < 12
        plotTrendLine = false;
        identStiffness = 0;
    end
    
    % Create a struct of formatting properties for the figure
    plotProps = func_initPlotPropsStruct(plotParams.formatType);

    % Specify the strain variable to use
    if ssCurveMethod == 2
        strainVar = strain.xnyAvg;
        strainLabel = '$\epsilon_{xx}+\nu\epsilon_{yy}$ [$mm.m^{-1}$]';
    else
        strainVar = strain.xAvg;
        strainLabel = '$\epsilon_{xx}$ [$mm.m^{-1}$]';
    end
    
    % Specify ranges for cropping pixels on image maps
    if  plotParams.cutEdgePx
        plotRangeX = (plotParams.cutPxX+1):size(strain.x,2)-plotParams.cutPxX;
        plotRangeY = (plotParams.cutPxY+1):size(strain.x,1)-plotParams.cutPxY;
    else
        plotRangeX = 1:size(strain.x,2);
        plotRangeY = 1:size(strain.x,1);
    end
    plotXPos = pos.x(plotParams.cutPxX+1:end-plotParams.cutPxX)*10^3;
    plotYPos = pos.y(plotParams.cutPxY+1:end-plotParams.cutPxY)*10^3;
    
    % Create variables for plotting the trendline on the SS curve
    if plotTrendLine
        trendLineStrain = linspace(min(strainVar(ssCurveLoc,plotParams.tRange)),...
            max(strainVar(ssCurveLoc,plotParams.tRange)),10);
        trendLineStress = identStiffness.linearFitCoeffs{ssCurveLoc}(1).*trendLineStrain +...
            identStiffness.linearFitCoeffs{ssCurveLoc}(2);
    end
    
    % Setup vars for location lines on the image
    locLineX(1:size(strain.x,1)) = pos.x(ssCurveLoc)*10^3;
    locLineY = pos.y(plotRangeY)*10^3;
    
    % Calculate the colorbar range for each variable to be plotted
    if isfield(plotParams,'cAxisType')
        plotParams.cRange{1} = func_calcColourBarRange(plotParams.cAxisType,accel.x(plotRangeY,plotRangeX,:));
        plotParams.cRange{2} = func_calcColourBarRange(plotParams.cAxisType,strain.x(plotRangeY,plotRangeX,:));
        plotParams.cRange{3} = func_calcColourBarRange(plotParams.cAxisType,strainRate.x(plotRangeY,plotRangeX,:));
    end
    
    % Create and size the figure
    hf = func_createFigure(plotProps,plotParams);
    
    % Loop over each frame and plot it
    for f = plotParams.tRange
        %------------------------------------------------------------------
        % SUBPLOT 1: Plot the stress strain curve 
        subplot(plotParams.Rows,plotParams.Cols,1)
        hold on
        plot(strainVar(ssCurveLoc,plotParams.tRange)*10^3,stress.xAvg(ssCurveLoc,plotParams.tRange)*10^-6,...
                '-xb','linewidth',plotProps.lw*1.5,'markersize',plotProps.ms*1.5)
        if plotTrendLine
            plot(trendLineStrain*10^3,trendLineStress*10^-6,...
                '-k','linewidth',plotProps.lw*1.5)    
        end
        hp = plot(strainVar(ssCurveLoc,f)*10^3,stress.xAvg(ssCurveLoc,f)*10^-6,...
                    'or','markersize',plotProps.ms*1.5);
        set(hp, 'MarkerFaceColor', get(hp, 'Color'));
        hold off
        title(['x = ',num2str(pos.x(ssCurveLoc)*10^3),'mm'],'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
        xlabel(strainLabel,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
        ylabel('$\sigma_{xx}$ [$MPa$]','fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
        xlim([round(1.25*min(strainVar(ssCurveLoc,plotParams.tRange)),2,'significant'),...
            round(1.25*max(strainVar(ssCurveLoc,plotParams.tRange)),2,'significant')]*10^3);
        ylim([round(1.1*min(stress.xAvg(ssCurveLoc,plotParams.tRange)),2,'significant'),...
            round(1.1*max(stress.xAvg(ssCurveLoc,plotParams.tRange)),2,'significant')]*10^-6);
        set(gca, 'fontsize', plotProps.fs,'fontname',plotProps.ft)
        set(gca,'XMinorTick','on','YMinorTick','on')
        box on
        grid on
        
        %--------------------------------------------------------------------------
        % SUBPLOT 2: Acceleration
        subplot(plotParams.Rows,plotParams.Cols,2)
        imagesc(plotXPos,plotYPos,accel.x(plotRangeY,plotRangeX,f))
        hold on
        plot(locLineX,locLineY,'--k','linewidth',plotProps.lw)
        hold off
        xlabel('X [mm]','fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
        ylabel('Y [mm]','fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
        titleStr = ['$a_{x}$ [$m.s^{-2}$], ',sprintf('t = %.2f',time.vec(f)*10^6),'$\mu s$'];
        title(titleStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
        set(gca,'fontsize',plotProps.fs,'fontname',plotProps.ft,'linewidth',plotProps.lw)
        colorbar
        colormap(jet)
        caxis(plotParams.cRange{1})
        axis image
        
        %--------------------------------------------------------------------------
        % SUBPLOT 3: Strain
        subplot(plotParams.Rows,plotParams.Cols,3)
        imagesc(plotXPos,plotYPos,strain.x(plotRangeY,plotRangeX,f)*1e03)
        hold on
        plot(locLineX,locLineY,'--k','linewidth',plotProps.lw)
        hold off
        xlabel('X [mm]','fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
        ylabel('Y [mm]','fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
        titleStr = ['$\epsilon_{xx}$ [$mm.m^{-1}$], ',sprintf('t = %.2f',time.vec(f)*10^6),'$\mu s$'];
        title(titleStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
        set(gca,'fontsize',plotProps.fs,'fontname',plotProps.ft,'linewidth',plotProps.lw)
        colorbar
        colormap(jet)
        caxis(plotParams.cRange{2}*10^3)       
        axis image
        
        %--------------------------------------------------------------------------
        % SUBPLOT4: Strain Rate
        subplot(plotParams.Rows,plotParams.Cols,4)
        imagesc(plotXPos,plotYPos,strainRate.x(plotRangeY,plotRangeX,f))
        hold on
        plot(locLineX,locLineY,'--k','linewidth',plotProps.lw)
        hold off
        xlabel('X [mm]','fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
        ylabel('Y [mm]','fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
        titleStr = ['$d\epsilon_{xx}/dt$ [$s^{-1}$], ',sprintf('t = %.2f',time.vec(f)*10^6),'$\mu s$'];
        title(titleStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
        set(gca,'fontsize',plotProps.fs,'fontname',plotProps.ft,'linewidth',plotProps.lw)
        colorbar
        caxis(plotParams.cRange{3})
        axis image
        
        % Save this frame to file and clear the figure
        print(hf,[imageSeqSavePath,'\SSMapVideo_Frame',num2str(f)],plotProps.format,'-r0')
        clf(hf)
    end
end
