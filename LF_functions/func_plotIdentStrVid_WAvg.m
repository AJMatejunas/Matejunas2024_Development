function func_plotIdentStrVid_WAvg(plotParams,savePath,pos,time,disp,accel,...
    strain,strainRate,stress,fracture)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 19/2/2018
% Date Edited: 25/4/2019
%
% Plots image sequence of the various full-field maps and the width based
% stress average at the fracture location.
    
    %----------------------------------------------------------------------
    % Initialise Variables for Plotting

    % Create a struct of formatting properties for the figure
    plotProps = func_initPlotPropsStruct(plotParams.formatType);
    plotParams.Cols = 3;
    plotParams.Rows = 2;

    % Create cropped position vectors for plotting
    plotRangeX = plotParams.cutPxX+1:length(pos.x)-plotParams.cutPxX;
    plotRangeY = plotParams.cutPxY+1:length(pos.y)-plotParams.cutPxY;
    plotXPos = pos.x(plotParams.cutPxX+1:end-plotParams.cutPxX)*10^3;
    plotYPos = pos.y(plotParams.cutPxY+1:end-plotParams.cutPxY)*10^3;

    % Create the vertical line to show the fracture plane on each image
    fractLineX(1:length(plotYPos)) = pos.x(fracture.locX)*10^3;
    fractLineY(1:length(plotYPos)) = plotYPos;

    % Calculate the RAW strains for plotting
    [rawStrain.x,~,~] = gradient(disp.x,pos.xStep,pos.yStep,time.step);

    if strcmp(plotParams.cAxisType,'Specified')
        plotParams.cRange{1} = plotParams.cAxisRawStrain;
        plotParams.cRange{2} = plotParams.cAxisStrain;
        plotParams.cRange{3} = plotParams.cAxisAccel;
        plotParams.cRange{4} = plotParams.cAxisStrainRate;
    else
        plotParams.cRange{1} = func_calcColourBarRange(plotParams.cAxisType,rawStrain.x(plotRangeY,plotRangeX,:))*10^3;
        plotParams.cRange{2} = func_calcColourBarRange(plotParams.cAxisType,strain.x(plotRangeY,plotRangeX,:))*10^3;
        plotParams.cRange{3} = func_calcColourBarRange(plotParams.cAxisType,accel.x(plotRangeY,plotRangeX,:));
        plotParams.cRange{4} = func_calcColourBarRange(plotParams.cAxisType,strainRate.x(plotRangeY,plotRangeX,:));
    end

    %----------------------------------------------------------------------
    % Create and size the figure
    plotProps.sizePerFigXcm = 1.2*plotProps.sizePerFigXcm;
    hf = func_createFigure(plotProps,plotParams);

    % Only plot a few frames after the fracture
    fStart = 1;
    fEnd = fracture.strengthFrame+7;
    if fEnd > time.numFrames
        fEnd = time.numFrames;
    end

    for ff = fStart:fEnd
        %-----------------------------------------------------------------------
        % Stress at the width section [takes up two subplots]
        subplot(plotParams.Rows,plotParams.Cols,[1,4])
        % 1.1) Stress-Gauge: Average
        plot(plotYPos,stress.xAvg(fracture.locX,ff)*ones(1,length(plotYPos))*10^-6,'-k','linewidth',plotProps.lw,'markersize',plotProps.ms)
        hold on
        % 1.2) Stress-Gauge: Linear
        plot(plotYPos,stress.xLinearGauge(plotRangeY,fracture.locX,ff)*10^-6,'--k','linewidth',plotProps.lw,'markersize',plotProps.ms)
        % 2.1) Stress from Strain: Full-Field
        plot(plotYPos,stress.QRecon.x(plotRangeY,fracture.locX,ff)*10^-6,'-ob','linewidth',plotProps.lw,'markersize',plotProps.ms/2)
        % 2.2) Stress from Strain: Average
        stressStrainMean = nanmean(stress.QRecon.x(plotRangeY,fracture.locX,ff)*10^-6);
        plot(plotYPos,stressStrainMean*ones(1,length(plotYPos)),'-b','linewidth',plotProps.lw,'markersize',plotProps.ms)
        hold off

        % Set the axis labels and the plot title using correct formatting
        th = xlabel('Y [mm]');
        set(th,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
        th = ylabel('$\sigma_{xx}$ [MPa]');
        set(th,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
        titleStr = {['(a) ',sprintf('x = %.2fmm, t = %.1f',pos.x(fracture.locX)*10^3,time.vec(ff)*10^6),'$\mu s$']};
        th = title(titleStr);
        set(th,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
        lh = legend({'$\overline{\sigma_{xx}}^{y}(SG)$','$\overline{\sigma_{xx}}^{y}(LSG)$','$\sigma_{xx}(\epsilon)$','$\overline{\sigma_{xx}}^{y}(\epsilon)$'},...
            'location','southeast');
        set(lh,'fontsize',plotProps.lfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText);
        set(gca,'fontsize',plotProps.fs,'fontname',plotProps.ft,'linewidth',plotProps.lw)

        xlim([0 pos.y(end)*10^3])
        if isfield(plotParams,'stressLims')
            ylim(plotParams.stressLims)
        else
            ylim([round(1.5*min(stress.xAvg(fracture.locX,:)),2,'significant')*10^-6,1.5*round(max(stress.xAvg(fracture.locX,:)),2,'significant')*10^-6])
        end
        grid on   

        %--------------------------------------------------------------------------
        % RAW Strain
        subplot(plotParams.Rows,plotParams.Cols,2)
        imagesc(plotXPos,plotYPos,rawStrain.x(plotRangeY,plotRangeX,ff)*10^3)
        hold on
        plot(fractLineX,fractLineY,'--k','linewidth',plotProps.lw)
        hold off
        xlabel('X [mm]','fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText)
        ylabel('Y [mm]','fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText)
        %titleStr = {'(b) $\epsilon_{xx}$ [$mm.m^{-1}$]',[sprintf('t = %.2f',time.vec(ff)*10^6),'$\mu s$']};
        titleStr = ['(b) $\epsilon_{xx}$ [$mm.m^{-1}$], ',sprintf('t = %.2f',time.vec(ff)*10^6),'$\mu s$'];
        title(titleStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText)
        set(gca,'fontsize',plotProps.fs,'fontname',plotProps.ft,'linewidth',plotProps.lw)
        colorbar
        colormap(jet)
        caxis(plotParams.cRange{1})
        axis image

        %--------------------------------------------------------------------------
        % Smoothed Strain
        subplot(plotParams.Rows,plotParams.Cols,3)
        imagesc(plotXPos,plotYPos,strain.x(plotRangeY,plotRangeX,ff)*10^3)
        hold on
        plot(fractLineX,fractLineY,'--k','linewidth',plotProps.lw)
        hold off
        xlabel('X [mm]','fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText)
        ylabel('Y [mm]','fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText)
        %titleStr = {'(c) $\epsilon_{xx}$ [$mm.m^{-1}$]',[sprintf('t = %.2f',time.vec(ff)*10^6),'$\mu s$']};
        titleStr = ['(c) $\epsilon_{xx}$ [$mm.m^{-1}$], ',sprintf('t = %.2f',time.vec(ff)*10^6),'$\mu s$'];
        title(titleStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText)
        set(gca,'fontsize',plotProps.fs,'fontname',plotProps.ft,'linewidth',plotProps.lw)
        colorbar
        colormap(jet)
        caxis(plotParams.cRange{2})
        axis image

        %--------------------------------------------------------------------------
        % Acceleration
        subplot(plotParams.Rows,plotParams.Cols,5)
        imagesc(plotXPos,plotYPos,accel.x(plotRangeY,plotRangeX,ff))
        hold on
        plot(fractLineX,fractLineY,'--k','linewidth',plotProps.lw)
        hold off
        xlabel('X [mm]','fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText)
        ylabel('Y [mm]','fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText)
        %titleStr = {'(d) $a_{x}$ [$m.s^{-2}$]',[sprintf('t = %.2f',time.vec(ff)*10^6),'$\mu s$']};
        titleStr = ['(d) $a_{x}$ [$m.s^{-2}$], ',sprintf('t = %.2f',time.vec(ff)*10^6),'$\mu s$'];
        title(titleStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText)
        set(gca,'fontsize',plotProps.fs,'fontname',plotProps.ft,'linewidth',plotProps.lw)
        colormap(jet)
        colorbar
        caxis(plotParams.cRange{3})
        axis image

        %--------------------------------------------------------------------------
        % Strain Rate
        subplot(plotParams.Rows,plotParams.Cols,6)
        imagesc(plotXPos,plotYPos,strainRate.x(plotRangeY,plotRangeX,ff))
        hold on
        plot(fractLineX,fractLineY,'--k','linewidth',plotProps.lw)
        hold off
        xlabel('X [mm]','fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText)
        ylabel('Y [mm]','fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText)
        %titleStr = {'(e) $\dot{\epsilon_{xx}}$ [$s^{-1}$]',[sprintf('t = %.2f',time.vec(ff)*10^6,'$\mu s$']};
        titleStr = ['(e) $\dot{\epsilon_{xx}}$ [$s^{-1}$], ',sprintf('t = %.2f',time.vec(ff)*10^6),'$\mu s$'];
        title(titleStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText)
        set(gca,'fontsize',plotProps.fs,'fontname',plotProps.ft,'linewidth',plotProps.lw)
        colorbar
        caxis(plotParams.cRange{4})
        axis image

        % Save this frame to file
        saveFile = [savePath,'\Frame_',num2str(ff)];
        func_saveFigureMultFormat(hf,saveFile,plotProps,plotParams)
        clf(hf);
    end
end



