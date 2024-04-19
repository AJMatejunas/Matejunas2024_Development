function func_plotIdentStrengthVideo_v2(plotParams,savePath,pos,time,disp,accel,...
    strain,strainRate,stress,fracture,stiffness)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 29/3/2017
%
% Plots image of fracture location and strength identification

% Check if the image sequence path exits, if it does ask the user if it
% needs to be plotted again
if exist(savePath,'file') == 7
    choice = questdlg('Strength identification image sequence folder found, plot again?', ...
    'Plot Image Sequence?', 'Yes','No','No');
    switch choice
        case 'Yes'
            plotImageSeq = true;
        case 'No'
            plotImageSeq = false;
    end
else
    % If the directory does not exist, create it and plot the image seq
    mkdir(savePath);
    plotImageSeq = true;
end
      
% Only plot the image sequence if it doesn't exist in file 
if plotImageSeq
    % Reconstruct the stress from the strains and identified stiffness
    if stiffness.method == 1
        stress.xQRecon = stiffness.Exx*strain.x;
    else
        stress.xQRecon = stiffness.Exx/(1-stiffness.nuxy^2)*(strain.x+stiffness.nuxy.*strain.y);
    end
    
    % Create a struct of formatting properties for the figure
    plotProps = func_initPlotPropsStruct(plotParams.formatType);
    plotParams.cAxisType = 'MaxQuantile';
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
    xStep = abs(pos.x(2) - pos.x(1));
    yStep = abs(pos.y(2) - pos.y(1));
    [rawStrain.x,~,~] = gradient(disp.x,xStep,yStep,time.step);

    % Create and size the figure
    hf = figure('Visible','On');
    set(hf,'color',[1 1 1])
    set(hf,'Position', [plotProps.locX,plotProps.locY,...
         1.5*plotParams.Cols*plotProps.sizePerFigX,1.5*plotParams.Rows*plotProps.sizePerFigY])
    set(hf,'PaperPositionMode','auto')

    calcType = 'MaxQuantile';
    %plotParams.cRange{1} =
    %func_calcColourBarRange(calcType,disp.xDef(plotRangeY,plotRangeX,:));
    plotParams.cRange{1} = func_calcColourBarRange(calcType,accel.x(plotRangeY,plotRangeX,:));
    plotParams.cRange{2} = func_calcColourBarRange(calcType,rawStrain.x(plotRangeY,plotRangeX,:));
    plotParams.cRange{3} = func_calcColourBarRange(calcType,strainRate.x(plotRangeY,plotRangeX,:));

    for f = 1:time.numFrames
        %--------------------------------------------------------------------------
        % Stress at the width section [takes up two subplots]
        subplot(plotParams.Rows,plotParams.Cols,[1,4])
        % 1.1) Stress-Gauge: Average
        plot(plotYPos,stress.xAvg(fracture.locX,f)*ones(1,length(pos.y))*10^-6,'-k','linewidth',plotProps.lw,'markersize',plotProps.ms)
        hold on
        % 1.2) Stress-Gauge: Linear
        plot(plotYPos,stress.xLinearGauge(plotRangeY,fracture.locX,f)*10^-6,'--k','linewidth',plotProps.lw,'markersize',plotProps.ms)
        
        % 2.1) Stress from Strain: Full-Field
        plot(plotYPos,stress.xQRecon(plotRangeY,fracture.locX,f)*10^-6,'-ob','linewidth',plotProps.lw,'markersize',plotProps.ms/2)
        
        % 2.2) Stress from Strain: Average
        stressStrainMean = mean(stress.xQRecon(plotRangeY,fracture.locX,f)*10^-6);
        plot(plotYPos,stressStrainMean*ones(1,length(pos.y)),'--b','linewidth',plotProps.lw,'markersize',plotProps.ms)
        hold off
        lh = legend({'SG','SG,L','\sigma(\epsilon)','\sigma(\epsilon)_{avg}'},'location','southeast');
        xlabel('Y (mm)','fontsize',plotProps.hfs,'fontname',plotProps.ft)
        xlim([0 pos.y(end)*10^3])
        ylabel('\sigma_{x} (MPa)','fontsize',plotProps.hfs,'fontname',plotProps.ft)
        if isfield(plotParams,'stressLims')
            ylim(plotParams.stressLims)
        else
            ylim([round(1.5*min(stress.xAvg(fracture.locX,:)),2,'significant')*10^-6,1.5*round(max(stress.xAvg(fracture.locX,:)),2,'significant')*10^-6])
        end
        titleStr = {[sprintf('x = %.2fmm, t = %.1f',pos.x(fracture.locX)*10^3,time.vec(f)*10^6),'\mus']};
        title(titleStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft)
        set(gca,'fontsize',plotProps.fs,'fontname',plotProps.ft,'linewidth',plotProps.lw)
        set(gca,'fontsize',plotProps.fs,'fontname',plotProps.ft,'linewidth',plotProps.lw)
        set(lh,'fontsize',round(plotProps.fs*0.75))
        grid on   
        
        %--------------------------------------------------------------------------
        % RAW Strain
        subplot(plotParams.Rows,plotParams.Cols,2)
        imagesc(plotXPos,plotYPos,rawStrain.x(plotRangeY,plotRangeX,f)*10^3)
        hold on
        plot(fractLineX,fractLineY,'--k','linewidth',plotProps.lw)
        hold off
        xlabel('X (mm)','fontsize',plotProps.hfs,'fontname',plotProps.ft)
        ylabel('Y (mm)','fontsize',plotProps.hfs,'fontname',plotProps.ft)
        titleStr = {'\it{\epsilon_{x}}\rm [mm.m^{-1}]',[sprintf('t = %.2f',time.vec(f)*10^6),'\mus']};
        title(titleStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft)
        set(gca,'fontsize',plotProps.fs,'fontname',plotProps.ft,'linewidth',plotProps.lw)
        colorbar
        colormap(jet)
        caxis(plotParams.cRange{2}*10^3)
        
        %--------------------------------------------------------------------------
        % Smoothed Strain
        subplot(plotParams.Rows,plotParams.Cols,3)
        imagesc(plotXPos,plotYPos,strain.x(plotRangeY,plotRangeX,f)*10^3)
        hold on
        plot(fractLineX,fractLineY,'--k','linewidth',plotProps.lw)
        hold off
        xlabel('X (mm)','fontsize',plotProps.hfs,'fontname',plotProps.ft)
        ylabel('Y (mm)','fontsize',plotProps.hfs,'fontname',plotProps.ft)
        titleStr = {'\it{\epsilon_{x}}\rm [mm.m^{-1}]',[sprintf('t = %.2f',time.vec(f)*10^6),'\mus']};
        title(titleStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft)
        set(gca,'fontsize',plotProps.fs,'fontname',plotProps.ft,'linewidth',plotProps.lw)
        colorbar
        colormap(jet)
        caxis(plotParams.cRange{2}*10^3)
        
        %--------------------------------------------------------------------------
        % Acceleration
        subplot(plotParams.Rows,plotParams.Cols,5)
        imagesc(plotXPos,plotYPos,accel.x(plotRangeY,plotRangeX,f))
        hold on
        plot(fractLineX,fractLineY,'--k','linewidth',plotProps.lw)
        hold off
        xlabel('X (mm)','fontsize',plotProps.hfs,'fontname',plotProps.ft)
        ylabel('Y (mm)','fontsize',plotProps.hfs,'fontname',plotProps.ft)
        titleStr = {'\it{a_{x}}\rm [m.s^{-2}]',[sprintf('t = %.2f',time.vec(f)*10^6),'\mus']};
        title(titleStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft)
        set(gca,'fontsize',plotProps.fs,'fontname',plotProps.ft,'linewidth',plotProps.lw)
        colormap(jet)
        colorbar
        caxis(plotParams.cRange{1})

        %--------------------------------------------------------------------------
        % Strain Rate
        subplot(plotParams.Rows,plotParams.Cols,6)
        imagesc(plotXPos,plotYPos,strainRate.x(plotRangeY,plotRangeX,f))
        hold on
        plot(fractLineX,fractLineY,'--k','linewidth',plotProps.lw)
        hold off
        xlabel('X (mm)','fontsize',plotProps.hfs,'fontname',plotProps.ft)
        ylabel('Y (mm)','fontsize',plotProps.hfs,'fontname',plotProps.ft)
        titleStr = {'\it{d\epsilon_{x}/dt}\rm [s^{-1}]',[sprintf('t = %.2f',time.vec(f)*10^6),'\mus']};
        title(titleStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft)
        set(gca,'fontsize',plotProps.fs,'fontname',plotProps.ft,'linewidth',plotProps.lw)
        colorbar
        caxis(plotParams.cRange{3})

        % Save this frame to file
        print(hf,[savePath,'\StrengthIdentification_Frame',num2str(f)],plotProps.format,'-r0')
        clf(hf);
    end
end

end

