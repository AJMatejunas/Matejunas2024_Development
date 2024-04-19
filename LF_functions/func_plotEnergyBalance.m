function hf = func_plotEnergyBalance(plotParams,EBalOpts,energy,time)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 3/10/2017
% Plot the energy balance figure

    % Based on the options structure create the title and legend strings
    sigFigs = 3;
    if EBalOpts.calcExtWorkFromAccel
        titleStr = {['Error_{max}(Strain): ',num2str(energy.balStrainErrPc,sigFigs),'%'],...
                    ['Error_{max}(Accel): ',num2str(energy.balAccelErrPc,sigFigs),'%']};
        legStrs = {'Impact(Strain)','Impact(Accel)','Kinetic','Strain',...
            'Balance(Strain)','Balance(Accel)','KE+SE'};
    else
        titleStr = {['Error_{max}(Strain):',num2str(energy.balStrainErrPc,sigFigs),'%']};
        legStrs = {'Impact(Strain)','Kinetic','Strain',...
            'Balance(Strain)','KE+SE'};
    end
    
    % Get some plot parameters for nice formatting
    plotProps = func_initPlotPropsStruct(plotParams.formatType);
    
    % Create and size the figure
    hf = func_createFigure(plotParams,plotProps);

    % Compare the impact energy measures
    hold on
    plot(time.vec.*10^6,energy.impactStrain,'-k','linewidth',plotProps.lw,'markersize',plotProps.ms)
    if EBalOpts.calcExtWorkFromAccel
        plot(time.vec.*10^6,energy.impactAccel,'--k','linewidth',plotProps.lw,'markersize',plotProps.ms)
    end
    plot(time.vec.*10^6,energy.kinetic,'-g','linewidth',plotProps.lw,'markersize',plotProps.ms)
    plot(time.vec.*10^6,energy.strain,'-b','linewidth',plotProps.lw,'markersize',plotProps.ms)
    plot(time.vec.*10^6,energy.balStrain,'-r','linewidth',plotProps.lw,'markersize',plotProps.ms)
    if EBalOpts.calcExtWorkFromAccel
        plot(time.vec.*10^6,energy.balAccel,'--r','linewidth',plotProps.lw,'markersize',plotProps.ms)
    end
    plot(time.vec.*10^6,energy.kineticPlusStrain,'-.k','linewidth',plotProps.lw,'markersize',plotProps.ms)
    hold off
     
    title(titleStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft)
    xlabel('Time, \mus','fontsize',plotProps.hfs,'fontname',plotProps.ft)
    ylabel('Energy, J','fontsize',plotProps.hfs,'fontname',plotProps.ft)
    if plotParams.specifyAxisLims
        xlim(plotParams.xLim)
        ylim(plotParams.yLim)
    end
    legend(legStrs,'location','eastoutside')
    set(gca, 'fontsize', plotProps.fs,'fontname',plotProps.ft)
    set(gca,'XMinorTick','on','YMinorTick','on')
    box on

end

