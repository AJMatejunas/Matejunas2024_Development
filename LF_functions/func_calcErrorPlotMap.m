function [hf,err] = func_calcErrorPlotMap(labelStrs,errRange,compVar,baseVar)
% Plots error map in space and time between two variables

    % Grab the plot properties for nice formatting
    plotProps = func_initPlotPropsStruct();
    
    % Calculate the errors in time and space
    % RAW error
    err.map = (compVar(errRange{1},errRange{2}) - baseVar(errRange{1},errRange{2}));
    % Percent Error
    err.mapPc = (err.map./max(max(abs(baseVar(errRange{1},errRange{2})))))*100;
    % Max Error
    err.absMax = sqrt(max(max(err.map.^2)));
    err.absMaxPc = sqrt(max(max(err.mapPc.^2)));
    % Median Error
    err.absMed = median(abs(err.map(:)));
    err.absMedPc = median(abs(err.mapPc(:)));
    % Max error for each slice
    err.absMaxVsL = max(abs(err.map),[],2);
    err.absMaxVsLPc = max(abs(err.mapPc),[],2);
    
    % Setup and plot the error map figures
    hf = figure;
    set(hf,'Position', [plotProps.locX,plotProps.locY,plotProps.sizeX*2,plotProps.sizeY])
    set(hf,'PaperPositionMode','auto')

    % Plot the raw error map
    subplot(1,2,1)
    hold on
    imagesc(err.map);
    title(labelStrs.t,'fontsize',plotProps.hfs,'fontname',plotProps.ft)
    xlabel(labelStrs.x,'fontsize',plotProps.hfs,'fontname',plotProps.ft)
    ylabel(labelStrs.y,'fontsize',plotProps.hfs,'fontname',plotProps.ft)
    set(gca, 'fontsize', plotProps.fs,'fontname',plotProps.ft)
    set(gca,'XMinorTick','on','YMinorTick','on')
    cbh = colorbar;
    set(cbh, 'fontsize', plotProps.fs,'fontname',plotProps.ft)
    colormap(jet)
    box on
    axis([1,size(err.map,2),1,size(err.map,1)])
    hold off 

    % Plot the percentage error
    subplot(1,2,2)
    hold on
    imagesc(err.mapPc);
    title(labelStrs.t,'fontsize',plotProps.hfs,'fontname',plotProps.ft)
    xlabel(labelStrs.x,'fontsize',plotProps.hfs,'fontname',plotProps.ft)
    ylabel(labelStrs.y,'fontsize',plotProps.hfs,'fontname',plotProps.ft)
    cbh = colorbar;
    set(cbh, 'fontsize', plotProps.fs,'fontname',plotProps.ft)
    colormap(jet)
    box on
    axis([1,size(err.map,2),1,size(err.map,1)])
    hold off 

end

