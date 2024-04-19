function [hf,Q] = func_plotIdentQvsLength(labelStrs,plotParams,time,grid,material,stressVars,strainVars,targPc)
% Fits a linear function to stress gauge data to calculate Q
    
    % Grab the plot properties for nice formatting
    plotProps = func_initPlotPropsStruct();

    % Make sure there are the same number of stress and strain vars
    check = [length(stressVars),length(strainVars)];
    if range(check) ~= 0;
        error('The lengths of the stress and strain variables must be the same.')
    end
    numVars = check(1);
    
    % Specify algorithm for linear fit
    options = fitoptions('Method','LinearLeastSquares','Robust','Off');
    % Loop through and fit stress strain curves
    for v = 1:numVars
       sizeVec(v) = size(strainVars{v}.xnyAvg,1);
    end  
    numXvals = min(sizeVec);
    for x = 1:numXvals
       for v = 1:numVars 
           % Fit Q for the image deformation data
           [Q{v}.xxFit{x},Q{v}.xxgof{x}] = fit(strainVars{v}.xnyAvg(x,1:end-time.cutFrames)',...
               stressVars{v}.xAvg(x,:)','poly1',options); 
           Q{v}.xx(x) = Q{v}.xxFit{x}.p1;
       end
    end

    % Calculate target value for Q and set axis lims based on this
    target = [];
    upperLim = [];
    lowerLim = [];
    target(1:numXvals) = material.Qxx;
    upperLim(1:numXvals) = material.Qxx*(1+targPc);
    lowerLim(1:numXvals) = material.Qxx*(1-targPc);
    xVec = grid.mmPerPx*(0:numXvals-1)*10^3;
    axisLims = [0,max(xVec),material.Qxx*(1-targPc-targPc/2),material.Qxx*(1+targPc+targPc/2)];
    
    % Create the cell array for the legend entry including the target value
    legendStrs = labelStrs.legStrs;
    legendStrs{end+1} = 'Target';
    legendStrs{end+1} = ['Target \pm',num2str(targPc*100),'%'];
    
    % Setup and plot the figure
    hf = figure;
    set(hf,'Position', [plotProps.locX,plotProps.locY,round(plotProps.ratio*plotProps.size),plotProps.size])
    set(hf,'PaperPositionMode','auto')
    hold on
    
    for v = 1:numVars 
        plot(xVec,Q{v}.xx,plotParams.markStrs{v},'linewidth',plotProps.lw,'markersize',plotProps.ms)
    end
    
    plot(xVec,target,'--k','linewidth',plotProps.lw,'markersize',plotProps.ms)
    plot(xVec,upperLim,'-k',xVec,lowerLim,'-k','linewidth',plotProps.lw,'markersize',plotProps.ms)
    xlabel(labelStrs.x,'fontsize',plotProps.hfs,'fontname',plotProps.ft)
    ylabel(labelStrs.y,'fontsize',plotProps.hfs,'fontname',plotProps.ft)    
    legend(legendStrs)
    axis(axisLims)
    set(gca,'fontsize', plotProps.fs,'fontname',plotProps.ft)
    set(gca,'XMinorTick','on','YMinorTick','on')
    box on
    hold off
end

