function plotProps = func_initPlotPropsStruct(typeStr)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 19/2/2018

    if nargin == 0
        typeStr = 'default';
    end
    
    plotProps.screenSize = get(groot,'Screensize');
    plotProps.interpText = 'latex';
    plotProps.units = 'centimeters';
    plotProps.ratio = 1.618; 
    plotProps.pageSizeXcm = 18;
    plotProps.singleColFigFactor = 1.5;
    plotProps.imageSeqSaveRes = '-r300';
    plotProps.imageSeqFormat = '-dpng';
    plotProps.titleAlign = 'center';
    plotProps.black = [0.2,0.2,0.2];
    plotProps.a4PageXcm = 21;
    plotProps.a4PageYcm = 29.7;
    plotProps.marginXcm = 2.5;
    plotProps.marginYcm = 2.5;
    plotProps.a4UsableXcm = plotProps.a4PageXcm - 2*plotProps.marginXcm;
    plotProps.a4UsableYcm = plotProps.a4PageYcm - 2*plotProps.marginYcm;
    %,'color',plotProps.black
    
    % Creates a data struct for plot properties
    if strcmp(typeStr,'article')
        plotProps.format = '-dpng';       % File format for image saving
        plotProps.saveRes = '-r300';
        
        plotProps.locXcm = 0;            % Image location in cm
        plotProps.locYcm = 0;            % Image location in cm

        plotProps.sizePerFigXcm = 9;
        plotProps.sizePerFigYcm = plotProps.sizePerFigXcm/plotProps.ratio;
        
        plotProps.lw = 0.6;             % plot line width 
        plotProps.ms = 4;               % marker size for plots 
        
        plotProps.ft = 'times';         % font for figures 
        plotProps.fs = 9;
        plotProps.hfs = 10;
        plotProps.lfs = 9;
        plotProps.txtFs = 8;
        
        plotProps.colormap = jet;
        
        plotProps.lLoc = 'eastoutside';
    elseif strcmp(typeStr,'article_v2')
        plotProps.format = '-dpng';       % File format for image saving
        plotProps.saveRes = '-r300';
        
        plotProps.locXcm = 0;            % Image location in cm
        plotProps.locYcm = 0;            % Image location in cm

        plotProps.sizePerFigXcm = 10;
        plotProps.sizePerFigYcm = plotProps.sizePerFigXcm/2;
        
        plotProps.lw = 0.5;             % plot line width 
        plotProps.ms = 4;               % marker size for plots 
        
        plotProps.ft = 'times';         % font for figures 
        plotProps.fs = 9;
        plotProps.hfs = 10;
        plotProps.lfs = 9;
        plotProps.txtFs = 8;
        
        plotProps.colormap = jet;
        
        plotProps.lLoc = 'eastoutside';
    
    elseif strcmp(typeStr,'presentation')
        plotProps.format = '-dpng';       % File format for image saving
        plotProps.saveRes = '-r300';
        plotProps.imageSeqFormat = '-dpng';
        plotProps.imageSeqSaveRes = '-r300';
        
        plotProps.locXcm = 0;            % Image location in cm
        plotProps.locYcm = 0;            % Image location in cm

        plotProps.sizePerFigXcm = 10;
        plotProps.sizePerFigYcm = plotProps.sizePerFigXcm/2;
        
        plotProps.lw = 0.6;             % plot line width 
        plotProps.ms = 5;               % marker size for plots 
        
        plotProps.ft = 'times';         % font for figures 
        plotProps.fs = 9;
        plotProps.hfs = 10;
        plotProps.lfs = 9;
        plotProps.txtFs = 8;
        
        plotProps.colormap = jet;
        
        plotProps.lLoc = 'eastoutside';    
         
    elseif strcmp(typeStr,'poster')
        plotProps.format = '-dmeta';       % File format for image saving
        plotProps.imageSeqFormat = '-dmeta';
        plotProps.saveRes = '-r300';
        plotProps.singleColFigFactor = 1.3;
        
        plotProps.locXcm = 0;            % Image location in cm
        plotProps.locYcm = 0;            % Image location in cm

        plotProps.sizePerFigXcm = 9;
        plotProps.sizePerFigYcm = 6;
        
        plotProps.lw = 0.6;             % plot line width 
        plotProps.ms = 4;               % marker size for plots 
        
        plotProps.ft = 'times';         % font for figures 
        plotProps.fs = 10;
        plotProps.hfs = 11;
        plotProps.lfs = 10;
        plotProps.txtFs = 9;
        
        plotProps.colormap = jet;
        
        plotProps.lLoc = 'eastoutside';      
    else
        plotProps.format = '-dpng';       % File format for image saving
        plotProps.saveRes = '-r150';
        
        plotProps.locXcm = 0;            % Image location in cm
        plotProps.locYcm = 0;            % Image location in cm

        plotProps.sizePerFigXcm = 9;
        plotProps.sizePerFigYcm = plotProps.sizePerFigXcm/plotProps.ratio;
        
        plotProps.lw = 1;             % plot line width 
        plotProps.ms = 3;               % marker size for plots 
        
        plotProps.ft = 'times';         % font for figures 
        plotProps.fs = 10;
        plotProps.hfs = 11;
        plotProps.lfs = 10;
        plotProps.txtFs = 9;
        
        plotProps.colormap = jet;
        
        plotProps.lLoc = 'eastoutside';
    end
    
    plotProps.markerVec = {'ob','or','ok','og','ob','or','ok','og','ob','or','ok','og'};
    plotProps.letter = {'a','b','c','d','e','f','g','h','i','j','k','l','m','n'};
    plotProps.colour = {'b','r','g','k'};
    plotProps.marker = {'x','+','o','s','d','^'};
    plotProps.line = {'-','--','-.',':'};
    plotProps.singleLS = '-xb';
    plotProps.defColourMat = get(groot,'DefaultAxesColorOrder');
    
    %----------------------------------------------------------------------
    % Create a vector of line styles using only line, marker and colour
    ii = 1;
    for ll = 1:length(plotProps.line)
        for mm = 1:length(plotProps.marker)
            for cc = 1:length(plotProps.colour)
                plotProps.lineStyleVec{ii} =...
                    [plotProps.line{ll},plotProps.marker{mm},plotProps.colour{cc}];
                ii = ii+1;
            end
        end
    end
    
    %----------------------------------------------------------------------
    % Create a vector of line styles using only line and marker specifiers
    ii = 1;
    for ll = 1:length(plotProps.line)
        for mm = 1:length(plotProps.marker)
            plotProps.lineStyleVecLineMarkOnly{ii} =...
                [plotProps.line{ll},plotProps.marker{mm}];
            ii = ii+1;
        end
    end
    
end

