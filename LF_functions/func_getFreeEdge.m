function [freeEdge,specimen,disp] = func_getFreeEdge(hardCodeFlag,imagePath,...
    imageFile,specimen,disp,printToCons)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 8/12/2017
% Date Edited: 3/5/2019
%
% Determines which side of the reference image contains the free edge of
% the specimen for stress gauge processing
        
    % Allow printing to the console by default
    if nargin < 6
        printToCons = true;
    end

    % Create the plot properties structure
    plotProps = func_initPlotPropsStruct();
    
    if ~hardCodeFlag
        % Load the reference image into memory
        refImage = imread([imagePath,imageFile]);

        % Show the first image in the sequence
        fh = figure;
        set(fh,'Position', [plotProps.locX,plotProps.locY,...
            2*plotProps.sizePerFigX,2*plotProps.sizePerFigY])
        set(fh,'PaperPositionMode','auto')
        imshow(refImage, []);

        % Ask the user which side the free edge is on using question dialogue
        freeEdge = questdlg('Which edge is the free edge on?', ...
            'Select free edge','Left','Right','Left');

        % Close the image
        clf(fh);
        close(fh);
        specimen.freeEdge = freeEdge;
    else
        freeEdge = specimen.freeEdge;
    end

    if strcmp(freeEdge,'Right')
        % If the free edge is on the right hand side flip the displacement
        % matrix left to right
        if printToCons
            fprintf('Free edge is on the right hand side, flipping displacement matrices\n')
        end
        for f = 1:size(disp.x,3)
            disp.x(:,:,f) = -fliplr(disp.x(:,:,f));
            disp.y(:,:,f) = fliplr(disp.y(:,:,f));
        end
    elseif strcmp(freeEdge,'Left')
        if printToCons
            fprintf('Free edge is on the left hand side, no correction required.\n')
        end
    end
end