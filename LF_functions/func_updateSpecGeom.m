function [specimen,grid] = func_updateSpecGeom(specimen,grid,disp)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 24/1/2017
% Date Edited: 24/1/2019
%
% Updates the specimen geometry in the specimen and grid structs based on
% the selected FOV from the grid method image processing.

% Check that the nominal specimen dimensions and the dimensions of the
% selected grid window are the same, if they are not update them
if isfield(grid,'asymmPitch')
if grid.asymmPitch
    checkLength = grid.mPerPxX*size(disp.x,2);
    if checkLength ~= specimen.length
        specimen.length = checkLength;
        grid.length = checkLength;
    end
    checkHeight = grid.mPerPxY*size(disp.x,1);
    if checkHeight ~= specimen.height
        specimen.height = checkHeight;
        grid.height = checkHeight;
    end
else
    checkLength = grid.mPerPx*size(disp.x,2);
    if checkLength ~= specimen.length
        specimen.length = checkLength;
        grid.length = checkLength;
    end
    checkHeight = grid.mPerPx*size(disp.x,1);
    if checkHeight ~= specimen.height
        specimen.height = checkHeight;
        grid.height = checkHeight;
    end
end
else
    checkLength = grid.mPerPx*size(disp.x,2);
    if checkLength ~= specimen.length
        specimen.length = checkLength;
        grid.length = checkLength;
    end
    checkHeight = grid.mPerPx*size(disp.x,1);
    if checkHeight ~= specimen.height
        specimen.height = checkHeight;
        grid.height = checkHeight;
    end
end

end

