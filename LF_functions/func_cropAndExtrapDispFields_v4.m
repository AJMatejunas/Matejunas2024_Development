function [dispExtrapX,dispExtrapY,dispRangeX,dispRangeY] = ...
    func_cropAndExtrapDispFields_v4(pos,dispX,dispY,extrapOptsD,printToCons)
% Crops bad edge data from the grid method. The missing edge data is then
% extrapolated using the options specified in the data structure
% 'extrapOpts'. The final step is to check if there are any NaNs left in
% the fields from the iterative grid calculation (method 2) and replace
% this with a local average over a specified window.
%
% Author: Lloyd Fletcher
% Image-Based Mechanics Group (IBMG), University of Southampton
% Date Created: 28/7/2020 - v4 Major Update! Fully overhauled kinematic
% calculation code to include better edge extrapolation options.
% Date Edited: 28/7/2020

if nargin < 5
    % Print to console as default
    printToCons = true;
end
debug = false;

% NOTE: function inputs are specified to avoid issues with parfor, here we
% change back to struct form for readability
disp.x = dispX;
disp.y = dispY;
clear dispX dispY

% Store the original size of the displacement fields before we crop and
% extrapolate. Use this later to create a mask to cut extra extrapolation.
[origSizeY,origSizeX,~] = size(disp.x);

%--------------------------------------------------------------------------
% 1) Crop specified number of pixels of specified edges
if printToCons
    fprintf('\tCropping displacement data from edges to remove poor data from grid method processing.\n')
end
[dispCrop,posCropX,posCropY] = func_cropDispFields(extrapOptsD,pos,disp);

%--------------------------------------------------------------------------
% 2) Extrapolate displacement field back to original size
dispExtrap = func_extrapolateDispFields(...
    extrapOptsD,pos,dispCrop,posCropX,posCropY,debug);

%--------------------------------------------------------------------------
% 3) Replace NaNs with local average
% NOTE: using the grid method with 'method 2' (iterative calculation) can 
% pull some NaNs into the fields on the edges.
if sum(isnan(dispExtrap.x(:))) > 0
    if printToCons
        fprintf('\tWARNING: NaNs found in x displacement field after extrapolation!\n')
        fprintf('\tReplacing NaNs with local average over a distance of %i px.\n',extrapOptsD.fixNaNKernel)
    end
    dispExtrap.x = func_fixNaNsInDispX(extrapOptsD,dispExtrap.x);
    if printToCons
        if sum(isnan(dispExtrap.x(:))) > 0
            fprintf('\tWARNING: NaNs still in x displacement field, widen kernel.\n')
        end
    end
end

if sum(isnan(dispExtrap.y(:))) > 0
    if printToCons
        fprintf('\tWARNING: NaNs found in y displacement field after extrapolation!\n')
        fprintf('\tReplacing NaNs with local average over a distance of %i px.\n',extrapOptsD.fixNaNKernel)
    end
    dispExtrap.y = func_fixNaNsInDispY(extrapOptsD,dispExtrap.y);
    if printToCons
        if sum(isnan(dispExtrap.y(:))) > 0
            fprintf('\tWARNING: NaNs still in y displacement field, widen kernel.\n')
        end
    end
end

%--------------------------------------------------------------------------
% 4) Create a mask for selecting displacement within original sample area
% NOTE: this is used to crop back to original size if we extend the
% extrapolation region as a kind of 'padding' for later spatial filters
[extSizeY,extSizeX,~] = size(dispExtrap.x);
edgePxY = (extSizeY - origSizeY)/2;
edgePxX = (extSizeX - origSizeX)/2;
dispRangeX = 1+edgePxX:extSizeX-edgePxX;
dispRangeY = 1+edgePxY:extSizeY-edgePxY;

if debug
   figure; imagesc(flipud(dispExtrap.x(dispRangeY,dispRangeX,end)));
   title('dispExtrap.x, masked');
   colorbar; axis image; set(gca,'YDir','normal'); colormap jet;
end

%--------------------------------------------------------------------------
% 5) Change back to non-struct form to return vars
dispExtrapX = dispExtrap.x;
dispExtrapY = dispExtrap.y;
end

