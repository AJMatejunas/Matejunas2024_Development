function [extrapVar] = func_spatialExtrapWidthAvgVar(inputVar,pos,extrapPixels,extrapMethod)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 27/3/2017
% Extrapolate width averaged variable along its length for a given number
% of pixels.

    % Delete the edge pixels
    crop_posX = pos.x(extrapPixels+1:end-extrapPixels);
    crop_inputVar = inputVar(extrapPixels+1:end-extrapPixels,:);
    extrapVar = inputVar;
     
    for i = 1:size(inputVar,2)
        % Extrapolate the variable over the given number of pixels
        temp_extrapVar(:,i) = interp1(crop_posX,crop_inputVar(:,i),pos.x,extrapMethod,'extrap');
        % Add the extrapolated data back onto the edges
        extrapVar(1:extrapPixels,i) = temp_extrapVar(1:extrapPixels,i);
        extrapVar(end-extrapPixels+1:end,i) = temp_extrapVar(end-extrapPixels+1:end,i);
    end

end

