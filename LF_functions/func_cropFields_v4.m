function[cropInputVar] = func_cropFields_v4(inputVar,cropDim,cropPxX,cropPxY)
% Author: Jared Van Blitterswyk
% PhotoDyn Group, University of Southampton
% Date Created: 19/9/2017
% Date Edited: 27/7/2020 - LF, added option to crop both edges
% Dated Edited: 28/7/2020 - LF, addded option to have asymmetric crop
%
% Crop poor data around edges of displacement fields prior to spatial
% smoothing to avoid spatial leakage of poor data

    % Delete the edge pixels and assign the position variable based on the
    % dimension we are extrapolating
    if strcmp(cropDim,'X') % crop only left and right edges
        cropInputVar = inputVar(:,cropPxX+1:end-cropPxX,:);   
    elseif strcmp(cropDim,'Y') % crop only top and bottom edges
        cropInputVar = inputVar(cropPxY+1:end-cropPxY,:,:);   
    else
        cropInputVar = inputVar(cropPxY+1:end-cropPxY,...
            cropPxX+1:end-cropPxX,:);  
    end

end

