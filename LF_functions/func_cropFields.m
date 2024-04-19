function[cropInputVar] = func_cropFields(inputVar,extrapDim,extrapPixels)
% Author: Jared Van Blitterswyk
% PhotoDyn Group, University of Southampton
% Date: 19/9/2017
% Crop poor data around edges of displacement fields prior to spatial
% smoothing to avoid spatial leakage of poor data

% Delete the edge pixels and assign the position variable based on the
% dimension we are extrapolating
if extrapDim == 2 % crop only left and right edges
    cropInputVar = inputVar(:,extrapPixels+1:end-extrapPixels,:);
else % crop only top and bottom edges
    cropInputVar = inputVar(extrapPixels+1:end-extrapPixels,:,:);     
end

