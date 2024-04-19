function [strainX,strainY,strainS] = ...
    func_calcStrainFromDisp(dispX,dispY,xStep,yStep)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 9/2/2017
% Date Edited: 20/12/2019
%
% Calculates the strain field from the 2D displacement field

    [dux_dx,dux_dy,~] = gradient(dispX,xStep,yStep,1);
    [duy_dx,duy_dy,~] = gradient(dispY,xStep,yStep,1);
    strainX = dux_dx;
    strainY = duy_dy;
    % Engineering shear strain!
    strainS = dux_dy + duy_dx;
end

