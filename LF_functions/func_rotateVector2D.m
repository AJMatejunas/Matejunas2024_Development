function [xRot,yRot] = func_rotateVector2D(xComp,yComp,alpha)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 8/2/2017

xRot = cosd(alpha).*xComp + sind(alpha).*yComp;
yRot = -sind(alpha).*xComp + cosd(alpha).*yComp;

end

