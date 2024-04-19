function [ xSlice , ySlice ] = func_calcAngSliceCoords(xInd,pos,slice)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 12/2/2018
% Date Edited: 8/12/2018 

% Calculate the X co-ords of the angled line - Y co-ords got to spec height
ySlice = pos.y;
if (slice.angle < 0 && slice.angle > -90)|| (slice.angle > 90 && slice.angle < 180)
    xSlice = 1/tand(slice.angle)*pos.y+pos.x(xInd)+slice.lengthX;
else
    xSlice = 1/tand(slice.angle)*pos.y+pos.x(xInd);
end

end

