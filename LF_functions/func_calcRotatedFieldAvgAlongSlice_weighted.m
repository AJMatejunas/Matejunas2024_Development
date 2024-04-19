function [sliceAvg,sliceField] = func_calcRotatedFieldAvgAlongSlice_weighted(...
    slice,pos,inField,weightVec)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 9/2/2018
% Date Edited: 14/12/2018
%
% Calculates the field average along an angled slice at material.rotAngle 
% uses a scattered interpolant field to generate the field along the 
% profile. 

% Calculate the number of points in the field for vectorisation
nPts = length(pos.x)*length(pos.y);
[sy,~,numFrames] = size(inField);

% Pre-alloc for speed
sliceAvg = zeros(slice.xMaxInd,numFrames);
sliceField = zeros(sy,slice.xMaxInd,numFrames);

for ff = 1:numFrames
    % Create the scattered interpolants for the material co-ord field for
    % this frame
    fieldInterp = scatteredInterpolant(reshape(pos.xGrid,nPts,1),...
                reshape(pos.yGrid,nPts,1),reshape(inField(:,:,ff),nPts,1));
 
    for xx = 1:slice.xMaxInd
        % Calculate the x and y co-ords of the angled slice
        [xSlice,ySlice] = func_calcAngSliceCoords(xx,pos,slice);
        
        temp = fieldInterp(xSlice,ySlice);
        sliceAvg(xx,ff) = nanmean(squeeze(temp).*weightVec);
        sliceField(:,xx,ff) = temp.*weightVec;
    end
end

end

