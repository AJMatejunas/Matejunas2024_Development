function [fieldSurfAvg,surfArea] = func_calcRotatedFieldAreaAvgAlongSlice_weighted(...
    slice,specimen,time,pos,inField,weightField)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 11/12/2018
% Date Edited: 11/12/2018

% Prealloc for speed
fieldSurfAvg = zeros(slice.xMaxInd,time.numFrames);
surfArea = zeros(slice.xMaxInd,time.numFrames);

for ff = 1:time.numFrames
    for xx = 1:slice.xMaxInd
        % Calculate the x and y co-ords of the angled slice
        [xSlice,~] = func_calcAngSliceCoords(xx,pos,slice);

        % Create a logical mask for the angled slice - get all points within the
        % bounds
        sliceMask = zeros(length(pos.y),length(pos.x));
        for ii = 1:length(pos.y)
           sliceMask(ii,:) = pos.x < xSlice(ii);  
        end
        sliceMask = logical(sliceMask);

        % Calculate the acceleration field average from the masked field
        fieldTemp = inField(:,:,ff).*weightField{xx};
        
        % Set the initial frame to nan to ignore values with nanmean
        fieldFrame = nan(length(pos.y),length(pos.x));
        fieldFrame(sliceMask) = fieldTemp(sliceMask);
        
        % Calculate the average of the masked field
        fieldSurfAvg(xx,ff) = nanmean(fieldFrame(:));

        % Calculate the area of the 'sliced' surface
        surfArea(xx,ff) = specimen.height*(pos.x(xx)+slice.lengthX/2);
    end
end

end

