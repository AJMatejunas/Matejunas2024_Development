function [stress22Avg,stress12Avg] = func_stressGaugeProcess_Ang(slice,material,specimen,time,pos,accel)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 9/2/2018
% Date Edited: 13/7/2018
%
% Calculates stress averages along an angled slice at the material co-ord
% system angle. Angle is specified in the slice.angle parameter.
% NOTE: acceleration field needs to be in material co-ords. The structure
% should have the fields .mat11 and .mat22 in the material co-ords. Use the
% function func_rotateVector2D to do this.

% Prealloc for speed
accelSurfAvg11 = zeros(slice.xMaxInd,time.numFrames);
accelSurfAvg22 = zeros(slice.xMaxInd,time.numFrames);
stress12Avg = zeros(slice.xMaxInd,time.numFrames);
stress22Avg = zeros(slice.xMaxInd,time.numFrames);

% Only use these vars for debugging to view masked accel fields
%accel11Masked = zeros(size(accel.x)));
%accel22Masked = zeros(size(accel.x)));

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
        accel11Temp = accel.mat11(:,:,ff);
        accel22Temp = accel.mat22(:,:,ff);
        
        % Set the initial frame to nan to ignore values with nanmean
        accel11Frame = nan(length(pos.y),length(pos.x));
        accel22Frame = nan(length(pos.y),length(pos.x));
        accel11Frame(sliceMask) = accel11Temp(sliceMask);
        accel22Frame(sliceMask) = accel22Temp(sliceMask);
        
        % Used for debugging to view masked acceleration
        %accel11Masked(:,:,ff) = accel11Frame; 
        %accel22Masked(:,:,ff) = accel22Frame;
        
        % Calculate the average of the masked field
        accelSurfAvg11(xx,ff) = nanmean(accel11Frame(:));
        accelSurfAvg22(xx,ff) = nanmean(accel22Frame(:));

        % Calculate the area of the 'sliced' surface
        S = specimen.height*(pos.x(xx)+slice.lengthX/2);

        % Calculate the average stress using the angled stress gauge
        stress12Avg(xx,ff) = -material.rho*S/slice.length * accelSurfAvg11(xx,ff);
        stress22Avg(xx,ff) = -material.rho*S/slice.length * accelSurfAvg22(xx,ff);
    end
end

end

