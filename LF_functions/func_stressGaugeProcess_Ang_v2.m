function [stress22Avg,stress12Avg] = func_stressGaugeProcess_Ang_v2(slice,material,time,pos,accel)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 12/7/2018
%
% Calculates stress averages along an angled slice at the material co-ord
% system angle. Angle is specified in the material.rotAngle parameter.
% NOTE: acceleration field needs to be in material co-ords. The structure
% should have the fields .mat11 and .mat22 in the material co-ords. Use the
% function func_rotateVector2D to do this. The slice properties slice.area
% and slice.length must also be created to use this function.

for ff = 1:time.numFrames
    for ss = 1:length(pos.x)
        % Create a logical mask for the angled slice - get all points within the
        % bounds
        for ii = 1:length(pos.y)
           sliceMask(ii,:) = pos.x < slice.linesXCoords{ss};  
        end

        % Calculate the acceleration field average from the masked field
        accel11Temp = accel.mat11(:,:,ff);
        accel22Temp = accel.mat22(:,:,ff);
        
        % Initialise the fields to use with nanmean to ignore points
        accel11Frame = nan(length(pos.y),length(pos.x));
        accel22Frame = nan(length(pos.y),length(pos.x)); 
        accel11Frame(sliceMask) = accel11Temp(sliceMask);
        accel22Frame(sliceMask) = accel22Temp(sliceMask);
        accel11Masked(:,:,ff) = accel11Frame;
        accel22Masked(:,:,ff) = accel22Frame;
        
        % Calculate the mean over the selected surfac
        accelSurfAvg11(ss,ff) = nanmean(accel11Frame(:));
        accelSurfAvg22(ss,ff) = nanmean(accel22Frame(:));

        % Calculate the average stress
        stress12Avg(ss,ff) = -material.rho*slice.area(ss)/slice.length(ss)...
            *accelSurfAvg11(ss,ff);
        stress22Avg(ss,ff) = -material.rho*slice.area(ss)/slice.length(ss)...
            *accelSurfAvg22(ss,ff);

        % Reset this loop iteration
        clear sliceMask;
    end
end

end

