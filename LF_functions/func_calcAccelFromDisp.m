function [velXAvg,accelXAvg] = func_calcAccelFromDisp(time,disp,option)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 9/3/2017
% Calculates acceleration and velocity from the displacement data
    if nargin < 3
        option = 'gradient';
    end

    % Differentiate in time to get the velocity and acceleration
    if strcmp(option,'cDiff')
        for x = 1:size(disp.xAvg,1)
            velXAvg(x,:) = gradient(disp.xAvg(x,:),time.step);
            accelXAvg(x,:) = func_centralDiff2(disp.xAvg(x,:),time.step);
        end
    else
        for x = 1:size(disp.xAvg,1)
            velXAvg(x,:) = gradient(disp.xAvg(x,:),time.step);
            accelXAvg(x,:) = gradient(velXAvg(x,:),time.step);
        end
    end
end

