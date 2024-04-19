function [vel,accel] = func_calcFFAccelFromDisp(time,disp,option)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 16/3/2017
    
    if nargin < 3
        option = 'gradient';
    end
    
    [~,~,vel.x] = gradient(disp.x,1,1,time.step);
    [~,~,vel.y] = gradient(disp.y,1,1,time.step);
    if strcmp(option,'cDiff')
        for y = 1:size(disp.x,1)
            for x = 1:size(disp.x,2)
                accel.x(y,x,:) = func_centralDiff2(disp.x(y,x,:),time.step);
                accel.y(y,x,:) = func_centralDiff2(disp.y(y,x,:),time.step);
            end
        end   
    else    
        [~,~,accel.x] = gradient(vel.x,1,1,time.step);
        [~,~,accel.y] = gradient(vel.y,1,1,time.step);
    end
end

