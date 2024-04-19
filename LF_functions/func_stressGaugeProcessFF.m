function [stressXAvg,accelSurfAvg] = func_stressGaugeProcessFF(material,time,pos,accel)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 19/8/2018
% Date Edited: 19/8/2018
%
% Processes width averaged accel data to obtain stress 
% sigma(x,t),lineAvg = material.rho*x*a(x,t),areaAvg

% Pre-alloc for speed
stressXAvg = zeros(size(accel.x,2),size(accel.x,3));
accelSurfAvg = zeros(size(accel.x,2),size(accel.x,3));

% Loop over each frame and each transverse slice
for tt = 1:time.numFrames
    for xx = 1:length(pos.x)
        % average accel averaged from free end (x = 0) to x
        tempAccel = accel.x(:,1:xx,tt);
        accelSurfAvg(xx,tt) = mean(tempAccel(:));
        % reconstructed stress on a line at x
        stressXAvg(xx,tt) = material.rho*pos.x(xx)*accelSurfAvg(xx,tt); 
    end
end

end

