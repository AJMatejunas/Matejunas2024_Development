function [stressXAvg,accelSurfAvg] = func_stressGaugeProcess(material,time,pos,accelWAvg,shiftData)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 9/2/2017
% Date Edited: 10/6/2019
%
% Processes width averaged accel data to obtain stress 
% sigma(x,t),lineAvg = material.rho*x*a(x,t),areaAvg

% If type is unspecified assume this is a grid not nodal FE data
if nargin < 5
    shiftData = true;
end

% Pre-alloc for speed
stressXAvg = zeros(length(pos.x),time.numFrames);
accelSurfAvg = zeros(length(pos.x),time.numFrames);

% Loop over each frame and each transverse slice
for tt = 1:time.numFrames
    for xx = 1:length(pos.x)
        % average accel averaged from free end (x = 0) to x
        accelSurfAvg(xx,tt) = mean(accelWAvg(1:xx,tt));
        
        % reconstructed stress on a line at x
        if shiftData
            stressXAvg(xx,tt) = material.rho*(pos.x(xx)+pos.xStep/2)*accelSurfAvg(xx,tt); 
        else
            stressXAvg(xx,tt) = material.rho*pos.x(xx)*accelSurfAvg(xx,tt); 
        end
    end
end

end

