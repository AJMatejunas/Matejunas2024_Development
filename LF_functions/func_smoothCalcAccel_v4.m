function [accel,vel,disp] = func_smoothCalcAccel_v4(pos,time, ...
    disp,...
    smoothOptsA,extrapOptsA,diffOpts,printToCons)
% TODO
%
% Author: Lloyd Fletcher
% Image-Based Mechanics Group (IBMG), University of Southampton
% Date Created: 31/7/2020 - Major update (v4)! Removed old width averaging
% before smoothing option, added better temporal padding options and post
% smoothing spatial edge extrapolation
% Date Edited: 31/7/2020.

% Allow printing to the console by default
if nargin < 7
    printToCons = true;
end
debug = false;

%--------------------------------------------------------------------------
% 1) Temporally pre and post pad the displacement data if required
if extrapOptsA.tempPadOn
    if printToCons
    fprintf('\tTemporal padding enabled.\n')
    fprintf('\t\tTemporal pad method is: %s.\n',extrapOptsA.tempPadMethod)
    fprintf('\t\tPre-padding with %i frames.\n',extrapOptsA.tempPadFrames(1))
    fprintf('\t\tPost-padding with %i frames.\n',extrapOptsA.tempPadFrames(2))
    end
    [dispTempPad,timeExtrap] = func_temporalPadDisp(extrapOptsA,time,disp,debug);
else
    if printToCons
    fprintf('\tTemporal padding disabled.\n')
    end
    dispTempPad.x = disp.x;
    dispTempPad.y = disp.y;
    timeExtrap = time;
    timeExtrap.keepRange = 1:length(time.vec);
end
disp.rT = timeExtrap.keepRange;

%--------------------------------------------------------------------------
% 2) Temporally smooth the displacements
if smoothOptsA.temporalSmooth
    if printToCons
        fprintf('\tTemporal smoothing of displacements enabled.\n')
        fprintf('\t\tTemporal filter type is %s.\n',smoothOptsA.temporalType)
        fprintf('\t\tTemporal filter algorithm is %s.\n',smoothOptsA.temporalAlgorithm)
        fprintf('\t\tX field temporal kernel is %i with order %i.\n',...
            smoothOptsA.temporalKernelSize(1),smoothOptsA.temporalKernelOrder(1))
        fprintf('\t\tY field temporal kernel is %i with order %i.\n',...
            smoothOptsA.temporalKernelSize(2),smoothOptsA.temporalKernelOrder(2))
    end
    disp.tSmooth.x = func_reshapeAndGolayFilt3D(dispTempPad.x,...
        smoothOptsA.temporalKernelOrder(1),smoothOptsA.temporalKernelSize(1));
    disp.tSmooth.y = func_reshapeAndGolayFilt3D(dispTempPad.y,...
        smoothOptsA.temporalKernelOrder(2),smoothOptsA.temporalKernelSize(2));
else
    if printToCons
    fprintf('\tTemporal smoothing of displacements disabled.\n')
    end
    disp.tSmooth.x = dispTempPad.x;
    disp.tSmooth.y = dispTempPad.y;
end

%--------------------------------------------------------------------------
% 3) Temporally differentiate to obtain velocity and acceleration
[vel,accel] = func_calcFFAccelFromDisp(time,disp.tSmooth,diffOpts.method);

%--------------------------------------------------------------------------
% 4) Crop back fields in time to original number of frames if they were
% padded
if extrapOptsA.tempPadOn
    vel.x = vel.x(:,:,timeExtrap.keepRange);
    vel.y = vel.y(:,:,timeExtrap.keepRange);
    accel.x = accel.x(:,:,timeExtrap.keepRange);
    accel.y = accel.y(:,:,timeExtrap.keepRange);
end

%--------------------------------------------------------------------------
% 4.1) Optional spatial edge extrapolation directly from acceleration fields
if extrapOptsA.postSpatExtrapOn
    if printToCons
        fprintf('\tOptional accel spatial crop and extrapolation enabled.\n')
        fprintf('\tOptional accel spatial crop and extrapolation processing...\n')
    end
    
    %----------------------------------------------------------------------
    % CROPPING
    [accelCrop,posCropX,posCropY] = func_cropAccelFields(...
        extrapOptsA,pos,accel);
    
    %----------------------------------------------------------------------
    % EXTRAPOLATION
    
    accel = func_extrapolateAccelFields(...
        extrapOptsA,pos,accelCrop,posCropX,posCropY,debug);
    clear accelCrop
end

%--------------------------------------------------------------------------
% 5) Average the full-field acceleration over the width.
if printToCons
    fprintf('\tAveraging full-field acceleration over the sample height.\n')
end
accel.xAvg = func_avgFFVarOverWidth(accel.x);  
accel.yAvg = func_avgFFVarOverWidth(accel.y);  

if printToCons
    fprintf('\tACCELERATION calculation complete.\n')
end

end