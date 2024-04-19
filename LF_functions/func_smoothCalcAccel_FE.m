function [disp,vel,accel] = func_smoothCalcAccel_FE(pos,time,grid,disp,smoothingOpts,extrapOpts,diffOpts)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 27/3/2017
% Modified by: Jared Van Blitterswyk
% Date: 15 Oct. 2017
%
% Processes output displacement fields from the grid method to obtain width
% averaged variables for processing with the stress gauge approach. Output 
% variables are obtained via numerical differentiation using a central 
% difference method. Variables output include velocity, acceleration. 

%% FULL-FIELD: Temporal Smoothing of Full-Field Displacements for Accel 

%--------------------------------------------------------------------------
% TEMPORAL SMOOTHING of full-field displacements
% Option to pad the data to avoid edge effects from the filter
if smoothingOpts.FFTempSmooth
    fprintf('\tTemporally smoothing full-field displacements.\n')
    frameSize = smoothingOpts.FFTemporalKernal(1);
    polyOrder = smoothingOpts.FFTemporalKernal(2);
    fprintf('\tTemporal full-field filter is %s of order %i and frame width %i.\n',...
        smoothingOpts.FFTemporalFilt,polyOrder,frameSize)

    if smoothingOpts.FFTemporalPad
        % Pad the displacement matrices in time for the golay filter
        padFrames = ceil(frameSize/2)+smoothingOpts.FFTemporalPadFrames;
        dispTemporalPadX = padarray(disp.x,[0,0,padFrames],smoothingOpts.FFTemporalPadMethod,'pre'); 
        dispTemporalPadY = padarray(disp.y,[0,0,padFrames],smoothingOpts.FFTemporalPadMethod,'pre');
        
        disp.tSmooth.x = func_reshapeAndGolayFilt3D(dispTemporalPadX,smoothingOpts);
        disp.tSmooth.y = func_reshapeAndGolayFilt3D(dispTemporalPadY,smoothingOpts);
        
        % Crop the back to the original size to remove padding
        disp.tSmooth.x = disp.tSmooth.x(:,:,(padFrames+1):end);
        disp.tSmooth.y = disp.tSmooth.y(:,:,(padFrames+1):end);
    else
        disp.tSmooth.x = func_reshapeAndGolayFilt3D(disp.x,smoothingOpts);
        disp.tSmooth.y = func_reshapeAndGolayFilt3D(disp.y,smoothingOpts);
    end
else
    disp.tSmooth.x = disp.x;
    disp.tSmooth.y = disp.y;
end

%--------------------------------------------------------------------------
% TEMPORAL DIFFERENTIATION of full-field displacements
% Option to pad the data to avoid edge effects from central diff
fprintf('\tTemporally differentiating to obtain full-field acceleration maps.\n')
if diffOpts.temporalPad
    padFrames = diffOpts.temporalPadFrames;
    
    % Pad frames at the start of the array
    disp.tSmooth.x = padarray(disp.tSmooth.x,[0,0,padFrames],diffOpts.temporalPadMethod,'pre'); 
    disp.tSmooth.y = padarray(disp.tSmooth.y,[0,0,padFrames],diffOpts.temporalPadMethod,'pre');
    
    % Calculate velocity and acceleration
    [vel,accel] = func_calcFFAccelFromDisp(time,disp.tSmooth,diffOpts.method);
    
    % Crop back to original size
    vel.x = vel.x(:,:,padFrames+1:end);
    vel.y = vel.y(:,:,padFrames+1:end);
    accel.x = accel.x(:,:,padFrames+1:end);
    accel.y = accel.y(:,:,padFrames+1:end);
    disp.tSmooth.x = disp.tSmooth.x(:,:,1+padFrames:end);
    disp.tSmooth.y = disp.tSmooth.y(:,:,1+padFrames:end);
else
    [vel,accel] = func_calcFFAccelFromDisp(time,disp.tSmooth,diffOpts.method);
end

%% WIDTH AVG: Average displacement and strain over the width
fprintf('\tAveraging displacement fields over the sample height.\n')
disp.xAvg = func_avgFFVarOverWidth(disp.tSmooth.x);         

%--------------------------------------------------------------------------
% Temporal Smoothing of Width Averaged Data
fprintf('\tTemporally smoothing width averaged displacement.\n')
if smoothingOpts.WATempSmooth == true
    frameSize = smoothingOpts.WATemporalKernal(1);
    polyOrder = smoothingOpts.WATemporalKernal(2);
    
    % Option to pad the array in time to mitigate filter edge effects
    if smoothingOpts.WATemporalPad
        padFrames = ceil(frameSize/2)+smoothingOpts.WATemporalPadFrames;
        dispTemporalPadX = padarray(disp.xAvg,[0,padFrames],smoothingOpts.WATemporalPadMethod,'pre');
        
        % Temporal smoothing
        disp.tSmooth.xAvg = sgolayfilt(dispTemporalPadX,polyOrder,frameSize,[],2);

        % Crop back to the original size
        disp.tSmooth.xAvg = disp.tSmooth.xAvg(:,(padFrames+1):end);
    else
        disp.tSmooth.xAvg = sgolayfilt(disp.xAvg,polyOrder,frameSize,[],2);
    end
end

%% WIDTH AVG: Calculate the Width Averaged Acceleration
fprintf('\tTemporally differentiating width averaged data to obtain acceleration.\n')
if smoothingOpts.WATempSmooth == true
    [~,accel.xAvg] = func_calcAccelFromDisp(time,disp.tSmooth,diffOpts);
else
    [~,accel.xAvg] = func_calcAccelFromDisp(time,disp,diffOpts);
end

end