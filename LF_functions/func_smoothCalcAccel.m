function [disp,vel,accel] = func_smoothCalcAccel(pos,time,grid,disp,...
    smoothingOpts,extrapOpts,diffOpts,printToCons)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 27 Mar. 2017
% Modified by: Jared Van Blitterswyk 
% Date: 15 Oct. 2017
%
% Processes output displacement fields from the grid method to obtain width
% averaged variables for processing with the stress gauge approach. Output 
% variables are obtained via numerical differentiation using a central 
% difference method. Variables output include velocity, acceleration.

% Allow printing to the console by default
if nargin < 8
    printToCons = true;
end

%% FULL-FIELD: Temporal Smoothing of Full-Field Displacements for Accel 

% Crop one pitch from the edges to remove intial rubbish data
dispCrop.x = func_cropFields(disp.x,2,extrapOpts.FFDispPx);
dispCrop.y = func_cropFields(disp.y,1,extrapOpts.FFDispPx);

%--------------------------------------------------------------------------
% Fix the NaNs in the displacement fields if required
if sum(isnan(dispCrop.x(:))) > 0
    if printToCons
        fprintf('\tWARNING: NaNs found in x displacement field.\n')
        fprintf('\tReplacing NaNs with local average over a distance of %i px.\n',extrapOpts.fixNaNKernal)
    end
    dispCrop.x = func_fixNaNsInDispX(extrapOpts,dispCrop.x);
    if printToCons
        if sum(isnan(dispCrop.x(:))) > 0
            fprintf('\tWARNING: NaNs still in x displacement field, widen kernel.\n')
        end
    end
end

if sum(isnan(dispCrop.y(:))) > 0
    if printToCons
        fprintf('\tWARNING: NaNs found in y displacement field.\n')
        fprintf('\tReplacing NaNs with local average over a distance of %i px.\n',extrapOpts.fixNaNKernal)
    end
    dispCrop.y = func_fixNaNsInDispY(extrapOpts,dispCrop.y);
    if printToCons
        if sum(isnan(dispCrop.y(:))) > 0
            fprintf('\tWARNING: NaNs still in y displacement field, widen kernel.\n')
        end
    end
end

% Pad out the array to reduce the edge effects
padRadius = 3*max([smoothingOpts.spatialKernal(1),smoothingOpts.spatialKernal(2),grid.pxPerPeriod]);
dispTemp = func_padDispArrayExt(pos,dispCrop,grid,padRadius,...
    extrapOpts.padEdgeMethod,extrapOpts.dispPx,extrapOpts.padFitWindow,...
    extrapOpts); 

% Crop back to the original FOV
dispTemp.x = dispTemp.x(:,...
    padRadius-extrapOpts.FFDispPx+1:end-padRadius+extrapOpts.FFDispPx,:);
dispTemp.y = dispTemp.y(padRadius-extrapOpts.FFDispPx+1:...
    end-padRadius+extrapOpts.FFDispPx,:,:);

%--------------------------------------------------------------------------
% TEMPORAL SMOOTHING of full-field displacements
% Option to pad the data to avoid edge effects from the filter
if smoothingOpts.FFTempSmooth
    frameSize = smoothingOpts.FFTemporalKernal(1);
    polyOrder = smoothingOpts.FFTemporalKernal(2);
    if printToCons
        fprintf('\tTemporally smoothing full-field displacements.\n')
        fprintf('\tTemporal full-field filter is %s of order %i and frame width %i.\n',...
            smoothingOpts.FFTemporalFilt,polyOrder,frameSize)
    end

    if smoothingOpts.FFTemporalPad
        % Pad the displacement matrices in time for the golay filter
        padFrames = ceil(frameSize/2)+smoothingOpts.FFTemporalPadFrames;
        dispTemporalPadX = padarray(dispTemp.x,[0,0,padFrames],smoothingOpts.FFTemporalPadMethod,'pre'); 
        dispTemporalPadY = padarray(dispTemp.y,[0,0,padFrames],smoothingOpts.FFTemporalPadMethod,'pre');
        
        disp.tSmooth.x = func_reshapeAndGolayFilt3D(dispTemporalPadX,smoothingOpts);
        disp.tSmooth.y = func_reshapeAndGolayFilt3D(dispTemporalPadY,smoothingOpts);
        
        % Crop the back to the original size to remove padding
        disp.tSmooth.x = disp.tSmooth.x(:,:,(padFrames+1):end);
        disp.tSmooth.y = disp.tSmooth.y(:,:,(padFrames+1):end);
    else
        disp.tSmooth.x = func_reshapeAndGolayFilt3D(dispTemp.x,smoothingOpts);
        disp.tSmooth.y = func_reshapeAndGolayFilt3D(dispTemp.y,smoothingOpts);
    end
else
    disp.tSmooth.x = dispTemp.x;
    disp.tSmooth.y = dispTemp.y;
end

%--------------------------------------------------------------------------
% TEMPORAL DIFFERENTIATION of full-field displacements
% Option to pad the data to avoid edge effects from central diff
if printToCons
    fprintf('\tTemporally differentiating to obtain full-field acceleration maps.\n')
end
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

%% WIDTH AVERAGED: Smoothing and differentiation of width averaged accel
% Set variables for backwards compatibility
if ~isfield(smoothingOpts,'WATemporalAvgFirst')
    % Default to averaging over the width first
    smoothingOpts.WATemporalAvgFirst = false; 
end

% Average displacement over the width
if printToCons
    fprintf('\tAveraging displacement fields over the sample height.\n')
end
disp.xAvg = func_avgFFVarOverWidth(disp.tSmooth.x);  

if smoothingOpts.WATemporalAvgFirst
    %--------------------------------------------------------------------------
    % Temporal Smoothing of Width Averaged Data
    if printToCons
        fprintf('\tTemporally smoothing width averaged displacement.\n')
    end
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

    % Calculate the Width Averaged Acceleration
    if printToCons
        fprintf('\tTemporally differentiating width averaged data to obtain acceleration.\n')
    end
    if smoothingOpts.WATempSmooth == true
        [~,accel.xAvg] = func_calcAccelFromDisp(time,disp.tSmooth,diffOpts);
    else
        [~,accel.xAvg] = func_calcAccelFromDisp(time,disp,diffOpts);
    end
else
    % Average displacement over the width
    if printToCons
        fprintf('\tAveraging full-field acceleration over the sample height.\n')
    end
    accel.xAvg = func_avgFFVarOverWidth(accel.x);  
end

end