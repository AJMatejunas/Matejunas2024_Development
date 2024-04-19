function [disp,vel,accel,strain,strainRate] = func_smoothCalcAccelStrain_Ebal(...
    pos,time,material,grid,disp,smoothingOpts,extrapOpts,fieldCalcOpts)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 27/3/2017
%
% Processes output displacement fields from the grid method to obtain width
% averaged variables for processing with the stress gauge approach. Output 
% variables are obtained via numerical differentiation using a central 
% difference method. Variables output include velocity, acceleration, strain 
% and strain rate. Width averaging then smoothing is performed before taking 
% derivatives. Following differentiation variables are extrapolated spatially 
% to regain the missing pitch from the grid method.

% Selects method for numerical differentiation
diffOpt = fieldCalcOpts.diffOpt;

% Initialise the variables we are not returning
if ~fieldCalcOpts.calcVel
    vel = -1;
end
if ~fieldCalcOpts.calcAccel
    accel = -1;
end
if ~fieldCalcOpts.calcStrainRate
    strainRate = -1;
end

%--------------------------------------------------------------------------
% Check for NaNs in the displacement fields
nanFlag = false;
if squeeze(sum(sum(sum(isnan(disp.x))))) > 0
    fprintf('WARNING: X displacement field contains NaNs.\n')
    nanFlag = true;
end
if squeeze(sum(sum(sum(isnan(disp.y))))) > 0
    fprintf('WARNING: Y displacement field contains NaNs.\n')
    nanFlag = true;
end

%--------------------------------------------------------------------------
% Fix the NaNs in the displacement fields if required
if extrapOpts.fixNaNs == true && nanFlag == true
    fprintf('\tReplacing NaNs with local average over a distance of %i px.\n',extrapOpts.fixNaNKernal)
    fprintf('\tReplacing NaNs in X displacement.\n')
    disp.x = func_fixNaNsInDisp(extrapOpts,disp.x);
    fprintf('\tReplacing NaNs in Y displacement.\n')
    disp.y = func_fixNaNsInDisp(extrapOpts,disp.y);
end

%--------------------------------------------------------------------------
% FULL-FIELD:  Temporal Smoothing of Full-Field Displacements to Create 
% Accel/Strain Rate Maps
if smoothingOpts.FFTempSmooth
    fprintf('\tTemporally smoothing raw full-field displacements.\n')
    frameSize = smoothingOpts.FFTemporalKernal(1);
    polyOrder = smoothingOpts.FFTemporalKernal(2);
    fprintf('\tTemporal full-field filter is %s of order %i and frame width %i.\n',...
        smoothingOpts.FFTemporalFilt,polyOrder,frameSize)
    
    if smoothingOpts.FFTemporalPad
        % Pad the displacement matrices
        padFrames = ceil(frameSize/2);
        dispPadX = padarray(disp.x,[0,0,padFrames],'replicate','pre'); 
        dispPadY = padarray(disp.y,[0,0,padFrames],'replicate','pre');
        
        % Smooth the padded array with the golay filter
        disp.tSmooth.x = sgolayfilt(dispPadX,polyOrder,frameSize,[],3);
        disp.tSmooth.y = sgolayfilt(dispPadY,polyOrder,frameSize,[],3);
        clear dispPadX dispPadY
        
        % Crop the vector back to the original size to remove padding
        disp.tSmooth.x = disp.tSmooth.x(:,:,(padFrames+1):end);
        disp.tSmooth.y = disp.tSmooth.y(:,:,(padFrames+1):end);
    else
        disp.tSmooth.x = sgolayfilt(disp.x,polyOrder,frameSize,[],3);
        disp.tSmooth.y = sgolayfilt(disp.y,polyOrder,frameSize,[],3);
    end
end

%--------------------------------------------------------------------------
% FULL-FIELD: Calculate full-field acceleration maps with no spatial
% smoothing to remove edge effects
fprintf('\tTemporally differentiating to obtain full-field velocity and acceleration maps.\n')
if fieldCalcOpts.calcAccel && fieldCalcOpts.calcVel
    if smoothingOpts.FFTempSmooth == true 
        [vel,accel] = func_calcFFAccelFromDisp(time,disp.tSmooth,diffOpt);
    else
        [vel,accel] = func_calcFFAccelFromDisp(time,disp,diffOpt);
    end
elseif fieldCalcOpts.calcVel
    if smoothingOpts.FFTempSmooth == true 
        [vel,~] = func_calcFFAccelFromDisp(time,disp.tSmooth,diffOpt);
    else
        [vel,~] = func_calcFFAccelFromDisp(time,disp,diffOpt);
    end
elseif fieldCalcOpts.calcAccel
    if smoothingOpts.FFTempSmooth == true 
        [~,accel] = func_calcFFAccelFromDisp(time,disp.tSmooth,diffOpt);
    else
        [~,accel] = func_calcFFAccelFromDisp(time,disp,diffOpt);
    end
end 

if fieldCalcOpts.reduceRAM
    if isfield(disp,'tSmooth')
        % Remove all smoothed displacement fields
        disp = rmfield(disp,'tSmooth');
    end
end

%--------------------------------------------------------------------------
% FULL-FIELD:  Spatial Smoothing Over the Full-Field
if smoothingOpts.spatSmooth == true
    fprintf('\tSpatially smoothing raw full-field displacements.\n')
    if strcmp(smoothingOpts.spatialFilt,'median')
        medFiltKernel = [smoothingOpts.spatialKernal(1),smoothingOpts.spatialKernal(1)];
        fprintf('\tSpatial full-field filter is %s with a [%i,%i] kernal.\n',...
            smoothingOpts.spatialFilt,medFiltKernel(1),medFiltKernel(2))
        for f = 1:time.numFrames 
           disp.sSmooth.x(:,:,f) = medfilt2(disp.x(:,:,f),medFiltKernel,smoothingOpts.spatialEdgeMode); 
           disp.sSmooth.y(:,:,f) = medfilt2(disp.y(:,:,f),medFiltKernel,smoothingOpts.spatialEdgeMode);
        end
    else
        kernelLength = smoothingOpts.spatialKernal(1);
        kernelRadius = ceil((kernelLength-1)/2);
        fprintf('\tSpatial full-field filter is %s with radius of %ipx and length of %ipx.\n',...
            smoothingOpts.spatialFilt,kernelRadius,kernelLength)
        gaussFilter = fspecial('gaussian',[kernelLength,kernelLength],kernelRadius);
        for f = 1:time.numFrames 
            disp.sSmooth.x(:,:,f) = imfilter(disp.x(:,:,f),gaussFilter,smoothingOpts.spatialEdgeMode);
            disp.sSmooth.y(:,:,f) = imfilter(disp.y(:,:,f),gaussFilter,smoothingOpts.spatialEdgeMode);
        end   
    end
end

%--------------------------------------------------------------------------
% FULL-FIELD: Calculate Strains from Full-Field Smoothed Displacement
fprintf('\tSpatially differentiating to obtain strains from smoothed displacement fields.\n')
if smoothingOpts.spatSmooth == true
    strain = func_calcStrainFromDisp(disp.sSmooth,grid.mPerPx,grid.mPerPx);
else
    strain = func_calcStrainFromDisp(disp,grid.mPerPx,grid.mPerPx);
end

%--------------------------------------------------------------------------
% FULL-FIELD: Calculate the Full-Field Strain Rate
if fieldCalcOpts.calcStrainRate
    fprintf('\tTemporally smoothing strains and calculating strain rate.\n')
    [strainRate,~] = func_calculateStrainRate(strain,time,smoothingOpts);
end

if fieldCalcOpts.reduceRAM
    if isfield(disp,'sSmooth')
        % Remove all smoothed displacement fields
        disp = rmfield(disp,'sSmooth');
    end
end

%--------------------------------------------------------------------------
% WIDTH AVG: Average displacement and strain over the width
if fieldCalcOpts.calcWidthAvgs
    fprintf('\tAveraging smoothed strain and raw displacement over the width.\n')
    disp.xAvg = func_avgFFVarOverWidth(disp.x);             % Raw (Noisy) Data
    strain.xAvg = func_avgFFVarOverWidth(strain.x);         % Smoothed
    strain.yAvg = func_avgFFVarOverWidth(strain.y);         % Smoothed
    strainRate.xAvg = func_avgFFVarOverWidth(strainRate.x); % Smoothed 
    strain.xnyAvg = strain.xAvg + material.nuxy*strain.yAvg;
end

%--------------------------------------------------------------------------
% WIDTH AVG: get the peak axial avg strain rates
if fieldCalcOpts.calcStrainRate && fieldCalcOpts.calcWidthAvgs
    if smoothingOpts.spatSmooth == true
        edgePx = grid.pxPerPeriod+round(smoothingOpts.spatialKernal(1)/2)+1;
    else
        edgePx = grid.pxPerPeriod+1;
    end
    xRange = edgePx:size(strainRate.xAvg,1) - edgePx;
    % Ignore the edges where the data is corrupted by smoothing effects
    strainRate.xAvgMin = min(min(strainRate.xAvg(xRange,:)));
    strainRate.xAvgMax = max(max(strainRate.xAvg(xRange,:)));
end

%--------------------------------------------------------------------------
% WIDTH AVG: Temporal Smoothing of Width Averaged Data
if fieldCalcOpts.calcWidthAvgs
    fprintf('\tTemporally smoothing width averaged displacement.\n')
    if smoothingOpts.WATempSmooth == true
        frameSize = smoothingOpts.WATemporalKernal(1);
        polyOrder = smoothingOpts.WATemporalKernal(2);

        if smoothingOpts.WATemporalPad
            padFrames = ceil(frameSize/2);
            dispPadX = padarray(disp.xAvg,[0,padFrames],'replicate','pre');
            disp.tSmooth.xAvg = sgolayfilt(dispPadX,polyOrder,frameSize,[],2);
            clear dispPadX
            %disp.tSmooth.xAvg = disp.tSmooth.xAvg(:,(padFrames+1):(end-padFrames));
            disp.tSmooth.xAvg = disp.tSmooth.xAvg(:,(padFrames+1):end);
        else
            disp.tSmooth.xAvg = sgolayfilt(disp.xAvg,polyOrder,frameSize,[],2);
        end
    end
end

%--------------------------------------------------------------------------
% WIDTH AVG: Calculate the Width Averaged Acceleration
if fieldCalcOpts.calcWidthAvgs
    fprintf('\tTemporally differentiating width averaged data to obtain acceleration.\n')
    if smoothingOpts.WATempSmooth == true
        [vel.xAvg,accel.xAvg] = func_calcAccelFromDisp(time,disp.tSmooth,diffOpt);
    else
        [vel.xAvg,accel.xAvg] = func_calcAccelFromDisp(time,disp,diffOpt);
    end
end

%--------------------------------------------------------------------------
% WIDTH AVG: Extrapolate the width averaged data to regain the missing pitch
if fieldCalcOpts.calcWidthAvgs
    if extrapOpts.extrapWidthAvgData == true
        fprintf('\tSpatially extrapolating width averaged data to regain the missing pitch.\n')
        disp.xAvg = func_spatialExtrapWidthAvgVar(disp.xAvg,pos,extrapOpts.dispPx,extrapOpts.dispMethod);
        accel.xAvg = func_spatialExtrapWidthAvgVar(accel.xAvg,pos,extrapOpts.dispPx,extrapOpts.dispMethod);
        strain.xAvg = func_spatialExtrapWidthAvgVar(strain.xAvg,pos,extrapOpts.strainPx,extrapOpts.strainMethod);
    end
end

%--------------------------------------------------------------------------
% FULL-FIELD: Spatially extrapolate full field data for use with the VFM
if extrapOpts.extrapFullFieldData == true
    fprintf('\tSpatially extrapolating full-field accel and strain for use with the VFM.\n')
    if fieldCalcOpts.calcAccel
        % Account for the size of the spatial smoothing kernal
        fprintf('\t\tExtrapolating X acceleration.\n')
        accel.x = func_spatialExtrapFullFieldVar(pos,accel.x,2,extrapOpts.FFDispPx,extrapOpts.FFDispMethod);
        fprintf('\t\tExtrapolating Y acceleration.\n')
        accel.y = func_spatialExtrapFullFieldVar(pos,accel.y,1,extrapOpts.FFDispPx,extrapOpts.FFDispMethod);
    end
    fprintf('\t\tExtrapolating X strain.\n')
    strain.x = func_spatialExtrapFullFieldVar(pos,strain.x,2,extrapOpts.FFStrainPx,extrapOpts.FFStrainMethod);
    fprintf('\t\tExtrapolating Y strain.\n')
    strain.y = func_spatialExtrapFullFieldVar(pos,strain.y,1,extrapOpts.FFStrainPx,extrapOpts.FFStrainMethod);
    fprintf('\t\tExtrapolating XY strain.\n')
    strain.s = func_spatialExtrapFullFieldVar(pos,strain.s,1,extrapOpts.FFStrainPx,extrapOpts.FFStrainMethod);
    strain.s = func_spatialExtrapFullFieldVar(pos,strain.s,2,extrapOpts.FFStrainPx,extrapOpts.FFStrainMethod);
    
    % Check if this field exists for backwards compatibility
    if isfield(extrapOpts,'extrapForEBal')
        fprintf('\tSpatially extrapolating full-field disp and vel for use with the E balance.\n')
        if extrapOpts.extrapForEBal
            fprintf('\t\tExtrapolating X displacement.\n')
            disp.x = func_spatialExtrapFullFieldVar(pos,disp.x,2,extrapOpts.FFDispPx,extrapOpts.FFDispMethod);
            fprintf('\t\tExtrapolating Y displacement.\n')
            disp.y = func_spatialExtrapFullFieldVar(pos,disp.y,1,extrapOpts.FFDispPx,extrapOpts.FFDispMethod);
            fprintf('\t\tExtrapolating X velocity.\n')
            vel.x = func_spatialExtrapFullFieldVar(pos,vel.x,2,extrapOpts.FFDispPx,extrapOpts.FFDispMethod);
            vel.x = func_spatialExtrapFullFieldVar(pos,vel.x,1,extrapOpts.FFDispPx,extrapOpts.FFDispMethod);
            fprintf('\t\tExtrapolating Y velocity.\n')
            vel.y = func_spatialExtrapFullFieldVar(pos,vel.y,1,extrapOpts.FFDispPx,extrapOpts.FFDispMethod);
            vel.y = func_spatialExtrapFullFieldVar(pos,vel.y,2,extrapOpts.FFDispPx,extrapOpts.FFDispMethod);
            
            fprintf('\t\tExtrapolating X strain in Y dim.\n')
            strain.x = func_spatialExtrapFullFieldVar(pos,strain.x,1,extrapOpts.FFStrainPx,extrapOpts.FFStrainMethod);
            fprintf('\t\tExtrapolating Y strain in X dim.\n')
            strain.y = func_spatialExtrapFullFieldVar(pos,strain.y,2,extrapOpts.FFStrainPx,extrapOpts.FFStrainMethod);
        end
    end
    fprintf('\tFull-field extrapolation complete.\n')
end

%--------------------------------------------------------------------------
% Remove Fields from Structures to Save RAM
if fieldCalcOpts.reduceRAM
    fprintf('\tRemoving some fields from data structs to reduce RAM usage.\n')
    % Remove all smoothed displacement fields
    if isfield(disp,'tSmooth')
        disp = rmfield(disp,'tSmooth');
    end
    if isfield(disp,'sSmooth')
        disp = rmfield(disp,'sSmooth');
    end
    if isfield(disp,'tsSmooth')
        disp = rmfield(disp,'tsSmooth');
    end
end

fprintf('\tAcceleration and Strain Calculation Complete.\n')

end