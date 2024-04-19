function [disp,vel,accel,strain,strainRate] = func_smoothCalcAccelStrain_WAvg(...
    pos,time,material,grid,disp,smoothingOpts,extrapOpts,saveRAM)
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

% if this argument is not specified clear fields in the structures to save
% space
if nargin < 7
    saveRAM = true;
end

%--------------------------------------------------------------------------
% WIDTH AVG: Pad Data to Remove NaNs over the width
if extrapOpts.padWidthNans == true
    disp = func_padWidthPixels(pos,disp,extrapOpts.padWidthPixels,extrapOpts.padWidthMethod);
end

%--------------------------------------------------------------------------
% FULL-FIELD:  Temporal Smoothing of Full-Field Displacements to Create 
% Accel/Strain Rate Maps
if smoothingOpts.ffTempSmooth
    frameSize = smoothingOpts.ffTemporalKernal(1);
    polyOrder = smoothingOpts.ffTemporalKernal(2);
    disp.tSmooth.x = sgolayfilt(disp.x,polyOrder,frameSize,[],3);
    disp.tSmooth.y = sgolayfilt(disp.y,polyOrder,frameSize,[],3);
else
    disp.tSmooth.x = disp.x;
    disp.tSmooth.y = disp.y;
end

%--------------------------------------------------------------------------
% FULL-FIELD:  Spatial Smoothing Over the Full-Field
if smoothingOpts.spatSmooth == true
    if strcmp(smoothingOpts.spatialFilt,'median')
        medFiltKernel = [smoothingOpts.spatialKernal(1),smoothingOpts.spatialKernal(1)];
        for f = 1:time.numFrames 
           disp.sSmooth.x(:,:,f) = medfilt2(disp.x(:,:,f),medFiltKernel,smoothingOpts.spatialEdgeMode); 
           disp.sSmooth.y(:,:,f) = medfilt2(disp.y(:,:,f),medFiltKernel,smoothingOpts.spatialEdgeMode);
           disp.tsSmooth.x(:,:,f) = medfilt2(disp.tSmooth.x(:,:,f),medFiltKernel,smoothingOpts.spatialEdgeMode);
           disp.tsSmooth.y(:,:,f) = medfilt2(disp.tSmooth.y(:,:,f),medFiltKernel,smoothingOpts.spatialEdgeMode);
        end
    else   
        kernelRadius = ceil(2*smoothingOpts.spatialKernal(1));
        kernelLength = (2*kernelRadius)+1;    
        gaussFilter = fspecial('gaussian',[kernelLength,kernelLength],kernelRadius);
        for f = 1:time.numFrames 
            disp.sSmooth.x(:,:,f) = imfilter(disp.x(:,:,f),gaussFilter,smoothingOpts.spatialEdgeMode);
            disp.sSmooth.y(:,:,f) = imfilter(disp.y(:,:,f),gaussFilter,smoothingOpts.spatialEdgeMode);
            disp.tsSmooth.x(:,:,f) = imfilter(disp.tSmooth.x(:,:,f),gaussFilter,smoothingOpts.spatialEdgeMode);
            disp.tsSmooth.y(:,:,f) = imfilter(disp.y(:,:,f),gaussFilter,smoothingOpts.spatialEdgeMode);
        end   
    end
end

%--------------------------------------------------------------------------
% FULL-FIELD: Calculate Full-Field Acceleration Maps
if strcmp(smoothingOpts.ffAccelSmooth,'temporal')  
    [~,accel] = func_calcFFAccelFromDisp(time,disp.tSmooth);
elseif strcmp(smoothingOpts.ffAccelSmooth,'temporalspatial')
    [~,accel] = func_calcFFAccelFromDisp(time,disp.tsSmooth);
else
    [~,accel] = func_calcFFAccelFromDisp(time,disp);
end

%--------------------------------------------------------------------------
% FULL-FIELD: Calculate Strains from Full-Field Smoothed Displacement
if smoothingOpts.spatSmooth == true
    strain = func_calcStrainFromDisp(disp.sSmooth,grid.mPerPx,grid.mPerPx);
else
    strain = func_calcStrainFromDisp(disp,grid.mPerPx,grid.mPerPx);
end

%--------------------------------------------------------------------------
% WIDTH AVG: Average Data Over the Width
% Pass the raw (unsmoothed) displacement data and spatially smoothed strains
[disp,strain] = func_avgFFDataOverWidth(time,disp,strain);

%--------------------------------------------------------------------------
% WIDTH AVG: Temporal Smoothing of Width Averaged Data
if smoothingOpts.tempSmooth == true
    frameSize = smoothingOpts.temporalKernal(1);
    polyOrder = smoothingOpts.temporalKernal(2);
    disp.tSmooth.xAvg = sgolayfilt(disp.xAvg,polyOrder,frameSize,[],2);
end

%--------------------------------------------------------------------------
% WIDTH AVG: Calculate the Width Averaged Acceleration
if smoothingOpts.tempSmooth == true
    [vel.xAvg,accel.xAvg] = func_calcAccelFromDisp(time,disp.tSmooth);
else
    [vel.xAvg,accel.xAvg] = func_calcAccelFromDisp(time,disp);
end

%--------------------------------------------------------------------------
% WIDTH AVG: Extrapolate the Data to Regain the Missing Pitch
if extrapOpts.extrapLength == true
    [disp,~,accel,strain] = func_spatialExtrapAvgDataOpts(extrapOpts,time,...
        pos,disp,vel,accel,strain);
end

%--------------------------------------------------------------------------
% Remove Fields from Structures to Save RAM
if saveRAM
    % Remove the smoothed 'y' displacement fields as they are unsused here
    disp.tSmooth = rmfield(disp.tSmooth,'y');
    disp.sSmooth = rmfield(disp.sSmooth,'y');
    disp.tsSmooth = rmfield(disp.tsSmooth,'y');
    % Remove the surface average accel as this is basically the stress
    accel = rmfield(accel,'surfAvg');
    % Remove the width averaged shear strain 
    strain = rmfield(strain,{'sAvg'});
end

end