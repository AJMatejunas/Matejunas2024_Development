function [stress,strain] = func_postProcessData(time,material,grid,disp,smoothingOpts,extrapOpts)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 9/2/2017
%
% Processes output data from the grid method to determine the stress

%--------------------------------------------------------------------------
% Create Position and Time Vectors
pos.x = grid.mmPerPx/2:grid.mmPerPx:(grid.mmPerPx*size(disp.x,2));
pos.y = grid.mmPerPx/2:grid.mmPerPx:(grid.mmPerPx*size(disp.x,1));
time.numFrames = size(disp.x,3);

%--------------------------------------------------------------------------
% Pad Data to Remove NaNs over the width
if extrapOpts.padWidthNans == true
    disp = func_padWidthPixels(pos,disp,1,'nearest');
end

%--------------------------------------------------------------------------
% Spatial Smoothing
if smoothingOpts.flag == true
    if strcmp(smoothingOpts.spatialFilt,'median')
        medFiltKernel = [smoothingOpts.spatialKernal(1),smoothingOpts.spatialKernal(1)];
        for f = 1:time.numFrames 
           disp.x(:,:,f) = medfilt2(disp.x(:,:,f),medFiltKernel,'symmetric'); 
           disp.y(:,:,f) = medfilt2(disp.y(:,:,f),medFiltKernel,'symmetric'); 
        end
    else   
        kernelRadius = ceil(2*smoothingOpts.spatialKernal(1));
        kernelLength = (2*kernelRadius)+1;    
        gaussFilter = fspecial('gaussian',[kernelLength,kernelLength],kernelRadius);
        for f = 1:time.numFrames 
            disp.x(:,:,f) = imfilter(disp.x(:,:,f),gaussFilter,'replicate');
            disp.y(:,:,f) = imfilter(disp.y(:,:,f),gaussFilter,'replicate'); 
        end   
    end
end

%--------------------------------------------------------------------------
% Calculate Strains
strain = func_calcStrainFromDisp(disp,grid.mmPerPx,grid.mmPerPx);

%--------------------------------------------------------------------------
% Average Smoothed Data Over the Width
[disp,strain] = func_avgFFDataOverWidth(time,disp,strain,false);

%--------------------------------------------------------------------------
% Filter in time
if smoothingOpts.flag == true
    frameSize = smoothingOpts.temporalKernal(1);
    polyOrder = smoothingOpts.temporalKernal(2);
    disp.xAvg = sgolayfilt(disp.xAvg,polyOrder,frameSize,[],2);
end

%--------------------------------------------------------------------------
% Calculate the acceleration
[disp,vel,accel] = func_calcAccelFromDisp(time,disp);

%--------------------------------------------------------------------------
% Extrapolate the Data to Regain the Missing Pitch
if extrapOpts.flag == true
    [disp,vel,accel,strain] = func_spatialExtrapAvgData(time,pos,disp,vel,accel,...
        strain,extrapOpts.method,extrapOpts.dispPx,extrapOpts.strainPx );
end

%--------------------------------------------------------------------------
% Process Accel Using the Stress-Gauge Approach
[stress,~] = func_stressGaugeProcess(material,grid,time,pos,accel);

end