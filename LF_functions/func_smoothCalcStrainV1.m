function [disp,strain,strainRate] = func_smoothCalcStrainV1(globalOpts,pos,...
    time,grid,disp,smoothingOpts,extrapOpts,printToCons)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 27 Mar. 2017
% Modified by: Jared Van Blitterswyk 
% Date: 15 Oct. 2017
%Modified by: Andrew Matejunas
    %Date 02 March 2023- changed to V1 to use the same version of
    %func_calcStrainFromDisp
% Processes output displacement fields from the grid method. Output 
% variables are obtained via numerical differentiation using a central 
% difference method. Variables output include strain and strain rate.

% Allow printing to the console by default
if nargin < 8
    printToCons = true;
end

%% Crop Data To Remove Edge Effects and Replace NaNs
if printToCons
    fprintf('\tCropping one pitch from edges to remove poor data from grid method processing.\n')
end
dispCrop.x = func_cropFields(disp.x,2,extrapOpts.FFDispPx);
dispCrop.y = func_cropFields(disp.y,1,extrapOpts.FFDispPx);

%--------------------------------------------------------------------------g100
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

%% FULL-FIELD:  Spatial Smoothing Over the Full-Field
if smoothingOpts.spatialSmooth == true
    if printToCons
        fprintf('\tSpatially smoothing full-field displacements.\n')
    end
    if strcmp(smoothingOpts.spatialFilt,'median')
        medFiltKernel = [smoothingOpts.spatialKernal(1),smoothingOpts.spatialKernal(1)];
        if printToCons
            fprintf('\tSpatial full-field filter is %s with a [%i,%i] kernel.\n',...
                smoothingOpts.spatialFilt,medFiltKernel(1),medFiltKernel(2))
        end
        
        dispTemp.x = padarray(dispCrop.x, [3*medFiltKernel,3*medFiltKernel],'replicate');
        dispTemp.y = padarray(dispCrop.y, [3*medFiltKernel,3*medFiltKernel],'replicate');
        
        for ff = 1:time.numFrames 
           disp.sSmooth.x(:,:,ff) = medfilt2(dispTemp.x(:,:,ff),medFiltKernel,smoothingOpts.spatialEdgeMode); 
           disp.sSmooth.y(:,:,ff) = medfilt2(dispTemp.y(:,:,ff),medFiltKernel,smoothingOpts.spatialEdgeMode);
        end
    else
        % Calculate and display the gaussian kernel size
        kernelLengthX = smoothingOpts.spatialKernal(1);
        kernelRadiusX = ceil((kernelLengthX-1)/2);
		kernelSigmaX = kernelRadiusX/2;
        
        kernelLengthY = smoothingOpts.spatialKernal(2);
        kernelRadiusY = ceil((kernelLengthY-1)/2);
		kernelSigmaY = kernelRadiusY/2;
        
        if printToCons
        fprintf('\tSpatial full-field filter is %s.\n',smoothingOpts.spatialFilt)
        fprintf('\tSpatial filter SD is [x,y] = [%.1f,%.1f].\n',kernelSigmaX,kernelSigmaY)
        fprintf('\tSpatial kernel length is [x,y] = [%i,%i].\n',kernelLengthX,kernelLengthY)
        end
         
        % Pad out the array to reconstruct the missing edge data, option to
        % extend the pad region instead of using the inbuilt padding
        % options in the matlab gaussian filter.
        padRadius = grid.pxPerPeriod;
        if isfield(extrapOpts,'extendPadBeforeSmooth')
            if extrapOpts.extendPadBeforeSmooth
                padRadius = 3*max([smoothingOpts.spatialKernal(1),smoothingOpts.spatialKernal(2),grid.pxPerPeriod]);
            end
        end
        dispTemp = func_padDispArrayExt(pos,dispCrop,grid,padRadius,extrapOpts.padEdgeMethod,extrapOpts.dispPx,extrapOpts.padFitWindow,extrapOpts); 

        % Perform spatial smoothing for each frame using the Gaussian filter 
        for ff = 1:time.numFrames 
            disp.sSmooth.x(:,:,ff) = imgaussfilt(dispTemp.x(:,:,ff),[kernelSigmaX,kernelSigmaY],'Padding',smoothingOpts.spatialEdgeMode);
            disp.sSmooth.y(:,:,ff) = imgaussfilt(dispTemp.y(:,:,ff),[kernelSigmaX,kernelSigmaY],'Padding',smoothingOpts.spatialEdgeMode);
        end

        % Crop fields back to the original size if required
        % If padRadius == extrapOpts.FFDispPx then the fields are unchanged
        disp.extrap.x = dispTemp.x(:,padRadius-extrapOpts.FFDispPx+1:end-padRadius+extrapOpts.FFDispPx,:);
        disp.extrap.y = dispTemp.y(padRadius-extrapOpts.FFDispPx+1:end-padRadius+extrapOpts.FFDispPx,:,:);
        disp.sSmooth.x = disp.sSmooth.x(:,padRadius-extrapOpts.FFDispPx+1:end-padRadius+extrapOpts.FFDispPx,:);
        disp.sSmooth.y = disp.sSmooth.y(padRadius-extrapOpts.FFDispPx+1:end-padRadius+extrapOpts.FFDispPx,:,:);
    end
else
    padRadius = extrapOpts.dispPx;
    dispTemp = func_padDispArrayExt(pos,dispCrop,grid,padRadius,extrapOpts.padEdgeMethod,extrapOpts.dispPx,extrapOpts.padFitWindow,extrapOpts); % pad out poor data at edges
    disp = dispTemp;
    disp.extrap.x = dispTemp.x(:,padRadius-extrapOpts.FFDispPx+1:end-padRadius+extrapOpts.FFDispPx,:);
    disp.extrap.y = dispTemp.y(padRadius-extrapOpts.FFDispPx+1:end-padRadius+extrapOpts.FFDispPx,:,:);
end

%% FULL-FIELD: Calculate Strains from Full-Field Smoothed Displacement
if printToCons
    fprintf('\tSpatially differentiating to obtain strains from smoothed displacement fields.\n')
end
if smoothingOpts.spatialSmooth == true
    [strainX,strainY,strainS] = func_calcStrainFromDisp(disp.sSmooth.x,disp.sSmooth.y, ...
        pos.xStep,pos.yStep);
else
    [strainX,strainY,strainS] = func_calcStrainFromDisp(dispTemp.x,dispTemp.y, ...
        pos.xStep,pos.yStep);
end
strain.x=strainX;
strain.y=strainY;
strain.s=strainS;
%% FULL-FIELD: Calculate the Full-Field Strain Rate
if strcmp(globalOpts.fieldComponents,'xOnly')
    if printToCons
        fprintf('\tCalculating [x] strain rate.\n')
    end
    [strainRate,~] = func_calculateStrainRate(strain,time,smoothingOpts);
else
    if printToCons
        fprintf('\tCalculating [x,y,s] strain rate.\n')
    end
    [strainRate,~] = func_calculateStrainRate(strain,time,smoothingOpts,true);
end

%% WIDTH AVG: Average displacement and strain over the width
if printToCons
    fprintf('\tAveraging smoothed strain and strain rate over the width.\n')
end
strain.xAvg = func_avgFFVarOverWidth(strain.x);         
strain.yAvg = func_avgFFVarOverWidth(strain.y);

if strcmp(globalOpts.fieldComponents,'xOnly')
    strainRate.xAvg = func_avgFFVarOverWidth(strainRate.x); 
else
    strainRate.xAvg = func_avgFFVarOverWidth(strainRate.x);
    strainRate.yAvg = func_avgFFVarOverWidth(strainRate.y);
    strainRate.sAvg = func_avgFFVarOverWidth(strainRate.s);
end

%% WIDTH AVG: Get the peak axial avg strain rates
strainRate = func_calcMaxStrainRate(smoothingOpts,grid,strainRate);

end