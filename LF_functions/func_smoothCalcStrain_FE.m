function [disp,strain,strainRate] = func_smoothCalcStrain_FE(...
    globalOpts,pos,time,grid,disp,smoothingOpts)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 27/3/2017
% Modified by: Jared Van Blitterswyk 
% Date: 15 Oct. 2017
%
% Processes output displacement fields from the grid method. Output 
% variables are obtained via numerical differentiation using a central 
% difference method. Variables output include strain and strain rate. 

%% FULL-FIELD: Spatial Smoothing Over the Full-Field
if smoothingOpts.spatialSmooth == true
    fprintf('\tSpatially smoothing full-field displacements.\n')
    if strcmp(smoothingOpts.spatialFilt,'median')
        medFiltKernel = [smoothingOpts.spatialKernal(1),smoothingOpts.spatialKernal(1)];
        fprintf('\tSpatial full-field filter is %s with a [%i,%i] kernel.\n',...
            smoothingOpts.spatialFilt,medFiltKernel(1),medFiltKernel(2))
        
        dispTemp.x = padarray(disp.x, [3*medFiltKernel,3*medFiltKernel],'replicate');
        dispTemp.y = padarray(disp.y, [3*medFiltKernel,3*medFiltKernel],'replicate');
        
        for f = 1:time.numFrames 
           disp.sSmooth.x(:,:,f) = medfilt2(dispTemp.x(:,:,f),medFiltKernel,smoothingOpts.spatialEdgeMode); 
           disp.sSmooth.y(:,:,f) = medfilt2(dispTemp.y(:,:,f),medFiltKernel,smoothingOpts.spatialEdgeMode);
        end
    else
        % Calculate and display the gaussian kernel size
        kernelLengthX = smoothingOpts.spatialKernal(1);
        kernelRadiusX = ceil((kernelLengthX-1)/2);
		kernelSigmaX = kernelRadiusX/2;
        
        kernelLengthY = smoothingOpts.spatialKernal(2);
        kernelRadiusY = ceil((kernelLengthY-1)/2);
		kernelSigmaY = kernelRadiusY/2;
        
        fprintf('\tSpatial full-field filter is %s.\n',smoothingOpts.spatialFilt)
        fprintf('\tSpatial filter SD is [x,y] = [%.1f,%.1f].\n',kernelSigmaX,kernelSigmaY)
        fprintf('\tSpatial kernel length is [x,y] = [%i,%i].\n',kernelLengthX,kernelLengthY)
        
        for f = 1:time.numFrames 
         	disp.sSmooth.x(:,:,f) = imgaussfilt(disp.x(:,:,f),[kernelSigmaX,kernelSigmaY],'Padding',smoothingOpts.spatialEdgeMode);
			disp.sSmooth.y(:,:,f) = imgaussfilt(disp.y(:,:,f),[kernelSigmaX,kernelSigmaY],'Padding',smoothingOpts.spatialEdgeMode);
		end

    end
end

%% FULL-FIELD: Calculate Strains from Full-Field Smoothed Displacement
fprintf('\tSpatially differentiating to obtain strains from smoothed displacement fields.\n')
if smoothingOpts.spatialSmooth == true
    strain = func_calcStrainFromDisp(disp.sSmooth,pos.xStep,pos.yStep);
else
    strain = func_calcStrainFromDisp(disp,pos.xStep,pos.yStep);
end

%% FULL-FIELD: Calculate the Full-Field Strain Rate
if strcmp(globalOpts.fieldComponents,'xOnly')
    fprintf('\tCalculating [x] strain rate.\n')
    [strainRate,~] = func_calculateStrainRate(strain,time,smoothingOpts);
else
    fprintf('\tCalculating [x,y,s] strain rate.\n')
    [strainRate,~] = func_calculateStrainRate(strain,time,smoothingOpts,true);
end

%% WIDTH AVG: Average displacement and strain over the width
fprintf('\tAveraging smoothed strain and strain rate over the width.\n')
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
if smoothingOpts.spatialSmooth == true
    edgePx = round(smoothingOpts.spatialKernal(1)/2)+1;
else
    edgePx = 1;
end

xRange = edgePx:size(strainRate.xAvg,1) - edgePx;

% Ignore the edges where the data is corrupted by smoothing effects
strainRate.xAvgMin = min(min(strainRate.xAvg(xRange,:)));
strainRate.xAvgMax = max(max(strainRate.xAvg(xRange,:)));

end