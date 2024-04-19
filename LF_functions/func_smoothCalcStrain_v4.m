function [strain,disp] = func_smoothCalcStrain_v4(pos,...
    time,disp,smoothOptsS,extrapOptsS,printToCons)
% Spatially differentiates displacement fields to obtain strain fields.
% Various smoothin filters can be applied to the displacement fields prior
% to numerical differentiation. Spatial differentiation is performed using
% the 'gradient' function with a centred finite difference. Also supports
% various options for directly cropping and extrapolating the strain fields
% to replace 'bad' edge data from smoothing edge effects.
%
% Author: Lloyd Fletcher
% Image-Based Mechanics Group (IBMG), University of Southampton
% Date Created: 29/7/2020 - Major update (v4)! fully overhauled edge 
% extrapolation and spatial smoothing code.
% Date Edited: 29/7/2020

% Allow printing to the console by default
if nargin < 6
    printToCons = true;
end
debug = false;

%--------------------------------------------------------------------------
% 1) Spatially smooth the displacements
if smoothOptsS.spatialSmooth == true
    if printToCons
    fprintf('\tSpatial smoothing is enabled.\n')
    fprintf('\t\tSpatial smoothing type is: %s.\n',smoothOptsS.spatialType)
    fprintf('\t\tSpatial smoothing algorithm is: %s.\n', ...
        smoothOptsS.spatialAlgorithm)
    fprintf('\t\tSpatial smoothing kernel is [X,Y]: [%i,%i].\n',...
        smoothOptsS.spatialKernelSize(1),smoothOptsS.spatialKernelSize(2))
    end
    
    if strcmp(smoothOptsS.spatialAlgorithm,'average')
        spatFilt = fspecial('average',smoothOptsS.spatialKernelSize);
        for ff = 1:time.numFrames 
           disp.sSmooth.x(:,:,ff) = imfilter(disp.x(:,:,ff),...
               spatFilt,smoothOptsS.spatialEdgeMode); 
           disp.sSmooth.y(:,:,ff) = imfilter(disp.y(:,:,ff),...
               spatFilt,smoothOptsS.spatialEdgeMode);
        end
    elseif strcmp(smoothOptsS.spatialAlgorithm,'median')
        for ff = 1:time.numFrames 
           disp.sSmooth.x(:,:,ff) = medfilt2(disp.x(:,:,ff),...
               smoothOptsS.spatialKernelSize,smoothOptsS.spatialEdgeMode); 
           disp.sSmooth.y(:,:,ff) = medfilt2(disp.y(:,:,ff),...
               smoothOptsS.spatialKernelSize,smoothOptsS.spatialEdgeMode);
        end
    else % Gaussian filter using imgaussfilt
        for ff = 1:time.numFrames 
           disp.sSmooth.x(:,:,ff) = imgaussfilt(disp.x(:,:,ff),...
               smoothOptsS.spatialKernelStd,...
               'FilterSize',smoothOptsS.spatialKernelSize,...
               'Padding',smoothOptsS.spatialEdgeMode); 
           disp.sSmooth.y(:,:,ff) = imgaussfilt(disp.y(:,:,ff),...
               smoothOptsS.spatialKernelStd,...
               'FilterSize',smoothOptsS.spatialKernelSize,...
               'Padding',smoothOptsS.spatialEdgeMode);
        end
    end
else
    if printToCons
        fprintf('\tSpatial smoothing is disabled.\n')
        %disp.sSmooth.x = disp.x;
        %disp.sSmooth.y = disp.y;
    end
end

%--------------------------------------------------------------------------
% 1.1) Temporally smooth the displacements
% TODO: implement optional temporal smoothing of displacements prior to 
% strain calculation - useful for plasticity because of the use of the
% strain increment

%--------------------------------------------------------------------------
% 2) Differentiate in space to calculate strains
if printToCons
    fprintf('\tSpatially differentiating to obtain strains from smoothed displacement fields.\n')
end
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% CO-ORDINATE SYSTEM UPDATE:
% Changed to a negative y step below because of how the displacement field
% is held in memory and gradient acts down along columns which is not
% consistent with a y axis pointing upwards. So a -y step is required

%NOTE Andrew removed the coordinate system change because strains were
%negative
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if smoothOptsS.spatialSmooth == true
    [strain.x,strain.y,strain.s] = func_calcStrainFromDisp(...
        disp.sSmooth.x,disp.sSmooth.y,pos.xStep,pos.yStep);
else
    [strain.x,strain.y,strain.s] = func_calcStrainFromDisp(...
        disp.x,disp.y,pos.xStep,pos.yStep);
end

%--------------------------------------------------------------------------
% 2.1) Crop large extrapolation area back to original size
strain.x = strain.x(disp.rY,disp.rX,:);
strain.y = strain.y(disp.rY,disp.rX,:);
strain.s = strain.s(disp.rY,disp.rX,:);
% Store the strain 'as is' for later comparison
% strain0 = strain;

%--------------------------------------------------------------------------
% 3) Optional cropping and direct extrapolation of strain data post-smooth
if extrapOptsS.postExtrapOn
    if printToCons
        fprintf('\tOptional strain crop and extrapolation enabled.\n')
        fprintf('\tOptional strain crop and extrapolation processing...\n')
    end
    
    %----------------------------------------------------------------------
    % CROPPING
    [strainCrop,posCropX,posCropY,posCropS] = func_cropStrainFields(...
        extrapOptsS,pos,strain);
    
    %----------------------------------------------------------------------
    % EXTRAPOLATION
    strain = func_extrapolateStrainFields(...
        extrapOptsS,pos,strainCrop,posCropX,posCropY,posCropS,debug);
    clear strainCrop
end

%--------------------------------------------------------------------------
% 3.1) Optional extrapolation of shear strain to zero on free edges
if strcmp(extrapOptsS.enforceGlobBCs,'shear')
    if printToCons
        fprintf('\tExtrapolating shear strain components to zero to enforce BCs.\n')
    end
    % Extrapolate shear strains to be zero on the top and bottom edges to
    % enforce the free edge boundary condition
    strain.s = func_extrapFieldYEdgesToBCs(pos,...
        extrapOptsS.enforceGlobBCsPx,strain.s);
    
    % Extrapolate the shear strains to be zero at the free edge to
    % enforce the free edge boundary condition
    strain.s = func_extrapFieldFreeEdgeToBCs(pos,...
        extrapOptsS.enforceGlobBCsPx,strain.s);
end

%--------------------------------------------------------------------------
% 4) CALC WIDTH AVG: Average strain over the width
if printToCons
    fprintf('\tAveraging smoothed strains over the width.\n')
end
strain.xAvg = func_avgFFVarOverWidth(strain.x);         
strain.yAvg = func_avgFFVarOverWidth(strain.y);
strain.sAvg = func_avgFFVarOverWidth(strain.s);

if printToCons
    fprintf('\tSTRAIN calculation complete.\n')
end

end