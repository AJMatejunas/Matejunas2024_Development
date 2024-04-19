function strainRate = func_smoothCalcStrainRate_v4(...
    time,grid,strain,smoothOpts,printToCons)
% TODO
%
% Author: Lloyd Fletcher
% Image-Based Mechanics Group (IBMG), University of Southampton
% Date Created: 29/7/2020 - Major update (v4)! fully overhauled edge 
% extrapolation and spatial smoothing code.
% Date Edited: 31/7/2020

% Allow printing to the console by default
if nargin < 5
    printToCons = true;
end

smoothOptsSR = smoothOpts.strainRate;
%--------------------------------------------------------------------------
% 1) Temporally smooth the strains if needed
if smoothOpts.strainRate.temporalSmooth
    strain.xtSmooth = func_reshapeAndGolayFilt3D(strain.x,...
        smoothOptsSR.temporalKernelOrder(1),smoothOptsSR.temporalKernelSize(1));
    strain.ytSmooth = func_reshapeAndGolayFilt3D(strain.y,...
        smoothOptsSR.temporalKernelOrder(2),smoothOptsSR.temporalKernelSize(2));
    strain.stSmooth = func_reshapeAndGolayFilt3D(strain.s,...
        smoothOptsSR.temporalKernelOrder(3),smoothOptsSR.temporalKernelSize(3));
else
    strain.xtSmooth = strain.x;
    strain.ytSmooth = strain.y;
    strain.stSmooth = strain.s;
end

%--------------------------------------------------------------------------
% 2) Temporally differentiate the strains to get the strain rate
[~,~,strainRate.x] = gradient(strain.xtSmooth,1,1,time.step);
[~,~,strainRate.y] = gradient(strain.ytSmooth,1,1,time.step);
[~,~,strainRate.s] = gradient(strain.stSmooth,1,1,time.step);

%--------------------------------------------------------------------------
% 3) Calculate Width Averaged Strain Rate
if printToCons
    fprintf('\tAveraging strain rate over the width.\n')
end
strainRate.xAvg = func_avgFFVarOverWidth(strainRate.x);
strainRate.sAvg = func_avgFFVarOverWidth(strainRate.s);

%--------------------------------------------------------------------------
% 4) Calculate Strain Weighted Strain Rate and Max Strain Rate
% TODO: check the functions below and update for v4

strainRate = func_calcMaxStrainRate_v4(smoothOpts,grid,strainRate);

[yR,xR,tR] = func_calcBoundsAvoidEdgeEffects_v4(strainRate.x,grid,smoothOpts);
strainRate.xWx = func_calcStrainWeightedStrainRate(strain.x,strainRate.x,yR,xR,tR);
strainRate.sWs = func_calcStrainWeightedStrainRate(strain.s,strainRate.s,yR,xR,tR);

if printToCons
    fprintf('\tSTRAIN RATE calculation complete.\n')
end

end