function strainRate = func_calcMaxStrainRate(smoothingOpts,grid,strainRate)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 21/2/2019
% Date Edited: 24/5/2019
%
% Calculates the max and min width average strain rate ignoring the
% smoothing edge effects.
    
    % Only take the max strain rate within a half smoothing kernel of the
    % edge
    if smoothingOpts.spatialSmooth == true
        edgePx = ceil(smoothingOpts.spatialKernal(1)/2)+grid.pxPerPeriod;
    else
        edgePx = grid.pxPerPeriod+1;
    end
    xRange = edgePx:size(strainRate.xAvg,1) - edgePx;
    
    % Only take the max strain rate avoiding temporal edge effects from the
    % filter
    if smoothingOpts.FFTempSmooth == true
        edgeFrames = ceil(smoothingOpts.FFTemporalKernal(1));
    else
        edgeFrames = 1;
    end
    tRange = edgeFrames:size(strainRate.xAvg,2) - edgeFrames;

    % Ignore the edges where the data is corrupted by smoothing effects
    strainRate.xAvgMin = min(min(strainRate.xAvg(xRange,tRange)));
    strainRate.xAvgMax = max(max(strainRate.xAvg(xRange,tRange)));
end

