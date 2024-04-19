function identStrength = func_getMaxStressOverFOV(grid,smoothingOpts,stress,identStrength)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 1/5/2019
% Date Edited: 1/5/2019
%
% Calculates the maximum stress values 'seen' by the specimen over the
% whole field of view and time period. Excludes data within a half
% grid pitch of the sample edge and within a half smoothing kernel in time.
    
    %----------------------------------------------------------------------
    % Specify the limits over which the maximum values can be taken
    %identStrength.FOVlims.rangeY = grid.pxPerPeriod+1:size(stress.xLinearGauge,1)-grid.pxPerPeriod-1;
    identStrength.FOVlims.rangeX = grid.pxPerPeriod+1:size(stress.xLinearGauge,2)-grid.pxPerPeriod-1;
    identStrength.FOVlims.rangeT = ceil(smoothingOpts.WATemporalKernal(1)/2):...
        size(stress.xLinearGauge,3)-ceil(smoothingOpts.WATemporalKernal(1)/2);
    
    %----------------------------------------------------------------------
    % MAX/MIN OF STANDARD STRESS GAUGE
    % Get the maximum and minimum stress guage values and their locs in
    % space and time over the whole FOV
    temp = stress.xAvg(identStrength.FOVlims.rangeX,identStrength.FOVlims.rangeT);
    
    [identStrength.FOVlims.SGMax,ind] = max(temp(:));
    [identStrength.FOVlims.SGMaxLocX,identStrength.FOVlims.SGMaxLocT]...
        = ind2sub(size(temp),ind);
    [identStrength.FOVlims.SGMin,ind] = min(temp(:));
    [identStrength.FOVlims.SGMinLocX,identStrength.FOVlims.SGMinLocT]...
        = ind2sub(size(temp),ind);
    
    % Correct for the pixels/frames cut from the edges
    identStrength.FOVlims.SGMaxLocX = identStrength.FOVlims.SGMaxLocX+identStrength.FOVlims.rangeX(1)-1;
    identStrength.FOVlims.SGMaxLocT = identStrength.FOVlims.SGMaxLocT+identStrength.FOVlims.rangeT(1)-1;
    identStrength.FOVlims.SGMinLocX = identStrength.FOVlims.SGMinLocX+identStrength.FOVlims.rangeX(1)-1;
    identStrength.FOVlims.SGMinLocT = identStrength.FOVlims.SGMinLocT+identStrength.FOVlims.rangeT(1)-1;
    
    %----------------------------------------------------------------------
    % MAX/MIN OF LINEAR STRESS GAUGE
    % Get the maximum and minimum stress guage values and their locs in
    % space and time over the whole FOV
    temp = stress.xLinearGauge(:,identStrength.FOVlims.rangeX,identStrength.FOVlims.rangeT);
    
    [identStrength.FOVlims.LSGMax,ind] = max(temp(:));
    [identStrength.FOVlims.LSGMaxLocY,identStrength.FOVlims.LSGMaxLocX,...
        identStrength.FOVlims.LSGMaxLocT] = ind2sub(size(temp),ind);
    [identStrength.FOVlims.LSGMin,ind] = min(temp(:));
    [identStrength.FOVlims.LSGMinLocY,identStrength.FOVlims.LSGMinLocX,...
        identStrength.FOVlims.LSGMinLocT] = ind2sub(size(temp),ind);
    
    % Correct for the pixels/frames cut from the edges
    identStrength.FOVlims.LSGMaxLocX = identStrength.FOVlims.LSGMaxLocX+identStrength.FOVlims.rangeX(1)-1;
    identStrength.FOVlims.LSGMaxLocT = identStrength.FOVlims.LSGMaxLocT+identStrength.FOVlims.rangeT(1)-1;
    identStrength.FOVlims.LSGMinLocX = identStrength.FOVlims.LSGMinLocX+identStrength.FOVlims.rangeX(1)-1;
    identStrength.FOVlims.LSGMinLocT = identStrength.FOVlims.LSGMinLocT+identStrength.FOVlims.rangeT(1)-1;

end

