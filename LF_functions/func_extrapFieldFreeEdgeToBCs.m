function outField = func_extrapFieldFreeEdgeToBCs(pos,...
    extrapPx,inField)
% Authors: Lloyd Fletcher, Jared Van-Blitterswyk
% PhotoDyn Group, University of Southampton
% Date Created:  6/8/2019
% Date Edited: 6/8/2019
% Date Edited: 29/7/2020 - updated to conform to new (v4) strain calculation
% method.
%
% Takes a 3D kinematic field with dimensions [y,x,t] and extrapolates the
% free edge to zero. Assumes an IBII like sample configuration with the 
% left edge being free and the right edge being impacted. Useful for 
% enforcing the free edge conditions on the 's' (shear) strains.

% Get the size of the field we are going to process
[sy,~,st] = size(inField);

% Initiliase the field we are going to return
outField = inField;

% Get the index of the last point of 'good' data in the field that is not
% affected by the smoothing kernel
goodEdgeIndX = extrapPx+1;

% Step through all frames in time and all slices in y
for ff = 1:st
    for yy = 1:sy
        % Calculate the slope of the extrapolation based on the last good
        % data point and zero at the free edge (LHS)
        slopeEdge = inField(yy,goodEdgeIndX,ff)/pos.x(goodEdgeIndX);
        for ee = 1:goodEdgeIndX
            outField(yy,ee,ff) = slopeEdge*pos.x(ee)+0;
        end
    end 
end

end

