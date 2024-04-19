function outField = func_extrapFieldYEdgesToBCs(pos,...
    extrapPx,inField)
% Authors: Lloyd Fletcher, Jared Van-Blitterswyk
% PhotoDyn Group, University of Southampton
% Date Created:  6/8/2019
% Date Edited: 6/8/2019
% Date Edited: 29/7/2020 - updated to conform to new (v4) strain calculation
% method.
%
% Takes a 3D kinematic field with dimensions [y,x,t] and extrapolates the
% top and bottom edges to zero. Assumes an IBII like sample configuration
% with the left edge being free and the right edge being impacted, top and
% bottom edges are assumed to be free. Useful for enforcing the free edge
% conditions on the 's' (shear) strains.

% Get the size of the field we are going to process
[sy,sx,st] = size(inField);

% Initiliase the field we are going to return
outField = inField;
% Get the index of the last point of 'good' data
goodEdgeIndY1 = 1+extrapPx;
goodEdgeIndY2 = sy-extrapPx-1;

% Step through all frames in time and all slices in x
for ff = 1:st
    for xx = 1:sx
        % Calculate the slope of the extrapolation based on the last good
        % data point and zero at the top and bottom edge
        slopeEdge1 = inField(goodEdgeIndY1,xx,ff)/pos.y(goodEdgeIndY1);
        slopeEdge2 = inField(goodEdgeIndY2,xx,ff)/(pos.y(goodEdgeIndY2)-pos.lengthY);
        % Use the calculated slope to replace the values in the kinematic
        % field array that will be returned
        for e1 = 1:goodEdgeIndY1
            outField(e1,xx,ff) = slopeEdge1*pos.y(e1) + 0;
        end
        for e2 = goodEdgeIndY2:sy
            outField(e2,xx,ff) = slopeEdge2*pos.y(e2) + (-slopeEdge2*pos.lengthY);
        end
    end
end

end

