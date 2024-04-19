function [slice11Avg,slice22Avg,slice12Avg,slice11Field,slice22Field,slice12Field] =...
    func_calcRotatedFieldAvgAlongSlice_AllComps(slice,pos,inField)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 9/2/2018
%
% Calculates the field average along an angled slice at material.rotAngle uses
% a scattered interpolant field to generate the field along the profile.
% returnComps = vector of booleans for returning 11,22,12 components
% NOTE: inField must contain components in the rotated co-ords e.g.
% inField.mat11,inField.mat22,inField.mat12. Use func_rotateStrain2D or 
% func_rotateStress2D for this purpose.

% Calculate the number of points in the field for vectorisation
nPts = length(pos.x)*length(pos.y);
[sy,~,numFrames] = size(inField.x);

% Pre-alloc for speed
slice11Avg = zeros(slice.xMaxInd,numFrames);
slice22Avg = zeros(slice.xMaxInd,numFrames);
slice12Avg = zeros(slice.xMaxInd,numFrames);
slice11Field = zeros(sy,slice.xMaxInd,numFrames);
slice22Field = zeros(sy,slice.xMaxInd,numFrames);
slice12Field = zeros(sy,slice.xMaxInd,numFrames);

for ff = 1:numFrames
    % Create the scattered interpolants for the material co-ord field for
    % this frame
    if slice.returnComps(1) % Return the 11 component averages
        field11Interp = scatteredInterpolant(reshape(pos.xGrid,nPts,1),...
                    reshape(pos.yGrid,nPts,1),reshape(inField.mat11(:,:,ff),nPts,1));
    end
    if slice.returnComps(2) % Return the 22 component averages
        field22Interp = scatteredInterpolant(reshape(pos.xGrid,nPts,1),...
                reshape(pos.yGrid,nPts,1),reshape(inField.mat22(:,:,ff),nPts,1));
    end
    if slice.returnComps(3) % Return the 12 component averages
        field12Interp = scatteredInterpolant(reshape(pos.xGrid,nPts,1),...
                reshape(pos.yGrid,nPts,1),reshape(inField.mat12(:,:,ff),nPts,1));
    end
 
    for xx = 1:slice.xMaxInd
        % Calculate the x and y co-ords of the angled slice
        [xSlice,ySlice] = func_calcAngSliceCoords(xx,pos,slice);
        
        % Calculate the field components along each angled slice
        if slice.returnComps(1) % Return the 11 component averages
            temp11 = field11Interp(xSlice,ySlice);
            slice11Avg(xx,ff) = nanmean(squeeze(temp11));
            slice11Field(:,xx,ff) = temp11;
        end
        
        if slice.returnComps(2) % Return the 22 component averages
            temp22 = field22Interp(xSlice,ySlice);
            slice22Avg(xx,ff) = nanmean(squeeze(temp22));
            slice22Field(:,xx,ff) = temp22;
        end
        
        if slice.returnComps(3) % Return the 12 component averages
            temp12 = field12Interp(xSlice,ySlice);
            slice12Avg(xx,ff) = nanmean(squeeze(temp12));
            slice12Field(:,xx,ff) = temp12;
        end
    end
end

end

