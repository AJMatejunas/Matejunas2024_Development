function inField = func_calcRotatedFieldAvgAlongSlice_v2(slice,pos,inField)
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

nPts = length(pos.x)*length(pos.y);
numFrames = size(inField.x,3);
for ff = 1:numFrames
    fprintf('Calculating slice averages for frame %i.\n',ff)
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
    
    for ss = 1:slice.numSlices
         % Calculate the field components along each angled slice
        if slice.returnComps(1)
            temp11 = field11Interp(slice.linesXCoords{ss},slice.linesYCoords{ss});
            inField.slice11Avg(ss,ff) = nanmean(squeeze(temp11));
        end
        if slice.returnComps(2)
            temp22 = field22Interp(slice.linesXCoords{ss},slice.linesYCoords{ss});
            inField.slice22Avg(ss,ff) = nanmean(squeeze(temp22));
        end
        if slice.returnComps(3)
            temp12 = field12Interp(slice.linesXCoords{ss},slice.linesYCoords{ss});
            inField.slice12Avg(ss,ff) = nanmean(squeeze(temp12));
        end 
    end 
end

end

