function [interpVar,Xq,Yq] = func_interp2DFEToCentroid(inputVar,X,Y,elemSize,interpMethod)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 9/8/2017
%
% Interpolates nodal FE data to the element centroid. Can also be used to
% interpolate FE data to the pixel centroid for comparison with image
% deformation data.

    % Get the size of the input arrays
    [sy,sx,st] = size(inputVar.x);
    
    % Create the mesh grid of query points
    xql = elemSize/2:elemSize:(sx-1)*elemSize;
    yql = elemSize/2:elemSize:(sy-1)*elemSize;
    [Xq,Yq] = meshgrid(xql,yql);
    
    % Interpolate frame by frame, if the shear component exists interpolate
    % this as well
    if isfield(inputVar,'s')
        for i = 1:st
            interpVar.x(:,:,i) = interp2(X,Y,squeeze(inputVar.x(:,:,i)),Xq,Yq,interpMethod);
            interpVar.y(:,:,i) = interp2(X,Y,squeeze(inputVar.y(:,:,i)),Xq,Yq,interpMethod);
            interpVar.s(:,:,i) = interp2(X,Y,squeeze(inputVar.s(:,:,i)),Xq,Yq,interpMethod);
        end
    else
        for i = 1:st
            interpVar.x(:,:,i) = interp2(X,Y,squeeze(inputVar.x(:,:,i)),Xq,Yq,interpMethod);
            interpVar.y(:,:,i) = interp2(X,Y,squeeze(inputVar.y(:,:,i)),Xq,Yq,interpMethod);
        end
    end
end

