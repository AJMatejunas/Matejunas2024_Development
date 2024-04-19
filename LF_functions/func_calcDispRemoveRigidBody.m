function [dispDeform,dispRigid] = func_calcDispRemoveRigidBody(dispField)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 9/3/2017
%
% Calculates the rigid body motion using the mean over the field and
% subtracts this to work out the deformation component of the displacement
% field

dispRigid = squeeze(nanmean(nanmean(dispField)));
dispDeform = zeros(size(dispField));
for f = 1:size(dispField,3)
    dispDeform(:,:,f) = dispField(:,:,f) - dispRigid(f);
end

end