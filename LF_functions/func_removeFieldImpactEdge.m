function [pos,accel,strain] = func_removeFieldImpactEdge(VFOpts,pos,accel,strain)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 11/12/2017

% Crop all fields - assumes impact edge is on the far right hand side of
% the field in the +'ve x direction
pos.x = pos.x(1:end-VFOpts.cutEdgePx);
pos.xGrid = pos.xGrid(:,1:end-VFOpts.cutEdgePx);
pos.yGrid = pos.yGrid(:,1:end-VFOpts.cutEdgePx);
accel.x = accel.x(:,1:end-VFOpts.cutEdgePx,:);
accel.y = accel.y(:,1:end-VFOpts.cutEdgePx,:);
strain.x = strain.x(:,1:end-VFOpts.cutEdgePx,:);
strain.y = strain.y(:,1:end-VFOpts.cutEdgePx,:);
strain.s = strain.s(:,1:end-VFOpts.cutEdgePx,:);

end

