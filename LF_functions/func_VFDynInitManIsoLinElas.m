function VFs = func_VFDynInitManIsoLinElas(VFOpts,pos,accel,strain)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 13/9/2017
%
% Creates two virtual fields for stiffness identification in a linear
% elastic isotropic material

% Remove the extrapolated data on the impact edge if needed
if VFOpts.cutImpactEdge
    [pos,~,~] = func_removeFieldImpactEdge(VFOpts,pos,accel,strain);
end

% Geometric Variables to create the fields
L = pos.x(end)+pos.xStep/2;
H = pos.y(end)+pos.yStep/2;
x = pos.xGrid;
y = pos.yGrid;

% First Virtual Field
VFs{1}.uX = L-x;
VFs{1}.uY = zeros(size(x));
VFs{1}.epsXX = -1.*ones(size(x));
VFs{1}.epsYY = zeros(size(x));
VFs{1}.epsXY = zeros(size(x));

% Second Virtual Field
VFs{2}.uX = zeros(size(x));
VFs{2}.uY = (L-x).*(y/H);
VFs{2}.epsXX = zeros(size(x));
VFs{2}.epsYY = (L-x)./H;
VFs{2}.epsXY = -y./H;


end

