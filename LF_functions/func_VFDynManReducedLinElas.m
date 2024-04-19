function Exx = func_VFDynManReducedLinElas(VFOpts,pos,material,accel,strain)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 27/3/2017
%
% Uses manual virtual fields to calculate Qxx for a nominally uniaxial
% test. This ignores any effects coming from non-axial strains

    % Remove the extrapolated data on the impact edge if needed
    if VFOpts.cutImpactEdge
        [pos,accel,strain] = func_removeFieldImpactEdge(VFOpts,pos,accel,strain);
    end
    
    % Define reduced manual virtual fields
    LL = pos.x(end)+pos.xStep/2;
    VFOpts.ux_star = @(x) x - LL;
    VFOpts.epsx_star = @(x) 1;

    % Pre-allocation for the weighted acceleration 
    ax_fux = zeros(size(accel.x));
    epsx_fepsx = zeros(size(strain.x));
    
    % Calculate the acceleration weighted virtual field
    for tt = VFOpts.startFrame:VFOpts.endFrame
        ax_fux(:,:,tt) = accel.x(:,:,tt).*VFOpts.ux_star(pos.xGrid);
        epsx_fepsx(:,:,tt) = strain.x(:,:,tt).*VFOpts.epsx_star(pos.xGrid);
    end
    
    % Calculate the modulus
    Exx = -material.rho*squeeze(mean(mean(ax_fux)))./squeeze(mean(mean(epsx_fepsx)));
    
end

