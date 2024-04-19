function identKVMan = func_VFIdentManKV(VFOpts,pos,time,material,VFs,strain,accel)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 12/9/2017
%
% Uses two manually defined virtual fields to identify the stiffness
% components Qxx and Qxy for a linear elastic isotropic material in a state
% of plane stress

% Remove the extrapolated data on the impact edge if needed
if VFOpts.cutImpactEdge
    [pos,accel,strain] = func_removeFieldImpactEdge(VFOpts,pos,accel,strain);
end

% Vectorise the virtual fields for speed
for v = 1:length(VFs)
    VFs{v}.uX = reshape(VFs{v}.uX, [], 1);
    VFs{v}.uY = reshape(VFs{v}.uY, [], 1);
    VFs{v}.epsXX = reshape(VFs{v}.epsXX, [], 1);
    VFs{v}.epsYY = reshape(VFs{v}.epsYY, [], 1);
    VFs{v}.epsXY = reshape(VFs{v}.epsXY, [], 1);
end

% Loop over each frame and calculate the modulus
for t = VFOpts.startFrame:VFOpts.endFrame
    % Vectorise the fields for speed.
    eXXvec = reshape(strain.x(:,:,t), [], 1);
    eYYvec = reshape(strain.y(:,:,t), [], 1);
    eXYvec = reshape(strain.s(:,:,t), [], 1);
    aXvec = reshape(accel.x(:,:,t), [], 1);
    aYvec = reshape(accel.y(:,:,t), [], 1);

    % Create the A and B matrix
    for v = 1:length(VFs)
        A(v,1) = -(nanmean(eXXvec.*VFs{v}.epsXX)+...
                 nanmean(eYYvec.*VFs{v}.epsYY)+...
                 nanmean(0.5*eXYvec.*VFs{v}.epsXY));
        A(v,2) = -(nanmean(eXXvec.*VFs{v}.epsYY)+...
                 nanmean(eYYvec.*VFs{v}.epsXX)-...
                 nanmean(0.5*eXYvec.*VFs{v}.epsXY));
        B(v,1) = material.rho*(nanmean(aXvec.*VFs{v}.uX)+...
                 nanmean(aYvec.*VFs{v}.uY));
            
    end

    % Invert the A matrix and solve the system for Qxx and Qxy
    stabilityCheck(t) = cond(A);
    Q = A\B;
    Qxx(t) = Q(1);
    Qxy(t) = Q(2);
end

% Calculate modulus and poissons ratio from the stiffness components
nuxy = Qxy./Qxx;
Exx = Qxx.*(1-nuxy.^2);

% Push everything into the identified properties struct and return it
identStiffVFMan.QxxVsT = Qxx;
identStiffVFMan.QxyVsT = Qxy;
identStiffVFMan.ExxVsT = Exx;
identStiffVFMan.NuxyVsT = nuxy;
identStiffVFMan.matCond = stabilityCheck;

end

