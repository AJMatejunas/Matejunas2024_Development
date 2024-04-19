function [ genSSCurves ] = func_genSSCurvesIso(pos,material,strain,accel)
% Author: L. Fletcher
% PhotoDyn Research Group
% Date Created: 7/12/2018 
% Date Edited: 7/12/2018
% 
% Calculates the generalised stress-strain curves for a  linear elastic
% isotropic constituitive law. These equations are based on the virtual
% fields method and rigid body virtual fields. The equations give a
% formulation that allows one to plot some weighted average of strain
% against a weighted average of the acceleration to obtain Q11 and Q12.
%
% Full details can be found in the following paper:
% Pierron, F. and Fletcher, L., "Image-based generalized stress-strain 
% curves for inertial high strain rate tests on isotropic and orthotropic
% materials" 
%
% The function takes input data structs describing the kinematic fields in
% specimen and outputs a data structure containing the values of each
% equation set allowing the linear relationship for each Q component to be 
% plotted and fitted to obtain the relevant stiffness component.

%--------------------------------------------------------------------------
% Calculate field averages and weighted averages
% NOTATION: 
% 'a' acceleration, 'e' strain, following subscript notatation gives direction
% '_y' or '_x' weighted average by given coordinate
% '_l' suffix average over a line, '_S' suffix average over the surface
numFrames = size(strain.x,3);
pos.xGridF = padarray(pos.xGrid,[0,0,numFrames-1],'replicate','post');
pos.yGridF = padarray(pos.yGrid,[0,0,numFrames-1],'replicate','post');
pos.x0F = squeeze(padarray(pos.x,[0,0,numFrames-1],'replicate','post'));
ax_S = func_calcSurfAvgFromFreeEdge(accel.x);
ay_S = func_calcSurfAvgFromFreeEdge(accel.y);
ax_y_S = func_calcSurfAvgFromFreeEdge(accel.x.*pos.yGridF);
ay_x_S = func_calcSurfAvgFromFreeEdge(accel.y.*pos.xGridF);
exx_l = func_avgFFVarOverWidth(strain.x);
eyy_l = func_avgFFVarOverWidth(strain.y);
exy_l = func_avgFFVarOverWidth(strain.s/2);
exx_y_l = func_avgFFVarOverWidth(strain.x.*pos.yGridF);
eyy_y_l = func_avgFFVarOverWidth(strain.y.*pos.yGridF);

%--------------------------------------------------------------------------
% Equations 1 + 2

% Q11
genSSCurves.Q11_strainAvg_Eq12 = (exx_l+eyy_l).*exy_l;
genSSCurves.Q11_accelAvg_Eq12 = material.rho.*pos.x0F.*(ax_S.*exy_l + ay_S.*eyy_l);

% Q12
genSSCurves.Q12_strainAvg_Eq12 = (exx_l+eyy_l).*exy_l;
genSSCurves.Q12_accelAvg_Eq12 = material.rho.*pos.x0F.*(ax_S.*exy_l - ay_S.*exx_l);

%--------------------------------------------------------------------------
% Equations 1 + 3

% Q11
genSSCurves.Q11_strainAvg_Eq13 = (-exx_y_l.*eyy_l + exx_l.*eyy_y_l) + ...
    (pos.x0F.*exy_l.*(exx_l+eyy_l));
genSSCurves.Q11_accelAvg_Eq13 = material.rho.*pos.x0F.*...
        (eyy_l.*(ay_x_S - ax_y_S) + ax_S.*(pos.x0F.*exy_l+eyy_y_l));
    
% Q12
genSSCurves.Q12_strainAvg_Eq13 = (-exx_y_l.*eyy_l + exx_l.*eyy_y_l) + ...
    (pos.x0F.*exy_l.*(exx_l+eyy_l));
genSSCurves.Q12_accelAvg_Eq13 = material.rho.*pos.x0F.*...
        (exx_l.*(ax_y_S - ay_x_S) + ax_S.*(pos.x0F.*exy_l-exx_y_l));

%--------------------------------------------------------------------------
% Equations 2 + 3

% Q11
genSSCurves.Q11_strainAvg_Eq23 = (exx_y_l + eyy_y_l).*exy_l;
genSSCurves.Q11_accelAvg_Eq23 = material.rho.*pos.x0F.*...
        (exy_l.*(ax_y_S - ay_x_S + pos.x0F.*ay_S) + ay_S.*eyy_y_l);
    
% Q12
genSSCurves.Q12_strainAvg_Eq23 = (exx_y_l + eyy_y_l).*exy_l;
genSSCurves.Q12_accelAvg_Eq23 = material.rho.*pos.x0F.*...
        (exy_l.*(ax_y_S - ay_x_S + pos.x0F.*ay_S) - ay_S.*exx_y_l);

end

