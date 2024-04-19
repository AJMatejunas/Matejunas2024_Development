function [SxRot,SyRot,SsRot] = func_rotateStrain2D(Sx,Sy,Ss,alpha)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 7/2/2018
%
% Rotates strain components by angle 'alpha' 
        
    SxRot= (cosd(alpha))^2.*Sx...
            +(sind(alpha))^2.*Sy...
            +sind(alpha)*cosd(alpha).*Ss;
    SyRot=(sind(alpha))^2.*Sx...
            +(cosd(alpha))^2.*Sy...
            -sind(alpha)*cosd(alpha).*Ss; 
    SsRot=-2*sind(alpha)*cosd(alpha).*Sx...
            +2*sind(alpha)*cosd(alpha).*Sy...
            +((cosd(alpha))^2-(sind(alpha))^2).*Ss;

end
