function VFs = VFs=func_VFManDynKV(VFOpts,pos,accel,strain,material)
%Author: Andrew Matejunas
%date 22/06/2020

%adapept from 
% PhotoDyn Group, University of Southampton
% Date: 13/9/2017
%
% Creates two virtual fields for identification of parameters for a
% Kelvin-voigt viscoelastic solid


% the funtional form of the Kelvin-Voigt model i 3 dimensions is
%Stress_ij=LameModulus*(Strain_kk+theta_LameModulus*StrainRate_kk)*KroneckerDelta_ij+2*ShearModulus(Strain_ij+theta_shear*StrainRate_ij)

% assuming rate sensitivity parameters of the lame and shear moduli are the
% same (theta_Lame=Theta_shear)that leaves three parameters that must be
% identified
    %Lame Modulus
    %Shear Modulus
    %Rate Sensitivity parameter (Theta)

%This will require 3 linearly independent virtual fields


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

% Third virtual field
VFs{3}.uX=L^2-x.^2
VFs{3}.uY=zeros(size(x));
VFs{3}.epsXX=2*x;
VFs{3}.epsYY=0;
VFs{3}.epsXY=0;

% %alternate  third field
% % Do we lose too much data here? Do we need a field that results in shear?
% % Maybe add a virrtual field on the uY?
% VFs{3}.uX=sin((x-L)*pi/L):
% VFs{3}.uY=zeros(size(x));
% VFs{3}.epsXX=pi/L*cos((x-L)*pi/L);
% VFs{3}.epsYY=0;
% VFs{3}.epsXY=0;

end

