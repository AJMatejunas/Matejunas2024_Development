function [ConstAVxx,ConstAVxy,xxSG,ShearSG] =...
    func_TempDS(ConstAVxx,ConstAVxy,xxSG,ShearSG,DSFact)
%This function is written to Downsample stress gage and average
    %constitutive stresses in time, prior to evaluating the cost function

%Author: Andrew Matejunas

%Date Created: 2022/07/14

%Change Log/Version history:
    
%function input arguments
    %ConstAVxx- Normal Stresses in the X direction at every x point
        %averaged over y slices
    %ConstAVxy- in-plane shear stresses at every x point averaged over y
        %slices    
    %xxSG- Average stress at each x coordinate calculated using the normal
        %stress gage equation
    %ShearSG- Average stress at each coordinated calculated with the shear
        %stress gage equation

%function output arguments
    %Same as the input arguments, but temporally downsampled
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Downsample Temporally
%Constitutive normal stress
ConstAVxx_DS=downsample(ConstAVxx',DSFact)';
ConstAVxx=ConstAVxx_DS;

%Constitutive shear stress
ConstAVxy_DS=downsample(ConstAVxy',DSFact)';
ConstAVxy=ConstAVxy_DS;

%Stress Gage Normal stress
xxSG_DS=downsample(xxSG',DSFact)';
xxSG=xxSG_DS;

%Stress Gage Shear Stress
ShearSG_DS=downsample(ShearSG',DSFact)';
ShearSG=ShearSG_DS;


end

