function [ShearErr,ShearLogErr,XXErr,XXLogErr]=func_SGerror(Shear_SG,SG,...
    avgXY_stress,stressX_ref,X_vec,time)
% this script is written to compare errors between reference stresses and
    %stresses measured with the shear stress gauge
    
    %% Initialize
    clear all; close all, clc
    
    %% Load stress information
    SGfile=uigetfile('Choose file with stress info','.mat');
    SGdesig=char(inputdlg('test designation'));
    load(SGfile);
    testdeg=SGdesig;
    clear SGdesig
    
    %% Compute errors in stresses in space and time
    ShearErr=abs(avgXY_stress-Shear_SG)./avgXY_stress*100;
    ShearErr(avgXY_stress==0)=0;
    
    XXErr=(stressX_ref-SG)./stressX_ref*100;
    XXErr(stressX_ref==0)=0;
   
 %% compute log of the error
ShearLogErr=log10(abs(ShearErr));
ShearLogErr(avgXY_stress==0)=0;

XXLogErr=log10(abs(ShearErr));
XXLogErr(stressX_ref==0)=0;
   
end
