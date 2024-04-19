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
   
   %% Plot heat maps of the error in stresses in space and time
   
figure
contourf(X_vec,time.vec,ShearErr')
xlabel('X coordinate')
ylabel('time (s)')
c=colorbar;
c.Label.String='Shear error (\%)';
saveas(gcf,strcat(testdeg,'ShearError_map'))
saveas(gcf,strcat(testdeg,'ShearError_map.png'))

figure
contourf(X_vec,time.vec,XXErr')
xlabel('X coordinate')
ylabel('time (s)')
c=colorbar;
c.Label.String='X_stress error (\%)';
saveas(gcf,strcat(testdeg,'XXError_map'))
saveas(gcf,strcat(testdeg,'XXError_map.png'))
   
   %% Plot heat maps of the log_10 of the stress error in space and time
   figure
contourf(X_vec,time.vec,ShearLogErr')
xlabel('X coordinate')
ylabel('time (s)')
c=colorbar;
c.Label.String='log(Shear error (%))';
saveas(gcf,strcat(testdeg,'LogXYError_map'))
saveas(gcf,strcat(testdeg,'LogXYError_map.png'))

figure
contourf(X_vec,time.vec,XXLogErr')
xlabel('X coordinate')
ylabel('time (s)')
c=colorbar;
c.Label.String='log(X stress error (%))';
saveas(gcf,strcat(testdeg,'LogXXError_map'))
saveas(gcf,strcat(testdeg,'LogXXError_map.png'))


