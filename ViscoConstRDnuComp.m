% This script is written to compare shear stress in a rate dependent
% Poisson's ratio with a rate independent Poisson's ratio

%% initialize
clear all
close all
clc

%% Obtain files of constitutive model stresses
[GeoFile,GeoPath]=uigetfile('*.mat','File containing geometric info');
[ConstFile,ConstPath]=uigetfile('*.mat','Constant Nu file');
[RdFile,RdPath]=uigetfile('*.mat','Rate dependent Nu file');

%% load files
fprintf('Loading geometry file \n')
load(strcat(GeoPath,'\',GeoFile),'pos','X_vec','Y_vec');

fprintf('Loading constant nu data file \n')
Cnu=load(strcat(ConstPath,'/',ConstFile),...
    'time','strain','stress','strainxy','stressxyFE')

fprintf('Loading rate dependent nu data file \n')
RDnu=load(strcat(RdPath,'/',RdFile),...
    'time','strain','stress','strainxy','stressxyFE','StressModel')

fprintf('Files Loaded \n')

%% Define three indexes to plot
ind1=140;
ind2=420;
ind3=560;

lgd1=num2str(X_vec(ind1)*1000);
lgd2=num2str(X_vec(ind2)*1000);
lgd3=num2str(X_vec(ind3)*1000);

%% plot stress time curves
figure(1)
plot(Cnu.time.vec,Cnu.stressxyFE(ind1,:),'-k')
hold on
plot(RDnu.time.vec,RDnu.stressxyFE(ind1,:),'-b')
plot(RDnu.time.vec,RDnu.StressModel.Avxy(ind1,:),'--')
hold off

title(strcat('Shear Stress at x=',lgd1,'mm'))
xlabel('time (s)')
ylabel('Shear Stress (Pa)')
legend('ABQ constant \nu','ABQ RD \nu','Model RD \nu')

saveas(gcf,strcat('ConstRDnu_CompTime_',lgd1,'mm.fig'))
saveas(gcf,strcat('ConstRDnu_CompTime_',lgd1,'mm.png'))

%% plot stress-strain curves
figure(2)
plot(Cnu.strainxy(ind1,:),Cnu.stressxyFE(ind1,:),'-k')
hold on
plot(RDnu.strainxy(ind1,:),RDnu.stressxyFE(ind1,:),'-b')
plot(RDnu.strainxy(ind1,:),RDnu.StressModel.Avxy(ind1,:),'--')
hold off

title(strcat('Shear Stress at x=',lgd1,'mm'))
xlabel('Shear Strain')
ylabel('Shear Stress (Pa)')
legend('ABQ constant \nu','ABQ RD \nu','Model RD \nu')

saveas(gcf,strcat('ConstRDnu_CompStrain_',lgd1,'mm.fig'))
saveas(gcf,strcat('ConstRDnu_CompStrain_',lgd1,'mm.png'))


%% plot stress time curves
figure(1)
plot(Cnu.time.vec,Cnu.stressxyFE(ind3,:),'-k')
hold on
plot(RDnu.time.vec,RDnu.stressxyFE(ind3,:),'-b')
plot(RDnu.time.vec,RDnu.StressModel.Avxy(ind3,:),'--')
hold off

title(strcat('Shear Stress at x=',lgd3,'mm'))
xlabel('time (s)')
ylabel('Shear Stress (Pa)')
legend('ABQ constant \nu','ABQ RD \nu','Model RD \nu')

saveas(gcf,strcat('ConstRDnu_CompTime_',lgd3,'mm.fig'))
saveas(gcf,strcat('ConstRDnu_CompTime_',lgd3,'mm.png'))

%% plot stress-strain curves
figure(2)
plot(Cnu.strainxy(ind3,:),Cnu.stressxyFE(ind3,:),'-k')
hold on
plot(RDnu.strainxy(ind3,:),RDnu.stressxyFE(ind3,:),'-b')
plot(RDnu.strainxy(ind3,:),RDnu.StressModel.Avxy(ind3,:),'--')
hold off

title(strcat('Shear Stress at x=',lgd3,'mm'))
xlabel('Shear Strain')
ylabel('Shear Stress (Pa)')
legend('ABQ constant \nu','ABQ RD \nu','Model RD \nu')

saveas(gcf,strcat('ConstRDnu_CompStrain_',lgd3,'mm.fig'))
saveas(gcf,strcat('ConstRDnu_CompStrain_',lgd3,'mm.png'))