% This script is written to compare shear stress calculated with
% Visco_const_verificationV1.m and  ViscoConstVer_KGdiffV1ser

%% initialize
clear all
close all
clc

%% Obtain files of constitutive model stresses
[GeoFile,GeoPath]=uigetfile('*.mat','File containing geometric info');

[ConstFile,ConstPath]=uigetfile('*.mat',...
    'calculated with Constant Nu props file');
[RdFile,RdPath]=uigetfile('*.mat','Calculated with original code and rd nu');
[KGdiffFile,KGdiffPath]=uigetfile('*.mat','Calculated with KGdiff code');
    
%% load files
fprintf('Loading geometry file \n')
load(strcat(GeoPath,'\',GeoFile),'pos','X_vec','Y_vec');

fprintf('Loading constant nu data file \n')
Cnu=load(strcat(ConstPath,'/',ConstFile),...
    'time','strain','stress','strainxy','stressxyFE','StressModel')

fprintf('Loading rate dependent nu data file with old code \n')
RDnu=load(strcat(RdPath,'/',RdFile),...
    'time','strain','stress','strainxy','stressxyFE','StressModel')

fprintf('Loading data calculated with KGdiff code \n')
KGdiff=load(strcat(KGdiffPath,'/',KGdiffFile),...
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
plot(Cnu.time.vec,Cnu.StressModel.Avxy(ind1,:),'-b')
plot(RDnu.time.vec,RDnu.StressModel.Avxy(ind1,:),'--')
plot(KGdiff.time.vec,KGdiff.StressModel.Avxy(ind1,:),':')
hold off

title(strcat('Shear Stress at x=',lgd1,'mm'))
xlabel('time (s)')
ylabel('Shear Stress (Pa)')
legend('ABQ','calculated the k=g','RD \nu old code','KGdiff',...
    'location','northwest')

saveas(gcf,strcat('OldNew_CompTime_',lgd1,'mm.fig'))
saveas(gcf,strcat('OldNew_CompTime_',lgd1,'mm.png'))

%% plot stress-strain curves
figure(2)
plot(Cnu.strainxy(ind1,:),Cnu.stressxyFE(ind1,:),'-k')
hold on
plot(Cnu.strainxy(ind1,:),Cnu.StressModel.Avxy(ind1,:),'-b')
plot(RDnu.strainxy(ind1,:),RDnu.StressModel.Avxy(ind1,:),'--')
plot(KGdiff.strainxy(ind1,:),KGdiff.StressModel.Avxy(ind1,:),':')
hold off


title(strcat('Shear Stress at x=',lgd1,'mm'))
xlabel('Shear Strain')
ylabel('Shear Stress (Pa)')
legend('ABQ','calculated the k=g','RD \nu old code','KGdiff',...
    'location','northwest')

saveas(gcf,strcat('OldNew_CompStrain_',lgd1,'mm.fig'))
saveas(gcf,strcat('OldNew_CompStrain_',lgd1,'mm.png'))


%% plot stress time curves
figure(3)
plot(Cnu.time.vec,Cnu.stressxyFE(ind3,:),'-k')
hold on
plot(Cnu.time.vec,Cnu.StressModel.Avxy(ind3,:),'-b')
plot(RDnu.time.vec,RDnu.StressModel.Avxy(ind3,:),'--')
plot(KGdiff.time.vec,KGdiff.StressModel.Avxy(ind3,:),':')
hold off

title(strcat('Shear Stress at x=',lgd3,'mm'))
xlabel('time (s)')
ylabel('Shear Stress (Pa)')
legend('ABQ','calculated the k=g','RD \nu old code','KGdiff',...
    'location','northwest')


saveas(gcf,strcat('OldNew_CompTime_',lgd3,'mm.fig'))
saveas(gcf,strcat('OldNew_CompTime_',lgd3,'mm.png'))

%% plot stress-strain curves
figure(4)
plot(Cnu.strainxy(ind3,:),Cnu.stressxyFE(ind3,:),'-k')
hold on
plot(Cnu.strainxy(ind3,:),Cnu.StressModel.Avxy(ind3,:),'-b')
plot(RDnu.strainxy(ind3,:),RDnu.StressModel.Avxy(ind3,:),'--')
plot(KGdiff.strainxy(ind3,:),KGdiff.StressModel.Avxy(ind3,:),':')
hold off


title(strcat('Shear Stress at x=',lgd3,'mm'))
xlabel('Shear Strain')
ylabel('Shear Stress (Pa)')
legend('ABQ','calculated the k=g','RD \nu old code','KGdiff',...
    'location','northwest')


saveas(gcf,strcat('OldNew_CompStrain_',lgd3,'mm.fig'))
saveas(gcf,strcat('OldNew_CompStrain_',lgd3,'mm.png'))



%% plot stress time curves
figure(5)
plot(Cnu.time.vec,Cnu.stressxyFE(ind2,:),'-k')
hold on
plot(Cnu.time.vec,Cnu.StressModel.Avxy(ind2,:),'-b')
plot(RDnu.time.vec,RDnu.StressModel.Avxy(ind2,:),'--')
plot(KGdiff.time.vec,KGdiff.StressModel.Avxy(ind2,:),':')
hold off

title(strcat('Shear Stress at x=',lgd3,'mm'))
xlabel('time (s)')
ylabel('Shear Stress (Pa)')
legend('ABQ','calculated the k=g','RD \nu old code','KGdiff',...
    'location','northwest')


saveas(gcf,strcat('OldNew_CompTime_',lgd2,'mm.fig'))
saveas(gcf,strcat('OldNew_CompTime_',lgd2,'mm.png'))

%% plot stress-strain curves
figure(6)
plot(Cnu.strainxy(ind2,:),Cnu.stressxyFE(ind2,:),'-k')
hold on
plot(Cnu.strainxy(ind2,:),Cnu.StressModel.Avxy(ind2,:),'-b')
plot(RDnu.strainxy(ind2,:),RDnu.StressModel.Avxy(ind2,:),'--')
plot(KGdiff.strainxy(ind2,:),KGdiff.StressModel.Avxy(ind2,:),':')
hold off


title(strcat('Shear Stress at x=',lgd3,'mm'))
xlabel('Shear Strain')
ylabel('Shear Stress (Pa)')
legend('ABQ','calculated the k=g','RD \nu old code','KGdiff',...
    'location','northwest')


saveas(gcf,strcat('OldNew_CompStrain_',lgd2,'mm.fig'))
saveas(gcf,strcat('OldNew_CompStrain_',lgd2,'mm.png'))