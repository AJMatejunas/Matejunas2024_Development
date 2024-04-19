%this code is written to evaluate the effect of pretrigger frames

%% initialize
clear all; close all; clc

%% Select data file
[SGfile,SGpath]=uigetfile('*.mat',...
    'select SG data file containing FE stresses');
SGname=strcat(SGpath,'/',SGfile);

[Pre5file,Pre5path]=uigetfile('*.mat',...
    'select bdata file containing 5 pretrigger model stresses');
Pre5name=strcat(Pre5path,'/',Pre5file);

[Pre60file,Pre60path]=uigetfile('*.mat',...
    'select bdata file containing 60 pretrigger model stresses');
Pre60name=strcat(Pre60path,'/',Pre60file);

%% load data
fprintf('Loading FE data \n')
SG=load(SGname);

fprintf('Loading 5 pretrigger model data /n')
Pre5=load(Pre5name);

fprintf('Loading 60 Pretigger model data \n')
Pre60=load(Pre60name);

%% plot 
ind4=SG.ind4;
lgd4=SG.lgd4;

figure
plot(SG.avgXY_strain(ind4,:),SG.avgXY_stress(ind4,:))
hold on
plot(Pre60.strainxy(ind4,:),Pre60.StressModel.Avxy(ind4,:),'--')
plot(Pre5.strainxy(ind4,:),Pre5.StressModel.Avxy(ind4,:),'-.')
hold off

legend('FE data','60 Pretrigger Frames','5 Frames','location','northwest')
xlabel('Strain XY')
ylabel('Stress XY (MPa)')

saveas(gcf,'5_60Pretrigger_comparison.fig')
saveas(gcf,'5_60Pretrigger_comparison.png')