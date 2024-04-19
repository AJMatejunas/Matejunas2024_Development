% This script is written to compare the stress strain curves at different
% locations for two different numbers of empty frames

%% initialize
clear all, close all, clc

%% find data file
[P5file,P5path]=uigetfile('*.mat',...
    'Select 5 pretrigger frame SG data file');
P5name=strcat(P5path,'/',P5file);
[P60file,P60path]=uigetfile('*.mat',...
    'Select 60 pretrigger frame SG data file');
P60name=strcat(P60path,'/',P60file);

%% load data
fprintf('Loading 5 pretrigger frame data \n')
Pre5=load(P5name);

fprintf('Loading 60 pretrigger frame data \n')
Pre60=load(P60name);

fprintf('Loading complete \n')

%% Plot two locations

ind3=Pre5.ind3;
ind4=Pre5.ind4;

figure()
plot(Pre5.avgXY_strain(ind3,:),Pre5.Shear_SG(ind3,:));
hold on
plot(Pre5.avgXY_strain(ind4,:),Pre5.Shear_SG(ind4,:));
plot(Pre60.avgXY_strain(ind3,:),Pre60.Shear_SG(ind3,:),'--');
plot(Pre60.avgXY_strain(ind4,:),Pre60.Shear_SG(ind4,:),'--');

hold off
xlabel('average strain_{xy}')
ylabel('\sigma_{xy}')

legend('5Pre 27mm','5Pre 42mm',...
    '60pre 27mm','60Pre 42mm',...
    'location','southeast');

saveas(gcf,'Pretrigger_SG_Comp.fig');
saveas(gcf,'Pretrigger_SG_Comp.png');
