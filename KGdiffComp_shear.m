% This script is written to compare matlab results for different shear
% stress-strain curves for different dynamic shear moduli

%% Initialize
clear all, close all, clc

%% Choose the SG data files for getting FE stresses

%k=g
[KGS.file,KGS.path]=uigetfile('*.mat','Choose k=g SGdata file');
KGS.name=strcat(KGS.path,'/',KGS.file);

%g=.9k
[G09.file,G09.path]=uigetfile('*.mat','Choose g=.9k SGdata file');
G09.name=strcat(G09.path,'/',G09.file);
%g=.8k
[G08.file,G08.path]=uigetfile('*.mat','Choose g=.8k SGdata file');
G08.name=strcat(G08.path,'/',G08.file);
%g=.75k
[G075.file,G075.path]=uigetfile('*.mat','Choose g=.75k SGdata file');
G075.name=strcat(G075.path,'/',G075.file);

%% Choose Constitutive Model Data Files
%k=g
[KGS.Constfile,KGS.Constpath]=uigetfile('*.mat',...
    'Choose k=g Constdata file');
    KGS.Constname=strcat(KGS.Constpath,'/',KGS.Constfile);
%g=.9k
[G09.Constfile,G09.Constpath]=uigetfile('*.mat',...
    'Choose k=g Constdata file');
    G09.Constname=strcat(G09.Constpath,'/',G09.Constfile);
%g=.8k
[G08.Constfile,G08.Constpath]=uigetfile('*.mat',...
    'Choose k=g Constdata file');
    G08.Constname=strcat(G08.Constpath,'/',G08.Constfile);
%g=.75k
[G075.Constfile,G075.Constpath]=uigetfile('*.mat',...
    'Choose k=g Constdata file');
    G075.Constname=strcat(G075.Constpath,'/',G075.Constfile);

%% load the  data
fprintf('Loading constant nu data \n')
KGS.FE=load(KGS.name,'X_vec','avgXY_strain','avgXY_stress',...
    'ind1','ind2','ind3','ind4','ind5','ind6',...
    'lgd1','lgd2','lgd3','lgd4','lgd5','lgd6');
KGS.Const=load(KGS.Constname,'StressModel');    

fprintf('Loading g0.9 data \n')
G09.FE=load(G09.name,'X_vec','avgXY_strain','avgXY_stress',...
    'ind1','ind2','ind3','ind4','ind5','ind6',...
    'lgd1','lgd2','lgd3','lgd4','lgd5','lgd6');
G09.Const=load(G09.Constname,'StressModel');    

fprintf('Loading g0.8 data \n')
G08.FE=load(G08.name,'X_vec','avgXY_strain','avgXY_stress',...
    'ind1','ind2','ind3','ind4','ind5','ind6',...
    'lgd1','lgd2','lgd3','lgd4','lgd5','lgd6');
G08.Const=load(G08.Constname,'StressModel');    

fprintf('Loading g0.75 data \n')
G075.FE=load(G075.name,'X_vec','avgXY_strain','avgXY_stress',...
    'ind1','ind2','ind3','ind4','ind5','ind6',...
    'lgd1','lgd2','lgd3','lgd4','lgd5','lgd6');
G075.Const=load(G075.Constname,'StressModel');    

fprintf('Data Loading Complete \n')


%% plot FE comparison at a specific slice
figure(1)
plot(KGS.FE.avgXY_strain(KGS.FE.ind4,:),KGS.FE.avgXY_stress(KGS.FE.ind4,:))
hold on
plot(G09.FE.avgXY_strain(G09.FE.ind4,:),G09.FE.avgXY_stress(G09.FE.ind4,:),'-')
plot(G08.FE.avgXY_strain(G08.FE.ind4,:),G08.FE.avgXY_stress(G08.FE.ind4,:),'-.')
plot(G075.FE.avgXY_strain(G075.FE.ind4,:),G075.FE.avgXY_stress(G075.FE.ind4,:),':')
hold off

title(strcat('Comparison of Shear stress at',KGS.FE.lgd4,...
    'm from impact edge'))
xlabel('shear strain')
ylabel('shear stress (Pa)')
legend('k=g','g=0.9k','g=0.8k','g=0.75k','location','northwest')

saveas(gcf,'KGdiff_ShearCompFE.fig')
saveas(gcf,'KGdiff_ShearCompFE.png')

%% Plot Constitutive Comparison
figure(1)
plot(KGS.FE.avgXY_strain(KGS.FE.ind4,:),KGS.Const.StressModel.Avxy(KGS.FE.ind4,:))
hold on
plot(G09.FE.avgXY_strain(G09.FE.ind4,:),G09.Const.StressModel.Avxy(G09.FE.ind4,:),'-')
plot(G08.FE.avgXY_strain(G08.FE.ind4,:),G08.Const.StressModel.Avxy(G08.FE.ind4,:),'-.')
plot(G075.FE.avgXY_strain(G075.FE.ind4,:),G075.Const.StressModel.Avxy(G075.FE.ind4,:),':')
hold off

title(strcat('Comparison of Model Shear stress at',KGS.FE.lgd4,...
    'm from impact edge'))
xlabel('shear strain')
ylabel('shear stress (Pa)')
legend('k=g','g=0.9k','g=0.8k','g=0.75k','location','northwest')

saveas(gcf,'KGdiff_ShearCompConst.fig')
saveas(gcf,'KGdiff_ShearCompConst.png')