%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This code is written to extract and plot the stress and strain data from 
    %QS baseline tension experiments for viscoelastic IBII experiments

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize
clear variables; close all; clc

%% Load Representative stress strain data
  
%Load Frame
[dataName,dataPath]=uigetfile('*.xlsx', ...
    'Choose excel database containing stress and strain information');
dataFile=strcat(dataPath,'/',dataName);
Rep.LFdata=xlsread(dataFile,'S04_analysis','E4:G895');

Rep.LFstrain=Rep.LFdata(:,1);
Rep.LFstress=Rep.LFdata(:,2);
Rep.LFstrainS=Rep.LFdata(:,3);

%DIC Data
Rep.DICdata=xlsread(dataFile,'S04_analysis','H4:L325');

Rep.DICstrainY=Rep.DICdata(:,1);
Rep.DICstrainX=Rep.DICdata(:,3);
Rep.DICstress=Rep.DICdata(:,4);

%% Generate Plot

figure('Units','inches','InnerPosition',[1,1,5,5])
scatter(Rep.LFstrain,Rep.LFstress*10^-6)
hold on
scatter(Rep.LFstrainS,Rep.LFstress*10^-6)
scatter(Rep.DICstrainY,Rep.DICstress*10^-6)
hold off
xlabel('strain_{yy}')
ylabel('\sigma_{yy} (MPa)')
title('Representative Stress Strain Curve')
legend('Raw LF','Toe Corrected','DIC')
ylim([0,70])
xlim([0,0.025])

%% Generate Representative Plot for Poisson's Ration
nu=0.33;
figure('Units','inches','InnerPosition',[1,1,5,5])
scatter(Rep.DICstrainY,Rep.DICstrainX)
hold on
plot(Rep.DICstrainY,-nu*Rep.DICstrainY)
hold off
xlabel('strain_{yy}')
ylabel('strain_{xx}')
title('Match ID \nu Calculation')
legend('Raw LF','-\nu*strain_{yy}')
