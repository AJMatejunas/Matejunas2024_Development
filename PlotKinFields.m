% This code is written to plot stress strain and strain rates for processed
    % grid method data
    
    clear all; close all; clc
 %% Choose data and test designation
 GMdatafile=uigetfile('.mat','choose GM data file');
 FEdata=uigetfile('.mat','choose FE data file');
 desig=char(inputdlg('input test designation'));
 
 %% Load data file and define test deg
 GM=load(GMdatafile);
 FE=load(FEdata);
 
testdeg=desig;
clear desig

%% Create X and Y vectors
X_vecGM=GM.pos.x*10^3;
Y_vecGM=GM.pos.y'*10^3;

%% Generate limits 
XXStressLim=[min(FE.stress.x,[],'all'),max(FE.stress.x,[],'all')]*10^-6;
XYStressLim=[min(FE.stress.s,[],'all'),max(FE.stress.s,[],'all')]*10^-6;

XXstrainLim=[min(FE.strain.x,[],'all'),max(FE.strain.x,[],'all')];
XYstrainLim=[min(FE.strain.s,[],'all'),max(FE.strain.s,[],'all')];

% XXRateLim=[min(GM.strainRate.x,[],'all'),max(GM.strainRate.x,[],'all')];
% XYRateLim=[min(GM.strainRate.s,[],'all'),max(GM.strainRate.s,[],'all')];

XXRateLim=[-2*10^3,2*10^3];
XYRateLim=[-2*10^3,2*10^3];

%% generate plots
timemic=GM.time.vec*10^6; %us

Xcoord=FE.pos.x*10^3;
Ycoord=FE.pos.y'*10^3;

figure('units','normalized','outerposition',[0 0 1 1])
for n=1:length(timemic)
   timecount=num2str(timemic(n)); 
   framecount=num2str(n); 
   
   Xstress=squeeze(FE.stress.x(:,:,n))/10^6;
   XYstress=squeeze(FE.stress.s(:,:,n))/10^6;
   
   Xstrain=squeeze(FE.strain.x(:,:,n));
   XYstrain=squeeze(FE.strain.s(:,:,n));
   
   Xrate=squeeze(GM.strainRate.x(:,:,n));
   XYrate=squeeze(GM.strainRate.s(:,:,n));
   
   
   %% Plot stresses
    subplot(3,2,1)
    imagesc(Xcoord,Ycoord,Xstress)
    title(strcat('$\sigma{}_{xx} $ ',timecount,' $ \mu{}s$'),'interpreter','latex')
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='$\sigma_{xx}$(MPa)';
    cx.Label.Interpreter='latex';
    caxis(XXStressLim)
    
    subplot(3,2,2)
    imagesc(Xcoord,Ycoord,XYstress)
    title(strcat('$\sigma{}_{xy}$ ',timecount,' $ \mu{}s$'),'interpreter','latex')
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='$\sigma_{xy}$(MPa)';
    cx.Label.Interpreter='latex';
    caxis(XYStressLim)
    
  %% plot strains
    subplot(3,2,3)
    imagesc(Xcoord,Ycoord,Xstrain)
    title(strcat('$\varepsilon{}_{xx} $ ',timecount,' $\mu{}s$'),'interpreter','latex')
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='$\varepsilon{}_{xx}$';
    cx.Label.Interpreter='latex';
    caxis(XXstrainLim)
    
    subplot(3,2,4)
    imagesc(Xcoord,Ycoord,XYstrain)
    title(strcat('$\varepsilon{}_{xy} $ ',timecount,' $\mu{}s$'),'interpreter','latex')
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='$\varepsilon{}_{xy}$';
    cx.Label.Interpreter='latex';
    caxis(XYstrainLim)
    
  %% plot rates
  subplot(3,2,5)
    imagesc(X_vecGM,Y_vecGM,Xrate)
    title(strcat('$\dot{\varepsilon{}}_{xx} $ ',timecount,' $\mu{}s$'),...
        'interpreter','latex')
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='$\dot{\varepsilon{}}_{xx} (s^{-1})$';
    cx.Label.Interpreter='latex';
    caxis(XXRateLim)
    
    subplot(3,2,6)
    imagesc(X_vecGM,Y_vecGM,XYrate)
    title(strcat('$\dot{\varepsilon{}}_{xy} $ ',timecount,' $\mu{}s$'),...
        'interpreter','latex')
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='$\dot{\varepsilon{}}_{xx} (s^{-1})$';
    cx.Label.Interpreter='latex';
    caxis(XYRateLim)
    
a=findall(gcf,'Type','axes');
h=findall(gcf,'Type','line');
t=findall(gcf,'Type','text');

set(t, ...
    'FontName','Helvetica', ...
    'FontSize',18, ...
    'FontWeight','bold' );

set(a                ,               ...
    'FontName'       , 'Helvetica' , ...
    'FontSize'       , 18     , ...
    'FontWeight'     , 'normal');

    
 saveas(gcf,strcat(testdeg,'_Fields_',framecount,'.fig'))
 saveas(gcf,strcat(testdeg,'_Fields_',framecount,'.png'))
end
