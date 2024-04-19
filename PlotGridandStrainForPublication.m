% This code is written to plot stress strain and strain krates for processed
    % grid method data
    
    clear all; close all; clc
 %% Choose data and test designation
GMdatafile=uigetfile('.mat','choose GM data file');
FEdatafile=uigetfile('.mat','Choose Finite Element Data File');
[GridImageFile,GridImagePath]=uigetfile('.tiff','Choose Grid Image With Scale Bar');
desig=char(inputdlg('input test designation'));
 
 %% Load data file and define test deg
GM=load(GMdatafile);
FE=load(FEdatafile);
GridImage=strcat(GridImagePath,'/',GridImageFile);
 
testdeg=desig;
clear desig

%% Create X and Y vectors
X_vecGM=GM.pos.x*10^3;
Y_vecGM=GM.pos.y'*10^3;

%% Generate limits 

XXstrainLim=[min(FE.strain.x,[],'all'),max(FE.strain.x,[],'all')];
YYstrainLim=[min(FE.strain.y,[],'all'),max(FE.strain.y,[],'all')];
XYstrainLim=[min(FE.strain.s,[],'all'),max(FE.strain.s,[],'all')];


%% generate plots
timemic=GM.time.vec*10^6; %us

Xcoord=GM.pos.x*10^3;
Ycoord=GM.pos.y'*10^3;

figure('units','Centimeters','outerposition',[0 0 27.4 17.2])
for n=1:length(timemic)
   timecount=num2str(timemic(n)); 
   framecount=num2str(n); 
      
   Xstrain=squeeze(GM.strain.x(:,:,n));
   XYstrain=squeeze(GM.strain.s(:,:,n));
   Ystrain=squeeze(GM.strain.y(:,:,n));
   
   Xrate=squeeze(GM.strainRate.x(:,:,n));
   XYrate=squeeze(GM.strainRate.s(:,:,n));
   
   
   %% Plot Synthetic Image
    subplot(2,2,1)
    imshow(GridImage)
    title('Synthetic Grid','interpreter','latex')
   
    
  %% plot Shear strain
    subplot(2,2,2)
    imagesc(Xcoord,Ycoord,XYstrain)
    title(strcat('$$\varepsilon_{xy}$$~',timecount,'$$\mu{}s$$'),'interpreter','latex')
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='$\varepsilon$';
    cx.Label.Interpreter='latex';
    caxis(XYstrainLim)
   
  %% Plot XX strain
    subplot(2,2,3)
    imagesc(Xcoord,Ycoord,Xstrain)
    title(strcat('$$\varepsilon_{xx}$$~',timecount,'$$\mu{}s$$'),'interpreter','latex')
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='$\varepsilon$';
    cx.Label.Interpreter='latex';
    caxis(XXstrainLim)
    
 %% Plot YY strain
    subplot(2,2,4)
    imagesc(Xcoord,Ycoord,Ystrain)
    title(strcat('$$\varepsilon_{yy}$$~',timecount,'$$\mu{}s$$'),'interpreter','latex')
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='$\varepsilon$';
    cx.Label.Interpreter='latex';
    caxis(YYstrainLim)  
  
a=findall(gcf,'Type','axes');
h=findall(gcf,'Type','line');
t=findall(gcf,'Type','text');

set(t, ...
    'FontName','Helvetica', ...
    'FontSize',12, ...
    'FontWeight','bold' );

set(a                ,               ...
    'FontName'       , 'Helvetica' , ...
    'FontSize'       , 12     , ...
    'FontWeight'     , 'normal');

    
 saveas(gcf,strcat(testdeg,'_PubStrains_',framecount,'.fig'))
 saveas(gcf,strcat(testdeg,'_PubStrains_',framecount,'.png'))
end
