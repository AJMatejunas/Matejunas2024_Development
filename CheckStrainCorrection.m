% This script is written to check the correction of Grid method strains at
% the egdges of the specimen

%Change log
    %2022-03-07: Added capabilities to deal with files that are outside of
        %the local path 
%% Initialize
clear all, close all, clc

%% Find files to load

desig=char(inputdlg('input test designation'));
[Rawfile,RawPath]=uigetfile('*.mat','Choose Raw GM data');
[Corrfile,CorrPath]=uigetfile('*.mat','Choose Corrected GM data');
[FEfile,FEpath]=uigetfile('*.mat','Choose FE data');

fprintf('Loading Raw data file \n')
Raw=load(strcat(RawPath,'/',Rawfile));

fprintf('Loading corrected data \n')
Cordata=load(strcat(CorrPath,'/',Corrfile));

%% Load FE strains
fprintf('loading FE strains \n')
FEstrain=load(strcat(FEpath,'/',FEfile),'strain');

testdeg=desig;
clear desig;

%% subtract shear strain fields
fprintf('Calculating differences \n')

ShearDiff=Raw.strain.s-Cordata.strain.s;
XXdiff=Raw.strain.x-Cordata.strain.x;

%% Determine plot limits
ShearLim=[min(FEstrain.strain.s,[],'all'),max(FEstrain.strain.s,[],'all')];
XXLim=[min(FEstrain.strain.x,[],'all'),max(FEstrain.strain.x,[],'all')];
XYdiffLim=[min(ShearDiff,[],'all'),max(ShearDiff,[],'all')];
XXdiffLim=[min(XXdiff,[],'all'),max(XXdiff,[],'all')];


%% Create X and Y and time vectors
X_vec=Raw.pos.x*10^3;
Y_vec=Raw.pos.y*10^3;

timemic=Raw.time.vec*10^6;

%% plot all fields
fprintf('Plotting \n')

figure('units','normalized','outerposition',[0 0 1 1])

for n=1:length(timemic)
   timecount=num2str(timemic(n)); 
   framecount=num2str(n); 
  
   
   XYstrainRaw=squeeze(Raw.strain.s(:,:,n));
   XYstrainCor=squeeze(Cordata.strain.s(:,:,n));
   XYstrainDiff=squeeze(ShearDiff(:,:,n));
   
   XXstrainRaw=squeeze(Raw.strain.x(:,:,n));
   XXstrainDiff=squeeze(XXdiff(:,:,n));
   XXstrainCor=squeeze(Cordata.strain.x(:,:,n));
  
   
   
   %% Plot FE strains
    subplot(3,2,1)
    imagesc(X_vec,Y_vec,XXstrainRaw)
    title(strcat('Raw GM strain_{xx} ',timecount,'\mus'))
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='strain';
    caxis(XXLim)
    
    subplot(3,2,2)
    imagesc(X_vec,Y_vec,XYstrainRaw)
    title(strcat('Raw shear Strain ',timecount,'\mus'))
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='shear strain';
    caxis(ShearLim)
    
  %% plot GM strains
    subplot(3,2,3)
    imagesc(X_vec,Y_vec,XXstrainCor)
    title(strcat('Corrected ',timecount,'\mus'))
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='strain';
    caxis(XXLim)
    
    subplot(3,2,4)
    imagesc(X_vec,Y_vec,XYstrainCor)
    title(strcat('Corrected shear Strain ',timecount,'\mus'))
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='shear strain';
    caxis(ShearLim)
    
  %% plot rates
  subplot(3,2,5)
    imagesc(X_vec,Y_vec,XXstrainDiff)
    title(strcat('Raw-Corrected ',timecount,'\mus'))
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='strain';
    caxis(XXdiffLim)
    
    subplot(3,2,6)
    imagesc(X_vec,Y_vec,XYstrainDiff)
    title(strcat('Raw-Corrected ',timecount,'\mus'))
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='shear strain';
    caxis(XYdiffLim)
    
 saveas(gcf,strcat(testdeg,'_CorrDiff_',framecount,'.fig'))
 saveas(gcf,strcat(testdeg,'_CorrDiff_',framecount,'.png'))
end

%% Plot Shear difference only
fprintf('Plotting Detailed shear difference \n')
figure('units','normalized','outerposition',[0 0 1 1])
for n=1:length(timemic)
   timecount=num2str(timemic(n)); 
   framecount=num2str(n); 
   
   
   XYstrainDiff=squeeze(ShearDiff(:,:,n));
    
   %plot shear
    imagesc(X_vec,Y_vec,XYstrainDiff)
    title(strcat('Raw-Corrected shear strain',timecount,'\mus'))
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('default')
    cx=colorbar;
    cx.Label.String='shear strain';
    caxis(XYdiffLim)
 
    saveas(gcf,strcat(testdeg,'_ShearDiff_',framecount,'.fig'))
    saveas(gcf,strcat(testdeg,'_ShearDiff_',framecount,'.png'))
end

%% Plot abssolute value Shear difference only
fprintf('Plotting absolute Detailed shear difference \n')
figure('units','normalized','outerposition',[0 0 1 1])
for n=1:length(timemic)
   timecount=num2str(timemic(n)); 
   framecount=num2str(n); 
   
   
   absXYstrainDiff=abs(squeeze(ShearDiff(:,:,n)));
    
   %plot shear
    imagesc(X_vec,Y_vec,absXYstrainDiff)
    title(strcat('abs(Raw-Corrected) shear strain',timecount,'\mus'))
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('default')
    cx=colorbar;
    cx.Label.String='shear strain';
    caxis([0,XYdiffLim(2)])
 
    saveas(gcf,strcat(testdeg,'_AbsShearDiff_',framecount,'.fig'))
    saveas(gcf,strcat(testdeg,'_AbsShearDiff_',framecount,'.png'))
end

%% Plot zoomed abssolute value Shear difference only
fprintf('Plotting zoomed absolute Detailed shear difference \n')
figure('units','normalized','outerposition',[0 0 1 1])
for n=1:length(timemic)
   timecount=num2str(timemic(n)); 
   framecount=num2str(n); 
   
   
   absXYstrainDiff=abs(squeeze(ShearDiff(:,:,n)));
    
   %plot shear
    imagesc(X_vec,Y_vec,absXYstrainDiff)
    title(strcat('abs(Raw-Corrected) shear strain',timecount,'\mus'))
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('default')
    cx=colorbar;
    cx.Label.String='shear strain';
    caxis([0,XYdiffLim(2)/5])
 
    saveas(gcf,strcat(testdeg,'_zoomShearDiff_',framecount,'.fig'))
    saveas(gcf,strcat(testdeg,'_zoomShearDiff_',framecount,'.png'))
end
%% Plot XX difference only
fprintf('Plotting Detailed XXstrain difference \n')
figure('units','normalized','outerposition',[0 0 1 1])

for n=1:length(timemic)
   timecount=num2str(timemic(n)); 
   framecount=num2str(n); 
   
   XXstrainDiff=squeeze(XXdiff(:,:,n));
   
      
   %plot XX strain difference
    imagesc(X_vec,Y_vec,XXstrainDiff)
    title(strcat('Raw-Corrected XXstrain',timecount,'\mus'))
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('default')
    cx=colorbar;
    cx.Label.String='shear strain';
    caxis(XXdiffLim)
 
    saveas(gcf,strcat(testdeg,'_XXDiff_',framecount,'.fig'))
    saveas(gcf,strcat(testdeg,'_XXDiff_',framecount,'.png'))
end

fprintf('Done \n')