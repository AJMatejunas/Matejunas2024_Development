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
FEdata=load(strcat(FEpath,'/',FEfile),'strain','disp','accel');
FEstrain=FEdata.strain;
FEdisp=FEdata.disp;
FEaccel=FEdata.accel;

%% Save Test Designation
testdeg=desig;

clear desig;

%% subtract shear strain fields
fprintf('Calculating differences \n')

ShearDiff=Raw.strain.s-Cordata.strain.s;
XXdiff=Raw.strain.x-Cordata.strain.x;

%% Determine plot limits
ShearLim=[min(FEstrain.s,[],'all'),max(FEstrain.s,[],'all')];
XXLim=[min(FEstrain.x,[],'all'),max(FEstrain.x,[],'all')];
XYdiffLim=[min(ShearDiff,[],'all'),max(ShearDiff,[],'all')];
XXdiffLim=[min(XXdiff,[],'all'),max(XXdiff,[],'all')];

%% Subtract displacement fields
XdispDiff=Raw.disp.x-Cordata.disp.x;
YdispDiff=Raw.disp.y-Cordata.disp.y;

%% Subtract acceleration fields
XaccelDiff=Raw.accel.x-Cordata.accel.x;
YaccelDiff=Raw.accel.y-Cordata.accel.y;

%% Determine Plot Limits
XdispLim=[min(FEdisp.x,[],'all'),max(FEdisp.x,[],'all')];
XdispDiffLim=[min(XdispDiff,[],'all'),max(XdispDiff,[],'all')];

YdispLim=[min(FEdisp.y,[],'all'),max(FEdisp.y,[],'all')];
YdispDiffLim=[min(YdispDiff,[],'all'),max(YdispDiff,[],'all')];

XaccelLim=[min(FEaccel.x,[],'all'),max(FEaccel.y,[],'all')];
XaccelDiffLim=[min(XaccelDiff,[],'all'),max(XaccelDiff,[],'all')];

YaccelLim=[min(FEaccel.y,[],'all'),max(FEaccel.y,[],'all')];
YaccelDiffLim=[min(YaccelDiff,[],'all'),max(YaccelDiff,[],'all')];


%% Create X and Y and time vectors
X_vec=Raw.pos.x*10^3;
Y_vec=Raw.pos.y*10^3;

timemic=Raw.time.vec*10^6;

%% plot all strain fields
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
    
  %% plot differences
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

%% Plot Displacement Fields
fprintf('Plotting displacement fileds \n')

figure('units','normalized','outerposition',[0 0 1 1])

for n=1:length(timemic)
   timecount=num2str(timemic(n)); 
   framecount=num2str(n); 
  
   
   YYdispRaw=squeeze(Raw.disp.y(:,:,n));
   YYdispCor=squeeze(Cordata.disp.y(:,:,n));
   YYdispDiff=squeeze(YdispDiff(:,:,n));
   
   XXdispRaw=squeeze(Raw.disp.x(:,:,n));
   XXdispDiff=squeeze(XdispDiff(:,:,n));
   XXdispCor=squeeze(Cordata.disp.x(:,:,n));
  
   
   
   %% Plot Raw displacements
    subplot(3,2,1)
    imagesc(X_vec,Y_vec,XXdispRaw)
    title(strcat('Raw GM displacement_{x} ',timecount,'\mus'))
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='u_x';
    caxis(XdispLim)
    
    subplot(3,2,2)
    imagesc(X_vec,Y_vec,YYdispRaw)
    title(strcat('Raw GM diplacement_y ',timecount,'\mus'))
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='u_y';
    caxis(YdispLim)
    
  %% plot corrected GM displacements
    subplot(3,2,3)
    imagesc(X_vec,Y_vec,XXdispCor)
    title(strcat('Corrected ',timecount,'\mus'))
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='u_x';
    caxis(XdispLim)
    
    subplot(3,2,4)
    imagesc(X_vec,Y_vec,YYdispCor)
    title(strcat('Corrected ',timecount,'\mus'))
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='u_y';
    caxis(YdispLim)
    
  %% plot differences
  subplot(3,2,5)
    imagesc(X_vec,Y_vec,XXdispDiff)
    title(strcat('Raw-Corrected ',timecount,'\mus'))
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='u_x';
    caxis(XdispDiffLim)
    
    subplot(3,2,6)
    imagesc(X_vec,Y_vec,YYdispDiff)
    title(strcat('Raw-Corrected ',timecount,'\mus'))
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='u_y';
    caxis(YdispDiffLim)
    
 saveas(gcf,strcat(testdeg,'_DispDiff_',framecount,'.fig'))
 saveas(gcf,strcat(testdeg,'_DispDiff_',framecount,'.png'))
end


%% Plot Acceleration Fields
fprintf('Plotting acceleration fileds \n')

figure('units','normalized','outerposition',[0 0 1 1])

for n=1:length(timemic)
   timecount=num2str(timemic(n)); 
   framecount=num2str(n); 
  
   
   YYaccelRaw=squeeze(Raw.accel.y(:,:,n));
   YYaccelCor=squeeze(Cordata.accel.y(:,:,n));
   YYaccelDiff=squeeze(YaccelDiff(:,:,n));
   
   XXaccelRaw=squeeze(Raw.accel.x(:,:,n));
   XXaccelDiff=squeeze(XaccelDiff(:,:,n));
   XXaccelCor=squeeze(Cordata.accel.x(:,:,n));
  
   
   
   %% Plot Raw accelerations
    subplot(3,2,1)
    imagesc(X_vec,Y_vec,XXaccelRaw)
    title(strcat('Raw GM acceleration_{x} ',timecount,'\mus'))
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='a_x';
    caxis(XaccelLim)
    
    subplot(3,2,2)
    imagesc(X_vec,Y_vec,YYaccelRaw)
    title(strcat('Raw GM acceleration_y ',timecount,'\mus'))
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='a_y';
    caxis(YaccelLim)
    
  %% plot corrected GM accelaerations
    subplot(3,2,3)
    imagesc(X_vec,Y_vec,XXaccelCor)
    title(strcat('Corrected ',timecount,'\mus'))
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='a_x';
    caxis(XaccelLim)
    
    subplot(3,2,4)
    imagesc(X_vec,Y_vec,YYaccelCor)
    title(strcat('Corrected ',timecount,'\mus'))
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='a_y';
    caxis(YaccelLim)
    
  %% plot differences
  subplot(3,2,5)
    imagesc(X_vec,Y_vec,XXaccelDiff)
    title(strcat('Raw-Corrected ',timecount,'\mus'))
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='a_x';
    caxis(XaccelDiffLim)
    
    subplot(3,2,6)
    imagesc(X_vec,Y_vec,YYaccelDiff)
    title(strcat('Raw-Corrected ',timecount,'\mus'))
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='a_y';
    caxis(YaccelDiffLim)
    
 saveas(gcf,strcat(testdeg,'_AccelDiff_',framecount,'.fig'))
 saveas(gcf,strcat(testdeg,'_AccelDiff_',framecount,'.png'))
end



% %% Plot Shear difference only
% fprintf('Plotting Detailed shear difference \n')
% figure('units','normalized','outerposition',[0 0 1 1])
% for n=1:length(timemic)
%    timecount=num2str(timemic(n)); 
%    framecount=num2str(n); 
%    
%    
%    XYstrainDiff=squeeze(ShearDiff(:,:,n));
%     
%    %plot shear
%     imagesc(X_vec,Y_vec,XYstrainDiff)
%     title(strcat('Raw-Corrected shear strain',timecount,'\mus'))
%     xlabel('X coordinate (mm)')
%     ylabel('Y coordinate (mm)')
%     colormap('default')
%     cx=colorbar;
%     cx.Label.String='shear strain';
%     caxis(XYdiffLim)
%  
%     saveas(gcf,strcat(testdeg,'_ShearDiff_',framecount,'.fig'))
%     saveas(gcf,strcat(testdeg,'_ShearDiff_',framecount,'.png'))
% end
% 
% %% Plot abssolute value Shear difference only
% fprintf('Plotting absolute Detailed shear difference \n')
% figure('units','normalized','outerposition',[0 0 1 1])
% for n=1:length(timemic)
%    timecount=num2str(timemic(n)); 
%    framecount=num2str(n); 
%    
%    
%    absXYstrainDiff=abs(squeeze(ShearDiff(:,:,n)));
%     
%    %plot shear
%     imagesc(X_vec,Y_vec,absXYstrainDiff)
%     title(strcat('abs(Raw-Corrected) shear strain',timecount,'\mus'))
%     xlabel('X coordinate (mm)')
%     ylabel('Y coordinate (mm)')
%     colormap('default')
%     cx=colorbar;
%     cx.Label.String='shear strain';
%     caxis([0,XYdiffLim(2)])
%  
%     saveas(gcf,strcat(testdeg,'_AbsShearDiff_',framecount,'.fig'))
%     saveas(gcf,strcat(testdeg,'_AbsShearDiff_',framecount,'.png'))
% end
% 
% %% Plot zoomed abssolute value Shear difference only
% fprintf('Plotting zoomed absolute Detailed shear difference \n')
% figure('units','normalized','outerposition',[0 0 1 1])
% for n=1:length(timemic)
%    timecount=num2str(timemic(n)); 
%    framecount=num2str(n); 
%    
%    
%    absXYstrainDiff=abs(squeeze(ShearDiff(:,:,n)));
%     
%    %plot shear
%     imagesc(X_vec,Y_vec,absXYstrainDiff)
%     title(strcat('abs(Raw-Corrected) shear strain',timecount,'\mus'))
%     xlabel('X coordinate (mm)')
%     ylabel('Y coordinate (mm)')
%     colormap('default')
%     cx=colorbar;
%     cx.Label.String='shear strain';
%     caxis([0,XYdiffLim(2)/5])
%  
%     saveas(gcf,strcat(testdeg,'_zoomShearDiff_',framecount,'.fig'))
%     saveas(gcf,strcat(testdeg,'_zoomShearDiff_',framecount,'.png'))
% end
% %% Plot XX difference only
% fprintf('Plotting Detailed XXstrain difference \n')
% figure('units','normalized','outerposition',[0 0 1 1])
% 
% for n=1:length(timemic)
%    timecount=num2str(timemic(n)); 
%    framecount=num2str(n); 
%    
%    XXstrainDiff=squeeze(XXdiff(:,:,n));
%    
%       
%    %plot XX strain difference
%     imagesc(X_vec,Y_vec,XXstrainDiff)
%     title(strcat('Raw-Corrected XXstrain',timecount,'\mus'))
%     xlabel('X coordinate (mm)')
%     ylabel('Y coordinate (mm)')
%     colormap('default')
%     cx=colorbar;
%     cx.Label.String='shear strain';
%     caxis(XXdiffLim)
%  
%     saveas(gcf,strcat(testdeg,'_XXDiff_',framecount,'.fig'))
%     saveas(gcf,strcat(testdeg,'_XXDiff_',framecount,'.png'))
% end

fprintf('Done \n')