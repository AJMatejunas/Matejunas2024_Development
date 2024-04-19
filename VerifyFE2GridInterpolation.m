% This script is written to verify the interpolation of the kinematic
    %fields from the Finite element coordinate system to the grid method
    %coordinate system

 %Author: Andrew Matejunas

 %% Initialize
 clear variables; close all; clc

 %% Test designation
SimDesig=char(cell2mat(inputdlg('Test Designation')));

 %% load the finite element data
[FEname,FEpath]=uigetfile('*.mat', ...
    'Choose file containing finite element kinematic fields');
FEfile=strcat(FEpath,'/',FEname);
FE=load(FEfile,'pos','disp','accel','strain','stress','time');

 %% Load grid method pos file
[PosName,PosPath]=uigetfile('*.mat', ...
    'choose file containing grid method coordinates');
PosFile=strcat(PosPath,'/',PosName); 
GridPos=load(PosFile,'pos');


%% Perform the interpolation
[pos,disp,accel,strain,stress]=func_interpFEt2Grid(FE.pos,GridPos.pos,...
    FE.disp,FE.accel,FE.strain,FE.stress);

%% Choose directory to save results
SaveDir=uigetdir('','Choose folder in which to save results');

%% Save the kinematic fields
save(strcat(SaveDir,'/',SimDesig,'_interpolatedFields.mat'));

%% Plot the displacement results
dispDir=strcat(SaveDir,'/Disp');
figDir=strcat(dispDir,'/Fig');
pngDir=strcat(dispDir,'/Png');
mkdir(figDir);
mkdir(pngDir);
dispLim.x=[min(FE.disp.x,[],'all'),max(FE.disp.x,[],'all')];
dispLim.y=[min(FE.disp.y,[],'all'),max(FE.disp.y,[],'all')];
figure('units','normalized','OuterPosition',[0,0,1,1])
for t=1:length(FE.time.vec)
    frameNum =num2str(t);
    dispRawX=squeeze(FE.disp.x(:,:,t));
    dispIntX=squeeze(disp.x(:,:,t));
    dispRawY=squeeze(FE.disp.x(:,:,t));
    dispIntY=squeeze(disp.x(:,:,t));

    %% plot raw displacements in the x direction
    subplot(2,2,1)
    imagesc(FE.pos.x,FE.pos.y,dispRawX);
    xlabel('X coordinate')
    ylabel('Y coordinate')
    title(strcat('u_x FE coordinates frame~',frameNum))
    colorbar
    caxis(dispLim.x)

    %% plot raw displacements in the y direction
    subplot(2,2,2)
    imagesc(FE.pos.x,FE.pos.y,dispRawY);
    xlabel('X coordinate')
    ylabel('Y coordinate')
    title(strcat('u_y FE coordinates frame~',frameNum))
    colorbar
    caxis(dispLim.y)

    %% plot interpolated displacements in the x direction
    subplot(2,2,3)
    imagesc(pos.x,pos.y,dispIntX);
    xlabel('X coordinate')
    ylabel('Y coordinate')
    title(strcat('u_x Interpolated frame~',frameNum))
    colorbar
    caxis(dispLim.x)

    %% plot Interpolated displacements in the y direction
    subplot(2,2,4)
    imagesc(pos.x,pos.y,dispIntY);
    xlabel('X coordinate')
    ylabel('Y coordinate')
    title(strcat('u_y Interpolated frame~',frameNum))
    colorbar
    caxis(dispLim.y)

    %% save figure
    pngName=strcat(pngDir,'/',SimDesig,'InterpDisp_Frame',frameNum);
    figName=strcat(figDir,'/',SimDesig,'InterpDisp_Frame',frameNum);
    saveas(gcf,pngName,'png')
    saveas(gcf,figName,'fig')

end

%% Plot the acceleration results
accelDir=strcat(SaveDir,'/accel');
figDir=strcat(accelDir,'/Fig');
pngDir=strcat(accelDir,'/Png');
mkdir(figDir);
mkdir(pngDir);
accelLim.x=[min(FE.accel.x,[],'all'),max(FE.accel.x,[],'all')];
accelLim.y=[min(FE.accel.y,[],'all'),max(FE.accel.y,[],'all')];
figure('units','normalized','OuterPosition',[0,0,1,1])
for t=1:length(FE.time.vec)
    frameNum =num2str(t);
    accelRawX=squeeze(FE.accel.x(:,:,t));
    accelIntX=squeeze(accel.x(:,:,t));
    accelRawY=squeeze(FE.accel.x(:,:,t));
    accelIntY=squeeze(accel.x(:,:,t));

    %% plot raw accelerations in the x direction
    subplot(2,2,1)
    imagesc(FE.pos.x,FE.pos.y,accelRawX);
    xlabel('X coordinate')
    ylabel('Y coordinate')
    title(strcat('a_x FE coordinates frame~',frameNum))
    colorbar
    caxis(accelLim.x)

    %% plot raw accelerations in the y direction
    subplot(2,2,2)
    imagesc(FE.pos.x,FE.pos.y,accelRawY);
    xlabel('X coordinate')
    ylabel('Y coordinate')
    title(strcat('a_y FE coordinates frame~',frameNum))
    colorbar
    caxis(accelLim.y)

    %% plot interpolated accelerations in the x direction
    subplot(2,2,3)
    imagesc(pos.x,pos.y,accelIntX);
    xlabel('X coordinate')
    ylabel('Y coordinate')
    title(strcat('a_x Interpolated frame~',frameNum))
    colorbar
    caxis(accelLim.x)

    %% plot Interpolated accelerations in the y direction
    subplot(2,2,4)
    imagesc(pos.x,pos.y,accelIntY);
    xlabel('X coordinate')
    ylabel('Y coordinate')
    title(strcat('a_y Interpolated frame~',frameNum))
    colorbar
    caxis(accelLim.y)

    %% save figure
    pngName=strcat(pngDir,'/',SimDesig,'InterpAccel_Frame',frameNum);
    figName=strcat(figDir,'/',SimDesig,'InterpAccel_Frame',frameNum);
    saveas(gcf,pngName,'png')
    saveas(gcf,figName,'fig')

end

%% Plot the strain results
strainDir=strcat(SaveDir,'/strain');
figDir=strcat(strainDir,'/Fig');
pngDir=strcat(strainDir,'/Png');
mkdir(figDir);
mkdir(pngDir);
strainLim.x=[min(FE.strain.x,[],'all'),max(FE.strain.x,[],'all')];
strainLim.y=[min(FE.strain.y,[],'all'),max(FE.strain.y,[],'all')];
strainLim.s=[min(FE.strain.s,[],'all'),max(FE.strain.s,[],'all')];
figure('units','normalized','OuterPosition',[0,0,1,1])
for t=1:length(FE.time.vec)
    frameNum =num2str(t);
    strainRawX=squeeze(FE.strain.x(:,:,t));
    strainIntX=squeeze(strain.x(:,:,t));
    strainRawY=squeeze(FE.strain.x(:,:,t));
    strainIntY=squeeze(strain.x(:,:,t));
    strainRawS=squeeze(FE.strain.s(:,:,t));
    strainIntS=squeeze(strain.s(:,:,t));
    %% plot raw strain in the x direction
    subplot(2,3,1)
    imagesc(FE.pos.x,FE.pos.y,strainRawX);
    xlabel('X coordinate')
    ylabel('Y coordinate')
    title(strcat('strain_{xx} FE coordinates frame~',frameNum))
    colorbar
    caxis(strainLim.x)

    %% plot raw strain in the y direction
    subplot(2,3,2)
    imagesc(FE.pos.x,FE.pos.y,strainRawY);
    xlabel('X coordinate')
    ylabel('Y coordinate')
    title(strcat('strain_{yy} FE coordinates frame~',frameNum))
    colorbar
    caxis(strainLim.y)

    %% plot raw strain in shear
    subplot(2,3,3)
    imagesc(FE.pos.x,FE.pos.y,strainRawS);
    xlabel('X coordinate')
    ylabel('Y coordinate')
    title(strcat('strain_{xy} FE coordinates frame~',frameNum))
    colorbar
    caxis(strainLim.s)

    %% plot interpolated strain in the x direction
    subplot(2,3,4)
    imagesc(pos.x,pos.y,strainIntX);
    xlabel('X coordinate')
    ylabel('Y coordinate')
    title(strcat('strain_{xx} Interpolated frame~',frameNum))
    colorbar
    caxis(strainLim.x)

    %% plot Interpolated strains in the y direction
    subplot(2,3,5)
    imagesc(pos.x,pos.y,strainIntY);
    xlabel('X coordinate')
    ylabel('Y coordinate')
    title(strcat('strain_{yy} Interpolated frame~',frameNum))
    colorbar
    caxis(strainLim.y)

    %% plot Interpolated strains in shear
    subplot(2,3,6)
    imagesc(pos.x,pos.y,strainIntS);
    xlabel('X coordinate')
    ylabel('Y coordinate')
    title(strcat('strain_{xy} Interpolated frame~',frameNum))
    colorbar
    caxis(strainLim.s)

    %% save figure
    pngName=strcat(pngDir,'/',SimDesig,'Interpstrain_Frame',frameNum);
    figName=strcat(figDir,'/',SimDesig,'Interpstrain_Frame',frameNum);
    saveas(gcf,pngName,'png')
    saveas(gcf,figName,'fig')
end

%% Plot the stress results
stressDir=strcat(SaveDir,'/stress');
figDir=strcat(stressDir,'/Fig');
pngDir=strcat(stressDir,'/Png');
mkdir(figDir);
mkdir(pngDir);
stressLim.x=[min(FE.stress.x,[],'all'),max(FE.stress.x,[],'all')];
stressLim.y=[min(FE.stress.y,[],'all'),max(FE.stress.y,[],'all')];
stressLim.s=[min(FE.stress.s,[],'all'),max(FE.stress.s,[],'all')];
figure('units','normalized','OuterPosition',[0,0,1,1])
for t=1:length(FE.time.vec)
    frameNum =num2str(t);
    stressRawX=squeeze(FE.stress.x(:,:,t));
    stressIntX=squeeze(stress.x(:,:,t));
    stressRawY=squeeze(FE.stress.x(:,:,t));
    stressIntY=squeeze(stress.x(:,:,t));
    stressRawS=squeeze(FE.stress.s(:,:,t));
    stressIntS=squeeze(stress.s(:,:,t));
    %% plot raw stress in the x direction
    subplot(2,3,1)
    imagesc(FE.pos.x,FE.pos.y,stressRawX);
    xlabel('X coordinate')
    ylabel('Y coordinate')
    title(strcat('stress_{xx} FE coordinates frame~',frameNum))
    colorbar
    caxis(stressLim.x)

    %% plot raw stress in the y direction
    subplot(2,3,2)
    imagesc(FE.pos.x,FE.pos.y,stressRawY);
    xlabel('X coordinate')
    ylabel('Y coordinate')
    title(strcat('stress_{yy} FE coordinates frame~',frameNum))
    colorbar
    caxis(stressLim.y)

    %% plot raw stress in shear
    subplot(2,3,3)
    imagesc(FE.pos.x,FE.pos.y,stressRawS);
    xlabel('X coordinate')
    ylabel('Y coordinate')
    title(strcat('stress_{xy} FE coordinates frame~',frameNum))
    colorbar
    caxis(stressLim.s)

    %% plot interpolated stress in the x direction
    subplot(2,3,4)
    imagesc(pos.x,pos.y,stressIntX);
    xlabel('X coordinate')
    ylabel('Y coordinate')
    title(strcat('stress_{xx} Interpolated frame~',frameNum))
    colorbar
    caxis(stressLim.x)

    %% plot Interpolated stresss in the y direction
    subplot(2,3,5)
    imagesc(pos.x,pos.y,stressIntY);
    xlabel('X coordinate')
    ylabel('Y coordinate')
    title(strcat('stress_{yy} Interpolated frame~',frameNum))
    colorbar
    caxis(stressLim.y)

    %% plot Interpolated stresss in shear
    subplot(2,3,6)
    imagesc(pos.x,pos.y,stressIntS);
    xlabel('X coordinate')
    ylabel('Y coordinate')
    title(strcat('stress_{xy} Interpolated frame~',frameNum))
    colorbar
    caxis(stressLim.s)

    %% save figure
    pngName=strcat(pngDir,'/',SimDesig,'Interpstress_Frame',frameNum);
    figName=strcat(figDir,'/',SimDesig,'Interpstress_Frame',frameNum);
    saveas(gcf,pngName,'png')
    saveas(gcf,figName,'fig')
end