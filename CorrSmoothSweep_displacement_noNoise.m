%% This script is written to produce heatmaps that determine the optimal 
    %Displacement cropping and correction kernal to produce the lowest
    %RMSE error in displacement measurement 

%Author: Andrew Matejunas

%Date: 2023/03/08

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize
close all; clear variables; clc

%% Choose directory to save results
SweepDir=uigetdir('','Choose Directory to save results');

%% Load colormaps
[CmapName,CmapPath]=uigetfile('*.mat', ...
    'Choose file containing Hot-Cold color Maps');
CmapFile=strcat(CmapPath,'/',CmapName);
Cmap=load(CmapFile);

%% Choose finite element data
[FEname,FEpath]=uigetfile('*.mat', ...
    'Choose file containing finite element kinematic fields');
FEfile=strcat(FEpath,'/',FEname);
%% Remove some options from original IBII codes

FEValidMode=0;
calcKinFieldsFromDisp=1;
saveGridData=0;

%% INITIALISE: Add path for processing functions
%(this section taken from main_IBIIProcessing_v1_0r
fprintf('Adding pathes for processing functions.\n')
% Add the path for the grid method code and other useful functions:
funcPath = [pwd,'\Functions\'];

% If the default path is not found we should find it
if exist(funcPath,'file') ~= 7
    hWarn = warndlg('Folder for processing functions not found.',...
        'Function folder not found');
    waitfor(hWarn);
    funcPath = uigetdir(pwd,'Locate Processing Function Folder');
end
addpath(funcPath);
addpath([funcPath,'GridMethodToolbox\']);


%% LOAD RAW DATA FILES: Raw .tiff images
    hardCodePath=0;
    fprintf('Loading reference image from the selected test data folder.\n')
    [imageFile,imagePath] = uigetfile({'*.*','All Files'},'Select the first image in the sequence');
   
%% Load Processing Parameter Data Structures
fprintf('Loading processing parameters file.\n')
initPath = imagePath;
initFile = 'processingParameters.mat';
if exist([initPath,initFile],'file') ~= 2
    hWarn = warndlg('Processing parameter file does not exist.','Processing parameters not found');
    waitfor(hWarn);
    [initFile,initPath,~] = uigetfile('*.mat','Locate processing parameter file');
end
% Load the processing parameters from file


load([initPath,initFile])

% Store the FE valid mode parameter
globalOpts.FEValidMode = FEValidMode;
globalOpts.calcKinFieldsFromDisp = calcKinFieldsFromDisp;


%% Define Test Dseignation
TestDeg=char(cell2mat(inputdlg('Input Test Designation')));

%% load the FE data
fprintf('Loading FE data \n')
FE=load(FEfile,'pos','disp','accel','strain','stress','time');

%% Calculate the Raw Displacement Fields 
% Process the raw tiff images using the grid method code developed by
% Grediac et al.
fprintf('\n--------------------------------------------------------------\n')
fprintf('GRID METHOD PROCESSING\n')
fprintf('--------------------------------------------------------------\n')

%--------------------------------------------------------------------------
% GRID IMAGE PROCESSING

% Check if the images have already been processed otherwise process the data
gridDataSavePath = imagePath;
gridDataFile = 'GridMethod_ProcessedData.mat';

%fprintf('Checking for existing processed data file.\n')
processGridImages = true;


% fprintf('Processing images using the grid method toolbox.\n')
% Process the image squence with the grid method toolbox
[grid,pos,dispRaw] = func_gridMethodImageProcessing_AJM(imagePath,...
    imageFile,...
    grid,gridMethodOpts,imageNoise);



%--------------------------------------------------------------------------
% Update Geometry and Number of Frames Based on Displacement Matrix Size
%fprintf('Updating geometry variables based on the size of displacement field matrices.\n')
[specimen,grid] = func_updateSpecGeom(specimen,grid,dispRaw);

% Currently the rotations are unused so remove them to save RAM
%disp = rmfield(disp,'rot');


%--------------------------------------------------------------------------
% Create the time vector based on the number of frames in the disp struct
time.numFrames = size(dispRaw.x,3);
time.vec = 0:time.step:(size(dispRaw.x,3)-1)*time.step;

%Create arrays of x and y vectors
X_vec=pos.x;
Y_vec=pos.y;
pos.lengthX = pos.x(end)+pos.xStep/2;
pos.lengthY = pos.y(end)+pos.yStep/2;
pos.xGridF = padarray(pos.xGrid,[0,0,time.numFrames-1], ...
    'replicate','post');
pos.yGridF = padarray(pos.yGrid,[0,0,time.numFrames-1], ...
    'replicate','post');
pos.x0F = squeeze(padarray(pos.x,[0,0,time.numFrames-1], ...
    'replicate','post'));

%--------------------------------------------------------------------------
% Load the Reference Image and Determine Where the Free Edge is
%fprintf('Obtaining and setting the free edge location.\n')
[freeEdge,specimen,dispRaw] =...
    func_getFreeEdge(globalOpts.hardCodeFreeEdge,...
    imagePath,imageFile,specimen,dispRaw);

%% Interpolate the FE data 
fprintf('Interpolating FE fields to grid coordinates \n')
interpMethod='nearest';
[Int.pos,Int.disp,~,~,~] =...
    func_interpFEt2Grid(FE.pos,pos,FE.disp,FE.accel,FE.strain,FE.stress, ...
    interpMethod);

%% Remove finite Element fields
clear FE

%% Set up displacement conditioning sweep parameters
cropVec=1:20;
corrVec=1:20;

%% Initialize RMSE variables for speed
RMSEx=zeros(length(corrVec),length(cropVec));
RMSEy=zeros(length(corrVec),length(cropVec));

% Create cell array of processing structs for Parfor loop
extrapArray=cell(length(cropVec),length(corrVec));
tempOpts=extrapOpts.disp;
for k=1:length(cropVec)
    for m=1:length(corrVec)
    tempOpts.cropPx1st=cropVec(k);
    tempOpts.cropPx2nd=tempOpts.cropPx1st;
    if corrVec(m)>=tempOpts.cropPx1st
        tempOpts.extrapPx1st=corrVec(m);
        tempOpts.extrapPx2nd=tempOpts.extrapPx1st;
    else
        tempOpts.extrapPx1st=tempOpts.cropPx1st;
        tempOpts.extrapPx2nd=tempOpts.cropPx1st;
    end
    extrapArray{k,m}=tempOpts;
    end
end
%% Perform the parametric sweep
extrapOptsD=extrapOpts.disp;
RawDispX=dispRaw.x;
RawDispY=dispRaw.y;
IntDispX=Int.disp.x;
IntDispY=Int.disp.y;

for k=1:length(cropVec)
    parfor m=1:length(corrVec)
          ExtrapOptsD=cell2mat(extrapArray(k,m));
          [Xdisp,Ydisp,dispRx,dispRy]=...
            func_cropAndExtrapDispFields_v4(pos,RawDispX,RawDispY, ...
            ExtrapOptsD,true);
          %crop back to original size
          Xdisp=func_crop2size(Xdisp,dispRx,dispRy);
          Ydisp=func_crop2size(Ydisp,dispRx,dispRy);
          %calculate RMSE
          RMSEx(m,k)=func_calcFFRMSE(IntDispX,Xdisp);
          RMSEy(m,k)=func_calcFFRMSE(IntDispY,Ydisp);
                    
    end
end

%% Plot results
figure('units','centimeters','InnerPosition',[1,1,18.2,9.1])
subplot(1,2,1)
RMSExLim=[min(RMSEx,[],'all'),median(RMSEx,'all')];
RMSEyLim=[min(RMSEy,[],'all'),median(RMSEy,'all')];
imagesc(cropVec,corrVec,RMSEx)
title('RMSE(u_x)')
xlabel('u_x Crop Kernel (px)')
ylabel('u Correction Kernel (px)')
colormap(Cmap.PBW_map)
cx=colorbar;
cx.Label.String='RMSE(u_x) [m]';
set(gca,'YDir','normal')
set(gca,'FontSize',12)
% clim(RMSExLim)

subplot(1,2,2)
imagesc(cropVec,corrVec,RMSEy)
title('RMSE(u_y)')
xlabel('u Crop Kernel (px)')
colormap(Cmap.PBW_map)
cx=colorbar;
cx.Label.String='RMSE(u_y) [m]';
set(gca,'YDir','normal')
set(gca,'FontSize',12)

saveName=strcat(SweepDir,'/',TestDeg,'_DispCorrRMSE');
saveas(gcf,saveName,'fig')
saveas(gcf,saveName,'png')

%% Plot results
figure('units','centimeters','InnerPosition',[1,1,18.2,9.1])
subplot(1,2,1)
imagesc(cropVec,corrVec,log10(RMSEx))
title('RMSE(u_x)')
xlabel('u Crop Kernel (px)')
ylabel('u Correction Kernel (px)')
colormap(Cmap.PBW_map)
cx=colorbar;
cx.Label.String='log_{10}(RMSE(u_x))';
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(1,2,2)
imagesc(cropVec,corrVec,log10(RMSEy))
title('RMSE(u_y)')
xlabel('u Crop Kernel (px)')
colormap(Cmap.PBW_map)
cx=colorbar;
cx.Label.String='log_{10}(RMSE(u_y))';
set(gca,'YDir','normal')
set(gca,'FontSize',12)

saveNameLog=strcat(SweepDir,'/',TestDeg,'_DispCorrLogRMSE');
saveas(gcf,saveNameLog,'fig')
saveas(gcf,saveNameLog,'png')

%%
save(strcat(saveName,'_Results.mat'))