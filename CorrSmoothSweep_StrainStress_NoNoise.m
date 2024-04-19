% This function is written to plot the errors in the kinematic fields as a
    % function of extrapolation kernal and smoothing kernal. Strain is a
    % function spatial smoothing kernal and stress gauge stresses is a
    % funciton of temporal smoothing kernal

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
globalOpts.FEValidMode = false;
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

RawDispX=dispRaw.x;
RawDispY=dispRaw.y;
%% Interpolate the FE data 
fprintf('Interpolating FE fields to grid coordinates \n')
interpMethod='nearest';
[Int.pos,Int.disp,Int.accel,Int.strain,Int.stress] =...
    func_interpFEt2Grid(FE.pos,pos,FE.disp,FE.accel,FE.strain,FE.stress, ...
    interpMethod);
%% Add reference parameters file 
[RefPar.file,RefPar.path]=uigetfile('*.mat',...
    'Choose Mat file containing reference constitutive pamrameters');
load(strcat(RefPar.path,'/',RefPar.file),'MatProps');


%% Calcualte Constitutive model stresses
fprintf('Calculating out of plane strains for FE data \n')
Int.Model=func_ViscoConstitutiveV6(...
    Int.strain.x,Int.strain.y,Int.strain.s,time.vec,MatProps, ...
    0,0,0);
Int.strain.z=Int.Model.zzstrain;
% Remove constitutive model stresses for memory
Int=rmfield(Int,'Model');

%% Calculate Average stresses
Int.stress.xAvg=squeeze(mean(Int.stress.x));
Int.stress.sAvg=squeeze(mean(Int.stress.s));

%% Remove finite Element fields;
clear FE


%% Generate smoothing and extrapolation kernal vectrors
CorrVec=7:41; %Vector of extrapolation kernals
TkernVec=[1,5:2:41]; %Vector of temporal extrapolation kernalsclc
SkernVec=[1,5:2:41]; %Vector of spatial smoothing kernals

%% Generate Array of Spatial Correction interval vs spatial smoothing 
    %Kernal
StrainCorrOpts=cell(length(CorrVec),1);
StrainSmooth=cell(length(SkernVec),1);

DispOpts=cell(length(CorrVec),1);
smoothOpts.strain.spatialSmooth=0;
StrainSmooth{1}=smoothOpts.strain;

for k=2:length(SkernVec)
    tempStrainSmooth=smoothOpts.strain;
    tempStrainSmooth.spatialSmooth=1;
    tempStrainSmooth.spatialKernelSize=SkernVec(k)*[1,1];
    StrainSmooth{k}=tempStrainSmooth;
end
    
 for m=1:length(CorrVec)
    tempStrainOpts=extrapOpts.strain;
    tempDispOpts=extrapOpts.disp;
    tempDispOpts.extrapPx1st=CorrVec(m);
    tempDispOpts.extrapPx2nd=tempDispOpts.extrapPx1st;
    tempStrainOpts.extrapPx1st=tempDispOpts.extrapPx1st;
    tempStrainOpts.extrapPx2nd=tempDispOpts.extrapPx1st;
    tempStrainOpts.enforceGlobBCsPx=CorrVec(m);  
    DispOpts{m}=tempDispOpts;
    StrainCorrOpts{m}=tempStrainOpts;
 end

 clear tempDispOpts tempAccelOpts tempStrainSmooth tempStrainOpts


%% set window of interest
RpitchImp=10;
RpitchFree=10;
RpxImp=RpitchImp*grid.pxPerPeriod;
RpxFree=RpitchFree*grid.pxPerPeriod;
Xind=RpxFree:(length(pos.x)-RpxImp);

RframeStart=5;
RframeEnd=5;
Tind=RframeStart:(length(time.vec)-RframeEnd);

%% Set reference strains
StrainXref=Int.strain.x(:,Xind,Tind);
StrainYref=Int.strain.y(:,Xind,Tind);
StrainSref=Int.strain.s(:,Xind,Tind);
StrainZref=Int.strain.z(:,Xind,Tind);

%% Set reference Stresses
StressXref=Int.stress.xAvg(Xind,Tind);
StressSref=Int.stress.sAvg(Xind,Tind);

%% Loop over processed data
fprintf('Running Strain smooth-corr sweep \n')
RMSEstrainX=zeros(length(CorrVec),length(SkernVec));
RMSEstressX=zeros(length(CorrVec),length(SkernVec));
RMSEstrainY=zeros(length(CorrVec),length(SkernVec));
RMSEstrainS=zeros(length(CorrVec),length(SkernVec));
RMSEstressS=zeros(length(CorrVec),length(SkernVec));
RMSEstrainZ=zeros(length(CorrVec),length(SkernVec));
%% Loop over processing parameters

for m=1:length(CorrVec)
    extrapD=cell2mat(DispOpts(m));
    extrapS=cell2mat(StrainCorrOpts(m));
    parfor k=1:length(SkernVec)
       SmoothS=cell2mat(StrainSmooth(k));
       %Correct Displacements
       tempDisp=[];
       [tempDisp.x,tempDisp.y,tempDisp.rX,tempDisp.rY]=...
           func_cropAndExtrapDispFields_v4(pos,RawDispX,RawDispY, ...
           extrapD,true);
       %Calculate Strains
       [strain,~]=func_smoothCalcStrain_v4(pos,time,tempDisp,...
           SmoothS,extrapS,true);
       %Calculate constituitive models stresses
       tempModel=func_ViscoConstitutiveV6(strain.x,strain.y,strain.s, ...
           time.vec,MatProps,0,0,0);

       %crop Strain to the ROI
       strain.x=strain.x(:,Xind,Tind);
       strain.y=strain.y(:,Xind,Tind);
       strain.s=strain.s(:,Xind,Tind);
       strain.z=tempModel.zzstrain(:,Xind,Tind);
       
       %Crop stress to ROI
       stressX=tempModel.Avxx(Xind,Tind);
       stressS=tempModel.Avxy(Xind,Tind);
       
       %Calculate RMSE for strains
       RMSEstrainX(m,k)=func_calcFFRMSE(StrainXref,strain.x);
       RMSEstrainY(m,k)=func_calcFFRMSE(StrainYref,strain.y);
       RMSEstrainS(m,k)=func_calcFFRMSE(StrainSref,strain.s);
       RMSEstrainZ(m,k)=func_calcFFRMSE(StrainZref,strain.z);
       % Calculate RMSE for Stress
       RMSEstressX(m,k)=func_calcFFRMSE(StressXref,stressX);
       RMSEstressS(m,k)=func_calcFFRMSE(StressSref,stressS);
    end
end

clear strain stressX stressS tempModel tempDisp 
Int=rmfield(Int,{'strain','disp','accel'});

%% Plot Heat Maps of strain RMSE errors
figure('units','centimeters','InnerPosition',[1,1,18.2,18.2])
subplot(2,2,1)
imagesc(SkernVec,CorrVec,RMSEstrainX)
title('RMSE(strain_{xx}')
ylabel('Extrapolation Interval (px)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(2,2,2)
imagesc(SkernVec,CorrVec,RMSEstrainY)
title('RMSE(strain_{yy}')
ylabel('Extrapolation Interval (px)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(2,2,3)
imagesc(SkernVec,CorrVec,RMSEstrainS)
title('RMSE(strain_{xy}')
xlabel('SK (px)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(2,2,4)
imagesc(SkernVec,CorrVec,RMSEstrainZ)
title('RMSE(strain_{zz}')
xlabel('SK (px)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

saveName=strcat(SweepDir,'/',TestDeg,'_SpatialCorrSmoothSweep_strainRMSE');
saveas(gcf,saveName,'fig')
saveas(gcf,saveName,'png')


%% Plot Heat Maps of constitutive Model stress RMSE
figure('units','centimeters','InnerPosition',[1,1,18.2,9.1])
subplot(1,2,1)
imagesc(SkernVec,CorrVec,RMSEstressX)
title('RMSE(\sigma^{model}_{xx}')
xlabel('SK (px)')
ylabel('Extrapolation Interval (px)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(1,2,2)
imagesc(SkernVec,CorrVec,RMSEstressS)
title('RMSE(\sigma^{model}_{xy}')
xlabel('SK (px)')

cx=colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

saveName=strcat(SweepDir,'/',TestDeg,'_SpatialCorrS,ppthSweep_modelRMSE');
saveas(gcf,saveName,'fig')
saveas(gcf,saveName,'png')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up smoothing Parameters for the acceleration sweep
extrapOpts.accel.tempPadMethod='replicate';
smoothOpts.accel.temporalSmooth=1;
AccelCorrOpts=cell(length(TkernVec),length(CorrVec));
smoothOpts.accel.temporalSmooth=0;
AccelSmooth=cell(length(TkernVec),1);
AccelSmooth{1}=smoothOpts.accel;

 for m=1:length(CorrVec)
    tempAccelOpts=extrapOpts.accel;
    tempAccelOpts.extrapPx1st=CorrVec(m);
    tempAccelOpts.extrapPx2nd=tempAccelOpts.extrapPx1st;
    AccelCorrOpts{1,m}=tempAccelOpts;
    for k=2:length(TkernVec)
        tempAccelOpts.tempPadFrames=round(TkernVec(k)/2)*[1,1];
        tempSmoothA=smoothOpts.accel;
        tempSmoothA.temporalSmooth=true;
        tempSmoothA.temporalKernelSize=TkernVec(k)*[1,1];
        AccelSmooth{k}=tempSmoothA;
        AccelCorrOpts{k,m}=tempAccelOpts;
    end
 end

%% Run acceleration sweep
fprintf('Running Acceleration smooth-corr sweep \n')
RMSEsgX=zeros(length(CorrVec),length(SkernVec));
RMSEsgS=zeros(length(CorrVec),length(SkernVec));


rho=MatProps.rho;
for m=1:length(CorrVec)
    extrapD=cell2mat(DispOpts(m));
    parfor k=1:length(TkernVec)
       extrapA=cell2mat(AccelCorrOpts(k,m));
       smoothA=cell2mat(AccelSmooth(k));
       %correct displacements
       tempDisp=[];
       [tempDisp.x,tempDisp.y,tempDisp.Rx,tempDisp.Ry]=...
           func_cropAndExtrapDispFields_v4(pos,RawDispX,RawDispY, ...
           extrapD,true);
       %calaculate Accelerations
       [accel,~,~]=func_smoothCalcAccel_v4(pos,time,tempDisp,...
       smoothA,extrapA,diffOpts,true);
       %Calculate stress gauge stresses
       SG=func_Full_SG(accel,X_vec,time,rho);

       %Crop to ROI
       SGx=SG.x(Xind,Tind);
       SGs=SG.s(Xind,Tind);

       %% Calculate RMSE
       RMSEsgX(m,k)=func_calcFFRMSE(StressXref,SGx);
       RMSEsgS(m,k)=func_calcFFRMSE(StressSref,SGs);
    end
end

%% Plot Heat Maps of stress gauge RMSE
figure('units','centimeters','InnerPosition',[1,1,18.2,9.1])
subplot(1,2,1)
imagesc(TkernVec,CorrVec,RMSEsgX)
title('RMSE(\sigma^{SG}_{xx})')
ylabel('Extrapolation Interval (px)')
xlabel('TK (frames)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(1,2,2)
imagesc(TkernVec,CorrVec,RMSEsgS)
title('RMSE(\sigma^{SG}_{xy})')
xlabel('TK (frames)')
cx=colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

saveName=strcat(SweepDir,'/',TestDeg, ...
    '_SpatialCorr_TemPSmooth_Sweep_SG_RMSE');
saveas(gcf,saveName,'fig')
saveas(gcf,saveName,'png')

%% Calculate Normnlized RMSE
RMSENorm.strainX=RMSEstrainX/(max(abs(StrainXref),[],'all'))*100;
RMSENorm.strainY=RMSEstrainY/(max(abs(StrainYref),[],'all'))*100;
RMSENorm.strainS=RMSEstrainS/(max(abs(StrainSref),[],'all'))*100;
RMSENorm.strainZ=RMSEstrainZ/(max(abs(StrainZref),[],'all'))*100;

RMSENorm.modelX=RMSEstressX/max(abs(StressXref),[],'all')*100;
RMSENorm.modelS=RMSEstressS/max(abs(StressSref),[],'all')*100;

RMSENorm.SGX=RMSEsgX/max(abs(StressXref),[],'all')*100;
RMSENorm.SGS=RMSEsgS/max(abs(StressSref),[],'all')*100;

%% Plot Heat Maps of stress gauge RMSE
figure('units','centimeters','InnerPosition',[1,1,18.2,9.1])
subplot(1,2,1)
imagesc(TkernVec,CorrVec,RMSENorm.SGX)
title('RMSE(\sigma^{SG}_{xx})/max(\sigma_{xx})')
ylabel('Extrapolation Interval (px)')
xlabel('TK (frames)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)
%caxis([0,25])

subplot(1,2,2)
imagesc(TkernVec,CorrVec,RMSENorm.SGS)
title('RMSE(\sigma^{SG}_{xy})/max(\sigma_{xy})')
xlabel('TK (frames)')
cx=colorbar;
cx.Label.String='%';
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)
%caxis([0,25])

saveName=strcat(SweepDir,'/',TestDeg, ...
    '_SpatialCorr_TemPSmooth_Sweep_SG_RMSE_Norm');
saveas(gcf,saveName,'fig')
saveas(gcf,saveName,'png')
%% Plot Heat Maps of stress gauge RMSE
figure('units','centimeters','InnerPosition',[1,1,18.2,9.1])
subplot(1,2,1)
imagesc(TkernVec,CorrVec,RMSENorm.SGX)
title('RMSE(\sigma^{SG}_{xx})/max(\sigma_{xx})')
ylabel('Extrapolation Interval (px)')
xlabel('TK (frames)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)
caxis([0,25])

subplot(1,2,2)
imagesc(TkernVec,CorrVec,RMSENorm.SGS)
title('RMSE(\sigma^{SG}_{xy})/max(\sigma_{xy})')
xlabel('TK (frames)')
cx=colorbar;
cx.Label.String='%';
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)
caxis([0,25])

saveName=strcat(SweepDir,'/',TestDeg, ...
    '_SpatialCorr_TemPSmooth_Sweep_SG_RMSE_NormZoom');
saveas(gcf,saveName,'fig')
saveas(gcf,saveName,'png')

%% Plot Heat Maps of stress RMSE
figure('units','centimeters','InnerPosition',[1,1,18.2,9.1])
subplot(1,2,1)
imagesc(SkernVec,CorrVec,RMSENorm.modelX)
title('RMSE(\sigma^{Model}_{xx})/max(\sigma_{xx})')
ylabel('Extrapolation Interval (px)')
xlabel('SK (pixels)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)
%caxis([0,25])

subplot(1,2,2)
imagesc(SkernVec,CorrVec,RMSENorm.modelS)
title('RMSE(\sigma^{Model}_{xy})/max(\sigma_{xy})')
xlabel('SK (pixels)')
cx=colorbar;
cx.Label.String='%';
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)
%caxis([0,25])

saveName=strcat(SweepDir,'/',TestDeg, ...
    '_SpatialCorr_TemPSmooth_Sweep_Stress_RMSE_Norm');
saveas(gcf,saveName,'fig')
saveas(gcf,saveName,'png')

%% Plot Heat Maps of strain RMSE errors
figure('units','centimeters','InnerPosition',[1,1,18.2,18.2])
subplot(2,2,1)
imagesc(SkernVec,CorrVec,RMSENorm.strainX)
title('RMSE(strain_{xx})/max(strain_{xx})')
ylabel('Extrapolation Interval (px)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(2,2,2)
imagesc(SkernVec,CorrVec,RMSENorm.strainY)
title('RMSE(strain_{yy})/max(strain_{yy})')
cx=colorbar;
cx.Label.String='%';
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(2,2,3)
imagesc(SkernVec,CorrVec,RMSENorm.strainS)
title('RMSE(strain_{xy})/max(strain_{xy})')
xlabel('SK (px)')
ylabel('Extrapolation Interval (px)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(2,2,4)
imagesc(SkernVec,CorrVec,RMSENorm.strainZ)
title('RMSE(strain_{zz})/max(strain_{zz})')
xlabel('SK (px)')
cx=colorbar;
colormap(Cmap.PBW_map)
cx.Label.String='%';
set(gca,'YDir','normal')
set(gca,'FontSize',12)

saveName=strcat(SweepDir,'/',TestDeg, ...
    '_SpatialCorrSmoothSweep_strainRMSE_Norm');
saveas(gcf,saveName,'fig')
saveas(gcf,saveName,'png')

%% Zoom IN
figure('units','centimeters','InnerPosition',[1,1,18.2,18.2])
subplot(2,2,1)
imagesc(SkernVec,CorrVec,RMSENorm.strainX)
title('RMSE(strain_{xx})/max(strain_{xx})')
ylabel('Extrapolation Interval (px)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)
clim([min(RMSENorm.strainX,[],'all'),3.6])

subplot(2,2,2)
imagesc(SkernVec,CorrVec,RMSENorm.strainY)
title('RMSE(strain_{yy})/max(strain_{yy})')
cx=colorbar;
cx.Label.String='%';
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)
clim([min(RMSENorm.strainY,[],'all'),3.5])

subplot(2,2,3)
imagesc(SkernVec,CorrVec,RMSENorm.strainS)
title('RMSE(strain_{xy})/max(strain_{xy})')
xlabel('SK (px)')
ylabel('Extrapolation Interval (px)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)
clim([min(RMSENorm.strainS,[],'all'),1])

subplot(2,2,4)
imagesc(SkernVec,CorrVec,RMSENorm.strainZ)
title('RMSE(strain_{zz})/max(strain_{zz})')
xlabel('SK (px)')
cx=colorbar;
colormap(Cmap.PBW_map)
cx.Label.String='%';
set(gca,'YDir','normal')
set(gca,'FontSize',12)
clim([min(RMSENorm.strainZ,[],'all'),3.7])
saveName=strcat(SweepDir,'/',TestDeg, ...
    '_SpatialCorrSmoothSweep_strainRMSE_NormZoom');
saveas(gcf,saveName,'fig')
saveas(gcf,saveName,'png')

%% Save results
saveName=strcat(SweepDir,'/',TestDeg, ...
    '_SpatialCorr_TempSpaSmooth_Sweep_results');
save(saveName)