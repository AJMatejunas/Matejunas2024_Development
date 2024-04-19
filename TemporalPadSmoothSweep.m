% This function is written to plot the errors in the acceleration and 
    %stress gauge kinematic fields as a function of the temporal padding
    %and smoothing kernals

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


%% Calculate Average stresses
Int.stress.xAvg=squeeze(mean(Int.stress.x));
Int.stress.sAvg=squeeze(mean(Int.stress.s));

%% Remove finite Element fields;
clear FE


%% Generate smoothing and extrapolation kernal vectrors
PadVec=0:41; %Vector of extrapolation kernals
TkernVec=[1,5:2:41]; %Vector of temporal extrapolation kernals


%% set window of interest
RpitchImp=10;
RpitchFree=10;
RpxImp=RpitchImp*grid.pxPerPeriod;
RpxFree=RpitchFree*grid.pxPerPeriod;
Xind=RpxFree:(length(pos.x)-RpxImp);

RframeStart=5;
RframeEnd=5;
Tind=RframeStart:(length(time.vec)-RframeEnd);

%% Set reference Stresses
StressXref=Int.stress.xAvg(Xind,Tind);
StressSref=Int.stress.sAvg(Xind,Tind);

%% Set up smoothing Parameters for the acceleration sweep
extrapOpts.accel.tempPadOn=false;
extrapOptsR=extrapOpts;
extrapOptsL=extrapOpts;
extrapOptsQ=extrapOpts;
extrapOptsR.accel.tempPadMethod='replicate';
extrapOptsL.accel.tempPadMethod='linear';
extrapOptsQ.accel.tempPadMethod='quadratic';

AccelPadOptsR=cell(length(PadVec),1);
AccelPadOptsL=cell(length(PadVec),1);
AccelPadOptsQ=cell(length(PadVec),1);

smoothOpts.accel.temporalSmooth=0;
AccelSmooth=cell(length(TkernVec),1);
AccelSmooth{1}=smoothOpts.accel;

PadOn=ones(length(PadVec),1);
PadOn(1)=0;

 for m=1:length(PadVec)
     %replicate padding
     TempAccelOptsR=extrapOptsR.accel;
     TempAccelOptsR.tempPadOn=PadOn(m);
     TempAccelOptsR.tempPadFrames=PadVec(m)*[1,1];
     AccelPadOptsR{m}=TempAccelOptsR;
          
     %linear padding
     TempAccelOptsL=extrapOptsL.accel;
     TempAccelOptsL.tempPadOn=PadOn(m);
     TempAccelOptsL.tempPadFrames=PadVec(m)*[1,1];
     AccelPadOptsL{m}=TempAccelOptsL;

     %Quadratic Padding
     TempAccelOptsQ=extrapOptsQ.accel;
     TempAccelOptsQ.tempPadOn=PadOn(m);
     TempAccelOptsQ.tempPadFrames=PadVec(m)*[1,1];
     AccelPadOptsQ{m}=TempAccelOptsQ;
 end
      
for k=2:length(TkernVec)
    tempSmoothA=smoothOpts.accel;
    tempSmoothA.temporalSmooth=true;
    tempSmoothA.temporalKernelSize=TkernVec(k)*[1,1];
    AccelSmooth{k}=tempSmoothA;
end


%% Calculate reference accelerations
AccelXref=Int.accel.x(:,Xind,Tind);
AccelYref=Int.accel.x(:,Xind,Tind);

%% set up progress bar
PadTot=length(PadVec);
PadNum=0;
padStr=num2str(PadTot);

progmsg=strcat('0/',padStr,'~Padding Iterations Complete');
Progress=waitbar(0,progmsg);

%% Run acceleration sweep
fprintf('Running Acceleration smooth-corr sweep \n')
RMSELsgX=zeros(length(PadVec),length(TkernVec));
RMSELsgS=zeros(length(PadVec),length(TkernVec));
RMSERsgX=zeros(length(PadVec),length(TkernVec));
RMSERsgS=zeros(length(PadVec),length(TkernVec));
RMSEQsgX=zeros(length(PadVec),length(TkernVec));
RMSEQsgS=zeros(length(PadVec),length(TkernVec));

RMSELaX=zeros(length(PadVec),length(TkernVec));
RMSELaY=zeros(length(PadVec),length(TkernVec));
RMSERaX=zeros(length(PadVec),length(TkernVec));
RMSERaY=zeros(length(PadVec),length(TkernVec));
RMSEQaX=zeros(length(PadVec),length(TkernVec));
RMSEQaY=zeros(length(PadVec),length(TkernVec));

extrapD=extrapOpts.disp;
rho=MatProps.rho;

for m=1:length(PadVec)       
       extrapR=cell2mat(AccelPadOptsR(m));
       extrapL=cell2mat(AccelPadOptsL(m));
       extrapQ=cell2mat(AccelPadOptsQ(m));
    parfor k=1:length(TkernVec)
       smoothA=cell2mat(AccelSmooth(k));
       %correct displacements
       tempDisp=[];
       [tempDisp.x,tempDisp.y,tempDisp.Rx,tempDisp.Ry]=...
           func_cropAndExtrapDispFields_v4(pos,RawDispX,RawDispY, ...
           extrapD,true);
       %calaculate Accelerations for replicate padding
       [accelR,~,~]=func_smoothCalcAccel_v4(pos,time,tempDisp,...
       smoothA,extrapR,diffOpts,true);
       %calaculate Accelerations for linear padding
       [accelL,~,~]=func_smoothCalcAccel_v4(pos,time,tempDisp,...
       smoothA,extrapL,diffOpts,true);
       %calaculate Accelerations for quadratic padding
       [accelQ,~,~]=func_smoothCalcAccel_v4(pos,time,tempDisp,...
       smoothA,extrapQ,diffOpts,true);

       %Calculate stress gauge stresses
       SGR=func_Full_SG(accelR,X_vec,time,rho);
       SGL=func_Full_SG(accelL,X_vec,time,rho);
       SGQ=func_Full_SG(accelQ,X_vec,time,rho);

       %Crop to ROI
       SGxR=SGR.x(Xind,Tind);
       SGsR=SGR.s(Xind,Tind);

       SGxL=SGL.x(Xind,Tind);
       SGsL=SGL.s(Xind,Tind);
 
       SGxQ=SGQ.x(Xind,Tind);
       SGsQ=SGQ.s(Xind,Tind);
       
       accelR.x=accelR.x(:,Xind,Tind);
       accelR.y=accelR.y(:,Xind,Tind);
       accelL.x=accelL.x(:,Xind,Tind);
       accelL.y=accelL.y(:,Xind,Tind);
       accelQ.x=accelQ.x(:,Xind,Tind);
       accelQ.y=accelQ.y(:,Xind,Tind);


       % Calculate RMSE
       RMSERsgX(m,k)=func_calcFFRMSE(StressXref,SGxR);
       RMSERsgS(m,k)=func_calcFFRMSE(StressSref,SGsR);
       RMSELsgX(m,k)=func_calcFFRMSE(StressXref,SGxL);
       RMSELsgS(m,k)=func_calcFFRMSE(StressSref,SGsL);
       RMSEQsgX(m,k)=func_calcFFRMSE(StressXref,SGxQ);
       RMSEQsgS(m,k)=func_calcFFRMSE(StressSref,SGsQ);

       RMSERaX(m,k)=func_calcFFRMSE(AccelXref,accelR.x);
       RMSERaY(m,k)=func_calcFFRMSE(AccelYref,accelR.y);
       RMSELaX(m,k)=func_calcFFRMSE(AccelXref,accelL.x);
       RMSELaY(m,k)=func_calcFFRMSE(AccelYref,accelL.y);
       RMSEQaX(m,k)=func_calcFFRMSE(AccelXref,accelQ.x);
       RMSEQaY(m,k)=func_calcFFRMSE(AccelYref,accelQ.y);

    end
    PadNum=PadNum+1;
    progmsg=strcat(num2str(PadNum),'/',padStr,...
        '~Padding Iterations Complete');
    Progress=waitbar(PadNum/PadTot,Progress,progmsg);
end

%% Plot Heat Maps of stress gauge RMSE
figure('units','centimeters','InnerPosition',[1,1,18.2,27.5])

subplot(3,2,1)
imagesc(TkernVec,PadVec,RMSERsgX)
title('RMSE(\sigma^{SG}_{xx}) Replicate')
ylabel('Temporal Pad (frames)')
%xlabel('TK (frames)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(3,2,2)
imagesc(TkernVec,PadVec,RMSERsgS)
title('RMSE(\sigma^{SG}_{xy}) Replicate')
% xlabel('TK (frames)')
cx=colorbar;
cx.Label.String='RMSE [Pa]';
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(3,2,3)
imagesc(TkernVec,PadVec,RMSELsgX)
title('RMSE(\sigma^{SG}_{xx}) Linear')
ylabel('Temporal Pad (frames)')
%xlabel('TK (frames)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(3,2,4)
imagesc(TkernVec,PadVec,RMSELsgS)
title('RMSE(\sigma^{SG}_{xy}) Linear')
% xlabel('TK (frames)')
cx=colorbar;
cx.Label.String='RMSE [Pa]';
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(3,2,5)
imagesc(TkernVec,PadVec,RMSEQsgX)
title('RMSE(\sigma^{SG}_{xx}) Quad')
ylabel('Temporal Pad (frames)')
xlabel('TK (frames)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(3,2,6)
imagesc(TkernVec,PadVec,RMSEQsgS)
title('RMSE(\sigma^{SG}_{xy}) Quad')
xlabel('TK (frames)')
cx=colorbar;
cx.Label.String='RMSE [Pa]';
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

saveName=strcat(SweepDir,'/',TestDeg, ...
    '_TemporalSmoothPadSweep_SG_RMSE');
saveas(gcf,saveName,'fig')
saveas(gcf,saveName,'png')

%% Calculate Normnlized RMSE

RMSENorm.R.SGX=RMSERsgX/max(abs(StressXref),[],'all')*100;
RMSENorm.R.SGS=RMSERsgS/max(abs(StressSref),[],'all')*100;
RMSENorm.R.aX=RMSERaX/max(abs(AccelXref),[],'all')*100;
RMSENorm.R.aY=RMSERaY/max(abs(AccelYref),[],'all')*100;

RMSENorm.L.SGX=RMSELsgX/max(abs(StressXref),[],'all')*100;
RMSENorm.L.SGS=RMSELsgS/max(abs(StressSref),[],'all')*100;
RMSENorm.L.aX=RMSELaX/max(abs(AccelXref),[],'all')*100;
RMSENorm.L.aY=RMSELaY/max(abs(AccelYref),[],'all')*100;

RMSENorm.Q.SGX=RMSEQsgX/max(abs(StressXref),[],'all')*100;
RMSENorm.Q.SGS=RMSEQsgS/max(abs(StressSref),[],'all')*100;
RMSENorm.Q.aX=RMSEQaX/max(abs(AccelXref),[],'all')*100;
RMSENorm.Q.aY=RMSEQaY/max(abs(AccelYref),[],'all')*100;

%% Plot Normalized stress gauge sweeps
figure('units','centimeters','InnerPosition',[1,1,18.2,27.5])

subplot(3,2,1)
imagesc(TkernVec,PadVec,RMSENorm.R.SGX)
title('RMSE(\sigma^{SG}_{xx})/max(\sigma_{xx}) Replicate')
ylabel('Temporal Pad (frames)')
%xlabel('TK (frames)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(3,2,2)
imagesc(TkernVec,PadVec,RMSENorm.R.SGS)
title('RMSE(\sigma^{SG}_{xy})/max(\sigma_{xy}) Replicate')
% xlabel('TK (frames)')
cx=colorbar;
colormap(Cmap.PBW_map)
cx.Label.String='%';
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(3,2,3)
imagesc(TkernVec,PadVec,RMSENorm.L.SGX)
title('RMSE(\sigma^{SG}_{xx})/max(\sigma_{xx}) Linear')
ylabel('Temporal Pad (frames)')
%xlabel('TK (frames)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(3,2,4)
imagesc(TkernVec,PadVec,RMSENorm.L.SGS)
title('RMSE(\sigma^{SG}_{xy})/max(\sigma_{xy}) Linear')
% xlabel('TK (frames)')
cx=colorbar;
cx.Label.String='%';
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(3,2,5)
imagesc(TkernVec,PadVec,RMSENorm.Q.SGX)
title('RMSE(\sigma^{SG}_{xx})/max(\sigma_{xx}) Quad')
ylabel('Temporal Pad (frames)')
xlabel('TK (frames)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(3,2,6)
imagesc(TkernVec,PadVec,RMSENorm.Q.SGS)
title('RMSE(\sigma^{SG}_{xy})/max(\sigma_{xy}) Quad')
xlabel('TK (frames)')
cx=colorbar;
cx.Label.String='%';
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

saveName=strcat(SweepDir,'/',TestDeg, ...
    '_TemporalSmoothPadSweep_SG_RMSE_Norm');
saveas(gcf,saveName,'fig')
saveas(gcf,saveName,'png')


%% Plot Normalized stress gauge sweeps Zoom IN
figure('units','centimeters','InnerPosition',[1,1,18.2,27.5])

subplot(3,2,1)
imagesc(TkernVec,PadVec,RMSENorm.R.SGX)
title('RMSE(\sigma^{SG}_{xx})/max(\sigma_{xx}) Replicate')
ylabel('Temporal Pad (frames)')
%xlabel('TK (frames)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)
caxis([0.5,3])

subplot(3,2,2)
imagesc(TkernVec,PadVec,RMSENorm.R.SGS)
title('RMSE(\sigma^{SG}_{xy})/max(\sigma_{xy}) Replicate')
% xlabel('TK (frames)')
cx=colorbar;
colormap(Cmap.PBW_map)
cx.Label.String='%';
set(gca,'YDir','normal')
set(gca,'FontSize',12)
caxis([0.5,3])

subplot(3,2,3)
imagesc(TkernVec,PadVec,RMSENorm.L.SGX)
title('RMSE(\sigma^{SG}_{xx})/max(\sigma_{xx}) Linear')
ylabel('Temporal Pad (frames)')
%xlabel('TK (frames)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)
caxis([0.5,3])

subplot(3,2,4)
imagesc(TkernVec,PadVec,RMSENorm.L.SGS)
title('RMSE(\sigma^{SG}_{xy})/max(\sigma_{xy}) Linear')
% xlabel('TK (frames)')
cx=colorbar;
cx.Label.String='%';
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)
caxis([0.5,3])

subplot(3,2,5)
imagesc(TkernVec,PadVec,RMSENorm.Q.SGX)
title('RMSE(\sigma^{SG}_{xx})/max(\sigma_{xx}) Quad')
ylabel('Temporal Pad (frames)')
xlabel('TK (frames)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)
caxis([0.5,3])

subplot(3,2,6)
imagesc(TkernVec,PadVec,RMSENorm.Q.SGS)
title('RMSE(\sigma^{SG}_{xy})/max(\sigma_{xy}) Quad')
xlabel('TK (frames)')
cx=colorbar;
cx.Label.String='%';
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)
caxis([0.5,3])

saveName=strcat(SweepDir,'/',TestDeg, ...
    '_TemporalSmoothPadSweep_SG_RMSE_NormZoom');
saveas(gcf,saveName,'fig')
saveas(gcf,saveName,'png')

%% Plot Heat Maps of Acceleration RMSE
figure('units','centimeters','InnerPosition',[1,1,18.2,27.5])

subplot(3,2,1)
imagesc(TkernVec,PadVec,RMSERaX)
title('RMSE(a_x) Replicate')
ylabel('Temporal Pad (frames)')
%xlabel('TK (frames)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(3,2,2)
imagesc(TkernVec,PadVec,RMSERaY)
title('RMSE(a_y) Replicate')
% xlabel('TK (frames)')
cx=colorbar;
cx.Label.String='RMSE (m/s^2)';
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(3,2,3)
imagesc(TkernVec,PadVec,RMSELaX)
title('RMSE(a_x) Linear')
ylabel('Temporal Pad (frames)')
%xlabel('TK (frames)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(3,2,4)
imagesc(TkernVec,PadVec,RMSELaY)
title('RMSE(a_y) Linear')
% xlabel('TK (frames)')
cx=colorbar;
cx.Label.String='RMSE (m/s^2)';
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(3,2,5)
imagesc(TkernVec,PadVec,RMSEQaX)
title('RMSE(a_x) Quad')
ylabel('Temporal Pad (frames)')
xlabel('TK (frames)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(3,2,6)
imagesc(TkernVec,PadVec,RMSEQaY)
title('RMSE(a_y) Quad')
xlabel('TK (frames)')
cx=colorbar;
cx.Label.String='RMSE (m/s^2)';
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

saveName=strcat(SweepDir,'/',TestDeg, ...
    '_TemporalSmoothPadSweep_Accel_RMSE');
saveas(gcf,saveName,'fig')
saveas(gcf,saveName,'png')

%% Plot Heat Maps of Normalized Acceleration RMSE
figure('units','centimeters','InnerPosition',[1,1,18.2,27.5])

subplot(3,2,1)
imagesc(TkernVec,PadVec,RMSENorm.R.aX)
title('RMSE(a_x)/max(a_x) Replicate')
ylabel('Temporal Pad (frames)')
%xlabel('TK (frames)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(3,2,2)
imagesc(TkernVec,PadVec,RMSENorm.R.aY)
title('RMSE(a_y) Replicate')
% xlabel('TK (frames)')
cx=colorbar;
cx.Label.String='%';
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(3,2,3)
imagesc(TkernVec,PadVec,RMSENorm.L.aX)
title('RMSE(a_x)/max(a_x) Linear')
ylabel('Temporal Pad (frames)')
%xlabel('TK (frames)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(3,2,4)
imagesc(TkernVec,PadVec,RMSENorm.L.aY)
title('RMSE(a_y) Linear')
% xlabel('TK (frames)')
cx=colorbar;
cx.Label.String='%';
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(3,2,5)
imagesc(TkernVec,PadVec,RMSENorm.Q.aX)
title('RMSE(a_x)/max(a_x) Quad')
ylabel('Temporal Pad (frames)')
xlabel('TK (frames)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(3,2,6)
imagesc(TkernVec,PadVec,RMSENorm.Q.aY)
title('RMSE(a_y) Quad')
xlabel('TK (frames)')
cx=colorbar;
cx.Label.String='%';
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

saveName=strcat(SweepDir,'/',TestDeg, ...
    '_TemporalSmoothPadSweep_Accel_RMSE_Norm');
saveas(gcf,saveName,'fig')
saveas(gcf,saveName,'png')

%% Save results
saveName=strcat(SweepDir,'/',TestDeg, ...
    '_SpatialCorr_TempSpaSmooth_Sweep_results');
save(saveName)