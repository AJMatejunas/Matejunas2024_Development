% This function is written to plot the errors in the kinematic fields as a
    % function of extrapolation kernal and smoothing kernal. Strain is a
    % function spatial smoothing kernal and stress gauge stresses is a
    % funciton of temporal smoothing kernal. In this case 30 copis of noise
    % will be evaluated

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

%% Choose finite element data interpolated to grid coordinates
[IntName,IntPath]=uigetfile('*.mat', ...
    'Choose file containing interpolated finite element kinematic fields');
IntFile=strcat(IntPath,'/',IntName);
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
Int=load(IntFile,'pos','disp','accel','strain','stress','time');

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

 %% obtain Grid Coordinates
 fprintf('Obtaining Grid Coordinates \n')
 [~,pos,~] = func_gridMethodImageProcessing_AJM(imagePath,...
                imageFile,...
                grid,gridMethodOpts,imageNoise);
 

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
RMSEstrainX=zeros(length(CorrVec),length(SkernVec),noiseNum);
RMSEstressX=zeros(length(CorrVec),length(SkernVec),noiseNum);
RMSEstrainY=zeros(length(CorrVec),length(SkernVec),noiseNum);
RMSEstrainS=zeros(length(CorrVec),length(SkernVec),noiseNum);
RMSEstressS=zeros(length(CorrVec),length(SkernVec),noiseNum);
RMSEstrainZ=zeros(length(CorrVec),length(SkernVec),noiseNum);

%% set up noise sweep
    imageNoise.addNoise=true;
    imageNoise.numCopies=30;
    noiseNum=imageNoise.numCopies;
%%
hardCodeFreeEdge=globalOpts.hardCodeFreeEdge;
%% Loop over processing parameters

for m=1:length(CorrVec)
    extrapD=cell2mat(DispOpts(m));
    extrapS=cell2mat(StrainCorrOpts(m));
    for k=1:length(SkernVec)
        SmoothS=cell2mat(StrainSmooth(k));
        parfor n=1:imageNoise.numCopies
            timeSpa=time;
            % Grediac et al.
            fprintf('\n--------------------------------------------------------------\n')
            fprintf('GRID METHOD PROCESSING\n')
            fprintf('--------------------------------------------------------------\n')

            %--------------------------------------------------------------------------
            % GRID IMAGE PROCESSING


            %fprintf('Checking for existing processed data file.\n')
            processGridImages = true;


            % fprintf('Processing images using the grid method toolbox.\n')
            % Process the image squence with the grid method toolbox
            [gridSpa,posSpa,dispRaw] = func_gridMethodImageProcessing_AJM(imagePath,...
                imageFile,...
                grid,gridMethodOpts,imageNoise);
            %--------------------------------------------------------------------------
            % Update Geometry and Number of Frames Based on Displacement Matrix Size
            %fprintf('Updating geometry variables based on the size of displacement field matrices.\n')
            [specimenSpa,gridSpa] = func_updateSpecGeom(specimen,gridSpa,dispRaw);

            % Currently the rotations are unused so remove them to save RAM
            %disp = rmfield(disp,'rot');


            %--------------------------------------------------------------------------
            % Create the time vector based on the number of frames in the disp struct
            timeSpa.numFrames = size(dispRaw.x,3);
            timeSpa.vec = 0:timeSpa.step:(size(dispRaw.x,3)-1)*timeSpa.step;

            %Create arrays of x and y vectors
            XSpa=posSpa.x;
            YSpac=posSpa.y;
            posSpa.lengthX = posSpa.x(end)+posSpa.xStep/2;
            posSpa.lengthY = posSpa.y(end)+posSpa.yStep/2;
            posSpa.xGridF = padarray(posSpa.xGrid,[0,0,timeSpa.numFrames-1], ...
                'replicate','post');
            posSpa.yGridF = padarray(posSpa.yGrid,[0,0,timeSpa.numFrames-1], ...
                'replicate','post');
            posSpa.x0F = squeeze(padarray(posSpa.x,[0,0,timeSpa.numFrames-1], ...
                'replicate','post'));

            [freeEdge,specimenSpa,dispRaw] =...
                func_getFreeEdge(hardCodeFreeEdge,...
                imagePath,imageFile,specimenSpa,dispRaw);

            SpaDispX=dispRaw.x;
            SpaDispY=dispRaw.y;

            %Correct Displacements
            tempDisp=[];
            [tempDisp.x,tempDisp.y,tempDisp.rX,tempDisp.rY]=...
                func_cropAndExtrapDispFields_v4(posSpa,SpaDispX,SpaDispY, ...
                extrapD,true);
            %Calculate Strains
            [strain,~]=func_smoothCalcStrain_v4(posSpa,timeSpa,tempDisp,...
                SmoothS,extrapS,true);
            %Calculate constituitive models stresses
            tempModel=func_ViscoConstitutiveV6(strain.x,strain.y,strain.s, ...
                timeSpa.vec,MatProps,0,0,0);

            %crop Strain to the ROI
            strain.x=strain.x(:,Xind,Tind);
            strain.y=strain.y(:,Xind,Tind);
            strain.s=strain.s(:,Xind,Tind);
            strain.z=tempModel.zzstrain(:,Xind,Tind);

            %Crop stress to ROI
            stressX=tempModel.Avxx(Xind,Tind);
            stressS=tempModel.Avxy(Xind,Tind);

            %Calculate RMSE for strains
            RMSEstrainX(m,k,n)=func_calcFFRMSE(StrainXref,strain.x);
            RMSEstrainY(m,k,n)=func_calcFFRMSE(StrainYref,strain.y);
            RMSEstrainS(m,k,n)=func_calcFFRMSE(StrainSref,strain.s);
            RMSEstrainZ(m,k,n)=func_calcFFRMSE(StrainZref,strain.z);
            % Calculate RMSE for Stress
            RMSEstressX(m,k,n)=func_calcFFRMSE(StressXref,stressX);
            RMSEstressS(m,k,n)=func_calcFFRMSE(StressSref,stressS);
        end
    end
end

clear strain stressX stressS tempModel tempDisp 
Int=rmfield(Int,{'strain','disp','accel'});

%% Calculate Strain Error statistics for each smoothing kernal
RMSEsys.strainX=mean(RMSEstrainX,3);
RMSEran.strainX=std(RMSEstrainX,0,3);
RMSEtot.strainX=RMSEsys.strainX+2*RMSEran.strainX;

RMSEsys.strainY=mean(RMSEstrainY,3);
RMSEran.strainY=std(RMSEstrainY,0,3);
RMSEtot.strainY=RMSEsys.strainY+2*RMSEran.strainY;

RMSEsys.strainS=mean(RMSEstrainS,3);
RMSEran.strainS=std(RMSEstrainS,0,3);
RMSEtot.strainS=RMSEsys.strainS+2*RMSEran.strainS;

RMSEsys.strainZ=mean(RMSEstrainZ,3);
RMSEran.strainZ=std(RMSEstrainZ,0,3);
RMSEtot.strainZ=RMSEsys.strainZ+2*RMSEran.strainZ;

%Normalized
RMSEsysNorm.strainX=mean(RMSEstrainX,3)/...
    (max(abs(StrainXref),[],'all'))*100;
RMSEranNorm.strainX=std(RMSEstrainX,0,3)/...
(max(abs(StrainXref),[],'all'))*100;
RMSEtotNorm.strainX=RMSEtot.strainX/...
    (max(abs(StrainXref),[],'all'))*100;

RMSEsysNorm.strainY=mean(RMSEstrainY,3)/...
    (max(abs(StrainYref),[],'all'))*100;
RMSEranNorm.strainY=std(RMSEstrainY,0,3)/...
    (max(abs(StrainYref),[],'all'))*100;
RMSEtotNorm.strainY=RMSEtot.strainY/...
    (max(abs(StrainYref),[],'all'))*100;

RMSEsysNorm.strainS=mean(RMSEstrainS,3)/...
    (max(abs(StrainSref),[],'all'))*100;
RMSEranNorm.strainS=std(RMSEstrainS,0,3)/...
    (max(abs(StrainSref),[],'all'))*100;
RMSEtotNorm.strainS=RMSEtot.strainS/...
    (max(abs(StrainSref),[],'all'))*100;

RMSEsysNorm.strainZ=mean(RMSEstrainZ,3)/...
    (max(abs(StrainZref),[],'all'))*100;
RMSEranNorm.strainZ=std(RMSEstrainZ,0,3)/...
    (max(abs(StrainZref),[],'all'))*100;
RMSEtotNorm.strainZ=RMSEtot.strainZ/...
    (max(abs(StrainZref),[],'all'))*100;

%% Stress Error Statistics
RMSEsys.stressX=mean(RMSEstressX,3);
RMSEran.stressX=std(RMSEstressX,0,3);
RMSEtot.stressX=RMSEsys.stressX+2*RMSEran.stressX;

RMSEsys.stressS=mean(RMSEstressS,3);
RMSEran.stressS=std(RMSEstressS,0,3);
RMSEtot.stressS=RMSEsys.stressS+2*RMSEran.stressS;

%Normalized
RMSEsysNorm.stressX=mean(RMSEstressX,3)/...
    (max(abs(StressXref),[],'all'))*100;
RMSEranNorm.stressX=std(RMSEstressX,0,3)/...
    (max(abs(StressXref),[],'all'))*100;
RMSEtotNorm.stressX=RMSEtot.stressX/...
    (max(abs(StressXref),[],'all'))*100;

RMSEsysNorm.stressS=mean(RMSEstressS,3)/...
    (max(abs(StressSref),[],'all'))*100;
RMSEranNorm.stressS=std(RMSEstressS,0,3)/...
    (max(abs(StressSref),[],'all'))*100;
RMSEtotNorm.stressS=RMSEtot.stressS/...
    (max(abs(StressSref),[],'all'))*100;

%% Plot Heat Maps of strain RMSE errors
figure('units','centimeters','InnerPosition',[1,1,18.2,18.2])
subplot(2,2,1)
imagesc(SkernVec,CorrVec,RMSEsys.strainX)
title('Random RMSE(strain_{xx}')
ylabel('Extrapolation Interval (px)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(2,2,2)
imagesc(SkernVec,CorrVec,RMSEsys.strainY)
title('Random RMSE(strain_{yy}')
ylabel('Extrapolation Interval (px)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(2,2,3)
imagesc(SkernVec,CorrVec,RMSEsys.strainS)
title('Systematic RMSE(strain_{xy}')
xlabel('SK (px)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(2,2,4)
imagesc(SkernVec,CorrVec,RMSEsys.strainZ)
title('Systematic RMSE(strain_{zz}')
xlabel('SK (px)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

saveName=strcat(SweepDir,'/',TestDeg, ...
    '_SpatialCorrSmoothSweep30Noise_strainRMSEsys');
saveas(gcf,saveName,'fig')
saveas(gcf,saveName,'png')


%% Plot Heat Maps of constitutive Model stress RMSE
figure('units','centimeters','InnerPosition',[1,1,18.2,9.1])
subplot(1,2,1)
imagesc(SkernVec,CorrVec,RMSEsys.stressX)
title('Systematic RMSE(\sigma^{model}_{xx}')
xlabel('SK (px)')
ylabel('Extrapolation Interval (px)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(1,2,2)
imagesc(SkernVec,CorrVec,RMSEsys.stressS)
title('Systematic RMSE(\sigma^{model}_{xy}')
xlabel('SK (px)')

colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

saveName=strcat(SweepDir,'/',TestDeg, ...
    '_SpatialCorrSmoothSweep30Noise_modelRMSEsys');
saveas(gcf,saveName,'fig')
saveas(gcf,saveName,'png')



%% Plot Heat Maps of strain RMSE random errors
figure('units','centimeters','InnerPosition',[1,1,18.2,18.2])
subplot(2,2,1)
imagesc(SkernVec,CorrVec,RMSEran.strainX)
title('Random RMSE(strain_{xx}')
ylabel('Extrapolation Interval (px)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(2,2,2)
imagesc(SkernVec,CorrVec,RMSEran.strainY)
title('Random RMSE(strain_{yy}')
ylabel('Extrapolation Interval (px)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(2,2,3)
imagesc(SkernVec,CorrVec,RMSEran.strainS)
title('Random RMSE(strain_{xy}')
xlabel('SK (px)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(2,2,4)
imagesc(SkernVec,CorrVec,RMSEran.strainZ)
title('Random RMSE(strain_{zz}')
xlabel('SK (px)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

saveName=strcat(SweepDir,'/',TestDeg, ...
    '_SpatialCorrSmoothSweep30Noise_strainRMSEran');
saveas(gcf,saveName,'fig')
saveas(gcf,saveName,'png')


%% Plot Heat Maps of constitutive Model stress RMSE
figure('units','centimeters','InnerPosition',[1,1,18.2,9.1])
subplot(1,2,1)
imagesc(SkernVec,CorrVec,RMSEran.stressX)
title('Random RMSE(\sigma^{model}_{xx}')
xlabel('SK (px)')
ylabel('Extrapolation Interval (px)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(1,2,2)
imagesc(SkernVec,CorrVec,RMSEran.stressS)
title('Random RMSE(\sigma^{model}_{xy}')
xlabel('SK (px)')

colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

saveName=strcat(SweepDir,'/',TestDeg, ...
    '_SpatialCorrSmoothSweep30Noise_modelRMSEran');
saveas(gcf,saveName,'fig')
saveas(gcf,saveName,'png')


%% Plot Heat Maps of strain RMSE random errors
figure('units','centimeters','InnerPosition',[1,1,18.2,18.2])
subplot(2,2,1)
imagesc(SkernVec,CorrVec,RMSEtot.strainX)
title('Total RMSE(strain_{xx}')
ylabel('Extrapolation Interval (px)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(2,2,2)
imagesc(SkernVec,CorrVec,RMSEtot.strainY)
title('Total RMSE(strain_{yy}')
ylabel('Extrapolation Interval (px)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(2,2,3)
imagesc(SkernVec,CorrVec,RMSEtot.strainS)
title('Total RMSE(strain_{xy}')
xlabel('SK (px)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(2,2,4)
imagesc(SkernVec,CorrVec,RMSEtot.strainZ)
title('Total RMSE(strain_{zz}')
xlabel('SK (px)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

saveName=strcat(SweepDir,'/',TestDeg, ...
    '_SpatialCorrSmoothSweep30Noise_strainRMSEtot');
saveas(gcf,saveName,'fig')
saveas(gcf,saveName,'png')


%% Plot Heat Maps of constitutive Model stress RMSE
figure('units','centimeters','InnerPosition',[1,1,18.2,9.1])
subplot(1,2,1)
imagesc(SkernVec,CorrVec,RMSEtot.stressX)
title('Total RMSE(\sigma^{model}_{xx}')
xlabel('SK (px)')
ylabel('Extrapolation Interval (px)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(1,2,2)
imagesc(SkernVec,CorrVec,RMSEtot.stressS)
title('Total RMSE(\sigma^{model}_{xy}')
xlabel('SK (px)')

colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

saveName=strcat(SweepDir,'/',TestDeg, ...
    '_SpatialCorrSmoothSweep30Noise_modelRMSEtot');
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
RMSEsgX=zeros(length(CorrVec),length(SkernVec),noiseNum);
RMSEsgS=zeros(length(CorrVec),length(SkernVec),noiseNum);
rho=MatProps.rho;
X_vec=pos.x;
for m=1:length(CorrVec)
    extrapD=cell2mat(DispOpts(m));
    for k=1:length(TkernVec)
        extrapA=cell2mat(AccelCorrOpts(k,m));
       smoothA=cell2mat(AccelSmooth(k)); 
       parfor n=1:noiseNum
            timeTemp=time;
            % Grediac et al.
            fprintf('\n--------------------------------------------------------------\n')
            fprintf('GRID METHOD PROCESSING\n')
            fprintf('--------------------------------------------------------------\n')

            %--------------------------------------------------------------------------
            % GRID IMAGE PROCESSING


            %fprintf('Checking for existing processed data file.\n')
            processGridImages = true;


            % fprintf('Processing images using the grid method toolbox.\n')
            % Process the image squence with the grid method toolbox
            [gridTemp,posTemp,dispRaw] = func_gridMethodImageProcessing_AJM(imagePath,...
                imageFile,...
                grid,gridMethodOpts,imageNoise);
            %--------------------------------------------------------------------------
            % Update Geometry and Number of Frames Based on Displacement Matrix Size
            %fprintf('Updating geometry variables based on the size of displacement field matrices.\n')
            [specimenTemp,gridTemp] = func_updateSpecGeom(specimen,gridTemp,dispRaw);

            % Currently the rotations are unused so remove them to save RAM
            %disp = rmfield(disp,'rot');


            %--------------------------------------------------------------------------
            % Create the time vector based on the number of frames in the disp struct
            timeTemp.numFrames = size(dispRaw.x,3);
            timeTemp.vec = 0:timeTemp.step:(size(dispRaw.x,3)-1)*timeTemp.step;

            %Create arrays of x and y vectors
            XTemp=posTemp.x;
            YTempc=posTemp.y;
            posTemp.lengthX = posTemp.x(end)+posTemp.xStep/2;
            posTemp.lengthY = posTemp.y(end)+posTemp.yStep/2;
            posTemp.xGridF = padarray(posTemp.xGrid,[0,0,timeTemp.numFrames-1], ...
                'replicate','post');
            posTemp.yGridF = padarray(posTemp.yGrid,[0,0,timeTemp.numFrames-1], ...
                'replicate','post');
            posTemp.x0F = squeeze(padarray(posTemp.x,[0,0,timeTemp.numFrames-1], ...
                'replicate','post'));

            [freeEdge,specimenTemp,dispRaw] =...
                func_getFreeEdge(hardCodeFreeEdge,...
                imagePath,imageFile,specimenTemp,dispRaw);
       
       %correct displacements
       tempDisp=[];
       [tempDisp.x,tempDisp.y,tempDisp.Rx,tempDisp.Ry]=...
           func_cropAndExtrapDispFields_v4(pos,dispRaw.x,dispRaw.y, ...
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
       RMSEsgX(m,k,n)=func_calcFFRMSE(StressXref,SGx);
       RMSEsgS(m,k,n)=func_calcFFRMSE(StressSref,SGs);
        end
    end
end

%% Calculate statistics on the SG error
RMSEsys.sgX=mean(RMSEsgX,3);
RMSEran.sgX=std(RMSEsgX,0,3);
RMSEsys.sgS=mean(RMSEsgS,3);
RMSEran.sgS=std(RMSEsgS,0,3);

RMSEsysNorm.sgX=mean(RMSEsgX,3)/max(abs(StressXref),[],'all')*100;
RMSEranNorm.sgX=std(RMSEsgX,0,3)/max(abs(StressXref),[],'all')*100;
RMSEsysNorm.sgS=mean(RMSEsgS,3)/max(abs(StressSref),[],'all')*100;
RMSEranNorm.sgS=std(RMSEsgS,0,3)/max(abs(StressSref),[],'all')*100;

RMSEtot.sgX=RMSEsys.sgX+2*RMSEran.sgX;
RMSEtotNorm.sgX=RMSEtot.sgX/max(abs(StressXref),[],'all')*100;
RMSEtot.sgS=RMSEsys.sgS+2*RMSEran.sgS;
RMSEtotNorm.sgS=RMSEtot.sgS/max(abs(StressSref),[],'all')*100;

%% Plot Heat Maps of stress gauge RMSE
figure('units','centimeters','InnerPosition',[1,1,18.2,9.1])
subplot(1,2,1)
imagesc(TkernVec,CorrVec,RMSEsys.sgX)
title('Systematic RMSE(\sigma^{SG}_{xx})')
ylabel('Extrapolation Interval (px)')
xlabel('TK (frames)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(1,2,2)
imagesc(TkernVec,CorrVec,RMSEsys.sgS)
title('Systematic RMSE(\sigma^{SG}_{xy})')
xlabel('TK (frames)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

saveName=strcat(SweepDir,'/',TestDeg, ...
    '_SpatialCorr_TemPSmooth_Sweep_SG_RMSEsys');
saveas(gcf,saveName,'fig')
saveas(gcf,saveName,'png')


%% Plot Heat Maps of Normalized systematic stress gauge RMSE
figure('units','centimeters','InnerPosition',[1,1,18.2,9.1])
subplot(1,2,1)
imagesc(TkernVec,CorrVec,RMSEsysNorm.sgX)
title('Systematic RMSE(\sigma^{SG}_{xx})/max(\sigma_{xx})')
ylabel('Extrapolation Interval (px)')
xlabel('TK (frames)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)
%caxis([0,25])

subplot(1,2,2)
imagesc(TkernVec,CorrVec,RMSEsysNorm.sgS)
title('Systematic RMSE(\sigma^{SG}_{xy})/max(\sigma_{xy})')
xlabel('TK (frames)')
cx=colorbar;
cx.Label.String='%';
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)
%caxis([0,25])

saveName=strcat(SweepDir,'/',TestDeg, ...
    '_SpatialCorr_TemPSmooth_Sweep_SG_RMSE_NormSys');
saveas(gcf,saveName,'fig')
saveas(gcf,saveName,'png')
%% Plot Zoomed in Heat Maps of systematic stress gauge RMSE
figure('units','centimeters','InnerPosition',[1,1,18.2,9.1])
subplot(1,2,1)
imagesc(TkernVec,CorrVec,RMSEsysNorm.sgX)
title('Systematic RMSE(\sigma^{SG}_{xx})/max(\sigma_{xx})')
ylabel('Extrapolation Interval (px)')
xlabel('TK (frames)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)
caxis([0,25])

subplot(1,2,2)
imagesc(TkernVec,CorrVec,RMSEsysNorm.sgS)
title('Systematic RMSE(\sigma^{SG}_{xy})/max(\sigma_{xy})')
xlabel('TK (frames)')
cx=colorbar;
cx.Label.String='%';
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)
caxis([0,25])

saveName=strcat(SweepDir,'/',TestDeg, ...
    '_SpatialCorr_TemPSmooth_Sweep_SG_RMSEsys_NormZoom');
saveas(gcf,saveName,'fig')
saveas(gcf,saveName,'png')

%% Plot Heat Maps of Ramdom  stress gauge RMSE
figure('units','centimeters','InnerPosition',[1,1,18.2,9.1])
subplot(1,2,1)
imagesc(TkernVec,CorrVec,RMSEran.sgX)
title('Random RMSE(\sigma^{SG}_{xx})')
ylabel('Extrapolation Interval (px)')
xlabel('TK (frames)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(1,2,2)
imagesc(TkernVec,CorrVec,RMSEran.sgS)
title('Random RMSE(\sigma^{SG}_{xy})')
xlabel('TK (frames)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

saveName=strcat(SweepDir,'/',TestDeg, ...
    '_SpatialCorr_TemPSmooth_Sweep_SG_RMSERan');
saveas(gcf,saveName,'fig')
saveas(gcf,saveName,'png')


%% Plot Heat Maps of stress gauge normalized Random RMSE
figure('units','centimeters','InnerPosition',[1,1,18.2,9.1])
subplot(1,2,1)
imagesc(TkernVec,CorrVec,RMSEranNorm.sgX)
title('Random RMSE(\sigma^{SG}_{xx})/max(\sigma_{xx})')
ylabel('Extrapolation Interval (px)')
xlabel('TK (frames)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)
%caxis([0,25])

subplot(1,2,2)
imagesc(TkernVec,CorrVec,RMSEranNorm.sgS)
title('Random RMSE(\sigma^{SG}_{xy})/max(\sigma_{xy})')
xlabel('TK (frames)')
cx=colorbar;
cx.Label.String='%';
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)
%caxis([0,25])

saveName=strcat(SweepDir,'/',TestDeg, ...
    '_SpatialCorr_TemPSmooth_Sweep_SG_RMSE_NormRan');
saveas(gcf,saveName,'fig')
saveas(gcf,saveName,'png')
%% Plot zoomed Heat Maps of stress gauge normalized random RMSE
figure('units','centimeters','InnerPosition',[1,1,18.2,9.1])
subplot(1,2,1)
imagesc(TkernVec,CorrVec,RMSEranNorm.sgX)
title('Random RMSE(\sigma^{SG}_{xx})/max(\sigma_{xx})')
ylabel('Extrapolation Interval (px)')
xlabel('TK (frames)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)
caxis([0,25])

subplot(1,2,2)
imagesc(TkernVec,CorrVec,RMSEranNorm.sgS)
title('RMSE(\sigma^{SG}_{xy})/max(\sigma_{xy})')
xlabel('TK (frames)')
cx=colorbar;
cx.Label.String='%';
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)
caxis([0,25])

saveName=strcat(SweepDir,'/',TestDeg, ...
    '_SpatialCorr_TemPSmooth_Sweep_SG_RMSEran_NormZoom');
saveas(gcf,saveName,'fig')
saveas(gcf,saveName,'png')


%% Plot Heat Maps of Total  stress gauge RMSE
figure('units','centimeters','InnerPosition',[1,1,18.2,9.1])
subplot(1,2,1)
imagesc(TkernVec,CorrVec,RMSEtot.sgX)
title('Total RMSE(\sigma^{SG}_{xx})')
ylabel('Extrapolation Interval (px)')
xlabel('TK (frames)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(1,2,2)
imagesc(TkernVec,CorrVec,RMSEtot.sgS)
title('Total RMSE(\sigma^{SG}_{xy})')
xlabel('TK (frames)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

saveName=strcat(SweepDir,'/',TestDeg, ...
    '_SpatialCorr_TemPSmooth_Sweep_SG_RMSEtot');
saveas(gcf,saveName,'fig')
saveas(gcf,saveName,'png')


%% Plot Heat Maps of stress gauge normalized Total RMSE
figure('units','centimeters','InnerPosition',[1,1,18.2,9.1])
subplot(1,2,1)
imagesc(TkernVec,CorrVec,RMSEtotNorm.sgX)
title('Total RMSE(\sigma^{SG}_{xx})/max(\sigma_{xx})')
ylabel('Extrapolation Interval (px)')
xlabel('TK (frames)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)
%caxis([0,25])

subplot(1,2,2)
imagesc(TkernVec,CorrVec,RMSEtotNorm.sgS)
title('Total RMSE(\sigma^{SG}_{xy})/max(\sigma_{xy})')
xlabel('TK (frames)')
cx=colorbar;
cx.Label.String='%';
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)
%caxis([0,25])

saveName=strcat(SweepDir,'/',TestDeg, ...
    '_SpatialCorr_TemPSmooth_Sweep_SG_RMSE_Normtot');
saveas(gcf,saveName,'fig')
saveas(gcf,saveName,'png')

%% Plot zoomed Heat Maps of stress gauge normalized total RMSE
figure('units','centimeters','InnerPosition',[1,1,18.2,9.1])
subplot(1,2,1)
imagesc(TkernVec,CorrVec,RMSEtotNorm.sgX)
title('Total RMSE(\sigma^{SG}_{xx})/max(\sigma_{xx})')
ylabel('Extrapolation Interval (px)')
xlabel('TK (frames)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)
caxis([0,50])

subplot(1,2,2)
imagesc(TkernVec,CorrVec,RMSEtotNorm.sgS)
title(' Total RMSE(\sigma^{SG}_{xy})/max(\sigma_{xy})')
xlabel('TK (frames)')
cx=colorbar;
cx.Label.String='%';
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)
caxis([0,50])

saveName=strcat(SweepDir,'/',TestDeg, ...
    '_SpatialCorr_TemPSmooth_Sweep_SG_RMSEtot_NormZoom');
saveas(gcf,saveName,'fig')
saveas(gcf,saveName,'png')

%% Plot Heat Maps of systematic stress RMSE
figure('units','centimeters','InnerPosition',[1,1,18.2,9.1])
subplot(1,2,1)
imagesc(SkernVec,CorrVec,RMSEsysNorm.stressX)
title('Systematic RMSE(\sigma^{Model}_{xx})/max(\sigma_{xx})')
ylabel('Extrapolation Interval (px)')
xlabel('SK (pixels)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)
%caxis([0,25])

subplot(1,2,2)
imagesc(SkernVec,CorrVec,RMSEsysNorm.stressS)
title('Systematic RMSE(\sigma^{Model}_{xy})/max(\sigma_{xy})')
xlabel('SK (pixels)')
cx=colorbar;
cx.Label.String='%';
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)
%caxis([0,25])

saveName=strcat(SweepDir,'/',TestDeg, ...
    '_SpatialCorr_TemPSmooth_Sweep_Stress_RMSEsys_Norm');
saveas(gcf,saveName,'fig')
saveas(gcf,saveName,'png')

%% Plot Heat Maps of Random Normalized stress RMSE
figure('units','centimeters','InnerPosition',[1,1,18.2,9.1])
subplot(1,2,1)
imagesc(SkernVec,CorrVec,RMSEranNorm.stressX)
title('Random RMSE(\sigma^{Model}_{xx})/max(\sigma_{xx})')
ylabel('Extrapolation Interval (px)')
xlabel('SK (pixels)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)
%caxis([0,25])

subplot(1,2,2)
imagesc(SkernVec,CorrVec,RMSEranNorm.stressS)
title('Random RMSE(\sigma^{Model}_{xy})/max(\sigma_{xy})')
xlabel('SK (pixels)')
cx=colorbar;
cx.Label.String='%';
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)
%caxis([0,25])

saveName=strcat(SweepDir,'/',TestDeg, ...
    '_SpatialCorr_TemPSmooth_Sweep_Stress_RMSEran_Norm');
saveas(gcf,saveName,'fig')
saveas(gcf,saveName,'png')

%% Plot Heat Maps of Total Normalized stress RMSE
figure('units','centimeters','InnerPosition',[1,1,18.2,9.1])
subplot(1,2,1)
imagesc(SkernVec,CorrVec,RMSEtotNorm.stressX)
title('Total RMSE(\sigma^{Model}_{xx})/max(\sigma_{xx})')
ylabel('Extrapolation Interval (px)')
xlabel('SK (pixels)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)
%caxis([0,25])

subplot(1,2,2)
imagesc(SkernVec,CorrVec,RMSEtotNorm.stressS)
title('Total RMSE(\sigma^{Model}_{xy})/max(\sigma_{xy})')
xlabel('SK (pixels)')
cx=colorbar;
cx.Label.String='%';
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)
%caxis([0,25])

saveName=strcat(SweepDir,'/',TestDeg, ...
    '_SpatialCorr_TemPSmooth_Sweep_Stress_RMSEtot_Norm');
saveas(gcf,saveName,'fig')
saveas(gcf,saveName,'png')


%% Plot Heat Maps of systematic strain Normalized RMSE errors
figure('units','centimeters','InnerPosition',[1,1,18.2,18.2])
subplot(2,2,1)
imagesc(SkernVec,CorrVec,RMSEsysNorm.strainX)
title('Systematic RMSE(strain_{xx})/max(strain_{xx})')
ylabel('Extrapolation Interval (px)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(2,2,2)
imagesc(SkernVec,CorrVec,RMSEsysNorm.strainY)
title('Systematic RMSE(strain_{yy})/max(strain_{yy})')
cx=colorbar;
cx.Label.String='%';
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(2,2,3)
imagesc(SkernVec,CorrVec,RMSEsysNorm.strainS)
title('Systematic RMSE(strain_{xy})/max(strain_{xy})')
xlabel('SK (px)')
ylabel('Extrapolation Interval (px)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(2,2,4)
imagesc(SkernVec,CorrVec,RMSEsysNorm.strainZ)
title('Systematic RMSE(strain_{zz})/max(strain_{zz})')
xlabel('SK (px)')
cx=colorbar;
colormap(Cmap.PBW_map)
cx.Label.String='%';
set(gca,'YDir','normal')
set(gca,'FontSize',12)

saveName=strcat(SweepDir,'/',TestDeg, ...
    '_SpatialCorrSmoothSweep_strainRMSEsys_Norm');
saveas(gcf,saveName,'fig')
saveas(gcf,saveName,'png')

%% Plot Heat Maps of strain RMSE errors
figure('units','centimeters','InnerPosition',[1,1,18.2,18.2])
subplot(2,2,1)
imagesc(SkernVec,CorrVec,RMSEranNorm.strainX)
title('Random RMSE(strain_{xx})/max(strain_{xx})')
ylabel('Extrapolation Interval (px)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(2,2,2)
imagesc(SkernVec,CorrVec,RMSEranNorm.strainY)
title('Random RMSE(strain_{yy})/max(strain_{yy})')
cx=colorbar;
cx.Label.String='%';
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(2,2,3)
imagesc(SkernVec,CorrVec,RMSEranNorm.strainS)
title('Random RMSE(strain_{xy})/max(strain_{xy})')
xlabel('SK (px)')
ylabel('Extrapolation Interval (px)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(2,2,4)
imagesc(SkernVec,CorrVec,RMSEranNorm.strainZ)
title('Random RMSE(strain_{zz})/max(strain_{zz})')
xlabel('SK (px)')
cx=colorbar;
colormap(Cmap.PBW_map)
cx.Label.String='%';
set(gca,'YDir','normal')
set(gca,'FontSize',12)

saveName=strcat(SweepDir,'/',TestDeg, ...
    '_SpatialCorrSmoothSweep_strainRMSEran_Norm');
saveas(gcf,saveName,'fig')
saveas(gcf,saveName,'png')

%% Plot Heat Maps of Total Normalized strain RMSE errors
figure('units','centimeters','InnerPosition',[1,1,18.2,18.2])
subplot(2,2,1)
imagesc(SkernVec,CorrVec,RMSEtotNorm.strainX)
title('Total RMSE(strain_{xx})/max(strain_{xx})')
ylabel('Extrapolation Interval (px)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(2,2,2)
imagesc(SkernVec,CorrVec,RMSEtotNorm.strainY)
title('Total RMSE(strain_{yy})/max(strain_{yy})')
cx=colorbar;
cx.Label.String='%';
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(2,2,3)
imagesc(SkernVec,CorrVec,RMSEtotNorm.strainS)
title('Total RMSE(strain_{xy})/max(strain_{xy})')
xlabel('SK (px)')
ylabel('Extrapolation Interval (px)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)

subplot(2,2,4)
imagesc(SkernVec,CorrVec,RMSEtotNorm.strainZ)
title('Total RMSE(strain_{zz})/max(strain_{zz})')
xlabel('SK (px)')
cx=colorbar;
colormap(Cmap.PBW_map)
cx.Label.String='%';
set(gca,'YDir','normal')
set(gca,'FontSize',12)

saveName=strcat(SweepDir,'/',TestDeg, ...
    '_SpatialCorrSmoothSweep_strainRMSEtot_Norm');
saveas(gcf,saveName,'fig')
saveas(gcf,saveName,'png')
%% Zoom IN
figure('units','centimeters','InnerPosition',[1,1,18.2,18.2])
subplot(2,2,1)
imagesc(SkernVec,CorrVec,RMSEsysNorm.strainX)
title('Systematic RMSE(strain_{xx})/max(strain_{xx})')
ylabel('Extrapolation Interval (px)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)
clim([min(RMSEsysNorm.strainX,[],'all'),10])

subplot(2,2,2)
imagesc(SkernVec,CorrVec,RMSEsysNorm.strainY)
title('Systematic RMSE(strain_{yy})/max(strain_{yy})')
cx=colorbar;
cx.Label.String='%';
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)
clim([min(RMSEsysNorm.strainY,[],'all'),10])

subplot(2,2,3)
imagesc(SkernVec,CorrVec,RMSEsysNorm.strainS)
title('Systematic RMSE(strain_{xy})/max(strain_{xy})')
xlabel('SK (px)')
ylabel('Extrapolation Interval (px)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)
clim([min(RMSEsysNorm.strainS,[],'all'),10])

subplot(2,2,4)
imagesc(SkernVec,CorrVec,RMSEsysNorm.strainZ)
title('Systematic RMSE(strain_{zz})/max(strain_{zz})')
xlabel('SK (px)')
cx=colorbar;
colormap(Cmap.PBW_map)
cx.Label.String='%';
set(gca,'YDir','normal')
set(gca,'FontSize',12)
clim([min(RMSEsysNorm.strainZ,[],'all'),10])


saveName=strcat(SweepDir,'/',TestDeg, ...
    '_SpatialCorrSmoothSweep_strainRMSEsys_NormZoom');
saveas(gcf,saveName,'fig')
saveas(gcf,saveName,'png')

%% Zoom IN
figure('units','centimeters','InnerPosition',[1,1,18.2,18.2])
subplot(2,2,1)
imagesc(SkernVec,CorrVec,RMSEranNorm.strainX)
title('Random RMSE(strain_{xx})/max(strain_{xx})')
ylabel('Extrapolation Interval (px)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)
clim([min(RMSEranNorm.strainX,[],'all'),1])

subplot(2,2,2)
imagesc(SkernVec,CorrVec,RMSEranNorm.strainY)
title('Random RMSE(strain_{yy})/max(strain_{yy})')
cx=colorbar;
cx.Label.String='%';
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)
clim([min(RMSEranNorm.strainY,[],'all'),1])

subplot(2,2,3)
imagesc(SkernVec,CorrVec,RMSEranNorm.strainS)
title('Random RMSE(strain_{xy})/max(strain_{xy})')
xlabel('SK (px)')
ylabel('Extrapolation Interval (px)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)
clim([min(RMSEranNorm.strainS,[],'all'),0.5])

subplot(2,2,4)
imagesc(SkernVec,CorrVec,RMSEranNorm.strainZ)
title('Random RMSE(strain_{zz})/max(strain_{zz})')
xlabel('SK (px)')
cx=colorbar;
colormap(Cmap.PBW_map)
cx.Label.String='%';
set(gca,'YDir','normal')
set(gca,'FontSize',12)
clim([min(RMSEranNorm.strainZ,[],'all'),1])


saveName=strcat(SweepDir,'/',TestDeg, ...
    '_SpatialCorrSmoothSweep_strainRMSEran_NormZoom');
saveas(gcf,saveName,'fig')
saveas(gcf,saveName,'png')

%% Zoom IN
figure('units','centimeters','InnerPosition',[1,1,18.2,18.2])
subplot(2,2,1)
imagesc(SkernVec,CorrVec,RMSEtotNorm.strainX)
title('Total RMSE(strain_{xx})/max(strain_{xx})')
ylabel('Extrapolation Interval (px)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)
clim([min(RMSEtotNorm.strainX,[],'all'),4])

subplot(2,2,2)
imagesc(SkernVec,CorrVec,RMSEtotNorm.strainY)
title('Total RMSE(strain_{yy})/max(strain_{yy})')
cx=colorbar;
cx.Label.String='%';
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)
clim([min(RMSEtotNorm.strainY,[],'all'),10])

subplot(2,2,3)
imagesc(SkernVec,CorrVec,RMSEtotNorm.strainS)
title('Total RMSE(strain_{xy})/max(strain_{xy})')
xlabel('SK (px)')
ylabel('Extrapolation Interval (px)')
colorbar;
colormap(Cmap.PBW_map)
set(gca,'YDir','normal')
set(gca,'FontSize',12)
clim([min(RMSEtotNorm.strainS,[],'all'),2.5])

subplot(2,2,4)
imagesc(SkernVec,CorrVec,RMSEtotNorm.strainZ)
title('Total RMSE(strain_{zz})/max(strain_{zz})')
xlabel('SK (px)')
cx=colorbar;
colormap(Cmap.PBW_map)
cx.Label.String='%';
set(gca,'YDir','normal')
set(gca,'FontSize',12)
clim([min(RMSEtotNorm.strainZ,[],'all'),6])


saveName=strcat(SweepDir,'/',TestDeg, ...
    '_SpatialCorrSmoothSweep_strainRMSEtot_NormZoom');
saveas(gcf,saveName,'fig')
saveas(gcf,saveName,'png')
%% Save results
saveName=strcat(SweepDir,'/',TestDeg, ...
    '_SpatialCorr_TempSpaSmooth_Sweep_results');
save(saveName)