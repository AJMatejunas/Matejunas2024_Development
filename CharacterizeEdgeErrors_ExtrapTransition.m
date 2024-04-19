%This script is written to quantify the error and produce error metric
    %plots for the corrected kinematic fields. 

%Author: Andrew Matejunas

%Date Created: 2023/03/06

%% initialize
clear variables; close all; clc

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
hardCodePath=0;
FEValidMode=0;
calcKinFieldsFromDisp=1;
saveGridData=0;

%% Add reference parameters file 
[RefPar.file,RefPar.path]=uigetfile('*.mat',...
    'Choose Mat file containing reference constitutive pamrameters');
load(strcat(RefPar.path,'/',RefPar.file),'MatProps');

%% Choose Directory to save results
savePath=uigetdir({},'Choose Folder to Save Results');

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
    
    fprintf('Loading reference image from the selected test data folder.\n')
    if ~hardCodePath
        [imageFile,imagePath] = uigetfile({'*.*','All Files'},'Select the first image in the sequence');
    else   
        imageFile = 'DefGridImage_001.tiff';
        imagePath = 'E:\Matlab_WorkingDirectory\1_IBIITest_Data\TestData_ID_CFIPQIso_HalfPulse\'; 
    end

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

%% Calculate the kinematic fields with grid method
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
[grid,pos,disp] = func_gridMethodImageProcessing_AJM(imagePath,...
    imageFile,...
    grid,gridMethodOpts,imageNoise);



%--------------------------------------------------------------------------
% Update Geometry and Number of Frames Based on Displacement Matrix Size
%fprintf('Updating geometry variables based on the size of displacement field matrices.\n')
[specimen,grid] = func_updateSpecGeom(specimen,grid,disp);

% Currently the rotations are unused so remove them to save RAM
%disp = rmfield(disp,'rot');


%--------------------------------------------------------------------------
% Create the time vector based on the number of frames in the disp struct
time.numFrames = size(disp.x,3);
time.vec = 0:time.step:(size(disp.x,3)-1)*time.step;

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

%% POST-PROCESSING: Smoothing and Kinematic Field Derivation
% Smooth the displacement data and then calculate acceleration and strains
% Extrapolate the data to account for the missing pitch on the edges
fprintf('\n--------------------------------------------------------------\n')
fprintf('POST PROCESSING: Smoothing and Kinematic Field Calculation\n')
fprintf('--------------------------------------------------------------\n')
%% Crop and Extrapolate displacement fields
[dispCorr.x,dispCorr.y,dispCorr.rX,dispCorr.rY]=...
    func_cropAndExtrapDispFields_v4(pos,disp.x,disp.y, ...
    extrapOpts.disp,true);


% Calculate the kinematic fields from displacement fields using
% displacements from images or displacements from FE data

%--------------------------------------------------------------------------
% Load the Reference Image and Determine Where the Free Edge is
%fprintf('Obtaining and setting the free edge location.\n')
[freeEdge,specimen,dispCorr] =...
    func_getFreeEdge(globalOpts.hardCodeFreeEdge,...
    imagePath,imageFile,specimen,dispCorr);

%--------------------------------------------------------------------------
% Smooth and Calculate Strain
fprintf('Calculating strain from the displacement fields.\n')

[strain,disp]=func_smoothCalcStrain_v4(pos,time,dispCorr,...
    smoothOpts.strain,extrapOpts.strain,true);



%--------------------------------------------------------------------------
% Smooth and Calculate Acceleration
fprintf('Calculating acceleration from the displacement fields.\n')
[accel,~,~] = func_smoothCalcAccel_v4(pos,time,dispCorr, ...
    smoothOpts.accel,...
    extrapOpts.accel,diffOpts,true);

% Remove some 3D fields from the structs to save memory.

%% Record spatially smoothed displacements
disp.x=disp.sSmooth.x;
disp.y=disp.sSmooth.y;
%% Crop disp back to oringinal size
disp.x = disp.x(dispCorr.rY,dispCorr.rX,:);
disp.y = disp.y(dispCorr.rY,dispCorr.rX,:);

%% Crop acceleration fields back to origninal size
accel.x = accel.x(dispCorr.rY,dispCorr.rX,:);
accel.y = accel.y(dispCorr.rY,dispCorr.rX,:);
%% Interpolate the FE data 
fprintf('Interpolating FE fields to grid coordinates \n')
interpMethod='nearest';
[Int.pos,Int.disp,Int.accel,Int.strain,Int.stress] =...
    func_interpFEt2Grid(FE.pos,pos,FE.disp,FE.accel,FE.strain,FE.stress, ...
    interpMethod);


%% Calculate Stress Gage Stresses

fprintf('Calculating Stress Gage Stresses \n')
SG=func_Full_SG(accel,X_vec,time,material.rho);
Int.SG=func_Full_SG(Int.accel,pos.x,time,material.rho);
%% Calculate constitutive model stresses
fprintf('Calculating Constitutive Model Stresses \n')
stressModel=func_ViscoConstitutiveV6(strain.x,strain.y,strain.s, ...
    time.vec,MatProps,0,0,0);
Int.StressModel=func_ViscoConstitutiveV6(Int.strain.x,Int.strain.y, ...
    Int.strain.s, ...
    time.vec,MatProps,0,0,0);

%% Calculate vertical averages for disp, accel strain, and stress fields
fprintf('Calculating Vertical averages of Kinematic Fields \n')
disp.xAvg=squeeze(mean(disp.x));
disp.yAvg=squeeze(mean(disp.y));
Int.disp.xAvg=squeeze(mean(Int.disp.x));
Int.disp.yAvg=squeeze(mean(Int.disp.y));

accel.xAvg=squeeze(mean(accel.x));
accel.yAvg=squeeze(mean(accel.y));
Int.accel.xAvg=squeeze(mean(Int.accel.x));
Int.accel.yAvg=squeeze(mean(Int.accel.y));

strain.xAvg=squeeze(mean(strain.x));
strain.yAvg=squeeze(mean(strain.y));
strain.sAvg=squeeze(mean(strain.s));
Int.strain.xAvg=squeeze(mean(Int.strain.x));
Int.strain.yAvg=squeeze(mean(Int.strain.y));
Int.strain.sAvg=squeeze(mean(Int.strain.s));


Int.stress.xAvg=squeeze(mean(Int.stress.x));
Int.stress.yAvg=squeeze(mean(Int.stress.y));
Int.stress.sAvg=squeeze(mean(Int.stress.s));

%% Calculate Errors
fprintf('Calculating Errors in the kinematic Fields\n')
Error.disp.x=(disp.x-Int.disp.x)./Int.disp.x*100;
Error.disp.y=(disp.y-Int.disp.y)./Int.disp.y*100;
Error.disp.xAvg=(disp.xAvg-Int.disp.xAvg)./Int.disp.xAvg*100;
Error.disp.yAvg=(disp.yAvg-Int.disp.yAvg)./Int.disp.yAvg*100;
Error.disp.x(abs(disp.x)<1e-9 & abs(Int.disp.x)<1e-9)=0;
Error.disp.y(abs(disp.y)<1e-9 & abs(Int.disp.y)<1e-9)=0;
Error.disp.xAvg(abs(disp.xAvg)<1e-9 & abs(Int.disp.xAvg)<1e-9)=0;
Error.disp.yAvg(abs(disp.yAvg)<1e-9 & abs(Int.disp.yAvg)<1e-9)=0;

Error.accel.x=(accel.x-Int.accel.x)./Int.accel.x*100;
Error.accel.y=(accel.y-Int.accel.y)./Int.accel.y*100;
Error.accel.xAvg=(accel.xAvg-Int.accel.xAvg)./Int.accel.xAvg*100;
Error.accel.yAvg=(accel.yAvg-Int.accel.yAvg)./Int.accel.yAvg*100;
% Error.accel.x(abs(accel.x)<1e-9 & abs(Int.accel.x)<1e-9)=0;
% Error.accel.y(abs(accel.y)<1e-9 & abs(Int.accel.y)<1e-9)=0;
% Error.accel.xAvg(abs(accel.xAvg)<1e-9 & abs(Int.accel.xAvg)<1e-9)=0;
% Error.accel.yAvg(abs(accel.yAvg)<1e-9 & abs(Int.accel.yAvg)<1e-9)=0;


Error.strain.x=(strain.x-Int.strain.x)./Int.strain.x*100;
Error.strain.y=(strain.y-Int.strain.y)./Int.strain.y*100;
Error.strain.s=(strain.s-Int.strain.s)./Int.strain.s*100;
Error.strain.xAvg=(strain.xAvg-Int.strain.xAvg)./Int.strain.xAvg*100;
Error.strain.yAvg=(strain.yAvg-Int.strain.yAvg)./Int.strain.yAvg*100;
Error.strain.sAvg=(strain.sAvg-Int.strain.sAvg)./Int.strain.sAvg*100;
Error.strain.x(abs(strain.x)<1e-6 & abs(Int.strain.x)<1e-6)=0;
Error.strain.y(abs(strain.y)<1e-6 & abs(Int.strain.y)<1e-6)=0;
Error.strain.s(abs(strain.s)<1e-6 & abs(Int.strain.s)<1e-6)=0;

Max.strain.xAvg=max(abs(Int.strain.xAvg),[],'all');
Error.strain.xAvg(abs(strain.xAvg)<.01*Max.strain.xAvg &...
    abs(Int.strain.xAvg)<.01*Max.strain.xAvg)=0;
Max.strain.yAvg=max(abs(Int.strain.yAvg),[],'all');
Error.strain.yAvg(abs(strain.yAvg)<.01*Max.strain.yAvg &...
    abs(Int.strain.yAvg)<.01*Max.strain.yAvg)=0;
Max.strain.sAvg=max(abs(Int.strain.sAvg),[],'all');
Error.strain.sAvg(abs(strain.sAvg)<.01*Max.strain.sAvg &...
    abs(Int.strain.sAvg)<.01*Max.strain.sAvg)=0;

Error.stress.x=(stressModel.xx-Int.stress.x)./Int.stress.x*100;
Error.stress.y=(stressModel.yy-Int.stress.y)./Int.stress.y*100;
Error.stress.s=(stressModel.xy-Int.stress.s)./Int.stress.s*100;
Error.stress.xAvg=(stressModel.Avxx-Int.stress.xAvg)./Int.stress.xAvg*100;
Error.stress.yAvg=(stressModel.Avyy-Int.stress.yAvg)./Int.stress.yAvg*100;
Error.stress.sAvg=(stressModel.Avxy-Int.stress.sAvg)./Int.stress.sAvg*100;
Error.stress.x(abs(stressModel.xx)<1e3 & abs(Int.stress.x)<1e3)=0;
Error.stress.y(abs(stressModel.yy)<1e3 & abs(Int.stress.y)<1e3)=0;
Error.stress.s(abs(stressModel.xy)<1e3 & abs(Int.stress.s)<1e3)=0;

Med.stress.xAvg=median(abs(Int.stress.xAvg),'all');
Med.stress.sAvg=median(abs(Int.stress.sAvg),'all');
Error.stress.xAvg(abs(stressModel.Avxx)<0.01*Med.stress.xAvg &...
    abs(Int.stress.xAvg)<0.01*Med.stress.xAvg)=0;
Error.stress.sAvg(abs(stressModel.Avxy)<0.01*Med.stress.sAvg & ...
    abs(Int.stress.sAvg)<0.01*Med.stress.sAvg)=0;
%% clear Raw FE data to save 
clear FE

Xmm=pos.x*10^3;
timems=time.vec*10^6;
ErrLim.stress.xAvg=10*[-1,1];
ErrLim.stress.sAvg=10*[-1,1];


ConstDir=strcat(savePath,'/ConstError');
mkdir(ConstDir)

figure('units','normalized','OuterPosition',[0,0,1,1])
subplot(2,1,1)
imagesc(Xmm,timems,Error.stress.xAvg')
xlabel('x_0 (mm)')
ylabel('t (\mus)')
title('\sigma_x Error')
cx=colorbar;
cx.Label.String='$\overline{\sigma_x}^y$ Error';
cx.Label.Interpreter='latex';
colormap(Cmap.PBWYRK_map)
caxis (ErrLim.stress.xAvg);
set(gca,'YDir','normal')

subplot(2,1,2)
imagesc(Xmm,timems,Error.stress.sAvg')
title('\sigma_s Error')
xlabel('x_0 (mm)')
ylabel('t (\mus)')
cx=colorbar;
cx.Label.String='$\overline{\sigma_y}^y$ Error';
cx.Label.Interpreter='latex';
colormap(Cmap.PBWYRK_map)
caxis (ErrLim.stress.sAvg);
set(gca,'YDir','normal')

figName=(strcat(ConstDir,'/',TestDeg,'_ConstErrorXt'));
saveas(gcf,figName,'png')
saveas(gcf,figName,'fig')
%% Verify interpolation
figure('units','normalized','OuterPosition',[0,0,1,1])
Xmm=pos.x*10^3;
Ymm=pos.y*10^3;
intDir=strcat(savePath,'/IntVer/disp');
mkdir(intDir);
ErrLim.dispY=[-5,5];
DyLim=max(abs(Int.disp.y),[],'all')*[-1,1];
Error.disp.x(abs(disp.x)<1e-9 & abs(Int.disp.x)<1e-9)=0;
Error.disp.y(abs(disp.y)<1e-9 & abs(Int.disp.y)<1e-9)=0;

for t=1:length(time.vec)
    fnum=num2str(t);
    IDy=squeeze(Int.disp.y(:,:,t));
    Dy=squeeze(disp.y(:,:,t));
    ErrDy=squeeze(Error.disp.y(:,:,t));
    
    
    subplot(3,1,1)
    imagesc(Xmm,Ymm,IDy)
    xlabel('x')
    ylabel('y')
    colormap(Cmap.PBWYRK_map)
    caxis (DyLim)
    title(strcat('Interpolated disp_y Frame~',fnum))
    cx=colorbar;
    cx.Label.String='u_y';

    subplot(3,1,2)
    imagesc(Xmm,Ymm,Dy)
    xlabel('x')
    ylabel('y')
    colormap(Cmap.PBWYRK_map)
    caxis (DyLim)
    title(strcat('GM disp_y Frame~',fnum))
    colorbar

    subplot(3,1,3)
    imagesc(Xmm,Ymm,ErrDy)
    xlabel('x')
    ylabel('y')
    colormap(Cmap.PBWYRK_map)
    caxis (ErrLim.dispY)
    title(strcat('Error disp_y Frame~',fnum))
    cx=colorbar;
    cx.Label.String='Error (%)';

    saveas(gcf,strcat(intDir,'/',TestDeg,'_DispIntVer_Frame',fnum),'png')
end
%% Verify strain_yy interpolation
figure('units','normalized','OuterPosition',[0,0,1,1])
Xmm=pos.x*10^3;
Ymm=pos.y*10^3;
intDir=strcat(savePath,'/IntVer/strain');
mkdir(intDir);
ErrLim.strainY=[-5,5];
StrainYLim=max(abs(Int.strain.y),[],'all')*[-1,1];
Error.strain.x(abs(strain.x)<1e-9 & abs(Int.strain.x)<1e-9)=0;
Error.strain.y(abs(strain.y)<1e-9 & abs(Int.strain.y)<1e-9)=0;

for t=1:length(time.vec)
    fnum=num2str(t);
    IStrainY=squeeze(Int.strain.y(:,:,t));
    StrainY=squeeze(strain.y(:,:,t));
    ErrStrainY=squeeze(Error.strain.y(:,:,t));
    
    
    subplot(3,1,1)
    imagesc(Xmm,Ymm,IStrainY)
    xlabel('x')
    ylabel('y')
    colormap(Cmap.PBWYRK_map)
    caxis (StrainYLim)
    title(strcat('Interpolated strain_y Frame~',fnum))
    cx=colorbar;
    cx.Label.String='u_y';

    subplot(3,1,2)
    imagesc(Xmm,Ymm,StrainY)
    xlabel('x')
    ylabel('y')
    colormap(Cmap.PBWYRK_map)
    caxis (StrainYLim)
    title(strcat('GM strain_y Frame~',fnum))
    colorbar

    subplot(3,1,3)
    imagesc(Xmm,Ymm,ErrStrainY)
    xlabel('x')
    ylabel('y')
    colormap(Cmap.PBWYRK_map)
    caxis (ErrLim.strainY)
    title(strcat('Error strain_y Frame~',fnum))
    cx=colorbar;
    cx.Label.String='Error (%)';

    saveas(gcf,strcat(intDir,'/',TestDeg,'_StrainYYIntVer_Frame',fnum),'png')
end
%% Verify strain_xx interpolation
figure('units','normalized','OuterPosition',[0,0,1,1])
Xmm=pos.x*10^3;
Ymm=pos.y*10^3;
intDir=strcat(savePath,'/IntVer/strainXX');
mkdir(intDir);
ErrLim.strainX=[-5,5];
StrainXLim=max(abs(Int.strain.x),[],'all')*[-1,1];
Error.strain.x(abs(strain.x)<1e-6 & abs(Int.strain.x)<1e-6)=0;
Error.strain.y(abs(strain.y)<1e-6 & abs(Int.strain.y)<1e-6)=0;

for t=1:length(time.vec)
    fnum=num2str(t);
    IStrainX=squeeze(Int.strain.x(:,:,t));
    StrainX=squeeze(strain.x(:,:,t));
    ErrStrainX=squeeze(Error.strain.x(:,:,t));
    
    
    subplot(3,1,1)
    imagesc(Xmm,Ymm,IStrainX)
    xlabel('x')
    ylabel('y')
    colormap(Cmap.PBWYRK_map)
    caxis (StrainXLim)
    title(strcat('Interpolated strain_x Frame~',fnum))
    cx=colorbar;
    cx.Label.String='u_y';

    subplot(3,1,2)
    imagesc(Xmm,Ymm,StrainX)
    xlabel('x')
    ylabel('y')
    colormap(Cmap.PBWYRK_map)
    caxis (StrainXLim)
    title(strcat('GM strain_x Frame~',fnum))
    colorbar

    subplot(3,1,3)
    imagesc(Xmm,Ymm,ErrStrainX)
    xlabel('x')
    ylabel('y')
    colormap(Cmap.PBWYRK_map)
    caxis (ErrLim.strainX)
    title(strcat('Error strain_x Frame~',fnum))
    cx=colorbar;
    cx.Label.String='Error (%)';

    saveas(gcf,strcat(intDir,'/',TestDeg,'_StrainxxIntVer_Frame',fnum),'png')
end

%% Verify Shear interpolation
figure('units','normalized','OuterPosition',[0,0,1,1])
Xmm=pos.x*10^3;
Ymm=pos.y*10^3;
intDir=strcat(savePath,'/IntVer/strainXY');
mkdir(intDir);
ErrLim.strainS=[-5,5];
StrainSLim=max(abs(Int.strain.s),[],'all')*[-1,1];
Error.strain.s(abs(strain.s)<1e-6 & abs(Int.strain.s)<1e-6)=0;


for t=1:length(time.vec)
    fnum=num2str(t);
    IStrainS=squeeze(Int.strain.s(:,:,t));
    StrainS=squeeze(strain.s(:,:,t));
    ErrStrainS=squeeze(Error.strain.s(:,:,t));
    
    
    subplot(3,1,1)
    imagesc(Xmm,Ymm,IStrainS)
    xlabel('x')
    ylabel('y')
    colormap(Cmap.PBWYRK_map)
    caxis (StrainXLim)
    title(strcat('Interpolated strain_s Frame~',fnum))
    cx=colorbar;
   

    subplot(3,1,2)
    imagesc(Xmm,Ymm,StrainS)
    xlabel('x')
    ylabel('y')
    colormap(Cmap.PBWYRK_map)
    caxis (StrainSLim)
    title(strcat('GM strain_s Frame~',fnum))
    colorbar

    subplot(3,1,3)
    imagesc(Xmm,Ymm,ErrStrainS)
    xlabel('x')
    ylabel('y')
    colormap(Cmap.PBWYRK_map)
    caxis (ErrLim.strainS)
    title(strcat('Error strain_s Frame~',fnum))
    cx=colorbar;
    cx.Label.String='Error (%)';

    saveas(gcf,strcat(intDir,'/',TestDeg,'_StrainSIntVer_Frame',fnum),'png')
end

%% Displacement Error  Heat Maps
Xmm=pos.x*10^3;
timems=time.vec*10^6;
ErrLim.disp.xAvg=10*[-1,1];
ErrLim.disp.yAvg=10*[-1,1];


Error.disp.xAvg(abs(disp.xAvg)<1e-9 & abs(Int.disp.xAvg)<1e-9)=0;
Error.disp.yAvg(abs(disp.xAvg)<1e-9 & abs(Int.disp.yAvg)<1e-9)=0;

dispDir=strcat(savePath,'/dispError');
mkdir(dispDir)

figure('units','normalized','OuterPosition',[0,0,1,1])
subplot(2,1,1)
imagesc(Xmm,timems,Error.disp.xAvg')
xlabel('x_0 (mm)')
ylabel('t (\mus)')
cx=colorbar;
cx.Label.String='$\overline{u_x}^y$ Error';
cx.Label.Interpreter='latex';
colormap(Cmap.PBWYRK_map)
caxis (ErrLim.disp.xAvg);
set(gca,'YDir','normal')

subplot(2,1,2)
imagesc(Xmm,timems,Error.disp.yAvg')
xlabel('x_0 (mm)')
ylabel('t (\mus)')
cx=colorbar;
cx.Label.String='$\overline{u_y}^y$ Error';
cx.Label.Interpreter='latex';
colormap(Cmap.PBWYRK_map)
caxis (ErrLim.disp.yAvg);
set(gca,'YDir','normal')

figName=(strcat(dispDir,'/',TestDeg,'_dispErrorXt'));
saveas(gcf,figName,'png')
saveas(gcf,figName,'fig')


%% Get indexes for middle and 4 grid pitches from free and impact edges

GMPos.NumPoints=size(disp.y,2);
GMPos.Free=20; %4 pitches from free edge
GMPos.Mid=round(GMPos.NumPoints/2);
GMPos.Imp=GMPos.NumPoints-20; %4 pitches from impact edge
% Get cooridinates of position indexes
GMPos.Xfree=pos.x(GMPos.Free);
GMPos.Xmid=pos.x(GMPos.Mid);
GMPos.Ximp=pos.x(GMPos.Imp);


%% plot Y displacement curves
Error.disp.x(abs(disp.x)<1e-9 & abs(Int.disp.x)<1e-9)=0;
Error.disp.y(abs(disp.y)<1e-9 & abs(Int.disp.y)<1e-9)=0;

MovDirFig=strcat(savePath,'/dispYMov/fig');
MovDirPng=strcat(savePath,'/dispYMov/png');
mkdir(MovDirPng);
mkdir(MovDirFig);

figure('Units','normalized','OuterPosition',[0,0,1,1])
timems=time.vec*10^6;
Ymm=pos.y*10^3;
% Set plot limits
YdispLim=1.1*[min(Int.disp.y,[],'all'),max(Int.disp.y,[],'all')];
YErrLim=[-10,10];
for k=1:length(time.vec)
    Fnum=num2str(k);
    Fstring=strcat('Frame~',Fnum);

    % Extract vector
    %FE Model Displacement
    FEDyImp=squeeze(Int.disp.y(:,GMPos.Imp,k));
    FEDyMid=squeeze(Int.disp.y(:,GMPos.Mid,k));
    FEDyFree=squeeze(Int.disp.y(:,GMPos.Free,k));

    %Uncorrected Grid Method Displacement
    DyImp=squeeze(disp.y(:,GMPos.Imp,k));
    DyMid=squeeze(disp.y(:,GMPos.Mid,k));
    DyFree=squeeze(disp.y(:,GMPos.Free,k));

    ErrYImp=squeeze(Error.disp.y(:,GMPos.Imp,k));
    ErrYMid=squeeze(Error.disp.y(:,GMPos.Mid,k));
    ErrYFree=squeeze(Error.disp.y(:,GMPos.Free,k));

    %  Y Displacement at the impact edge
    subplot(3,1,1)
    yyaxis left
    plot(Ymm,FEDyImp,'k')
    hold on
    plot(Ymm,DyImp,'--b')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('u_y (m)')
    ylim(YdispLim)
    yyaxis right
    plot(Ymm,ErrYImp,'r')
    ylabel('Error (%)')
    ylim(YErrLim)
    title(strcat('4P from Impact ',Fstring))
    legend('FE','RawGM','Error','location','southwest')


    % Y Displacement in middle
    subplot(3,1,2)
    yyaxis left
    plot(Ymm,FEDyMid,'k')
    hold on
    plot(Ymm,DyMid,'--b')
    hold off
    ylim(YdispLim)
    xlabel('Y Coordinate (mm)')
    ylabel('u_y (m)')
    yyaxis right
    plot(Ymm,ErrYMid,'r')
    ylim(YErrLim)
    title(strcat('Middle ',Fstring))


    % Y displacement Free Surface
    subplot(3,1,3)
    plot(Ymm,FEDyFree,'k')
    hold on
    plot(Ymm,DyFree,'--b')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('u_y (m)')
    ylim(YdispLim)
    
    yyaxis right
    plot(Ymm,ErrYFree,'r')
    ylabel('Error (%)')
    ylim(YErrLim)
    title(strcat('4P from Free Edge ',Fstring))
    % SaveFiles
    FigName=strcat(MovDirFig,'/',TestDeg,'_1DyD_',Fnum);
    saveas(gcf,FigName,'fig');
   PngName=strcat(MovDirPng,'/',TestDeg,'_1DyD_',Fnum);
    saveas(gcf,PngName,'fig');
end

%% plot x displacement curves
Error.disp.x(abs(disp.x)<1e-9 & abs(Int.disp.x)<1e-9)=0;
Error.disp.y(abs(disp.y)<1e-9 & abs(Int.disp.y)<1e-9)=0;

MovDirFig=strcat(savePath,'/dispXMov/fig');
MovDirPng=strcat(savePath,'/dispXMov/png');
mkdir(MovDirPng);
mkdir(MovDirFig);

figure('Units','normalized','OuterPosition',[0,0,1,1])
timems=time.vec*10^6;

Ymm=pos.y*10^3;
% Set plot limits
XdispLim=1.1*[min(Int.disp.x,[],'all'),max(Int.disp.x,[],'all')];
XErrLim=[-10,10];
for k=1:length(time.vec)
    Fnum=num2str(k);
    Fstring=strcat('Frame~',Fnum);

    % Extract vector
    %FE Model Displacement
    FEDxImp=squeeze(Int.disp.x(:,GMPos.Imp,k));
    FEDxMid=squeeze(Int.disp.x(:,GMPos.Mid,k));
    FEDxFree=squeeze(Int.disp.x(:,GMPos.Free,k));

    %Uncorrected Grid Method Displacement
    DxImp=squeeze(disp.x(:,GMPos.Imp,k));
    DxMid=squeeze(disp.x(:,GMPos.Mid,k));
    DxFree=squeeze(disp.x(:,GMPos.Free,k));

    ErrXImp=squeeze(Error.disp.x(:,GMPos.Imp,k));
    ErrXMid=squeeze(Error.disp.x(:,GMPos.Mid,k));
    ErrXFree=squeeze(Error.disp.x(:,GMPos.Free,k));

    %  X Displacement at the impact edge
    subplot(3,1,1)
    yyaxis left
    plot(Ymm,FEDxImp,'k')
    hold on
    plot(Ymm,DxImp,'--b')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('u_x (m)')
    ylim(XdispLim)
    yyaxis right
    plot(Ymm,ErrXImp,'r')
    ylabel('Error (%)')
    ylim(XErrLim)
    title(strcat('4P from Impact ',Fstring))
    legend('FE','RawGM','Error','location','southwest')


    % X Displacement in middle
    subplot(3,1,2)
    yyaxis left
    plot(Ymm,FEDxMid,'k')
    hold on
    plot(Ymm,DxMid,'--b')
    hold off
    ylim(XdispLim)
    xlabel('Y Coordinate (mm)')
    ylabel('u_x (m)')
    yyaxis right
    plot(Ymm,ErrXMid,'r')
    ylim(XErrLim)
    title(strcat('Middle ',Fstring))


    % X displacement Free Surface
    subplot(3,1,3)
    plot(Ymm,FEDxFree,'k')
    hold on
    plot(Ymm,DxFree,'--b')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('u_x (m)')
    ylim(XdispLim)
    
    yyaxis right
    plot(Ymm,ErrXFree,'r')
    ylabel('Error (%)')
    ylim(XErrLim)
    title(strcat('4P from Free Edge ',Fstring))
    % SaveFiles
    FigName=strcat(MovDirFig,'/',TestDeg,'_1DxD_',Fnum);
    saveas(gcf,FigName,'fig');
   PngName=strcat(MovDirPng,'/',TestDeg,'_1DxD_',Fnum);
    saveas(gcf,PngName,'fig');
end

%% Acceleration Error  Heat Maps
Xmm=pos.x*10^3;
timems=time.vec*10^6;
ErrLim.accel.xAvg=10*[-1,1];
ErrLim.accel.yAvg=10*[-1,1];
Med.accel.xAvg=median(abs(accel.xAvg),'all');
Med.accel.yAvg=median(abs(accel.yAvg),'all');

Error.accel.xAvg(abs(accel.xAvg)<.001*Med.accel.xAvg ...
    & abs(Int.accel.xAvg)<Med.accel.xAvg)=0;
Error.accel.yAvg(abs(accel.xAvg)<Med.accel.yAvg ...
    & abs(Int.accel.yAvg)<Med.accel.yAvg)=0;

accelDir=strcat(savePath,'/accelError');
mkdir(accelDir)

figure('units','normalized','OuterPosition',[0,0,1,1])
subplot(2,1,1)
imagesc(Xmm,timems,Error.accel.xAvg')
xlabel('x_0 (mm)')
ylabel('t (\mus)')
title('a_x Error')
cx=colorbar;
cx.Label.String='$\overline{a_x}^y$ Error';
cx.Label.Interpreter='latex';
colormap(Cmap.PBWYRK_map)
caxis (ErrLim.accel.xAvg);
set(gca,'YDir','normal')

subplot(2,1,2)
imagesc(Xmm,timems,Error.accel.yAvg')
title('a_y Error')
xlabel('x_0 (mm)')
ylabel('t (\mus)')
cx=colorbar;
cx.Label.String='$\overline{a_y}^y$ Error';
cx.Label.Interpreter='latex';
colormap(Cmap.PBWYRK_map)
caxis (ErrLim.accel.yAvg);
set(gca,'YDir','normal')

figName=(strcat(accelDir,'/',TestDeg,'_accelErrorXt'));
saveas(gcf,figName,'png')
saveas(gcf,figName,'fig')

%% Calculate surface average acceleration errors
Error.accel.xSurf=(SG.Ax_surfx-Int.SG.Ax_surfx)./Int.SG.Ax_surfx*100;
Error.accel.ySurf=(SG.AY_surfx-Int.SG.AY_surfx)./Int.SG.AY_surfx*100;
Med.accel.xSurf=median(abs(Int.SG.Ax_surfx),'all');
Med.accel.ySurf=median(abs(Int.SG.AY_surfx),'all');
Error.accel.xSurf(abs(SG.Ax_surfx)<.001*Med.accel.xSurf ...
    & abs(Int.SG.Ax_surfx)<Med.accel.xSurf)=0;
Error.accel.ySurf(abs(SG.AY_surfx)<Med.accel.ySurf ...
    & abs(Int.SG.AY_surfx)<Med.accel.ySurf)=0;

%% Surface Acceleration Error  Heat Maps
Xmm=pos.x*10^3;
timems=time.vec*10^6;
ErrLim.accel.xSurf=10*[-1,1];
ErrLim.accel.ySurf=10*[-1,1];


accelDir=strcat(savePath,'/SurfaccelError');
mkdir(accelDir)

figure('units','normalized','OuterPosition',[0,0,1,1])
subplot(2,1,1)
imagesc(Xmm,timems,Error.accel.xSurf')
xlabel('x_0 (mm)')
ylabel('t (\mus)')
title('Surface a_x Error')
cx=colorbar;
cx.Label.String='$\overline{a_x}^S$ Error';
cx.Label.Interpreter='latex';
colormap(Cmap.PBWYRK_map)
caxis (ErrLim.accel.xSurf);
set(gca,'YDir','normal')

subplot(2,1,2)
imagesc(Xmm,timems,Error.accel.ySurf')
title('Surface a_y Error')
xlabel('x_0 (mm)')
ylabel('t (\mus)')
cx=colorbar;
cx.Label.String='$\overline{a_y}^S$ Error';
cx.Label.Interpreter='latex';
colormap(Cmap.PBWYRK_map)
caxis (ErrLim.accel.ySurf);
set(gca,'YDir','normal')

figName=(strcat(accelDir,'/',TestDeg,'_SurfaceAccelErrorXt'));
saveas(gcf,figName,'png')
saveas(gcf,figName,'fig')


%% Calculate Strain33 Errors
Error.strain.z=(stressModel.zzstrain-Int.StressModel.zzstrain)./...
    Int.StressModel.zzstrain*100;
Error.strain.zAvg=(stressModel.ZZAvstrain-Int.StressModel.ZZAvstrain)./...
    Int.StressModel.ZZAvstrain*100;
Max.strain.zAvg=max(abs(Int.StressModel.ZZAvstrain),[],'all');
Error.strain.zAvg(abs(stressModel.ZZAvstrain)<.01*Max.strain.zAvg &...
    abs(Int.StressModel.ZZAvstrain)<.01*Max.strain.zAvg)=0;

%% Strain Error Heat Maps
Xmm=pos.x*10^3;
timems=time.vec*10^6;
ErrLim.strain.xAvg=5*[-1,1];
ErrLim.strain.yAvg=10*[-1,1];
ErrLim.strain.sAvg=5*[-1,1];
ErrLim.strain.zAvg=10*[-1,1];

strainDir=strcat(savePath,'/AvStrainError');
mkdir(accelDir)

figure('units','normalized','OuterPosition',[0,0,1,1])
subplot(2,2,1)
imagesc(Xmm,timems,Error.strain.xAvg')
xlabel('x_0 (mm)')
ylabel('t (\mus)')
title('strain_x Error')
cx=colorbar;
cx.Label.String='$\overline{\varepsilon_x}^y$ Error';
cx.Label.Interpreter='latex';
colormap(Cmap.PBWYRK_map)
caxis (ErrLim.strain.xAvg);
set(gca,'YDir','normal')

subplot(2,2,2)
imagesc(Xmm,timems,Error.strain.yAvg')
title('strain_y Error')
xlabel('x_0 (mm)')
ylabel('t (\mus)')
cx=colorbar;
cx.Label.String='$\overline{\varepsilon_y}^y$ Error';
cx.Label.Interpreter='latex';
colormap(Cmap.PBWYRK_map)
caxis (ErrLim.strain.yAvg);
set(gca,'YDir','normal')

subplot(2,2,3)
imagesc(Xmm,timems,Error.strain.sAvg')
xlabel('x_0 (mm)')
ylabel('t (\mus)')
title('strain_s Error')
cx=colorbar;
cx.Label.String='$\overline{\varepsilon_s}^y$ Error';
cx.Label.Interpreter='latex';
colormap(Cmap.PBWYRK_map)
caxis (ErrLim.strain.sAvg);
set(gca,'YDir','normal')

subplot(2,2,4)
imagesc(Xmm,timems,Error.strain.zAvg')
title('strain_z Error')
xlabel('x_0 (mm)')
ylabel('t (\mus)')
cx=colorbar;
cx.Label.String='$\overline{\varepsilon_z}^y$ Error';
cx.Label.Interpreter='latex';
colormap(Cmap.PBWYRK_map)
caxis (ErrLim.strain.zAvg);
set(gca,'YDir','normal')


figName=(strcat(accelDir,'/',TestDeg,'_SurfaceAccelErrorXt'));
saveas(gcf,figName,'png')
saveas(gcf,figName,'fig')

%% Constitutive Model Error
Xmm=pos.x*10^3;
timems=time.vec*10^6;
ErrLim.stress.xAvg=10*[-1,1];
ErrLim.stress.sAvg=10*[-1,1];


ConstDir=strcat(savePath,'/ConstError');
mkdir(ConstDir)

figure('units','normalized','OuterPosition',[0,0,1,1])
subplot(2,1,1)
imagesc(Xmm,timems,Error.stress.xAvg')
xlabel('x_0 (mm)')
ylabel('t (\mus)')
title('\sigma_x Error')
cx=colorbar;
cx.Label.String='$\overline{\sigma_x}^y$ Error';
cx.Label.Interpreter='latex';
colormap(Cmap.PBWYRK_map)
caxis (ErrLim.stress.xAvg);
set(gca,'YDir','normal')

subplot(2,1,2)
imagesc(Xmm,timems,Error.stress.sAvg')
title('\sigma_s Error')
xlabel('x_0 (mm)')
ylabel('t (\mus)')
cx=colorbar;
cx.Label.String='$\overline{\sigma_s}^y$ Error';
cx.Label.Interpreter='latex';
colormap(Cmap.PBWYRK_map)
caxis (ErrLim.stress.sAvg);
set(gca,'YDir','normal')

figName=(strcat(ConstDir,'/',TestDeg,'_ConstErrorXt'));
saveas(gcf,figName,'png')
saveas(gcf,figName,'fig')
%% Calculate Stress Gauge Error
Error.SG.x=(SG.x-Int.stress.xAvg)./Int.stress.xAvg*100;
Error.SG.s=(SG.s-Int.stress.sAvg)./Int.stress.sAvg*100;
Error.SG.x(abs(SG.x)<0.01*Med.stress.xAvg & ...
    abs(Int.stress.xAvg)<0.01*Med.stress.xAvg)=0;
Error.SG.s(abs(SG.s)<0.01*Med.stress.sAvg & ...
    abs(Int.stress.sAvg)<0.01*Med.stress.sAvg)=0;
%% SG error
Xmm=pos.x*10^3;
timems=time.vec*10^6;
ErrLim.stress.xAvg=10*[-1,1];
ErrLim.stress.sAvg=10*[-1,1];


SGDir=strcat(savePath,'/SGError');
mkdir(SGDir)

figure('units','normalized','OuterPosition',[0,0,1,1])
subplot(2,1,1)
imagesc(Xmm,timems,Error.SG.x')
xlabel('x_0 (mm)')
ylabel('t (\mus)')
title('SG_x Error')
cx=colorbar;
cx.Label.String='$\overline{SG_x}^y$ Error';
cx.Label.Interpreter='latex';
colormap(Cmap.PBWYRK_map)
caxis (ErrLim.stress.xAvg);
set(gca,'YDir','normal')

subplot(2,1,2)
imagesc(Xmm,timems,Error.SG.s')
title('SG_s Error')
xlabel('x_0 (mm)')
ylabel('t (\mus)')
cx=colorbar;
cx.Label.String='$\overline{SG_s}^y$ Error';
cx.Label.Interpreter='latex';
colormap(Cmap.PBWYRK_map)
caxis (ErrLim.stress.sAvg);
set(gca,'YDir','normal')

figName=(strcat(SGDir,'/',TestDeg,'_SGErrorXt'));
saveas(gcf,figName,'png')
saveas(gcf,figName,'fig')
