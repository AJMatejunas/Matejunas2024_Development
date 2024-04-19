% This code is written to calculate stress gauge and constitutive model
    % stresses from images with and without shear correction and with exact
    % and identified constitutive parameter inputs


%Author: Andrew Matejunas

%Date: 2022-10-26





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialize
clear variables; close all; clc

%% load reference parameter file
[refPar.Name,refPar.Path]=uigetfile('*.mat', ...
    'Select file containing reference parameters');
refPar.File=strcat(refPar.Path,'/',refPar.Name);
refPar=load(refPar.File);

%% Define identified parameters
identPar=refPar;
identPar.MatProps.Gi=8.8737e+08;
identPar.MatProps.Ki=1.4135e9;
identPar.MatProps.tau=9.8378e-6;

%% Choose first image frame
[imageFile,imagePath] = uigetfile({'*.*','All Files'}, ...
    'Select the first image in the sequence');
%% Add Process function path
funcPath = uigetdir(pwd,'Locate Processing Function Folder');
addpath(funcPath);
addpath([funcPath,'GridMethodToolbox\']);

%% Find processing parameter file
[initFile,initPath,~] = uigetfile('*.mat', ...
    'Locate processing parameter file');
load([initPath,initFile])

%% Choose Directory to Save Results
MainSaveDir=uigetdir('','Choose Directory to Save Results');

%% Define Parent Designation
ParentDesig=char(inputdlg('Test Designation'));

%% Define smoothing opts for no smoothing
smoothingOpts.FFTempSmooth=0;
smoothingOpts.WATempSmooth=0;
smoothingOpts.spatialSmooth=0; 
imageNoise.addNoise=false;


%% Define Sample Condtioning Options
CondOpts.ImpCens=15;
CondOpts.FreeCens=20;
CondOpts.Xds=3;
CondOpts.Yds=1;
CondOpts.Tds=1;


if CondOpts.Tds <=1
    CondOpts.TempDS=false;
else
    CondOpts.TempDS=true;
end
   
%% Determine weather to run interpolated correction on the displacement 
    %fields
    
fprintf('Determining dispacement filed correction options \n')
% quest='Correct Disp fields?';
% dlgtitle='Displacement correction';
% DispCorr.Opt=questdlg(quest,dlgtitle,'Yes');
% clear quest dlgtitle
% 
%Displacement corrections are hardcoded
DispCorr.Opt='Yes';
DispCorr.int=10;
DispCorr.Method='LinGrad';
DispCorr.PitchFitKern=2;
DispCorr.strainMethod='GridPeriod';
DispCorr.strainPitchNum=2;


%% Run Grid method displacement measurement
fprintf('\n--------------------------------------------------------------\n')
fprintf('GRID METHOD PROCESSING\n')
fprintf('--------------------------------------------------------------\n')

%--------------------------------------------------------------------------
% GRID IMAGE PROCESSING

gridDataSavePath = imagePath;
gridDataFile = 'GridMethod_ProcessedData.mat';

fprintf('Checking for existing processed data file.\n')
processGridImages = true;
fprintf('Processing images using the grid method toolbox.\n')
% Process the image squence with the grid method toolbox
[grid,pos,disp] = func_gridMethodImageProcessing_AJM(imagePath,...
    imageFile,...
    grid,gridMethodOpts,imageNoise);

% Update Geometry and Number of Frames Based on Displacement Matrix Size
fprintf('Updating geometry variables based on the size of displacement field matrices.\n')
[specimen,grid] = func_updateSpecGeom(specimen,grid,disp);

% Create the time vector based on the number of frames in the disp struct
time.numFrames = size(disp.x,3);
time.vec = 0:time.step:(size(disp.x,3)-1)*time.step;

%Create arrays of x and y vectors
X_vec=pos.x;
Y_vec=pos.y;

%% POST-PROCESSING: Smoothing and Kinematic Field Derivation
% Smooth the displacement data and then calculate acceleration and strains
% Extrapolate the data to account for the missing pitch on the edges
fprintf('\n--------------------------------------------------------------\n')
fprintf('POST PROCESSING: Smoothing and Kinematic Field Calculation\n')
fprintf('--------------------------------------------------------------\n')

%% Perfom corrections if needed
switch    DispCorr.Opt
    case 'Yes'
        fprintf('Saving Raw Displacement Fields \n')
        RawDisp=disp;
        fprintf('Correcting Grid Method Displacements along specimen edges')
        [disp,DispCorr,grid,ProgramVersions]=func_CorrectGMDisp(disp,...
            DispCorr,grid,pos);

    case 'No'
        fprintf('No Displacement corrections performed \n')
end


% Calculate the kinematic fields from displacement fields using
% displacements from images or displacements from FE data

%--------------------------------------------------------------------------
% Load the Reference Image and Determine Where the Free Edge is
fprintf('Obtaining and setting the free edge location.\n')
[freeEdge,specimen,disp] = func_getFreeEdge(globalOpts.hardCodeFreeEdge,...
    imagePath,imageFile,specimen,disp);

% Smooth and Calculate Strain
fprintf('Calculating strain from the displacement fields.\n')
[disp,strain,strainRate] = func_smoothCalcStrain(globalOpts,pos,time,...
    grid,disp,smoothingOpts,extrapOpts);

% Smooth and Calculate Acceleration without strain correction
fprintf('Calculating acceleration from the displacement fields.\n')
[disp,~,accel] = func_smoothCalcAccel(pos,time,grid,disp,smoothingOpts,...
    extrapOpts,diffOpts);

% Remove some 3D fields from the structs to save memory.
    if globalOpts.reduceMemory
            disp.x = disp.extrap.x;
            disp.y = disp.extrap.y;
            disp = rmfield(disp,'extrap');
            disp = rmfield(disp,'tSmooth');
    end
%% Set Constitutive parameters to reference parameters
MatProps=refPar.MatProps;

%% Calculate the stress gauge stresses without strain correction
SG=func_Full_SG(accel,X_vec,time,MatProps.rho);

%% Calculate Constitutive Model Stresses
StressModel=func_ViscoConstitutiveV6(strain.x,strain.y,strain.s,time.vec, ...
    MatProps,0,0,0);

%% Save Data
saveName=strcat(ParentDesig, ...
    '_exactParam_NoShearCorr_StressStrainData.mat');
saveFile=strcat(MainSaveDir,'/',saveName);
save(saveFile,'disp','accel','pos','time','SG','X_vec',"Y_vec", ...
    "strainRate",'StressModel','strain','CondOpts','smoothingOpts', ...
    'DispCorr','MatProps','grid','imageNoise')

%% clear Unneccessary Data to Prepare for the next run
clear StressModel saveName saveFile MatProps 

%% Set Material Properties to identified constitutive parameters
MatProps=identPar.MatProps;

%% Calculate stress gauge stresses 
SG=func_Full_SG(accel,X_vec,time,MatProps.rho);

%% Calculate constitutive model stresses
StressModel=func_ViscoConstitutiveV6(strain.x,strain.y,strain.s,time.vec, ...
    MatProps,0,0,0);

%% Save Data
saveName=strcat(ParentDesig, ...
    '_identParam_NoShearCorr_StressStrainData.mat');
saveFile=strcat(MainSaveDir,'/',saveName);
save(saveFile,'disp','accel','pos','time','SG','X_vec',"Y_vec", ...
    "strainRate",'StressModel','strain','CondOpts','smoothingOpts', ...
    'DispCorr','MatProps','grid','imageNoise')

%% clear unneccessary data for next run
clear StressModel saveName saveFile MatProps disp accel strain SG 


%% Reprocess grid images to make sure there is no carryover of kinematic
    %fields
fprintf('\n--------------------------------------------------------------\n')
fprintf('GRID METHOD PROCESSING\n')
fprintf('--------------------------------------------------------------\n')

%--------------------------------------------------------------------------
% GRID IMAGE PROCESSING

gridDataSavePath = imagePath;
gridDataFile = 'GridMethod_ProcessedData.mat';

fprintf('Checking for existing processed data file.\n')
processGridImages = true;
fprintf('Processing images using the grid method toolbox.\n')
% Process the image squence with the grid method toolbox
[grid,pos,disp] = func_gridMethodImageProcessing_AJM(imagePath,...
    imageFile,...
    grid,gridMethodOpts,imageNoise);

% Update Geometry and Number of Frames Based on Displacement Matrix Size
fprintf('Updating geometry variables based on the size of displacement field matrices.\n')
[specimen,grid] = func_updateSpecGeom(specimen,grid,disp);

% Create the time vector based on the number of frames in the disp struct
time.numFrames = size(disp.x,3);
time.vec = 0:time.step:(size(disp.x,3)-1)*time.step;

%Create arrays of x and y vectors
X_vec=pos.x;
Y_vec=pos.y;

%% POST-PROCESSING: Smoothing and Kinematic Field Derivation
% Smooth the displacement data and then calculate acceleration and strains
% Extrapolate the data to account for the missing pitch on the edges
fprintf('\n--------------------------------------------------------------\n')
fprintf('POST PROCESSING: Smoothing and Kinematic Field Calculation\n')
fprintf('--------------------------------------------------------------\n')

%% Perfom corrections if needed
switch    DispCorr.Opt
    case 'Yes'
        fprintf('Saving Raw Displacement Fields \n')
        RawDisp=disp;
        fprintf('Correcting Grid Method Displacements along specimen edges')
        [disp,DispCorr,grid,ProgramVersions]=func_CorrectGMDisp(disp,...
            DispCorr,grid,pos);

    case 'No'
        fprintf('No Displacement corrections performed \n')
end


% Calculate the kinematic fields from displacement fields using
% displacements from images or displacements from FE data

%--------------------------------------------------------------------------
% Load the Reference Image and Determine Where the Free Edge is
fprintf('Obtaining and setting the free edge location.\n')
[freeEdge,specimen,disp] = func_getFreeEdge(globalOpts.hardCodeFreeEdge,...
    imagePath,imageFile,specimen,disp);

% Smooth and Calculate Strain
fprintf('Calculating strain from the displacement fields.\n')
[disp,strain,strainRate] = func_smoothCalcStrain(globalOpts,pos,time,...
    grid,disp,smoothingOpts,extrapOpts);
%Post Correct Strain
strain=func_PostCorrectShear(grid,strain,pos,time,DispCorr,...
    smoothingOpts);

% Smooth and Calculate Acceleration without strain correction
fprintf('Calculating acceleration from the displacement fields.\n')
[disp,~,accel] = func_smoothCalcAccel(pos,time,grid,disp,smoothingOpts,...
    extrapOpts,diffOpts);

% Remove some 3D fields from the structs to save memory.
    if globalOpts.reduceMemory
            disp.x = disp.extrap.x;
            disp.y = disp.extrap.y;
            disp = rmfield(disp,'extrap');
            disp = rmfield(disp,'tSmooth');
    end

%% Set Constitutive parameters to reference parameters
MatProps=refPar.MatProps;

%% Calculate the stress gauge stresses without strain correction
SG=func_Full_SG(accel,X_vec,time,MatProps.rho);

%% Calculate Constitutive Model Stresses
StressModel=func_ViscoConstitutiveV6(strain.x,strain.y,strain.s,time.vec, ...
    MatProps,0,0,0);

%% Save Data
saveName=strcat(ParentDesig, ...
    '_exactParam_ShearCorr_StressStrainData.mat');
saveFile=strcat(MainSaveDir,'/',saveName);
save(saveFile,'disp','accel','pos','time','SG','X_vec',"Y_vec", ...
    "strainRate",'StressModel','strain','CondOpts','smoothingOpts', ...
    'DispCorr','MatProps','grid','imageNoise')


%% clear Unneccessary Data to Prepare for the next run
clear StressModel saveName saveFile MatProps 

%% Set Material Properties to identified constitutive parameters
MatProps=identPar.MatProps;

%% Calculate stress gauge stresses 
SG=func_Full_SG(accel,X_vec,time,MatProps.rho);

%% Calculate constitutive model stresses
StressModel=func_ViscoConstitutiveV6(strain.x,strain.y,strain.s,time.vec, ...
    MatProps,0,0,0);

%% Save Data
saveName=strcat(ParentDesig, ...
    '_identParam_ShearCorr_StressStrainData.mat');
saveFile=strcat(MainSaveDir,'/',saveName);
save(saveFile,'disp','accel','pos','time','SG','X_vec',"Y_vec", ...
    "strainRate",'StressModel','strain','CondOpts','smoothingOpts', ...
    'DispCorr','MatProps','grid','imageNoise')
