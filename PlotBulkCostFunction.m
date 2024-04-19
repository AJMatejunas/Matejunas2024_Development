%This script is written to plot the cost funciton in the bulk modulus for a
    %large range of identified parameters and certain inputs of tau and G

 %Author: Andrew Matejunas

 %Date: 2022-10-13

 %Change Log:


 %% initialize
clear vars; close all; clc

%% Define test designation and save path
MainSaveDir=uigetdir({},'Choose Folder to Save Results');
ParentDesig=char(cell2mat(inputdlg('Input Test Designation')));

%% Choose constitutive parameter input file
[RefPar.Name,RefPar.Path]=uigetfile('*.mat',...
    'Choose file containing reference parameters');
RefPar.File=strcat(RefPar.Path,'/',RefPar.Name);
RefPar=load(RefPar.File);

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

        %--------------------------------------------------------------------------
        % Smooth and Calculate Strain
        fprintf('Calculating strain from the displacement fields.\n')
        [disp,strain,strainRate] = func_smoothCalcStrain(globalOpts,pos,time,...
            grid,disp,smoothingOpts,extrapOpts);
        %Post Correct Strain
        strain=func_PostCorrectShear(grid,strain,pos,time,DispCorr,...
            smoothingOpts);
        %--------------------------------------------------------------------------
        % Smooth and Calculate Acceleration
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
        %     else
        %         fprintf('Using kinematic fields directly from FE data.\n')
        %         % For pure FE data the kinematic fields can be taken directly from the
        %         % FE calculation
        %         disp.xAvg = func_avgFFVarOverWidth(disp.x);
        %         accel.xAvg = func_avgFFVarOverWidth(accel.x);
        %         strain.xAvg = func_avgFFVarOverWidth(strain.x);
        %         strain.yAvg = func_avgFFVarOverWidth(strain.y);
        %         [strainRate,~] = func_calculateStrainRate(strain,time,smoothingOpts,true);
        %         strainRate.xAvg = func_avgFFVarOverWidth(strainRate.x);
        %         strainRate = func_calcMaxStrainRate(smoothingOpts,grid,strainRate);
        %         disp.extrap.x = disp.x;
        %         disp.extrap.y = disp.y;
  

 
%% Create Vectors of Material Parameters

%Upper and lower limit of bulk modulus. 
Kul=RefPar.MatProps.Ki*1.2;
Kll=RefPar.MatProps.Ki*.8;
% Vector of K to plot cost function for
K_vec=linspace(Kll,Kul,100);

%Upper and lower limits for shear modulus
Gul=RefPar.MatProps.Gi*1.1;
Gll=RefPar.MatProps.Gi*.9;
G_vec=linspace(Gll,Gul,100);

%Identified shear modulus with no smoothing
Gident=8.8737e+08;

%Exact shear modulus input
Gexact=8.7698e8;

%Upper and lower limits for tau
tau_ul=10e-6*1.1;
tau_ll=10e-6*.9;

tau_vec=linspace(tau_ll,tau_ul,100);
tau_ident=9.8378e-6;
tau_exact=10e-6;

%% Calculate Stress Gage Stresses
 fprintf('Calculating Stress Gage Stresses \n')
 SG=func_Full_SG(accel,X_vec,time,RefPar.MatProps.rho);
%% Plot Heat Maps of the cost function with identified inputs
costDesig=strcat(ParentDesig,'_identParam');
[phi_ctau_Ident,phi_cG_ident,phi3D_ident] = func_PlotBulkCost(K_vec,G_vec,tau_vec, ...
    Gident,tau_ident,RefPar.MatProps,...
    strain, time,SG, ...
    CondOpts,MainSaveDir,costDesig);

%% Plot Heat Maps with exact inputs
costDesig=strcat(ParentDesig,'_exactParam');
[phi_ctau_exact,phi_cG_exact,phi3D_exact] = func_PlotBulkCost(K_vec,G_vec,tau_vec, ...
    Gident,tau_exact,RefPar.MatProps,...
    strain, time,SG, ...
    CondOpts,MainSaveDir,costDesig);




