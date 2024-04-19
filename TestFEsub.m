%This script is written to test the finite element substitution of the edge
    %fields
 

 %% Initialize
 close all; clear variables; clc

 %%
 ParentDesig='TESTING';
 
 %% Identify Directory to save analysis data

MainSaveDir=uigetdir('','Select Folder to save analysis');
%% Open file containing the processing parameters
[ProcessParam.Name,ProcessParam.Path]=uigetfile('*.mat',...
    'Choose Process Parameter Data File');
ProcessParam.File=strcat(ProcessParam.Path,'/',ProcessParam.Name);
fprintf('Loading Process Parameter Data File \n')
load(ProcessParam.File);

%% Choose Files Containing Kinematic Fields
        [Raw.Name,Raw.Path]=uigetfile('*.mat',...
            'Choose File Containing Kinematic Fields');
        Raw.File=strcat(Raw.Path,'\',Raw.Name);

        fprintf('Loading Kinematic Field Without Noise \n')
        RawFields=load(Raw.File);
        RawFields.File=Raw.File;

 %% Define Smoothing options
   smoothingOpts.spatialSmooth=false;
   smoothingOpts.spatialKernal=[0,0];
   smoothingOpts.FFTempSmooth=false;
   smoothingOpts.FFTemporalKernal=[0,3];
   RawFields.smoothingOpts=smoothingOpts;
   RawFields.globalOpts=globalOpts;
   RawFields.extrapOpts=extrapOpts;
   RawFields.diffOpts=diffOpts;


   %% Load Image Def File
   [ImDef.Name,ImDef.Path]=uigetfile('*.mat',...
       'Choose Grid Image Deformation File');
   ImDef.File=strcat(ImDef.Path,'/',ImDef.Name);
   %% Interpolate Finite element data onto grid coordinates
   subset=1;
   fprintf('Interpolating Finite Element Data to Grid Coordinates \n')
   FEinterp=func_InterpFEtoGM(ImDef.File,MainSaveDir,ParentDesig,subset);

%% Run the substitution
SubDesig='EDGEsubTest';
subOpts.method='GridPeriod';
subOpts.PitchNum=2;
   FESub=func_subGMEdgeFields(RawFields.time,FEinterp.disp, ...
            RawFields.disp,RawFields.strain,RawFields.accel, ...
            RawFields.pos,RawFields.grid,... 
            MainSaveDir,SubDesig, ...
            subOpts,RawFields.smoothingOpts,RawFields.globalOpts, ...
            RawFields.extrapOpts,RawFields.diffOpts);
       