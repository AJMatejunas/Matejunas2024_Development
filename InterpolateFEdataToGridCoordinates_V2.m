% This script is written to interpolate and save the 
 %Author: Andrew Matejunas

 %Version History/Change Log
    %2023-04-25-V2: Added the option to choose between a nearest neightbor
        %interpolant and a linear interpolant

 %% Initialize
 clear variables; close all; clc

 %% Test designation
SimDesig=char(cell2mat(inputdlg('Test Designation')));

%% Load the reference parameter file
[refName,refPath]=uigetfile('*.mat', ...
    'Choose file containing reference material parameters');
refFile=strcat(refPath,'/',refName);

 %% load the finite element data
[FEname,FEpath]=uigetfile('*.mat', ...
    'Choose file containing finite element kinematic fields');
FEfile=strcat(FEpath,'/',FEname);

%% Load grid method pos file
[PosName,PosPath]=uigetfile('*.mat', ...
    'choose file containing grid method coordinates');
PosFile=strcat(PosPath,'/',PosName);

%% Choose directory to save results
SaveDir=uigetdir('','Choose folder in which to save results');

%% load selected filed
fprintf('Loading finite element data \n')
FE=load(FEfile,'pos','disp','accel','strain','stress','time');
time=FE.time;
fprintf('Loading material properties \n')
load(refFile);

fprintf('Loading the grid method coordinate system \n')
GridPos=load(PosFile,'pos');

%% Choose Interpolation Method
InterpMethod=questdlg('Choose Interpolation Method','Interp Method',...
    'nearest','linear','nearest');
%% Perform the interpolation
fprintf('Interpolating finite eleement deformations to grid coordinates \n')
[pos,disp,accel,strain,stress]=func_interpFEt2Grid(FE.pos,GridPos.pos,...
    FE.disp,FE.accel,FE.strain,FE.stress,InterpMethod);

%% Save the kinematic fields
fprintf('Saving Results \n')
save(strcat(SaveDir,'/',SimDesig,'_',InterpMethod,...
    '_interpolatedFields.mat'));

