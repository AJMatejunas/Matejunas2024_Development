% This function is written to run a correction to correct for the errors in
    %the strain field along the top and bottom free edge caused by grid
    %method processing of FE data

%Author: Andrew Matejunas
%date completed:
%date verified:

    
%% initialize
clear all, close all, clc;

%% select grid method data file and define test designation
GMdatafile=uigetfile('.mat','choose GM data file');
desig=char(inputdlg('input test designation'));

%% define corection options
corOpts.corDisp=false;
corOpts.corStrain=true;
corOpts.corAccel=false;

prompt='How many grid pitches to correct?';
corOpts.CorNum=str2double(cell2mat(inputdlg(prompt)));
clear prompt


%% load the file
fprintf('Loading Data \n')
load(GMdatafile);
testdeg=desig;
clear desig

%% run correction
fprintf('Correct Strain Fields \n')
[disp,strain,accel]=func_StrainEdgeCorrection(grid,pos,...
    disp,strain,accel,time,...
    corOpts);
    
%% save data
fprintf('saving data \n')
save(strcat(testdeg,'_1PstrainCorr.mat'));







