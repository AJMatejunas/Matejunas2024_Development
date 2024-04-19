% This script is written to perform a  parametric sweep of different spatial
% and temporal smoothing techniques on processed grid method data. The
% smoothed data is then fed into the viscoelastic VFM parameter
% identification algorithm in order to determine optimal smoothing
% parameters.

%Author: Andrew Matejunas
%Date Completed:

%Version History:


%% Initialize
clear all
close all
clc

%% Define test designation and choose either processed grid data.
desig=char(inputdlg('input test designation'));

quest='data source?';
dlgtitle='source';
btn1='Images';
btn2='Processed GM Data';
defbtn='Images';
DataType=questdlg(quest,dlgtitle,btn1,btn2,defbtn);

clear quest dlgtitle btn1 btn2 defbtn
switch DataType
    case 'Processed GM Data'
[datafile,datapath]=uigetfile('*.mat','Choose GM SG data file');
end

%% Define an existing processing parameter file
[paramfile,parampath]=uigetfile('*.mat','Choose existing process parameter file');

%% Define material properties file
[propsfile,propspath]=uigetfile('*.mat','Choose file containing reference parameters');

%% Load Data
fprintf('Loading Data \n')

switch DataType
    case 'Processed GM Data'
load(strcat(datapath,'/',datafile));
end

load(strcat(parampath,'/',paramfile));
load(strcat(propspath,'/',propsfile));

fprintf('loading complete')
testdeg=desig;
clear desig

%% Set up image processing if using image data
switch DataType
    case 'Images'
ImgNoise=questdlg('Add Noise?');

switch ImgNoise
    case 'Yes'
        globalOpts.imageDefAddNoise=true;
% Image Noise Parameters
imageNoiseSweep.addNoise = true;
imageNoiseSweep.pcNoise = 0.4; % Normally 0.35-0.5% for HPVX
imageNoiseSweep.bits = 16;
imageNoiseSweep.convToUInt16 = true;    
    
    case 'No'
globalOpts.imageDefAddNoise=false;


       
end


if globalOpts.imageDefAddNoise==true

%Define number of noise copies
NumCopies=str2num(cell2mat(inputdlg('How many iterations of noise?',{},...
    [1,50],'30')));

else
imageNoiseSweep.addNoise=false;
NumCopies=1;

end
    
end

%% Set Up parallel pool
%Determine number of available workers and assign
simOpts.MaxWorkers=parcluster('local').NumWorkers;
simOPts.NumWorkers=simOpts.MaxWorkers;


%% define spatial smoothing parameters to try
spatialScaleFactor = 1; %used to compare different pixel arrays with
                            %similar smoothing parameters
sKern=inputdlg('Input vector of spatial smoothing kernels',...
    'spatial kernels',[1,50],'11,21,31,41,51,61');
sKernels=str2num(cell2mat(sKern));
clear sKern


%% Define temporal smoothing parameters to try
temporalScaleFactor = 1; %used to compare different frame rates with
                            %similar smoothing parameters
tKern=inputdlg('Input vector of temporal smoothing kernels',...
    'spatial kernels',[1,50],',11,15,21,25,31');
tKernels=str2num(cell2mat(tKern));
tOrder = 3; 

clear tKern


%% Initialize image Mask and Kernel Worker Assignment

switch DataType
    case 'Images'
fprintf('Loading reference image from the selected test data folder.\n')
hardCodePath = true;

[imageFile,imagePath] = uigetfile({'*.*','All Files'},...
    'Select the first image in the sequence');

% Specimen Location for Masking Grid Images
locPath = imagePath;
locFile = 'specimenLocation.mat';
    if exist([locPath,locFile],'file') == 2
        load([locPath,locFile]);
    else
        waitfor(warndlg({'Specimen location file not found.',...
            'Using default specimen location.'},'Warning!'))
        specimenLoc.bottomLeft = [14,5];
        specimenLoc.topRight = [398,246];
    end

end

