%This script is written to compare the performance of the  original
%constitutive model code calculating all in plane stresses
%(func_ViscoConstitutiveV6) and the constitutive model only dealing with 
%shear (func_ViscoShearConstitutive).

%It will be used to compare run time and the output of shear stress in the
%code


%Author: Andrew Matejunas

%Date completed: 2023/01/10

%Version History/Changelog:


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize (clear memory and figures)
close all
clear variables
clc

%% Remove some options from original IBII codes
hardCodePath=0;
FEValidMode=0;
calcKinFieldsFromDisp=1;
saveGridData=0;
%% Define Test Dseignation
TestDeg=char(cell2mat(inputdlg('Input Test Designation')));


 %% Create structute of known parameters
    Einf=2.98E9;
    nu=0.26;

    exactProps.Kinf=Einf/(3*(1-2*nu));
    exactProps.Ginf=Einf/(2*(1+nu));
    exactProps.nu=0;

%% Add reference parameters file 
[RefPar.file,RefPar.path]=uigetfile('*.mat',...
    'Choose Mat file containing reference constitutive pamrameters');
RefPar=load(strcat(RefPar.path,'/',RefPar.file),'MatProps');
MatProps=RefPar.MatProps;

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

%% Define Smoothing Opts

smoothingOpts.spatialSmooth=0;
smoothingOpts.FFTempSmooth=0;
smoothingOpts.WATempSmooth=0;

%% Define whether Noise should be added
quest='Add Noise?';
noiseChoice='No';

switch noiseChoice
    case 'Yes'
        %Addition by Andrew Matejunas
        imageNoise.addNoise=1;
        imageNoise.NumCopies=1;
        %% Choose noise magnitude
        noiseMag=str2double(inputdlg('Noise Magnitude',...
            'noise',[1,20],... %Dlgtitle, dims
            {'0.4000'})); %Default magnitude

    case 'No'
        imageNoise.addNoise=0;
end
    
   
%% Determine weather to run interpolated correction on the displacement 
    %fields
    
fprintf('Determining dispacement filed correction options \n')
DispCorr.Opt='Yes';
DispCorr.int=10;
DispCorr.Method='LinGrad';
DispCorr.PitchFitKern=2;
DispCorr.strainMethod='GridPeriod';
DispCorr.strainPitchNum=1;

%% Define Sample Condtioning Options

CondOpts.ImpCens=20;
CondOpts.FreeCens=20;
CondOpts.Xds=1;
CondOpts.Yds=1;
CondOpts.Tds=1;


if CondOpts.Tds <=1
    CondOpts.TempDS=false;
else
    CondOpts.TempDS=true;
end


%% Obtain strain fields
%             fprintf('------------------------------------------------------------- \n')
%             fprintf(strcat('Processing iteration',itNum,' \n'))
%             fprintf('------------------------------------------------------------- \n')
%% IMAGE PROCESSING: Use the Grid Method to extract displacement fields
if ~FEValidMode
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
end

%--------------------------------------------------------------------------
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
%fprintf('Obtaining and setting the free edge location.\n')
[freeEdge,specimen,disp] = func_getFreeEdge(globalOpts.hardCodeFreeEdge,...
    imagePath,imageFile,specimen,disp);

%--------------------------------------------------------------------------
% Smooth and Calculate Strain
fprintf('Calculating strain from the displacement fields.\n')
[disp,strain,strainRate] = func_smoothCalcStrain(globalOpts,pos,time,...
    grid,disp,smoothingOpts,extrapOpts);


fprintf('Correcting Y strain \n')
strain=func_PostCorrectStrain(grid,strain,pos,time,DispCorr,...
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

%% Calculate Stress Gage Stresses

fprintf('Calculating Stress Gage Stresses \n')
SG=func_Full_SG(accel,X_vec,time,material.rho);

%% Condition the data (censoring/downsampling/smoothing)
%% Define Censorship Parameter. (Standard is 2 grid pitches from
%impact and free edges)
fprintf('Conditioning IBII Fields \n')
[SG,accel,strain,X_vec]=func_conditionIBIIData(SG,accel,...
    strain,X_vec,time,CondOpts);

fprintf('Calculating the full in-plane stresses using ViscoConstitutiveV6 \n')
ConV6=func_ViscoConstitutiveV6(strain.x,strain.y,strain.s,...
    time.vec, MatProps,0,0,0);

fprintf('Calculating shear stress only using ViscoShearConstitutive \n')
ShearOnly=func_ViscoShearConstitutive(strain.s,...
    time.vec,MatProps);

%% Calculate error
fprintf('Calculating Error \n')
Error=(ShearOnly.xy-ConV6.xy)./ConV6.xy*100;

%%
fprintf('Done \n')

%% SaveResults
SaveName=strcat(savePath,'\',TestDeg,'_ConstV6vsShearComp.mat');
save(SaveName)

%% Plot Results

figure('units','normalized','outerposition',[0 0 1 1])

x_coord=pos.x*10^3;
y_coord=pos.y*10^3;
for k=1:length(time.vec)
    TConV6=squeeze(ConV6.xy(:,:,k))*10^-6;
    TShear=squeeze(ShearOnly.xy(:,:,k))*10^-6;
    TErr=squeeze(Error(:,:,k));
    timeinc=num2str(time.vec(k)*10^6);

    %%plot constitutive model v6 shear stress
    subplot(3,1,1)
    imagesc(x_coord,y_coord,TConV6)
    xlabel('X (mm)')
    ylabel('Y (mm)')
    title(strcat('ViscoConstitutiveV6 ',timeinc,'\mu s'))


    colormap('hot')
    cx=colorbar;
    cx.Label.String='\sigma_{xy} (MPa)';
    
    %%plot shear stress calculated with the new shear only algorithm
    subplot(3,1,2)
    imagesc(x_coord,y_coord,TShear)
    xlabel('X (mm)')
    ylabel('Y (mm)')
    title(strcat('Shear Only ',timeinc,'\mu s'))
    colormap('hot')
    cx=colorbar;
    cx.Label.String='\sigma_{xy} (MPa)';


    %Plot Errors
    subplot(3,1,3)
    imagesc(x_coord,y_coord,TErr)
    xlabel('X (mm)')
    ylabel('Y (mm)')
    title(strcat('Error',timeinc,'\mu s'))
    colormap('hot')
    cx=colorbar;
    cx.Label.String='Error (percent)';

  FigSaveName=strcat(savePath,'/',TestDeg,'ConstV6vsShearOnlyComp');
  saveas(gcf,FigSaveName,'fig')
   saveas(gcf,FigSaveName,'png')
    saveas(gcf,FigSaveName,'svg')
end