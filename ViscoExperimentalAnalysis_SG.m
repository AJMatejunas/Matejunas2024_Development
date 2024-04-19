%% IMAGE BASED INERTIAL IMPACT (IBII) TEST : Data Processing - v1.0r
% Authors: Lloyd Fletcher, Jared Van-Blitterswyk, Andrew Matejunas
% PhotoDyn Group, University of Southampton
% http://photodyn.org/
% Date Created: 12th Dec. 2017 - v0.1
% Date Edited: 6th Jun. 2019 - v1.0r release
% Date Edited: 24th Feb. 2022- ExtracFieldsGrid-V1
               %Added capabilities to correct displacement fields to
                %account for errors from grid method extraction within one
                %pitch of the free surfaces and by shear induced rotation of
                %the free surfaces (these errors can be more than one grid
                %pitch)
 %
% Takes input IBII test images processes them with the grid method code 
% to obtain displacement fields. Further kinematic fields are derived from 
% the displacement fields (acceleration, strain) after applying smoothing 
% to the displacement fields. Kinematic fields are then used to identify 
% material properties including: stiffness and strength. This version
% includes two material models: reduced orthotropic and isotropic.
%
% This code uses the grid processing tool box that can be found at:
% www.thegridmethod.net, developed by Grediac et al.
%
% The following papers describe the IBII method:
%[1] L. Fletcher, J. Van-Blitterswyk, F. Pierron, A Novel Image-Based 
%   Inertial Impact Test (IBII) for the Transverse Properties of 
%   Composites at High Strain Rates, J. Dynamic Behavior Mater. (2019). 
%   doi:10.1007/s40870-019-00186-y.
%[2] L. Fletcher, F. Pierron, An image-based inertial impact (IBII) test 
%   for tungsten carbide cermets, J. Dynamic Behavior Mater. 4 (2018) 
%   481–504. doi:10.1007/s40870-018-0172-4.
%[3] J. Van Blitterswyk, L. Fletcher, F. Pierron, Image-based inertial 
%   impact test for composite interlaminar tensile properties, J. 
%   Dynamic Behavior Mater. 4 (2018) 543–572. doi:10.1007/s40870-018-0175-1.
%
% This work is licensed under the Creative Commons Attribution-
% NonCommercial-ShareAlike 4.0 International License. To view a copy of 
% this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or 
% send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, 
% USA.
%
% If there are any issues with this code please contact: 
% Lloyd Fletcher: l.c.fletcher@soton.ac.uk / lloydcolinfletcher@gmail.com 

clc
clear all
close all

fprintf('--------------------------------------------------------------\n')
fprintf('IMAGE-BASED INERTIAL IMPACT (IBII): Data Processing v1.0r\n')
fprintf('--------------------------------------------------------------\n')

%% Set test designation
testdeg=char(inputdlg('Test designation'));

%%
% Flag for processing pure FE data for validation purposes, bypassess the
% grid method processing module and directly loads FE kinematic data
FEValidMode =false;
calcKinFieldsFromDisp = ~FEValidMode;

%% INITIALISE: Add path for processing functions
fprintf('Adding pathes for processing functions.\n')
% Add the path for the grid method code and other useful functions:
funcPath = [pwd,'\Functions\'];

% If the default path is not found we should find it
if exist(funcPath,'file') ~= 7
    hWarn = warndlg('Folder for processing functions not found.','Function folder not found');
    waitfor(hWarn);
    funcPath = uigetdir(pwd,'Locate Processing Function Folder');
end
addpath(funcPath);
addpath([funcPath,'GridMethodToolbox\']);

%% LOAD RAW DATA FILES: Images or FE
hardCodePath = false; % Useful for running a single file repeatedly for debugging
if ~FEValidMode    
    % LOAD RAW DATA FILES: Raw .tiff images
    fprintf('Loading reference image from the selected test data folder.\n')
    if ~hardCodePath
        [imageFile,imagePath] = uigetfile({'*.*','All Files'},'Select the first image in the sequence');
    else   
        imageFile = 'DefGridImage_001.tiff';
        imagePath = 'E:\Matlab_WorkingDirectory\1_IBIITest_Data\TestData_ID_CFIPQIso_HalfPulse\'; 
    end
else
    % LOAD RAW DATA FILES: FE Validation Data
    fprintf('Loading FE data from selected .mat file.\n')
    if ~hardCodePath
        [FEDataFile,FEDataPath] = uigetfile({'*.*','All Files'},'Select the FE data file');
    else   
        FEDataFile = '';
        FEDataPath = [pwd,'\'];    
    end
    load([FEDataPath,FEDataFile])
    imagePath = FEDataPath;
    imageFile = FEDataFile;   
end

%% INITIALISE: Load Processing Parameter Data Structures
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

%% Determine weather to run 
quest='Correct Disp fields?';
dlgtitle='Displacement correction';
DispCorr.Opt=questdlg(quest,dlgtitle,'Yes');
clear quest dlgtitle

DispCorr.Opt='Yes';
DispCorr.int=10;
DispCorr.Method='LinGrad';
DispCorr.PitchFitKern=2;
DispCorr.strainMethod='GridPeriod';
DispCorr.strainPitchNum=1;

        
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

    fprintf('Checking for existing processed data file.\n')
    processGridImages = true;
    if exist([gridDataSavePath,gridDataFile],'file') == 2
        if gridMethodOpts.autoLoadProccessedDataFile
            processGridImages = false;
        else   
            choice = questdlg('Processed grid data file found, process again?', ...
            'Process Grid Images', 'Yes','No','No');
            switch choice
                case 'Yes'
                    processGridImages = true;
                case 'No'
                    processGridImages = false;
            end
        end
    end

    if processGridImages
        fprintf('Processing images using the grid method toolbox.\n')
        % Process the image squence with the grid method toolbox
        [grid,pos,disp] = func_gridMethodImageProcessing(imagePath,imageFile,...
            grid,gridMethodOpts,imageNoise);    

        % Save the data to file so a reprocess is not necessary
        fprintf('Saving processed grid data to: %s.\n',gridDataFile)
        save([gridDataSavePath,gridDataFile],'grid','pos','disp') 

        fprintf('Grid Method Processing Complete.\n')
    else
        fprintf('Loading processed grid data from: %s.\n',gridDataFile)
            load([gridDataSavePath,gridDataFile]);
    end

    %--------------------------------------------------------------------------
    % Update Geometry and Number of Frames Based on Displacement Matrix Size
    fprintf('Updating geometry variables based on the size of displacement field matrices.\n')
    [specimen,grid] = func_updateSpecGeom(specimen,grid,disp);

    % Currently the rotations are unused so remove them to save RAM
    disp = rmfield(disp,'rot');
end

%--------------------------------------------------------------------------
% Create the time vector based on the number of frames in the disp struct
time.numFrames = size(disp.x,3);
time.vec = 0:time.step:(size(disp.x,3)-1)*time.step;


%% POST-PROCESSING: Smoothing and Kinematic Field Derivation
% Smooth the displacement data and then calculate acceleration and strains
% Extrapolate the data to account for the missing pitch on the edges
fprintf('\n--------------------------------------------------------------\n')
fprintf('POST PROCESSING: Smoothing and Kinematic Field Calculation\n')
fprintf('--------------------------------------------------------------\n')
%% Perfom corrections if needed
fprintf('Saving Raw Displacement Fields \n')
RawDisp=disp;
fprintf('Correcting Grid Method Displacements along specimen edges')
[disp,DispCorr,grid,ProgramVersions]=func_CorrectGMDisp(disp,...
    DispCorr,grid,pos);

% Calculate the kinematic fields from displacement fields using
% displacements from images or displacements from FE data
if calcKinFieldsFromDisp
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
else
    fprintf('Using kinematic fields directly from FE data.\n')
    % For pure FE data the kinematic fields can be taken directly from the
    % FE calculation
    disp.xAvg = func_avgFFVarOverWidth(disp.x);
    accel.xAvg = func_avgFFVarOverWidth(accel.x);
    strain.xAvg = func_avgFFVarOverWidth(strain.x);
    strain.yAvg = func_avgFFVarOverWidth(strain.y);
    [strainRate,~] = func_calculateStrainRate(strain,time,smoothingOpts,true);
    strainRate.xAvg = func_avgFFVarOverWidth(strainRate.x);
    strainRate = func_calcMaxStrainRate(smoothingOpts,grid,strainRate);
    disp.extrap.x = disp.x;
    disp.extrap.y = disp.y;
end

%% POST-PROCESSING - Image Sequences of Kinematic Fields
fprintf('\n--------------------------------------------------------------\n')
fprintf('POST PROCESSING: Image Sequences\n')
fprintf('--------------------------------------------------------------\n')
plotProps = func_initPlotPropsStruct(plotParams.formatType);

%--------------------------------------------------------------------------
% Create Video Sequences of Disp, Strain, Accel and Strain Rate
% NOTE: these images are required for the crack location identification
if strcmp(globalOpts.plotImageSeqs,'no')
    fprintf('Image sequence plotting disabled.\n')
else
    fprintf('Plotting image sequence of disp, accel, strain and strain rate.\n')
    fprintf('NOTE: image sequence used to identify crack location for failure stress identification.\n')
    func_plotAllFullFieldVideos(globalOpts,plotParams,imagePath,pos,time,...
        disp,strain,strainRate,accel);
end
close all

%% Save data
fprintf('Saving Data \n')
[GMsaveName,GMsavePath]=uiputfile('*.mat','Save GM data',strcat(testdeg,'_GMdata.mat'));
save(strcat(GMsavePath,'/',GMsaveName))

fprintf('GM data extraction complete \n')

%% Calculate the stress guage stresses
X_vec=pos.x;
SGstress=func_Full_SG(accel,X_vec,time,material.rho);

%% Plot SG stresses heat map
figure('units','inches','innerposition',[5,5,7,3.5])
imagesc(X_vec*10^3,time.vec*10^6,SGstress.x'*10^-6)
xlabel('X position (mm)')
ylabel('time (\mu s)')
cx=colorbar;
cx.Label.String='Stress (MPa)';
set(gca,'YDir','normal')
set(gca,'FontSize',14)

figSaveName=strcat(testdeg,'_SGxHeat');
saveas(gcf,figSaveName,'fig');
saveas(gcf,figSaveName,'svg');
saveas(gcf,figSaveName,'png')

%% plot shear stresses in a heat map
figure('units','inches','innerposition',[5,5,7,3.5])
imagesc(X_vec*10^3,time.vec*10^6,SGstress.s'*10^-6)
xlabel('X position (mm)')
ylabel('time (\mu s)')
cx=colorbar;
cx.Label.String=['Shear Stress (MPa)'];
set(gca,'YDir','normal')
set(gca,'FontSize',14)

figSaveName=strcat(testdeg,'_SGsHeat');
saveas(gcf,figSaveName,'fig');
saveas(gcf,figSaveName,'svg');
saveas(gcf,figSaveName,'png')

%% calculate average shear strains
strain.sAvg=squeeze(mean(strain.s));

%% Choose indexes to pull data from
SR=grid.pxPerPeriod;
freeInd=4*SR;
impInd=length(pos.x)-4*SR;
midInd=round(length(pos.x)/2);



%% Plot 1-D strain-time plots
timemic=time.vec*10^6;
figure('units','inches','innerposition',[5,5,7,3.5])
subplot(1,2,1)
plot(timemic,strain.xAvg(impInd,:),'k')
hold on
plot(timemic,strain.xAvg(midInd,:),'--R')
plot(timemic,strain.xAvg(freeInd,:),'.-B')
hold off
xlabel('time (\mu s)')
ylabel('strain_{xx}')
legend('4P from Imp','Mid','4P from Free','location','southwest')

subplot(1,2,2)
plot(timemic,strain.sAvg(impInd,:),'k')
hold on
plot(timemic,strain.sAvg(midInd,:),'--R')
plot(timemic,strain.sAvg(freeInd,:),'-.B')
hold off
xlabel('time (\mu s)')
ylabel('strain_{xy}')
legend('4P from Imp','Mid','4P from Free','location','southwest')

figSaveName=strcat(testdeg,'_1DStrainTime');
saveas(gcf,figSaveName,'fig')
saveas(gcf,figSaveName,'png')
saveas(gcf,figSaveName,'svg')

%% Plot 1D SG-time plots
timemic=time.vec*10^6;
figure('units','inches','innerposition',[5,5,7,3.5])
subplot(1,2,1)
plot(timemic,SGstress.x(impInd,:)*10^-6,'k')
hold on
plot(timemic,SGstress.x(midInd,:)*10^-6,'--R')
plot(timemic,SGstress.x(freeInd,:)*10^-6,'.-B')
hold off
xlabel('time (\mu s)')
ylabel('SG_{xx} (MPa)')
legend('4P from Imp','Mid','4P from Free','location','southwest')

subplot(1,2,2)
plot(timemic,SGstress.s(impInd,:)*10^-6,'k')
hold on
plot(timemic,SGstress.s(midInd,:)*10^-6,'--R')
plot(timemic,SGstress.s(freeInd,:)*10^-6,'-.B')
hold off
xlabel('time (\mu s)')
ylabel('SG_{xy} (MPa)')
legend('4P from Imp','Mid','4P from Free','location','southwest')

figSaveName=strcat(testdeg,'_1DsgTime');
saveas(gcf,figSaveName,'fig')
saveas(gcf,figSaveName,'png')
saveas(gcf,figSaveName,'svg')

%% Plot stress strain
figure('units','inches','innerposition',[5,5,7,3.5])
subplot(1,2,1)
plot(strain.xAvg(impInd,:),SGstress.x(impInd,:)*10^-6,'k')
hold on
plot(strain.xAvg(midInd,:),SGstress.x(midInd,:)*10^-6,'--R')
plot(strain.xAvg(freeInd,:),SGstress.x(freeInd,:)*10^-6,'.-B')
hold off
xlabel('strain_{xx}')
ylabel('SG_{xx} (MPa)')
legend('4P from Imp','Mid','4P from Free','location','southwest')

subplot(1,2,2)
plot(strain.sAvg(impInd,:),SGstress.s(impInd,:)*10^-6,'k')
hold on
plot(strain.sAvg(midInd,:),SGstress.s(midInd,:)*10^-6,'--R')
plot(strain.sAvg(freeInd,:),SGstress.s(freeInd,:)*10^-6,'-.B')
hold off
xlabel('strain_{xy}')
ylabel('SG_{xy} (MPa)')
legend('4P from Imp','Mid','4P from Free','location','southwest')

figSaveName=strcat(testdeg,'_1DSGStrain');
saveas(gcf,figSaveName,'fig')
saveas(gcf,figSaveName,'png')
saveas(gcf,figSaveName,'svg')