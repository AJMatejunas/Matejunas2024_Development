%This script is written to plot the kinematic fields for a set of grid
%images both with and without significant smoothing


%% Choose folder to save data
savePath=uigetdir({},'Choose Folder to Save Results');

%% Load properties file
[refpar.name,refpar.path]=uigetfile('*.mat',...
    'Load File containing constitutive parameters');
refpar.File=strcat(refpar.path,'/',refpar.name);
load(refpar.File);

%% Remove some options from original IBII codes
hardCodePath=0;
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

%% Define Sample Condtioning Options
% prompt={'number of points to censor on impact edge',...
%     'censor on free edge',...
%     'Downsample factor in X','DS factor in Y','Temporal DS factor'};
% dims=[1, 50; 1, 50; 1, 50; 1, 50; 1,50];
% definputs={'10','10','1','1','1'};
% dlgtitle='define Specimen Conditioning options';
% CondParam=inputdlg(prompt,dlgtitle,dims,definputs);
% 
% CondOpts.ImpCens=str2num(CondParam{1});
% CondOpts.FreeCens=str2num(CondParam{2});
% CondOpts.Xds=str2num(CondParam{3});
% CondOpts.Yds=str2num(CondParam{4});
% CondOpts.Tds=str2num(CondParam{5});


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

clear prompt dims definputs dlgtitle CondParam

%% Create smoothing opts data structure without smoothing
smoothingOpts.spatialSmooth=0;
smoothingOpts.WATempSmooth=0;
smoothingOpts.FFTempSmooth=0;
smoothingOpts.spatialFilt='gauss';
smoothingOpts.spatialKernal=[3,3];
smoothingOpts.spatialEdgeMode='symmetric';
smoothingOpts.FFTemporalFilt='sgolay';
smoothingOpts.FFTemporalPad=0;
smoothingOpts.FFTemporalPadFrames=3;
smoothingOpts.FFTemporalPadMethod='replicate';
smoothingOpts.WATemporalAvgFirst=0;
smoothingOpts.WATemporalFilt='sgolay';
smoothingOpts.WATemporalKernal=[11,3];
smoothingOpts.WATemporalPad=0;
smoothingOpts.WATemporalPadFrames=3;
smoothingOpts.WATemporalPadMethod='replicate';

%% Generate full field data with No Noise or smoothing
fprintf('Calculating Fields without Noise or smoothing\n')
imageNoise.addNoise=0;

    
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
DispCorr.strainPitchNum=1;


%% %Preallocate identified constitutive parameter
% %     fprintf(strcat('Processing iteration',itNum,' \n'))
%     fprintf('------------------------------------------------------------- \n')
    %% IMAGE PROCESSING: Use the Grid Method to extract displacement fields
    if ~FEValidMode
        % Process the raw tiff images using the grid method code developed by
        % Grediac et al.
        fprintf('\n--------------------------------------------------------------\n')
        fprintf('GRID METHOD PROCESSING\n')
        fprintf('--------------------------------------------------------------\n')

        %--------------------------------------------------------------------------
        % GRID IMAGE PROCESSING



              
%             fprintf('Processing images using the grid method toolbox.\n')
            % Process the image squence with the grid method toolbox
            [grid,pos,disp] = func_gridMethodImageProcessing_AJM(imagePath,...
                imageFile,...
                grid,gridMethodOpts,imageNoise);    
             fprintf('Grid Method Processing Complete.\n')
       

        %--------------------------------------------------------------------------
        % Update Geometry and Number of Frames Based on Displacement Matrix Size
%         fprintf('Updating geometry variables based on the size of displacement field matrices.\n')
        [specimen,grid] = func_updateSpecGeom(specimen,grid,disp);

        % Currently the rotations are unused so remove them to save RAM
        disp = rmfield(disp,'rot');
    end

    %--------------------------------------------------------------------------
    % Create the time vector based on the number of frames in the disp struct
    time.numFrames = size(disp.x,3);
    time.vec = 0:time.step:(size(disp.x,3)-1)*time.step;

    %Create arrays of x and y vectors
    X_vec=pos.x;
    %Y_vec=pos.y;

    %% POST-PROCESSING: Smoothing and Kinematic Field Derivation
    % Smooth the displacement data and then calculate acceleration and strains
    % Extrapolate the data to account for the missing pitch on the edges
    fprintf('\n--------------------------------------------------------------\n')
    fprintf('POST PROCESSING: Smoothing and Kinematic Field Calculation\n')
    fprintf('--------------------------------------------------------------\n')
    %% Perfom corrections if needed
%     
%             fprintf('Saving Raw Displacement Fields \n')
%             RawDisp=disp;
            fprintf('Correcting Grid Method Displacements along specimen edges')
            [disp,DispCorr,grid,~]=func_CorrectGMDisp(disp,...
                DispCorr,grid,pos);
            
       


    % Calculate the kinematic fields from displacement fields using
    % displacements from images or displacements from FE data
    if calcKinFieldsFromDisp
        %--------------------------------------------------------------------------
        % Load the Reference Image and Determine Where the Free Edge is
%         fprintf('Obtaining and setting the free edge location.\n')
        [~,~,disp] = func_getFreeEdge(globalOpts.hardCodeFreeEdge,...
            imagePath,imageFile,specimen,disp);

        %--------------------------------------------------------------------------
        % Smooth and Calculate Strain
        fprintf('Calculating strain from the displacement fields.\n')
        [disp,strain,~] = func_smoothCalcStrain(globalOpts,pos,time,...
            grid,disp,smoothingOpts,extrapOpts);
        %Post Correct Strain
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
    end

    NSNN.disp=disp;
    NSNN.strain=strain;
    NSNN.accel=accel;
    NSNN.smoothingOpts=smoothingOpts;

%% Calculate with no smoothing with correction and noise
fprintf('Calculating Fields without smoothing with noise\n')
imageNoise.addNoise=1;
imageNoise.pcNoise=0.4;

    
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
DispCorr.strainPitchNum=1;


%% %Preallocate identified constitutive parameter
% %     fprintf(strcat('Processing iteration',itNum,' \n'))
%     fprintf('------------------------------------------------------------- \n')
    %% IMAGE PROCESSING: Use the Grid Method to extract displacement fields
    if ~FEValidMode
        % Process the raw tiff images using the grid method code developed by
        % Grediac et al.
        fprintf('\n--------------------------------------------------------------\n')
        fprintf('GRID METHOD PROCESSING\n')
        fprintf('--------------------------------------------------------------\n')

        %--------------------------------------------------------------------------
        % GRID IMAGE PROCESSING



              
%             fprintf('Processing images using the grid method toolbox.\n')
            % Process the image squence with the grid method toolbox
            [grid,pos,disp] = func_gridMethodImageProcessing_AJM(imagePath,...
                imageFile,...
                grid,gridMethodOpts,imageNoise);    
             fprintf('Grid Method Processing Complete.\n')
       

        %--------------------------------------------------------------------------
        % Update Geometry and Number of Frames Based on Displacement Matrix Size
%         fprintf('Updating geometry variables based on the size of displacement field matrices.\n')
        [specimen,grid] = func_updateSpecGeom(specimen,grid,disp);

        % Currently the rotations are unused so remove them to save RAM
        disp = rmfield(disp,'rot');
    end

    %--------------------------------------------------------------------------
    % Create the time vector based on the number of frames in the disp struct
    time.numFrames = size(disp.x,3);
    time.vec = 0:time.step:(size(disp.x,3)-1)*time.step;

    %Create arrays of x and y vectors
    X_vec=pos.x;
    %Y_vec=pos.y;

    %% POST-PROCESSING: Smoothing and Kinematic Field Derivation
    % Smooth the displacement data and then calculate acceleration and strains
    % Extrapolate the data to account for the missing pitch on the edges
    fprintf('\n--------------------------------------------------------------\n')
    fprintf('POST PROCESSING: Smoothing and Kinematic Field Calculation\n')
    fprintf('--------------------------------------------------------------\n')
    %% Perfom corrections if needed
%     
%             fprintf('Saving Raw Displacement Fields \n')
%             RawDisp=disp;
            fprintf('Correcting Grid Method Displacements along specimen edges')
            [disp,DispCorr,grid,~]=func_CorrectGMDisp(disp,...
                DispCorr,grid,pos);
            
       


    % Calculate the kinematic fields from displacement fields using
    % displacements from images or displacements from FE data
    if calcKinFieldsFromDisp
        %--------------------------------------------------------------------------
        % Load the Reference Image and Determine Where the Free Edge is
%         fprintf('Obtaining and setting the free edge location.\n')
        [~,~,disp] = func_getFreeEdge(globalOpts.hardCodeFreeEdge,...
            imagePath,imageFile,specimen,disp);

        %--------------------------------------------------------------------------
        % Smooth and Calculate Strain
        fprintf('Calculating strain from the displacement fields.\n')
        [disp,strain,~] = func_smoothCalcStrain(globalOpts,pos,time,...
            grid,disp,smoothingOpts,extrapOpts);
        %Post Correct Strain
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
    end

    NSN.disp=disp;
    NSN.strain=strain;
    NSN.accel=accel;
    NSN.smoothingOpts=smoothingOpts;

  %% Calculate with noise and smoothing 
  fprintf('Calculating fields with noise and smoothing 13px spatial 11 frames \n')
  %% Create smoothing opts data structure without smoothing
smoothingOpts.spatialSmooth=1;
smoothingOpts.WATempSmooth=0;
smoothingOpts.FFTempSmooth=1;
smoothingOpts.spatialFilt='gauss';
smoothingOpts.spatialKernal=[13,13];
smoothingOpts.spatialEdgeMode='symmetric';
smoothingOpts.FFTemporalFilt='sgolay';
smoothingOpts.FFTemporalPad=0;
smoothingOpts.FFTemporalPadFrames=3;
smoothingOpts.FFTemporalPadMethod='replicate';
smoothingOpts.WATemporalAvgFirst=0;
smoothingOpts.WATemporalFilt='sgolay';
smoothingOpts.WATemporalKernal=[11,3];
smoothingOpts.WATemporalPad=0;
smoothingOpts.WATemporalPadFrames=3;
smoothingOpts.WATemporalPadMethod='replicate';



fprintf('Calculating Fields without Noise or smoothing\n')
imageNoise.addNoise=1;
imageNoise.pcNoise=0.4;

    
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
DispCorr.strainPitchNum=1;


%% %Preallocate identified constitutive parameter
% %     fprintf(strcat('Processing iteration',itNum,' \n'))
%     fprintf('------------------------------------------------------------- \n')
    %% IMAGE PROCESSING: Use the Grid Method to extract displacement fields
    if ~FEValidMode
        % Process the raw tiff images using the grid method code developed by
        % Grediac et al.
        fprintf('\n--------------------------------------------------------------\n')
        fprintf('GRID METHOD PROCESSING\n')
        fprintf('--------------------------------------------------------------\n')

        %--------------------------------------------------------------------------
        % GRID IMAGE PROCESSING



              
%             fprintf('Processing images using the grid method toolbox.\n')
            % Process the image squence with the grid method toolbox
            [grid,pos,disp] = func_gridMethodImageProcessing_AJM(imagePath,...
                imageFile,...
                grid,gridMethodOpts,imageNoise);    
             fprintf('Grid Method Processing Complete.\n')
       

        %--------------------------------------------------------------------------
        % Update Geometry and Number of Frames Based on Displacement Matrix Size
%         fprintf('Updating geometry variables based on the size of displacement field matrices.\n')
        [specimen,grid] = func_updateSpecGeom(specimen,grid,disp);

        % Currently the rotations are unused so remove them to save RAM
        disp = rmfield(disp,'rot');
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
%     
%             fprintf('Saving Raw Displacement Fields \n')
%             RawDisp=disp;
            fprintf('Correcting Grid Method Displacements along specimen edges')
            [disp,DispCorr,grid,~]=func_CorrectGMDisp(disp,...
                DispCorr,grid,pos);
            
       


    % Calculate the kinematic fields from displacement fields using
    % displacements from images or displacements from FE data
    if calcKinFieldsFromDisp
        %--------------------------------------------------------------------------
        % Load the Reference Image and Determine Where the Free Edge is
%         fprintf('Obtaining and setting the free edge location.\n')
        [~,~,disp] = func_getFreeEdge(globalOpts.hardCodeFreeEdge,...
            imagePath,imageFile,specimen,disp);

        %--------------------------------------------------------------------------
        % Smooth and Calculate Strain
        fprintf('Calculating strain from the displacement fields.\n')
        [disp,strain,~] = func_smoothCalcStrain(globalOpts,pos,time,...
            grid,disp,smoothingOpts,extrapOpts);
        %Post Correct Strain
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
    end

    SN.disp=disp;
    SN.strain=strain;
    SN.accel=accel;
    SN.smoothingOpts=smoothingOpts;


 %% Plot kinematic Fields
for k=1:length(time.vec)
timePoint=num2str(time.vec(k)*10^6);
figure('units','normalized','outerposition',[0 0 1 1])

subplot(3,3,1)
colormap('hot');
tstrain=squeeze(NSNN.strain.x(:,:,k));
imagesc(X_vec,Y_vec,tstrain)
xlabel('X (mm)')
ylabel('Y (mm)')
cx=colorbar;
cx.Label.String='strain';
tstring=strcat('Strain_{xx} No Smoothing No Noise', timePoint, '\mnu s');
title(tstring);

subplot(3,3,2)
tstrain=squeeze(NSN.strain.x(:,:,k));
imagesc(X_vec,Y_vec,tstrain)
xlabel('X (mm)')
ylabel('Y (mm)')
cx=colorbar;
cx.Label.String='strain';
title('Strain_{xx} No Smoothing With Noise')

subplot(3,3,3)
tstrain=squeeze(SN.strain.x(:,:,k));
imagesc(X_vec,Y_vec,tstrain)
xlabel('X (mm)')
ylabel('Y (mm)')
cx=colorbar;
cx.Label.String='strain';
title('Strain_{xx} with Smoothing and Noise')


subplot(3,3,4)
tstrain=squeeze(NSNN.strain.y(:,:,k));
imagesc(X_vec,Y_vec,tstrain)
xlabel('X (mm)')
ylabel('Y (mm)')
cx=colorbar;
cx.Label.String='strain';
title('Strain_{yy} No Smoothing No Noise')

subplot(3,3,5)
tstrain=squeeze(NSN.strain.y(:,:,k));
imagesc(X_vec,Y_vec,tstrain)
xlabel('X (mm)')
ylabel('Y (mm)')
cx=colorbar;
cx.Label.String='strain';
title('Strain_{yy} No Smoothing With Noise')

subplot(3,3,6)
tstrain=squeeze(SN.strain.y(:,:,k));
imagesc(X_vec,Y_vec,tstrain)
xlabel('X (mm)')
ylabel('Y (mm)')
cx=colorbar;
cx.Label.String='strain';

title('Strain_{yy} with Smoothing and Noise')


subplot(3,3,7)
tstrain=squeeze(NSNN.strain.s(:,:,k));
imagesc(X_vec,Y_vec,tstrain)
xlabel('X (mm)')
ylabel('Y (mm)')
cx=colorbar;
cx.Label.String='strain';
title('Strain_{xy} No Smoothing No Noise')

subplot(3,3,8)
tstrain=squeeze(NSN.strain.s(:,:,k));
imagesc(X_vec,Y_vec,tstrain)
xlabel('X (mm)')
ylabel('Y (mm)')
cx=colorbar;
cx.Label.String='strain';
title('Strain_{xy} No Smoothing With Noise')

subplot(3,3,9)
tstrain=squeeze(SN.strain.s(:,:,k));
imagesc(X_vec,Y_vec,tstrain)
xlabel('X (mm)')
ylabel('Y (mm)')
cx=colorbar;
cx.Label.String='strain';

title('Strain_{xy} with Smoothing and Noise')

it=k;
FrameNum=num2str(it);
FrameSaveName=strcat(savePath,'NoiseAndSmoothStrainCompPlots_Frame', ...
    FrameNum,'.png');
saveas(gcf,FrameSaveName)
end

%% Save Console
SaveName=strcat(savePath,'FullFieldNoiseAndSmooth.mat');
save(SaveName);