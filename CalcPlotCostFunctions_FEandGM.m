%This script is written to calculate and plot the weighted cost fucntions
    %with both finite element and grid method deformation data.
        %The code will plot the cost function using exact values along with
            %5% errors in K, G, and tau

%Author: Andrew Matejunas

%Date created: 2023/03/05

%Version history/changelog



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize
clear variables; close all; clc

%% Select Directory to save files in
MainDir=uigetdir('','Choose folder to save results');
%% Add reference parameters file 
[RefPar.file,RefPar.path]=uigetfile('*.mat',...
    'Choose Mat file containing reference constitutive pamrameters');

%% Select finite element data (note only loading certain values to save
    %memory)
[FEname,FEpath]=uigetfile('*.mat', ...
    'Select file containing Finite Element data');
FEfile=strcat(FEpath,'/',FEname);

%% Load Processing Parameter Data Structures
fprintf('Loading processing parameters file.\n')
 [initFile,initPath,~] = uigetfile('*.mat', ...
        'Locate processing parameter file');

% Load the processing parameters from file


load([initPath,initFile])

% Store the FE valid mode parameter
globalOpts.FEValidMode =false;
globalOpts.calcKinFieldsFromDisp = true;

%% Load selected files
load(strcat(RefPar.path,'/',RefPar.file),'MatProps');
fprintf('Loading the finite element fields \n')
FE=load(FEfile,'disp','strain','accel','pos');

%% Main designation
MainDesig=char(cell2mat(inputdlg('Input Test Designation')));
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
          [imageFile,imagePath] = uigetfile({'*.*','All Files'},'Select the first image in the sequence');
   
%% Generate known parameter data structure
knownParam=[MatProps.Kinf,MatProps.Ginf,MatProps.nu,0];

%% Process Grid Data
           
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
            [disp.x,disp.y,disp.rX,disp.rY]=...
                func_cropAndExtrapDispFields_v4(pos,disp.x,disp.y, ...
                extrapOpts.disp,true);
            

            % Calculate the kinematic fields from displacement fields using
            % displacements from images or displacements from FE data
           
                %--------------------------------------------------------------------------
                % Load the Reference Image and Determine Where the Free Edge is
                %fprintf('Obtaining and setting the free edge location.\n')
                [freeEdge,specimen,disp] =...
                    func_getFreeEdge(globalOpts.hardCodeFreeEdge,...
                    imagePath,imageFile,specimen,disp);

                %--------------------------------------------------------------------------
                % Smooth and Calculate Strain
                fprintf('Calculating strain from the displacement fields.\n')

                [strain,~]=func_smoothCalcStrain_v4(pos,time,disp,...
                    smoothOpts.strain,extrapOpts.strain,true);



                %--------------------------------------------------------------------------
                % Smooth and Calculate Acceleration
                fprintf('Calculating acceleration from the displacement fields.\n')
                [accel,~,~] = func_smoothCalcAccel_v4(pos,time,disp, ...
                    smoothOpts.accel,...
                    extrapOpts.accel,diffOpts,true);

                % Remove some 3D fields from the structs to save memory.
               
            %% Calculate Stress Gage Stresses

            fprintf('Calculating Stress Gage Stresses \n')
            SG=func_Full_SG(accel,X_vec,time,material.rho);
            FE.SG=func_Full_SG(FE.accel,FE.pos.x,time,material.rho);

            %% Condition the data (censoring/downsampling/smoothing)
            fprintf('Conditioning IBII Fields \n')
            timefull=time;
            [SG,accel,strain,X_vec,time]=...
            func_conditionIBIIDataV2(SG,accel,...
                strain,X_vec,time,CondOpts);

%% Condition the Finite element fields to match that of the grid method data
FE.CondOpts=CondOpts;
FE.CondOpts.ImpCens=length(...
    FE.pos.x(FE.pos.x>=pos.x(end-CondOpts.ImpCens)));
FE.CondOpts.FreeCens=length(...
    FE.pos.x(FE.pos.x<=pos.x(CondOpts.FreeCens)));

[FE.SG,FE.accel,FE.strain,~,FE.time]=func_conditionIBIIDataV2(...
    FE.SG,FE.accel,FE.strain,FE.pos.x,timefull,FE.CondOpts);


%% Ident Params structure for exact property inputs
identParams.K=MatProps.Ki;
identParams.G=MatProps.Gi;
identParams.tau=MatProps.tau;

%% reference Parameter data structure
refParams.K=MatProps.Ki;
refParams.G=MatProps.Gi;
refParams.tau=MatProps.tau;

%% Define bounds of parameter ratios to plot
ub=1.25*[1,1,1];
lb=0.75*[1,1,1];
%% Calculate and plot the FE cost fucntion for exact inputs
fprintf('Plotting FE cost functions with exact inputs \n')
CFdesig=strcat(MainDesig,'_exact');
[FE.phiKTax.exact,FE.phiKTtot.exact,FE.phiKGax.exact,...
    FE.phiKGtot.exact,FE.phiGTax.exact,...
    FE.phiGTtot.exact,FE.phiGTs.exact] =...
    func_plotCostFunction(FE.SG,knownParam,FE.strain,time.vec,FE.CondOpts, ...
    SolveOpts.wK,SolveOpts.wG,...
    refParams,identParams,ub,lb,MainDir,CFdesig,'FE');

%% Calculate and plot the cost function for the grid method with exact 
% inputs
fprintf('Plotting Grid Method cost function with exact inputs \n')
[GM.phiKTax.exact,GM.phiKTtot.exact,GM.phiKGax.exact,...
    GM.phiKGtot.exact,GM.phiGTax.exact,...
    GM.phiGTtot.exact,GM.phiGTs.exact] =...
    func_plotCostFunction(SG,knownParam,strain,time.vec,CondOpts, ...
    SolveOpts.wK,SolveOpts.wG,...
    refParams,identParams,ub,lb,MainDir,CFdesig,'GM');

%% Ident Params structure for 3% error in G property inputs
identParams.K=MatProps.Ki;
identParams.G=MatProps.Gi*(1+.03);
identParams.tau=MatProps.tau;

%% Calculate and plot the FE cost fucntion for exact inputs
fprintf('Plotting FE cost functions with 3percent error in G \n')
CFdesig=strcat(MainDesig,'_3percentG');
[FE.phiKTax.G3p,FE.phiKTtot.G3p,FE.phiKGax.G3p,...
    FE.phiKGtot.G3p,FE.phiGTax.G3p,...
    FE.phiGTtot.G3p,FE.phiGTs.G3p] =...
    func_plotCostFunction(FE.SG,knownParam,FE.strain,time.vec,FE.CondOpts, ...
    SolveOpts.wK,SolveOpts.wG,...
    refParams,identParams,ub,lb,MainDir,CFdesig,'FE');

%% Calculate and plot the cost function for the grid method with exact 
% inputs
fprintf('Plotting Grid Method cost function with 3percent error in G \n')
[GM.phiKTax.G3p,GM.phiKTtot.G3p,GM.phiKGax.G3p,...
    GM.phiKGtot.G3p,GM.phiGTax.G3p,...
    GM.phiGTtot.G3p,GM.phiGTs.G3p] =...
    func_plotCostFunction(SG,knownParam,strain,time.vec,CondOpts, ...
    SolveOpts.wK,SolveOpts.wG,...
    refParams,identParams,ub,lb,MainDir,CFdesig,'GM');

%% Ident Params structure for 3% error in G property inputs
identParams.K=MatProps.Ki;
identParams.G=MatProps.Gi*(1+.03);
identParams.tau=MatProps.tau*(1+.03);

%% Calculate and plot the FE cost fucntion for exact inputs
fprintf('Plotting FE cost functions with 3percent error in G and tau \n')
CFdesig=strcat(MainDesig,'_3percentGandTau');
[FE.phiKTax.GT3p,FE.phiKTtot.GT3p,FE.phiKGax.GT3p,...
    FE.phiKGtot.GT3p,FE.phiGTax.GT3p,...
    FE.phiGTtot.GT3p,FE.phiGTs.GT3p] =...
    func_plotCostFunction(FE.SG,knownParam,FE.strain,time.vec,FE.CondOpts, ...
    SolveOpts.wK,SolveOpts.wG,...
    refParams,identParams,ub,lb,MainDir,CFdesig,'FE');

%% Calculate and plot the cost function for the grid method with exact 
% inputs
fprintf('Plotting Grid Method cost function with 3percent error in G and tau \n')
[GM.phiKTax.GT3p,GM.phiKTtot.GT3p,GM.phiKGax.GT3p,...
    GM.phiKGtot.GT3p,GM.phiGTax.GT3p,...
    GM.phiGTtot.GT3p,GM.phiGTs.GT3p] =...
    func_plotCostFunction(SG,knownParam,strain,time.vec,CondOpts, ...
    SolveOpts.wK,SolveOpts.wG,...
    refParams,identParams,ub,lb,MainDir,CFdesig,'GM');

%% Ident Params structure for 3% error in G property inputs
identParams.K=MatProps.Ki;
identParams.G=MatProps.Gi;
identParams.tau=MatProps.tau*(1+.03);

%% Calculate and plot the FE cost fucntion for exact inputs
fprintf('Plotting FE cost functions with 3percent error in tau \n')
CFdesig=strcat(MainDesig,'_3percentTau');
[FE.phiKTax.T3p,FE.phiKTtot.T3p,FE.phiKGax.T3p,...
    FE.phiKGtot.T3p,FE.phiGTax.T3p,...
    FE.phiGTtot.T3p,FE.phiGTs.T3p] =...
    func_plotCostFunction(FE.SG,knownParam,FE.strain,time.vec,FE.CondOpts, ...
    SolveOpts.wK,SolveOpts.wG,...
    refParams,identParams,ub,lb,MainDir,CFdesig,'FE');

%% Calculate and plot the cost function for the grid method with exact 
% inputs
fprintf('Plotting Grid Method cost function with 3percent error in tau \n')
[GM.phiKTax.T3p,GM.phiKTtot.T3p,GM.phiKGax.T3p,...
    GM.phiKGtot.T3p,GM.phiGTax.T3p,...
    GM.phiGTtot.T3p,GM.phiGTs.T3p] =...
    func_plotCostFunction(SG,knownParam,strain,time.vec,CondOpts, ...
    SolveOpts.wK,SolveOpts.wG,...
    refParams,identParams,ub,lb,MainDir,CFdesig,'GM');


