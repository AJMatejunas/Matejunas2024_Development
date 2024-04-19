% This code is written to compare identification results betweeen raw grid
    %Method data, extrapolated grid method data, and grid method data where
    %the edge displacements are replaced with finite element data


%% Initialize
clear variables; close all; clc

%% Choose Save Directory
SaveDir=uigetfile('','Choose Directory to save results file');

%% Load interpolated data file
[INTname,INTpath]=uigetfile('*.mat', ...
    'Choose file containing interpolated kinematic fields');
INTfile=strcat(INTpath,'/',INTname);
%% Add reference parameters file 
[RefPar.file,RefPar.path]=uigetfile('*.mat',...
    'Choose Mat file containing reference constitutive pamrameters');
fprintf('Loading Reference Parameter File \n')
load(strcat(RefPar.path,'/',RefPar.file),'MatProps');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%3

%% Remove some options from original IBII codes
hardCodePath=0;
FEValidMode=0;
calcKinFieldsFromDisp=1;
saveGridData=0;

%% Create reasonable upper bounds
%bounds on the initial guess are +/- 25%
%% Create reasonable upper bounds
%bounds on the initial guess are +/- 25%
ub=[1.875e9;
    1.0962E9;... G1
    12.5E-6]; %tau_1
ubG=ub(2:3);

%% Create reasonable lower bounds
lb=[1.125e9;
    657.7E6;... G1
    7.5E-6]; %tau_1
lbG=lb(2:3);

%% Create Matrix of initial guesses
    %Eint=2.5E9;
    %Choose a random 3x3 matrix of initial guress multipliers. Multipliers
        %will be uniformly distributed random numbers in the interval [0,1]
 
    intMult=rand(3,3);
    intMultK=rand(1,3);
  
    %Create 3X3 matrix of initial guess parameters
     intMatrix=lb+(ub-lb).*intMult;
     intMatrix=intMatrix(2:3,:);

%     intGuess=[Eint/(3*(1-2*nu));... K1
%               Eint/(2*(1+nu));.... G1
%               10E-6];... tau1
   


%% Create structute of known parameters
    Einf=2.98E9;
    nu=0.26;

    exactProps.Kinf=Einf/(3*(1-2*nu));
    exactProps.Ginf=Einf/(2*(1+nu));
    exactProps.nu=0;
    exactProps.Ki=RefPar.MatProps.Ki;
   
%% Calculate the the kinematic fields from images
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n')
fprintf('Calculating Kinematic Fields From Images without edge corrections \n')
fprintf('--------------------------------------------------------------- \n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
fprintf('Loading processing parameters file without edge extrapolation.\n')
initPath = imagePath;
initFile = 'processingParameters.mat';
if exist([initPath,initFile],'file') ~= 2
    hWarn = warndlg('Processing parameter file does not exist.','Processing parameters not found');
    waitfor(hWarn);
    [initFile,initPath,~] = uigetfile('*.mat','Locate processing parameter file');
end
% Load the processing parameters from file


load([initPath,initFile])

% Store the FE valid mode parametervccvv.MatProps.Gi*100;
globalOpts.FEValidMode = FEValidMode;
globalOpts.calcKinFieldsFromDisp = calcKinFieldsFromDisp;

%% Load Interpolated data file
INT=load(INTfile,'disp','pos');
%% Define Test Dseignation
TestDeg=char(cell2mat(inputdlg('Input Test Designation')));

%% Generate Kinematic Fields for Images Without Corrections
        timefull=time;
        %             fprintf('------------------------------------------------------------- \n')
        %             fprintf(strcat('Processing iteration',itNum,' \n'))
        %             fprintf('------------------------------------------------------------- \n')
        % IMAGE PROCESSING: Use the Grid Method to extract displacement fields

        % Process the raw tiff images using the grid method code developed by
        % Grediac et al.
        %                 fprintf('\n--------------------------------------------------------------\n')
        %                 fprintf('GRID METHOD PROCESSING\n')
        %                 fprintf('--------------------------------------------------------------\n')
        %
        %                 %--------------------------------------------------------------------------
        %                 % GRID IMAGE PROCESSING

        % fprintf('Processing images using the grid method toolbox.\n')
        % Process the image squence with the grid method toolbox
        [grid,pos,disp] = func_gridMethodImageProcessing_AJM(imagePath,...
            imageFile,...
            grid,gridMethodOpts,imageNoise);

        %--------------------------------------------------------------------------
        % Update Geometry and Number of Frames Based on Displacement Matrix Size
        %fprintf('Updating geometry variables based on the size of displacement field matrices.\n')
        [specimen,~] = func_updateSpecGeom(specimen,grid,disp);

        % Currently the rotations are unused so remove them to save RAM
        %disp = rmfield(disp,'rot');


        %--------------------------------------------------------------------------
        % Create the time vector based on the number of frames in the disp struct
        % time.numFrames = size(disp.x,3);
        %time.vec = 0:time.step:(size(disp.x,3)-1)*time.step;

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

        % POST-PROCESSING: Smoothing and Kinematic Field Derivation
        % Smooth the displacement data and then calculate acceleration and strains
        % Extrapolate the data to account for the missing pitch on the edges
        %             fprintf('\n--------------------------------------------------------------\n')
        %             fprintf('POST PROCESSING: Smoothing and Kinematic Field Calculation\n')
        %             fprintf('--------------------------------------------------------------\n')
        % Crop and Extrapolate displacement fields
        [disp.x,disp.y,disp.rX,disp.rY]=...
            func_cropAndExtrapDispFields_v4(pos,disp.x,disp.y, ...
            extrapOpts.disp,false);


        % Calculate the kinematic fields from displacement fields using
        % displacements from images or displacements from FE data

        %--------------------------------------------------------------------------
        % Load the Reference Image and Determine Where the Free Edge is
        %fprintf('Obtaining and setting the free edge location.\n')
        [freeEdge,specimenN,disp] =...
            func_getFreeEdge(globalEdge,...
            imagePath,imageFile,specimen,disp);

        %--------------------------------------------------------------------------
        % Smooth and Calculate Strain
        % fprintf('Calculating strain from the displacement fields.\n')

        [strain,~]=func_smoothCalcStrain_v4(pos,time,disp,...
            smoothOpts.strain,extrapOpts.strain,false);

        %--------------------------------------------------------------------------
        % Smooth and Calculate Acceleration
        %fprintf('Calculating acceleration from the displacement fields.\n')
        [accel,~,~] = func_smoothCalcAccel_v4(pos,time,disp, ...
            smoothOpts.accel,...
            extrapOpts.accel,diffOpts,false);

        % Remove some 3D fields from the structs to save memory.

        % Calculate Stress Gage Stresses

        %fprintf('Calculating Stress Gage Stresses \n')
        SG=func_Full_SG(accel,X_vec,time,Rho);

        %% Save Full Resolution kinematic fields
        fprintf('Recording Kinematic Fields with no edge extrapolation \n')
        NoCorr.disp=disp;
        NoCorr.accel=accel;
        NoCorr.SG=SG;
        NoCorr.strain=strain;

%% Perform Identification
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n')
fprintf('Identifying Parameters from Images Without Edge Extrapolation \n')
fprintf('-------------------------------------------------------------- \n')


%% Condition fields
[SG,~,strain,X_vec,time]=...
func_conditionIBIIDataV2(SG,accel,...
strain,X_vec,time,CondOpts);

%% Run the minimization
[TempG,TempTau,TempPhi]=deal(zeros(3,1));
fprintf('Running the Shear Identification on Imaging Data Without Edge Extrapolation \n')
for ig=1:3
intGuess=squeeze(intMatrix(:,ig)); %#ok<PFBNS>

[TempParam,TempPhi(ig)]=func_PronyShearSGvfm(strain.s, ...
    time.vec,SG.s,exactProps,...
    intGuess,ubG,lbG,SolveOpts,minOpts,CondOpts);
TempG(ig)=TempParam(1);
TempTau(ig)=TempParam(2);


end

minPhi=min(TempPhi);
NoCorr.G=TempG(TempPhi==minPhi);
NoCorr.tauG=TempTau(TempPhi==minPhi);
NoCorr.phi(:)=TempPhi';
NoCorr.ConstitutiveParamG(k,m,:)=[G(1);tau(1)];

fprintf('Running Bulk Identification on Imaging Data Without Edge Extrapolation \n')

%% Create reasonable upper bounds
%bounds on the initial guess are +/- 25%
ubK=[1.875e9;... K1
    12.5E-6]; %tau_1


%% Create reasonable lower bounds
lbK=[1.125e9;...
    7.5E-6]; %tau_1

%% Create Matrix of initial guesses
    %Eint=2.5E9;
    %Choose a random 3x3 matrix of initial guress multipliers. Multipliers
        %will be uniformly distributed random numbers in the interval [0,1]
 
    intMultK=rand(2,3);
      
    %Create 3X3 matrix of initial guess parameters
     intMatrixK=lbK+(ubK-lbK).*intMultK;
   
%     intGuess=[Eint/(3*(1-2*nu));... K1
%               Eint/(2*(1+nu));.... G1
%               10E-6];... tau1
   
%% Initialize bulk identification variables
[TempK,TempPhiK,FETempTau]=deal(zeros([10,1]));

%% Run Bulk Minimization
for ig=1:3
    [TempParam,TempPhi]=func_PronyBulkTauSGvfm(FEstrain, ...
        FEtime,FESGx,exactPropsKFE,...
        intGuess,ubK,lbK,SolveOpts,minOpts,CondOptsFE);
    TempK(ig)=TempParam(1);
    TempTau(ig)=TempParam(2);
    TempPhiK(ig)=TempPhi;
end