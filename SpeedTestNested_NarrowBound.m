%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NOTE: THIS SCRIPT IS ONLY FOR USE WITH GRID IMAGES. USE OTHER CODES FOR
%FINITE ELEMENT VALIDATION.


%This code is written to evaluate the effects of smoothing and noise on the
%identification of viscoelastic constitutitve parameters using the image
%based inertial impact test. 

%Authors: Andrew Matejunas, 
% THIS CODE IS HEAVILY ADAPTED FROM main_IBIIProcessing_v1_0r 
%Date completed: 2022-06-24

%Version History/changelog
 %2022-08-24: Began optimizing code for HPC if needed
 %2022-10-04: V2 Added in post correction of shear strains
 %2022-12-05: V3 Added in 3 initial guesses
              %- Added in post correction of yy strains
              %- Moved plotting of errors to a separate file
              %- Increased spatial downsampling to 6 pixels

 
                    
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

%% Add reference parameters file 
[RefPar.file,RefPar.path]=uigetfile('*.mat',...
    'Choose Mat file containing reference constitutive pamrameters');
RefPar=load(strcat(RefPar.path,'/',RefPar.file),'MatProps');

 %% Create structute of known parameters
    Einf=2.98E9;
    nu=0.26;

    exactProps.Kinf=Einf/(3*(1-2*nu));
    exactProps.Ginf=Einf/(2*(1+nu));
    exactProps.nu=0;

    %% create vector of initial guesses
%     Eint=2.5E9;   
%     intGuess=[Eint/(3*(1-2*nu));... K1
%               Eint/(2*(1+nu));.... G1
%               10E-6];... tau1

    %% Create reasonable upper bounds
    ub=[5E9;...
        2E9;...
        100E-6];



    %% Create reasonable lower bounds
    lb=[0.9E9;...
        100E6;...
        1E-6];

    %% Create Solve Opts Structure
    SolveOpts.constnu=0; %No longer constrain Poisson's ratio
    SolveOpts.KGsame=false; %simultaneously identify K and G
    SolveOpts.identEinf=false;
    SolveOpts.identForm='KG';
    SolveOpts.minFunc='fmincon';

%% Create Matrix of initial guresses
    %Eint=2.5E9;
    %Choose a random 3x3 matrix of initial guress multipliers. Multipliers
        %will be uniformly distributed random numbers in the interval [0,1]
 
    intMult=rand(3);
  
    %Create 3X3 matrix of initial guess parameters

     intMatrix=lb+(ub-lb).*intMult;
%     intGuess=[Eint/(3*(1-2*nu));... K1
%               Eint/(2*(1+nu));.... G1
%               10E-6];... tau1
    %% options for the fmincon algorithm
    minOpts.Algorithm='interior-point';
    %minOpts.Algorithm='sqp';
    minOpts.UseParallel=true;
    minOpts.StepTolerance=1e-3;

%% Choose Directory to save results
savePath=uigetdir({},'Choose Folder to Save Results');

%% Generate bounds of reasonable accuracy (10%)
%K
RefPar.Kp10=RefPar.MatProps.Ki*(1+.1);
RefPar.Km10=RefPar.MatProps.Ki*(1-.1);

%G
RefPar.Gp10=RefPar.MatProps.Gi*(1+.1);
RefPar.Gm10=RefPar.MatProps.Gi*(1-.1);

%tau
RefPar.taup10=RefPar.MatProps.tau*(1+.1);
RefPar.taum10=RefPar.MatProps.tau*(1-.1);

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


quest='Choose Temporal Smoothing Method';
TempOpts=questdlg(quest,'Temporal Method','FF only','FF and WA','None',...
    'FF only');

switch TempOpts
    case 'FF only'
        smoothingOpts.FFTempSmooth=1;
        smoothingOpts.WATempSmooth=0;
        
        % Sweep of Temporal frame lengths
        %TempKernVec=[0,5:2:13];
        TempKernVec=5;
        TempKernNum=length(TempKernVec);
        
    case 'FF and WA'
        smoothingOpts.FFTempSmooth=1;
        smoothingOpts.WATempSmooth=1;
        %TempKernVec=[0,5:2:13];
        TempKernVec=5;
        TempKernNum=length(TempKernVec);
        
    case 'None'
        smoothingOpts.WATempSmooth=0;
        smoothingOpts.FFTempSmooth=0;
end

if smoothingOpts.FFTempSmooth==1
    
    quest='Pad Temporal Filter?';
    TempPadOpts=questdlg(quest,'Temporal Padding','Yes','No','Yes');

    switch TempPadOpts
        case 'Yes'
            smoothingOpts.FFTemporalPad=1;
            smoothingOPts.FFTemporalPadFrames=3;
        case 'No'
            smoothingOpts.FFTemporalPad=0;
    end
end

smoothingOpts.spatialSmooth=0;      
%Spatial Smoothing disabled
%Define sweep of spatial Kernels
%SpaKernVec=[0,5:2:15];
SpaKernVec=5;
%record Number of Kernels
SpaKernNum=length(SpaKernVec);

%% Define whether Noise should be added
%Noise Is hardcoded

% quest='Add Noise?';
% noiseChoice=questdlg(quest,'Noise Choice','Yes');

noiseChoice='Yes';
% 
% switch noiseChoice
%     case 'Yes'
%         %Addition by Andrew Matejunas
%         imageNoise.addNoise=1;    
%         imageNoise.NumCopies=1;
%         %% Choose noise magnitude
%         noiseMag=str2double(inputdlg('Noise Magnitude',...
%             'noise',[1,20],... %Dlgtitle, dims
%             {'0.4000'})); %Default magnitude
%         
%     case 'No'
%         imageNoise.addNoise=0;
% end
    
    imageNoise.addNoise=1;    
    imageNoise.NumCopies=30; 
   
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
DispCorr.strainPitchNum=1;

% switch DispCorr.Opt
%     case 'Yes'
%         quest='How many pixels to correct?';
%         %Correction interverval
%         DispCorr.int=str2num(cell2mat(inputdlg(quest)));
%         clear quest
%         
%         quest='Correction Method';
%         DispCorr.Method=questdlg(quest,'Method',...
%             'Direct',... %sets displacement equal to last valid value
%             'LinGrad',... %follows the average linear gradient of the last
%             ...           %Valid pitch (additional interpolation functions
%             ...               %will be added if needed)
%            'Cancel',...    %Cancels the correction
%            'Direct');    %Defaults to direct method
% 
% 
%     switch DispCorr.Method
%         case 'Cancel'
%             DispCorr.Opt='No';
%     
%         case 'LinGrad'
%         prompt='Number of Grid Pirches to Calculate Gradient';
%         DispCorr.PitchFitKern=str2num(cell2mat(inputdlg(prompt)));
% 
%        
%     end
% end
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
CondOpts.Xds=6;
CondOpts.Yds=1;
CondOpts.Tds=1;


if CondOpts.Tds <=1
    CondOpts.TempDS=false;
else
    CondOpts.TempDS=true;
end

clear prompt dims definputs dlgtitle CondParam

%% Generate Progress bar
itCount=0;
itTotal=TempKernNum*SpaKernNum*imageNoise.NumCopies*3;

progmsg=strcat('0/',num2str(itTotal),...
    ' Smoothing Iterations Complete');
Progress=waitbar(0,progmsg);

%% Run the identification for every copy smoothing Kernal

%Initialize the constitutive model identifications
Ident.K=zeros(TempKernNum,SpaKernNum,imageNoise.NumCopies);
Ident.G=zeros(TempKernNum,SpaKernNum,imageNoise.NumCopies);
Ident.tau=zeros(TempKernNum,SpaKernNum,imageNoise.NumCopies);
phi=zeros(TempKernNum,SpaKernNum,imageNoise.NumCopies);

TempParam=zeros(1,3);
TempK=zeros(1,3);
TempG=zeros(1,3);
TempTau=zeros(1,3);
TempPhi=zeros(3,1);
%Preallocate identified constitutive parameter
 for k=1:TempKernNum
    for m=1:SpaKernNum
        for n=1:imageNoise.NumCopies
               

    if TempKernVec(k)==0
        smoothingOpts.FFTempSmooth=0;
        smoothingOpts.WATempSmooth=0;
    else
        switch TempOpts
        case 'FF only'
            smoothingOpts.FFTempSmooth=1;
            smoothingOpts.WATempSmooth=0;


        case 'FF and WA'
            smoothingOpts.FFTempSmooth=1;
            smoothingOpts.WATempSmooth=1;

        case 'None'
            smoothingOpts.WATempSmooth=0;
            smoothingOpts.FFTempSmooth=0;
        end
    smoothingOpts.FFTemporalKernal=[TempKernVec(k),3];
    smoothingOpts.WATemporalKernal=[TempKernVec(k),3];
    end
    
    if SpaKernVec(m)==0
        smoothingOpts.spatialSmooth=0;
    else
        smoothingOpts.spatialSmooth=1;
        smoothingOpts.spatialKernal=SpaKernVec(m)*[1,1];
    end
    
    fprintf('------------------------------------------------------------- \n')
    fprintf(strcat('Processing iteration',itNum,' \n'))
    fprintf('------------------------------------------------------------- \n')
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

%         fprintf('Checking for existing processed data file.\n')
        processGridImages = true;
        
%             fprintf('Processing images using the grid method toolbox.\n')
            % Process the image squence with the grid method toolbox
            [grid,pos,disp] = func_gridMethodImageProcessing_AJM(imagePath,...
                imageFile,...
                grid,gridMethodOpts,imageNoise);    
%             fprintf('Grid Method Processing Complete.\n')
       

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
            [disp,DispCorr,grid,ProgramVersions]=func_CorrectGMDisp(disp,...
                DispCorr,grid,pos);
            
       


    % Calculate the kinematic fields from displacement fields using
    % displacements from images or displacements from FE data
    if calcKinFieldsFromDisp
        %--------------------------------------------------------------------------
        % Load the Reference Image and Determine Where the Free Edge is
%         fprintf('Obtaining and setting the free edge location.\n')
        [freeEdge,specimen,disp] = func_getFreeEdge(globalOpts.hardCodeFreeEdge,...
            imagePath,imageFile,specimen,disp);

        %--------------------------------------------------------------------------
        % Smooth and Calculate Strain
        fprintf('Calculating strain from the displacement fields.\n')
        [disp,strain,strainRate] = func_smoothCalcStrain(globalOpts,pos,time,...
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
  
   
   %% Run the minimization
            for ig=1:3
                 itCount=itCount+1;
                 itNum=num2str(itCount);
                intGuess=squeeze(intMatrix(:,ig));
    fprintf('Running Parameter Identification \n')
    [TempParam,TempPhi(ig)]=func_PronySGvfmV5(strain,time.vec,SG.x,SG.s,...
                        exactProps,...
                        intGuess,ub,lb,SolveOpts,minOpts,CondOpts);
           TempK(n)=TempParam(1);
           TempG(n)=TempParam(2);
           TempTau(n)=TempParam(3);
   fprintf(strcat('Parameter Identification complete for iteration number',...
       itNum,'\n'))
   
   progmsg=strcat(itNum,'/',num2str(itTotal),...
    ' Smoothing Iterations complete');
   Progress=waitbar(itCount/itTotal,Progress,progmsg);
        progmsg=strcat(itNum,'/',num2str(itTotal),...
    ' Smoothing Iterations complete');
   Progress=waitbar(itCount/itTotal,Progress,progmsg);         
            end
    minPhi=min(TempPhi);
    Ident.K(k,m,n)=TempK(TempPhi==minPhi);
    Ident.G(k,m,n)=TempG(TempPhi==minPhi);
    Ident.Tau(k,m,n)=TempTau(TempPhi==minPhi);

        end
    end
end
close(Progress)

%% Calculate Identification Errors
Errors.K=(Ident.K-...
    RefPar.MatProps.Ki)/RefPar.MatProps.Ki*100;

Errors.G=(Ident.G...
    -RefPar.MatProps.Gi)/RefPar.MatProps.Gi*100;

Errors.tau=(Ident.tau...
    -RefPar.MatProps.tau)/RefPar.MatProps.tau*100;


%% Save Results
 save(strcat(savePath,'/',TestDeg,...
        '_NoisySpeedTestZ_NestedFor','_KGident.mat'),...
        'intGuess','ub','lb','ConstitutiveParam','SolveOpts','imageNoise',...
        'CondOpts','Errors','smoothingOpts','TempKernVec','SpaKernVec',...
        'Ident');
    

