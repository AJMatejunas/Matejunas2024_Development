%This script is written to investigate the causes of noise and smoothing
    %induced systematic error
    
    %Author: Andrew Matejunas
    
    %Date Created: 2022-09-21

    %Change Log:

%% initialize
clear variables; close all; clc

%% Deine Parent Designation
ParentDesig=char(inputdlg('Define Test Designation'));

%% Identify Directory to save analysis data

MainSaveDir=uigetdir('','Select Folder to save analysis');

%% Record Smoothing Vectors
%Decide whether to manually input smoothing vectors or extract from a file
quest='Source of Smoothing Vectors';
SmoothSource=questdlg(quest,'smooth Source','Input','File','File');

switch SmoothSource
    case 'File'
 [SmoothVecName,SmoothVecPath]=uigetfile('*.mat',...
     'Select File Containing Smoothing vectors');
  load(strcat(SmoothVecPath,'/',SmoothVecName));
    case 'Input'
    prompt={'Spatial Smoothing Vector',...
            'Temporal Smoothing Vector'};
        SweepVars=inputdlg(prompt);
        
        SpaKernVec=eval(char(SweepVars{1}));
        SpaSmoothChar=num2str(SpaKernVec');

        TempKernVec=eval(char(SweepVars{2}));
        TempSmoothChar=num2str(TempKernVec');

       save(strcat(MainSaveDir,'\',ParentDesig,'_SmoothingVectors.mat'))
end


%% Define displacement correction options
%Displacement corrections are hardcoded
fprintf('Hardcoding Displacement Correction \n')
DispCorr.Opt='Yes';
DispCorr.int=10;
DispCorr.Method='LinGrad';
DispCorr.PitchFitKern=2;

%% Hardcode Sample conditioning options
fprintf('Hardcoding Data conditoning \n')
CondOpts.ImpCens=15;
CondOpts.FreeCens=20;
CondOpts.Xds=3;
CondOpts.Yds=1;
CondOpts.Tds=1;

%% Identify files containing the Noise and NoNoise analysis

% Selcect file without noise
[NoNoiseIdent.FileName,NoNoiseIdent.FilePath]=uigetfile('*.mat', ...
    'Choose file Containing Identified Parameters without noise');
NoNoiseIdent.FullFile=strcat(NoNoiseIdent.FilePath,'\',NoNoiseIdent.FileName);

% Select file with noise
[NoiseIdent.FileName,NoiseIdent.FilePath]=uigetfile('*.mat', ...
    'Choose File Containing Identified Parameters With Noise');
NoiseIdent.FullFile=strcat(NoiseIdent.FilePath,'\',NoiseIdent.FileName);

%% Identify FE data file
[FE.DataName,FE.DataPath]=uigetfile('*.mat', ...
    'Choose Raw FE SG data file');
FE.FullFile=strcat(FE.DataPath,'\',FE.DataName);

%% Identify First Image Frame
 [imageFile,imagePath] = uigetfile({'*.*','All Files'},...
     'Select the first image in the sequence');
%% identify Processing parameters file
[ProcessParam.Name,ProcessParam.Path]=uigetfile('*.mat',...
    'Choose Process Parameter Data File');
ProcessParam.File=strcat(ProcessParam.Path,'/',ProcessParam.Name);
fprintf('Loading Process Parameter Data File \n')
load(ProcessParam.File);

%% Load indentified parameters files
fprintf('Loading identifications without noise \n')
load(NoNoiseIdent.FullFile);
fprintf('Loading identifications with noise \n')
load(NoiseIdent.FullFile);

%% Calculate Kinematic Fields on single Copy
quest='Source of Kinematic Fields';
KinFieldSource=questdlg(quest,'Kinematic Field data source','Files', ...
    'images','Files');
switch KinFieldSource
    case 'Files'
        %% Choose Files Containing Kinematic Fields
        [Raw.Name,Raw.Path]=uigetfile('*.mat',...
            'Choose File Containing Kinematic Fields');
        Raw.File=strcat(Raw.Path,'\',Raw.Name);
        %% noisy
        [NoiseInfo.Name,NoiseInfo.Path]=uigetfile('.mat', ...
            'Choose file containing fields with noise');
        NoiseInfo.File=strcat(NoiseInfo.Path,'\',NoiseInfo.Name);

        %% Load
        fprintf('Loading Kinematic Field Without Noise \n')
        RawFields=load(Raw.File);
        RawFields.File=Raw.File;
        clear Raw

        fprintf('Loading Kinematic Fields With Noise \n')
        NoiseFields=load(NoiseInfo.File);
        NoiseFields.File=NoiseInfo.File;
        clear NoiseInfo

        fprintf('Kinematic Fields Loaded \n')
    case 'images'
    %% Calculate displacemnent, strain, acceleration without noise or smoothing
    fprintf('Calculating Kinematic Fields Without Noise Or Smoothing \n')
    
   NoiseCall='NoNoise';
   imageNoise.addNoise=false;
   imageNoise.bits=16;
   imageNoise.pcNoise=0.4;
   imageNoise.convToUInt16=true;
   
   %Define Smoothing options
   smoothingOpts.spatialSmooth=false;
   smoothingOpts.spatialKernal=[0,0];
   smoothingOpts.FFTempSmooth=false;
   smoothingOpts.FFTemporalKernal=[0,3];

    RawFields= func_CalcSaveKinFieldsFromImagesNoSmooth(MainSaveDir,...
    imageFile,imagePath,ParentDesig,NoiseCall, ...
        grid,specimen,time, ...
        gridMethodOpts,imageNoise,DispCorr,globalOpts,extrapOpts,...
        smoothingOpts,diffOpts);

    clear NoiseCall
    RawFields.smoothingOpts=smoothingOpts;
    fprintf('Calculating Noisy Fields without smoothing \n')
    NoiseCall='WithNoise1copy';
    imageNoise.addNoise=true;
    NoiseFields= func_CalcSaveKinFieldsFromImagesNoSmooth(MainSaveDir,...
    imageFile,imagePath,ParentDesig,NoiseCall, ...
        grid,specimen,time, ...
        gridMethodOpts,imageNoise,DispCorr,globalOpts,extrapOpts,...
        smoothingOpts,diffOpts);

    

        
end

%% Obtain max smoothing strain, and acceleration fields
quest='Source of Maximum Smoothing Fields \n';
MaxSmoothSource=questdlg(quest,'Max Smoothing data source','Files', ...
    'images','Files');

switch MaxSmoothSource
    case 'Files'
        [MaxSmooth.Name,MaxSmooth.Path]=uigetfile('*.mat', ...
   'choose File Containing No Noise Data Smoothed with the Maximum kernal');
        MaxSmooth.File=strcat(MaxSmooth.Path,'\',MaxSmooth.Name);
        

        [MaxSmoothNoise.Name,MaxSmoothNoise.Path]=uigetfile('*.mat', ...
          'Choose File Containing Noisy Data smoothed with Maximum Kernal');
        MaxSmoothNoise.File=strcat(MaxSmoothNoise.Path,'/', ...
            MaxSmoothNoise.Name);
        %% Load
        fprintf('Loading Max smooth without noisy \n')
        MaxSmooth=load(MaxSmooth.File);
        fprintf('Loading Max Noisy max smoothed \n')
        MaxSmoothNoise=load(MaxSmoothNoise.File);
    
    case 'images'
        %% Set Smoothing Options
        smoothingOpts.spatialSmooth=true;
        smoothingOpts.FFTempSmooth=true;
        smoothingOpts.WATempSmooth=false;
        smoothingOpts.FFTemporalPad=1;
        smoothingOpts.FFTemporalPadFrames=3;
        smoothingOpts.spatialFilt='gauss';
        smoothingOpts.spatialEdgeMode='symmetric';
        smoothingOpts.FFTemporalFilt='sgolay';
        smoothingOpts.FFTemporalPadMethod='replicate';


        smoothingOpts.spatialKernal=SpaKernVec(end)*[1,1];
        smoothingOpts.FFTemporalKernal=[TempKernVec(end),3];

       %% Calculate and smooth strain and acceeleration fileds No Noise
       fprintf('Smoothing No Noise Fields with Max Kernal \n')
       SmoothLevel='Max';
       NoiseCall='NoNoise';
       
       MaxSmooth=func_SmoothSaveKinFields(RawFields.globalOpts, ...
            RawFields.extrapOpts,smoothingOpts,RawFields.diffOpts, ...
            RawFields.pos,RawFields.time,RawFields.grid,RawFields.disp, ...
            MainSaveDir,ParentDesig,NoiseCall,SmoothLevel, ...
            RawFields.DispCorr,RawFields.specimen,RawFields.imageNoise);
      
      %% Calculate and smooth strain and accelerations fields with noise
      fprintf('Smoothing Noisy Fields with Max kernal \n')
      NoiseCall='withNoise1copy';
      
      MaxSmoothNoise=func_SmoothSaveKinFields(NoiseFields.globalOpts, ...
            NoiseFields.extrapOpts,smoothingOpts,NoiseFields.diffOpts, ...
            NoiseFields.pos,NoiseFields.time,NoiseFields.grid, ...
            NoiseFields.disp, ...
            MainSaveDir,ParentDesig,NoiseCall,SmoothLevel, ...
            NoiseFields.DispCorr,NoiseFields.specimen, ...
            NoiseFields.imageNoise);



end

%% Obtain Kinematic Fields with the optimal smoothing kernal
quest='Source of Optimal Smoothing Fields \n';
OptSmoothSource=questdlg(quest,'Opt Smooth data source source','Files', ...
    'images','Files');

switch OptSmoothSource
    case 'Files'
    [OptSmooth.Name,OptSmooth.Path]=uigetfile('*.mat', ...
        'Choose file containing optimally smoothed fields');
    OptSmooth.File=strcat(OptSmooth.Path,'/',OptSmooth.Name);
   

    [OptSmoothNoise.Name,OptSmoothNoise.Path]=uigetfile('*.mat',...
        'Choose file containing optimally smoothed noisy data');
    OptSmoothNoise.File=strcat(OptSmoothNoise.Path,'/', ...
        OptSmoothNoise.Name);
    
    %% Load
    fprintf('Loading Optimal smoothed No Noise \n')
    OptSmooth=load(OptSmooth.File);
    fprintf('Loading Optimal Smoothed Noisy \n')
    OptSmoothNoise=load(OptSmoothNoise.File);

    case 'images'
        %% Input Optimum Smoothing Kernal
        prompt={'Input Optimum Spatial Smoothing Kernal','Temporal Kernal'};
        OptKerns=inputdlg(prompt,'Define Optimum Smoothing',[1,50;1,50]);
        OptSpa=str2double(OptKerns{1});
        OptTemp=str2double(OptKerns{2});
        
        
        %% Set smoothing Options
        smoothingOpts.spatialKernal=OptSpa*[1,1];
        smoothingOpts.FFTemporalKernal=[OptTemp,3];

        fprintf('Smoothing Raw Fields with Optimal Smoothing Kernal \n')
        NoiseCall='NoNoise';
        SmoothLevel='Opt';
        OptSmooth=func_SmoothSaveKinFields(RawFields.globalOpts, ...
            RawFields.extrapOpts,smoothingOpts,RawFields.diffOpts, ...
            RawFields.pos,RawFields.time,RawFields.grid,RawFields.disp, ...
            MainSaveDir,ParentDesig,NoiseCall,SmoothLevel, ...
            RawFields.DispCorr,RawFields.specimen,RawFields.imageNoise);

        fprintf('Smoothing Noisy Fields with Optimal Smdoothing Kernal \n')
        NoiseCall='WithNoise1copy';
        OptSmoothNoise=func_SmoothSaveKinFields(NoiseFields.globalOpts, ...
            NoiseFields.extrapOpts,smoothingOpts,NoiseFields.diffOpts, ...
            NoiseFields.pos,NoiseFields.time,NoiseFields.grid, ...
            NoiseFields.disp, ...
            MainSaveDir,ParentDesig,NoiseCall,SmoothLevel, ...
            NoiseFields.DispCorr,NoiseFields.specimen, ...
            NoiseFields.imageNoise);      
end

%% Load FE data
fprintf('Loading Finite Element Reference Data \n')
FEfields=load(FE.FullFile);
FEfields.File=FE.FullFile;

%%  Generate Strain Time Plots

fprintf('Generating Strain-Time Diagnostic Plots With No Noise \n')
    % Set Position Indexes for Grid image data
    GridPos.NumPoints=size(RawFields.disp.x,2);
    GridPos.Free=20; %4 pitches from free edge
    GridPos.Mid=round(GridPos.NumPoints/2);
    GridPos.Imp=GridPos.NumPoints-20; %4 pitches from impact edge
    % Get cooridinates of position indexes
    GridPos.Xfree=RawFields.pos.x(GridPos.Free);
    GridPos.Xmid=RawFields.pos.x(GridPos.Mid);
    GridPos.Ximp=RawFields.pos.x(GridPos.Imp);

    %Set Position indexes and coordintes for finite element data
    FEPos.NumPoints=size(FEfields.disp.x,2);
    FEPos.GridRat=FEPos.NumPoints/GridPos.NumPoints;
    FEPos.Free=round(GridPos.Free*FEPos.GridRat); %4 pitches from free edge
    FEPos.Mid=round(FEPos.NumPoints/2);
    FEPos.Imp=round(GridPos.Imp*FEPos.GridRat); %4 pitches from impact edge
    
    FEPos.Xfree=FEfields.pos.x(FEPos.Free);
    FEPos.Xmid=FEfields.pos.x(FEPos.Mid);
    FEPos.Ximp=FEfields.pos.x(FEPos.Imp);
    
  %Print Plot selections for evaluation purposes
  FEPos.freeString=num2str(FEPos.Xfree);
  GridPos.freeString=num2str(GridPos.Xfree);
  fprintf(strcat('Free Coordinate FE:',FEPos.freeString,' Grid: ',...
      GridPos.freeString,'\n'));
  
  FEPos.midString=num2str(FEPos.Xmid);
  GridPos.midString=num2str(GridPos.Xmid);
  fprintf(strcat('Mid Coordinate FE:',FEPos.midString,' Grid: ',...
      GridPos.midString,'\n'));

  FEPos.impString=num2str(FEPos.Ximp);
  GridPos.impString=num2str(GridPos.Ximp);
  fprintf(strcat('Impact Coordinate FE:',FEPos.impString,' Grid: ',...
      GridPos.impString,'\n'));
  
    %Generate plot variables for easier coding
        %Finite Element
        FETime=FEfields.time.vec*10^6;
        FEfields.avgX_strain=squeeze(FEfields.avgX_strain);
        FEfields.avgXY_strain=squeeze(FEfields.avgXY_strain);
        FEFreeStrainX=squeeze(FEfields.avgX_strain(FEPos.Free,:));
        FEimpStrainX=squeeze(FEfields.avgX_strain(FEPos.Imp,:));
        FEmidStrainX=squeeze(FEfields.avgX_strain(FEPos.Mid,:));
    
        FEFreeStrainS=squeeze(FEfields.avgXY_strain(FEPos.Free,:));
        FEimpStrainS=squeeze(FEfields.avgXY_strain(FEPos.Imp,:));
        FEmidStrainS=squeeze(FEfields.avgXY_strain(FEPos.Mid,:));

        %Grid No Smoothing
        GridTime=RawFields.time.vec*10^6;
        RawFreeStrainX=squeeze(RawFields.avgX_strain(GridPos.Free,:));
        RawImpStrainX=squeeze(RawFields.avgX_strain(GridPos.Imp,:));
        RawMidStrainX=squeeze(RawFields.avgX_strain(GridPos.Mid,:));

        RawFreeStrainS=squeeze(RawFields.avgXY_strain(GridPos.Free,:));
        RawImpStrainS=squeeze(RawFields.avgXY_strain(GridPos.Imp,:));
        RawMidStrainS=squeeze(RawFields.avgXY_strain(GridPos.Mid,:));

        %Grid Optimal Smoothing
        OptFreeStrainX=squeeze(OptSmooth.avgX_strain(GridPos.Free,:));
        OptImpStrainX=squeeze(OptSmooth.avgX_strain(GridPos.Imp,:));
        OptMidStrainX=squeeze(OptSmooth.avgX_strain(GridPos.Mid,:));

        OptFreeStrainS=squeeze(OptSmooth.avgXY_strain(GridPos.Free,:));
        OptImpStrainS=squeeze(OptSmooth.avgXY_strain(GridPos.Imp,:));
        OptMidStrainS=squeeze(OptSmooth.avgXY_strain(GridPos.Mid,:));

        %Grid Maximum Smoothing
        MaxFreeStrainX=squeeze(MaxSmooth.avgX_strain(GridPos.Free,:));
        MaxImpStrainX=squeeze(MaxSmooth.avgX_strain(GridPos.Imp,:));
        MaxMidStrainX=squeeze(MaxSmooth.avgX_strain(GridPos.Mid,:));

        MaxFreeStrainS=squeeze(MaxSmooth.avgXY_strain(GridPos.Free,:));
        MaxImpStrainS=squeeze(MaxSmooth.avgXY_strain(GridPos.Imp,:));
        MaxMidStrainS=squeeze(MaxSmooth.avgXY_strain(GridPos.Mid,:));


%Set figure size to full screen for diagnostic purposes
figure('units','normalized','outerposition',[0 0 1 1])
title('Systematic Error No Noise')
% X strain impact edge
subplot(2,3,1)
%FE
plot(FETime,FEimpStrainX,'--k')
hold on
%No Smooth
plot(GridTime,RawImpStrainX,'k')
%OptSmooth
plot(GridTime,OptImpStrainX,'b')
%Max Smooth
plot(GridTime,MaxImpStrainX,'r')

hold off
title('4P from Impact')
xlabel('Time ($\mu$s)','Interpreter','latex')
ylabel('$\varepsilon_{xx}$','Interpreter','latex')
legend('FE','No Smoothing','Optimum Smoothing','Max Smoothing')
% X strain Middle 
subplot(2,3,2)
%FE
plot(FETime,FEmidStrainX,'--k')
hold on
%No Smooth
plot(GridTime,RawMidStrainX,'k')
%OptSmooth
plot(GridTime,OptMidStrainX,'b')
%Max Smooth
plot(GridTime,MaxMidStrainX,'r')

hold off

title('Specimen Middle')
xlabel('Time ($\mu$s)','Interpreter','latex')
ylabel('$\varepsilon_{xx}$','Interpreter','latex')
% X Strain Free edge
subplot(2,3,3)
%FE
plot(FETime,FEFreeStrainX,'--k')
hold on
%No Smooth
plot(GridTime,RawFreeStrainX,'k')
%OptSmooth
plot(GridTime,OptFreeStrainX,'b')
%Max Smooth
plot(GridTime,MaxFreeStrainX,'r')

hold off

title('4P from Free Edge')
xlabel('Time ($\mu$s)','Interpreter','latex')
ylabel('$\varepsilon_{xx}$','Interpreter','latex')


% Shear Strain Impact
subplot(2,3,4)
%FE
plot(FETime,FEimpStrainS,'--k')
hold on
%No Smooth
plot(GridTime,RawImpStrainS,'k')
%OptSmooth
plot(GridTime,OptImpStrainS,'b')
%Max Smooth
plot(GridTime,MaxImpStrainS,'r')

hold off

title('4P from Impact')
xlabel('Time ($\mu$s)','Interpreter','latex')
ylabel('$\varepsilon_{xy}$','Interpreter','latex')

% Shear Strain Middle
subplot(2,3,5)
%FE
plot(FETime,FEmidStrainS,'--k')
hold on
%No Smooth
plot(GridTime,RawMidStrainS,'k')
%OptSmooth
plot(GridTime,OptMidStrainS,'b')
%Max Smooth
plot(GridTime,MaxMidStrainS,'r')

hold off
title('Specimen Middle')
xlabel('Time ($\mu$s)','Interpreter','latex')
ylabel('$\varepsilon_{xy}$','Interpreter','latex')

% Shear strain free edge
subplot(2,3,6)
%FE
plot(FETime,FEFreeStrainS,'--k')
hold on
%No Smooth
plot(GridTime,RawFreeStrainS,'k')
%OptSmooth
plot(GridTime,OptFreeStrainS,'b')
%Max Smooth
plot(GridTime,MaxFreeStrainS,'r')


hold off
title('Free Edge')
xlabel('Time ($\mu$s)','Interpreter','latex')
ylabel('$\varepsilon_{xy}$','Interpreter','latex')

sgtitle('Systematic Error No Noise')

STSaveName=strcat(MainSaveDir,'/',ParentDesig, ...
    '_SystematicErrorNoNoise_StrainTime');
saveas(gcf,strcat(STSaveName,'.fig'))
saveas(gcf,strcat(STSaveName,'.svg'))

%% Generate strain time plots with 1 Noise copy
 
        %Noisy No Smoothing
        GridTime=NoiseFields.time.vec*10^6;
        NoiseFreeStrainX=squeeze(NoiseFields.avgX_strain(GridPos.Free,:));
        NoiseImpStrainX=squeeze(NoiseFields.avgX_strain(GridPos.Imp,:));
        NoiseMidStrainX=squeeze(NoiseFields.avgX_strain(GridPos.Mid,:));

        NoiseFreeStrainS=squeeze(NoiseFields.avgXY_strain(GridPos.Free,:));
        NoiseImpStrainS=squeeze(NoiseFields.avgXY_strain(GridPos.Imp,:));
        NoiseMidStrainS=squeeze(NoiseFields.avgXY_strain(GridPos.Mid,:));

        %Noisy Optimal Smoothing
        OptFreeStrainX=squeeze(OptSmoothNoise.avgX_strain(GridPos.Free,:));
        OptImpStrainX=squeeze(OptSmoothNoise.avgX_strain(GridPos.Imp,:));
        OptMidStrainX=squeeze(OptSmoothNoise.avgX_strain(GridPos.Mid,:));

        OptFreeStrainS=squeeze(OptSmoothNoise.avgXY_strain(GridPos.Free,:));
        OptImpStrainS=squeeze(OptSmoothNoise.avgXY_strain(GridPos.Imp,:));
        OptMidStrainS=squeeze(OptSmoothNoise.avgXY_strain(GridPos.Mid,:));

        %Grid Maximum Smoothing
        MaxFreeStrainX=squeeze(MaxSmoothNoise.avgX_strain(GridPos.Free,:));
        MaxImpStrainX=squeeze(MaxSmoothNoise.avgX_strain(GridPos.Imp,:));
        MaxMidStrainX=squeeze(MaxSmoothNoise.avgX_strain(GridPos.Mid,:));

        MaxFreeStrainS=squeeze(MaxSmoothNoise.avgXY_strain(GridPos.Free,:));
        MaxImpStrainS=squeeze(MaxSmoothNoise.avgXY_strain(GridPos.Imp,:));
        MaxMidStrainS=squeeze(MaxSmoothNoise.avgXY_strain(GridPos.Mid,:));


    %Set figure size to full screen for diagnostic purposes
    figure('units','normalized','outerposition',[0 0 1 1])
    % X strain impact edge
    subplot(2,3,1)
    %FE
    plot(FETime,FEimpStrainX,'--k')
    hold on
    %No Smooth
    plot(GridTime,NoiseImpStrainX,'k')
    %OptSmooth
    plot(GridTime,OptImpStrainX,'b')
    %Max Smooth
    plot(GridTime,MaxImpStrainX,'r')

    hold off
    title('4P from Impact')
    xlabel('Time ($\mu$s)','Interpreter','latex')
    ylabel('$\varepsilon_{xx}$','Interpreter','latex')
    legend('FE','No Smoothing','Optimum Smoothing','Max Smoothing')
    % X strain Middle 
    subplot(2,3,2)
    %FE
    plot(FETime,FEmidStrainX,'--k')
    hold on
    %No Smooth
    plot(GridTime,NoiseMidStrainX,'k')
    %OptSmooth
    plot(GridTime,OptMidStrainX,'b')
    %Max Smooth
    plot(GridTime,MaxMidStrainX,'r')

    hold off

    title('Specimen Middle')
    xlabel('Time ($\mu$s)','Interpreter','latex')
    ylabel('$\varepsilon_{xx}$','Interpreter','latex')
    % X Strain Free edge
    subplot(2,3,3)
    %FE
    plot(FETime,FEFreeStrainX,'--k')
    hold on
    %No Smooth
    plot(GridTime,NoiseFreeStrainX,'k')
    %OptSmooth
    plot(GridTime,OptFreeStrainX,'b')
    %Max Smooth
    plot(GridTime,MaxFreeStrainX,'r')

    hold off

    title('4P from Free Edge')
    xlabel('Time ($\mu$s)','Interpreter','latex')
    ylabel('$\varepsilon_{xx}$','Interpreter','latex')


    % Shear Strain Impact
    subplot(2,3,4)
    %FE
    plot(FETime,FEimpStrainS,'--k')
    hold on
    %No Smooth
    plot(GridTime,NoiseImpStrainS,'k')
    %OptSmooth
    plot(GridTime,OptImpStrainS,'b')
    %Max Smooth
    plot(GridTime,MaxImpStrainS,'r')

    hold off

    title('4P from Impact')
    xlabel('Time ($\mu$s)','Interpreter','latex')
    ylabel('$\varepsilon_{xy}$','Interpreter','latex')

    % Shear Strain Middle
    subplot(2,3,5)
    %FE
    plot(FETime,FEmidStrainS,'--k')
    hold on
    %No Smooth
    plot(GridTime,NoiseMidStrainS,'k')
    %OptSmooth
    plot(GridTime,OptMidStrainS,'b')
    %Max Smooth
    plot(GridTime,MaxMidStrainS,'r')

    hold off
    title('Specimen Middle')
    xlabel('Time ($\mu$s)','Interpreter','latex')
    ylabel('$\varepsilon_{xy}$','Interpreter','latex')

    % Shear strain free edge
    subplot(2,3,6)
    %FE
    plot(FETime,FEFreeStrainS,'--k')
    hold on
    %No Smooth
    plot(GridTime,NoiseFreeStrainS,'k')
    %OptSmooth
    plot(GridTime,OptFreeStrainS,'b')
    %Max Smooth
    plot(GridTime,MaxFreeStrainS,'r')


    hold off
    title('Free Edge')
    xlabel('Time ($\mu$s)','Interpreter','latex')
    ylabel('$\varepsilon_{xy}$','Interpreter','latex')

    sgtitle('1 Noise Copy')
    STSaveName=strcat(MainSaveDir,'/',ParentDesig, ...
        '_SystematicError1Noise_StrainTime');
    saveas(gcf,strcat(STSaveName,'.fig'))
    saveas(gcf,strcat(STSaveName,'.svg'))

%% Calculate Stress Gauge Stresses
switch KinFieldSource
    case 'images'
       %% if smoothing is calculated from images recalculate stress gauge 
        SGsource='Calculate';

    case 'Files'
        %%  If smoothing is taken from files determine whether stress gauge
            %stresses must be calculated
        SGquest='Calculate Stress Gage?';
        SGsource=questdlg(SGquest,'SG source','Calculate','File','Files');
end

switch SGsource
    case 'Calculate'
        %% Add material field to all kinematic field structures
        RawFields.material=material;
        NoiseFields.material=material;
        MaxSmooth.material=material;
        MaxSmoothNoise.material=material;
        OptSmooth.material=material;
        OptSmoothNoise.material=material;

        %% Calculate stress gauge stresses for NonSmoothed data
        %No Noise
        fprintf('Calculating stress gauge stresses without noise or smoothing \n')
        RawFields.SG=func_Full_SG(RawFields.accel,RawFields.pos.x,...
            RawFields.time,material.rho);
        
        % With noise
        fprintf('Calculating Noisy SG without smoothing \n')
        NoiseFields.SG=func_Full_SG(NoiseFields.accel,NoiseFields.pos.x,...
            NoiseFields.time,material.rho);
        fprintf('Saving Unsmoothed stress gauge \n')
        save(RawFields.File,'-struct','RawFields') 
        save(NoiseFields.File,'-struct','NoiseFields')

        %% SG for Maximum smoothing
         %No Noise
        fprintf('Calculating No Noise SG max smoothing \n')
        MaxSmooth.SG=func_Full_SG(MaxSmooth.accel,MaxSmooth.pos.x,...
            MaxSmooth.time,material.rho);
        
        % With noise
        fprintf('Calculating Noisy SG Max smoothing \n')
        MaxSmoothNoise.SG=func_Full_SG(MaxSmoothNoise.accel,MaxSmoothNoise.pos.x,...
            MaxSmoothNoise.time,material.rho);
        fprintf('Saving Max Smoothed stress gauge \n')
        save(MaxSmooth.File,'-struct','MaxSmooth') 
        save(MaxSmoothNoise.File,'-struct','MaxSmoothNoise')

        %% SG for Optimum smoothing
         %No Noise
        fprintf('Calculating No Noise SG Optimum smoothing \n')
        OptSmooth.SG=func_Full_SG(OptSmooth.accel,OptSmooth.pos.x,...
            OptSmooth.time,material.rho);
        
        % With noise
        fprintf('Calculating Noisy SG Opt smoothing \n')
        OptSmoothNoise.SG=func_Full_SG(OptSmoothNoise.accel,OptSmoothNoise.pos.x,...
            OptSmoothNoise.time,material.rho);
        fprintf('Saving Opt Smoothed stress gauge \n')
        save(OptSmooth.File,'-struct','OptSmooth') 
        save(OptSmoothNoise.File,'-struct','OptSmoothNoise')


end

%% Plot SG time Plots Without Noise
fprintf('Generating SG-Time Diagnostic Plots With No Noise \n')
%Reformat SG in FEfields

FEfields=rmfield(FEfields,'SG');
FEfields.SG.x=FEfields.Full_SG.x;
FEfields.SG.s=FEfields.Full_SG.s;
    
   
    %Generate plot variables for easier coding
        %Finite Element

        FEFreeSGX=squeeze(FEfields.SG.x(FEPos.Free,:))*10^-6;
        FEimpSGX=squeeze(FEfields.SG.x(FEPos.Imp,:))*10^-6;
        FEmidSGX=squeeze(FEfields.SG.x(FEPos.Mid,:))*10^-6;
    
        FEFreeSGS=squeeze(FEfields.SG.s(FEPos.Free,:))*10^-6;
        FEimpSGS=squeeze(FEfields.SG.s(FEPos.Imp,:))*10^-6;
        FEmidSGS=squeeze(FEfields.SG.s(FEPos.Mid,:))*10^-6;

        %Grid No Smoothing
        GridTime=RawFields.time.vec*10^6;
        RawFreeSGX=squeeze(RawFields.SG.x(GridPos.Free,:))*10^-6;
        RawImpSGX=squeeze(RawFields.SG.x(GridPos.Imp,:))*10^-6;
        RawMidSGX=squeeze(RawFields.SG.x(GridPos.Mid,:))*10^-6;

        RawFreeSGS=squeeze(RawFields.SG.s(GridPos.Free,:))*10^-6;
        RawImpSGS=squeeze(RawFields.SG.s(GridPos.Imp,:))*10^-6;
        RawMidSGS=squeeze(RawFields.SG.s(GridPos.Mid,:))*10^-6;

        %Grid Optimal Smoothing
        OptFreeSGX=squeeze(OptSmooth.SG.x(GridPos.Free,:))*10^-6;
        OptImpSGX=squeeze(OptSmooth.SG.x(GridPos.Imp,:))*10^-6;
        OptMidSGX=squeeze(OptSmooth.SG.x(GridPos.Mid,:))*10^-6;

        OptFreeSGS=squeeze(OptSmooth.SG.s(GridPos.Free,:))*10^-6;
        OptImpSGS=squeeze(OptSmooth.SG.s(GridPos.Imp,:))*10^-6;
        OptMidSGS=squeeze(OptSmooth.SG.s(GridPos.Mid,:))*10^-6;

        %Grid Maximum Smoothing
        MaxFreeSGX=squeeze(MaxSmooth.SG.x(GridPos.Free,:))*10^-6;
        MaxImpSGX=squeeze(MaxSmooth.SG.x(GridPos.Imp,:))*10^-6;
        MaxMidSGX=squeeze(MaxSmooth.SG.x(GridPos.Mid,:))*10^-6;

        MaxFreeSGS=squeeze(MaxSmooth.SG.s(GridPos.Free,:))*10^-6;
        MaxImpSGS=squeeze(MaxSmooth.SG.s(GridPos.Imp,:))*10^-6;
        MaxMidSGS=squeeze(MaxSmooth.SG.s(GridPos.Mid,:))*10^-6;


%Set figure size to full screen for diagnostic purposes
figure('units','normalized','outerposition',[0 0 1 1])
title('Systematic Error No Noise')
% X SG impact edge
subplot(2,3,1)
%FE
plot(FETime,FEimpSGX,'--k')
hold on
%No Smooth
plot(GridTime,RawImpSGX,'k')
%OptSmooth
plot(GridTime,OptImpSGX,'b')
%Max Smooth
plot(GridTime,MaxImpSGX,'r')

hold off
title('4P from Impact')
xlabel('Time ($\mu$s)','Interpreter','latex')
ylabel('$\sigma_{xx}$ (MPa)','Interpreter','latex')
legend('FE','No Smoothing','Optimum Smoothing','Max Smoothing')
% X SG Middle 
subplot(2,3,2)
%FE
plot(FETime,FEmidSGX,'--k')
hold on
%No Smooth
plot(GridTime,RawMidSGX,'k')
%OptSmooth
plot(GridTime,OptMidSGX,'b')
%Max Smooth
plot(GridTime,MaxMidSGX,'r')

hold off

title('Specimen Middle')
xlabel('Time ($\mu$s)','Interpreter','latex')
ylabel('$\sigma_{xx}$ (MPa)','Interpreter','latex')
% X SG Free edge
subplot(2,3,3)
%FE
plot(FETime,FEFreeSGX,'--k')
hold on
%No Smooth
plot(GridTime,RawFreeSGX,'k')
%OptSmooth
plot(GridTime,OptFreeSGX,'b')
%Max Smooth
plot(GridTime,MaxFreeSGX,'r')

hold off

title('4P from Free Edge')
xlabel('Time ($\mu$s)','Interpreter','latex')
ylabel('$\sigma_{xx}$ (MPa)','Interpreter','latex')


% Shear SG Impact
subplot(2,3,4)
%FE
plot(FETime,FEimpSGS,'--k')
hold on
%No Smooth
plot(GridTime,RawImpSGS,'k')
%OptSmooth
plot(GridTime,OptImpSGS,'b')
%Max Smooth
plot(GridTime,MaxImpSGS,'r')

hold off

title('4P from Impact')
xlabel('Time ($\mu$s)','Interpreter','latex')
ylabel('$\sigma_{xy}$ (MPa)','Interpreter','latex')

% Shear SG Middle
subplot(2,3,5)
%FE
plot(FETime,FEmidSGS,'--k')
hold on
%No Smooth
plot(GridTime,RawMidSGS,'k')
%OptSmooth
plot(GridTime,OptMidSGS,'b')
%Max Smooth
plot(GridTime,MaxMidSGS,'r')

hold off
title('Specimen Middle')
xlabel('Time ($\mu$s)','Interpreter','latex')
ylabel('$\sigma_{xy}$ (MPa)','Interpreter','latex')

% Shear SG free edge
subplot(2,3,6)
%FE
plot(FETime,FEFreeSGS,'--k')
hold on
%No Smooth
plot(GridTime,RawFreeSGS,'k')
%OptSmooth
plot(GridTime,OptFreeSGS,'b')
%Max Smooth
plot(GridTime,MaxFreeSGS,'r')


hold off
title('Free Edge')
xlabel('Time ($\mu$s)','Interpreter','latex')
ylabel('$\sigma_{xy}$ (MPa)','Interpreter','latex')

sgtitle('Systematic Error No Noise')

STSaveName=strcat(MainSaveDir,'/',ParentDesig, ...
    '_SystematicErrorNoNoise_SGTime');
saveas(gcf,strcat(STSaveName,'.fig'))
saveas(gcf,strcat(STSaveName,'.svg'))


%% Plot SG time Plots With Noise
fprintf('Generating SG-Time Diagnostic Plots With added Noise \n')
%Reformat SG in FEfields

FEfields=rmfield(FEfields,'SG');
FEfields.SG.x=FEfields.Full_SG.x;
FEfields.SG.s=FEfields.Full_SG.s;
    
   
    %Generate plot variables for easier coding
        %Finite Element

        FEFreeSGX=squeeze(FEfields.SG.x(FEPos.Free,:))*10^-6;
        FEimpSGX=squeeze(FEfields.SG.x(FEPos.Imp,:))*10^-6;
        FEmidSGX=squeeze(FEfields.SG.x(FEPos.Mid,:))*10^-6;
    
        FEFreeSGS=squeeze(FEfields.SG.s(FEPos.Free,:))*10^-6;
        FEimpSGS=squeeze(FEfields.SG.s(FEPos.Imp,:))*10^-6;
        FEmidSGS=squeeze(FEfields.SG.s(FEPos.Mid,:))*10^-6;

        %Grid No Smoothing
        GridTime=RawFields.time.vec*10^6;
        NoiseFreeSGX=squeeze(NoiseFields.SG.x(GridPos.Free,:))*10^-6;
        NoiseImpSGX=squeeze(NoiseFields.SG.x(GridPos.Imp,:))*10^-6;
        NoiseMidSGX=squeeze(NoiseFields.SG.x(GridPos.Mid,:))*10^-6;

        NoiseFreeSGS=squeeze(NoiseFields.SG.s(GridPos.Free,:))*10^-6;
        NoiseImpSGS=squeeze(NoiseFields.SG.s(GridPos.Imp,:))*10^-6;
        NoiseMidSGS=squeeze(NoiseFields.SG.s(GridPos.Mid,:))*10^-6;

        %Grid Optimal Smoothing
        OptFreeSGX=squeeze(OptSmoothNoise.SG.x(GridPos.Free,:))*10^-6;
        OptImpSGX=squeeze(OptSmoothNoise.SG.x(GridPos.Imp,:))*10^-6;
        OptMidSGX=squeeze(OptSmoothNoise.SG.x(GridPos.Mid,:))*10^-6;

        OptFreeSGS=squeeze(OptSmoothNoise.SG.s(GridPos.Free,:))*10^-6;
        OptImpSGS=squeeze(OptSmoothNoise.SG.s(GridPos.Imp,:))*10^-6;
        OptMidSGS=squeeze(OptSmoothNoise.SG.s(GridPos.Mid,:))*10^-6;

        %Grid Maximum Smoothing
        MaxFreeSGX=squeeze(MaxSmoothNoise.SG.x(GridPos.Free,:))*10^-6;
        MaxImpSGX=squeeze(MaxSmoothNoise.SG.x(GridPos.Imp,:))*10^-6;
        MaxMidSGX=squeeze(MaxSmoothNoise.SG.x(GridPos.Mid,:))*10^-6;

        MaxFreeSGS=squeeze(MaxSmoothNoise.SG.s(GridPos.Free,:))*10^-6;
        MaxImpSGS=squeeze(MaxSmoothNoise.SG.s(GridPos.Imp,:))*10^-6;
        MaxMidSGS=squeeze(MaxSmoothNoise.SG.s(GridPos.Mid,:))*10^-6;


%Set figure size to full screen for diagnostic purposes
figure('units','normalized','outerposition',[0 0 1 1])
title('Systematic Error No Noise')
% X SG impact edge
subplot(2,3,1)
%FE
plot(FETime,FEimpSGX,'--k')
hold on
%No Smooth
plot(GridTime,NoiseImpSGX,'k')
%OptSmooth
plot(GridTime,OptImpSGX,'b')
%Max Smooth
plot(GridTime,MaxImpSGX,'r')

hold off
title('4P from Impact')
xlabel('Time ($\mu$s)','Interpreter','latex')
ylabel('$\sigma_{xx}$ (MPa)','Interpreter','latex')
legend('FE','No Smoothing','Optimum Smoothing','Max Smoothing')
% X SG Middle 
subplot(2,3,2)
%FE
plot(FETime,FEmidSGX,'--k')
hold on
%No Smooth
plot(GridTime,NoiseMidSGX,'k')
%OptSmooth
plot(GridTime,OptMidSGX,'b')
%Max Smooth
plot(GridTime,MaxMidSGX,'r')

hold off

title('Specimen Middle')
xlabel('Time ($\mu$s)','Interpreter','latex')
ylabel('$\sigma_{xx}$ (MPa)','Interpreter','latex')
% X SG Free edge
subplot(2,3,3)
%FE
plot(FETime,FEFreeSGX,'--k')
hold on
%No Smooth
plot(GridTime,NoiseFreeSGX,'k')
%OptSmooth
plot(GridTime,OptFreeSGX,'b')
%Max Smooth
plot(GridTime,MaxFreeSGX,'r')

hold off

title('4P from Free Edge')
xlabel('Time ($\mu$s)','Interpreter','latex')
ylabel('$\sigma_{xx}$ (MPa)','Interpreter','latex')


% Shear SG Impact
subplot(2,3,4)
%FE
plot(FETime,FEimpSGS,'--k')
hold on
%No Smooth
plot(GridTime,NoiseImpSGS,'k')
%OptSmooth
plot(GridTime,OptImpSGS,'b')
%Max Smooth
plot(GridTime,MaxImpSGS,'r')

hold off

title('4P from Impact')
xlabel('Time ($\mu$s)','Interpreter','latex')
ylabel('$\sigma_{xy}$ (MPa)','Interpreter','latex')

% Shear SG Middle
subplot(2,3,5)
%FE
plot(FETime,FEmidSGS,'--k')
hold on
%No Smooth
plot(GridTime,NoiseMidSGS,'k')
%OptSmooth
plot(GridTime,OptMidSGS,'b')
%Max Smooth
plot(GridTime,MaxMidSGS,'r')

hold off
title('Specimen Middle')
xlabel('Time ($\mu$s)','Interpreter','latex')
ylabel('$\sigma_{xy}$ (MPa)','Interpreter','latex')

% Shear SG free edge
subplot(2,3,6)
%FE
plot(FETime,FEFreeSGS,'--k')
hold on
%No Smooth
plot(GridTime,NoiseFreeSGS,'k')
%OptSmooth
plot(GridTime,OptFreeSGS,'b')
%Max Smooth
plot(GridTime,MaxFreeSGS,'r')


hold off
title('Free Edge')
xlabel('Time ($\mu$s)','Interpreter','latex')
ylabel('$\sigma_{xy}$ (MPa)','Interpreter','latex')

sgtitle('Systematic Error 1 Noise Copy')

STSaveName=strcat(MainSaveDir,'/',ParentDesig, ...
    '_SystematicError1NoiseCopy_SGTime');
saveas(gcf,strcat(STSaveName,'.fig'))
saveas(gcf,strcat(STSaveName,'.svg'))

%% Calculate Smoothed Fields For Spatial Sweep
quest='Source of spatial field sweep?';
SpaSource=questdlg(quest,'Spatial Sweep Source','Files', ...
    'Calculate','Files');
switch SpaSource
    case 'Files'
    %% Choose File without noise
    [SpaSweep.Name,SpaSweep.Path]=uigetfile('*.mat', ...
        'Choose File Containing Spatial Sweep Strains without Noise');
    SpaSweep.File=strcat(SpaSweep.Path,'/',SpaSweep.Name);
    %% Choose file with noise
    [SpaSweepNoise.Name,SpaSweepNoise.Path]=uigetfile('*.mat', ...
        'Choose Noisy Spatial Sweep Data');
    SpaSweepNoise.File=strcat(SpaSweepNoise.Path,'/',SpaSweepNoise.Name);
    %% Load files
    fprintf('Loading Spatial Sweep Data Without Noise \n')
    SpaSweep=load(SpaSweep.File);
    fprintf('Loading Noisy Spatial Sweep Data \n')
    SpaSweepNoise=load(SpaSweepNoise.File);
    fprintf('Spatial Sweep Data Loaded \n')
    
    case 'Calculate'
    %% Calculate Without Noise
    fprintf('Calculating Spatial Sweep Data Without Noise \n')
    NoiseCall='NoNoise';

    SpaSweep=func_SpatialSmoothingSweepStrain(SpaKernVec,RawFields.disp,...
        RawFields.time,RawFields.strain,RawFields.pos,RawFields.grid, ...
        smoothingOpts,RawFields.extrapOpts, ...
        RawFields.globalOpts, ...
        MainSaveDir,ParentDesig,NoiseCall);

    %% Calculate with Noise
    fprintf('Calculateing Spatial Sweep Data with one Noise Copy \n')
    NoiseCall='1NoiseCopy';
    SpaSweepNoise=func_SpatialSmoothingSweepStrain(SpaKernVec,NoiseFields.disp,...
        NoiseFields.time,NoiseFields.strain,NoiseFields.pos,NoiseFields.grid, ...
        smoothingOpts,NoiseFields.extrapOpts, ...
        NoiseFields.globalOpts, ...
        MainSaveDir,ParentDesig,NoiseCall);

    fprintf('spatial smoothing Sweep complete \n')

   

end
%% Calculate average strains
SpaSweep.strain.xAvg=squeeze(mean(SpaSweep.strain.x,1));
SpaSweep.strain.sAvg=squeeze(mean(SpaSweep.strain.s,1));
SpaSweepNoise.strain.xAvg=squeeze(mean(SpaSweepNoise.strain.x,1));
SpaSweepNoise.strain.sAvg=squeeze(mean(SmoothNoise.strain.s,1));
%% Plot Spatial smoothing sweep results at specimen center without noise


%Generate Figure Window
 figure('units','normalized','outerposition',[0 0 1 1])

%Set Persitent variables
GridTime=RawFields.time.vec*10^6;
FETime=FEfields.time.vec*10^6;
    % Set Position Indexes for Grid image data
    GridPos.NumPoints=size(RawFields.disp.x,2);
    GridPos.Free=20; %4 pitches from free edge
    GridPos.Mid=round(GridPos.NumPoints/2);
    GridPos.Imp=GridPos.NumPoints-20; %4 pitches from impact edge
    % Get cooridinates of position indexes
    GridPos.Xfree=RawFields.pos.x(GridPos.Free);
    GridPos.Xmid=RawFields.pos.x(GridPos.Mid);
    GridPos.Ximp=RawFields.pos.x(GridPos.Imp);

    %Set Position indexes and coordintes for finite element data
    FEPos.NumPoints=size(FEfields.disp.x,2);
    FEPos.GridRat=FEPos.NumPoints/GridPos.NumPoints;
    FEPos.Free=round(GridPos.Free*FEPos.GridRat); %4 pitches from free edge
    FEPos.Mid=round(FEPos.NumPoints/2);
    FEPos.Imp=round(GridPos.Imp*FEPos.GridRat); %4 pitches from impact edge
    
    FEPos.Xfree=FEfields.pos.x(FEPos.Free);
    FEPos.Xmid=FEfields.pos.x(FEPos.Mid);
    FEPos.Ximp=FEfields.pos.x(FEPos.Imp);

    FEfields.avgX_strain=squeeze(FEfields.avgX_strain);
    FEfields.avgXY_strain=squeeze(FEfields.avgXY_strain);
    FEmidStrainX=squeeze(FEfields.avgX_strain(FEPos.Mid,:));
    FEmidStrainS=squeeze(FEfields.avgXY_strain(FEPos.Mid,:));

    subplot(1,2,1)
plot(FETime,FEmidStrainX,'--k')
hold on
for k=1:length(SpaKernVec)
    TempStrainX=squeeze(SpaSweep.strain.xAvg(GridPos.Mid,:,k));
    plot(GridTime,TempStrainX)
   
end
hold off
xlabel('time (\mu s)','FontSize',14)
ylabel('$\varepsilon_{XX}$','Interpreter','latex','FontSize',14)
title('Normal Strain','FontSize',14)
legend(['FE';SpaSmoothChar],'FontSize',12)

 subplot(1,2,2)
plot(FETime,FEmidStrainS,'--k')

hold on
for k=1:length(SpaKernVec)
    TempStrainS=squeeze(SpaSweep.strain.sAvg(GridPos.Mid,:,k));
    plot(GridTime,TempStrainS)
   
end
hold off
xlabel('time (\mu s)','FontSize',14)
ylabel('$\varepsilon_{XY}$','Interpreter','latex','FontSize',14)
title('Shear Strain','FontSize',14)
legend(['FE';SpaSmoothChar],'FontSize',12)

sgtitle('Specimen Middle No Noise Spatial Smoothing Sweep','FontSize',18)

SaveName=strcat(MainSaveDir,'/',ParentDesig,...
    '_Spatial Smoothng Sweep_NoNoise_Middle');
saveas(gcf,strcat(SaveName,'.fig'));
saveas(gcf,strcat(SaveName,'.svg'));


%% Plot Spatial smoothing sweep results at specimen center with one noise


%Generate Figure Window
 figure('units','normalized','outerposition',[0 0 1 1])

%Set Persitent variables
GridTime=RawFields.time.vec*10^6;
FETime=FEfields.time.vec*10^6;
    % Set Position Indexes for Grid image data
    GridPos.NumPoints=size(RawFields.disp.x,2);
    GridPos.Free=20; %4 pitches from free edge
    GridPos.Mid=round(GridPos.NumPoints/2);
    GridPos.Imp=GridPos.NumPoints-20; %4 pitches from impact edge
    % Get cooridinates of position indexes
    GridPos.Xfree=RawFields.pos.x(GridPos.Free);
    GridPos.Xmid=RawFields.pos.x(GridPos.Mid);
    GridPos.Ximp=RawFields.pos.x(GridPos.Imp);

    %Set Position indexes and coordintes for finite element data
    FEPos.NumPoints=size(FEfields.disp.x,2);
    FEPos.GridRat=FEPos.NumPoints/GridPos.NumPoints;
    FEPos.Free=round(GridPos.Free*FEPos.GridRat); %4 pitches from free edge
    FEPos.Mid=round(FEPos.NumPoints/2);
    FEPos.Imp=round(GridPos.Imp*FEPos.GridRat); %4 pitches from impact edge
    
    FEPos.Xfree=FEfields.pos.x(FEPos.Free);
    FEPos.Xmid=FEfields.pos.x(FEPos.Mid);
    FEPos.Ximp=FEfields.pos.x(FEPos.Imp);

    FEfields.avgX_strain=squeeze(FEfields.avgX_strain);
    FEfields.avgXY_strain=squeeze(FEfields.avgXY_strain);
    FEmidStrainX=squeeze(FEfields.avgX_strain(FEPos.Mid,:));
    FEmidStrainS=squeeze(FEfields.avgXY_strain(FEPos.Mid,:));

    subplot(1,2,1)
plot(FETime,FEmidStrainX,'--k')
hold on
for k=1:length(SpaKernVec)
    TempStrainX=squeeze(SmoothNoise.strain.xAvg(GridPos.Mid,:,k));
    plot(GridTime,TempStrainX)
   
end
hold off
xlabel('time (\mu s)','FontSize',14)
ylabel('$\varepsilon_{XX}$','Interpreter','latex','FontSize',14)
title('Normal Strain','FontSize',14)
legend(['FE';SpaSmoothChar],'FontSize',12)

 subplot(1,2,2)
plot(FETime,FEmidStrainS,'--k')
TempStrainS=squeeze(SmoothNoise.strain.x(GridPos.Mid,:,k));
hold on
for k=1:length(SpaKernVec)
    TempStrainS=squeeze(SmoothNoise.strain.sAvg(GridPos.Mid,:,k));
    plot(GridTime,TempStrainS)
   
end
hold off
xlabel('time (\mu s)','FontSize',14)
ylabel('$\varepsilon_{XY}$','Interpreter','latex','FontSize',14)
title('Shear Strain','FontSize',14)
legend(['FE';SpaSmoothChar],'FontSize',12)

sgtitle('Specimen Middle 1 Noise copy Spatial Smoothing Sweep','FontSize',18)

SaveName=strcat(MainSaveDir,'/',ParentDesig,...
    '_Spatial Smoothng Sweep_1NoiseCopy_Middle');
saveas(gcf,strcat(SaveName,'.fig'));
saveas(gcf,strcat(SaveName,'.svg'));


%% Substitute FE displacements on specimen edges
%Choose source of interpolated data
quest='Source of the Substiuted FE data';
FESub.Source=questdlg(quest,'Sub Source','File','Calculate','File');

switch FESub.Source
    case 'File'
       [FESub.Name,FESub.Path]=uigetfile('*.mat',...
           'Choose File containing interpolated FE displacements');
       FESub.File=strcat(FESub.Path,'/',FESub.Name);
       FESub=load(FESub.File);

        [FESubOpt.Name,FESubOpt.Path]=uigetfile('*.mat',...
           'Choose Optimal Smoothed interpolated FE displacements');
       FESubOpt.File=strcat(FESubOpt.Path,'/',FESubOpt.Name);
       FESubOpt=load(FESubOpt.File);

    case 'Calculate'
        %% Load Image Def File
        [ImDef.Name,ImDef.Path]=uigetfile('*.mat',...
            'Choose Grid Image Deformation File');
        ImDef.File=strcat(ImDef.Path,'/',ImDef.Name);
        %% Interpolate Finite element data onto grid coordinates
        fprintf('Interpolating Finite Element Data to Grid Coordinates \n')
        
        subset=1;
        FEinterp=func_InterpFEtoGM(ImDef.File,MainSaveDir,ParentDesig,...
            subset);
        
        %% Choose substitution method parameters
        quest='Edge Substitution Method';
        subOpts.method=questdlg(quest,'sub Method','GridPeriod', ...
            'SpatialKernal','GridPeriod');
        switch subOpts.method
            case 'GridPeriod'
                subOpts.PitchNum=str2double(...
                    cell2mat(inputdlg('Number of Periods to subtitute')));
        end

        %% Substitute Interpolated finite element data along the edges and 
            %Calculate strain and acceleration fields
            fprintf('Substituting FE data on grid edges \n')
            SubDesig=strcat(ParentDesig,'NoSmooth_NoNoise');
            RawFields.globalOpts=globalOpts;
            smoothingOpts.spatialSmooth=false;
            smoothingOpts.spatialKernal=[0,0];
            smoothingOpts.FFTempSmooth=false;
            smoothingOpts.FFTemporalKernal=[0,3];
            RawFields.smoothingOpts=smoothingOpts;
            RawFields.globalOpts=globalOpts;
            RawFields.extrapOpts=extrapOpts;
            RawFields.diffOpts=diffOpts;
            
            
        FESub=func_subGMEdgeFields(RawFields.time,FEinterp.disp, ...
            RawFields.disp,RawFields.strain,RawFields.accel, ...
            RawFields.pos,RawFields.grid,... 
            MainSaveDir,SubDesig, ...
            subOpts,RawFields.smoothingOpts,RawFields.globalOpts, ...
            RawFields.extrapOpts,RawFields.diffOpts);
       
        %% Substitute with smooothing
        fprintf('Optimally smoothing substituted data \n')
        SubDesig='OptSmooth_NoNoise';

        FESubOpt=func_subGMEdgeFields(RawFields.time,FEinterp.disp, ...
            RawFields.disp,RawFields.strain,RawFields.accel, ...
            RawFields.pos,RawFields.grid,... 
            MainSaveDir,SubDesig, ...
            subOpts,OptSmooth.smoothingOpts,OptSmooth.globalOpts, ...
            OptSmooth.extrapOpts,OptSmooth.diffOpts);

        fprintf('Substituting with maximum smoothing')
        SubDesig='MaxSmooth_NoNoise';
        FESubMax=func_subGMEdgeFields(RawFields.time,FEinterp.disp, ...
            RawFields.disp,RawFields.strain,RawFields.accel, ...
            RawFields.pos,RawFields.grid,... 
            MainSaveDir,SubDesig, ...
            subOpts,MaxSmooth.smoothingOpts,MaxSmooth.globalOpts, ...
            MaxSmooth.extrapOpts,MaxSmooth.diffOpts);
end
%% Calculate average strains
FESub.strain.xAvg=squeeze(mean(FESub.strain.x,1));
FESub.strain.sAvg=squeeze(mean(FESub.strain.s,1));
FESubOpt.strain.xAvg=squeeze(mean(FESubOpt.strain.x,1));
FESubOpt.strain.sAvg=squeeze(mean(FESubOpt.strain.s,1));
FESubMax.strain.xAvg=squeeze(mean(FESubMax.strain.x,1));
FESubMax.strain.sAvg=squeeze(mean(FESubMax.strain.s,1));
%% Plot the substituted edge data
        FEfields.avgX_strain=squeeze(FEfields.avgX_strain);
        FEfields.avgXY_strain=squeeze(FEfields.avgXY_strain);
        FEFreeStrainX=squeeze(FEfields.avgX_strain(FEPos.Free,:));
        FEimpStrainX=squeeze(FEfields.avgX_strain(FEPos.Imp,:));
        FEmidStrainX=squeeze(FEfields.avgX_strain(FEPos.Mid,:));
    
        FEFreeStrainS=squeeze(FEfields.avgXY_strain(FEPos.Free,:));
        FEimpStrainS=squeeze(FEfields.avgXY_strain(FEPos.Imp,:));
        FEmidStrainS=squeeze(FEfields.avgXY_strain(FEPos.Mid,:));

        %Grid No Smoothing
        GridTime=RawFields.time.vec*10^6;
        SubFreeStrainX=squeeze(FESub.strain.xAvg(GridPos.Free,:));
        SubImpStrainX=squeeze(FESub.strain.xAvg(GridPos.Imp,:));
        SubMidStrainX=squeeze(FESub.strain.xAvg(GridPos.Mid,:));

        SubFreeStrainS=squeeze(FESub.strain.sAvg(GridPos.Free,:));
        SubImpStrainS=squeeze(FESub.strain.sAvg(GridPos.Imp,:));
        SubMidStrainS=squeeze(FESub.strain.sAvg(GridPos.Mid,:));

        %Grid Optimal Smoothing
        OptFreeStrainX=squeeze(FESubOpt.strain.xAvg(GridPos.Free,:));
        OptImpStrainX=squeeze(FESubOpt.strain.xAvg(GridPos.Imp,:));
        OptMidStrainX=squeeze(FESubOpt.strain.xAvg(GridPos.Mid,:));

        OptFreeStrainS=squeeze(FESubOpt.strain.sAvg(GridPos.Free,:));
        OptImpStrainS=squeeze(FESubOpt.strain.sAvg(GridPos.Imp,:));
        OptMidStrainS=squeeze(FESubOpt.strain.sAvg(GridPos.Mid,:));

        %Grid Maximum Smoothing
        MaxFreeStrainX=squeeze(FESubMax.strain.xAvg(GridPos.Free,:));
        MaxImpStrainX=squeeze(FESubMax.strain.xAvg(GridPos.Imp,:));
        MaxMidStrainX=squeeze(FESubMax.strain.xAvg(GridPos.Mid,:));

        MaxFreeStrainS=squeeze(FESubMax.strain.sAvg(GridPos.Free,:));
        MaxImpStrainS=squeeze(FESubMax.strain.sAvg(GridPos.Imp,:));
        MaxMidStrainS=squeeze(FESubMax.strain.sAvg(GridPos.Mid,:));


%Set figure size to full screen for diagnostic purposes
figure('units','normalized','outerposition',[0 0 1 1])
title('Systematic Error No Noise','FontSize',14)
% X strain impact edge
subplot(2,3,1)
%FE
plot(FETime,FEimpStrainX,'--k')
hold on
%No Smooth
plot(GridTime,SubImpStrainX,'k')
%OptSmooth
plot(GridTime,OptImpStrainX,'b')
%Max Smooth
plot(GridTime,MaxImpStrainX,'r')

hold off
title('4P from Impact','FontSize',14)
xlabel('Time ($\mu$s)','Interpreter','latex','FontSize',14)
ylabel('$\varepsilon_{xx}$','Interpreter','latex','FontSize',14)
legend('FE','No Smoothing','Optimum Smoothing','Max Smoothing','FontSize',14)
% X strain Middle 
subplot(2,3,2)
%FE
plot(FETime,FEmidStrainX,'--k')
hold on
%No Smooth
plot(GridTime,SubMidStrainX,'k')
%OptSmooth
plot(GridTime,OptMidStrainX,'b')
%Max Smooth
plot(GridTime,MaxMidStrainX,'r')

hold off

title('Specimen Middle','FontSize',14)
xlabel('Time ($\mu$s)','Interpreter','latex','FontSize',14)
ylabel('$\varepsilon_{xx}$','Interpreter','latex','FontSize',14)
% X Strain Free edge
subplot(2,3,3)
%FE
plot(FETime,FEFreeStrainX,'--k')
hold on
%No Smooth
plot(GridTime,SubFreeStrainX,'k')
%OptSmooth
plot(GridTime,OptFreeStrainX,'b')
%Max Smooth
plot(GridTime,MaxFreeStrainX,'r')

hold off

title('4P from Free Edge','FontSize',14)
xlabel('Time ($\mu$s)','Interpreter','latex','FontSize',14)
ylabel('$\varepsilon_{xx}$','Interpreter','latex','FontSize',14)


% Shear Strain Impact
subplot(2,3,4)
%FE
plot(FETime,FEimpStrainS,'--k')
hold on
%No Smooth
plot(GridTime,SubImpStrainS,'k')
%OptSmooth
plot(GridTime,OptImpStrainS,'b')
%Max Smooth
plot(GridTime,MaxImpStrainS,'r')

hold off

title('4P from Impact','FontSize',14)
xlabel('Time ($\mu$s)','Interpreter','latex','FontSize',14)
ylabel('$\varepsilon_{xy}$','Interpreter','latex','FontSize',14)

% Shear Strain Middle
subplot(2,3,5)
%FE
plot(FETime,FEmidStrainS,'--k')
hold on
%No Smooth
plot(GridTime,SubMidStrainS,'k')
%OptSmooth
plot(GridTime,OptMidStrainS,'b')
%Max Smooth
plot(GridTime,MaxMidStrainS,'r')

hold off
title('Specimen Middle','FontSize',14)
xlabel('Time ($\mu$s)','Interpreter','latex','FontSize',14)
ylabel('$\varepsilon_{xy}$','Interpreter','latex','FontSize',14)

% Shear strain free edge
subplot(2,3,6)
%FE
plot(FETime,FEFreeStrainS,'--k')
hold on
%No Smooth
plot(GridTime,SubFreeStrainS,'k')
%OptSmooth
plot(GridTime,OptFreeStrainS,'b')
%Max Smooth
plot(GridTime,MaxFreeStrainS,'r')


hold off
title('Free Edge','FontSize',14)
xlabel('Time ($\mu$s)','Interpreter','latex','FontSize',14)
ylabel('$\varepsilon_{xy}$','Interpreter','latex','FontSize',14)

sgtitle('Substituted Edge Fields','FontSize',14)

STSaveName=strcat(MainSaveDir,'/',ParentDesig, ...
    '_FESubNoNoise_StrainTime');
saveas(gcf,strcat(STSaveName,'.fig'))
saveas(gcf,strcat(STSaveName,'.svg'))

%% Investigate Post correcting shear strain
RawFields.DispCorr.strainMethod='GridPeriod';
RawFields.DispCorr.strainPitchNum=1;
            SubDesig=strcat(ParentDesig,'NoSmooth_NoNoise');
            RawFields.globalOpts=globalOpts;
            smoothingOpts.spatialSmooth=false;
            smoothingOpts.spatialKernal=[0,0];
            smoothingOpts.FFTempSmooth=false;
            smoothingOpts.FFTemporalKernal=[0,3];
            RawFields.smoothingOpts=smoothingOpts;
            RawFields.globalOpts=globalOpts;
            RawFields.extrapOpts=extrapOpts;
            RawFields.diffOpts=diffOpts;
            
fprintf('Correcting Shear strains in RawFields \n')
TempShear=func_PostCorrectShear(RawFields.grid,RawFields.strain, ...
    RawFields.pos,RawFields.time,RawFields.DispCorr, ...
    RawFields.smoothingOpts);
RawFields.strain.s_corr=TempShear.s;
clear TempShear

fprintf('Correcting Max Smoothed \n')
MaxSmooth.smoothingOpts=MaxSmooth.smoothingOPpts;
TempShear=func_PostCorrectShear(MaxSmooth.grid,MaxSmooth.strain, ...
    MaxSmooth.pos,MaxSmooth.time,RawFields.DispCorr, ...
    MaxSmooth.smoothingOpts);
MaxSmooth.strain.s_corr=TempShear.s;
clear TempShear

fprintf('Correcting Optimal smoothed \n')
OptSmooth.smoothingOpts=OptSmooth.smoothingOPpts;
TempShear=func_PostCorrectShear(OptSmooth.grid,OptSmooth.strain, ...
    OptSmooth.pos,OptSmooth.time,RawFields.DispCorr, ...
    OptSmooth.smoothingOpts);
OptSmooth.strain.s_corr=TempShear.s;

fprintf('strain correction complete \n')

%% Calculate average corrected strains
RawFields.strain.s_corrAvg=squeeze(mean(RawFields.strain.s_corr));
MaxSmooth.strain.s_corrAvg=squeeze(mean(MaxSmooth.strain.s_corr));
OptSmooth.strain.s_corrAvg=squeeze(mean(OptSmooth.strain.s_corr));

fprintf('Generating Strain-Time Diagnostic Plots With No Noise \n')
    % Set Position Indexes for Grid image data
    GridPos.NumPoints=size(RawFields.disp.x,2);
    GridPos.Free=20; %4 pitches from free edge
    GridPos.Mid=round(GridPos.NumPoints/2);
    GridPos.Imp=GridPos.NumPoints-20; %4 pitches from impact edge
    % Get cooridinates of position indexes
    GridPos.Xfree=RawFields.pos.x(GridPos.Free);
    GridPos.Xmid=RawFields.pos.x(GridPos.Mid);
    GridPos.Ximp=RawFields.pos.x(GridPos.Imp);

    %Set Position indexes and coordintes for finite element data
    FEPos.NumPoints=size(FEfields.disp.x,2);
    FEPos.GridRat=FEPos.NumPoints/GridPos.NumPoints;
    FEPos.Free=round(GridPos.Free*FEPos.GridRat); %4 pitches from free edge
    FEPos.Mid=round(FEPos.NumPoints/2);
    FEPos.Imp=round(GridPos.Imp*FEPos.GridRat); %4 pitches from impact edge
    
    FEPos.Xfree=FEfields.pos.x(FEPos.Free);
    FEPos.Xmid=FEfields.pos.x(FEPos.Mid);
    FEPos.Ximp=FEfields.pos.x(FEPos.Imp);
    
  %Print Plot selections for evaluation purposes
  FEPos.freeString=num2str(FEPos.Xfree);
  GridPos.freeString=num2str(GridPos.Xfree);
  fprintf(strcat('Free Coordinate FE:',FEPos.freeString,' Grid: ',...
      GridPos.freeString,'\n'));
  
  FEPos.midString=num2str(FEPos.Xmid);
  GridPos.midString=num2str(GridPos.Xmid);
  fprintf(strcat('Mid Coordinate FE:',FEPos.midString,' Grid: ',...
      GridPos.midString,'\n'));

  FEPos.impString=num2str(FEPos.Ximp);
  GridPos.impString=num2str(GridPos.Ximp);
  fprintf(strcat('Impact Coordinate FE:',FEPos.impString,' Grid: ',...
      GridPos.impString,'\n'));
  
    %Generate plot variables for easier coding
        %Finite Element
        FETime=FEfields.time.vec*10^6;
        FEfields.avgX_strain=squeeze(FEfields.avgX_strain);
        FEfields.avgXY_strain=squeeze(FEfields.avgXY_strain);
        FEFreeStrainX=squeeze(FEfields.avgX_strain(FEPos.Free,:));
        FEimpStrainX=squeeze(FEfields.avgX_strain(FEPos.Imp,:));
        FEmidStrainX=squeeze(FEfields.avgX_strain(FEPos.Mid,:));
    
        FEFreeStrainS=squeeze(FEfields.avgXY_strain(FEPos.Free,:));
        FEimpStrainS=squeeze(FEfields.avgXY_strain(FEPos.Imp,:));
        FEmidStrainS=squeeze(FEfields.avgXY_strain(FEPos.Mid,:));

        %Grid No Smoothing
        GridTime=RawFields.time.vec*10^6;
        RawFreeStrainX=squeeze(RawFields.avgX_strain(GridPos.Free,:));
        RawImpStrainX=squeeze(RawFields.avgX_strain(GridPos.Imp,:));
        RawMidStrainX=squeeze(RawFields.avgX_strain(GridPos.Mid,:));

        RawFreeStrainS=squeeze(RawFields.strain.s_corrAvg(GridPos.Free,:));
        RawImpStrainS=squeeze(RawFields.strain.s_corrAvg(GridPos.Imp,:));
        RawMidStrainS=squeeze(RawFields.strain.s_corrAvg(GridPos.Mid,:));

        %Grid Optimal Smoothing
        OptFreeStrainX=squeeze(OptSmooth.avgX_strain(GridPos.Free,:));
        OptImpStrainX=squeeze(OptSmooth.avgX_strain(GridPos.Imp,:));
        OptMidStrainX=squeeze(OptSmooth.avgX_strain(GridPos.Mid,:));

        OptFreeStrainS=squeeze(OptSmooth.strain.s_corrAvg(GridPos.Free,:));
        OptImpStrainS=squeeze(OptSmooth.strain.s_corrAvg(GridPos.Imp,:));
        OptMidStrainS=squeeze(OptSmooth.strain.s_corrAvg(GridPos.Mid,:));

        %Grid Maximum Smoothing
        MaxFreeStrainX=squeeze(MaxSmooth.avgX_strain(GridPos.Free,:));
        MaxImpStrainX=squeeze(MaxSmooth.avgX_strain(GridPos.Imp,:));
        MaxMidStrainX=squeeze(MaxSmooth.avgX_strain(GridPos.Mid,:));

        MaxFreeStrainS=squeeze(MaxSmooth.strain.s_corrAvg(GridPos.Free,:));
        MaxImpStrainS=squeeze(MaxSmooth.strain.s_corrAvg(GridPos.Imp,:));
        MaxMidStrainS=squeeze(MaxSmooth.strain.s_corrAvg(GridPos.Mid,:));


%Set figure size to full screen for diagnostic purposes
figure('units','normalized','outerposition',[0 0 1 1])
title('Systematic Error No Noise')
% X strain impact edge
subplot(2,3,1)
%FE
plot(FETime,FEimpStrainX,'--k')
hold on
%No Smooth
plot(GridTime,RawImpStrainX,'k')
%OptSmooth
plot(GridTime,OptImpStrainX,'b')
%Max Smooth
plot(GridTime,MaxImpStrainX,'r')

hold off
title('4P from Impact','FontSize',14)
xlabel('Time ($\mu$s)','Interpreter','latex','FontSize',14)
ylabel('$\varepsilon_{xx}$','Interpreter','latex','FontSize',14)
legend('FE','No Smoothing','Optimum Smoothing','Max Smoothing','location','northwest')
% X strain Middle 
subplot(2,3,2)
%FE
plot(FETime,FEmidStrainX,'--k')
hold on
%No Smooth
plot(GridTime,RawMidStrainX,'k')
%OptSmooth
plot(GridTime,OptMidStrainX,'b')
%Max Smooth
plot(GridTime,MaxMidStrainX,'r')

hold off

title('Specimen Middle','FontSize',14)
xlabel('Time ($\mu$s)','Interpreter','latex','FontSize',14)
ylabel('$\varepsilon_{xx}$','Interpreter','latex','FontSize',14)
% X Strain Free edge
subplot(2,3,3)
%FE
plot(FETime,FEFreeStrainX,'--k')
hold on
%No Smooth
plot(GridTime,RawFreeStrainX,'k')
%OptSmooth
plot(GridTime,OptFreeStrainX,'b')
%Max Smooth
plot(GridTime,MaxFreeStrainX,'r')

hold off

title('4P from Free Edge')
xlabel('Time ($\mu$s)','Interpreter','latex','FontSize',14)
ylabel('$\varepsilon_{xx}$','Interpreter','latex','FontSize',14)


% Shear Strain Impact
subplot(2,3,4)
%FE
plot(FETime,FEimpStrainS,'--k')
hold on
%No Smooth
plot(GridTime,RawImpStrainS,'k')
%OptSmooth
plot(GridTime,OptImpStrainS,'b')
%Max Smooth
plot(GridTime,MaxImpStrainS,'r')

hold off

title('4P from Impact')
xlabel('Time ($\mu$s)','Interpreter','latex')
ylabel('$\varepsilon_{xy}$','Interpreter','latex')

% Shear Strain Middle
subplot(2,3,5)
%FE
plot(FETime,FEmidStrainS,'--k')
hold on
%No Smooth
plot(GridTime,RawMidStrainS,'k')
%OptSmooth
plot(GridTime,OptMidStrainS,'b')
%Max Smooth
plot(GridTime,MaxMidStrainS,'r')

hold off
title('Specimen Middle')
xlabel('Time ($\mu$s)','Interpreter','latex','FontSize',14)
ylabel('$\varepsilon_{xy}$','Interpreter','latex','FontSize',14)

% Shear strain free edge
subplot(2,3,6)
%FE
plot(FETime,FEFreeStrainS,'--k')
hold on
%No Smooth
plot(GridTime,RawFreeStrainS,'k')
%OptSmooth
plot(GridTime,OptFreeStrainS,'b')
%Max Smooth
plot(GridTime,MaxFreeStrainS,'r')


hold off
title('Free Edge')
xlabel('Time ($\mu$s)','Interpreter','latex','FontSize',14)
ylabel('$\varepsilon_{xy}$','Interpreter','latex','FontSize',14)

sgtitle('Systematic Error No Noise','FontSize',14)

STSaveName=strcat(MainSaveDir,'/',ParentDesig, ...
    '_ShearCorr_StrainTime');
saveas(gcf,strcat(STSaveName,'.fig'))
saveas(gcf,strcat(STSaveName,'.svg'))



%% Calculate Average Fields over 30 Noise Copies


%%
