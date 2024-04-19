% This script is written to identify the noise floor on synthetic grid
% images

%Author: Andrew Matejunas
%Date: 2022-12-18

%Change Log:
    %2023-01-04: Changed the code to create contour plots of the noise
        %floor and bias in displacement, acceleration, strain, stress gauge
        %stress and constitutive model stresses for a sweep of smoothing
        %parameters.
  

%%
close all; clear variables; clc
%% Choose folder to save data
savePath=uigetdir({},'Choose Folder to Save Results');

%% Load properties file
[refpar.name,refpar.path]=uigetfile('*.mat',...
    'Load File containing constitutive parameters');
refpar.File=strcat(refpar.path,'/',refpar.name);
load(refpar.File);
RefPar=MatProps;
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

%% Create smoothing kernal vecotors

%Vector of spatial smoothing kernal sizes in pixels
SpaKernVec=[0,5:2:21];
SpaKernNum=length(SpaKernVec);
%Vector of temporal smoothing kernal sizes in frames
TemKernVec=[0,5:2:21];
TemKernNum=length(TemKernVec);
TotNum=SpaKernNum*TemKernNum;

%% Create smoothing options that do not change
smoothingOpts.spatialFilt='gauss';
smoothingOpts.spatialEdgeMode='symmetric';
smoothingOpts.FFTemporalFilt='sgolay';
smoothingOpts.FFTemporalPad=1;
smoothingOpts.FFTemporalPadMethod='replicate';
smoothingOpts.WATemporalAvgFirst=0;
smoothingOpts.WATempSmooth=0;
smoothingOpts.WATemporalFilt='sgolay';
smoothingOpts.WATemporalKernal=[11,3];
smoothingOpts.WATemporalPad=0;
smoothingOpts.WATemporalPadFrames=3;
smoothingOpts.WATemporalPadMethod='replicate';

%% initialize output variables
[DispNF.x,DispNF.y,DispBias.x,DispBias.y,...
    StrainNF.x,StrainNF.y,StrainNF.s,...
    StrainBias.x,StrainBias.y,StrainBias.s,...
    AccelNF.x,AccelNF.y,AccelBias.x,AccelBias.y,...
    SGNF.x,SGNF.s,SGBias.x,SGBias.s,...
    StressNF.x,StressNF.y,StressNF.s,...
    StressBias.x,StressBias.y,StressBias.s]=...
    deal(zeros(SpaKernNum,TemKernNum));

%% Generate the progress bar
itNum=0;
progmsg=strcat('0/',num2str(TotNum),' Smoothing Iterations Complete');
Progress=waitbar(0,progmsg);

%% Calculate the Noise Floor/resolution
for SK=1:SpaKernNum
    for TK=1:TemKernNum
itNum=itNum+1;
itCount=num2str(itNum);

        %% Set smoothing kernal sizes
        if SpaKernVec(SK)==0
            smoothingOpts.spatialSmooth=0;
        else
            smoothingOpts.spatialSmooth=1;
            smoothingOpts.spatialKernal=SpaKernVec(SK)*[1,1];
        end

        if TemKernVec(TK)==0
            smoothingOpts.FFTempSmooth=0;
        else
            smoothingOpts.FFTempSmooth=1;
            smoothingOpts.FFTemporalKernal=[TemKernVec(TK),3];
            smoothingOpts.FFTemporalPadFrames=3;
        end

        NoiseMag=0.4;

        [DispRes,StrainRes,AccelRes,SGRes,StressRes,...
            ~,~,~,~,~] =...
            func_computeResolutionNoiseFloorV2(NoiseMag,grid,time, ...
            imagePath,imageFile,gridMethodOpts,specimen, ...
            DispCorr,CondOpts,smoothingOpts,globalOpts, ...
            extrapOpts,diffOpts,material,RefPar);

        %% Store the noise floors and biases in the output variables
        %Displacement noise floor and bias
        DispNF.x(SK,TK)=DispRes.x.NF;
        DispBias.x(SK,TK)=DispRes.x.avg;
        DispNF.y(SK,TK)=DispRes.y.NF;
        DispBias.y(SK,TK)=DispRes.y.avg;

        %Acceleration noise floor and bias
        AccelNF.x(SK,TK)=AccelRes.x.NF;
        AccelBias.x(SK,TK)=AccelRes.x.avg;
        AccelNF.y(SK,TK)=AccelRes.y.NF;
        AccelBias.y(SK,TK)=AccelRes.y.avg;

        %Strain noise floor and bias
        StrainNF.x(SK,TK)=StrainRes.x.NF;
        StrainBias.x(SK,TK)=StrainRes.x.avg;
        StrainNF.y(SK,TK)=StrainRes.y.NF;
        StrainBias.y(SK,TK)=StrainRes.y.avg;
        StrainNF.s(SK,TK)=StrainRes.s.NF;
        StrainBias.s(SK,TK)=StrainRes.s.avg;

        %Stress noise floor and bias
        StressNF.x(SK,TK)=StressRes.xAv.NF;
        StressBias.x(SK,TK)=StressRes.xAv.avg;
        StressNF.y(SK,TK)=StressRes.y.NF;
        StressBias.y(SK,TK)=StressRes.y.avg;
        StressNF.s(SK,TK)=StressRes.s.NF;
        StressBias.s(SK,TK)=StressRes.s.avg;

        %Stress gauge noise floor and bias
        SGNF.x(SK,TK)=SGRes.x.NF;
        SGBias.x(SK,TK)=SGRes.x.avg;
        SGNF.s(SK,TK)=SGRes.s.NF;
        SGBias.s(SK,TK)=SGRes.s.avg;

        %% Generate Progress bar
        progmsg=strcat(itCount,'/',num2str(TotNum),' Iterations Complete');
        Progress=waitbar(itNum/TotNum,Progress,progmsg);

    end
end

%% Set colormap
ColorScheme='hot';
%% Plot the displacement Noise floor and bias as a function 
figure('units','normalized','outerposition',[0 0 1 1])
%x direction noise floor
subplot(2,2,1)
contourf(SpaKernVec,TemKernVec,DispNF.x',50)
title('u_x Noise Floor')
xlabel('Spatial Kernal (px)')
ylabel('Temporal Kernal (frames)')

cx=colorbar;
colormap(ColorScheme);
cx.Label.String='stdev(u_x) (m)';

%x direction noise bias
subplot(2,2,2)
contourf(SpaKernVec,TemKernVec,DispBias.x',50)
title('u_x Noise Bias')
xlabel('Spatial Kernal (px)')
ylabel('Temporal Kernal (frames)')

cx=colorbar;
colormap(ColorScheme);
cx.Label.String='mean(u_x) (m)';

%y direction noise floor
subplot(2,2,3)
contourf(SpaKernVec,TemKernVec,DispNF.y',50)
title('u_y Noise Floor')
xlabel('Spatial Kernal (px)')
ylabel('Temporal Kernal (frames)')

cx=colorbar;
colormap(ColorScheme);
cx.Label.String='stdev(u_y) (m)';

%x direction noise bias
subplot(2,2,4)
contourf(SpaKernVec,TemKernVec,DispBias.y',50)
title('u_y Noise Bias')
xlabel('Spatial Kernal (px)')
ylabel('Temporal Kernal (frames)')

cx=colorbar;
colormap(ColorScheme);
cx.Label.String='mean(u_y) (m)';

FigSaveName=strcat(savePath,'/SyntheticImage_DispFloorBias');
saveas(gcf,FigSaveName,'fig');
saveas(gcf,FigSaveName,'svg');

%% Acceleration Noise Floor and bias
figure('units','normalized','outerposition',[0 0 1 1])
%x direction noise floor
subplot(2,2,1)
contourf(SpaKernVec,TemKernVec,AccelNF.x',50)
title('a_x Noise Floor')
xlabel('Spatial Kernal (px)')
ylabel('Temporal Kernal (frames)')

cx=colorbar;
colormap(ColorScheme);
cx.Label.String='stdev(a_x) (m/s^2)';

%x direction noise bias
subplot(2,2,2)
contourf(SpaKernVec,TemKernVec,AccelBias.x',50)
title('a_x Noise Bias')
xlabel('Spatial Kernal (px)')
ylabel('Temporal Kernal (frames)')

cx=colorbar;
colormap(ColorScheme);
cx.Label.String='mean(a_x) (m/s^2)';

%y direction noise floor
subplot(2,2,3)
contourf(SpaKernVec,TemKernVec,AccelNF.y',50)
title('a_y Noise Floor')
xlabel('Spatial Kernal (px)')
ylabel('Temporal Kernal (frames)')

cx=colorbar;
colormap(ColorScheme);
cx.Label.String='stdev(a_y) (m/s^2)';

%x direction noise bias
subplot(2,2,4)
contourf(SpaKernVec,TemKernVec,AccelBias.y',50)
title('a_y Noise Bias')
xlabel('Spatial Kernal (px)')
ylabel('Temporal Kernal (frames)')

cx=colorbar;
colormap(ColorScheme);
cx.Label.String='mean(a_y) (m/s^2)';

FigSaveName=strcat(savePath,'/SyntheticImage_AccelFloorBias');
saveas(gcf,FigSaveName,'fig');
saveas(gcf,FigSaveName,'svg');

%% Plot Strain Noise Floor and Bias
figure('units','normalized','outerposition',[0 0 1 1])
%x direction noise floor
subplot(3,2,1)
contourf(SpaKernVec,TemKernVec,StrainNF.x',50)
title('$\varepsilon_{xx}$ Noise Floor','Interpreter','latex')
xlabel('Spatial Kernal (px)')
ylabel('Temporal Kernal (frames)')

cx=colorbar;
colormap(ColorScheme);
cx.Label.String='$stdev(\varepsilon_{xx})$';
cx.Label.Interpreter='latex';

%x direction noise bias
subplot(3,2,2)
contourf(SpaKernVec,TemKernVec,StrainBias.x',50)
title('$\varepsilon_{xx}$ Noise Bias','Interpreter','latex')
xlabel('Spatial Kernal (px)')
ylabel('Temporal Kernal (frames)')

cx=colorbar;
colormap(ColorScheme);
cx.Label.String='$mean(\varepsilon_{xx})$';
cx.Label.Interpreter='latex';

%y direction noise floor
subplot(3,2,3)
contourf(SpaKernVec,TemKernVec,StrainNF.y',50)
title('$\varepsilon_{yy}$ Noise Floor','Interpreter','latex')
xlabel('Spatial Kernal (px)')
ylabel('Temporal Kernal (frames)')

cx=colorbar;
colormap(ColorScheme);
cx.Label.String='$stdev(\varepsilon_{yy})$';
cx.Label.Interpreter='latex';

%y direction noise bias
subplot(3,2,4)
contourf(SpaKernVec,TemKernVec,StrainBias.y',50)
title('$\varepsilon_{yy}$ Noise Bias','Interpreter','latex')
xlabel('Spatial Kernal (px)')
ylabel('Temporal Kernal (frames)')

cx=colorbar;
colormap(ColorScheme);
cx.Label.String='$mean(\varepsilon_{yy})$';
cx.Label.Interpreter='latex';

%xy direction noise floor
subplot(3,2,5)
contourf(SpaKernVec,TemKernVec,StrainNF.s',50)
title('$\varepsilon_{xy}$ Noise Floor','Interpreter','latex')
xlabel('Spatial Kernal (px)')
ylabel('Temporal Kernal (frames)')

cx=colorbar;
colormap(ColorScheme);
cx.Label.String='$stdev(\varepsilon_{xy})$';
cx.Label.Interpreter='latex';

%y direction noise bias
subplot(3,2,6)
contourf(SpaKernVec,TemKernVec,StrainBias.s',50)
title('$\varepsilon_{xy}$ Noise Bias','Interpreter','latex')
xlabel('Spatial Kernal (px)')
ylabel('Temporal Kernal (frames)')

cx=colorbar;
colormap(ColorScheme);
cx.Label.String='$mean(\varepsilon_{xy})$';
cx.Label.Interpreter='latex';

FigSaveName=strcat(savePath,'/SyntheticImage_StrainFloorBias');
saveas(gcf,FigSaveName,'fig');
saveas(gcf,FigSaveName,'svg');

%% Plot Constitutive Model Stress Noise Floor and Bias
figure('units','normalized','outerposition',[0 0 1 1])
%x direction noise floor
subplot(3,2,1)
contourf(SpaKernVec,TemKernVec,StressNF.x',50)
title('$\sigma_{xx}$ Noise Floor','Interpreter','latex')
xlabel('Spatial Kernal (px)')
ylabel('Temporal Kernal (frames)')

cx=colorbar;
colormap(ColorScheme);
cx.Label.String='$stdev(\sigma_{xx})$';
cx.Label.Interpreter='latex';

%x direction noise bias
subplot(3,2,2)
contourf(SpaKernVec,TemKernVec,StressBias.x',50)
title('$\sigma_{xx}$ Noise Bias','Interpreter','latex')
xlabel('Spatial Kernal (px)')
ylabel('Temporal Kernal (frames)')

cx=colorbar;
colormap(ColorScheme);
cx.Label.String='$mean(\sigma_{xx})$';
cx.Label.Interpreter='latex';

%y direction noise floor
subplot(3,2,3)
contourf(SpaKernVec,TemKernVec,StressNF.y',50)
title('$\sigma_{yy}$ Noise Floor','Interpreter','latex')
xlabel('Spatial Kernal (px)')
ylabel('Temporal Kernal (frames)')

cx=colorbar;
colormap(ColorScheme);
cx.Label.String='$stdev(\sigma_{yy})$';
cx.Label.Interpreter='latex';

%y direction noise bias
subplot(3,2,4)
contourf(SpaKernVec,TemKernVec,StressBias.y',50)
title('$\sigma_{yy}$ Noise Bias','Interpreter','latex')
xlabel('Spatial Kernal (px)')
ylabel('Temporal Kernal (frames)')

cx=colorbar;
colormap(ColorScheme);
cx.Label.String='$mean(\sigma_{yy})$';
cx.Label.Interpreter='latex';

%xy direction noise floor
subplot(3,2,5)
contourf(SpaKernVec,TemKernVec,StressNF.s',50)
title('$\sigma_{xy}$ Noise Floor','Interpreter','latex')
xlabel('Spatial Kernal (px)')
ylabel('Temporal Kernal (frames)')

cx=colorbar;
colormap(ColorScheme);
cx.Label.String='$stdev(\sigma_{xy})$';
cx.Label.Interpreter='latex';

%y direction noise bias
subplot(3,2,6)
contourf(SpaKernVec,TemKernVec,StressBias.s',50)
title('$\sigma_{xy}$ Noise Bias','Interpreter','latex')
xlabel('Spatial Kernal (px)')
ylabel('Temporal Kernal (frames)')

cx=colorbar;
colormap(ColorScheme);
cx.Label.String='$mean(\sigma_{xy})$';
cx.Label.Interpreter='latex';

FigSaveName=strcat(savePath,'/SyntheticImage_StressFloorBias');
saveas(gcf,FigSaveName,'fig');
saveas(gcf,FigSaveName,'svg');

%% Plot Stress Guage Noise Floor and Bias
figure('units','normalized','outerposition',[0 0 1 1])
%x direction noise floor
subplot(2,2,1)
contourf(SpaKernVec,TemKernVec,SGNF.x',50)
title('$\sigma_{xx_{sg}}$ Noise Floor','Interpreter','latex')
xlabel('Spatial Kernal (px)')
ylabel('Temporal Kernal (frames)')

cx=colorbar;
colormap(ColorScheme);
cx.Label.String='$stdev(\sigma_{xx})$';
cx.Label.Interpreter='latex';

%x direction noise bias
subplot(2,2,2)
contourf(SpaKernVec,TemKernVec,SGBias.x',50)
title('$\sigma_{xx_{sg}}$ Noise Bias','Interpreter','latex')
xlabel('Spatial Kernal (px)')
ylabel('Temporal Kernal (frames)')

cx=colorbar;
colormap(ColorScheme);
cx.Label.String='$mean(\sigma_{xx})$';
cx.Label.Interpreter='latex';

%xy direction noise floor
subplot(2,2,3)
contourf(SpaKernVec,TemKernVec,SGNF.s',50)
title('$\sigma_{xy_{sg}}$ Noise Floor','Interpreter','latex')
xlabel('Spatial Kernal (px)')
ylabel('Temporal Kernal (frames)')

cx=colorbar;
colormap(ColorScheme);
cx.Label.String='$stdev(\sigma_{xy})$';
cx.Label.Interpreter='latex';

%y direction noise bias
subplot(2,2,4)
contourf(SpaKernVec,TemKernVec,SGBias.s',50)
title('$\sigma_{xy_{sg}}$ Noise Bias','Interpreter','latex')
xlabel('Spatial Kernal (px)')
ylabel('Temporal Kernal (frames)')

cx=colorbar;
colormap(ColorScheme);
cx.Label.String='$mean(\sigma_{xy})$';
cx.Label.Interpreter='latex';

FigSaveName=strcat(savePath,'/SyntheticImage_SGFloorBias');
saveas(gcf,FigSaveName,'fig');
saveas(gcf,FigSaveName,'svg');

%% Generate noise floor comparison plots for stress guage and constitutive 
    %Model stresses

ConstLineX=squeeze(StressNF.x(:,2));
ConstLineY=squeeze(StressNF.y(:,2));
ConstLineS=squeeze(StressNF.s(:,2));

SGLineX=squeeze(SGNF.x(2,:));
SGLineS=squeeze(SGNF.s(2,:));

StressLim=[min(ConstLineX),max(ConstLineX)];

figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1)
plot(SpaKernVec,ConstLineX)
hold on
plot(SpaKernVec,ConstLineS,'--')
plot(SpaKernVec,ConstLineY,'-.')
hold off
xlabel('Spatial Kernal (px)')
ylabel(' Constitutive Stress Noise Floor (Pa)')
legend('\sigma_{xx}','\sigma_{xy}','\sigma_{yy}','location','northeast')


subplot(1,2,2)
plot(SpaKernVec,SGLineX)
hold on
plot(SpaKernVec,SGLineS,'--')
hold off
xlabel('Temporal Kernal (frames)')
ylabel('Stress Gauge Noise Floor (Pa)')
legend('\sigma_{xx}','\sigma_{xy}','location','northeast')

FigSaveName=strcat(savePath,'/SyntheticImage_SGConstFloorComp');
saveas(gcf,FigSaveName,'fig');
saveas(gcf,FigSaveName,'svg');
SaveName=strcat(savePath,'/SyntheticImageNoiseFloor_Mag04.mat');
save(SaveName);
