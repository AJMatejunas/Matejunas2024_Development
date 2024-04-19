% This script is written to identify the noise floor on synthetic grid
% images

%Author: Andrew Matejunas
%Date: 2022-12-18

%Change Log:
    %2023-01-04: Changed the code to create contour plots of the noise
        %floor and bias in displacement, acceleration, strain, stress gauge
        %stress and constitutive model stresses for a sweep of smoothing
        %parameters.
    %V2-2023-07-02: Changed code to calculate accelerations and strains using
        %V4 versions
  

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
globalOpts.FEValidMode =false;
globalOpts.calcKinFieldsFromDisp = true;

%% Remove some options from original IBII codes
hardCodePath=0;
saveGridData=0;

%% Define Smoothing Opts
TempKernVec=[0,5:2:27];
%TempKernVec=0;
TempKernNum=length(TempKernVec);  
SpaKernVec=[0,5:2:15];
%SpaKernVec=0;

%record Number of Kernels
SpaKernNum=length(SpaKernVec);
TotNum=TempKernNum*SpaKernNum;

%% set the extrapolation options that do not depend on 
extrapOpts.accel.tempPadOn=false;       


[strainCorr,strainSmooth,dispCorr,accelCorr]=deal(cell(SpaKernNum,1));

for m=1:length(SpaKernVec)
            if SpaKernVec(m)==0
            smoothOpts.strain.spatialSmooth=false;
            % Set default displacement cropping and extrapolation
            %parameters at 7 pixels
            extrapOpts.disp.extrapPx1st=7;
            extrapOpts.disp.extrapPx2nd=7;
            extrapOpts.disp.cropPx1st=7;
            extrapOpts.disp.cropPx2nd=7;
        else
            smoothOpts.strain.spatialSmooth=true;
            smoothOpts.strain.spatialKernelSize=SpaKernVec(m)*[1,1];
            smoothOpts.strain.spatialKernelStd=[10,10];

            % Scale displacement cropping and extrapolation parameters
                %with smoothing Kernel size
                if SpaKernVec(m)<=13
                   extrapOpts.disp.cropPx1st=7;
                   extrapOpts.disp.cropPx2nd=7;
                   extrapOpts.disp.extrapPx1st=extrapOpts.disp.cropPx1st;
                   extrapOpts.disp.extrapPx2nd=extrapOpts.disp.cropPx2nd;
                else

                   extrapOpts.disp.cropPx1st=round(SpaKernVec(m)/2)+1;
                   extrapOpts.disp.cropPx2nd=round(SpaKernVec(m)/2)+1;
                   extrapOpts.disp.extrapPx1st=extrapOpts.disp.cropPx1st;
                   extrapOpts.disp.extrapPx2nd=extrapOpts.disp.cropPx2nd;

                end
        
         % Adjust sample conditioning to remove at least one extrapolation
            %kernel from the cost function evaluation

        if extrapOpts.strain.extrapPx1st<=20
            CondOpts.FreeCens=20;
        else
            CondOpts.FreeCens=extrapOpts.strain.extrapPx1st;
        end

            end

        extrapOpts.strain.extrapPx1st=extrapOpts.disp.extrapPx1st;
        extrapOpts.strain.extrapPx2nd=extrapOpts.disp.extrapPx2nd;
        extrapOpts.strain.enforceGlobBCsPx=...
            extrapOpts.strain.extrapPx1st;
        extrapOpts.accel.extrapPx1st=extrapOpts.disp.extrapPx1st;
        etrapOpts.accel.extrapPx2nd=extrapOpts.disp.extrapPx2nd;

        accelCorr{m}=extrapOpts.accel;
        strainSmooth{m}=smoothOpts.strain;
        strainCorr{m}=extrapOpts.strain;
        dispCorr{m}=extrapOpts.disp;
end




%% initialize output variables
[DispNF.x,DispNF.y,DispBias.x,DispBias.y,...
    StrainNF.x,StrainNF.y,StrainNF.s,...
    StrainBias.x,StrainBias.y,StrainBias.s,...
    AccelNF.x,AccelNF.y,AccelBias.x,AccelBias.y,...
    SGNF.x,SGNF.s,SGBias.x,SGBias.s,...
    StressNF.x,StressNF.y,StressNF.s,...
    StressBias.x,StressBias.y,StressBias.s]=...
    deal(zeros(SpaKernNum,TempKernNum));

%% Generate the progress bar
itNum=0;
progmsg=strcat('0/',num2str(TotNum),' Smoothing Iterations Complete');
Progress=waitbar(0,progmsg);

%% Calculate the Noise Floor/resolution
for SK=1:SpaKernNum
    for TK=1:TempKernNum

    if TempKernVec(TK)==0
        smoothOpts.accel.temporalSmooth=false;

    else
        smoothOpts.accel.temporalSmooth=true;
        smoothOpts.accel.temporalKernelSize=TempKernVec(TK)*[1,1];
        smoothOpts.accel.temporalKernelOrder=[3,3];
    end   

    if smoothOpts.accel.temporalSmooth(1)==true
        CondOpts.cutStartFrames=...
            round(smoothOpts.accel.temporalKernelSize(1)/2);
        if CondOpts.cutStartFrames<=4
            CondOpts.cutEndFrames=4;
        else
            CondOpts.CutEndFrames=CondOpts.cutStartFrames;
        end
    end
     

    itNum=itNum+1;
    progmsg=strcat(num2str(itNum),'/',num2str(TotNum), ...
    ' Smoothing Iterations Complete');
        %% Set smoothing kernal sizes
        smoothOpts.strain=cell2mat(strainSmooth(SK));
        extrapOpts.disp=cell2mat(dispCorr(SK));
        extrapOpts.strain=cell2mat(strainCorr(SK));
        extrapOpts.accel=cell2mat(accelCorr(SK));



        NoiseMag=0.4;

        [DispRes,StrainRes,AccelRes,SGRes,StressRes,...
            ~,~,~,~,~] =...
            func_computeResolutionNoiseFloorV3(NoiseMag,grid,time, ...
            imagePath,imageFile,gridMethodOpts,specimen, ...
           CondOpts,smoothOpts,globalOpts, ...
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
        StressNF.x(SK,TK)=StressRes.x.NF;
        StressBias.x(SK,TK)=StressRes.x.avg;
        StressNF.y(SK,TK)=StressRes.y.NF;
        StressBias.y(SK,TK)=StressRes.y.avg;
        StressNF.s(SK,TK)=StressRes.s.NF;
        StressBias.s(SK,TK)=StressRes.s.avg;

        %Stress gauge noise floor and bias
        SGNF.x(SK,TK)=SGRes.x.NF;
        SGBias.x(SK,TK)=SGRes.x.avg;
        SGNF.s(SK,TK)=SGRes.s.NF;
        SGBias.s(SK,TK)=SGRes.s.avg;

        %Post progress bar
        Progress=waitbar(itNum/TotNum,Progress,progmsg);
    end
end

%% Set colormap
ColorScheme='hot';
%% Plot the displacement Noise floor and bias as a function 
figure('units','normalized','outerposition',[0 0 1 1])
%x direction noise floor
subplot(2,2,1)
contourf(SpaKernVec,TempKernVec,DispNF.x',50)
title('u_x Noise Floor')
xlabel('Spatial Kernal (px)')
ylabel('Temporal Kernal (frames)')

cx=colorbar;
colormap(ColorScheme);
cx.Label.String='stdev(u_x) (m)';

%x direction noise bias
subplot(2,2,2)
contourf(SpaKernVec,TempKernVec,DispBias.x',50)
title('u_x Noise Bias')
xlabel('Spatial Kernal (px)')
ylabel('Temporal Kernal (frames)')

cx=colorbar;
colormap(ColorScheme);
cx.Label.String='mean(u_x) (m)';

%y direction noise floor
subplot(2,2,3)
contourf(SpaKernVec,TempKernVec,DispNF.y',50)
title('u_y Noise Floor')
xlabel('Spatial Kernal (px)')
ylabel('Temporal Kernal (frames)')

cx=colorbar;
colormap(ColorScheme);
cx.Label.String='stdev(u_y) (m)';

%x direction noise bias
subplot(2,2,4)
contourf(SpaKernVec,TempKernVec,DispBias.y',50)
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
contourf(SpaKernVec,TempKernVec,AccelNF.x',50)
title('a_x Noise Floor')
xlabel('Spatial Kernal (px)')
ylabel('Temporal Kernal (frames)')

cx=colorbar;
colormap(ColorScheme);
cx.Label.String='stdev(a_x) (m/s^2)';

%x direction noise bias
subplot(2,2,2)
contourf(SpaKernVec,TempKernVec,AccelBias.x',50)
title('a_x Noise Bias')
xlabel('Spatial Kernal (px)')
ylabel('Temporal Kernal (frames)')

cx=colorbar;
colormap(ColorScheme);
cx.Label.String='mean(a_x) (m/s^2)';

%y direction noise floor
subplot(2,2,3)
contourf(SpaKernVec,TempKernVec,AccelNF.y',50)
title('a_y Noise Floor')
xlabel('Spatial Kernal (px)')
ylabel('Temporal Kernal (frames)')

cx=colorbar;
colormap(ColorScheme);
cx.Label.String='stdev(a_y) (m/s^2)';

%x direction noise bias
subplot(2,2,4)
contourf(SpaKernVec,TempKernVec,AccelBias.y',50)
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
contourf(SpaKernVec,TempKernVec,StrainNF.x',50)
title('$\varepsilon_{xx}$ Noise Floor','Interpreter','latex')
xlabel('Spatial Kernal (px)')
ylabel('Temporal Kernal (frames)')

cx=colorbar;
colormap(ColorScheme);
cx.Label.String='$stdev(\varepsilon_{xx})$';
cx.Label.Interpreter='latex';

%x direction noise bias
subplot(3,2,2)
contourf(SpaKernVec,TempKernVec,StrainBias.x',50)
title('$\varepsilon_{xx}$ Noise Bias','Interpreter','latex')
xlabel('Spatial Kernal (px)')
ylabel('Temporal Kernal (frames)')

cx=colorbar;
colormap(ColorScheme);
cx.Label.String='$mean(\varepsilon_{xx})$';
cx.Label.Interpreter='latex';

%y direction noise floor
subplot(3,2,3)
contourf(SpaKernVec,TempKernVec,StrainNF.y',50)
title('$\varepsilon_{yy}$ Noise Floor','Interpreter','latex')
xlabel('Spatial Kernal (px)')
ylabel('Temporal Kernal (frames)')

cx=colorbar;
colormap(ColorScheme);
cx.Label.String='$stdev(\varepsilon_{yy})$';
cx.Label.Interpreter='latex';

%y direction noise bias
subplot(3,2,4)
contourf(SpaKernVec,TempKernVec,StrainBias.y',50)
title('$\varepsilon_{yy}$ Noise Bias','Interpreter','latex')
xlabel('Spatial Kernal (px)')
ylabel('Temporal Kernal (frames)')

cx=colorbar;
colormap(ColorScheme);
cx.Label.String='$mean(\varepsilon_{yy})$';
cx.Label.Interpreter='latex';

%xy direction noise floor
subplot(3,2,5)
contourf(SpaKernVec,TempKernVec,StrainNF.s',50)
title('$\varepsilon_{xy}$ Noise Floor','Interpreter','latex')
xlabel('Spatial Kernal (px)')
ylabel('Temporal Kernal (frames)')

cx=colorbar;
colormap(ColorScheme);
cx.Label.String='$stdev(\varepsilon_{xy})$';
cx.Label.Interpreter='latex';

%y direction noise bias
subplot(3,2,6)
contourf(SpaKernVec,TempKernVec,StrainBias.s',50)
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
contourf(SpaKernVec,TempKernVec,StressNF.x',50)
title('$\sigma_{xx}$ Noise Floor','Interpreter','latex')
xlabel('Spatial Kernal (px)')
ylabel('Temporal Kernal (frames)')

cx=colorbar;
colormap(ColorScheme);
cx.Label.String='$stdev(\sigma_{xx})$';
cx.Label.Interpreter='latex';

%x direction noise bias
subplot(3,2,2)
contourf(SpaKernVec,TempKernVec,StressBias.x',50)
title('$\sigma_{xx}$ Noise Bias','Interpreter','latex')
xlabel('Spatial Kernal (px)')
ylabel('Temporal Kernal (frames)')

cx=colorbar;
colormap(ColorScheme);
cx.Label.String='$mean(\sigma_{xx})$';
cx.Label.Interpreter='latex';

%y direction noise floor
subplot(3,2,3)
contourf(SpaKernVec,TempKernVec,StressNF.y',50)
title('$\sigma_{yy}$ Noise Floor','Interpreter','latex')
xlabel('Spatial Kernal (px)')
ylabel('Temporal Kernal (frames)')

cx=colorbar;
colormap(ColorScheme);
cx.Label.String='$stdev(\sigma_{yy})$';
cx.Label.Interpreter='latex';

%y direction noise bias
subplot(3,2,4)
contourf(SpaKernVec,TempKernVec,StressBias.y',50)
title('$\sigma_{yy}$ Noise Bias','Interpreter','latex')
xlabel('Spatial Kernal (px)')
ylabel('Temporal Kernal (frames)')

cx=colorbar;
colormap(ColorScheme);
cx.Label.String='$mean(\sigma_{yy})$';
cx.Label.Interpreter='latex';

%xy direction noise floor
subplot(3,2,5)
contourf(SpaKernVec,TempKernVec,StressNF.s',50)
title('$\sigma_{xy}$ Noise Floor','Interpreter','latex')
xlabel('Spatial Kernal (px)')
ylabel('Temporal Kernal (frames)')

cx=colorbar;
colormap(ColorScheme);
cx.Label.String='$stdev(\sigma_{xy})$';
cx.Label.Interpreter='latex';

%y direction noise bias
subplot(3,2,6)
contourf(SpaKernVec,TempKernVec,StressBias.s',50)
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
contourf(SpaKernVec,TempKernVec,SGNF.x',50)
title('$\sigma_{xx_{sg}}$ Noise Floor','Interpreter','latex')
xlabel('Spatial Kernal (px)')
ylabel('Temporal Kernal (frames)')

cx=colorbar;
colormap(ColorScheme);
cx.Label.String='$stdev(\sigma_{xx})$';
cx.Label.Interpreter='latex';

%x direction noise bias
subplot(2,2,2)
contourf(SpaKernVec,TempKernVec,SGBias.x',50)
title('$\sigma_{xx_{sg}}$ Noise Bias','Interpreter','latex')
xlabel('Spatial Kernal (px)')
ylabel('Temporal Kernal (frames)')

cx=colorbar;
colormap(ColorScheme);
cx.Label.String='$mean(\sigma_{xx})$';
cx.Label.Interpreter='latex';

%xy direction noise floor
subplot(2,2,3)
contourf(SpaKernVec,TempKernVec,SGNF.s',50)
title('$\sigma_{xy_{sg}}$ Noise Floor','Interpreter','latex')
xlabel('Spatial Kernal (px)')
ylabel('Temporal Kernal (frames)')

cx=colorbar;
colormap(ColorScheme);
cx.Label.String='$stdev(\sigma_{xy})$';
cx.Label.Interpreter='latex';

%y direction noise bias
subplot(2,2,4)
contourf(SpaKernVec,TempKernVec,SGBias.s',50)
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
