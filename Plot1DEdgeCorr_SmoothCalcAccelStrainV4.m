% This script was created to evaluate the strain and displacement 
    %corrections at the edges for func_cropAndExtrapDispFields_v4
    %func_smoothCalcAccel_v4, func_smoothCalcStrain_v4, and their
    %associated subfunctions with and without smoothing

%%
clear variables
close all
clc

%% Load finite element data
[FEname,FEpath]=uigetfile('*.mat', ...
    'Choose file containing finite element kinematic fields');
FEfile=strcat(FEpath,'/',FEname);
FE=load(FEfile);

%% Choose folder to save results
SaveDir=uigetdir('','Choose Folder in which to save results');

%% Define Test Dseignation
TestDeg=char(cell2mat(inputdlg('Input Test Designation')));

%% Make subdirectories to save the image frames 
FigDir=strcat(SaveDir,'/FIG/');
PngDir=strcat(SaveDir,'/PNG/');

mkdir(FigDir);
mkdir(PngDir);

%% LOAD RAW DATA FILES: Raw .tiff images
    hardCodePath=false;
    fprintf('Loading reference image from the selected test data folder.\n')
    if ~hardCodePath
        [imageFile,imagePath] = uigetfile({'*.*','All Files'},'Select the first image in the sequence');
    else   
        imageFile = 'DefGridImage_001.tiff';
        imagePath = 'E:\Matlab_WorkingDirectory\1_IBIITest_Data\TestData_ID_CFIPQIso_HalfPulse\'; 
    end

%% Load Processing Parameter Data Structures
fprintf('Loading processing parameters file.\n')

    [initFile,initPath,~] = uigetfile('*.mat', ...
        'Locate processing parameter file for original calculations');

% Load the processing parameters from file


load([initPath,initFile])

% Store the FE valid mode parameter
globalOpts.FEValidMode = false;
globalOpts.calcKinFieldsFromDisp = true;

extrapOptsO=extrapOpts;
diffOptsO=diffOpts;

%% Load V4 Processing Parameters
fprintf('Loading processing parameters file.\n')

    [initFile,initPath,~] = uigetfile('*.mat', ...
        'Locate V4 processing parameter file for original calculations');

% Load the processing parameters from file


load([initPath,initFile])

% Store the FE valid mode parameter
globalOpts.FEValidMode = false;
globalOpts.calcKinFieldsFromDisp = true;

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


%% Determine weather to run interpolated correction on the displacement 
    %fields
    
fprintf('Determining dispacement filed correction options \n')
DispCorr.Opt='Yes';
DispCorr.int=10;
DispCorr.Method='Quad';
DispCorr.PitchFitKern=2;
DispCorr.strainMethod='GridPeriod';
DispCorr.strainPitchNum=2;

%% Smoothing opts for original strain and acceleration code
smoothingOpts.spatialSmooth=false;
smoothingOpts.FFTempSmooth=false;
smoothingOpts.WATempSmooth=false;
smoothingOpts.spatialFilt='gauss';
smoothingOpts.spatialKernal=[11,11];
smoothingOpts.spatialEdgeMode='symmetric';
smoothingOpts.FFTemporalFilt='sgolay';
smoothingOpts.FFTemporalKernal=[11,3];
smoothingOpts.FFTemporalPad=true;
smoothingOpts.FFTemporalPadFrames=3;
smoothingOpts.FFTemporalPadMethod='replicate';
smoothingOpts.WATemporalAvgFirst=false;
smoothingOpts.WATemporalFilt='sgolay';
smoothingOpts.WATemporalKernal=[11,3];
smoothingOpts.WATemporalPad=false;
smoothingOpts.WATemporalPadFrames=3;
smoothingOpts.WATemporalPadMethod='replicate';

%%
smoothOpts.strain.spatialSmooth=0;
smoothOpts.accel.temporalSmooth=0;

%% %% IMAGE PROCESSING: Use the Grid Method to extract displacement fields
FEValidMode=false;
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

            % Create 3D arrays of position to simplify generalised stress-strain curve
            % calculations
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

                    fprintf('Saving Raw Displacement Fields \n')
                    RawDisp=disp;
                    fprintf('Correcting Grid Method Displacements along specimen edges')
                                       
                    [disp.x,disp.y,disp.rX,disp.rY] = ...
                        func_cropAndExtrapDispFields_v4(pos,RawDisp.x, ...
                        RawDisp.y, ...
                        extrapOpts.disp, ...
                        1); %flag for printing to console
            
            %convert corrected displacements back to structure form
%             disp.x=dispExtrapX;
%             disp.y=dispExtrapY;
%             disp.rX=dispRangeX;
%             disp.rY=dispRangeY;
%             %clear extrapolated displacement variables 
%             clear dispExtrapY dispExtrapX dispRangeY dispRangeX
            % Calculate the kinematic fields from displacement fields using
            % displacements from images or displacements from FE data
           
                %--------------------------------------------------------------------------
                % Load the Reference Image and Determine Where the Free Edge is
                %fprintf('Obtaining and setting the free edge location.\n')
                [freeEdge,specimen,disp] = func_getFreeEdge(...
                    globalOpts.hardCodeFreeEdge,...
                    imagePath,imageFile,specimen,disp);

                %--------------------------------------------------------------------------
                %% Smooth and Calculate Strain
                fprintf('Calculating strain from the displacement fields.\n')
                [~,strainRaw,~] = func_smoothCalcStrainV1(globalOpts, ...
                    pos,time,...
                    grid,disp,smoothingOpts,extrapOptsO);

                [strain,disp] = func_smoothCalcStrain_v4(pos,...
                    time,disp,smoothOpts.strain,extrapOpts.strain,true);
                %% Calculate accelerations
                [~,~,RawAccel] = func_smoothCalcAccel(pos,time,grid, ...
                    RawDisp,smoothingOpts,...
                    extrapOptsO,diffOptsO);
                [accel,~,~] = func_smoothCalcAccel_v4(pos,time,disp,...
                    smoothOpts.accel,extrapOpts.accel,diffOpts,true);

                
  
%% Get indexes for middle and 4 grid pitches from free and impact edges

  GMPos.NumPoints=size(disp.x,2);
    GMPos.Free=20; %4 pitches from free edge
    GMPos.Mid=round(GMPos.NumPoints/2);
    GMPos.Imp=GMPos.NumPoints-20; %4 pitches from impact edge
    % Get cooridinates of position indexes
    GMPos.Xfree=pos.x(GMPos.Free);
    GMPos.Xmid=pos.x(GMPos.Mid);
    GMPos.Ximp=pos.x(GMPos.Imp);

    %Set Position indexes and coordintes for finite element data
    FEPos.NumPoints=size(FE.disp.x,2);
    FEPos.GridRat=FEPos.NumPoints/GMPos.NumPoints;
    FEPos.Free=round(GMPos.Free*FEPos.GridRat); %4 pitches from free edge
    FEPos.Mid=round(FEPos.NumPoints/2);
    FEPos.Imp=round(GMPos.Imp*FEPos.GridRat); %4 pitches from impact edge
    
    FEPos.Xfree=FE.pos.x(FEPos.Free);
    FEPos.Xmid=FE.pos.x(FEPos.Mid);
    FEPos.Ximp=FE.pos.x(FEPos.Imp);
    
  %Print Plot selections for evaluation purposes
  FEPos.freeString=num2str(FEPos.Xfree);
  GMPos.freeString=num2str(GMPos.Xfree);
  fprintf(strcat('Free Coordinate FE:',FEPos.freeString,' Grid: ',...
      GMPos.freeString,'\n'));
  
  FEPos.midString=num2str(FEPos.Xmid);
  GMPos.midString=num2str(GMPos.Xmid);
  fprintf(strcat('Mid Coordinate FE:',FEPos.midString,' Grid: ',...
      GMPos.midString,'\n'));

  FEPos.impString=num2str(FEPos.Ximp);
  GMPos.impString=num2str(GMPos.Ximp);
  fprintf(strcat('Impact Coordinate FE:',FEPos.impString,' Grid: ',...
      GMPos.impString,'\n'));

%% plot Y displacement curves
figure('Units','normalized','OuterPosition',[0,0,1,1])
timems=time.vec*10^6;
Ymm=pos.y*10^3;
YmmFE=FE.pos.y*10^3;
% Set plot limits
YdispLim=1.1*[min(FE.disp.y,[],'all'),max(FE.disp.y,[],'all')];

for k=1:length(time.vec)
    Fnum=num2str(k);
    Fstring=strcat('Frame~',Fnum);
    
    %% Extract vector
    %FE Model Displacement
    FEDyImp=squeeze(FE.disp.y(:,FEPos.Imp,k));
    FEDyMid=squeeze(FE.disp.y(:,FEPos.Mid,k));
    FEDyFree=squeeze(FE.disp.y(:,FEPos.Free,k));
    
    %Uncorrected Grid Method Displacement
    RawDyImp=squeeze(RawDisp.y(:,GMPos.Imp,k));
    RawDyMid=squeeze(RawDisp.y(:,GMPos.Mid,k));
    RawDyFree=squeeze(RawDisp.y(:,GMPos.Free,k));


    %Corrected Grid Method Dipslacement
    CorrDyImp=squeeze(disp.y(:,GMPos.Imp,k));
    CorrDyMid=squeeze(disp.y(:,GMPos.Mid,k));
    CorrDyFree=squeeze(disp.y(:,GMPos.Free,k));

    %%  Y Displacement at the impact edge
    subplot(3,1,1)
    plot(YmmFE,FEDyImp,'k')
    hold on
    plot(Ymm,RawDyImp,'b')
    plot(Ymm,CorrDyImp,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('u_y (m)')
    title(strcat('4P from Impact ',Fstring))
    legend('FE','RawGM','Corrected','location','southwest')
    ylim(YdispLim)

    %% Y Displacement in middle
    subplot(3,1,2)
    plot(YmmFE,FEDyMid,'k')
    hold on
    plot(Ymm,RawDyMid,'b')
    plot(Ymm,CorrDyMid,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('u_y (m)')
    title(strcat('Middle ',Fstring))
    ylim(YdispLim)

    %% Y displacement Free Surface
    subplot(3,1,3)
    plot(YmmFE,FEDyFree,'k')
    hold on
    plot(Ymm,RawDyFree,'b')
    plot(Ymm,CorrDyFree,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('u_y (m)')
    title(strcat('4P from Free Edge ',Fstring))
    ylim(YdispLim)

    %% SaveFiles
    FigName=strcat(FigDir,TestDeg,'_1DyD_',Fnum,'.fig');
    saveas(gcf,FigName);
    PngName=strcat(PngDir,TestDeg,'_1DyD_',Fnum,'.png');
    saveas(gcf,PngName);

end

%% Grab top and bottom indices
GMindVec=1:length(pos.y);
FEindVec=1:length(FE.pos.y);

%top
GMTind=GMindVec((end-5*grid.pxPerPeriod+1):end);
GMTcoord=pos.y(GMTind);
FETind=FEindVec(FE.pos.y>=GMTcoord(1));
FETcoord=FE.pos.y(FETind);

%bottom
GMBind=GMindVec(1:5*grid.pxPerPeriod);
GMBcoord=pos.y(GMBind);
FEBind=FEindVec(FE.pos.y<=GMBcoord(end));
FEBcoord=FE.pos.y(FEBind);

%% Create Directory
ZoomPngDir=strcat(PngDir,'/ZoomYdisp/');
ZoomFigDir=strcat(FigDir,'/ZoomYdisp/');
mkdir(ZoomPngDir);
mkdir(ZoomFigDir);
%% Plot the zoomed displacements

%Determine limits
TDyLim=1.1*[min(FE.disp.y(FETind,:,:),[],'all'),...
    max(FE.disp.y(FETind,:,:),[],'all')];
BDyLim=1.1*[min(FE.disp.y(FEBind,:,:),[],'all'),...
    max(FE.disp.y(FEBind,:,:),[],'all')];

%Initialize figure
figure('Units','normalized','outerposition',[0,0,1,1])

for k=1:length(time.vec)
    Fnum=num2str(k);
    Fstring=strcat('Frame~',Fnum);

        %% Extract vector
    %FE Model Displacement
    FEDyImpT=squeeze(FE.disp.y(FETind,FEPos.Imp,k));
    FEDyMidT=squeeze(FE.disp.y(FETind,FEPos.Mid,k));
    FEDyFreeT=squeeze(FE.disp.y(FETind,FEPos.Free,k));
    
    FEDyImpB=squeeze(FE.disp.y(FEBind,FEPos.Imp,k));
    FEDyMidB=squeeze(FE.disp.y(FEBind,FEPos.Mid,k));
    FEDyFreeB=squeeze(FE.disp.y(FEBind,FEPos.Free,k));
    
    %Uncorrected Grid Method Displacement
    RawDyImpT=squeeze(RawDisp.y(GMTind,GMPos.Imp,k));
    RawDyMidT=squeeze(RawDisp.y(GMTind,GMPos.Mid,k));
    RawDyFreeT=squeeze(RawDisp.y(GMTind,GMPos.Free,k));

    RawDyImpB=squeeze(RawDisp.y(GMBind,GMPos.Imp,k));
    RawDyMidB=squeeze(RawDisp.y(GMBind,GMPos.Mid,k));
    RawDyFreeB=squeeze(RawDisp.y(GMBind,GMPos.Free,k));

    %Corrected Grid Method Dipslacement
    CorrDyImpT=squeeze(disp.y(GMTind,GMPos.Imp,k));
    CorrDyMidT=squeeze(disp.y(GMTind,GMPos.Mid,k));
    CorrDyFreeT=squeeze(disp.y(GMTind,GMPos.Free,k));

    CorrDyImpB=squeeze(disp.y(GMBind,GMPos.Imp,k));
    CorrDyMidB=squeeze(disp.y(GMBind,GMPos.Mid,k));
    CorrDyFreeB=squeeze(disp.y(GMBind,GMPos.Free,k));

    %% Top Edges
    %%  Y Displacement at the impact edge
    subplot(2,3,1)
    plot(FETcoord*10^3,FEDyImpT,'k')
    hold on
    plot(GMTcoord*10^3,RawDyImpT,'b')
    plot(GMTcoord*10^3,CorrDyImpT,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('u_y (m)')
    title(strcat('Top 4P from Impact ',Fstring))
    legend('FE','RawGM','Corrected','location','southwest')
    ylim(TDyLim)

    %% Y Displacement in middle
    subplot(2,3,2)
    plot(FETcoord*10^3,FEDyMidT,'k')
    hold on
    plot(GMTcoord*10^3,RawDyMidT,'b')
    plot(GMTcoord*10^3,CorrDyMidT,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('u_y (m)')
    title(strcat('Top Middle ',Fstring))
    ylim(TDyLim)

    %% Y displacement Free Surface
    subplot(2,3,3)
    plot(FETcoord*10^3,FEDyFreeT,'k')
    hold on
    plot(GMTcoord*10^3,RawDyFreeT,'b')
    plot(GMTcoord*10^3,CorrDyFreeT,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('u_y (m)')
    title(strcat('Top 4P from Free Edge ',Fstring))
    ylim(TDyLim)

    %%Bottom Edges
    %%  Y Displacement at the impact edge
    subplot(2,3,4)
    plot(FEBcoord*10^3,FEDyImpB,'k')
    hold on
    plot(GMBcoord*10^3,RawDyImpB,'b')
    plot(GMBcoord*10^3,CorrDyImpB,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('u_y (m)')
    title(strcat('Bottom 4P from Impact ',Fstring))
    legend('FE','RawGM','Corrected','location','southwest')
    ylim(BDyLim)

    %% Y Displacement in middle
    subplot(2,3,5)
    plot(FEBcoord*10^3,FEDyMidB,'k')
    hold on
    plot(GMBcoord*10^3,RawDyMidB,'b')
    plot(GMBcoord*10^3,CorrDyMidB,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('u_y (m)')
    title(strcat('Bottom Middle ',Fstring))
    ylim(BDyLim)

    %% Y displacement Free Surface
    subplot(2,3,6)
    plot(FEBcoord*10^3,FEDyFreeB,'k')
    hold on
    plot(GMBcoord*10^3,RawDyFreeB,'b')
    plot(GMBcoord*10^3,CorrDyFreeB,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('u_y (m)')
    title(strcat('Bottom 4P from Free Edge ',Fstring))
    ylim(BDyLim)

    %% SaveFiles
    FigName=strcat(ZoomFigDir,TestDeg,'_1DyDZoom_',Fnum,'.fig');
    saveas(gcf,FigName);
    PngName=strcat(ZoomPngDir,TestDeg,'_1DyDZoom_',Fnum,'.png');
    saveas(gcf,PngName);
    
end

%% X displacement
%Make directories
PngDirXdisp=strcat(PngDir,'/Xdisp/');
FigDirXdisp=strcat(FigDir,'/Xdisp/');
mkdir(PngDirXdisp);
mkdir(FigDirXdisp);
%% plot X displacement curves
figure('Units','normalized','OuterPosition',[0,0,1,1])
timems=time.vec*10^6;
Ymm=pos.y*10^3;
YmmFE=FE.pos.y*10^3;
% Set plot limits
XdispLim=1.1*[min(FE.disp.x,[],'all'),max(FE.disp.x,[],'all')];

for k=1:length(time.vec)
    Fnum=num2str(k);
    Fstring=strcat('Frame~',Fnum);
    
    %% Extract vector
    %FE Model Displacement
    FEDxImp=squeeze(FE.disp.x(:,FEPos.Imp,k));
    FEDxMid=squeeze(FE.disp.x(:,FEPos.Mid,k));
    FEDxFree=squeeze(FE.disp.x(:,FEPos.Free,k));
    
    %Uncorrected Grid Method Displacement
    RawDxImp=squeeze(RawDisp.x(:,GMPos.Imp,k));
    RawDxMid=squeeze(RawDisp.x(:,GMPos.Mid,k));
    RawDxFree=squeeze(RawDisp.x(:,GMPos.Free,k));


    %Corrected Grid Method Dipslacement
    CorrDxImp=squeeze(disp.x(:,GMPos.Imp,k));
    CorrDxMid=squeeze(disp.x(:,GMPos.Mid,k));
    CorrDxFree=squeeze(disp.x(:,GMPos.Free,k));

    %%  x Displacement at the impact edge
    subplot(3,1,1)
    plot(YmmFE,FEDxImp,'k')
    hold on
    plot(Ymm,RawDxImp,'b')
    plot(Ymm,CorrDxImp,'r')
    hold off
    xlabel('y Coordinate (mm)')
    ylabel('u_x (m)')
    title(strcat('4P from Impact ',Fstring))
    legend('FE','RawGM','Corrected','location','southwest')
    ylim(XdispLim)

    %% X Displacement in middle
    subplot(3,1,2)
    plot(YmmFE,FEDxMid,'k')
    hold on
    plot(Ymm,RawDxMid,'b')
    plot(Ymm,CorrDxMid,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('u_x (m)')
    title(strcat('Middle ',Fstring))
    ylim(XdispLim)

    %% X displacement Free Surface
    subplot(3,1,3)
    plot(YmmFE,FEDxFree,'k')
    hold on
    plot(Ymm,RawDxFree,'b')
    plot(Ymm,CorrDxFree,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('u_x (m)')
    title(strcat('4P from Free Edge ',Fstring))
    ylim(XdispLim)

    %% SaveFiles
    FigName=strcat(FigDirXdisp,TestDeg,'_1DxD_',Fnum,'.fig');
    saveas(gcf,FigName);
    PngName=strcat(PngDirXdisp,TestDeg,'_1DxD_',Fnum,'.png');
    saveas(gcf,PngName);

end

%% Grab top and bottom indices
GMindVec=1:length(pos.y);
FEindVec=1:length(FE.pos.y);

%top
GMTind=GMindVec((end-5*grid.pxPerPeriod+1):end);
GMTcoord=pos.y(GMTind);
FETind=FEindVec(FE.pos.y>=GMTcoord(1));
FETcoord=FE.pos.y(FETind);

%bottom
GMBind=GMindVec(1:5*grid.pxPerPeriod);
GMBcoord=pos.y(GMBind);
FEBind=FEindVec(FE.pos.y<=GMBcoord(end));
FEBcoord=FE.pos.y(FEBind);

%% Create Directory
ZoomPngDirX=strcat(PngDir,'/ZoomXdisp/');
ZoomFigDirX=strcat(FigDir,'/ZoomXdisp/');
mkdir(ZoomPngDirX);
mkdir(ZoomFigDirX);
%% Plot the zoomed displacements

%Determine limits
TDxLim=1.1*[min(FE.disp.x(FETind,:,:),[],'all'),...
    max(FE.disp.x(FETind,:,:),[],'all')];
BDxLim=1.1*[min(FE.disp.x(FEBind,:,:),[],'all'),...
    max(FE.disp.x(FEBind,:,:),[],'all')];

%Initialize figure
figure('Units','normalized','outerposition',[0,0,1,1])

for k=1:length(time.vec)
    Fnum=num2str(k);
    Fstring=strcat('Frame~',Fnum);

        %% Extract vector
    %FE Model Displacement
    FEDxImpT=squeeze(FE.disp.x(FETind,FEPos.Imp,k));
    FEDxMidT=squeeze(FE.disp.x(FETind,FEPos.Mid,k));
    FEDxFreeT=squeeze(FE.disp.x(FETind,FEPos.Free,k));
    
    FEDxImpB=squeeze(FE.disp.x(FEBind,FEPos.Imp,k));
    FEDxMidB=squeeze(FE.disp.x(FEBind,FEPos.Mid,k));
    FEDxFreeB=squeeze(FE.disp.x(FEBind,FEPos.Free,k));
    
    %Uncorrected Grid Method Displacement
    RawDxImpT=squeeze(RawDisp.x(GMTind,GMPos.Imp,k));
    RawDxMidT=squeeze(RawDisp.x(GMTind,GMPos.Mid,k));
    RawDxFreeT=squeeze(RawDisp.x(GMTind,GMPos.Free,k));

    RawDxImpB=squeeze(RawDisp.x(GMBind,GMPos.Imp,k));
    RawDxMidB=squeeze(RawDisp.x(GMBind,GMPos.Mid,k));
    RawDxFreeB=squeeze(RawDisp.x(GMBind,GMPos.Free,k));

    %Corrected Grid Method Dipslacement
    CorrDxImpT=squeeze(disp.x(GMTind,GMPos.Imp,k));
    CorrDxMidT=squeeze(disp.x(GMTind,GMPos.Mid,k));
    CorrDxFreeT=squeeze(disp.x(GMTind,GMPos.Free,k));

    CorrDxImpB=squeeze(disp.x(GMBind,GMPos.Imp,k));
    CorrDxMidB=squeeze(disp.x(GMBind,GMPos.Mid,k));
    CorrDxFreeB=squeeze(disp.x(GMBind,GMPos.Free,k));

    %% Top Edges
    %%  Y Displacement at the impact edge
    subplot(2,3,1)
    plot(FETcoord*10^3,FEDxImpT,'k')
    hold on
    plot(GMTcoord*10^3,RawDxImpT,'b')
    plot(GMTcoord*10^3,CorrDxImpT,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('u_x (m)')
    title(strcat('Top 4P from Impact ',Fstring))
    legend('FE','RawGM','Corrected','location','southwest')
    ylim(TDxLim)

    %% Y Displacement in middle
    subplot(2,3,2)
    plot(FETcoord*10^3,FEDxMidT,'k')
    hold on
    plot(GMTcoord*10^3,RawDxMidT,'b')
    plot(GMTcoord*10^3,CorrDxMidT,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('u_x (m)')
    title(strcat('Top Middle ',Fstring))
    ylim(TDxLim)

    %% Y displacement Free Surface
    subplot(2,3,3)
    plot(FETcoord*10^3,FEDxFreeT,'k')
    hold on
    plot(GMTcoord*10^3,RawDxFreeT,'b')
    plot(GMTcoord*10^3,CorrDxFreeT,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('u_x (m)')
    title(strcat('Top 4P from Free Edge ',Fstring))
    ylim(TDxLim)

    %%Bottom Edges
    %%  Y Displacement at the impact edge
    subplot(2,3,4)
    plot(FEBcoord*10^3,FEDxImpB,'k')
    hold on
    plot(GMBcoord*10^3,RawDxImpB,'b')
    plot(GMBcoord*10^3,CorrDxImpB,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('u_x (m)')
    title(strcat('Bottom 4P from Impact ',Fstring))
    legend('FE','RawGM','Corrected','location','southwest')
    ylim(BDxLim)

    %% Y Displacement in middle
    subplot(2,3,5)
    plot(FEBcoord*10^3,FEDxMidB,'k')
    hold on
    plot(GMBcoord*10^3,RawDxMidB,'b')
    plot(GMBcoord*10^3,CorrDxMidB,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('u_x (m)')
    title(strcat('Bottom Middle ',Fstring))
    ylim(BDxLim)

    %% Y displacement Free Surface
    subplot(2,3,6)
    plot(FEBcoord*10^3,FEDxFreeB,'k')
    hold on
    plot(GMBcoord*10^3,RawDxFreeB,'b')
    plot(GMBcoord*10^3,CorrDxFreeB,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('u_x (m)')
    title(strcat('Bottom 4P from Free Edge ',Fstring))
    ylim(BDxLim)

    %% SaveFiles
    FigName=strcat(ZoomFigDirX,TestDeg,'_1DxDZoom_',Fnum,'.fig');
    saveas(gcf,FigName);
    PngName=strcat(ZoomPngDirX,TestDeg,'_1DxDZoom_',Fnum,'.png');
    saveas(gcf,PngName);
    
end

%% YY Strain directory
PngDirYY=strcat(PngDir,'/YYstrain/');
FigDirYY=strcat(FigDir,'/YYstrain/');
mkdir(PngDirYY);
mkdir(FigDirYY);
%% YY strain
figure('Units','normalized','OuterPosition',[0,0,1,1])
timems=time.vec*10^6;
Ymm=pos.y*10^3;
YmmFE=FE.pos.y*10^3;
% Set plot limits
YYstrainLim=1.1*[min(FE.strain.y,[],'all'),max(FE.strain.y,[],'all')];

for k=1:length(time.vec)
    Fnum=num2str(k);
    Fstring=strcat('Frame~',Fnum);
    
    %% Extract vector
    %FE Model Displacement
    FEYYImp=squeeze(FE.strain.y(:,FEPos.Imp,k));
    FEYYMid=squeeze(FE.strain.y(:,FEPos.Mid,k));
    FEYYFree=squeeze(FE.strain.y(:,FEPos.Free,k));
    
    %Uncorrected Grid Method Displacement
    RawYYImp=squeeze(strainRaw.y(:,GMPos.Imp,k));
    RawYYMid=squeeze(strainRaw.y(:,GMPos.Mid,k));
    RawYYFree=squeeze(strainRaw.y(:,GMPos.Free,k));


    %Corrected Grid Method Dipslacement
    CorrYYImp=squeeze(strain.y(:,GMPos.Imp,k));
    CorrYYMid=squeeze(strain.y(:,GMPos.Mid,k));
    CorrYYFree=squeeze(strain.y(:,GMPos.Free,k));

    %%  Y Displacement at the impact edge
    subplot(3,1,1)
    plot(YmmFE,FEYYImp,'k')
    hold on
    plot(Ymm,RawYYImp,'b')
    plot(Ymm,CorrYYImp,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('strain_{yy}')
    title(strcat('4P from Impact ',Fstring))
    legend('FE','RawGM','Corrected','location','southwest')
    ylim(YYstrainLim)

    %% Y Displacement in middle
    subplot(3,1,2)
    plot(YmmFE,FEYYMid,'k')
    hold on
    plot(Ymm,RawYYMid,'b')
    plot(Ymm,CorrYYMid,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('strain_{yy}')
    title(strcat('Middle ',Fstring))
    ylim(YYstrainLim)

    %% Y displacement Free Surface
    subplot(3,1,3)
    plot(YmmFE,FEYYFree,'k')
    hold on
    plot(Ymm,RawYYFree,'b')
    plot(Ymm,CorrYYFree,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('strain_{yy}')
    title(strcat('4P from Free Edge ',Fstring))
    ylim(YYstrainLim)

    %% SaveFiles
    FigName=strcat(FigDirYY,TestDeg,'_1DstrainY_',Fnum,'.fig');
    saveas(gcf,FigName);
    PngName=strcat(PngDirYY,TestDeg,'_1DstrainY_',Fnum,'.png');
    saveas(gcf,PngName);

end

%% XY Strain directory
PngDirXY=strcat(PngDir,'/Shearstrain/');
FigDirXY=strcat(FigDir,'/Shearstrain/');
mkdir(PngDirXY);
mkdir(FigDirXY);
%% XY strain
figure('Units','normalized','OuterPosition',[0,0,1,1])
timems=time.vec*10^6;
Ymm=pos.y*10^3;
YmmFE=FE.pos.y*10^3;
% Set plot limits
XYstrainLim=1.1*[min(FE.strain.s,[],'all'),max(FE.strain.s,[],'all')];

for k=1:length(time.vec)
    Fnum=num2str(k);
    Fstring=strcat('Frame~',Fnum);
    
    %% Extract vector
    %FE Model Displacement
    FEXYImp=squeeze(FE.strain.s(:,FEPos.Imp,k));
    FEXYMid=squeeze(FE.strain.s(:,FEPos.Mid,k));
    FEXYFree=squeeze(FE.strain.s(:,FEPos.Free,k));
    
    %Uncorrected Grid Method Displacement
    RawXYImp=squeeze(strainRaw.s(:,GMPos.Imp,k));
    RawXYMid=squeeze(strainRaw.s(:,GMPos.Mid,k));
    RawXYFree=squeeze(strainRaw.s(:,GMPos.Free,k));


    %Corrected Grid Method Dipslacement
    CorrXYImp=squeeze(strain.s(:,GMPos.Imp,k));
    CorrXYMid=squeeze(strain.s(:,GMPos.Mid,k));
    CorrXYFree=squeeze(strain.s(:,GMPos.Free,k));

    %%  Shear strain at the impact edge
    subplot(3,1,1)
    plot(YmmFE,FEXYImp,'k')
    hold on
    plot(Ymm,RawXYImp,'b')
    plot(Ymm,CorrXYImp,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('strain_{xy}')
    title(strcat('4P from Impact ',Fstring))
    legend('FE','RawGM','Corrected','location','southwest')
    ylim(XYstrainLim)

    %% Shear strain in middle
    subplot(3,1,2)
    plot(YmmFE,FEXYMid,'k')
    hold on
    plot(Ymm,RawXYMid,'b')
    plot(Ymm,CorrXYMid,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('strain_{xy}')
    title(strcat('Middle ',Fstring))
    ylim(XYstrainLim)

    %% Shear strain Free Surface
    subplot(3,1,3)
    plot(YmmFE,FEXYFree,'k')
    hold on
    plot(Ymm,RawXYFree,'b')
    plot(Ymm,CorrXYFree,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('strain_{xy}')
    title(strcat('4P from Free Edge ',Fstring))
    ylim(XYstrainLim)

    %% SaveFiles
    FigName=strcat(FigDirXY,TestDeg,'_1DstrainShear_',Fnum,'.fig');
    saveas(gcf,FigName);
    PngName=strcat(PngDirXY,TestDeg,'_1DstrainShear_',Fnum,'.png');
    saveas(gcf,PngName);

end

%% XX Strain directory
PngDirXX=strcat(PngDir,'/XXstrain/');
FigDirXX=strcat(FigDir,'/XXstrain/');
mkdir(PngDirXX);
mkdir(FigDirXX);
%% XX strain
figure('Units','normalized','OuterPosition',[0,0,1,1])
timems=time.vec*10^6;
Ymm=pos.y*10^3;
YmmFE=FE.pos.y*10^3;
% Set plot limits
XXstrainLim=1.1*[min(FE.strain.x,[],'all'),max(FE.strain.x,[],'all')];

for k=1:length(time.vec)
    Fnum=num2str(k);
    Fstring=strcat('Frame~',Fnum);
    
    %% Extract vector
    %FE Model Displacement
    FEXXImp=squeeze(FE.strain.x(:,FEPos.Imp,k));
    FEXXMid=squeeze(FE.strain.x(:,FEPos.Mid,k));
    FEXXFree=squeeze(FE.strain.x(:,FEPos.Free,k));
    
    %Uncorrected Grid Method Displacement
    RawXXImp=squeeze(strainRaw.x(:,GMPos.Imp,k));
    RawXXMid=squeeze(strainRaw.x(:,GMPos.Mid,k));
    RawXXFree=squeeze(strainRaw.x(:,GMPos.Free,k));


    %Corrected Grid Method Dipslacement
    CorrXXImp=squeeze(strain.x(:,GMPos.Imp,k));
    CorrXXMid=squeeze(strain.x(:,GMPos.Mid,k));
    CorrXXFree=squeeze(strain.x(:,GMPos.Free,k));

    %%  XX strain at the impact edge
    subplot(3,1,1)
    plot(YmmFE,FEXXImp,'k')
    hold on
    plot(Ymm,RawXXImp,'b')
    plot(Ymm,CorrXXImp,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('strain_{xx}')
    title(strcat('4P from Impact ',Fstring))
    legend('FE','RawGM','Corrected','location','southwest')
    ylim(XXstrainLim)

    %% Y Displacement in middle
    subplot(3,1,2)
    plot(YmmFE,FEXXMid,'k')
    hold on
    plot(Ymm,RawXXMid,'b')
    plot(Ymm,CorrXXMid,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('strain_{xx}')
    title(strcat('Middle ',Fstring))
    ylim(XXstrainLim)

    %% Y displacement Free Surface
    subplot(3,1,3)
    plot(YmmFE,FEXXFree,'k')
    hold on
    plot(Ymm,RawXXFree,'b')
    plot(Ymm,CorrXXFree,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('strain_{xx}')
    title(strcat('4P from Free Edge ',Fstring))
    ylim(XXstrainLim)

    %% SaveFiles
    FigName=strcat(FigDirXX,TestDeg,'_1DstrainX_',Fnum,'.fig');
    saveas(gcf,FigName);
    PngName=strcat(PngDirXX,TestDeg,'_1DstrainX_',Fnum,'.png');
    saveas(gcf,PngName);

end


%% Create Directory
ZoomPngDir=strcat(PngDir,'/ZoomYYstrain/');
ZoomFigDir=strcat(FigDir,'/ZoomYYstrain/');
mkdir(ZoomPngDir);
mkdir(ZoomFigDir);
%% Plot the zoomed YY strains

%Determine limits
TYYLim=1.1*[min(FE.strain.y(FETind,:,:),[],'all'),...
    max(FE.strain.y(FETind,:,:),[],'all')];
BYYLim=1.1*[min(FE.strain.y(FEBind,:,:),[],'all'),...
    max(FE.strain.y(FEBind,:,:),[],'all')];

%Initialize figure
figure('Units','normalized','outerposition',[0,0,1,1])

for k=1:length(time.vec)
    Fnum=num2str(k);
    Fstring=strcat('Frame~',Fnum);

        %% Extract vector
    %FE Model Displacement
    FEYYImpT=squeeze(FE.strain.y(FETind,FEPos.Imp,k));
    FEYYMidT=squeeze(FE.strain.y(FETind,FEPos.Mid,k));
    FEYYFreeT=squeeze(FE.strain.y(FETind,FEPos.Free,k));
    
    FEYYImpB=squeeze(FE.strain.y(FEBind,FEPos.Imp,k));
    FEYYMidB=squeeze(FE.strain.y(FEBind,FEPos.Mid,k));
    FEYYFreeB=squeeze(FE.strain.y(FEBind,FEPos.Free,k));
    
    %Uncorrected Grid Method Displacement
    RawYYImpT=squeeze(strainRaw.y(GMTind,GMPos.Imp,k));
    RawYYMidT=squeeze(strainRaw.y(GMTind,GMPos.Mid,k));
    RawYYFreeT=squeeze(strainRaw.y(GMTind,GMPos.Free,k));

    RawYYImpB=squeeze(strainRaw.y(GMBind,GMPos.Imp,k));
    RawYYMidB=squeeze(strainRaw.y(GMBind,GMPos.Mid,k));
    RawYYFreeB=squeeze(strainRaw.y(GMBind,GMPos.Free,k));

    %Corrected Grid Method Dipslacement
    CorrYYImpT=squeeze(strain.y(GMTind,GMPos.Imp,k));
    CorrYYMidT=squeeze(strain.y(GMTind,GMPos.Mid,k));
    CorrYYFreeT=squeeze(strain.y(GMTind,GMPos.Free,k));

    CorrYYImpB=squeeze(strain.y(GMBind,GMPos.Imp,k));
    CorrYYMidB=squeeze(strain.y(GMBind,GMPos.Mid,k));
    CorrYYFreeB=squeeze(strain.y(GMBind,GMPos.Free,k));

    %% Top Edges
    %%  Y Displacement at the impact edge
    subplot(2,3,1)
    plot(FETcoord*10^3,FEYYImpT,'k')
    hold on
    plot(GMTcoord*10^3,RawYYImpT,'b')
    plot(GMTcoord*10^3,CorrYYImpT,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('strain_{yy}')
    title(strcat('Top 4P from Impact ',Fstring))
    legend('FE','RawGM','Corrected','location','southwest')
    ylim(TYYLim)

    %% Y Displacement in middle
    subplot(2,3,2)
    plot(FETcoord*10^3,FEYYMidT,'k')
    hold on
    plot(GMTcoord*10^3,RawYYMidT,'b')
    plot(GMTcoord*10^3,CorrYYMidT,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('strain_{yy}')
    title(strcat('Top Middle ',Fstring))
    ylim(TYYLim)

    %% Y displacement Free Surface
    subplot(2,3,3)
    plot(FETcoord*10^3,FEYYFreeT,'k')
    hold on
    plot(GMTcoord*10^3,RawYYFreeT,'b')
    plot(GMTcoord*10^3,CorrYYFreeT,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('strain_{yy}')
    title(strcat('Top 4P from Free Edge ',Fstring))
    ylim(TYYLim)

    %%Bottom Edges
    %%  Y Displacement at the impact edge
    subplot(2,3,4)
    plot(FEBcoord*10^3,FEYYImpB,'k')
    hold on
    plot(GMBcoord*10^3,RawYYImpB,'b')
    plot(GMBcoord*10^3,CorrYYImpB,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('strain_{yy}')
    title(strcat('Bottom 4P from Impact ',Fstring))
    legend('FE','RawGM','Corrected','location','southwest')
    ylim(BYYLim)

    %% Y Displacement in middle
    subplot(2,3,5)
    plot(FEBcoord*10^3,FEYYMidB,'k')
    hold on
    plot(GMBcoord*10^3,RawYYMidB,'b')
    plot(GMBcoord*10^3,CorrYYMidB,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('strain_{yy}')
    title(strcat('Bottom Middle ',Fstring))
    ylim(BYYLim)

    %% Y displacement Free Surface
    subplot(2,3,6)
    plot(FEBcoord*10^3,FEYYFreeB,'k')
    hold on
    plot(GMBcoord*10^3,RawYYFreeB,'b')
    plot(GMBcoord*10^3,CorrYYFreeB,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('strain_{yy}')
    title(strcat('Bottom 4P from Free Edge ',Fstring))
    ylim(BYYLim)

    %% SaveFiles
    FigName=strcat(ZoomFigDir,TestDeg,'_1DStrainYYZoom_',Fnum,'.fig');
    saveas(gcf,FigName);
    PngName=strcat(ZoomPngDir,TestDeg,'_1DStrainYYZoom_',Fnum,'.png');
    saveas(gcf,PngName);
    
end

%% Create Directory
ZoomPngDir=strcat(PngDir,'/ZoomXYstrain/');
ZoomFigDir=strcat(FigDir,'/ZoomXYstrain/');
mkdir(ZoomPngDir);
mkdir(ZoomFigDir);
%% Plot the zoomed Shear strains

%Determine limits
TXYLim=1.1*[min(FE.strain.s(FETind,:,:),[],'all'),...
    max(FE.strain.s(FETind,:,:),[],'all')];
BXYLim=1.1*[min(FE.strain.s(FEBind,:,:),[],'all'),...
    max(FE.strain.s(FEBind,:,:),[],'all')];

%Initialize figure
figure('Units','normalized','outerposition',[0,0,1,1])

for k=1:length(time.vec)
    Fnum=num2str(k);
    Fstring=strcat('Frame~',Fnum);

        %% Extract vector
    %FE Model Displacement
    FEXYImpT=squeeze(FE.strain.s(FETind,FEPos.Imp,k));
    FEXYMidT=squeeze(FE.strain.s(FETind,FEPos.Mid,k));
    FEXYFreeT=squeeze(FE.strain.s(FETind,FEPos.Free,k));
    
    FEXYImpB=squeeze(FE.strain.s(FEBind,FEPos.Imp,k));
    FEXYMidB=squeeze(FE.strain.s(FEBind,FEPos.Mid,k));
    FEXYFreeB=squeeze(FE.strain.s(FEBind,FEPos.Free,k));
    
    %Uncorrected Grid Method Displacement
    RawXYImpT=squeeze(strainRaw.s(GMTind,GMPos.Imp,k));
    RawXYMidT=squeeze(strainRaw.s(GMTind,GMPos.Mid,k));
    RawXYFreeT=squeeze(strainRaw.s(GMTind,GMPos.Free,k));

    RawXYImpB=squeeze(strainRaw.s(GMBind,GMPos.Imp,k));
    RawXYMidB=squeeze(strainRaw.s(GMBind,GMPos.Mid,k));
    RawXYFreeB=squeeze(strainRaw.s(GMBind,GMPos.Free,k));

    %Corrected Grid Method Dipslacement
    CorrXYImpT=squeeze(strain.s(GMTind,GMPos.Imp,k));
    CorrXYMidT=squeeze(strain.s(GMTind,GMPos.Mid,k));
    CorrXYFreeT=squeeze(strain.s(GMTind,GMPos.Free,k));

    CorrXYImpB=squeeze(strain.s(GMBind,GMPos.Imp,k));
    CorrXYMidB=squeeze(strain.s(GMBind,GMPos.Mid,k));
    CorrXYFreeB=squeeze(strain.s(GMBind,GMPos.Free,k));

    %% Top Edges
    %%  Shear strain at the impact edge
    subplot(2,3,1)
    plot(FETcoord*10^3,FEXYImpT,'k')
    hold on
    plot(GMTcoord*10^3,RawXYImpT,'b')
    plot(GMTcoord*10^3,CorrXYImpT,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('strain_{xy}')
    title(strcat('Top 4P from Impact ',Fstring))
    legend('FE','RawGM','Corrected','location','southwest')
    ylim(TXYLim)

    %% Shear strain in middle
    subplot(2,3,2)
    plot(FETcoord*10^3,FEXYMidT,'k')
    hold on
    plot(GMTcoord*10^3,RawXYMidT,'b')
    plot(GMTcoord*10^3,CorrXYMidT,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('strain_{xy}')
    title(strcat('Top Middle ',Fstring))
    ylim(TXYLim)

    %% Shear strain Free Surface
    subplot(2,3,3)
    plot(FETcoord*10^3,FEXYFreeT,'k')
    hold on
    plot(GMTcoord*10^3,RawXYFreeT,'b')
    plot(GMTcoord*10^3,CorrXYFreeT,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('strain_{xy}')
    title(strcat('Top 4P from Free Edge ',Fstring))
    ylim(TXYLim)

    %%Bottom Edges
    %%  Shear strain at the impact edge
    subplot(2,3,4)
    plot(FEBcoord*10^3,FEXYImpB,'k')
    hold on
    plot(GMBcoord*10^3,RawXYImpB,'b')
    plot(GMBcoord*10^3,CorrXYImpB,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('strain_{xy}')
    title(strcat('Bottom 4P from Impact ',Fstring))
    legend('FE','RawGM','Corrected','location','southwest')
    ylim(BXYLim)

    %% Shear strain in middle
    subplot(2,3,5)
    plot(FEBcoord*10^3,FEXYMidB,'k')
    hold on
    plot(GMBcoord*10^3,RawXYMidB,'b')
    plot(GMBcoord*10^3,CorrXYMidB,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('strain_{xy}')
    title(strcat('Bottom Middle ',Fstring))
    ylim(BXYLim)

    %% Shear strain Free Surface
    subplot(2,3,6)
    plot(FEBcoord*10^3,FEXYFreeB,'k')
    hold on
    plot(GMBcoord*10^3,RawXYFreeB,'b')
    plot(GMBcoord*10^3,CorrXYFreeB,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('strain_{xy}')
    title(strcat('Bottom 4P from Free Edge ',Fstring))
    ylim(BXYLim)

    %% SaveFiles
    FigName=strcat(ZoomFigDir,TestDeg,'_1DStrainXYZoom_',Fnum,'.fig');
    saveas(gcf,FigName);
    PngName=strcat(ZoomPngDir,TestDeg,'_1DStrainXYZoom_',Fnum,'.png');
    saveas(gcf,PngName);
    
end

%% Create Directory
ZoomPngDir=strcat(PngDir,'/ZoomXXstrain/');
ZoomFigDir=strcat(FigDir,'/ZoomXXstrain/');
mkdir(ZoomPngDir);
mkdir(ZoomFigDir);
%% Plot the zoomed Shear strains

%Determine limits
TXXLim=1.1*[min(FE.strain.x(FETind,:,:),[],'all'),...
    max(FE.strain.x(FETind,:,:),[],'all')];
BXXLim=1.1*[min(FE.strain.x(FEBind,:,:),[],'all'),...
    max(FE.strain.x(FEBind,:,:),[],'all')];

%Initialize figure
figure('Units','normalized','outerposition',[0,0,1,1])

for k=1:length(time.vec)
    Fnum=num2str(k);
    Fstring=strcat('Frame~',Fnum);

        %% Extract vector
    %FE Model Displacement
    FEXXImpT=squeeze(FE.strain.x(FETind,FEPos.Imp,k));
    FEXXMidT=squeeze(FE.strain.x(FETind,FEPos.Mid,k));
    FEXXFreeT=squeeze(FE.strain.x(FETind,FEPos.Free,k));
    
    FEXXImpB=squeeze(FE.strain.x(FEBind,FEPos.Imp,k));
    FEXXMidB=squeeze(FE.strain.x(FEBind,FEPos.Mid,k));
    FEXXFreeB=squeeze(FE.strain.x(FEBind,FEPos.Free,k));
    
    %Uncorrected Grid Method Displacement
    RawXXImpT=squeeze(strainRaw.x(GMTind,GMPos.Imp,k));
    RawXXMidT=squeeze(strainRaw.x(GMTind,GMPos.Mid,k));
    RawXXFreeT=squeeze(strainRaw.x(GMTind,GMPos.Free,k));

    RawXXImpB=squeeze(strainRaw.x(GMBind,GMPos.Imp,k));
    RawXXMidB=squeeze(strainRaw.x(GMBind,GMPos.Mid,k));
    RawXXFreeB=squeeze(strainRaw.x(GMBind,GMPos.Free,k));

    %Corrected Grid Method Dipslacement
    CorrXXImpT=squeeze(strain.x(GMTind,GMPos.Imp,k));
    CorrXXMidT=squeeze(strain.x(GMTind,GMPos.Mid,k));
    CorrXXFreeT=squeeze(strain.x(GMTind,GMPos.Free,k));

    CorrXXImpB=squeeze(strain.x(GMBind,GMPos.Imp,k));
    CorrXXMidB=squeeze(strain.x(GMBind,GMPos.Mid,k));
    CorrXXFreeB=squeeze(strain.x(GMBind,GMPos.Free,k));

    %% Top Edges
    %%  XX strain at the impact edge
    subplot(2,3,1)
    plot(FETcoord*10^3,FEXXImpT,'k')
    hold on
    plot(GMTcoord*10^3,RawXXImpT,'b')
    plot(GMTcoord*10^3,CorrXXImpT,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('strain_{xx}')
    title(strcat('Top 4P from Impact ',Fstring))
    legend('FE','RawGM','Corrected','location','southwest')
    ylim(TXXLim)

    %% XX strain in middle
    subplot(2,3,2)
    plot(FETcoord*10^3,FEXXMidT,'k')
    hold on
    plot(GMTcoord*10^3,RawXXMidT,'b')
    plot(GMTcoord*10^3,CorrXXMidT,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('strain_{xx}')
    title(strcat('Top Middle ',Fstring))
    ylim(TXXLim)

    %% XX strain Free Surface
    subplot(2,3,3)
    plot(FETcoord*10^3,FEXXFreeT,'k')
    hold on
    plot(GMTcoord*10^3,RawXXFreeT,'b')
    plot(GMTcoord*10^3,CorrXXFreeT,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('strain_{xx}')
    title(strcat('Top 4P from Free Edge ',Fstring))
    ylim(TXXLim)

    %%Bottom Edges
    %%  XX strain at the impact edge
    subplot(2,3,4)
    plot(FEBcoord*10^3,FEXXImpB,'k')
    hold on
    plot(GMBcoord*10^3,RawXXImpB,'b')
    plot(GMBcoord*10^3,CorrXXImpB,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('strain_{xx}')
    title(strcat('Bottom 4P from Impact ',Fstring))
    legend('FE','RawGM','Corrected','location','southwest')
    ylim(BXXLim)

    %% XX strain in middle
    subplot(2,3,5)
    plot(FEBcoord*10^3,FEXXMidB,'k')
    hold on
    plot(GMBcoord*10^3,RawXXMidB,'b')
    plot(GMBcoord*10^3,CorrXXMidB,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('strain_{xx}')
    title(strcat('Bottom Middle ',Fstring))
    ylim(BXXLim)

    %% XX strain Free Surface
    subplot(2,3,6)
    plot(FEBcoord*10^3,FEXXFreeB,'k')
    hold on
    plot(GMBcoord*10^3,RawXXFreeB,'b')
    plot(GMBcoord*10^3,CorrXXFreeB,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('strain_{xx}')
    title(strcat('Bottom 4P from Free Edge ',Fstring))
    ylim(BXXLim)

    %% SaveFiles
    FigName=strcat(ZoomFigDir,TestDeg,'_1DStrainXXZoom_',Fnum,'.fig');
    saveas(gcf,FigName);
    PngName=strcat(ZoomPngDir,TestDeg,'_1DStrainXXZoom_',Fnum,'.png');
    saveas(gcf,PngName);
    
end

%% Add reference parameters 
[RefPar.file,RefPar.path]=uigetfile('*.mat',...
    'Choose Mat file containing reference constitutive pamrameters');
load(strcat(RefPar.path,'/',RefPar.file),'MatProps');

%% Calculate out of plane strains
FE.StressModel=func_ViscoConstitutiveV6(FE.strain.x,FE.strain.y,...
    FE.strain.s,time.vec,MatProps,0,0,0);
StressModel=func_ViscoConstitutiveV6(strain.x,strain.y,strain.s, ...
    time.vec,MatProps,0,0,0);
RawModel=func_ViscoConstitutiveV6(strainRaw.x,strainRaw.y,strainRaw.s, ...
    time.vec,MatProps,0,0,0);

%% Plot StrainZZ
%% ZZ Strain directory
PngDirZZ=strcat(PngDir,'/ZZstrain/');
FigDirZZ=strcat(FigDir,'/ZZstrain/');
mkdir(PngDirZZ);
mkdir(FigDirZZ);
%% ZZ strain
figure('Units','normalized','OuterPosition',[0,0,1,1])
timems=time.vec*10^6;
Ymm=pos.y*10^3;
YmmFE=FE.pos.y*10^3;
% Set plot limits
ZZstrainLim=1.1*[min(FE.StressModel.zzstrain,[],'all'),...
    max(FE.StressModel.zzstrain,[],'all')];

for k=1:length(time.vec)
    Fnum=num2str(k);
    Fstring=strcat('Frame~',Fnum);
    
    %% Extract vector
    %FE Model Displacement
    FEZZImp=squeeze(FE.StressModel.zzstrain(:,FEPos.Imp,k));
    FEZZMid=squeeze(FE.StressModel.zzstrain(:,FEPos.Mid,k));
    FEZZFree=squeeze(FE.StressModel.zzstrain(:,FEPos.Free,k));
    
    %Uncorrected Grid Method Displacement
    RawZZImp=squeeze(RawModel.zzstrain(:,GMPos.Imp,k));
    RawZZMid=squeeze(RawModel.zzstrain(:,GMPos.Mid,k));
    RawZZFree=squeeze(RawModel.zzstrain(:,GMPos.Free,k));


    %Corrected Grid Method Dipslacement
    CorrZZImp=squeeze(StressModel.zzstrain(:,GMPos.Imp,k));
    CorrZZMid=squeeze(StressModel.zzstrain(:,GMPos.Mid,k));
    CorrZZFree=squeeze(StressModel.zzstrain(:,GMPos.Free,k));

    %%  ZZ strain at the impact edge
    subplot(3,1,1)
    plot(YmmFE,FEZZImp,'k')
    hold on
    plot(Ymm,RawZZImp,'b')
    plot(Ymm,CorrZZImp,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('strain_{zz}')
    title(strcat('4P from Impact ',Fstring))
    legend('FE','RawGM','Corrected','location','southwest')
    ylim(ZZstrainLim)

    %% Y Displacement in middle
    subplot(3,1,2)
    plot(YmmFE,FEZZMid,'k')
    hold on
    plot(Ymm,RawZZMid,'b')
    plot(Ymm,CorrZZMid,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('strain_{zz}')
    title(strcat('Middle ',Fstring))
    ylim(ZZstrainLim)

    %% Y displacement Free Surface
    subplot(3,1,3)
    plot(YmmFE,FEZZFree,'k')
    hold on
    plot(Ymm,RawZZFree,'b')
    plot(Ymm,CorrZZFree,'r')
    hold off
    xlabel('Y Coordinate (mm)')
    ylabel('strain_{zz}')
    title(strcat('4P from Free Edge ',Fstring))
    ylim(ZZstrainLim)

    %% SaveFiles
    FigName=strcat(FigDirZZ,TestDeg,'_1DstrainZ_',Fnum,'.fig');
    saveas(gcf,FigName);
    PngName=strcat(PngDirZZ,TestDeg,'_1DstrainZ_',Fnum,'.png');
    saveas(gcf,PngName);

end

%% Get X indicies
    GMPos.NumPoints=size(disp.x,2);
    GMPos.Free=20; %4 pitches from free edge
    GMPos.Mid=round(GMPos.NumPoints/2);
    GMPos.Imp=GMPos.NumPoints-20; %4 pitches from impact edge
    % Get cooridinates of position indexes
    GMPos.Xfree=pos.x(GMPos.Free);
    GMPos.Xmid=pos.x(GMPos.Mid);
    GMPos.Ximp=pos.x(GMPos.Imp);

    %Set Position indexes and coordintes for finite element data
    FEPos.NumPoints=size(FE.disp.x,2);
    FEPos.GridRat=FEPos.NumPoints/GMPos.NumPoints;
    FEPos.Free=round(GMPos.Free*FEPos.GridRat); %4 pitches from free edge
    FEPos.Mid=round(FEPos.NumPoints/2);
    FEPos.Imp=round(GMPos.Imp*FEPos.GridRat); %4 pitches from impact edge
    
    FEPos.Xfree=FE.pos.x(FEPos.Free);
    FEPos.Xmid=FE.pos.x(FEPos.Mid);
    FEPos.Ximp=FE.pos.x(FEPos.Imp);
    
  %Print Plot selections for evaluation purposes
  FEPos.freeString=num2str(FEPos.Xfree);
  GMPos.freeString=num2str(GMPos.Xfree);
  fprintf(strcat('Free Coordinate FE:',FEPos.freeString,' Grid: ',...
      GMPos.freeString,'\n'));
  
  FEPos.midString=num2str(FEPos.Xmid);
  GMPos.midString=num2str(GMPos.Xmid);
  fprintf(strcat('Mid Coordinate FE:',FEPos.midString,' Grid: ',...
      GMPos.midString,'\n'));

  FEPos.impString=num2str(FEPos.Ximp);
  GMPos.impString=num2str(GMPos.Ximp);
  fprintf(strcat('Impact Coordinate FE:',FEPos.impString,' Grid: ',...
      GMPos.impString,'\n'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Investigate accelerations 


%% Find Y indexes and coordinates
    GMPos.NumPointsY=size(disp.x,1);
    GMPos.Bot=20; %4 pitches from free edge
    GMPos.MidY=round(GMPos.NumPointsY/2);
    GMPos.Top=GMPos.NumPointsY-20; %4 pitches from impact edge
    
    % Get cooridinates of position indexes
    GMPos.YBot=pos.y(GMPos.Bot);
    GMPos.Ymid=pos.y(GMPos.MidY);
    GMPos.YTop=pos.y(GMPos.Top);

    %Set Position indexes and coordintes for finite element data
    FEPos.NumPointsY=size(FE.disp.x,1);
    FEPos.GridRatY=FEPos.NumPointsY/GMPos.NumPointsY;
    FEPos.Bot=round(GMPos.Bot*FEPos.GridRatY); %4 pitches from free edge
    FEPos.MidY=round(FEPos.NumPointsY/2);
    FEPos.Top=round(GMPos.Top*FEPos.GridRatY); %4 pitches from impact edge
    
    FEPos.YBot=FE.pos.y(FEPos.Bot);
    FEPos.Ymid=FE.pos.y(FEPos.MidY);
    FEPos.YTop=FE.pos.y(FEPos.Top);
    
  %Print Plot selections for evaluation purposes
  FEPos.BotString=num2str(FEPos.YBot);
  GMPos.BotString=num2str(GMPos.YBot);
  fprintf(strcat('Bottom Coordinate FE:',FEPos.BotString,' Grid: ',...
      GMPos.BotString,'\n'));
  
  FEPos.YmidString=num2str(FEPos.Ymid);
  GMPos.YmidString=num2str(GMPos.Ymid);
  fprintf(strcat('Mid Coordinate FE:',FEPos.YmidString,' Grid: ',...
      GMPos.YmidString,'\n'));

  FEPos.TopString=num2str(FEPos.YTop);
  GMPos.TopString=num2str(GMPos.YTop);
  fprintf(strcat('Top Coordinate FE:',FEPos.TopString,' Grid: ',...
      GMPos.TopString,'\n'));

%% Plot X accelerations

%create directory
PngDirXa=strcat(PngDir,'/AccelX/');
FigDirXa=strcat(FigDir,'/AccelX/');
mkdir(PngDirXa);
mkdir(FigDirXa);

figure('units','normalized','outerposition',[0,0,1,1])

timems=time.vec*10^6;
Xmm=pos.x*10^3;
XmmFE=FE.pos.x*10^3;
XALim=[min(FE.accel.x,[],'all'),max(FE.accel.x,[],'all')];

for k=1:length(time.vec)
    Fnum=num2str(k);
    Fstring=strcat('Frame~',Fnum);    
    %% Extract vector
    %FE Model Displacement
    FEAxBot=squeeze(FE.accel.x(FEPos.Bot,:,k));
    FEAxMidY=squeeze(FE.accel.x(FEPos.MidY,:,k));
    FEAxTop=squeeze(FE.accel.x(FEPos.Top,:,k));
    
    %Uncorrected Grid Method Displacement
    RawAxBot=squeeze(RawAccel.x(GMPos.Bot,:,k));
    RawAxMidY=squeeze(RawAccel.x(GMPos.MidY,:,k));
    RawAxTop=squeeze(RawAccel.x(GMPos.Top,:,k));


    %Corrected Grid Method Dipslacement
    CorrAxBot=squeeze(accel.x(GMPos.Bot,:,k));
    CorrAxMidY=squeeze(accel.x(GMPos.MidY,:,k));
    CorrAxTop=squeeze(accel.x(GMPos.Top,:,k));

    %%  x Acceleration at the Bottom edge
    subplot(3,1,1)
    plot(XmmFE,FEAxBot,'k')
    hold on
    plot(Xmm,RawAxBot,'b')
    plot(Xmm,CorrAxBot,'r')
    hold off
    xlabel('X Coordinate (mm)')
    ylabel('a_x (m/s^2)')
    title(strcat('4P from Bottom Edge',Fstring))
    legend('FE','RawGM','Corrected','location','southwest')
    ylim(XALim)

    %% X Acceleration in middle
    subplot(3,1,2)
    plot(XmmFE,FEAxMidY,'k')
    hold on
    plot(Xmm,RawAxMidY,'b')
    plot(Xmm,CorrAxMidY,'r')
    hold off
    xlabel('X Coordinate (mm)')
    ylabel('a_x (m/s^2)')
    title(strcat('Middle ',Fstring))
    ylim(XALim)

    %% X Acceleration Top
    subplot(3,1,3)
    plot(XmmFE,FEAxTop,'k')
    hold on
    plot(Xmm,RawAxTop,'b')
    plot(Xmm,CorrAxTop,'r')
    hold off
    xlabel('X Coordinate (mm)')
    ylabel('a_x (m/s^2)')
    title(strcat('4P from Top Edge ',Fstring))
    ylim(XALim)

    %% SaveFiles
    FigName=strcat(FigDirXa,TestDeg,'_1DxA_',Fnum,'.fig');
    saveas(gcf,FigName);
    PngName=strcat(PngDirXa,TestDeg,'_1DxA_',Fnum,'.png');
    saveas(gcf,PngName);
end

%% Calculate surface average accelerations
FE.Ay_surfx=zeros(length(FE.pos.x),length(time.vec));
FE.Ax_surfx=zeros(length(FE.pos.x),length(time.vec));
RawAccel.Ax_surfx=zeros(length(pos.x),length(time.vec));
RawAccel.Ay_surfx=zeros(length(pos.x),length(time.vec));
accel.Ax_surfx=zeros(length(pos.x),length(time.vec));
accel.Ay_surfx=zeros(length(pos.x),length(time.vec));

Ay_avx=mean(accel.y);
Ax_avx=mean(accel.x);
RawAccel.Ay_avx=mean(RawAccel.y);
RawAccel.Ax_avx=mean(RawAccel.x);
FE.Ay_avx=mean(FE.accel.y);
FE.Ax_avx=mean(FE.accel.x);

for m=1:length(pos.x)
    for n=1:length(time.vec)
        accel.Ay_surfx(m,n)=mean(Ay_avx(1,1:m,n));
        accel.Ax_surfx(m,n)=mean(Ax_avx(1,1:m,n));
        RawAccel.Ax_surfx(m,n)=mean(RawAccel.Ax_avx(1,1:m,n));
        RawAccel.Ay_surfx(m,n)=mean(RawAccel.Ay_avx(1,1:m,n));
     end
end

for m=1:length(FE.pos.x)
    for n=1:length(time.vec)
        FE.Ay_surfx(m,n)=mean(FE.Ay_avx(1,1:m,n));
        FE.Ax_surfx(m,n)=mean(FE.Ax_avx(1,1:m,n));
     end
end

%% Plot X surface accelerations

%create directory
PngDirSurf=strcat(PngDir,'/AccelSurf/');
FigDirSurf=strcat(FigDir,'/AccelSurf/');
mkdir(PngDirSurf);
mkdir(FigDirSurf);

figure('units','normalized','outerposition',[0,0,1,1])

timems=time.vec*10^6;
Xmm=pos.x*10^3;
XmmFE=FE.pos.x*10^3;
XsurfLim=[min(FE.Ax_surfx,[],'all'),max(FE.Ax_surfx,[],'all')];
YsurfLim=[min(FE.Ay_surfx,[],'all'),max(FE.Ay_surfx,[],'all')];

for k=1:length(time.vec)
    Fnum=num2str(k);
    Fstring=strcat('Frame~',Fnum);    
    %% Extract vector
    %FE Model Displacement
    FESurfx=squeeze(FE.Ax_surfx(:,k));
    FESurfy=squeeze(FE.Ay_surfx(:,k));
    %Uncorrected Grid Method Displacement
    RawSurfx=squeeze(RawAccel.Ax_surfx(:,k));
    RawSurfy=squeeze(RawAccel.Ay_surfx(:,k));
    %Corrected Grid Method Dipslacement
    CorrSurfx=squeeze(accel.Ax_surfx(:,k));
    CorrSurfy=squeeze(accel.Ay_surfx(:,k));

    %%  x surface average accelerations
    subplot(2,1,1)
    plot(XmmFE,FESurfx,'k')
    hold on
    plot(Xmm,RawSurfx,'b')
    plot(Xmm,CorrSurfx,'r')
    hold off
    xlabel('X Coordinate (mm)')
    ylabel('a^s_x (m/s^2)')
    title(strcat('X surface accelerations',Fstring))
    legend('FE','RawGM','Corrected','location','southwest')
    ylim(XsurfLim)

    subplot(2,1,2)
    plot(XmmFE,FESurfy,'k')
    hold on
    plot(Xmm,RawSurfy,'b')
    plot(Xmm,CorrSurfy,'r')
    hold off
    xlabel('X Coordinate (mm)')
    ylabel('a^s_y (m/s^2)')
    title(strcat('Y surface accelerations',Fstring))
    %legend('FE','RawGM','Corrected','location','southwest')
    ylim(YsurfLim)

    %% SaveFiles
    FigName=strcat(FigDirSurf,TestDeg,'_1DSurf_',Fnum,'.fig');
    saveas(gcf,FigName);
    PngName=strcat(PngDirSurf,TestDeg,'_1DSurf_',Fnum,'.png');
    saveas(gcf,PngName);
end
%% Calculate FE accelerations from displacements in the same manner as for
    %The grid method
 diffOpts.method='gradient';
 fprintf( ...
 'Calculating accelerations in the same manner as for grid method data \n')
 [~,FE.aCalc] = func_calcFFAccelFromDisp(FE.time,FE.disp,diffOpts.method);


%% Calculate stress guage stresses
rho=MatProps.rho;
FE.SG=func_Full_SG(FE.accel,FE.pos.x,FE.time,rho);
FE.SGcalc=func_Full_SG(FE.aCalc,FE.pos.x,FE.time,rho);
SG=func_Full_SG(accel,pos.x,time,rho);
SGRaw=func_Full_SG(RawAccel,pos.x,time,rho);

%% Plot stress gauge stresses against FE stresses and constitutive model
%create directory
PngDirStress=strcat(PngDir,'/AvStress/');
FigDirStress=strcat(FigDir,'/AvStress/');
mkdir(PngDirStress);
mkdir(FigDirStress);


%Calculate average stresses with FE
FE.stress.xAvg=squeeze(mean(FE.stress.x));
FE.stress.sAvg=squeeze(mean(FE.stress.s));



figure('units','normalized','outerposition',[0,0,1,1])

timems=time.vec*10^6;
Xmm=pos.x*10^3;
XmmFE=FE.pos.x*10^3;
XStressLim=[min(FE.stress.xAvg,[],'all'),max(FE.stress.xAvg,[],'all')];
SStressLim=[min(FE.stress.sAvg,[],'all'),max(FE.stress.sAvg,[],'all')];

for k=1:length(time.vec)
    Fnum=num2str(k);
    Fstring=strcat('Frame~',Fnum);    
    %% Extract vector
    %FE average stresses
    FEStressX=squeeze(FE.stress.xAvg(:,k));
    FEStressS=squeeze(FE.stress.sAvg(:,k));
    
    %FE SG
    FEsgX=squeeze(FE.SG.x(:,k));
    FEsgS=squeeze(FE.SG.s(:,k));
    FEsgDX=squeeze(FE.SGcalc.x(:,k));
    FEsgDS=squeeze(FE.SGcalc.s(:,k));
  
    %GM stress gauge
    RawSGx=squeeze(SGRaw.x(:,k));
    RawSGs=squeeze(SGRaw.s(:,k));
    CorrSGx=squeeze(SG.x(:,k));
    CorrSGs=squeeze(SG.s(:,k));

    %Constitutive Model 
    RawConstX=squeeze(RawModel.Avxx(:,k));
    RawConstS=squeeze(RawModel.Avxy(:,k));
    CorrConstX=squeeze(StressModel.Avxx(:,k));
    CorrConstS=squeeze(StressModel.Avxy(:,k));

     %%  X average stresses
    subplot(2,1,1)
    plot(XmmFE,FEStressX,'k')
    hold on
    plot(XmmFE,FEsgX,'b')
    plot(XmmFE,FEsgDX,'r')
    plot(Xmm,RawSGx,'b:')
    plot(Xmm,CorrSGx,'r:')
    plot(Xmm,RawConstX,'b--')
    plot(Xmm,CorrConstX,'r--')
    hold off
    xlabel('X Coordinate (mm)')
    ylabel('\sigma_{xx} (Pa)')
    title(strcat('X surface accelerations',Fstring))
    legend('FE stress','FEsg','FEsgDisp','RawSG','CorrectedSG', ...
        'RawModel','Corrected Model','location','southwest')
    ylim(XStressLim)

    %% average shear stresses
    subplot(2,1,2)
    plot(XmmFE,FEStressS,'k')
    hold on
    plot(XmmFE,FEsgS,'b')
    plot(XmmFE,FEsgDS,'r')
    plot(Xmm,RawSGs,'b:')
    plot(Xmm,CorrSGs,'r:')
    plot(Xmm,RawConstS,'b--')
    plot(Xmm,CorrConstS,'r--')
    hold off
    xlabel('X Coordinate (mm)')
    ylabel('\sigma_{xy} (Pa)')
    title(strcat('X surface accelerations',Fstring))
    legend('FE stress','FEsg','FEsgDisp','RawSG','CorrectedSG', ...
        'RawModel','Corrected Model','location','southwest')
    ylim(XStressLim)
    %% SaveFiles
    FigName=strcat(FigDirStress,TestDeg,'_1DStress_',Fnum,'.fig');
    saveas(gcf,FigName);
    PngName=strcat(PngDirStress,TestDeg,'_1DStress_',Fnum,'.png');
    saveas(gcf,PngName);
end


