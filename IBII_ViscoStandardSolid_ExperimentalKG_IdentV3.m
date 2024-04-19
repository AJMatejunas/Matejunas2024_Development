%This script is written to extract Generalized Maxwell Parameters for the
%standard solid model from experimental images.

%Author: Andrew Matejunas (andrew.matejunas@gmail.com)

%Version History/Change log:
    %2023-06-22: Initial version
    %2023-08-03: V2 added capability to extract displacements from DIC data
    %2023-08-11: Increased upper bound on Tau to 1ms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialize
clear variables; close all; clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('-----------------------------------------------------------------\n')
fprintf('INITIALIZING \n')
fprintf('------------------------------------------------------------------ \n')
fprintf('Defining known properties and processing parameters \n')


%% Set default fonts for figures to arial for consitency with inkscape
figFont='Arial';
figFsize=10;

labelProps={'Interpreter','latex','FontName',figFont,'FontSize',figFsize};
MajorProps={'Interpreter','latex','FontName',figFont,'FontSize',14};
axprops={'Ydir','Normal','FontSize',10};

%% Set figure size variables
%Set figure widths (Note these come from strain author guidelines
SCW=8.3; %Single Column in cm
DCW=17.5; %Double Column Width in cm
CW15=12.5; %1.5 column width

AR11=1; %standard aspect ratio of 1 to 1
AR43=3/4; %4:3 aspect ration for line plots
AR23=2/3;
AR13=1/3;
AR32=3/2;
%single figure (defaults to 1.5 column width) 
    %ideal for single 1-D plots

SingFig={'Units','Centimeters','InnerPosition',[1,1,CW15,AR43*CW15]};

%Double Colum figure for Double column width
DcSrFig={'Units','Centimeters','InnerPosition',[1,1,DCW,AR11*SCW]};
%Double column double row
DcDrFig={'Units','Centimeters','InnerPosition',[1,1,DCW,AR11*DCW]};
%Triple column
TcSrFig={'Units','Centimeters','InnerPosition',[1,1,DCW,AR13*DCW]};
TcDrFig={'Units','Centimeters','InnerPosition',[1,1,DCW,AR23*DCW]};
DcTrFig={'Units','Centimeters','InnerPosition',[1,1,DCW,AR11*DCW]};
T3B2Fig={'Units','Centimeters','InnerPosition',[1,1,DCW,AR32*DCW]};
TcTrFig={'Units','Centimeters','InnerPosition',[1,1,DCW,3*AR13*DCW]};
%Full Page Portrait Figure
FPpFig={'Units','Centimeters','InnerPosition',[1,1,DCW,20]};
%Full Page Landscape
FPlFig={'Units','Centimeters','InnerPosition',[0.5,0.5,23,DCW]};

%% Presentation figures
PreFig={'units','inches','InnerPosition',[0.5,0.5,8,3]};
PreFig2={'units','inches','InnerPosition',[0.5,0.5,4,5]};
PreLabel={'Interpreter','latex','FontName',figFont,'FontSize',18};
PreMajor={'Interpreter','latex','FontName',figFont,'FontSize',20};
PreAx={'Ydir','Normal','FontSize',16};

%% Choose Folder to save data
ExpSavePath=uigetdir('','Choose Folder to Save Experimental Results');

%% Input test designation
ExpDesig=char(cell2mat(inputdlg('Input Experiment Designation \n')));

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

%% Load Processing Parameter Data Structures
fprintf('Loading processing parameters file.\n')
[initFile,initPath,~] = uigetfile('*.mat','Locate processing parameter file');

% Load the processing parameters from file
load([initPath,initFile])

% Store the FE valid mode parameter
globalOpts.FEValidMode = 0;
globalOpts.calcKinFieldsFromDisp = true;

%% LOAD RAW DATA FILES: Raw .tiff images
switch globalOpts.dispSource
    case 'GM'
    fprintf('Loading reference image from the selected test data folder.\n')

    [imageFile,imagePath] = uigetfile({'*.*','All Files'},...
        'Select the first image in the sequence');

    Static.Img1=strcat(imagePath,'/',imageFile);
    [Static.file2,Static.path2]=uigetfile({'*.*','All Files'}, ...
        'Select the second image in the sequence');
    Static.Img2=strcat(Static.path2,'/',Static.file2);
end


%% Set Guess Parameters
prompt={'Literature Elastic Modulus'; 'QS Poissons ratio';
    'Constant Nu (Yes or No)'};
defin={'5.19e9','0.36','No'};
dims=1;
dlgtitle='Literature properties for initial guess calculation';
dlgOpts.Resize='on';
dlgOpts.WindowStyle='normal';
tempProps=inputdlg(prompt,dlgtitle,dims,defin,dlgOpts);

material.E0_Lit=str2double(tempProps(1));
material.Nu_inf=str2double(tempProps(2));


material.G0_Lit=material.E0_Lit/(2*(1+material.Nu_inf));
%effective shear modulus (literature high rate -  measured QS)
material.G1_eff=material.G0_Lit-material.Ginf;

material.K0_Lit=material.E0_Lit/(3*(1-2*material.Nu_inf));
%effective Bulk modulus (literature high rate -  measured QS)
material.K1_eff=material.K0_Lit-material.Kinf;

switch cell2mat(tempProps(3))
    case 'Yes'
    material.nu=material.Nu_inf;
    case 'No'
    material.nu=0;
end


% Create reasonable upper bounds
%bounds on the moduli initial guess are effective literature parameters +/- 50%
ub=[material.K1_eff*(1.5);
    material.G1_eff*(1.5);... G1
    100E-3]; %tau_1
ubG=ub(2:3);
ubK=[ub(1);ub(3)];

% Create reasonable lower bounds
lb=[material.K1_eff*(0.5);
    material.G1_eff*(0.5);... G1
    1E-7]; %tau_1
lbG=lb(2:3);
lbK=[lb(1);lb(3)];


%% Create Matrix of initial guesses
    %Eint=2.5E9;
    %Choose a random 3x3 matrix of initial guress multipliers. Multipliers
        %will be uniformly distributed random numbers in the interval [0,1]
 
    intMult=rand(3,3);
    intMultK=rand(1,3);
  
    %Create 3X3 matrix of initial guess parameters
     intMatrix=lb+(ub-lb).*intMult;
     
     intMatrixG=intMatrix(2:3,:);
    

     intMatrixK=[intMatrix(1,:);intMatrix(3,:)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Noise Floor Evaluation
NoiseQuest='Evaluate Kinematic Noise Floors?';
NoiseChoice=questdlg(NoiseQuest);

switch NoiseChoice
    case 'Yes'
     fprintf('--------------------------------------------------------\n')
     fprintf('EVALUATING KINEMATIC NOISE FLOORS\n')
     fprintf('----------------------------------------------------------\n ')
     
    staticQuest='Use impact images for underformed?';
    Static.Choice=questdlg(staticQuest);
    switch Static.Choice
        case 'Yes'
     Static.Path=imagePath;
     Static.File=imageFile;
        case 'No'
       [Static.File,Static.Path]=uigetfile({'*.*','All Files'}, ...
            'Select the first image in the undeformed reference image');
    end    
       
    %Determine bit depth of the image
    Static.Img1=strcat(Static.Path,'/',Static.File);
    Static.info=imfinfo(Static.Img1);
    Static.BitD=Static.info.BitDepth;

     %% Calculate Noise Floor
    
     NF=func_CalcExpNoiseFloor(Static.Path,Static.File,Static.BitD,grid, ...
        specimen,time,material,Shear.extrapOpts,Shear.smoothOpts, ...
        Bulk.extrapOpts,Bulk.smoothOpts,diffOpts,globalOpts,gridMethodOpts);


    case 'No'
   fprintf('NOT EVALUATING NOISE FLOORS \n')
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Kinematic Field Extraction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ensure No noise is added
imageNoise.addNoise=0;
%% Extract  Displacements
switch globalOpts.dispSource
    case 'GM'
    fprintf('Grid Method Displacement Extraction \n')
    fprintf('\n--------------------------------------------------------------\n')

    %                 fprintf('GRID METHOD PROCESSING\n')
    %                 fprintf('--------------------------------------------------------------\n')
    %
    %                 %--------------------------------------------------------------------------
    %                 % GRID IMAGE PROCESSING

    % fprintf('Processing images using the grid method toolbox.\n')
    % Process the image squence with the grid method toolbox
    [grid,pos,disp] = func_gridMethodImageProcessing(imagePath,...
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
    % time.numFrames = size(disp.x,3);
    %time.vec = 0:time.step:(size(disp.x,3)-1)*time.step;

    %Create arrays of x and y vectors
    case'MatchID'
    if globalOpts.ExtractDICStrain==false
        [disp,pos,DIC,specimen,~] = func_ExtractMatchIDFields(DIC, ...
            specimen,material,globalOpts,DIC.ConvUnit);
    elseif globalOpts.ExtractDICStrain==true
        [disp,pos,DIC,specimen,strain] = func_ExtractMatchIDFields(DIC, ...
            specimen,material,globalOpts,DIC.ConvUnit);

    end
end

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

%% Save Raw Displacements
fprintf('Saving Raw Displacement fields \n')
RawName=strcat(ExpSavePath,'/',ExpDesig,'_RawDisplacementData.mat');
save(RawName)



%% SHEAR IDENTIFICATION
fprintf('---------------------------------------------------------------- \n')
fprintf('SHEAR KINEMATIC CALCULATIONS \n')
fprintf('---------------------------------------------------------------- \n')

% Extrapolate the data to account for the missing pitch on the edges

% Crop and Extrapolate displacement fields


fprintf('Croping and Extrapolating Displacement Fields \n')
[Shear.disp.x,Shear.disp.y,Shear.disp.rX,Shear.disp.rY]=...
    func_cropAndExtrapDispFields_v4(pos,disp.x,disp.y, ...
    Shear.extrapOpts.disp,false);
pos.Raw=pos;
switch globalOpts.dispSource
    case 'MatchID'
        fprintf('Extrapolating Pos to be full ROI \n')
        pos=func_ExtrapolatePos(pos,DIC.edgeCut);
end

switch globalOpts.dispSource
    case 'MatchID'
       
        %Steps to crop
        cST=Shear.extrapOpts.disp.cropPx1st;
        %Steps to extrapolate
        eST=Shear.extrapOpts.disp.extrapPx1st;
        
        %Data Points missing from edge
        missST=DIC.edgeCut;

        %Determine range to crop back to full specimen size
        Shear.disp.rX=(1+eST-cST-missST):(size(Shear.disp.x,2)...
            -eST+cST+missST);
        Shear.disp.rY=(1+eST-cST-missST):(size(Shear.disp.x,1)...
            -eST+cST+missST);
end

X_vec=pos.x;
Y_vec=pos.y;

fprintf('Obtaining and setting the free edge location.\n')
switch globalOpts.hardCodeFreeEdge
    case false
    [freeEdge,specimen,Shear.disp] =...
    func_getFreeEdge(globalOpts.hardCodeFreeEdge,...
    imagePath,imageFile,specimen,Shear.disp);
    case true
        [freeEdge,specimen,Shear.disp] =...
    func_CheckFreeEdge(specimen,Shear.disp);
end

%%
if globalOpts.ExtractDICStrain==false
    fprintf('Calculating Strains for Shear Identification \n')
    [Shear.strain,Shear.disp]=func_smoothCalcStrain_v4(pos,time,Shear.disp,...
        Shear.smoothOpts.strain,Shear.extrapOpts.strain,false);
    Shear.strain.sAv=squeeze(mean(Shear.strain.s));
end
%--------------------------------------------------------------------------
% Smooth and Calculate Acceleration
if isfield(Shear.smoothOpts.accel,'spatialSmooth')
    if Shear.smoothOpts.accel.spatialSmooth==false
        fprintf('Calculating acceleration from the unsmoothed displacement fields.\n')
        [Shear.accel,~,Shear.disp] = func_smoothCalcAccel_v4(pos,time,Shear.disp, ...
            Shear.smoothOpts.accel,Shear.extrapOpts.accel,...
            diffOpts,false);
    else
        fprintf('Calculating acceleration from the spatially smoothed displacement fields.\n')
        [Shear.accel,~,temp.disp] = func_smoothCalcAccel_v4(pos,time,Shear.disp.sSmooth, ...
            Shear.smoothOpts.accel,Shear.extrapOpts.accel,...
            diffOpts,false);
        Shear.disp.tSmooth=temp.disp.tSmooth;
        clear temp

    end
else
    fprintf('Calculating acceleration from the unsmoothed displacement fields.\n')
    [Shear.accel,~,Shear.disp] = func_smoothCalcAccel_v4(pos,time,Shear.disp, ...
        Shear.smoothOpts.accel,Shear.extrapOpts.accel,...
        diffOpts,false);
end



%% Calculate Stress Gauge Stresses
fprintf('Calculating Stress Gage Stresses For Shear Identification \n')
Shear.SG=func_Full_SG(Shear.accel,pos.x,time,material.rho);

% %% Smooth Stress Gauge Stresses
% Shear.SG.smoothS=sgolayfilt(Shear.SG.s,3,27);
%% Convert time to microsecond
timeMS=time.vec*10^6;

%% 1-D Shear Diagnostic Plots
    %1/3 from free surface
    FS.ind=round(size(Shear.strain.sAv,1)/3);
    FS.strain.sAv=squeeze(Shear.strain.sAv(FS.ind,:));
    FS.SG.s=squeeze(Shear.SG.s(FS.ind,:))*10^-6;

    %Impact Edge
    IMP.ind=size(Shear.strain.sAv,1)-FS.ind;
    IMP.strain.sAv=squeeze(Shear.strain.sAv(IMP.ind,:));
    IMP.SG.s=squeeze(Shear.SG.s(IMP.ind,:))*10^-6;

    %Middle
    MID.ind=round(size(Shear.strain.sAv,1)/2);
    MID.strain.sAv=squeeze(Shear.strain.sAv(MID.ind,:));
    MID.SG.s=squeeze(Shear.SG.s(MID.ind,:))*10^-6;

    %Plot
    figure(PreFig{:})
    t=tiledlayout(1,3);
    
    %Strain Time
    ax=nexttile;
    plot(timeMS,FS.strain.sAv,'k')
    hold on
    plot(timeMS,MID.strain.sAv,'b')
    plot(timeMS,IMP.strain.sAv,'r')
    hold off
    xlabel('$time~\mathrm{(\mu s)}$',PreLabel{:})
    ylabel('$\overline{\varepsilon_{12}}$',PreLabel{:})
    

    %Stress Time
    nexttile
    plot(timeMS,FS.SG.s,'k')
    hold on
    plot(timeMS,MID.SG.s,'b')
    plot(timeMS,IMP.SG.s,'r')
    hold off
    xlabel('$time~\mathrm{(\mu s)}$',PreLabel{:})
    ylabel('$\overline{\sigma_{12}^\mathrm{SG}}~\mathrm{(MPa)}$',PreLabel{:})

     %Stress Time
    nexttile
    plot(FS.strain.sAv,FS.SG.s,'k')
    hold on
    plot(MID.strain.sAv,MID.SG.s,'b')
    plot(IMP.strain.sAv,IMP.SG.s,'r')
    hold off
    xlabel('$\overline{\varepsilon_{12}}$',PreLabel{:})
    ylabel('$\overline{\sigma_{12}^\mathrm{SG}}~\mathrm{(MPa)}$',PreLabel{:})
    
    lg=legend(ax,'1/3 from FS','Mid','1/3 from Impact');
    lg.Layout.Tile='South';
    lg.Orientation='Horizontal';

    D1FigName=strcat(ExpSavePath,'/',ExpDesig,'_Shear1D_Diagnostic');
    saveas(gcf,D1FigName,'fig')
    saveas(gcf,D1FigName,'svg')
    saveas(gcf,D1FigName,'png')
   
  
%% Full Field Shear Movie
fprintf('Generating Full Field Movies using Shear Processing Parameters \n')
ShearDesig=strcat(ExpDesig,'_Shear');
ShearDir=strcat(ExpSavePath,'/ShearProc');

func_IBIIMakeFFMovies(pos,time,Shear.disp,Shear.strain,Shear.accel, ...
    Shear.smoothOpts,Shear.extrapOpts,diffOpts,ShearDesig,ShearDir);

%% Shear X-t Diagram (X datapoint vs frame)
 figure(PreFig{:})
t=tiledlayout(1,2);
time.Frames=1:length(time.vec);
pos.Xpt=1:length(pos.x);
Shear.strain.sLim=max(abs(Shear.strain.sAv),[],'all')*[-1,1];
Shear.SG.sLim=10*[-1,1];

%Shear Strain
nexttile
%Shear.strain.sAvZoom=Shear.strain.sAv*10^3;
imagesc(pos.Xpt,time.Frames,Shear.strain.sAv')
crameri('-roma')
c=colorbar;
c.Label.Interpreter='latex';
c.Label.String='$\varepsilon_{xy}$';
caxis(Shear.strain.sLim)
set(gca,PreAx{:})

%Shear Stress Gauge
nexttile
imagesc(pos.Xpt,time.Frames,Shear.SG.s'*10^-6)
crameri('-roma')
c=colorbar;
c.Label.Interpreter='latex';
c.Label.String='$\overline{\sigma_{xy}^\mathrm{SG}}~\mathrm{(MPa)}$';
set(gca,PreAx{:})
caxis(Shear.SG.sLim)

xlabel(t,'X pixel',PreLabel{:})
ylabel(t,'Frame Number',PreLabel{:})

FigName=strcat(ExpSavePath,'/',ExpDesig,'_ShearXT_dataPoints');
    saveas(gcf,FigName,'fig')
    saveas(gcf,FigName,'svg')
    saveas(gcf,FigName,'png')

%% Shear X-t Diagram 
 figure(PreFig{:})
t=tiledlayout(1,2);
time.us=time.vec*10^6;
pos.Xmm=pos.x*10^3;
Shear.strain.sLim=max(abs(Shear.strain.sAv),[],'all')*[-1,1];
Shear.SG.sLim=10*[-1,1];

%Shear Strain
nexttile
%Shear.strain.sAvZoom=Shear.strain.sAv*10^3;
imagesc(pos.Xmm,time.us,Shear.strain.sAv')
crameri('-roma')
c=colorbar;
title('$\overline{\varepsilon_{12}}$',PreLabel{:})
% c.Label.Interpreter='latex';
% c.Label.String='$\varepsilon_{xy}$';
caxis([-1,1]*15e-3)
yticks([0,25,50,75,100,125])
set(gca,PreAx{:})

%Shear Stress Gauge
nexttile
imagesc(pos.Xmm,time.us,Shear.SG.s'*10^-6)
crameri('-roma')
c=colorbar;
title('$\overline{\sigma_{xy}^\mathrm{SG}}~\mathrm{[MPa]}$',PreLabel{:})
% c.Label.Interpreter='latex';
% c.Label.String='$\overline{\sigma_{xy}^\mathrm{SG}}~\mathrm{(MPa)}$';
set(gca,PreAx{:})
caxis(Shear.SG.sLim)
yticks([0,25,50,75,100,125]);

xlabel(t,'$x~\mathrm{(mm)}$',PreLabel{:})
ylabel(t,'$t~\mathrm{(\mu s)}$',PreLabel{:})

FigName=strcat(ExpSavePath,'/',ExpDesig,'_ShearXT');
    saveas(gcf,FigName,'fig')
    saveas(gcf,FigName,'svg')
    saveas(gcf,FigName,'png')

%% Determine Time and Space Window to evaluate the cost function for
fprintf('Adjusting the censorship of the data before solving the cost function \n')
fprintf(strcat('Total number of X data Points=',...
    num2str(length(pos.x)),'\n'));
prompt={'Number of frames to remove from the beginning of the analysis';
    'Number of frames to remove from the end of the analysis';
    'Number of pixels to remove from the impact edge';
    'Number of pixels to remove from the free edge'};
defImp=num2str(50);
defFree=num2str(10);
if defFree<Shear.smoothOpts.strain.spatialKernelSize(1)/2
    defFree=num2str(Shear.smoothOpts.strain.spatialKernelSize(1)/2);
end

defin={'5';'5';defImp;defFree};
dims=1;
dlgtitle='Choose Data Censorship Parameters';
dlgOpts.Resize='on';
dlgOpts.WindowStyle='normal';
tempOpts=inputdlg(prompt,dlgtitle,dims,defin,dlgOpts);

CondOpts.CutStartFrames=str2double(tempOpts{1});
CondOpts.CutEndFrames=str2double(tempOpts{2});
CondOpts.ImpCens=str2double(tempOpts{3});
CondOpts.FreeCens=str2double(tempOpts{4});

%% Run the minimization 
fprintf('-----------------------------------------------------------------\n')
fprintf('RUNNING SHEAR MINIMIZATION \n')
fprintf('------------------------------------------------------------------\n')

%% Condition the data (censoring/downsampling/smoothing)
fprintf('Conditioning IBII Fields \n')
[Shear.SGcens,~,Shear.strainCens,X_vec,timeCens]=...
    func_conditionIBIIDataV2(Shear.SG,Shear.accel,...
    Shear.strain,pos.x,time,CondOpts);

%% Run the minimization
[TempG,TempTau,TempPhi]=deal(zeros(3,1));
exactProps=material;
exactProps.Ki=material.K1_eff;
for ig=1:3
    fprintf(strcat('Running Shear Iteration ',num2str(ig),'/3 \n'))
    intGuess=squeeze(intMatrixG(:,ig)); %#ok<PFBNS>
    [TempParam,TempPhi(ig)]=func_PronyShearSGvfm(Shear.strainCens.s, ...
        timeCens.vec,Shear.SGcens.s,exactProps,...
        intGuess,ubG,lbG,SolveOpts,minOpts,CondOpts);
    TempG(ig)=TempParam(1);
    TempTau(ig)=TempParam(2);
end

minPhi=min(TempPhi);
G=TempG(TempPhi==minPhi);
tauG=TempTau(TempPhi==minPhi);
phi=TempPhi';
ConstitutiveParam=[G(1);tauG(1)];

Ident.G1=G(1);
Ident.tauG=tauG(1);

fprintf('SHEAR MODULUS IDENTIFICATION COMPLETE \n')


%% Bulk Modulus 

fprintf('---------------------------------------------------------------- \n')
fprintf('STARTING BULK MODULUS IDENTIFICATION \n')
fprintf('------------------------------------------------------------------ \n')

%% Bulk Modulus Kinematic field calculation
% Extrapolate the data to account for the missing pitch on the edges
% Crop and Extrapolate displacement fields

fprintf('BULK MODULUS KINEMATIC FIELD CALCULATION \n')

fprintf('Croping and Extrapolating Displacement Fields \n')
[Bulk.disp.x,Bulk.disp.y,Bulk.disp.rX,Bulk.disp.rY]=...
    func_cropAndExtrapDispFields_v4(pos.Raw,disp.x,disp.y, ...
    Bulk.extrapOpts.disp,false);

switch globalOpts.dispSource
    case 'MatchID'
        fprintf('Extrapolating Pos to be the same size as disp \n')
        pos=func_ExtrapolatePos(pos.Raw,DIC.edgeCut);
end

fprintf('Obtaining and setting the free edge location.\n')
switch globalOpts.hardCodeFreeEdge
    case false
    [freeEdge,specimen,Bulk.disp] =...
    func_getFreeEdge(globalOpts.hardCodeFreeEdge,...
    imagePath,imageFile,specimen,Bulk.disp);
    case true
        [freeEdge,specimen,Bulk.disp] =...
    func_CheckFreeEdge(specimen,Bulk.disp);
end

switch globalOpts.dispSource
    case 'MatchID'
       
        %Steps to crop
        cST=Bulk.extrapOpts.disp.cropPx1st;
        %Steps to extrapolate
        eST=Bulk.extrapOpts.disp.extrapPx1st;
        
        %Data Points missing from edge
        missST=DIC.edgeCut;

        %Determine range to crop back to full specimen size
        Bulk.disp.rX=(1+eST-cST-missST):(size(Bulk.disp.x,2)...
            -eST+cST+missST);
        Bulk.disp.rY=(1+eST-cST-missST):(size(Bulk.disp.x,1)...
            -eST+cST+missST);
end

if globalOpts.ExtractDICStrain==false
    fprintf('Calculating Strains for Bulk Identification \n')
    [Bulk.strain,Bulk.disp]=func_smoothCalcStrain_v4(pos,time,Bulk.disp,...
        Bulk.smoothOpts.strain,Bulk.extrapOpts.strain,false);
end

%--------------------------------------------------------------------------
% Smooth and Calculate Acceleration
if isfield(Bulk.smoothOpts.accel,'spatialSmooth')
    if Bulk.smoothOpts.accel.spatialSmooth==false
        fprintf('Calculating acceleration from the unsmoothed displacement fields.\n')
        [Bulk.accel,~,Bulk.disp] = func_smoothCalcAccel_v4(pos,time, ...
            Bulk.disp, ...
            Bulk.smoothOpts.accel,Bulk.extrapOpts.accel,...
            diffOpts,false);
    else
        fprintf('Calculating acceleration from the spatially smoothed displacement fields.\n')
        [Bulk.accel,~,temp.disp] = func_smoothCalcAccel_v4(pos,time, ...
            Bulk.disp.sSmooth, ...
            Bulk.smoothOpts.accel,Bulk.extrapOpts.accel,...
            diffOpts,false);
            Bulk.disp.tSmooth=temp.disp.tSmooth;
            clear temp
    end
else
    fprintf('Calculating acceleration from the unsmoothed displacement fields.\n')
    [Bulk.accel,~,Bulk.disp] = func_smoothCalcAccel_v4(pos,time, ...
        Bulk.disp, ...
        Bulk.smoothOpts.accel,Bulk.extrapOpts.accel,...
        diffOpts,false);
end


%% Calculate Stress Gauge Stresses
fprintf('Calculating Stress Gage Stresses For Bulk Modulus Identification \n')
Bulk.SG=func_Full_SG(Bulk.accel,pos.x,time,material.rho);


%% Bulk 1-D Diagnostic Figures
fprintf('BULK DIAGNOTIC FIGURES \n')
 %1/3 from free surface
    FS.ind=round(size(Bulk.strain.xAvg,1)/3);
    FS.strain.xAvg=squeeze(Bulk.strain.xAvg(FS.ind,:));
    FS.strain.yAvg=squeeze(Bulk.strain.yAvg(FS.ind,:));
    FS.SG.x=squeeze(Bulk.SG.x(FS.ind,:))*10^-6;

    %Impact Edge
    IMP.ind=size(Bulk.strain.xAvg,1)-FS.ind;
    IMP.strain.xAvg=squeeze(Bulk.strain.xAvg(IMP.ind,:));
    IMP.strain.yAvg=squeeze(Bulk.strain.yAvg(IMP.ind,:));
    IMP.SG.x=squeeze(Bulk.SG.x(IMP.ind,:))*10^-6;

    %Middle
    MID.ind=round(size(Bulk.strain.xAvg,1)/2);
    MID.strain.xAvg=squeeze(Bulk.strain.xAvg(MID.ind,:));
    MID.strain.yAvg=squeeze(Bulk.strain.yAvg(MID.ind,:));
    MID.SG.x=squeeze(Shear.SG.x(MID.ind,:))*10^-6;

    %Plot
    figure(PreFig{:})
    t=tiledlayout(1,3);
    
    %Strain Time
    ax=nexttile;
    plot(timeMS,FS.strain.xAvg,'k')
    hold on
    plot(timeMS,MID.strain.xAvg,'b')
    plot(timeMS,IMP.strain.xAvg,'r')
    hold off
    xlabel('$time~\mathrm{(\mu s)}$',PreLabel{:})
    ylabel('$\overline{\varepsilon_{11}}$',PreLabel{:})
    

    %Stress Time
    nexttile
    plot(timeMS,FS.SG.x,'k')
    hold on
    plot(timeMS,MID.SG.x,'b')
    plot(timeMS,IMP.SG.x,'r')
    hold off
    ylabel('$\overline{\sigma_{11}^\mathrm{SG}}~\mathrm{(MPa)}$',PreLabel{:})

     %Stress Time
    nexttile
    plot(FS.strain.xAvg,FS.SG.x,'k')
    hold on
    plot(MID.strain.xAvg,MID.SG.x,'b')
    plot(IMP.strain.xAvg,IMP.SG.x,'r')
    hold off
    xlabel('$\overline{\varepsilon_{11}}$',PreLabel{:})
    ylabel('$\overline{\sigma_{11}^\mathrm{SG}}~\mathrm{(MPa)}$',PreLabel{:})
    
    lg=legend(ax,'1/3 from FS','Mid','1/3 from Impact');
    lg.Layout.Tile='South';
    lg.Orientation='Horizontal';

    D1FigName=strcat(ExpSavePath,'/',ExpDesig,'_Bulk1D_Diagnostic');
    saveas(gcf,D1FigName,'fig')
    saveas(gcf,D1FigName,'svg')
    saveas(gcf,D1FigName,'png')
   
    %% 1D with YY 

    %Plot
    figure(PreFig{:})
    t=tiledlayout(1,3);
    
    %Strain Time
    ax=nexttile;
    plot(timeMS,FS.strain.xAvg,'k')
    hold on
    plot(timeMS,MID.strain.xAvg,'b')
    plot(timeMS,IMP.strain.xAvg,'r')
    hold off
    xlabel('$time~\mathrm{(\mu s)}$',PreLabel{:})
    ylabel('$\overline{\varepsilon_{11}}$',PreLabel{:})
    

    %Stress Time
    nexttile
    plot(timeMS,FS.strain.yAvg,'k')
    hold on
    plot(timeMS,MID.strain.yAvg,'b')
    plot(timeMS,IMP.strain.yAvg,'r')
    hold off
    xlabel('$time~\mathrm{(\mu s)}$',PreLabel{:})
    ylabel('$\overline{\varepsilon_{22}}$',PreLabel{:})

     %Stress Time
    nexttile
    plot(FS.strain.xAvg,FS.strain.yAvg,'k')
    hold on
    plot(MID.strain.xAvg,MID.strain.yAvg,'b')
    plot(IMP.strain.xAvg,FS.strain.yAvg,'r')
    hold off
    xlabel('$\overline{\varepsilon_{11}}$',PreLabel{:})
    ylabel('$\overline{\varepsilon_{22}}$',PreLabel{:})
    
    lg=legend(ax,'1/3 from FS','Mid','1/3 from Impact');
    lg.Layout.Tile='South';
    lg.Orientation='Horizontal';

    D1FigName=strcat(ExpSavePath,'/',ExpDesig,'_Bulk1DStrains_Diagnostic');
    saveas(gcf,D1FigName,'fig')
    saveas(gcf,D1FigName,'svg')
    saveas(gcf,D1FigName,'png')

%% Bulk X-t Diagram (X datapoint vs frame)
 figure(PreFig2{:})
t=tiledlayout(4,4);
time.Frames=1:length(time.vec);
pos.Xpt=1:length(pos.x);
Bulk.strain.xLim=max(abs(Bulk.strain.xAvg),[],'all')*[-1,1];
Bulk.strain.yLim=max(abs(Bulk.strain.yAvg),[],'all')*[-1,1];
Bulk.SG.xLim=40*[-1,1];

%Bulk Axial Strain
nexttile([2,2])
%Shear.strain.sAvZoom=Shear.strain.sAv*10^3;
imagesc(pos.Xpt,time.Frames,Bulk.strain.xAvg')
crameri('-roma')
c=colorbar;
c.Label.Interpreter='latex';
c.Label.String='$\varepsilon_{11}$';
caxis(Bulk.strain.xLim)
set(gca,PreAx{:})

%Bulk Axial Stress
nexttile([2,2])
imagesc(pos.Xpt,time.Frames,Bulk.SG.x'*10^-6)
crameri('-roma')
c=colorbar;
c.Label.Interpreter='latex';
c.Label.String='$\overline{\sigma_{11}^\mathrm{SG}}~\mathrm{(MPa)}$';
set(gca,PreAx{:})
caxis(Bulk.SG.xLim)

%Bulk Transverse Setress
nexttile(10,[2,2])
%Shear.strain.sAvZoom=Shear.strain.sAv*10^3;
imagesc(pos.Xpt,time.Frames,Bulk.strain.yAvg')
crameri('-roma')
c=colorbar;
c.Label.Interpreter='latex';
c.Label.String='$\varepsilon_{22}$';
caxis(Bulk.strain.yLim)
set(gca,PreAx{:})
xlabel(t,'X pixel',PreLabel{:})
ylabel(t,'Frame Number',PreLabel{:})

FigName=strcat(ExpSavePath,'/',ExpDesig,'_BulkXT_DataPoints');
    saveas(gcf,FigName,'fig')
    saveas(gcf,FigName,'svg')
    saveas(gcf,FigName,'png')


 %% Bulk X-t Diagram 
figure(PreFig2{:})
t=tiledlayout(4,4);
time.Frames=1:length(time.vec);
pos.Xpt=1:length(pos.x);
Bulk.strain.xLim=max(abs(Bulk.strain.xAvg),[],'all')*[-1,1];
Bulk.strain.yLim=max(abs(Bulk.strain.yAvg),[],'all')*[-1,1];
Bulk.SG.xLim=40*[-1,1];

%Bulk Axial Strain
nexttile([2,2])
%Shear.strain.sAvZoom=Shear.strain.sAv*10^3;
imagesc(pos.Xmm,timeMS,Bulk.strain.xAvg')
crameri('-roma')
c=colorbar;
c.Label.Interpreter='latex';
c.Label.String='$\varepsilon_{11}$';
caxis(Bulk.strain.xLim)
set(gca,PreAx{:})

%Bulk Axial Stress
nexttile([2,2])
imagesc(pos.Xmm,timeMS,Bulk.SG.x'*10^-6)
crameri('-roma')
c=colorbar;
c.Label.Interpreter='latex';
c.Label.String='$\overline{\sigma_{11}^\mathrm{SG}}~\mathrm{(MPa)}$';
set(gca,PreAx{:})
caxis(Bulk.SG.xLim)

%Bulk Transverse Setress
nexttile(10,[2,2])
%Shear.strain.sAvZoom=Shear.strain.sAv*10^3;
imagesc(pos.Xmm,timeMS,Bulk.strain.yAvg')
crameri('-roma')
c=colorbar;
c.Label.Interpreter='latex';
c.Label.String='$\varepsilon_{22}$';
caxis(Bulk.strain.yLim)
set(gca,PreAx{:})

xlabel(t,'$x~\mathrm{(mm)}$',PreLabel{:})
ylabel(t,'$t~\mathrm{(\mu s)}$',PreLabel{:})

FigName=strcat(ExpSavePath,'/',ExpDesig,'_BulkXT');
    saveas(gcf,FigName,'fig')
    saveas(gcf,FigName,'svg')
    saveas(gcf,FigName,'png')

 %% Bulk X-t Diagram 
figure(PreFig{:})
t=tiledlayout(1,2);
time.Frames=1:length(time.vec);
pos.Xpt=1:length(pos.x);
Bulk.strain.xLim=max(abs(Bulk.strain.xAvg),[],'all')*[-1,1];
Bulk.strain.yLim=max(abs(Bulk.strain.yAvg),[],'all')*[-1,1];
Bulk.SG.xLim=40*[-1,1];

%Bulk Axial Strain
nexttile
%Shear.strain.sAvZoom=Shear.strain.sAv*10^3;
imagesc(pos.Xmm,timeMS,Bulk.strain.xAvg')
crameri('-roma')
title('$\overline{\varepsilon_{11}}$',PreLabel{:})
c=colorbar;
% c.Label.Interpreter='latex';
% c.Label.String='$\overline{\varepsilon_{11}}$';
caxis([-1,1]*50e-3)
set(gca,PreAx{:})
yticks([0,25,50,75,100,125])

%Bulk Axial Stress
nexttile
imagesc(pos.Xmm,timeMS,Bulk.SG.x'*10^-6)
crameri('-roma')
c=colorbar;
title('$\overline{\sigma_{11}^\mathrm{SG}}~\mathrm{[MPa]}$',PreLabel{:})
% c.Label.Interpreter='latex';
% c.Label.String='$\overline{\sigma_{11}^\mathrm{SG}}~\mathrm{[MPa]}$';
set(gca,PreAx{:})
caxis([-100,100])
yticks([0,25,50,75,100,125])

xlabel(t,'$x~\mathrm{(mm)}$',PreLabel{:})
ylabel(t,'$t~\mathrm{(\mu s)}$',PreLabel{:})

FigName=strcat(ExpSavePath,'/',ExpDesig,'_BulkXT');
    saveas(gcf,FigName,'fig')
    saveas(gcf,FigName,'svg')
    saveas(gcf,FigName,'png')
    %% Full Field Movies
    fprintf('FULL FIELD MOVIES \n')

fprintf('Generating Full Field Movies using Bulk Processing Parameters \n')
BulkDesig=strcat(ExpDesig,'_Bulk');
BulkDir=strcat(ExpSavePath,'/BulkProc');

func_IBIIMakeFFMovies(pos,time,Bulk.disp,Bulk.strain,Bulk.accel, ...
    Bulk.smoothOpts,Bulk.extrapOpts,diffOpts,BulkDesig,BulkDir);

%% Extract Bulk Modulus and Time Constant
fprintf('-------------------------------------------------------------- \n')
fprintf('BULK MODULUS IDENTIFICATION \n')

%% Condition the data (censoring/downsampling/smoothing)
fprintf('Conditioning Bulk IBII Fields \n')
[Bulk.SGcens,~,Bulk.strainCens,X_vec,timeCens]=...
    func_conditionIBIIDataV2(Bulk.SG,Bulk.accel,...
    Bulk.strain,pos.x,time,CondOpts);

% 
% %% Plot Axial and Shear stress gauges vertically
% 
% %Bulk Axial Stress
% nexttile
% imagesc(pos.Xmm,timeMS,Bulk.SG.x'*10^-6)
% crameri('-roma')
% c=colorbar;
% title('$\overline{\sigma_{11}^\mathrm{SG}}~\mathrm{[MPa]}$',PreLabel{:})
% % c.Label.Interpreter='latex';
% % c.Label.String='$\overline{\sigma_{11}^\mathrm{SG}}~\mathrm{[MPa]}$';
% set(gca,PreAx{:})
% caxis([-100,100])
% yticks([0,25,50,75,100,125])

%% Set the material properties
BulkProps=material;
BulkProps.Gi=Ident.G1;
 
[TempK,TempTau,TempPhi]=...
                deal(zeros(3,1));
            for ig=1:3
                fprintf(strcat('K ident iteration',num2str(ig),'/3 \n'))
                intGuess=squeeze(intMatrixK(:,ig)); 
                [TempParam,TempPhi(ig)]=func_PronyBulkTauSGvfm(Bulk.strainCens, ...
                    timeCens.vec,Bulk.SGcens.x,BulkProps,...
                    intGuess,ubK,lbK,SolveOpts,minOpts,CondOpts);
                TempK(ig)=TempParam(1);
                TempTau(ig)=TempParam(2);

            end

%             minPhiE=min(TempPhiE);
%             KE=TempKE(TempPhiE==minPhiE);
%             tauE=TempTauE(TempPhiE==minPhiE);
%             phiE(k,m,:)=TempPhiE';
%             ConstitutiveParamE(k,m,:)=[KE(1);tauE(1)];

            minPhi=min(TempPhi);
            K=TempK(TempPhi==minPhi);
            tauK=TempTau(TempPhi==minPhi);
            phiK(:)=TempPhi';
            ConstitutiveParam(:)=[K(1);tauK(1)];

            Ident.K1=K(1);
            Ident.tau=tauK(1);

fprintf('')            
%% Save Results
ExpSaveName=strcat(ExpSavePath,'/',ExpDesig,'KG_Ident_Results.mat');
save(ExpSaveName)
