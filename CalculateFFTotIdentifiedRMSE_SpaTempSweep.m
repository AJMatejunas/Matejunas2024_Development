%This script is written to calculate the total root mean square error using
    %identified constitutive parameters over a sweep of spatial and
    %temporal smoothing kernels on noisy images. 

% RMSE is calculated from the full-field finite element
    %stresses along with the full field constitutive model stresses
    %calculated from full field finite element strains at the full
    %resolution.

%Author: Andrew Matejunas

%Version History/Changelog

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize
clear variables; close all; clc

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
%Full Page Portrait Figure
FPpFig={'Units','Centimeters','InnerPosition',[1,1,DCW,20]};
%Full Page Landscape
FPlFig={'Units','Centimeters','InnerPosition',[0.5,0.5,230,DCW]};


%% Choose file with full field finite element time, strain, and stress arrays
[FEname,FEpath]=uigetfile('*.mat', ...
    'Choose file containing full-field FE kinematic fields');
FEfile=strcat(FEpath,'/',FEname);

%% Choose file containing identified constitutive paramters
[ShearName,ShearPath]=uigetfile('*.mat', ...
    'Choose file containing identified shear parameters');
ShearFile=strcat(ShearPath,'/',ShearName);

%% Choose file containing identified bulk parameters
[BulkName,BulkPath]=uigetfile('*.mat', ...
    'Choose file containing identified Bulk parameters');
BulkFile=strcat(BulkPath,'/',BulkName);

%% Choose Reference Parameter file
[RefName,RefPath]=uigetfile('*.mat', ...
    'Choose file containing reference parameters');
RefFile=strcat(RefPath,'/',RefName)
%% Choose location to save results
SaveDir=uigetdir('','Choose directory to save results in');

%% Input Test Designation
Desig=char(cell2mat(inputdlg('Input Test Designation')));

%% Load files;

fprintf('Loading FE fields \n')
FE=load(FEfile,'strain','stress','time');

fprintf('Loading Identified shear parameters \n')
Shear=load(ShearFile,'Ident','SpaKernVec','TempKernVec');

fprintf('Loading Identified bulk parameters \n')
Bulk=load(BulkFile,'Ident','SpaKernVec','TempKernVec');

fprintf('Loading Reference properties \n')
load(RefFile);

%% Input Identified shear
ShearIdent=str2double(cell2mat(inputdlg('Enter Identified Shear')));

%% Extract stresses, strains, and time to form that can be used in parfor
timevec=FE.time.vec;
strainX=FE.strain.x;
strainY=FE.strain.y;
strainS=FE.strain.s;

stressX=FE.stress.x;
stressY=FE.stress.y;
stressS=FE.stress.s;
fprintf('Calculating FE Hydrostatic Stress \n')
stressH=(FE.stress.x+FE.stress.y)/3;

fprintf('Extracting Maximum stresses \n')
MaxStressX=max(abs(FE.stress.x),[],"all");
MaxStressY=max(abs(FE.stress.y),[],'all');
MaxStressS=max(abs(FE.stress.s),[],"all");
MaxStressH=max(abs(stressH),[],'all');

%% Set up MatPropsFile
ShearPropsArray=cell(size(Shear.Ident.G));
ShearG=Shear.Ident.G;
ShearTau=Shear.Ident.tauG;

fprintf('Generating Shear Property Input File \n')
parfor k=1:numel(Shear.Ident.G)
TempProps=MatProps;
TempProps.Gi=ShearG(k);
TempProps.tau=ShearTau(k);
TempProps.G0=MatProps.Ginf+ShearG(k);
ShearPropsArray{k}=TempProps;
end

clear TempProps
%% calculate shear RSME 
ShearRMSE=zeros(size(Shear.Ident.G));
fprintf('Calculating Shear RMSE \n')
parfor m=1:numel(ShearRMSE)
    % Calculate Constitutive Model stresses
    ShearProps=cell2mat(ShearPropsArray(m));
    %Calculate Shear Stress
    ShearModel=func_ViscoShearConstitutive(strainS,timevec,ShearProps);

    %Shear RMSE
    ShearRMSE(m)=func_calcFFRMSE(stressS,ShearModel.xy);
end
clear ShearModel
%% Consolodate Shear RMSE to Single Variable
RMSE.Shear=ShearRMSE;
RMSE.sys.shear=squeeze(mean(ShearRMSE,3));
RMSE.ran.shear=squeeze(std(ShearRMSE,0,3));
RMSE.tot.shear=abs(RMSE.sys.shear)+2*RMSE.ran.shear;

RMSE.ShearNorm=RMSE.Shear/MaxStressS;
RMSE.sys.shearNorm=squeeze(mean(RMSE.ShearNorm,3));
RMSE.ran.shearNorm=squeeze(std(RMSE.ShearNorm,0,3));
RMSE.tot.shearNorm=abs(RMSE.sys.shearNorm)+2*RMSE.ran.shearNorm;

%% Plot Heat Map of Shear RMSE Errors
figure(TcSrFig{:})
t=tiledlayout(1,3);

Spaticks=[0,25,50];
% Systematic Errors
nexttile
imagesc(Shear.SpaKernVec,Shear.TempKernVec,RMSE.sys.shearNorm)
crameri bilbao
set(gca,axprops{:})
colorbar
%clim(5*[-1,1])
title('(a)~$||\mathrm{RMSE_{sys}}(\sigma_{12}^\mathrm{Model})||_{\mathrm{max}}$ [\%]', ...
    labelProps{:})
xticks(Spaticks);


% Random Error
nexttile
imagesc(Shear.SpaKernVec,Shear.TempKernVec,RMSE.ran.shearNorm)
crameri bilbao
set(gca,axprops{:})
colorbar
%clim(5*[-1,1])
title('(b)~$||\mathrm{RMSE_{ran}}(\sigma_{12}^\mathrm{Model})||_{\mathrm{max}}$ [\%]', ...
    labelProps{:})
xticks(Spaticks);



% Total Error
nexttile
imagesc(Shear.SpaKernVec,Shear.TempKernVec,RMSE.tot.shearNorm)
crameri bilbao
set(gca,axprops{:})
colorbar
%clim(5*[-1,1])
title('(c)~$||\mathrm{RMSE_{tot}}(\sigma_{12}^\mathrm{Model})||_{\mathrm{max}}$ [\%]', ...
    labelProps{:})
xticks(Spaticks);

ylabel(t,'$TK$~[frames]',MajorProps{:})
xlabel(t,'$SK$~[pixels]',MajorProps{:})

FigName=strcat(SaveDir,'/',Desig,'_shearRMSE');
saveas(gcf,FigName,'fig');
saveas(gcf,FigName,'svg');
saveas(gcf,FigName,'eps');


%% Plot Heat Map of Shear RMSE Errors
figure(TcSrFig{:})
t=tiledlayout(1,3);

Spaticks=[0,25,50];
% Systematic Errors
nexttile
imagesc(Shear.SpaKernVec,Shear.TempKernVec,RMSE.sys.shear)
crameri bilbao
set(gca,axprops{:})
colorbar
%clim(5*[-1,1])
title('(a)~$\mathrm{RMSE_{sys}}(\sigma_{12}^\mathrm{Model})$ [Pa]', ...
    labelProps{:})
xticks(Spaticks);


% Random Error
nexttile
imagesc(Shear.SpaKernVec,Shear.TempKernVec,RMSE.ran.shear)
crameri bilbao
set(gca,axprops{:})
colorbar
%clim(5*[-1,1])
title('(b)~$\mathrm{RMSE_{ran}}(\sigma_{12}^\mathrm{Model})$ [Pa]', ...
    labelProps{:})
xticks(Spaticks);



% Total Error
nexttile
imagesc(Shear.SpaKernVec,Shear.TempKernVec,RMSE.tot.shear)
crameri bilbao
set(gca,axprops{:})
colorbar
%clim(5*[-1,1])
title('(c)~$\mathrm{RMSE_{tot}}(\sigma_{12}^\mathrm{Model})$ [Pa]', ...
    labelProps{:})
xticks(Spaticks);

ylabel(t,'$TK$~[frames]',MajorProps{:})
xlabel(t,'$SK$~[pixels]',MajorProps{:})

FigName=strcat(SaveDir,'/',Desig,'_shearRMSE_NoNorm');
saveas(gcf,FigName,'fig');
saveas(gcf,FigName,'svg');
saveas(gcf,FigName,'eps');

%% Set up Bulk MatPropsFile
BulkPropsArray=cell(size(Shear.Ident.G));
BulkK=Bulk.Ident.K;
BulkTau=Bulk.Ident.tau;

fprintf('Generating Bulk Property Input File \n')
parfor k=1:numel(Shear.Ident.G)
TempProps=MatProps;
TempProps.Gi=ShearIdent;
TempProps.G0=MatProps.Ginf+ShearIdent;
TempProps.Ki=BulkK(k);
TempProps.K0=BulkK(k)+MatProps.Kinf;
TempProps.tau=BulkTau(k);
BulkPropsArray{k}=TempProps;
end

clear TempProps
%% calculate Bulk RMSE 
AxialRMSE=zeros(size(Bulk.Ident.K));
HydroRMSE=zeros(size(Bulk.Ident.K));

fprintf('Calculating Axial and Hydrostatic RMSE \n')
parfor m=1:numel(AxialRMSE)
    % Calculate Constitutive Model stresses
    BulkProps=cell2mat(BulkPropsArray(m));
    %Calculate COnsititutive Model Stresses
    BulkModel=func_ViscoConstitutiveV6(strainX,strainY,strainS,timevec, ...
        BulkProps,0,0,0);

    %Shear RMSE
    AxialRMSE(m)=func_calcFFRMSE(stressX,BulkModel.xx);
    HydroRMSE(m)=func_calcFFRMSE(stressH,BulkModel.hy);
end
clear BulkModel
%% Consolodate Axial RMSE to Single Variable
RMSE.Axial=AxialRMSE;
RMSE.sys.Axial=squeeze(mean(AxialRMSE,3));
RMSE.ran.Axial=squeeze(std(AxialRMSE,0,3));
RMSE.tot.Axial=abs(RMSE.sys.Axial)+2*RMSE.ran.Axial;

RMSE.AxialNorm=RMSE.Axial/MaxStressX;
RMSE.sys.AxialNorm=squeeze(mean(RMSE.AxialNorm,3));
RMSE.ran.AxialNorm=squeeze(std(RMSE.AxialNorm,0,3));
RMSE.tot.AxialNorm=abs(RMSE.sys.AxialNorm)+2*RMSE.ran.AxialNorm;

%% Consolodate Hydrostatic RMSE to Single Variable
RMSE.Hydro=HydroRMSE;
RMSE.sys.Hydro=squeeze(mean(HydroRMSE,3));
RMSE.ran.Hydro=squeeze(std(HydroRMSE,0,3));
RMSE.tot.Hydro=abs(RMSE.sys.Hydro)+2*RMSE.ran.Hydro;

RMSE.HydroNorm=RMSE.Hydro/MaxStressH;
RMSE.sys.HydroNorm=squeeze(mean(RMSE.HydroNorm,3));
RMSE.ran.HydroNorm=squeeze(std(RMSE.HydroNorm,0,3));
RMSE.tot.HydroNorm=abs(RMSE.sys.HydroNorm)+2*RMSE.ran.HydroNorm;

%% Plot Heat Map of Bulk RMSE Errors
figure(TcDrFig{:})
t=tiledlayout(2,3);

Spaticks=[0,25,50];
% Systematic Errors
nexttile
imagesc(Bulk.SpaKernVec,Bulk.TempKernVec,RMSE.sys.AxialNorm)
crameri bilbao
set(gca,axprops{:})
colorbar
%clim(5*[-1,1])
title('(a)~$||\mathrm{RMSE_{sys}}(\sigma_{11}^\mathrm{Model})||_{\mathrm{max}}$ [\%]', ...
    labelProps{:})
xticks(Spaticks);


% Random Error
nexttile
imagesc(Bulk.SpaKernVec,Bulk.TempKernVec,RMSE.ran.AxialNorm)
crameri bilbao
set(gca,axprops{:})
colorbar
%clim(5*[-1,1])
title('(b)~$||\mathrm{RMSE_{ran}}(\sigma_{11}^\mathrm{Model})||_{\mathrm{max}}$ [\%]', ...
    labelProps{:})
xticks(Spaticks);



% Total Error
nexttile
imagesc(Bulk.SpaKernVec,Bulk.TempKernVec,RMSE.tot.AxialNorm)
crameri bilbao
set(gca,axprops{:})
colorbar
%clim(5*[-1,1])
title('(c)~$||\mathrm{RMSE_{tot}}(\sigma_{11}^\mathrm{Model})||_{\mathrm{max}}$ [\%]', ...
    labelProps{:})
xticks(Spaticks);

%HYDROSTATIC STRESSES
% Systematic Errors
nexttile
imagesc(Bulk.SpaKernVec,Bulk.TempKernVec,RMSE.sys.HydroNorm)
crameri bilbao
set(gca,axprops{:})
colorbar
%clim(5*[-1,1])
title('(d)~$||\mathrm{RMSE_{sys}}(\sigma_{h}^\mathrm{Model})||_{\mathrm{max}}$ [\%]', ...
    labelProps{:})
xticks(Spaticks);


% Random Error
nexttile
imagesc(Bulk.SpaKernVec,Bulk.TempKernVec,RMSE.ran.HydroNorm)
crameri bilbao
set(gca,axprops{:})
colorbar
%clim(5*[-1,1])
title('(e)~$||\mathrm{RMSE_{ran}}(\sigma_{h}^\mathrm{Model})||_{\mathrm{max}}$ [\%]', ...
    labelProps{:})
xticks(Spaticks);



% Total Error
nexttile
imagesc(Bulk.SpaKernVec,Bulk.TempKernVec,RMSE.tot.HydroNorm)
crameri bilbao
set(gca,axprops{:})
colorbar
%clim(5*[-1,1])
title('(f)~$||\mathrm{RMSE_{tot}}(\sigma_{h}^\mathrm{Model})||_{\mathrm{max}}$ [\%]', ...
    labelProps{:})
xticks(Spaticks);


ylabel(t,'$TK$~[frames]',MajorProps{:})
xlabel(t,'$SK$~[pixels]',MajorProps{:})

FigName=strcat(SaveDir,'/',Desig,'_BulkRMSE_Norm');
saveas(gcf,FigName,'fig');
saveas(gcf,FigName,'svg');
saveas(gcf,FigName,'eps');


%% Plot Heat Map of Bulk RMSE Errors
figure(TcDrFig{:})
t=tiledlayout(2,3);

Spaticks=[0,25,50];
% Systematic Errors
nexttile
imagesc(Bulk.SpaKernVec,Bulk.TempKernVec,RMSE.sys.Axial)
crameri bilbao
set(gca,axprops{:})
colorbar
%clim(5*[-1,1])
title('(a)~$\mathrm{RMSE_{sys}}(\sigma_{11}^\mathrm{Model})$ [Pa]', ...
    labelProps{:})
xticks(Spaticks);


% Random Error
nexttile
imagesc(Bulk.SpaKernVec,Bulk.TempKernVec,RMSE.ran.Axial)
crameri bilbao
set(gca,axprops{:})
colorbar
%clim(5*[-1,1])
title('(b)~$\mathrm{RMSE_{ran}}(\sigma_{11}^\mathrm{Model})$ [Pa]', ...
    labelProps{:})
xticks(Spaticks);



% Total Error
nexttile
imagesc(Bulk.SpaKernVec,Bulk.TempKernVec,RMSE.tot.Axial)
crameri bilbao
set(gca,axprops{:})
colorbar
%clim(5*[-1,1])
title('(c)~$\mathrm{RMSE_{tot}}(\sigma_{11}^\mathrm{Model})$ [Pa]', ...
    labelProps{:})
xticks(Spaticks);

%HYDROSTATIC STRESSES
% Systematic Errors
nexttile
imagesc(Bulk.SpaKernVec,Bulk.TempKernVec,RMSE.sys.Hydro)
crameri bilbao
set(gca,axprops{:})
colorbar
%clim(5*[-1,1])
title('(d)~$\mathrm{RMSE_{sys}}(\sigma_{h}^\mathrm{Model})$ [Pa]', ...
    labelProps{:})
xticks(Spaticks);


% Random Error
nexttile
imagesc(Bulk.SpaKernVec,Bulk.TempKernVec,RMSE.ran.Hydro)
crameri bilbao
set(gca,axprops{:})
colorbar
%clim(5*[-1,1])
title('(e)~$\mathrm{RMSE_{ran}}(\sigma_{h}^\mathrm{Model})$ [Pa]', ...
    labelProps{:})
xticks(Spaticks);



% Total Error
nexttile
imagesc(Bulk.SpaKernVec,Bulk.TempKernVec,RMSE.tot.Hydro)
crameri bilbao
set(gca,axprops{:})
colorbar
%clim(5*[-1,1])
title('(f)~$\mathrm{RMSE_{tot}}(\sigma_{h}^\mathrm{Model})$ [Pa]', ...
    labelProps{:})
xticks(Spaticks);


ylabel(t,'$TK$~[frames]',MajorProps{:})
xlabel(t,'$SK$~[pixels]',MajorProps{:})

FigName=strcat(SaveDir,'/',Desig,'_BulkRMSE_NoNorm');
saveas(gcf,FigName,'fig');
saveas(gcf,FigName,'svg');
saveas(gcf,FigName,'png');

%% Save the results
fprintf('Saving Results \n')
DataSaveName=strcat(SaveDir,'/',Desig,'_RMSE_SmoothSweep.mat');
save(DataSaveName)

fprintf('COMPLETE \n')