% This code is written to check the FE cost fucntion for a variety of
    % initial constitutive parameter guesses

%% Initialize
clear variables; close all; clc

%% Load FE data
[FE.name,FE.path]=uigetfile('*.mat','Choose File containing FE fields');
FE.file=strcat(FE.path,'/',FE.name);

fprintf('Loading finite element fields \n')
load(FE.file);

fprintf('Kinematic fields loaded \n')

%% Choose directory in which to save files
MainSaveDir=uigetdir('','Choose directory to save results');

%% Define Parent Desig
ParentDesig=char(inputdlg('Test Designation'));

%% Choose constitutive parameter input file
[RefPar.Name,RefPar.Path]=uigetfile('*.mat',...
    'Choose file containing reference parameters');
RefPar.File=strcat(RefPar.Path,'/',RefPar.Name);
RefPar=load(RefPar.File,'MatProps');

%% Create Vectors of Material Parameters
%Upper and lower limit of bulk modulus.
Kexact=1.5347e9;
Kident=1.4273e9;
Kul=RefPar.MatProps.Ki*1.5;
Kll=RefPar.MatProps.Ki*.5;
%Vector of K to plot cost function for
K_vec=linspace(Kll,Kul,201);
%K_vec=[Kll,Kident,Kexact,Kul];
%Identified shear modulus with no smoothing
Gident=8.8737e+08;

%Exact shear modulus input
Gexact=8.7698e8;

%Upper and lower limits for shear modulus
Gul=RefPar.MatProps.Gi*1.5;
Gll=RefPar.MatProps.Gi*.5;
G_vec=linspace(Gll,Gul,200);
%Test vector for code testing/debugging pourposes;
%G_vec=[Gll,Gexact,Gul];


%Upper and lower limits for tau
tau_ul=10e-6*1.25;
tau_ll=10e-6*.75;

tau_vec=linspace(tau_ll,tau_ul,100);
tau_ident=9.8378e-6;
tau_exact=10e-6;

%% Set sample condtioning options
CondOpts.ImpCens=27;
CondOpts.FreeCens=36;
CondOpts.Xds=5;
CondOpts.Yds=1;
CondOpts.Tds=1;

if CondOpts.Tds <=1
    CondOpts.TempDS=false;
else
    CondOpts.TempDS=true;
end
   
%% Condition data
fprintf('Pre Conditioning data \n')
clear SG
SG=Full_SG;

[SG,accel,strain,X_vec]=func_conditionIBIIData(SG,accel,...
                strain,X_vec,time,CondOpts);

%% Calculate and plot cost function for exact tau
fprintf('Plotting Cost Function \n')
EtauDesig=strcat(ParentDesig,'_ExactTau');
knownParam=zeros(1,4);
knownParam(1)=RefPar.MatProps.Kinf;
knownParam(2)=RefPar.MatProps.Ginf;
 [Etau.Phi,Etau.Phi_K,Etau.Phi_G] = func_PlotPhiConstTau(K_vec,G_vec, ...
     tau_exact,...
    SG,knownParam,strain,time,CondOpts,EtauDesig,MainSaveDir, ...
    Kexact,Gexact);

%% Calculate Constituitve model stresses
fprintf('Calculating Constituitve Model Stresses \n')
StressModel=func_ViscoConstitutiveV6(strain.x,strain.y,strain.s,time.vec,...
    RefPar.MatProps,0,0,0);
fprintf('Model stresses calculated \n')

%% save
fprintf('Saving Workspace \n')
SaveName=strcat(MainSaveDir,'/',ParentDesig,'_FE_CostFunctions.mat');
save(SaveName)

fprintf('Complete \n')
