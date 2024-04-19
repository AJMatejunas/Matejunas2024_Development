%This script is written to perform initial evaluations of the the
    %minimization algorithm

close all
clear all
clc

%% Load stress gauge and strain data
SG_filename=uigetfile('.mat','Choose File Containing Stress Gage Data');


%% Get Test designation
prompt='input test designation';
desig=char(inputdlg(prompt));

fprintf('Loading Data file \n')
load(SG_filename);
testdeg=desig;
clear desig
%% Create structute of known parameters
Einf=2.98E9;
nu=0.26;

exactProps.Kinf=Einf/(3*(1-2*nu));
exactProps.Ginf=Einf/(2*(1+nu));
exactProps.nu=nu;

%% create vector of initial guesses
%Eint=2.21E9; 
Gint=0.75E9;
intGuess=[2*Gint*(1+nu)/(3*(1-2*nu));... K1
         Gint;.... G1
        12E-6];... tau1

%% Create reasonable upper bounds
ub=[11E9;...
    10E9;...
    1E-3];
    
  
 
%% Create reasonable lower bounds
lb=[0;...
    0;...
    1E-9];

%% Create Solve Opts Structure
SolveOpts.constnu=1;
SolveOpts.KGsame=true;
SolveOpts.identEinf=false;
SolveOpts.identForm='KG';
SolveOpts.minFunc='fmincon';

%% options for the fmincon algorithm
minOpts.Algorithm='interior-point';
%minOpts.Algorithm='sqp';
minOpts.UseParallel=true;

%% Run the minimization
fprintf('running minimization \n')
ConstitutiveParam=func_PronyXYSGvfm(strain,time.vec,Shear_SG,exactProps,...
    intGuess,ub,lb,SolveOpts,minOpts);

%% Print Results
fprintf('Identification finished, saving results \n')
save(strcat(testdeg,'_ident.mat'),'intGuess','ub','lb','ConstitutiveParam');
fprintf('done \n')