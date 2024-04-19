%This script is written to perform initial evaluations of the the
    %minimization algorithm

close all
clear all
clc

%% Load stress gauge and strain data
SG_filename=uigetfile('.mat','Choose File Containing Stress Gage Data');

%% Get Test designation
prompt='input test designation';
deg=char(inputdlg(prompt));

%% Stress gauge and strain data
fprintf('Loading Data \n')
load(SG_filename);
testdeg=deg;
clear deg

%% Create structute of known parameters
Einf=2.98E9;
nu=0.26;

exactProps.Kinf=Einf/(3*(1-2*nu));
exactProps.Ginf=Einf/(2*(1+nu));
exactProps.nu=0;

%% create vector of initial guesses
Eint=2E9;   
intGuess=[Eint/(3*(1-2*nu));... K1
          Eint/(2*(1+nu));.... G1
          10E-6];... tau1

%% Create reasonable upper bounds
ub=[5E9;...
    3E9;...
    100E-6];
    
 
 
%% Create reasonable lower bounds
lb=[1E8;...
    1E8;...
    1E-6];

%% Create Solve Opts Structure
SolveOpts.constnu=1; %No longer constrain Poisson's ratio
SolveOpts.KGsame=true; %simultaneously identify K and G
SolveOpts.identEinf=false;
SolveOpts.identForm='KG';
SolveOpts.minFunc='fmincon';

%% options for the fmincon algorithm
minOpts.Algorithm='interior-point';
%minOpts.Algorithm='sqp';
minOpts.UseParallel=true;

%% Run the minimization
fprintf('running minimization \b')
fprintf('Poissons ratio constrained \n') 
ConstitutiveParam=func_PronySGvfmV3(strain,time.vec,SG,Shear_SG,...
                    exactProps,...
                    intGuess,ub,lb,SolveOpts,minOpts);

%% Print Results
save(strcat(testdeg,'_Cnu_KGident.mat'),'intGuess','ub','lb','ConstitutiveParam','SolveOpts');

