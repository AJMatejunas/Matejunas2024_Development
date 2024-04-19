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
load(SG_filename);
testdeg=deg;
clear deg

%% Create structute of known parameters
fprintf('Saving known parameters \n')
Einf=2.98E9;
nu=0.26;

exactProps.Kinf=Einf/(3*(1-2*nu));
exactProps.Ginf=Einf/(2*(1+nu));
exactProps.K1=0.7673e9;
exactProps.K2=exactProps.K1;

exactProps.G1=0.4385e9;
exactProps.G2=exactProps.G1;

exactProps.nu=0;

%% create vector of initial guesses
fprintf('Generating initial guesses and bounds \n')
%Randomly assign a K guess
Kint=[10*rand(1)*exactProps.K1,10*rand(1)*exactProps.K2]; 
%Randomly assign G guess
Gint=[10*rand(1)*exactProps.G1,10*rand(1)*exactProps.G2];   
%

intGuess=[Kint;... K_i
          Gint;.... G_i
          50E-6,500E-6];... tau_i guesses are placed in the middle of the 
                            %order of magnitudes that are to be identified

%% Create reasonable upper bounds
ub=[5E9,5e9;...
    3E9,3e9;...
    100E-6,1000e-6];
    
 
 
%% Create reasonable lower bounds
lb=[1E6,1e6;...
    1E6,1e6;...
    10E-6,100e-6];

%% Create Solve Opts Structure
fprintf('creating solve options \n')
SolveOpts.constnu=0; %No longer constrain Poisson's ratio
SolveOpts.KGsame=false; %simultaneously identify K and G
SolveOpts.identEinf=false;
SolveOpts.identForm='KG';
SolveOpts.minFunc='fmincon';

%% options for the fmincon algorithm
minOpts.Algorithm='interior-point';
%minOpts.Algorithm='sqp';
minOpts.UseParallel=true;
%% Run the minimization
fprintf('Running the solution in parallel \n')
ConstitutiveParam=func_PronySGvfmV4(strain,time.vec,SG,Shear_SG,...
                    exactProps,...
                    intGuess,ub,lb,SolveOpts,minOpts);

%% Print Results
save(strcat(testdeg,'_KGident.mat'),'intGuess','ub','lb','ConstitutiveParam','SolveOpts');

