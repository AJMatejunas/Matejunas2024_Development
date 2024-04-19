%This script was written to validate the identification of a single maxwell
    %element for a viscoelastic material using the kinematic fields output
    %by a finite element model of an IBII experiment. 


%% Select file containing finite element kinematic fields
[FEname,FEpath]=uigetfile('*.mat', ...
    'Choose file containing finite element kinematic fields');

%% Choose folder to save results
SaveDir=uigetdir('','Choose Folder in which to save results');

%% Define Test Dseignation
TestDeg=char(cell2mat(inputdlg('Input Test Designation')));
%% Add reference parameters file 
[RefPar.file,RefPar.path]=uigetfile('*.mat',...
    'Choose Mat file containing reference constitutive pamrameters');


%% Load Finite element data
fprintf('Loading finite element kinematic fields \n')
FEfile=strcat(FEpath,'/',FEname);
load(FEfile);
fprintf('kinematic fields loaded \n')

%% load reference parameter file
fprintf('loading reference parameter file \n')
load(strcat(RefPar.path,'/',RefPar.file),'MatProps');
fprintf('Reference parameters loaded \n')

%% Define Sample Condtioning Options
CondOpts.ImpCens=20;
CondOpts.FreeCens=40;
CondOpts.Xds=3;
CondOpts.Yds=1;
CondOpts.Tds=1;
CondOpts.CutEndFrames=2;

if CondOpts.Tds <=1
    CondOpts.TempDS=false;
else
    CondOpts.TempDS=true;
end
%% Create structute of known parameters
    Einf=2.98E9;
    nu=0.26;

    exactProps.Kinf=Einf/(3*(1-2*nu));
    exactProps.Ginf=Einf/(2*(1+nu));
    exactProps.nu=0;

%% Create reasonable upper bounds
%bounds on the initial guess are +/- 25%
ub=[1.0962E9;... G1
    12.5E-6]; %tau_1
ubK=1.875e9;
%% Create reasonable lower bounds
lb=[657.7E6;... G1
    7.5E-6]; %tau_1
lbK=1.125e9;
%% Create Matrix of initial guesses
    %Eint=2.5E9;
    %Choose a random 3x3 matrix of initial guress multipliers. Multipliers
        %will be uniformly distributed random numbers in the interval [0,1]
 
    intMult=rand(2,3);
    intMultK=rand(1,3);
  
    %Create 3X3 matrix of initial guess parameters
     intMatrix=lb+(ub-lb).*intMult;
     intVecK=lbK+(ubK-lbK).*intMultK;
%     intGuess=[Eint/(3*(1-2*nu));... K1
%               Eint/(2*(1+nu));.... G1
%               10E-6];... tau1

    %% Create Solve Opts Structure
    SolveOpts.constnu=0; %No longer constrain Poisson's ratio
    SolveOpts.KGsame=false; %simultaneously identify K and G
    SolveOpts.identEinf=false;
    SolveOpts.identForm='KG';
    SolveOpts.minFunc='fmincon';

    %% options for the fmincon algorithm
    minOpts.Algorithm='interior-point';
    %minOpts.Algorithm='sqp';
    minOpts.UseParallel=true;
    minOpts.StepTolerance=1e-3;

%% Calculate Accelerations from disp
 diffOpts.method='gradient';
 fprintf('Calculating accelerations in the same manner as for grid method data \n')
 [~,accelD] = func_calcFFAccelFromDisp(time,disp,diffOpts.method);
%% Calculate stress gauge stresses
fprintf('Calculating Stress Gage Stresses using FE accelerations \n')
SG=func_Full_SG(accel,pos.x,time,material.rho);
fprintf('Calculating stress guage stresses with calculated accelerations \n')
SGD=func_Full_SG(accelD,pos.x,time,material.rho);

%% Initialize temporaty variables
[TempPhi,TempG,TempTau,TempPhiD,TempGD,TempTauD,...
    TempK,TempKD,TempPhiK,TempPhiKD]=deal(zeros([3,1]));
%% Run parameter identification
for n=1:3
    %% Shear identification
    Kguess=intVecK(n);
    exactProps.Ki=Kguess;
    intGuess=squeeze(intMatrix(:,n));
    fprintf('Running Shear Parameter Identification on FE accelerations \n')
    [TempParam,TempPhi(n)]=func_PronyShearSGvfm(strain.s, ...
        time.vec,SG.s,exactProps,...
        intGuess,ub,lb,SolveOpts,minOpts,CondOpts);
    TempG(n)=TempParam(1);
    TempTau(n)=TempParam(2);

    fprintf('Running parameter identification on calculated accelerations \n')
    [TempParamD,TempPhiD(n)]=func_PronyShearSGvfm(strain.s, ...
        time.vec,SGD.s,exactProps,...
        intGuess,ub,lb,SolveOpts,minOpts,CondOpts);
    TempGD(n)=TempParamD(1);
    TempTauD(n)=TempParamD(2);


    fprintf('Running Bulk Identification with FE accelerations \n')
    %% Bulk Identification
    %remove exactPropsK.Ki
    exactPropsK=rmfield(exactProps,'Ki');
    exactPropsK.Gi=TempG(n);
    exactPropsK.tau=TempTau(n);
    [TempK(n),TempPhiK(n)]=func_PronyBulkSGvfm(strain,...
        time.vec,SG.x,...
        exactPropsK,Kguess,ubK,lbK,SolveOpts,minOpts,CondOpts);

    fprintf('Running bulk parameter identification with calculated accelerations \n')
    exactPropsK.Gi=TempGD(n);
    exactPropsK.tau=TempTauD(n);
    [TempKD(n),TempPhiKD(n)]=func_PronyBulkSGvfm(strain,...
        time.vec,SGD.x,...
        exactPropsK,Kguess,ubK,lbK,SolveOpts,minOpts,CondOpts);
    
end

    minPhi=min(TempPhi);
    minPhiD=min(TempPhiD);
    minPhiK=min(TempPhiK);
    minPhiKD=min(TempPhiKD);
    Kident=TempK(TempPhiK==minPhiK)
    Gident=TempG(TempPhi==minPhi)
    Tauident=TempTau(TempPhi==minPhi)
    KDident=TempKD(TempPhiKD==minPhiKD)
    GDident=TempGD(TempPhiD==minPhiD)
    TauDident=TempTauD(TempPhiD==minPhiD)
    phi=TempPhi;
    phiD=TempPhiD;

    IdentParam=[Kident;Gident;Tauident];
    IdentParamD=[KDident,GDident,TauDident];

%% Calculate Errors
Errors.K=(Kident-MatProps.Ki)/MatProps.Ki*100;
Errors.G=(Gident-MatProps.Gi)/MatProps.Gi*100;
Errors.tau=(Tauident-MatProps.tau)/MatProps.tau*100;
Errors.KD=(KDident-MatProps.Ki)/MatProps.Ki*100;
Errors.GD=(GDident-MatProps.Gi)/MatProps.Gi*100;
Errors.tauD=(TauDident-MatProps.tau)/MatProps.tau*100

%% Save data
SaveName=strcat(SaveDir,'/',TestDeg,'_FE_SeparateKGIdent.mat');
save(SaveName,'accel','accelD','SG','SGD','pos','MatProps','strain', ...
    'IdentParam','IdentParamD','stress','Errors');