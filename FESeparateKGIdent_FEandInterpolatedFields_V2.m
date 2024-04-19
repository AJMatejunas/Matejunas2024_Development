% This script was created to run multiple identifications on finite element
    %data with different strain and acceleration inputs including
        %Accelerations directly from FE output
        %Accelerations calculated from FE displacements
        %Accelerations interpolated from FE output to  grid coordinates
        %Accelerations calculated from interpolated FE displacements
        %Strains output directly from the FE simulation
        %Strains interpolated onto the grid coordinates


 %Author: Andrew Matejunas

 %Date Created: 2023-04-22

 %Version History/Changelog
    %2023-04-24: No version change. Increased number of initial guesses to
        %10
        %Also changed step tolerance to 1e-4 for more robust identification
    %2023-04-28-V2: Changed number of initial guesses to 30
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize

clear variables; close all; clc

%% Choose folder to save results
SaveDir=uigetdir('','Choose folder in which to save results');

%% Define simulation designation
SimDesig=char(cell2mat(inputdlg('Input Simulation Designation')));

%% Choose files containing raw FE and Interpolated FE kinematic fields
[FEname,FEpath]=uigetfile('*.mat', ...
    'Choose file containing finite element kinematic fields');
FEfile=strcat(FEpath,'/',FEname);

[INTname,INTpath]=uigetfile('*.mat', ...
    'Choose file containing kinematic fields interpolated onto grid coordinates');
INTfile=strcat(INTpath,'/',INTname);

%% Load File containing reference constitutive parameters
[PropsName,PropsPath]=uigetfile('*.mat',...
    'Choose File Containing Reference Parameters');
PropsFile=strcat(PropsPath,'/',PropsName);
%% Load files
fprintf('Loading Finite Element File \n')
FE=load(FEfile,'disp','time','strain','accel','pos');

fprintf('Loading Interpolated File \n')
INT=load(INTfile,'disp','time','strain','accel','pos','GridPos');

fprintf('Loading Material Properties File \n')
load(PropsFile);

%% Define Sample Condtioning Options
CondOptsFE.ImpCens=90;
CondOptsFE.FreeCens=36;
CondOptsFE.Xds=1;
CondOptsFE.Yds=1;
CondOptsFE.Tds=1;
CondOptsFE.CutEndFrames=4;

if CondOptsFE.Tds <=1
    CondOptsFE.TempDS=false;
else
    CondOptsFE.TempDS=true;
end



%% Define Sample Condtioning Options
CondOptsINT.ImpCens=50;
CondOptsINT.FreeCens=20;
CondOptsINT.Xds=1;
CondOptsINT.Yds=1;
CondOptsINT.Tds=1;
CondOptsINT.CutEndFrames=4;

if CondOptsINT.Tds <=1
    CondOptsINT.TempDS=false;
else
    CondOptsINT.TempDS=true;
end
%% Create structute of known parameters
    Einf=2.98E9;
    nu=0.26;

    exactPropsG.Kinf=Einf/(3*(1-2*nu));
    exactPropsG.Ginf=Einf/(2*(1+nu));
    exactPropsG.nu=0;
    exactPropsG.Ki=MatProps.Ki;
   
%% Create reasonable upper bounds
%bounds on the initial guess are +/- 25%
ub=[1.875e9;
    1.0962E9;... G1
    12.5E-6]; %tau_1
ubG=ub(2:3);

%% Create reasonable lower bounds
lb=[1.125e9;
    100e6;... G1
    7.5E-6]; %tau_1
lbG=lb(2:3);

%% Create Matrix of initial guesses
    %Eint=2.5E9;
    %Choose a random 3x3 matrix of initial guress multipliers. Multipliers
        %will be uniformly distributed random numbers in the interval [0,1]
 
    intMult=rand(3,30);
   
  
    %Create 3X3 matrix of initial guess parameters
     intMatrix=lb+(ub-lb).*intMult;
     intMatrix=intMatrix(2:3,:);

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
    minOpts.UseParallel=false;
    minOpts.StepTolerance=1e-4;
%% Calcuate Accelerations 
 diffOpts.method='gradient';
 fprintf('Calculating FE accelerations in the same manner as for grid method data \n')
 [~,FE.accelD] = func_calcFFAccelFromDisp(FE.time,FE.disp,diffOpts.method);

 
 fprintf('Calculating Accelerations on Interpolated Kinematic Fields \n')
[~,INT.accelD] = func_calcFFAccelFromDisp(INT.time,INT.disp,diffOpts.method);

 
%% Calculate stress gauge stresses
%Finite Element Coordinate
fprintf('Calculating stress gauge stresses for finite element accelerations \n')
FE.SG=func_Full_SG(FE.accel,FE.pos.x,FE.time,MatProps.rho);

fprintf('Calculating SG stresses for accelerations calculated from finite element displacements \n')
FE.SGd=func_Full_SG(FE.accelD,FE.pos.x,FE.time,MatProps.rho);

%Interpolated to Grid Coordinates
fprintf('Calculating SG for Interpolated FE accelerations \n')
INT.SG=func_Full_SG(INT.accel,INT.pos.x,INT.time,MatProps.rho);

fprintf('Calculating SG for Accelerations calculated from interpoplated accelerations \n')
INT.SGd=func_Full_SG(INT.accelD,INT.pos.x,INT.time,MatProps.rho);


%% Record original kinematic fields
TimeFull=FE.time;
FE.strainFull=FE.strain;
FE.SGfull=FE.SG;
FE.X_VecFull=FE.pos.x;

INT.strainFull=INT.strain;
INT.SGfull=INT.SG;
INT.X_VecFull=INT.pos.x;

%% Condition Fields

fprintf('Conditioning FE kinematic Fields \n')
[FE.SG,FE.accel,FE.strain,FE.X_vec,FE.time]=func_conditionIBIIDataV2(...
    FE.SGfull,FE.accel,FE.strainFull,FE.pos.x,TimeFull,CondOptsFE); 

fprintf('Conditioning FE kinematic fields with calculated accelerations \n')
[FE.SGd,FE.accelD,~,~,~]=func_conditionIBIIDataV2(...
    FE.SGd,FE.accelD,FE.strainFull,FE.pos.x,TimeFull,CondOptsFE); 

fprintf('Conditioning interpolated fields \n')
[INT.SG,INT.accel,INT.strain,INT.X_vec,INT.time]=func_conditionIBIIDataV2(...
    INT.SGfull,INT.accel,INT.strainFull,INT.pos.x,TimeFull,CondOptsINT);

fprintf('Conditioning interpolated fields with calculated accelerations \n')
[INT.SGd,INT.accelD,~,~,~]=func_conditionIBIIDataV2(...
    INT.SGd,INT.accelD,INT.strainFull,INT.pos.x,TimeFull,CondOptsINT); 

%% Initialize identification variables
[FETempG,FETempPhiG,FETempGd,FETempPhiGd,...
    INTTempG,INTTempPhiG,INTTempGd,...
    INTTempPhiGd]=deal(zeros([30,1]));

%% Set up Par For loops
FEshear=FE.strain.s;
FEtime=FE.time.vec;
FESGs=FE.SG.s;
FESGds=FE.SGd.s;

INTshear=INT.strain.s;
INTtime=INT.time.vec;
INTSGs=INT.SG.s;
INTSGds=INT.SGd.s;
%% Run Shear Identifications
parfor ig=1:30
    intGuess=squeeze(intMatrix(:,ig));

    fprintf('Running Shear identification of raw FE data \n')
    [TempParam,TempPhi]=func_PronyShearSGvfm(FEshear, ...
        FEtime,FESGs,exactPropsG,...
        intGuess,ubG,lbG,SolveOpts,minOpts,CondOptsFE);
    FETempG(ig)=TempParam(1);
    FETempPhiG(ig)=TempPhi;

    fprintf('Running shear identfication on FE data with calculated accelerations \n')
    [TempParam,TempPhi]=func_PronyShearSGvfm(FEshear, ...
        FEtime,FESGds,exactPropsG,...
        intGuess,ubG,lbG,SolveOpts,minOpts,CondOptsFE);
    FETempGd(ig)=TempParam(1);
    FETempPhiGd(ig)=TempPhi;

    fprintf('Running Shear Identification on Interpolated FE data \n')
    [TempParam,TempPhi]=func_PronyShearSGvfm(INTshear, ...
        INTtime,INTSGs,exactPropsG,...
        intGuess,ubG,lbG,SolveOpts,minOpts,CondOptsINT);
   INTTempG(ig)=TempParam(1);
   INTTempPhiG(ig)=TempPhi;

    fprintf('Running shear identfication on interpolated FE data with calculated accelerations \n')
    [TempParam,TempPhi]=func_PronyShearSGvfm(INTshear, ...
        INTtime,INTSGds,exactPropsG,...
        intGuess,ubG,lbG,SolveOpts,minOpts,CondOptsINT);
    INTTempGd(ig)=TempParam(1);
    INTTempPhiGd(ig)=TempPhi;

end

%% Put Vectors into data structures
FE.Ident.TempG=FETempG;
FE.TempPhiG=FETempPhiG;
FE.Ident.TempGd=FETempGd;
FE.TempPhiGd=FETempPhiGd;
clear FETempG FETempPhiG FETempGd FETempPhiGd FEshear FESGs FESGds

INT.Ident.TempG=INTTempG;
INT.TempPhiG=INTTempPhiG;
INT.Ident.TempGd=INTTempGd;
INT.TempPhiGd=INTTempPhiGd;
clear INTTempG INTTempPhiG INTTempGd INTTempPhiGd

%% Pull out identified parameters
% Raw Finite Element Kinematic Fields
minPhi=min(FE.TempPhiG);
FE.Ident.G=FE.Ident.TempG(FE.TempPhiG==minPhi);
FE.Ident.G=FE.Ident.G(1);

%FE with Accelerations Calculated From Displacements
minPhi=min(FE.TempPhiGd);
FE.Ident.Gd=FE.Ident.TempGd(FE.TempPhiGd==minPhi);
FE.Ident.Gd=FE.Ident.Gd(1);

%FE interpolated onto grid coordinates
minPhi=min(INT.TempPhiG);
INT.Ident.G=INT.Ident.TempG(INT.TempPhiG==minPhi);
INT.Ident.G=INT.Ident.G(1);

%Accelerations calculated from interpolated displacements
minPhi=min(INT.TempPhiGd);
INT.Ident.Gd=INT.Ident.TempGd(INT.TempPhiGd==minPhi);
INT.Ident.Gd=INT.Ident.Gd(1);
%% Calculate Errors
FE.Error.G=(FE.Ident.G-MatProps.Gi)/MatProps.Gi*100;
FE.Error.Gd=(FE.Ident.Gd-MatProps.Gi)/MatProps.Gi*100;

INT.Error.G=(INT.Ident.G-MatProps.Gi)/MatProps.Gi*100;
INT.Error.Gd=(INT.Ident.Gd-MatProps.Gi)/MatProps.Gi*100;

%% Calculate Error statistics
FE.Ident.Gmean=mean(FE.Ident.TempG);
FE.Ident.Gstd=std(FE.Ident.TempG);
FE.Ident.Gdmean=mean(FE.Ident.TempGd);
FE.Ident.Gdstd=std(FE.Ident.TempGd);

FE.Error.TempG=(FE.Ident.TempG-MatProps.Gi)/MatProps.Gi*100;
FE.Error.TempGd=(FE.Ident.TempGd-MatProps.Gi)/MatProps.Gi*100;
FE.Error.Gmean=mean(FE.Error.TempG);
FE.Error.Gdmean=mean(FE.Error.TempGd);
FE.Error.Gstd=std(FE.Error.TempG);
FE.Error.Gdstd=mean(FE.Error.TempGd);

%Interpolated fields
INT.Ident.Gmean=mean(INT.Ident.TempG);
INT.Ident.Gstd=std(INT.Ident.TempG);
INT.Ident.Gdmean=mean(INT.Ident.TempGd);
INT.Ident.Gdstd=std(INT.Ident.TempGd);

INT.Error.TempG=(INT.Ident.TempG-MatProps.Gi)/MatProps.Gi*100;
INT.Error.TempGd=(INT.Ident.TempGd-MatProps.Gi)/MatProps.Gi*100;
INT.Error.Gmean=mean(INT.Error.TempG);
INT.Error.Gdmean=mean(INT.Error.TempGd);
INT.Error.Gstd=std(INT.Error.TempG);
INT.Error.Gdstd=mean(INT.Error.TempGd);
%% Set up variables for K and Tau identification
exactPropsK=exactPropsG;
exactPropsK=rmfield(exactPropsK,'Ki');

[exactPropsKFE,exactPropsKFEd,exactPropsKINT,exactPropsKINTd]=...
    deal(exactPropsK);

exactPropsKFE.Gi=FE.Ident.G;
exactPropsKFEd.Gi=FE.Ident.Gd;
exactPropsKINT.Gi=INT.Ident.G;
exactPropsKINTd.Gi=INT.Ident.Gd;

%% Create reasonable upper bounds
%bounds on the initial guess are +/- 25%
ubK=[1.875e9;... K1
    12.5E-6]; %tau_1


%% Create reasonable lower bounds
lbK=[1.125e9;...
    7.5E-6]; %tau_1

%% Create Matrix of initial guesses
    %Eint=2.5E9;
    %Choose a random 3x3 matrix of initial guress multipliers. Multipliers
        %will be uniformly distributed random numbers in the interval [0,1]
 
    intMultK=rand(2,30);
      
    %Create 3X3 matrix of initial guess parameters
     intMatrixK=lbK+(ubK-lbK).*intMultK;
   
%     intGuess=[Eint/(3*(1-2*nu));... K1
%               Eint/(2*(1+nu));.... G1
%               10E-6];... tau1
   
%% Initialize bulk identification variables
[FETempK,FETempPhiK,FETempKd,FETempPhiKd,...
    INTTempK,INTTempPhiK,INTTempKd,...
    INTTempPhiKd,...
    FETempTau,FETempTaud,...
    INTTempTau,INTTempTaud]=deal(zeros([10,1]));

%%
FEstrain=FE.strain;
FESGx=FE.SG.x;
FESGdx=FE.SGd.x;

INTstrain=INT.strain;
INTSGx=INT.SG.x;
INTSGdx=INT.SGd.x;

%% Run Bulk Identifications
parfor ig=1:30
    intGuess=squeeze(intMatrixK(:,ig));

    fprintf('Running Bulk identification of raw FE data \n')
    [TempParam,TempPhi]=func_PronyBulkTauSGvfm(FEstrain, ...
        FEtime,FESGx,exactPropsKFE,...
        intGuess,ubK,lbK,SolveOpts,minOpts,CondOptsFE);
    FETempK(ig)=TempParam(1);
    FETempTau(ig)=TempParam(2);
    FETempPhiK(ig)=TempPhi;

    fprintf('Running sBulk identfication on FE data with calculated accelerations \n')
    [TempParam,TempPhi]=func_PronyBulkTauSGvfm(FEstrain, ...
        FEtime,FESGdx,exactPropsKFEd,...
        intGuess,ubK,lbK,SolveOpts,minOpts,CondOptsFE);
    FETempKd(ig)=TempParam(1);
    FETempTaud(ig)=TempParam(2);
    FETempPhiKd(ig)=TempPhi;
    

    fprintf('Running Bulk Identification on Interpolated FE data \n')
    [TempParam,TempPhi]=func_PronyBulkTauSGvfm(INTstrain, ...
        INTtime,INTSGx,exactPropsKINT,...
        intGuess,ubK,lbK,SolveOpts,minOpts,CondOptsINT);
   INTTempK(ig)=TempParam(1);
   INTTempTau(ig)=TempParam(2);
   INTTempPhiK(ig)=TempPhi;

    fprintf('Running Bulk identfication on interpolated FE data with calculated accelerations \n')
    [TempParam,TempPhi]=func_PronyBulkTauSGvfm(INTstrain, ...
        INTtime,INTSGdx,exactPropsKINTd,...
        intGuess,ubK,lbK,SolveOpts,minOpts,CondOptsINT);
   INTTempKd(ig)=TempParam(1);
   INTTempTaud(ig)=TempParam(2);
   INTTempPhiKd(ig)=TempPhi;
end

%% Move into structures
FE.Ident.TempK=FETempK;
FE.TempPhiK=FETempPhiK;
FE.Ident.TempKd=FETempKd;
FE.TempPhiKd=FETempPhiKd;
clear FETempK FETempPhiK FETempKd FETempPhiKd FEshear FESGs FESGds

INT.Ident.TempK=INTTempK;
INT.TempPhiK=INTTempPhiK;
INT.Ident.TempKd=INTTempKd;
INT.TempPhiKd=INTTempPhiKd;
clear INTTempK INTTempPhiK INTTempKd INTTempPhiKd

FE.Ident.TempTau=FETempTau;
% FE.TempPhiTau=FETempPhiK;
FE.Ident.TempTaud=FETempTaud;
% FE.TempPhiTaud=FETempPhiKd;
clear FETempTau  FETempTaud  FEshear 

INT.Ident.TempTau=INTTempTau;
% INT.TempPhiTau=INTTempPhiTau;
INT.Ident.TempTaud=INTTempTaud;
% INT.TempPhiTaud=INTTempPhiTaud;
clear INTTempTau  INTTempTaud 

%% Extract
% Raw Finite Element Kinematic Fields
minPhi=min(FE.TempPhiK);
FE.Ident.K=FE.Ident.TempK(FE.TempPhiK==minPhi);
FE.Ident.K=FE.Ident.K(1);
FE.Ident.tau=FE.Ident.TempTau(FE.TempPhiK==minPhi);
FE.Ident.tau=FE.Ident.tau(1);

%FE with Accelerations Calculated From Displacements
minPhi=min(FE.TempPhiKd);
FE.Ident.Kd=FE.Ident.TempKd(FE.TempPhiKd==minPhi);
FE.Ident.Kd=FE.Ident.Kd(1);
FE.Ident.taud=FE.Ident.TempTaud(FE.TempPhiKd==minPhi);
FE.Ident.taud=FE.Ident.taud(1);

%FE interpolated onto grid coordinates
minPhi=min(INT.TempPhiK);
INT.Ident.K=INT.Ident.TempK(INT.TempPhiK==minPhi);
INT.Ident.K=INT.Ident.K(1);
INT.Ident.tau=INT.Ident.TempTau(INT.TempPhiK==minPhi);
INT.Ident.tau=INT.Ident.tau(1);

%Accelerations calculated from interpolated displacements
minPhi=min(INT.TempPhiKd);
INT.Ident.Kd=INT.Ident.TempKd(INT.TempPhiKd==minPhi);
INT.Ident.Kd=INT.Ident.Kd(1);
INT.Ident.taud=INT.Ident.TempTaud(INT.TempPhiKd==minPhi);
INT.Ident.taud=INT.Ident.taud(1);

%% Calculate bulk and tau errors
FE.Error.K=(FE.Ident.K-MatProps.Ki)/MatProps.Ki*100;
FE.Error.Kd=(FE.Ident.Kd-MatProps.Ki)/MatProps.Ki*100;
FE.Error.tau=(FE.Ident.tau-MatProps.tau)/MatProps.tau*100;
FE.Error.taud=(FE.Ident.taud-MatProps.tau)/MatProps.tau*100;

INT.Error.K=(INT.Ident.K-MatProps.Ki)/MatProps.Ki*100;
INT.Error.Kd=(INT.Ident.Kd-MatProps.Ki)/MatProps.Ki*100;
INT.Error.tau=(INT.Ident.tau-MatProps.tau)/MatProps.tau*100;
INT.Error.taud=(INT.Ident.taud-MatProps.tau)/MatProps.tau*100;

%% Calculate Error statistics
FE.Ident.Kmean=mean(FE.Ident.TempK);
FE.Ident.Kstd=std(FE.Ident.TempK);
FE.Ident.Kdmean=mean(FE.Ident.TempKd);
FE.Ident.Kdstd=std(FE.Ident.TempKd);

FE.Error.TempK=(FE.Ident.TempK-MatProps.Ki)/MatProps.Ki*100;
FE.Error.TempKd=(FE.Ident.TempKd-MatProps.Ki)/MatProps.Ki*100;
FE.Error.Kmean=mean(FE.Error.TempK);
FE.Error.Kdmean=mean(FE.Error.TempKd);
FE.Error.Kstd=std(FE.Error.TempK);
FE.Error.Kdstd=std(FE.Error.TempKd);

FE.Ident.taumean=mean(FE.Ident.TempTau);
FE.Ident.taustd=std(FE.Ident.TempTau);
FE.Ident.taudmean=mean(FE.Ident.TempTaud);
FE.Ident.taudstd=std(FE.Ident.TempTaud);

FE.Error.TempTau=(FE.Ident.TempTau-MatProps.tau)/MatProps.tau*100;
FE.Error.TempTaud=(FE.Ident.TempTaud-MatProps.tau)/MatProps.tau*100;
FE.Error.taumean=mean(FE.Error.TempTau);
FE.Error.taudmean=mean(FE.Error.TempTaud);
FE.Error.taustd=std(FE.Error.TempTau);
FE.Error.taudstd=std(FE.Error.TempTaud);

%Interpolated fields
INT.Ident.Kmean=mean(INT.Ident.TempK);
INT.Ident.Kstd=std(INT.Ident.TempK);
INT.Ident.Kdmean=mean(INT.Ident.TempKd);
INT.Ident.Kdstd=std(INT.Ident.TempKd);

INT.Error.TempK=(INT.Ident.TempK-MatProps.Ki)/MatProps.Ki*100;
INT.Error.TempKd=(INT.Ident.TempKd-MatProps.Ki)/MatProps.Ki*100;
INT.Error.Kmean=mean(INT.Error.TempK);
INT.Error.Kdmean=mean(INT.Error.TempKd);
INT.Error.Kstd=std(INT.Error.TempK);
INT.Error.Kdstd=std(INT.Error.TempKd);

INT.Ident.taumean=mean(INT.Ident.TempTau);
INT.Ident.taustd=std(INT.Ident.TempTau);
INT.Ident.taudmean=mean(INT.Ident.TempTaud);
INT.Ident.taudstd=std(INT.Ident.TempTaud);

INT.Error.TempTau=(INT.Ident.TempTau-MatProps.tau)/MatProps.tau*100;
INT.Error.TempTaud=(INT.Ident.TempTaud-MatProps.tau)/MatProps.tau*100;
INT.Error.taumean=mean(INT.Error.TempTau);
INT.Error.taudmean=mean(INT.Error.TempTaud);
INT.Error.taustd=std(INT.Error.TempTau);
INT.Error.taudstd=std(INT.Error.TempTaud);
%% Save Results
SaveName=strcat(SaveDir,'/',SimDesig,'_FE_INT_IdentValidation30runs.mat');
save(SaveName);
