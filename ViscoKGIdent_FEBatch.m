%This script is written to batch extract constitutive parameters from finite
    %element data from a parametric sweep
close all
clear all
clc


%% Create sweep variables
SweepDesig=char(inputdlg('Sweep Designation'));
MainDir=uigetdir('Choose Home Directory');
[DesigVec,SaveDirVec,FileNameVec]=...
func_BatchProcessVars(SweepDesig,MainDir);

%% Create Data Conditioning Options
CondOpts.ImpCens=10;
CondOpts.FreeCens=10;
CondOpts.Xds=1;
CondOpts.Yds=1;
CondOpts.TempDS=false;

%% Create Progress Bar
TotIter=length(FileNameVec);
ItNum=1;
prognum=ItNum/TotIter;
progmsg2=strcat('0/',num2str(TotIter),'Minimizations Complete');
Progress2=waitbar(prognum,progmsg2);

%% Run Sweep
for mm=1:length(FileNameVec)
    
    Desig=DesigVec{mm};
    FileName=FileNameVec{mm};  
    prognum=mm/TotIter;
    
  %% Load Data
    fprintf(strcat('Loading ',Desig,'\n'))
        load(FileName,'-regexp', '^(?!Progress)\w');
    
    Desig=DesigVec{mm};  

  %%
    TotIter=length(FileNameVec);
    ItStr=num2str(mm);
    SaveDir=SaveDirVec{mm};
    %% Generate Stress Gauge Data
    fprintf(strcat('Calculating Stress Gauge for ',Desig,'\n'))
    [X_vec,Y_vec,Full_SG]=...
    func_FEViscoSGCalc(Desig,SaveDir,Xq,Yq,accel,strain,stress,time,...
    material);
  
    %% Create structute of known parameters
    Einf=2.98E9;
    nu=0.26;

    exactProps.Kinf=Einf/(3*(1-2*nu));
    exactProps.Ginf=Einf/(2*(1+nu));
    exactProps.nu=0;

    %% create vector of initial guesses
    Eint=2.5E9;
    intGuess=[1E9;... K1
       1E8;.... G1
        6.785e-6];... tau1

    %% Create reasonable upper bounds
    ub=[5E9;...
        3E9;...
        1000E-3];



    %% Create reasonable lower bounds
    lb=[1E6;...
        1E6;...
        1E-9];

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
    
    
    %% Define Censorship options
    fprintf('Censoring Impact and Free Edges \n')
      
    [Full_SG,accel,strain,X_vec]=func_conditionIBIIData(Full_SG,accel,...
    strain,X_vec,time,CondOpts);
   
    SG=Full_SG.x;
    Shear_SG=Full_SG.s;

    %% Run the minimization
    fprintf(strcat('Running Minimization \n'))
    ConstitutiveParam=func_PronySGvfmV5(strain,time.vec,SG,Shear_SG,...
        exactProps,...
        intGuess,ub,lb,SolveOpts,minOpts,CondOpts);

    %% Print Results
    mkdir(SaveDir)
    save(strcat(SaveDir,Desig,'_KGident.mat'),'intGuess','ub','lb',...
        'ConstitutiveParam','SolveOpts');

    progmsg2=strcat(ItStr,'/',num2str(TotIter),'Minimizations Complete');
    Progress2=waitbar(prognum,Progress2,progmsg2);
    
    %% Ready Workspace for next iteration
    clearvars -except SweepDesig MainDir DesigVec SaveDirVec FileNameVec...
        progmsg2 Progress2 TotIter CondOpts

end
