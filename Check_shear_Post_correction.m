%This code is written to investigate the effects of correcting shear
%strains on the constitutive parameter identification for viscoelastic
%image based inertial impact tests

%Author: Andrew Matejunas

%Date: 2022-10-24

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% initialize
close all; clear variables; clc

%% Choose directory in which to save files
MainSaveDir=uigetdir('','Choose directory to save results');

%% Define Parent Desig
ParentDesig=char(inputdlg('Test Designation'));

%% Identify Kinematic Fields
%exact inputs shear without correction
[ENC.name,ENC.path]=uigetfile('*.mat', ...
    'Choose file containing kin fields with exact inputs without shear correction');
ENC.file=strcat(ENC.path,'/',ENC.name);

%exact input with shear correction
[ESC.name,ESC.path]=uigetfile('*.mat', ...
    'Choose file containing kin fields with exact inputs with shear correction');
ESC.file=strcat(ESC.path,'/',ESC.name);

%identified inputs No Correction
[INC.name,INC.path]=uigetfile('*.mat', ...
    'Choose file containing kin fields with identified inputs without shear correction');
INC.file=strcat(INC.path,'/',INC.name);

%identified inputs with shear correction
[ISC.name,ISC.path]=uigetfile('*.mat', ...
    'Choose file containing kin fields with exact inputs with shear correction');
ISC.file=strcat(ISC.path,'/',ISC.name);

%FE fields
[FE.name,FE.path]=uigetfile('*.mat','Choose File containing FE fields');
FE.file=strcat(FE.path,'/',FE.name);


%% Load kinematic field files
fprintf('Loading Kinematic fields with exact inputs without correction \n')
ENC=load(ENC.file);

fprintf('Loading Fields with exact inputs and shear correction \n')
ESC=load(ESC.file);

fprintf('Loading fields with identified inputs without correction \n')
INC=load(INC.file);

fprintf('Loading fields with identified inputs and shear correction \n')
ISC=load(ISC.file);

fprintf('Loading finite element fields \n')
FE=load(FE.file);

fprintf('Kinematic fields loaded \n')

%% Calculate average shear strain curves
ENC.strain.sAvg=squeeze(mean(ENC.strain.s));
ESC.strain.sAvg=squeeze(mean(ESC.strain.s));

INC.strain.sAvg=squeeze(mean(INC.strain.s));
ISC.strain.sAvg=squeeze(mean(ISC.strain.s));

FE.strain.sAvg=squeeze(mean(FE.strain.s));
FE.strain.xAvg=squeeze(mean(FE.strain.x));
FE.strain.yAvg=squeeze(mean(FE.strain.y));

%% Get indexes for 1-D plots
% Set Position Indexes for Grid image data
    GridPos.NumPoints=size(ENC.disp.x,2);
    GridPos.Free=20; %4 pitches from free edge
    GridPos.Mid=round(GridPos.NumPoints/2);
    GridPos.Imp=GridPos.NumPoints-20; %4 pitches from impact edge
    % Get cooridinates of position indexes
    GridPos.Xfree=ENC.pos.x(GridPos.Free);
    GridPos.Xmid=ENC.pos.x(GridPos.Mid);
    GridPos.Ximp=ENC.pos.x(GridPos.Imp);

    %Set Position indexes and coordintes for finite element data
    FEPos.NumPoints=size(FE.disp.x,2);
    FEPos.GridRat=FEPos.NumPoints/GridPos.NumPoints;
    FEPos.Free=round(GridPos.Free*FEPos.GridRat); %4 pitches from free edge
    FEPos.Mid=round(FEPos.NumPoints/2);
    FEPos.Imp=round(GridPos.Imp*FEPos.GridRat); %4 pitches from impact edge
    
    FEPos.Xfree=FE.pos.x(FEPos.Free);
    FEPos.Xmid=FE.pos.x(FEPos.Mid);
    FEPos.Ximp=FE.pos.x(FEPos.Imp);
    
  %Print Plot selections for evaluation purposes
  FEPos.freeString=num2str(FEPos.Xfree);
  GridPos.freeString=num2str(GridPos.Xfree);
  fprintf(strcat('Free Coordinate FE:',FEPos.freeString,' Grid: ',...
      GridPos.freeString,'\n'));
  
  FEPos.midString=num2str(FEPos.Xmid);
  GridPos.midString=num2str(GridPos.Xmid);
  fprintf(strcat('Mid Coordinate FE:',FEPos.midString,' Grid: ',...
      GridPos.midString,'\n'));

  FEPos.impString=num2str(FEPos.Ximp);
  GridPos.impString=num2str(GridPos.Ximp);
  fprintf(strcat('Impact Coordinate FE:',FEPos.impString,' Grid: ',...
      GridPos.impString,'\n'));
%% Plot 1-D strain time curves
time=ENC.time.vec*10^6;

figure('units','normalized','outerposition',[0 0 1 1])

%% strain in the x direction
%Free surface
subplot(3,3,1)
plot(time,FE.strain.xAvg(FEPos.Free,:),'k')
hold on
plot(time,ENC.strain.xAvg(GridPos.Free,:),'--b')
plot(time,ESC.strain.xAvg(GridPos.Free,:),':g')
plot(time,INC.strain.xAvg(GridPos.Free,:),'--r')
plot(time,ISC.strain.xAvg(GridPos.Free,:),':c')
hold off

title('4P from Free Surface')
xlabel('time (\mu s)','FontSize',14)
ylabel('$\varepsilon_{11}$','Interpreter','latex','FontSize',14)
legend('FE','Exact no Corr','Exact Shear Corr','ident No corr', ...
    'ident shear corr','Location','northwest')

%Middle 
subplot(3,3,2)
plot(time,FE.strain.xAvg(FEPos.Mid,:),'k')
hold on
plot(time,ENC.strain.xAvg(GridPos.Mid,:),'--b')
plot(time,ESC.strain.xAvg(GridPos.Mid,:),':g')
plot(time,INC.strain.xAvg(GridPos.Mid,:),'--r')
plot(time,ISC.strain.xAvg(GridPos.Mid,:),':c')
hold off

title('Specimen Middle')
xlabel('time (\mu s)','FontSize',14)
ylabel('$\varepsilon_{11}$','Interpreter','latex','FontSize',14)

%Impact Edge
subplot(3,3,3)
plot(time,FE.strain.xAvg(FEPos.Imp,:),'k')
hold on
plot(time,ENC.strain.xAvg(GridPos.Imp,:),'--b')
plot(time,ESC.strain.xAvg(GridPos.Imp,:),':g')
plot(time,INC.strain.xAvg(GridPos.Imp,:),'--r')
plot(time,ISC.strain.xAvg(GridPos.Imp,:),':c')
hold off

title('4P from Impact')
xlabel('time (\mu s)','FontSize',14)
ylabel('$\varepsilon_{11}$','Interpreter','latex','FontSize',14)


%% Normal Strain in the Y direction
%Free surface
subplot(3,3,4)
plot(time,FE.strain.yAvg(FEPos.Free,:),'k')
hold on
plot(time,ENC.strain.yAvg(GridPos.Free,:),'--b')
plot(time,ESC.strain.yAvg(GridPos.Free,:),':g')
plot(time,INC.strain.yAvg(GridPos.Free,:),'--r')
plot(time,ISC.strain.yAvg(GridPos.Free,:),':c')
hold off

title('4P from Free Surface')
xlabel('time (\mu s)','FontSize',14)
ylabel('$\varepsilon_{22}$','Interpreter','latex','FontSize',14)

%Middle 
subplot(3,3,5)
plot(time,FE.strain.yAvg(FEPos.Mid,:),'k')
hold on
plot(time,ENC.strain.yAvg(GridPos.Mid,:),'--b')
plot(time,ESC.strain.yAvg(GridPos.Mid,:),':g')
plot(time,INC.strain.yAvg(GridPos.Mid,:),'--r')
plot(time,ISC.strain.yAvg(GridPos.Mid,:),':c')
hold off

title('Specimen Middle')
xlabel('time (\mu s)','FontSize',14)
ylabel('$\varepsilon_{22}$','Interpreter','latex','FontSize',14)


%Impact Edge
subplot(3,3,6)
plot(time,FE.strain.yAvg(FEPos.Imp,:),'k')
hold on
plot(time,ENC.strain.yAvg(GridPos.Imp,:),'--b')
plot(time,ESC.strain.yAvg(GridPos.Imp,:),':g')
plot(time,INC.strain.yAvg(GridPos.Imp,:),'--r')
plot(time,ISC.strain.yAvg(GridPos.Imp,:),':c')
hold off

title('4P from Impact')
xlabel('time (\mu s)','FontSize',14)
ylabel('$\varepsilon_{22}$','Interpreter','latex','FontSize',14)


%% Shear strains
%Free surface
subplot(3,3,7)
plot(time,FE.strain.sAvg(FEPos.Free,:),'k')
hold on
plot(time,ENC.strain.sAvg(GridPos.Free,:),'--b')
plot(time,ESC.strain.sAvg(GridPos.Free,:),':g')
plot(time,INC.strain.sAvg(GridPos.Free,:),'--r')
plot(time,ISC.strain.sAvg(GridPos.Free,:),':c')
hold off

title('4P from Free Surface')
xlabel('time (\mu s)','FontSize',14)
ylabel('$\varepsilon_{12}$','Interpreter','latex','FontSize',14)

%Middle 
subplot(3,3,8)
plot(time,FE.strain.sAvg(FEPos.Mid,:),'k')
hold on
plot(time,ENC.strain.sAvg(GridPos.Mid,:),'--b')
plot(time,ESC.strain.sAvg(GridPos.Mid,:),':g')
plot(time,INC.strain.sAvg(GridPos.Mid,:),'--r')
plot(time,ISC.strain.sAvg(GridPos.Mid,:),':c')
hold off

title('Specimen Middle')
xlabel('time (\mu s)','FontSize',14)
ylabel('$\varepsilon_{12}$','Interpreter','latex','FontSize',14)

%Impact Edge
subplot(3,3,9)
plot(time,FE.strain.sAvg(FEPos.Imp,:),'k')
hold on
plot(time,ENC.strain.sAvg(GridPos.Imp,:),'--b')
plot(time,ESC.strain.sAvg(GridPos.Imp,:),':g')
plot(time,INC.strain.sAvg(GridPos.Imp,:),'--r')
plot(time,ISC.strain.sAvg(GridPos.Imp,:),':c')
hold off

title('$P from Impact')
xlabel('time (\mu s)','FontSize',14)
ylabel('$\varepsilon_{12}$','Interpreter','latex','FontSize',14)


%% Save 1-D strain figures
SaveName1D=strcat(MainSaveDir,'/',ParentDesig,'_ShearCorr_1DstrainComp');
saveas(gcf,strcat(SaveName1D,'.fig'))
saveas(gcf,strcat(SaveName1D,'.svg'))
saveas(gcf,strcat(SaveName1D,'.png'))

%% Plot shear only at free surface
figure('units','normalized','outerposition',[0 0 1 1])
plot(time,FE.strain.sAvg(FEPos.Free,:),'k')
hold on
plot(time,ENC.strain.sAvg(GridPos.Free,:),'--b')
plot(time,ESC.strain.sAvg(GridPos.Free,:),':g')
%plot(time,INC.strain.sAvg(GridPos.Free,:),'--r')
%plot(time,ISC.strain.sAvg(GridPos.Free,:),':c')
hold off

title('4P from Free Surface')
xlabel('time (\mu s)','FontSize',14)
ylabel('$\varepsilon_{12}$','Interpreter','latex','FontSize',14)


%% Plot Heat Maps of Strain Time

%% Calculate the constitutive model stresses for the finite element data
FE.StressModel=func_ViscoConstitutiveV6(FE.strain.x,FE.strain.y, ...
    FE.strain.s,FE.time.vec, ...
    ENC.MatProps,0,0,0);

FE=rmfield(FE,'SG');
FE.SG=FE.Full_SG;

%% Calculate average FE stresses
fprintf('Calculating average FE stresses \n')
FE.stress.xAvg=squeeze(mean(FE.stress.x));
FE.stress.yAvg=squeeze(mean(FE.stress.y));
FE.stress.sAvg=squeeze(mean(FE.stress.s));

FE.StressModel.xAvg=squeeze(mean(FE.StressModel.xx));
FE.StressModel.yAvg=squeeze(mean(FE.StressModel.yy));
FE.StressModel.sAvg=squeeze(mean(FE.StressModel.xy));

%% Check constitutive model
figure(figure('units','normalized','outerposition',[0 0 1 1]))


%X-direction free edge
subplot(3,3,1)
plot(time,FE.stress.xAvg(FEPos.Free,:),'K')
hold on
plot(time,FE.SG.x(FEPos.Free,:),'--b')
plot(time,FE.StressModel.xAvg(FEPos.Free,:),'--r')
hold off

xlabel('time (\mu{}s','FontSize',14)
ylabel('\sigma_{11} (Pa)','FontSize',14)
title('4P from Free Edge')
legend('FE','SG','Const','location','northwest')

%X-direction Middle
subplot(3,3,2)
plot(time,FE.stress.xAvg(FEPos.Mid,:),'K')
hold on
plot(time,FE.SG.x(FEPos.Mid,:),'--b')
plot(time,FE.StressModel.xAvg(FEPos.Mid,:),'--r')
hold off

xlabel('time (\mu{}s','FontSize',14)
ylabel('\sigma_{11} (Pa)','FontSize',14)
title('Specimen Middle')


%X-direction free edge
subplot(3,3,3)
plot(time,FE.stress.xAvg(FEPos.Imp,:),'K')
hold on
plot(time,FE.SG.x(FEPos.Imp,:),'--b')
plot(time,FE.StressModel.xAvg(FEPos.Imp,:),'--r')
hold off

xlabel('time (\mu{}s','FontSize',14)
ylabel('\sigma_{11} (Pa)','FontSize',14)
title('4P from Impact Edge')


%Y-direction free edge
subplot(3,3,4)
plot(time,FE.stress.yAvg(FEPos.Free,:),'K')
hold on
plot(time,FE.StressModel.yAvg(FEPos.Free,:),'--r')
hold off

xlabel('time (\mu{}s','FontSize',14)
ylabel('\sigma_{22} (Pa)','FontSize',14)
title('4P from Free Edge')
legend('FE','Const','location','northwest')

%Y-direction Middle
subplot(3,3,5)
plot(time,FE.stress.yAvg(FEPos.Mid,:),'K')
hold on
plot(time,FE.StressModel.yAvg(FEPos.Mid,:),'--r')
hold off

xlabel('time (\mu{}s','FontSize',14)
ylabel('\sigma_{22} (Pa)','FontSize',14)
title('Specimen Middle')


%Y-direction free edge
subplot(3,3,6)
plot(time,FE.stress.yAvg(FEPos.Imp,:),'K')
hold on
plot(time,FE.StressModel.yAvg(FEPos.Imp,:),'--r')
hold off

xlabel('time (\mu{}s','FontSize',14)
ylabel('\sigma_{22} (Pa)','FontSize',14)
title('4P from Impact Edge')

%Shear free edge
subplot(3,3,7)
plot(time,FE.stress.sAvg(FEPos.Free,:),'K')
hold on
plot(time,FE.SG.s(FEPos.Free,:),'--b')
plot(time,FE.StressModel.sAvg(FEPos.Free,:),'--r')
hold off

xlabel('time (\mu{}s','FontSize',14)
ylabel('\sigma_{11} (Pa)','FontSize',14)
title('4P from Free Edge')
legend('FE','SG','Const','location','northwest')

%Shear Middle
subplot(3,3,8)
plot(time,FE.stress.sAvg(FEPos.Mid,:),'K')
hold on
plot(time,FE.SG.s(FEPos.Mid,:),'--b')
plot(time,FE.StressModel.sAvg(FEPos.Mid,:),'--r')
hold off

xlabel('time (\mu{}s','FontSize',14)
ylabel('\sigma_{11} (Pa)','FontSize',14)
title('Specimen Middle')


%Shear free edge
subplot(3,3,9)
plot(time,FE.stress.sAvg(FEPos.Imp,:),'K')
hold on
plot(time,FE.SG.s(FEPos.Imp,:),'--b')
plot(time,FE.StressModel.sAvg(FEPos.Imp,:),'--r')
hold off

xlabel('time (\mu{}s','FontSize',14)
ylabel('\sigma_{11} (Pa)','FontSize',14)
title('4P from Impact Edge')

SaveName1Dstress=strcat(MainSaveDir,'/',ParentDesig,'_1DConstVer');
saveas(gcf,strcat(SaveName1Dstress,'.fig'))
saveas(gcf,strcat(SaveName1Dstress,'.svg'))
saveas(gcf,strcat(SaveName1Dstress,'.png'))

%% Change format of Average Stress Model Stresses
ENC.StressModel.xAvg=squeeze(mean(ENC.StressModel.xx));
ENC.StressModel.yAvg=squeeze(mean(ENC.StressModel.yy));
ENC.StressModel.sAvg=squeeze(mean(ENC.StressModel.xy));

ESC.StressModel.xAvg=squeeze(mean(ESC.StressModel.xx));
ESC.StressModel.yAvg=squeeze(mean(ESC.StressModel.yy));
ESC.StressModel.sAvg=squeeze(mean(ESC.StressModel.xy));

INC.StressModel.xAvg=squeeze(mean(INC.StressModel.xx));
INC.StressModel.yAvg=squeeze(mean(INC.StressModel.yy));
INC.StressModel.sAvg=squeeze(mean(INC.StressModel.xy));

ISC.StressModel.xAvg=squeeze(mean(ISC.StressModel.xx));
ISC.StressModel.yAvg=squeeze(mean(ISC.StressModel.yy));
ISC.StressModel.sAvg=squeeze(mean(ISC.StressModel.xy));

%% Plot Stress Gauge (1-D) and constitutive model stresses Exact Inputs
figure(figure('units','normalized','outerposition',[0 0 1 1]))


%X-direction free edge
subplot(3,3,1)
plot(time,FE.stress.xAvg(FEPos.Free,:),'K')
hold on
plot(time,ENC.SG.x(GridPos.Free,:),'--b')
plot(time,ENC.StressModel.xAvg(GridPos.Free,:),':b')
plot(time,ESC.SG.x(GridPos.Free,:),'--r')
plot(time,ESC.StressModel.xAvg(GridPos.Free,:),':r')
hold off

xlabel('time (\mu{}s','FontSize',14)
ylabel('\sigma_{11} (Pa)','FontSize',14)
title('4P from Free Edge')
legend('FE','SG','Const','SG_{corr}','Const_{corr}','location','northwest')

%X-direction Middle
subplot(3,3,2)
plot(time,FE.stress.xAvg(FEPos.Mid,:),'K')
hold on
plot(time,ENC.SG.x(GridPos.Mid,:),'--b')
plot(time,ENC.StressModel.xAvg(GridPos.Mid,:),':b')
plot(time,ESC.SG.x(GridPos.Mid,:),'--r')
plot(time,ESC.StressModel.xAvg(GridPos.Mid,:),':r')
hold off

xlabel('time (\mu{}s','FontSize',14)
ylabel('\sigma_{11} (Pa)','FontSize',14)
title('Specimen Middle')


%X-direction free edge
subplot(3,3,3)
plot(time,FE.stress.xAvg(FEPos.Imp,:),'K')
hold on
plot(time,ENC.SG.x(GridPos.Imp,:),'--b')
plot(time,ENC.StressModel.xAvg(GridPos.Imp,:),':b')
plot(time,ESC.SG.x(GridPos.Imp,:),'--r')
plot(time,ESC.StressModel.xAvg(GridPos.Imp,:),':r')
hold off

xlabel('time (\mu{}s','FontSize',14)
ylabel('\sigma_{11} (Pa)','FontSize',14)
title('4P from Impact Edge')


%Y-direction free edge
subplot(3,3,4)
plot(time,FE.stress.yAvg(FEPos.Free,:),'K')
hold on
plot(time,ENC.StressModel.yAvg(GridPos.Free,:),':b')
plot(time,ESC.StressModel.yAvg(GridPos.Free,:),':r')
hold off

xlabel('time (\mu{}s','FontSize',14)
ylabel('\sigma_{22} (Pa)','FontSize',14)
title('4P from Free Edge')
legend('FE','Const','Const_{corr}','location','northwest')

%Y-direction Middle
subplot(3,3,5)
plot(time,FE.stress.yAvg(FEPos.Mid,:),'K')
hold on
plot(time,ENC.StressModel.yAvg(GridPos.Mid,:),':b')
plot(time,ESC.StressModel.yAvg(GridPos.Mid,:),':r')
hold off

xlabel('time (\mu{}s','FontSize',14)
ylabel('\sigma_{22} (Pa)','FontSize',14)
title('Specimen Middle')


%Y-direction free edge
subplot(3,3,6)
plot(time,FE.stress.yAvg(FEPos.Imp,:),'K')
hold on
plot(time,ENC.StressModel.yAvg(GridPos.Imp,:),':b')
plot(time,ESC.StressModel.yAvg(GridPos.Imp,:),':r')
hold off

xlabel('time (\mu{}s','FontSize',14)
ylabel('\sigma_{22} (Pa)','FontSize',14)
title('4P from Impact Edge')

%Shear free edge
subplot(3,3,7)
plot(time,FE.stress.sAvg(FEPos.Free,:),'K')
hold on
plot(time,ENC.SG.s(GridPos.Free,:),'--b')
plot(time,ENC.StressModel.sAvg(GridPos.Free,:),':b')
plot(time,ESC.SG.s(GridPos.Free,:),'--r')
plot(time,ESC.StressModel.sAvg(GridPos.Free,:),':r')
hold off

xlabel('time (\mu{}s','FontSize',14)
ylabel('\sigma_{12} (Pa)','FontSize',14)
title('4P from Free Edge')
legend('FE','SG','Const','SG_{corr}','Const_{corr}','location','northwest')

%Shear Middle
subplot(3,3,8)
plot(time,FE.stress.sAvg(FEPos.Mid,:),'K')
hold on
plot(time,ENC.SG.s(GridPos.Mid,:),'--b')
plot(time,ENC.StressModel.sAvg(GridPos.Mid,:),':b')
plot(time,ESC.SG.s(GridPos.Mid,:),'--r')
plot(time,ESC.StressModel.sAvg(GridPos.Mid,:),':r')
hold off

xlabel('time (\mu{}s','FontSize',14)
ylabel('\sigma_{12} (Pa)','FontSize',14)
title('Specimen Middle')


%Shear free edge
subplot(3,3,9)
plot(time,FE.stress.sAvg(FEPos.Imp,:),'K')
hold on
plot(time,ENC.SG.s(GridPos.Imp,:),'--b')
plot(time,ENC.StressModel.sAvg(GridPos.Imp,:),':b')
plot(time,ESC.SG.s(GridPos.Imp,:),'--r')
plot(time,ESC.StressModel.sAvg(GridPos.Imp,:),':r')
hold off

xlabel('time (\mu{}s','FontSize',14)
ylabel('\sigma_{12} (Pa)','FontSize',14)
title('4P from Impact Edge')

SaveName1Dstress=strcat(MainSaveDir,'/',ParentDesig, ...
    '_1DStressTime_ExactProps');
saveas(gcf,strcat(SaveName1Dstress,'.fig'))
saveas(gcf,strcat(SaveName1Dstress,'.svg'))
saveas(gcf,strcat(SaveName1Dstress,'.png'))
 

%% Plot Stress Gauge (1-D) and constitutive model stresses Identified Inputs
figure(figure('units','normalized','outerposition',[0 0 1 1]))


%X-direction free edge
subplot(3,3,1)
plot(time,FE.stress.xAvg(FEPos.Free,:),'K')
hold on
plot(time,INC.SG.x(GridPos.Free,:),'--b')
plot(time,INC.StressModel.xAvg(GridPos.Free,:),':b')
plot(time,ISC.SG.x(GridPos.Free,:),'--r')
plot(time,ISC.StressModel.xAvg(GridPos.Free,:),':r')
hold off

xlabel('time (\mu{}s','FontSize',14)
ylabel('\sigma_{11} (Pa)','FontSize',14)
title('4P from Free Edge')
legend('FE','SG','Const','SG_{corr}','Const_{corr}','location','northwest')

%X-direction Middle
subplot(3,3,2)
plot(time,FE.stress.xAvg(FEPos.Mid,:),'K')
hold on
plot(time,INC.SG.x(GridPos.Mid,:),'--b')
plot(time,INC.StressModel.xAvg(GridPos.Mid,:),':b')
plot(time,ISC.SG.x(GridPos.Mid,:),'--r')
plot(time,ISC.StressModel.xAvg(GridPos.Mid,:),':r')
hold off

xlabel('time (\mu{}s','FontSize',14)
ylabel('\sigma_{11} (Pa)','FontSize',14)
title('Specimen Middle')


%X-direction free edge
subplot(3,3,3)
plot(time,FE.stress.xAvg(FEPos.Imp,:),'K')
hold on
plot(time,INC.SG.x(GridPos.Imp,:),'--b')
plot(time,INC.StressModel.xAvg(GridPos.Imp,:),':b')
plot(time,ISC.SG.x(GridPos.Imp,:),'--r')
plot(time,ISC.StressModel.xAvg(GridPos.Imp,:),':r')
hold off

xlabel('time (\mu{}s','FontSize',14)
ylabel('\sigma_{11} (Pa)','FontSize',14)
title('4P from Impact Edge')


%Y-direction free edge
subplot(3,3,4)
plot(time,FE.stress.yAvg(FEPos.Free,:),'K')
hold on
plot(time,INC.StressModel.yAvg(GridPos.Free,:),':b')
plot(time,ISC.StressModel.yAvg(GridPos.Free,:),':r')
hold off

xlabel('time (\mu{}s','FontSize',14)
ylabel('\sigma_{22} (Pa)','FontSize',14)
title('4P from Free Edge')
legend('FE','Const','Const_{corr}','location','northwest')

%Y-direction Middle
subplot(3,3,5)
plot(time,FE.stress.yAvg(FEPos.Mid,:),'K')
hold on
plot(time,INC.StressModel.yAvg(GridPos.Mid,:),':b')
plot(time,ISC.StressModel.yAvg(GridPos.Mid,:),':r')
hold off

xlabel('time (\mu{}s','FontSize',14)
ylabel('\sigma_{22} (Pa)','FontSize',14)
title('Specimen Middle')


%Y-direction free edge
subplot(3,3,6)
plot(time,FE.stress.yAvg(FEPos.Imp,:),'K')
hold on
plot(time,INC.StressModel.yAvg(GridPos.Imp,:),':b')
plot(time,ISC.StressModel.yAvg(GridPos.Imp,:),':r')
hold off

xlabel('time (\mu{}s','FontSize',14)
ylabel('\sigma_{22} (Pa)','FontSize',14)
title('4P from Impact Edge')

%Shear free edge
subplot(3,3,7)
plot(time,FE.stress.sAvg(FEPos.Free,:),'K')
hold on
plot(time,INC.SG.s(GridPos.Free,:),'--b')
plot(time,INC.StressModel.sAvg(GridPos.Free,:),':b')
plot(time,ISC.SG.s(GridPos.Free,:),'--r')
plot(time,ISC.StressModel.sAvg(GridPos.Free,:),':r')
hold off

xlabel('time (\mu{}s','FontSize',14)
ylabel('\sigma_{12} (Pa)','FontSize',14)
title('4P from Free Edge')
legend('FE','SG','Const','SG_{corr}','Const_{corr}','location','northwest')

%Shear Middle
subplot(3,3,8)
plot(time,FE.stress.sAvg(FEPos.Mid,:),'K')
hold on
plot(time,INC.SG.s(GridPos.Mid,:),'--b')
plot(time,INC.StressModel.sAvg(GridPos.Mid,:),':b')
plot(time,ISC.SG.s(GridPos.Mid,:),'--r')
plot(time,ISC.StressModel.sAvg(GridPos.Mid,:),':r')
hold off

xlabel('time (\mu{}s','FontSize',14)
ylabel('\sigma_{12} (Pa)','FontSize',14)
title('Specimen Middle')


%Shear free edge
subplot(3,3,9)
plot(time,FE.stress.sAvg(FEPos.Imp,:),'K')
hold on
plot(time,INC.SG.s(GridPos.Imp,:),'--b')
plot(time,INC.StressModel.sAvg(GridPos.Imp,:),':b')
plot(time,ISC.SG.s(GridPos.Imp,:),'--r')
plot(time,ISC.StressModel.sAvg(GridPos.Imp,:),':r')
hold off

xlabel('time (\mu{}s','FontSize',14)
ylabel('\sigma_{12} (Pa)','FontSize',14)
title('4P from Impact Edge')

SaveName1Dstress=strcat(MainSaveDir,'/',ParentDesig, ...
    '_1DStressTime_IdentProps');
saveas(gcf,strcat(SaveName1Dstress,'.fig'))
saveas(gcf,strcat(SaveName1Dstress,'.svg'))
saveas(gcf,strcat(SaveName1Dstress,'.png'))

%% Calculate average FE strains
FE.strain.xAvg=squeeze(mean(FE.strain.x));
FE.strain.yAvg=squeeze(mean(FE.strain.y));
FE.strain.sAvg=squeeze(mean(FE.strain.s));
%% Plot 1-D stress-strain identified parameters
figure(figure('units','normalized','outerposition',[0 0 1 1]))


%X-direction free edge
subplot(3,3,1)
plot(FE.strain.xAvg(FEPos.Free,:),FE.stress.xAvg(FEPos.Free,:),'K')
hold on
plot(INC.strain.xAvg(GridPos.Free,:),INC.SG.x(GridPos.Free,:),'--b')
plot(INC.strain.xAvg(GridPos.Free,:),INC.StressModel.xAvg(GridPos.Free,:),':b')
plot(ISC.strain.xAvg(GridPos.Free,:),ISC.SG.x(GridPos.Free,:),'--r')
plot(ISC.strain.xAvg(GridPos.Free,:),ISC.StressModel.xAvg(GridPos.Free,:),':r')
hold off

xlabel('$\varepsilon_{11}$','interpreter','latex','FontSize',14)
ylabel('\sigma_{11} (Pa)','FontSize',14)
title('4P from Free Edge')
legend('FE','SG','Const','SG_{corr}','Const_{corr}','location','northwest')

%X-direction Middle
subplot(3,3,2)
plot(FE.strain.xAvg(FEPos.Mid,:),FE.stress.xAvg(FEPos.Mid,:),'K')
hold on
plot(INC.strain.xAvg(GridPos.Mid,:),INC.SG.x(GridPos.Mid,:),'--b')
plot(INC.strain.xAvg(GridPos.Mid,:),INC.StressModel.xAvg(GridPos.Mid,:),':b')
plot(ISC.strain.xAvg(GridPos.Mid,:),ISC.SG.x(GridPos.Mid,:),'--r')
plot(ISC.strain.xAvg(GridPos.Mid,:),ISC.StressModel.xAvg(GridPos.Mid,:),':r')
hold off

xlabel('$\varepsilon_{11}$','interpreter','latex','FontSize',14)
ylabel('\sigma_{11} (Pa)','FontSize',14)
title('Specimen Middle')


%X-direction free edge
subplot(3,3,3)
plot(FE.strain.xAvg(FEPos.Imp,:),FE.stress.xAvg(FEPos.Imp,:),'K')
hold on
plot(INC.strain.xAvg(GridPos.Imp,:),INC.SG.x(GridPos.Imp,:),'--b')
plot(INC.strain.xAvg(GridPos.Imp,:),INC.StressModel.xAvg(GridPos.Imp,:),':b')
plot(ISC.strain.xAvg(GridPos.Imp,:),ISC.SG.x(GridPos.Imp,:),'--r')
plot(ISC.strain.xAvg(GridPos.Imp,:),ISC.StressModel.xAvg(GridPos.Imp,:),':r')
hold off

xlabel('$\varepsilon_{11}$','interpreter','latex','FontSize',14)
ylabel('\sigma_{11} (Pa)','FontSize',14)
title('4P from Impact Edge')


%Y-direction free edge
subplot(3,3,4)
plot(FE.strain.yAvg(FEPos.Free,:),FE.stress.yAvg(FEPos.Free,:),'K')
hold on
plot(INC.strain.yAvg(GridPos.Free,:),INC.StressModel.yAvg(GridPos.Free,:),':b')
plot(ISC.strain.yAvg(GridPos.Free,:),ISC.StressModel.yAvg(GridPos.Free,:),':r')
hold off

xlabel('$\varepsilon_{22}$','interpreter','latex','FontSize',14)
ylabel('\sigma_{22} (Pa)','FontSize',14)
title('4P from Free Edge')
legend('FE','Const','Const_{corr}','location','northwest')

%Y-direction Middle
subplot(3,3,5)
plot(FE.strain.yAvg(FEPos.Mid,:),FE.stress.yAvg(FEPos.Mid,:),'K')
hold on
plot(INC.strain.yAvg(GridPos.Mid,:),INC.StressModel.yAvg(GridPos.Mid,:),':b')
plot(ISC.strain.yAvg(GridPos.Mid,:),ISC.StressModel.yAvg(GridPos.Mid,:),':r')
hold off

xlabel('$\varepsilon_{22}$','interpreter','latex','FontSize',14)
ylabel('\sigma_{22} (Pa)','FontSize',14)
title('Specimen Middle')


%Y-direction free edge
subplot(3,3,6)
plot(FE.strain.yAvg(FEPos.Imp,:),FE.stress.yAvg(FEPos.Imp,:),'K')
hold on
plot(INC.strain.yAvg(GridPos.Imp,:),INC.StressModel.yAvg(GridPos.Imp,:),':b')
plot(ISC.strain.yAvg(GridPos.Imp,:),ISC.StressModel.yAvg(GridPos.Imp,:),':r')
hold off

xlabel('$\varepsilon_{22}$','interpreter','latex','FontSize',14)
ylabel('\sigma_{22} (Pa)','FontSize',14)
title('4P from Impact Edge')

%Shear free edge
subplot(3,3,7)
plot(FE.strain.sAvg(FEPos.Free,:),FE.stress.sAvg(FEPos.Free,:),'K')
hold on
plot(INC.strain.sAvg(GridPos.Free,:),INC.SG.s(GridPos.Free,:),'--b')
plot(INC.strain.sAvg(GridPos.Free,:),INC.StressModel.sAvg(GridPos.Free,:),':b')
plot(ISC.strain.sAvg(GridPos.Free,:),ISC.SG.s(GridPos.Free,:),'--r')
plot(ISC.strain.sAvg(GridPos.Free,:),ISC.StressModel.sAvg(GridPos.Free,:),':r')
hold off

xlabel('\gamma','FontSize',14)
ylabel('\sigma_{12} (Pa)','FontSize',14)
title('4P from Free Edge')
legend('FE','SG','Const','SG_{corr}','Const_{corr}','location','northwest')

%Shear Middle
subplot(3,3,8)
plot(FE.strain.sAvg(FEPos.Mid,:),FE.stress.sAvg(FEPos.Mid,:),'K')
hold on
plot(INC.strain.sAvg(GridPos.Mid,:),INC.SG.s(GridPos.Mid,:),'--b')
plot(INC.strain.sAvg(GridPos.Mid,:),INC.StressModel.sAvg(GridPos.Mid,:),':b')
plot(ISC.strain.sAvg(GridPos.Mid,:),ISC.SG.s(GridPos.Mid,:),'--r')
plot(ISC.strain.sAvg(GridPos.Mid,:),ISC.StressModel.sAvg(GridPos.Mid,:),':r')
hold off

xlabel('\gamma','FontSize',14)
ylabel('\sigma_{12} (Pa)','FontSize',14)
title('Specimen Middle')


%Shear free edge
subplot(3,3,9)
plot(FE.strain.sAvg(FEPos.Imp,:),FE.stress.sAvg(FEPos.Imp,:),'K')
hold on
plot(INC.strain.sAvg(GridPos.Imp,:),INC.SG.s(GridPos.Imp,:),'--b')
plot(INC.strain.sAvg(GridPos.Imp,:),INC.StressModel.sAvg(GridPos.Imp,:),':b')
plot(ISC.strain.sAvg(GridPos.Imp,:),ISC.SG.s(GridPos.Imp,:),'--r')
plot(ISC.strain.sAvg(GridPos.Imp,:),ISC.StressModel.sAvg(GridPos.Imp,:),':r')
hold off

xlabel('\gamma','FontSize',14)
ylabel('\sigma_{12} (Pa)','FontSize',14)
title('4P from Impact Edge')

SaveName1Dstress=strcat(MainSaveDir,'/',ParentDesig, ...
    '_1DStressStrain_IdentProps');
saveas(gcf,strcat(SaveName1Dstress,'.fig'))
saveas(gcf,strcat(SaveName1Dstress,'.svg'))
saveas(gcf,strcat(SaveName1Dstress,'.png'))
%% Plot 1-D stress-strain exact parameters
figure(figure('units','normalized','outerposition',[0 0 1 1]))


%X-direction free edge
subplot(3,3,1)
plot(FE.strain.xAvg(FEPos.Free,:),FE.stress.xAvg(FEPos.Free,:),'K')
hold on
plot(ENC.strain.xAvg(GridPos.Free,:),ENC.SG.x(GridPos.Free,:),'--b')
plot(ENC.strain.xAvg(GridPos.Free,:),ENC.StressModel.xAvg(GridPos.Free,:),':b')
plot(ESC.strain.xAvg(GridPos.Free,:),ESC.SG.x(GridPos.Free,:),'--r')
plot(ESC.strain.xAvg(GridPos.Free,:),ESC.StressModel.xAvg(GridPos.Free,:),':r')
hold off

xlabel('$\varepsilon_{11}$','interpreter','latex','FontSize',14)
ylabel('\sigma_{11} (Pa)','FontSize',14)
title('4P from Free Edge')
legend('FE','SG','Const','SG_{corr}','Const_{corr}','location','northwest')

%X-direction Middle
subplot(3,3,2)
plot(FE.strain.xAvg(FEPos.Mid,:),FE.stress.xAvg(FEPos.Mid,:),'K')
hold on
plot(ENC.strain.xAvg(GridPos.Mid,:),ENC.SG.x(GridPos.Mid,:),'--b')
plot(ENC.strain.xAvg(GridPos.Mid,:),ENC.StressModel.xAvg(GridPos.Mid,:),':b')
plot(ESC.strain.xAvg(GridPos.Mid,:),ESC.SG.x(GridPos.Mid,:),'--r')
plot(ESC.strain.xAvg(GridPos.Mid,:),ESC.StressModel.xAvg(GridPos.Mid,:),':r')
hold off

xlabel('$\varepsilon_{11}$','interpreter','latex','FontSize',14)
ylabel('\sigma_{11} (Pa)','FontSize',14)
title('Specimen Middle')


%X-direction free edge
subplot(3,3,3)
plot(FE.strain.xAvg(FEPos.Imp,:),FE.stress.xAvg(FEPos.Imp,:),'K')
hold on
plot(ENC.strain.xAvg(GridPos.Imp,:),ENC.SG.x(GridPos.Imp,:),'--b')
plot(ENC.strain.xAvg(GridPos.Imp,:),ENC.StressModel.xAvg(GridPos.Imp,:),':b')
plot(ESC.strain.xAvg(GridPos.Imp,:),ESC.SG.x(GridPos.Imp,:),'--r')
plot(ESC.strain.xAvg(GridPos.Imp,:),ESC.StressModel.xAvg(GridPos.Imp,:),':r')
hold off

xlabel('$\varepsilon_{11}$','interpreter','latex','FontSize',14)
ylabel('\sigma_{11} (Pa)','FontSize',14)
title('4P from Impact Edge')


%Y-direction free edge
subplot(3,3,4)
plot(FE.strain.yAvg(FEPos.Free,:),FE.stress.yAvg(FEPos.Free,:),'K')
hold on
plot(ENC.strain.yAvg(GridPos.Free,:),ENC.StressModel.yAvg(GridPos.Free,:),':b')
plot(ESC.strain.yAvg(GridPos.Free,:),ESC.StressModel.yAvg(GridPos.Free,:),':r')
hold off

xlabel('$\varepsilon_{22}$','interpreter','latex','FontSize',14)
ylabel('\sigma_{22} (Pa)','FontSize',14)
title('4P from Free Edge')
legend('FE','Const','Const_{corr}','location','northwest')

%Y-direction Middle
subplot(3,3,5)
plot(FE.strain.yAvg(FEPos.Mid,:),FE.stress.yAvg(FEPos.Mid,:),'K')
hold on
plot(ENC.strain.yAvg(GridPos.Mid,:),ENC.StressModel.yAvg(GridPos.Mid,:),':b')
plot(ESC.strain.yAvg(GridPos.Mid,:),ESC.StressModel.yAvg(GridPos.Mid,:),':r')
hold off

xlabel('$\varepsilon_{22}$','interpreter','latex','FontSize',14)
ylabel('\sigma_{22} (Pa)','FontSize',14)
title('Specimen Middle')


%Y-direction free edge
subplot(3,3,6)
plot(FE.strain.yAvg(FEPos.Imp,:),FE.stress.yAvg(FEPos.Imp,:),'K')
hold on
plot(ENC.strain.yAvg(GridPos.Imp,:),ENC.StressModel.yAvg(GridPos.Imp,:),':b')
plot(ESC.strain.yAvg(GridPos.Imp,:),ESC.StressModel.yAvg(GridPos.Imp,:),':r')
hold off

xlabel('$\varepsilon_{22}$','interpreter','latex','FontSize',14)
ylabel('\sigma_{22} (Pa)','FontSize',14)
title('4P from Impact Edge')

%Shear free edge
subplot(3,3,7)
plot(FE.strain.sAvg(FEPos.Free,:),FE.stress.sAvg(FEPos.Free,:),'K')
hold on
plot(ENC.strain.sAvg(GridPos.Free,:),ENC.SG.s(GridPos.Free,:),'--b')
plot(ENC.strain.sAvg(GridPos.Free,:),ENC.StressModel.sAvg(GridPos.Free,:),':b')
plot(ESC.strain.sAvg(GridPos.Free,:),ESC.SG.s(GridPos.Free,:),'--r')
plot(ESC.strain.sAvg(GridPos.Free,:),ESC.StressModel.sAvg(GridPos.Free,:),':r')
hold off

xlabel('\gamma','FontSize',14)
ylabel('\sigma_{12} (Pa)','FontSize',14)
title('4P from Free Edge')
legend('FE','SG','Const','SG_{corr}','Const_{corr}','location','northwest')

%Shear Middle
subplot(3,3,8)
plot(FE.strain.sAvg(FEPos.Mid,:),FE.stress.sAvg(FEPos.Mid,:),'K')
hold on
plot(ENC.strain.sAvg(GridPos.Mid,:),ENC.SG.s(GridPos.Mid,:),'--b')
plot(ENC.strain.sAvg(GridPos.Mid,:),ENC.StressModel.sAvg(GridPos.Mid,:),':b')
plot(ESC.strain.sAvg(GridPos.Mid,:),ESC.SG.s(GridPos.Mid,:),'--r')
plot(ESC.strain.sAvg(GridPos.Mid,:),ESC.StressModel.sAvg(GridPos.Mid,:),':r')
hold off

xlabel('\gamma','FontSize',14)
ylabel('\sigma_{12} (Pa)','FontSize',14)
title('Specimen Middle')


%Shear free edge
subplot(3,3,9)
plot(FE.strain.sAvg(FEPos.Imp,:),FE.stress.sAvg(FEPos.Imp,:),'K')
hold on
plot(ENC.strain.sAvg(GridPos.Imp,:),ENC.SG.s(GridPos.Imp,:),'--b')
plot(ENC.strain.sAvg(GridPos.Imp,:),ENC.StressModel.sAvg(GridPos.Imp,:),':b')
plot(ESC.strain.sAvg(GridPos.Imp,:),ESC.SG.s(GridPos.Imp,:),'--r')
plot(ESC.strain.sAvg(GridPos.Imp,:),ESC.StressModel.sAvg(GridPos.Imp,:),':r')
hold off

xlabel('\gamma','FontSize',14)
ylabel('\sigma_{12} (Pa)','FontSize',14)
title('4P from Impact Edge')

SaveName1Dstress=strcat(MainSaveDir,'/',ParentDesig, ...
    '_1DStressStrain_ExactProps');
saveas(gcf,strcat(SaveName1Dstress,'.fig'))
saveas(gcf,strcat(SaveName1Dstress,'.svg'))
saveas(gcf,strcat(SaveName1Dstress,'.png'))

%% Calculate Errors
fprintf('Calculating Stress Errors \n')
ENC.SGxErr.Free=(ENC.SG.x(GridPos.Free,:)-FE.stress.xAvg(FEPos.Free,:))/...
    FE.stress.xAvg(FEPos.Free,:)*100;
ENC.SGsErr.Free=(ENC.SG.s(GridPos.Free,:)-FE.stress.sAvg(FEPos.Free,:))/...
    FE.stress.sAvg(FEPos.Free,:)*100;
ENC.ConxErr.Free=(ENC.StressModel.xAvg(GridPos.Free,:)-...
    FE.stress.xAvg(FEPos.Free,:))/FE.stress.xAvg(FEPos.Free,:)*100;
ENC.ConsErr.Free=(ENC.StressModel.sAvg(GridPos.Free,:)-...
    FE.stress.sAvg(FEPos.Free,:))/FE.stress.sAvg(FEPos.Free,:)*100;
ENC.ConyErr.Free=(ENC.StressModel.yAvg(GridPos.Free,:)-...
    FE.stress.yAvg(FEPos.Free,:))/FE.stress.yAvg(FEPos.Free,:)*100;

INC.SGxErr.Free=(INC.SG.x(GridPos.Free,:)-FE.stress.xAvg(FEPos.Free,:))/...
    INC.StressModel.xAvg(FEPos.Free,:)*100;
INC.SGsErr.Free=(INC.SG.s(GridPos.Free,:)-FE.stress.sAvg(FEPos.Free,:))/...
    FE.stress.sAvg(FEPos.Free,:)*100;
INC.ConxErr.Free=(INC.StressModel.xAvg(GridPos.Free,:)-...
    FE.stress.xAvg(FEPos.Free,:))/FE.stress.xAvg(FEPos.Free,:)*100;
INC.ConsErr.Free=(INC.StressModel.sAvg(GridPos.Free,:)-...
    FE.stress.sAvg(FEPos.Free,:))/FE.stress.sAvg(FEPos.Free,:)*100;
INC.ConyErr.Free=(INC.StressModel.yAvg(GridPos.Free,:)-...
    FE.stress.yAvg(FEPos.Free,:))/FE.stress.yAvg(FEPos.Free,:)*100;

%Specimen Middle
ENC.SGxErr.Free=(ENC.SG.x(GridPos.Free,:)-FE.stress.xAvg(FEPos.Free,:))/...
    FE.stress.xAvg(FEPos.Free,:)*100;
ENC.SGsErr.Free=(ENC.SG.s(GridPos.Free,:)-FE.stress.sAvg(FEPos.Free,:))/...
    FE.stress.sAvg(FEPos.Free,:)*100;
ENC.ConxErr.Free=(ENC.StressModel.xAvg(GridPos.Free,:)-...
    FE.stress.xAvg(FEPos.Free,:))/FE.stress.xAvg(FEPos.Free,:)*100;
ENC.ConsErr.Free=(ENC.StressModel.sAvg(GridPos.Free,:)-...
    FE.stress.sAvg(FEPos.Free,:))/FE.stress.sAvg(FEPos.Free,:)*100;
ENC.ConyErr.Free=(ENC.StressModel.yAvg(GridPos.Free,:)-...
    FE.stress.yAvg(FEPos.Free,:))/FE.stress.yAvg(FEPos.Free,:)*100;

INC.SGxErr.Free=(INC.SG.x(GridPos.Free,:)-FE.stress.xAvg(FEPos.Free,:))/...
    INC.StressModel.xAvg(FEPos.Free,:)*100;
INC.SGsErr.Free=(INC.SG.s(GridPos.Free,:)-FE.stress.sAvg(FEPos.Free,:))/...
    FE.stress.sAvg(FEPos.Free,:)*100;
INC.ConxErr.Free=(INC.StressModel.xAvg(GridPos.Free,:)-...
    FE.stress.xAvg(FEPos.Free,:))/FE.stress.xAvg(FEPos.Free,:)*100;
INC.ConsErr.Free=(INC.StressModel.sAvg(GridPos.Free,:)-...
    FE.stress.sAvg(FEPos.Free,:))/FE.stress.sAvg(FEPos.Free,:)*100;
INC.ConyErr.Free=(INC.StressModel.yAvg(GridPos.Free,:)-...
    FE.stress.yAvg(FEPos.Free,:))/FE.stress.yAvg(FEPos.Free,:)*100;

%Impact edge
ENC.SGxErr.Free=(ENC.SG.x(GridPos.Free,:)-FE.stress.xAvg(FEPos.Free,:))/...
    FE.stress.xAvg(FEPos.Free,:)*100;
ENC.SGsErr.Free=(ENC.SG.s(GridPos.Free,:)-FE.stress.sAvg(FEPos.Free,:))/...
    FE.stress.sAvg(FEPos.Free,:)*100;
ENC.ConxErr.Free=(ENC.StressModel.xAvg(GridPos.Free,:)-...
    FE.stress.xAvg(FEPos.Free,:))/FE.stress.xAvg(FEPos.Free,:)*100;
ENC.ConsErr.Free=(ENC.StressModel.sAvg(GridPos.Free,:)-...
    FE.stress.sAvg(FEPos.Free,:))/FE.stress.sAvg(FEPos.Free,:)*100;
ENC.ConyErr.Free=(ENC.StressModel.yAvg(GridPos.Free,:)-...
    FE.stress.yAvg(FEPos.Free,:))/FE.stress.yAvg(FEPos.Free,:)*100;

INC.SGxErr.Free=(INC.SG.x(GridPos.Free,:)-FE.stress.xAvg(FEPos.Free,:))/...
    INC.StressModel.xAvg(FEPos.Free,:)*100;
INC.SGsErr.Free=(INC.SG.s(GridPos.Free,:)-FE.stress.sAvg(FEPos.Free,:))/...
    FE.stress.sAvg(FEPos.Free,:)*100;
INC.ConxErr.Free=(INC.StressModel.xAvg(GridPos.Free,:)-...
    FE.stress.xAvg(FEPos.Free,:))/FE.stress.xAvg(FEPos.Free,:)*100;
INC.ConsErr.Free=(INC.StressModel.sAvg(GridPos.Free,:)-...
    FE.stress.sAvg(FEPos.Free,:))/FE.stress.sAvg(FEPos.Free,:)*100;
INC.ConyErr.Free=(INC.StressModel.yAvg(GridPos.Free,:)-...
    FE.stress.yAvg(FEPos.Free,:))/FE.stress.yAvg(FEPos.Free,:)*100;

%Shear corrected

%% Plot shear cost function only

%% Plot Bulk Cost Function only

%% 