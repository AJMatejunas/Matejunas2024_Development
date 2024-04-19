% This code is written to compare the evolution of instantaneous shear
    % modulus and Poisson's ratio
    
    clear all; close all; clc
 %% Choose data and test designation
 desig=char(inputdlg('Input Test Designation'));
 
 [FEfile,FEpath]=uigetfile('*.mat','Load file containing FE stresses');
 FEname=strcat(FEpath,'/',FEfile);
 
 [Modelfile,Modelpath]=uigetfile('*.mat',...
     'Load file containing constitutive model data');
 Modelname=strcat(Modelpath,'/',Modelfile);
 
 
  
 %% Load data file and define test deg
fprintf('Loading FE data \n')
FE=load(FEname);

fprintf('Loading Constitutive Model Data \n')
Model=load(Modelname);
testdeg=desig;
clear desig

%% Create position an dtime vectors
X_vec=FE.X_vec;
Y_vec=FE.Y_vec;
timemic=FE.time.vec*10^6; %us

Xcoord=X_vec*10^3;
Ycoord=Y_vec*10^3;

%% Calculate shear moduli
fprintf('Calculating Instantneous Shear Modulus for FE data \n')
FE.Ginst=FE.stress.s./FE.strain.s;

fprintf('Calculating instantaneous shear modulus for Model data \n')
Model.Ginst=Model.StressModel.xy./FE.strain.s;

%ratio of instantaneous shear moduli
GinstRat=Model.Ginst./FE.Ginst;

%% Generate limits 
FE.Gmax=max(FE.Ginst,[],'all');
FE.Gmin=min(FE.Ginst,[],'all');

Model.Gmax=max(Model.Ginst,[],'all');
Model.Gmin=min(Model.Ginst,[],'all');

GLim=[1E9,2E9];

RatLim=[0,1];

%% generate plots
figure('units','normalized','outerposition',[0 0 1 1])
for n=1:length(timemic)
   timecount=num2str(timemic(n)); 
   framecount=num2str(n); 
   
  FEinst=squeeze(FE.Ginst(:,:,n));
  Modelinst=squeeze(Model.Ginst(:,:,n));
  Rat=squeeze(GinstRat(:,:,n));   
   
   %% Plot Instanteous Shear Modulus FE
    subplot(3,1,1)
    imagesc(Xcoord,Ycoord,FEinst)
    title(strcat(' FE G_{inst}',timecount,'\mus'))
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='G_{inst} (Pa)';
    caxis(GLim)
    
       
  %% plot Instantaneous shear modul from constitutive model
    subplot(3,1,2)
    imagesc(Xcoord,Ycoord,Modelinst)
    title(strcat('Model G_{inst}',timecount,'\mus'))
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='G_{inst} (Pa)';
    caxis(GLim)

 %% Plot instantaneous shear modulus ratio
    subplot(3,1,3)
    imagesc(Xcoord,Ycoord,Rat)
    title(strcat('G_{Model}/G_{FE} ',timecount,'\mus'))
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='G_{Model}/G_{FE}';
    caxis(RatLim)
    

a=findall(gcf,'Type','axes');
h=findall(gcf,'Type','line');
t=findall(gcf,'Type','text');

set(t, ...
    'FontName','Helvetica', ...
    'FontSize',18, ...
    'FontWeight','bold' );

set(a                ,               ...
    'FontName'       , 'Helvetica' , ...
    'FontSize'       , 18     , ...
    'FontWeight'     , 'normal');

    
 saveas(gcf,strcat(testdeg,'_GinstComp_',framecount,'.fig'))
 saveas(gcf,strcat(testdeg,'_GinstComp_',framecount,'.png'))
end
