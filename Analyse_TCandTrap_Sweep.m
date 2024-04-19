%This code is written to analyze the time sweep of viscoelatic constitutive
    %parameter identifications performed on a parametric sweep of time
    %constants and trapezoidal pulse durations

%% Initialize
clearvars; close all; clc;

%% Create Parameter Vectors
TCvec=[0.1e-6,1e-6,10e-6,100e-6,1000e-6];
RTvec=[1e-6,5e-6,10e-6,15e-6,20e-6,25e-6];
PWvec=RTvec*3;

%% Chose Save Directory
FigSaveDir=uigetdir('','Choose Directory to Save Evaluation Figure');

%% Add reference parameters file 
[RefPar.file,RefPar.path]=uigetfile('*.mat',...
    'Choose Mat file containing reference constitutive pamrameters');
RefPar=load(strcat(RefPar.path,'/',RefPar.file),'MatProps');

%% Create sweep variables
SweepDesig=char(inputdlg('Sweep Designation'));
MainDir=uigetdir('Choose Home Directory');
[DesigVec,SaveDirVec,FileNameVec]=...
func_BatchProcessVars(SweepDesig,MainDir);

%% Reshape Sweep Into A proper Heat Map
DesigArray=reshape(DesigVec,[length(TCvec),length(RTvec)]);
LocDirArray=reshape(SaveDirVec,[length(TCvec),length(RTvec)]);
FileNameArray=reshape(FileNameVec,[length(TCvec),length(RTvec)]);

%% Initiate the array of file names that are loaded
LoadFileArray=cell(length(TCvec),length(PWvec));
FileEndChar='_KGident.mat';

%% Initiate Constitutive Parameter arrays
Ident.K=zeros(length(TCvec),length(PWvec));
Ident.G=zeros(length(TCvec),length(PWvec));
Ident.tau=zeros(length(TCvec),length(PWvec));

Errors.tau=zeros(length(TCvec),length(PWvec));

%% Read the constitutive Parameter Files
LoadTotIter=length(TCvec)*length(PWvec);
LoadIt=0;
LoadProg=LoadIt/LoadTotIter;
LoadMsg=strcat(num2str(LoadIt),'/',num2str(LoadTotIter),...
'Identification Files Loaded');
LoadProgress=waitbar(LoadProg,LoadMsg);

for jj=1:length(TCvec)
    for kk=1:length(PWvec)
        %% Set Up Progress Bar
        LoadIt=LoadIt+1;
        LoadProg=LoadIt/LoadTotIter;

        %% Find File
        LoadDir=LocDirArray{jj,kk};
        LoadDesig=DesigArray{jj,kk};
        fprintf(strcat('Loading~',LoadDesig,'~Identification file \n'))
        LoadFileName=strcat(LoadDir,LoadDesig,FileEndChar);
        LoadFileArray{jj,kk}=LoadFileName;
        load(LoadFileName,'-mat','ConstitutiveParam')
        
        %% Save Constitutive Parameters
        Ident.K(jj,kk)=ConstitutiveParam(1);
        Ident.G(jj,kk)=ConstitutiveParam(2);
        Ident.tau(jj,kk)=ConstitutiveParam(3);
        
        %% Calculate time constant identification errors
        Errors.tau(jj,kk)=(Ident.tau(jj,kk)-TCvec(jj))/TCvec(jj)*100;

        %% Create Progress Bar
        LoadMsg=strcat(num2str(LoadIt),'/',num2str(LoadTotIter),...
            'Identification Files Loaded');
        LoadProgress=waitbar(LoadProg,LoadProgress,LoadMsg);
        %% Clear Constitutive Parameters file
        clearvars ConstitutiveParam
    end
end
fprintf('Data Loading Complete \n')
fprintf('Calculating and Plotting Errors \n')
clear jj kk LoadProgress LoadMsg LoadProg

%% Calculate Identification Errors
Errors.K=(Ident.K-RefPar.MatProps.Ki)/RefPar.MatProps.Ki*100;
Errors.G=(Ident.G-RefPar.MatProps.Gi)/RefPar.MatProps.Gi*100;

%% Save data
fprintf('Saving Identification Results \n')
SaveName=strcat(FigSaveDir,'/',SweepDesig,'_IdentErrors.mat');
save(SaveName)

%% Plot Identfication Errors
figure('units','Centimeters','InnerPosition',[10 10 18 12])
ColorScheme='jet';

% K errors
subplot(2,2,1)
imagesc(TCvec,PWvec*10^6,Errors.K')
title('K_1 identification Errors')
xlabel('Time Constant (s)')
ylabel('Pulse Duration (\mu{}s)')
colormap(ColorScheme)
cx=colorbar;
cx.Label.String='Error (%)';
set(gca,'YDir','normal')

% G errors
subplot(2,2,2)
imagesc(TCvec,PWvec*10^6,Errors.G')
title('G_1 identification Errors')
xlabel('Time Constant (s)')
ylabel('Pulse Duration (\mu{}s)')
colormap(ColorScheme)
cx=colorbar;
cx.Label.String='Error (%)';
set(gca,'YDir','normal')

% tau errors
subplot(2,2,3.5)
imagesc(TCvec,PWvec*10^6,Errors.tau')
title('\tau_1 identification Errors')
xlabel('Time Constant (s)')
ylabel('Pulse Duration (\mu{}s)')
colormap(ColorScheme)
cx=colorbar;
cx.Label.String='Error (%)';
set(gca,'YDir','normal')

saveas(gcf,strcat(FigSaveDir,'\TC_TrapSweepErr.fig'))
saveas(gcf,strcat(FigSaveDir,'\TC_TrapSweepErr.svg'))
saveas(gcf,strcat(FigSaveDir,'\TC_TrapSweepErr.eps'))

%% Plot as Contour Plots
figure('units','Centimeters','InnerPosition',[10 10 18 12])
ColorScheme='jet';

% K errors
subplot(2,2,1)
kc=contourf(TCvec,PWvec*10^6,Errors.K');
title('K_1 identification Errors')
xlabel('Time Constant (s)')
ylabel('Pulse Duration (\mu{}s)')
colormap(ColorScheme)
cx=colorbar;
cx.Label.String='Error (%)';
caxis([-5,5])
% set(gca,'YDir','normal')
% set(kc,'LineColor','none')

% G errors
subplot(2,2,2)
gc=contourf(TCvec,PWvec*10^6,Errors.G');
title('G_1 identification Errors')
xlabel('Time Constant (s)')
ylabel('Pulse Duration (\mu{}s)')
colormap(ColorScheme)
cx=colorbar;
cx.Label.String='Error (%)';
% set(gca,'YDir','normal')
% set(gc,'LineColor','none')

% tau errors
subplot(2,2,3.5)
tc=contourf(TCvec,PWvec*10^6,Errors.tau');
title('\tau_1 identification Errors')
xlabel('Time Constant (s)')
ylabel('Pulse Duration (\mu{}s)')
colormap(ColorScheme)
cx=colorbar;
cx.Label.String='Error (%)';
% set(gca,'YDir','normal')
% set(tc,'LineColor','none')

saveas(gcf,strcat(FigSaveDir,'\TC_TrapSweepErrContour.fig'))
saveas(gcf,strcat(FigSaveDir,'\TC_TrapSweepErrContour.svg'))
saveas(gcf,strcat(FigSaveDir,'\TC_TrapSweepErrContour.eps'))


%% Plot the identified parameters
figure('units','Centimeters','InnerPosition',[10 10 18 12])
ColorScheme='jet';

% K errors
subplot(2,2,1)
kc=contourf(TCvec,PWvec*10^6,Ident.K');
title('K_1 identification Errors')
xlabel('Time Constant (s)')
ylabel('Pulse Duration (\mu{}s)')
colormap(ColorScheme)
cx=colorbar;
cx.Label.String='K_1';
caxis([-5,5])
% set(gca,'YDir','normal')
% set(kc,'LineColor','none')

% G errors
subplot(2,2,2)
gc=contourf(TCvec,PWvec*10^6,Ident.G');
title('G_1 identification Errors')
xlabel('Time Constant (s)')
ylabel('Pulse Duration (\mu{}s)')
colormap(ColorScheme)
cx=colorbar;
cx.Label.String='G_1';
% set(gca,'YDir','normal')
% set(gc,'LineColor','none')

% tau errors
subplot(2,2,3.5)
tc=contourf(TCvec,PWvec*10^6,Ident.tau');
title('\tau_1 identification Errors')
xlabel('Time Constant (s)')
ylabel('Pulse Duration (\mu{}s)')
colormap(ColorScheme)
cx=colorbar;
cx.Label.String='\tau';
% set(gca,'YDir','normal')
% set(tc,'LineColor','none')

saveas(gcf,strcat(FigSaveDir,'\TC_TrapSweepErrContour.fig'))
saveas(gcf,strcat(FigSaveDir,'\TC_TrapSweepErrContour.svg'))
saveas(gcf,strcat(FigSaveDir,'\TC_TrapSweepErrContour.eps'))

