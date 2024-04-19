%% Initialize
clear all, close all, clc

%% Choose File Containing FE data
[FEfile,FEpath]=uigetfile('*.mat','Choose File Contating FE Data');
FEname=strcat(FEpath,'/',FEfile);

%% Choose File raw Containing Displacement Data
[GMfile,GMpath]=uigetfile('*.mat','Choose File Containing GM data');
GMname=strcat(GMpath,'/',GMfile);

%% Choose Folder to Save Results
SavePath=uigetdir({},'Choose Parent Folder to Save Results');

%% Define designation
Desig=char(inputdlg('Test Designation'));


%% Load
fprintf('Loading FE displacement data \n')
FE=load(FEname,'pos','disp');

fprintf('Loading Grid Method Displacement Data \n')
load(GMname,'pos','grid','disp','time');

fprintf('Data Loading Complete \n')
%% Extract X, Y, and time vectors
FE.X_vec=pos.x;
FE.Y_vec=pos.y;
FE.Xmm=FE.X_vec*10^3;
FE.Ymm=FE.Y_vec*10^3;

X_vec=pos.x;
Y_vec=pos.y;

Xmm=FE.X_vec*10^3;
Ymm=FE.Y_vec*10^3;
timemic=time.vec*10^6; %Time vector in microseconds

%% Calculate Limits for plot
FE.plotlim.x=[min(FE.disp.x,[],'all'),max(FE.disp.x,[],'all')];
FE.plotlim.y=[min(FE.disp.y,[],'all'),max(FE.disp.y,[],'all')];

plotlim.x=[min(disp.x,[],'all'),max(disp.x,[],'all')];
plotlim.y=[min(disp.y,[],'all'),max(disp.y,[],'all')];

%% Choose Correction Options
quest='Correct Disp fields?';
dlgtitle='Displacement correction';
DispCorr.Opt=questdlg(quest,dlgtitle,'Yes');
clear quest dlgtitle

switch DispCorr.Opt
    case 'Yes'
        quest='How many pixels to correct?';
        DispCorr.int=str2num(cell2mat(inputdlg(quest)));%Correction interverval
        clear quest
        
        quest='Correction Method';
        DispCorr.Method=questdlg(quest,'Method',...
            'Direct',... %sets displacement equal to last valid value
            'LinGrad',... %follows the average linear gradient of the last
            ...           %Valid pitch (additional interpolation functions
            ...               %will be added if needed
           'Cancel',...    %Cancels the correction
           'Direct');    %Defaults to direct method
end

switch DispCorr.Method
    case 'Cancel'
        DispCorr.Opt='No';
end

switch DispCorr.Method
    case 'LinGrad'
        prompt='Number of Grid Pirches to Calculate Gradient';
        DispCorr.PitchFitKern=str2num(cell2mat(inputdlg(prompt)));
end

fprintf('--------------------------------------------------------------\n')
%% Perfom corrections if needed
switch    DispCorr.Opt
    case 'Yes'
        fprintf('Saving Raw Displacement Fields \n')
        RawDisp=disp;
        fprintf('Correcting Grid Method Displacements along specimen edges \n')
        [Corrdisp,DispCorr,grid,ProgramVersions]=func_CorrectGMDispV2(disp,...
            DispCorr,grid,pos);
        
        
    case 'No'
        fprintf('No Displacement corrections performed \n')
end

%% Save data
SaveFile=strcat(SavePath,'/',Desig,'Disp_Corr_Eval.mat');
save(SaveFile);


%% Generate plots
PngPath=strcat(SavePath,'/','PNG');
FigPath=strcat(SavePath,'/','FIG');

if exist(PngPath,'dir')==0
    mkdir(PngPath);
end

if exist('FigPath','dir')==0
    mkdir(FigPath)
end

figure('units','normalized','outerposition',[0 0 1 1])
for k=1:length(timemic)
    timestr=num2str(timemic(k));
    
    FEx=squeeze(FE.disp.x(:,:,k));
    FEy=squeeze(FE.disp.y(:,:,k));

    GMx=squeeze(disp.x(:,:,k));
    GMy=squeeze(disp.y(:,:,k));

    Corrx=squeeze(Corrdisp.x(:,:,k));
    Corry=squeeze(Corrdisp.y(:,:,k));

    % FE x displacements
    subplot(3,2,1)
    imagesc(FE.Xmm,FE.Ymm,FEx)
    xlabel('X Coordinate (mm)')
    ylabel('Y Coordinate (mm)')
    title(strcat('Finite Element u_x at',timestr,'\mu s'))
    colormap('jet')
    cx=colorbar;
    cx.Label.String='u_x (m)';
    caxis(FE.plotlim.x)
    set(gca,'YDir','normal')

    %FE y displacements
    subplot(3,2,2)
    imagesc(FE.Xmm,FE.Ymm,FEy)
    xlabel('X Coordinate (mm)')
    ylabel('Y Coordinate (mm)')
    title(strcat('Finite Element u_y at',timestr,'\mu s'))
    colormap('jet')
    cx=colorbar;
    cx.Label.String='u_y (m)';
    caxis(FE.plotlim.y)
    set(gca,'YDir','normal')
    
    %Raw Grid Method x
    subplot(3,2,3)
    imagesc(FE.Xmm,FE.Ymm,GMx)
    xlabel('X Coordinate (mm)')
    ylabel('Y Coordinate (mm)')
    title('Raw Grid Method u_x')
    colormap('jet')
    cx=colorbar;
    cx.Label.String='u_x (m)';
    caxis(plotlim.x)
    set(gca,'YDir','normal')

    % Raw Grid y displacements
    subplot(3,2,4)
    imagesc(Xmm,Ymm,GMy)
    xlabel('X Coordinate (mm)')
    ylabel('Y Coordinate (mm)')
    title('Raw Grid Method u_y')
    colormap('jet')
    cx=colorbar;
    cx.Label.String='u_y (m)';
    caxis(plotlim.y)
    set(gca,'YDir','normal')
    
     %Corrected Grid Method x
    subplot(3,2,5)
    imagesc(Xmm,Ymm,Corrx)
    xlabel('X Coordinate (mm)')
    ylabel('Y Coordinate (mm)')
    title('Corrected u_x')
    colormap('jet')
    cx=colorbar;
    cx.Label.String='u_x (m)';
    caxis(plotlim.x)
    set(gca,'YDir','normal')

    %Corrected y displacements
    subplot(3,2,6)
    imagesc(Xmm,Ymm,Corry)
    xlabel('X Coordinate (mm)')
    ylabel('Y Coordinate (mm)')
    title('Corrected u_y')
    colormap('jet')
    cx=colorbar;
    cx.Label.String='u_y (m)';
    caxis(plotlim.y)
    set(gca,'YDir','normal')

    %Save
    FigName=strcat(FigPath,'/',Desig,'_',timestr,'_DispComp.fig');
    saveas(gcf,FigName);

    PngName=strcat(PngPath,'/',Desig,'_',timestr,'_DispComp.png');
    saveas(gcf,PngName);
     

end
%% Generate Difference plots
CorrDiff.x=disp.x-Corrdisp.x;
CorrDiff.y=disp.y-Corrdisp.y;
Difflim.x=[min(CorrDiff.x,[],'all'),max(CorrDiff.x,[],'all')];

if Difflim.x==[0,0]
    Difflim.x=[-1,1];
end

Difflim.y=[-mean(CorrDiff.y,'all'),mean(CorrDiff.y,'all')];


PngPath=strcat(SavePath,'/Diff/','PNG');
FigPath=strcat(SavePath,'/Diff/','FIG');

if exist(PngPath,'dir')==0
    mkdir(PngPath);
end

if exist('FigPath','dir')==0
    mkdir(FigPath)
end

figure('units','normalized','outerposition',[0 0 1 1])
for k=1:length(timemic)
    timestr=num2str(timemic(k));
    
    FEx=squeeze(FE.disp.x(:,:,k));
    FEy=squeeze(FE.disp.y(:,:,k));

    GMx=squeeze(disp.x(:,:,k));
    GMy=squeeze(disp.y(:,:,k));

    Corrx=squeeze(Corrdisp.x(:,:,k));
    Corry=squeeze(Corrdisp.y(:,:,k));

    Diffx=squeeze(CorrDiff.x(:,:,k));
    Diffy=squeeze(CorrDiff.y(:,:,k));

   
    %Raw Grid Method x
    subplot(3,2,1)
    imagesc(FE.Xmm,FE.Ymm,GMx)
    xlabel('X Coordinate (mm)')
    ylabel('Y Coordinate (mm)')
    title(strcat('Raw Grid Method u_x ',timestr,'\mu s'))
    colormap('jet')
    cx=colorbar;
    cx.Label.String='u_x (m)';
    caxis(plotlim.x)
    set(gca,'YDir','normal')

    % Raw Grid y displacements
    subplot(3,2,2)
    imagesc(Xmm,Ymm,GMy)
    xlabel('X Coordinate (mm)')
    ylabel('Y Coordinate (mm)')
    title(strcat('Raw Grid Method u_y ',timestr,'\mu s'))
    colormap('jet')
    cx=colorbar;
    cx.Label.String='u_y (m)';
    caxis(plotlim.y)
    set(gca,'YDir','normal')
    
    %Corrected Grid Method x
    subplot(3,2,3)
    imagesc(Xmm,Ymm,Corrx)
    xlabel('X Coordinate (mm)')
    ylabel('Y Coordinate (mm)')
    title('Corrected u_x')
    colormap('jet')
    cx=colorbar;
    cx.Label.String='u_x (m)';
    caxis(plotlim.x)
set(gca,'YDir','normal')

    %Corrected y displacements
    subplot(3,2,4)
    imagesc(Xmm,Ymm,Corry)
    xlabel('X Coordinate (mm)')
    ylabel('Y Coordinate (mm)')
    title('Corrected u_y')
    colormap('jet')
    cx=colorbar;
    cx.Label.String='u_y (m)';
    caxis(plotlim.y)
    set(gca,'YDir','normal')


    %Correction Difference Grid Method x
    subplot(3,2,5)
    imagesc(Xmm,Ymm,Diffx)
    xlabel('X Coordinate (mm)')
    ylabel('Y Coordinate (mm)')
    title('Raw - Corrected u_x')
    colormap('jet')
    cx=colorbar;
    cx.Label.String='u_x (m)';
    caxis(Difflim.x)
    set(gca,'YDir','normal')

    %Correction Difference y displacements
    subplot(3,2,6)
    imagesc(Xmm,Ymm,Diffy)
    xlabel('X Coordinate (mm)')
    ylabel('Y Coordinate (mm)')
    title('Raw - Corrected u_y')
    colormap('jet')
    cx=colorbar;
    cx.Label.String='u_y (m)';
    caxis(Difflim.y)
    set(gca,'YDir','normal')

    %Save
    FigName=strcat(FigPath,'/',Desig,'_',timestr,'_DispDiff.fig');
    saveas(gcf,FigName);

    PngName=strcat(PngPath,'/',Desig,'_',timestr,'_DispDiff.png');
    saveas(gcf,PngName);
     

end

%% Y Diff Zoomed

PngPath=strcat(SavePath,'/YDiff/','PNG');
FigPath=strcat(SavePath,'/YDiff/','FIG');

if exist(PngPath,'dir')==0
    mkdir(PngPath);
end

if exist('FigPath','dir')==0
    mkdir(FigPath)
end
figure('units','normalized','outerposition',[0 0 1 1])

for k=1:length(timemic)
    timestr=num2str(timemic(k));
    
    Diffy=squeeze(CorrDiff.y(:,:,k));

    imagesc(Xmm,Ymm,Diffy)
    xlabel('X Coordinate (mm)')
    ylabel('Y Coordinate (mm)')
   title(strcat('Raw - Corrected u_y ',timestr,'\mu s'))
    colormap('hot')
    cx=colorbar;
    cx.Label.String='u_y (m)';
    caxis(Difflim.y)
    set(gca,'YDir','normal')

    %Save
    FigName=strcat(FigPath,'/',Desig,'_',timestr,'_YDiff.fig');
    saveas(gcf,FigName);

    PngName=strcat(PngPath,'/',Desig,'_',timestr,'_YDiff.png');
    saveas(gcf,PngName);
     

end