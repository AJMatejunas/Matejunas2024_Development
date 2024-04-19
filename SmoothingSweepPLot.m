%% Load The Error information
[ErrorFile,ErrorPath]=uigetfile('*.mat',...
    'load File containing Identification error')
load(strcat(ErrorPath,'/',ErrorFile))
savePath=ErrorPath;

%% Define Test Designation
TestDeg=char(inputdlg('Choose Test Designation'));
%% Choose Colormap
defaultScheme={'Heat'};
prompt='Choose Color Map'
ColorScheme=char(inputdlg(prompt,'color scheme',[1,25],defaultScheme));

%% Calculate Minimun Errors

% minKerror=min(abs(squeeze(Errors.K(:,:,1)),[],'all'));
%[minKtempIn,MinKspaIn]=find(abs(E

%% Plot Heatmaps for K_1

figure('units','Centimeters','outerposition',[0 0 17 17])
imagesc(TempKernVec,SpaKernVec,squeeze(Errors.K(:,:,1)))
title('K_1 identification Errors')
xlabel('Temporal Smoothing Kernal')
ylabel('Spatial Smoothing Kernal')
colormap(ColorScheme)
cx=colorbar;
cx.Label.String='Error (%)';
set(gca,'YDir','normal')
% save
saveas(gcf,strcat(savePath,'/',TestDeg,...
        '_CombinedSpatialTemporal_Kerrors',ColorScheme,'.fig'));
    saveas(gcf,strcat(savePath,'/',TestDeg,...
        '_CombinedSpatialTemporal_Kerrors',ColorScheme,'.svg'));
   saveas(gcf,strcat(savePath,'/',TestDeg,...
        '_CombinedSpatialTemporal_Kerrors',ColorScheme,'.png'));
    


%% Plot Heatmaps for G_1
figure('units','Centimeters','outerposition',[0 0 17 17])
imagesc(TempKernVec,SpaKernVec,squeeze(Errors.G(:,:)))
title('G_1 identification Errors')
xlabel('Temporal Smoothing Kernal')
ylabel('Spatial Smoothing Kernal')
colormap(ColorScheme)
cx=colorbar;
cx.Label.String='Error (%)';
set(gca,'YDir','normal')
   
% save
saveas(gcf,strcat(savePath,'/',TestDeg,...
        '_CombinedSpatialTemporal_Gerrors',ColorScheme,'.fig'));
    saveas(gcf,strcat(savePath,'/',TestDeg,...
        '_CombinedSpatialTemporal_Gerrors',ColorScheme,'.svg'));
   saveas(gcf,strcat(savePath,'/',TestDeg,...
        '_CombinedSpatialTemporal_Gerrors',ColorScheme,'.png'));


%% Plot Heatmaps for tau_1

figure('units','Centimeters','outerposition',[0 0 17 17])
imagesc(TempKernVec,SpaKernVec,squeeze(Errors.tau(:,:,1)))
title('\tau_1 identification Errors')
xlabel('Temporal Smoothing Kernal')
ylabel('Spatial Smoothing Kernal')
colormap(ColorScheme)
cx=colorbar;
cx.Label.String='Error (%)';
set(gca,'YDir','normal')
% save
saveas(gcf,strcat(savePath,'/',TestDeg,...
        '_CombinedSpatialTemporal_Tauerrors',ColorScheme,'.fig'));
    saveas(gcf,strcat(savePath,'/',TestDeg,...
        '_CombinedSpatialTemporal_Tauerrors',ColorScheme,'.svg'));
   saveas(gcf,strcat(savePath,'/',TestDeg,...
        '_CombinedSpatialTemporal_Tauerrors',ColorScheme,'.png'));

%% calculate total Error 
Errors.E=sqrt(Errors.K.^2+Errors.G.^2+Errors.tau.^2);

figure('units','Centimeters','outerposition',[0 0 17 17])
imagesc(TempKernVec,SpaKernVec,squeeze(Errors.E(:,:,1)))
title('Resultant identification Errors')
xlabel('Temporal Smoothing Kernal')
ylabel('Spatial Smoothing Kernal')
colormap(ColorScheme)
cx=colorbar;
cx.Label.String='Error (%)';
set(gca,'YDir','normal')
% save
saveas(gcf,strcat(savePath,'/',TestDeg,...
        '_CombinedSpatialTemporal_TotalErrors',ColorScheme,'.fig'));
    saveas(gcf,strcat(savePath,'/',TestDeg,...
        '_CombinedSpatialTemporal_TotalErrors',ColorScheme,'.svg'));
   saveas(gcf,strcat(savePath,'/',TestDeg,...
        '_CombinedSpatialTemporal_TotalErrors',ColorScheme,'.png'));
