%% This script is written to check the difference between evaluations
    %Of func_ViscoConstitutiveV4 & V5
    
    %initialize
    fprintf('initializing \n')
    
    clear all, close all, clc
    
%% select data files to load

%Constitutive model V4 
[Oldfile,OldPath]=uigetfile('*.mat',...
    'Choose file with original evaluation of the constitutive model');

%constitutive model new file
[Newfile,NewPath]=uigetfile('*.mat',...
    'Choose file with new evaluation fo the constitutive model');

%Geometric data
[Gfile,Gpath]=uigetfile('*.mat',...
    'Choose file containing geometric information \n');
%% define test designation
testdeg=char(inputdlg('Input test designation'));

%% Record versions of the model code
Versions=inputdlg({'new version','0ld version'},'Version definition',...
    [1,50;1,50],{'V4','V5'});
Vold=cell2mat(Versions(1));
Vnew=cell2mat(Versions(2));

clear versions


%% load the files
fprintf('Loading old data file \n')
Old=load(strcat(OldPath,'\',Oldfile));

fprintf('Loading new data file \n')
New=load(strcat(NewPath,'\',Newfile));

fprintf('Loading Geometry File \n')
load(strcat(Gpath,'\',Gfile),'X_vec','Y_vec');

fprintf('Loading Complete \n')

%% Create plots of the the Poisson's Ratio
fprintf('Comparing Poisons Ratio evolution \n')
figure(1)
told=Old.time.vec*10^6;
tnew=New.time.vec*10^6;

NuOld=Old.ConstStress.nuR;
NuNew=New.ConstStress.nuR;

plot(told,NuOld,tnew,NuNew)
legend(Vold,Vnew)
xlabel('time (\mus)')
ylabel('\nu')

saveas(gcf,strcat(testdeg,'-',Vold,'Vs',Vnew,'_Nu.fig'))
saveas(gcf,strcat(testdeg,'-',Vold,'Vs',Vnew,'_Nu.png'))

%% Create time and geometry info for plotting
timemic=New.time.vec*10^6; %us

Xcoord=squeeze(X_vec*1000)';  %mm
Ycoord=squeeze(Y_vec*1000);   %mm

%% Compare the hydrostatic stresses and strains
HstrainDiff=Old.ConstStress.hstrain-New.ConstStress.hstrain;
HstressDiff=Old.ConstStress.hy-New.ConstStress.hy;

HstrainDlim=[min(HstrainDiff,[],'all'),max(HstrainDiff,[],'all')];
HstressDlim=[min(HstressDiff,[],'all'),max(HstressDiff,[],'all')]*10^-6;

%create limits
    %determine maximum and minimu stresses and strains
    MaxHstress=max([New.ConstStress.hy,Old.ConstStress.hy],[],'all');
    MinHstress=min([New.ConstStress.hy,Old.ConstStress.hy],[],'all');
    
    MaxHstrain=max([New.ConstStress.hstrain,Old.ConstStress.hstrain],...
        [],'all');
    MinHstrain=min([New.ConstStress.hstrain,Old.ConstStress.hstrain],...
        [],'all');
    
    HstrainLim=[MinHstrain,MaxHstrain];
    HstressLim=[MinHstress,MaxHstress]*10^-6;
    


% Plot
fprintf('Plotting Hydrostatic comparison \n')

%create directory for images
foldername=strcat('Hydro',Vold,Vnew);
mkdir(foldername);

figure('units','normalized','outerposition',[0 0 1 1])
for k=1:length(timemic)
   
 timecount=num2str(timemic(k)); 
 framecount=num2str(k);
 
 NHstress=squeeze(New.ConstStress.hy(:,:,k))*10^-6;
 NHstrain=squeeze(New.ConstStress.hstrain(:,:,k));
 OHstress=squeeze(Old.ConstStress.hy(:,:,k))*10^-6;
 OHstrain=squeeze(Old.ConstStress.hstrain(:,:,k));
    
 DHstress=squeeze(HstressDiff(:,:,k))*10^-6;
 DHstrain=squeeze(HstrainDiff(:,:,k));

% Original stress
subplot(3,2,1)
   imagesc(Xcoord,Ycoord,OHstress)
   title(strcat(Vold,' HStress ',timecount,'\mu s'))
   xlabel('X coordinate (mm)')
   ylabel('Y coordinate (mm)')
   colormap('hot')
   cx=colorbar;
   cx.Label.String='Hydrostatic Stress(MPa)';
   caxis(HstressLim)

   %Original Strain
subplot(3,2,2)
   imagesc(Xcoord,Ycoord,OHstrain)
   title(strcat(Vold,' HStrain ',timecount,'\mu s'))
   xlabel('X coordinate (mm)')
   ylabel('Y coordinate (mm)')
   colormap('hot')
   cx=colorbar;
   cx.Label.String='Hydrostatic Strain';
   caxis(HstrainLim)

%New Stress
subplot(3,2,3)
   imagesc(Xcoord,Ycoord,NHstress)
   title(strcat(Vnew,' HStress ',timecount,'\mu s'))
   xlabel('X coordinate (mm)')
   ylabel('Y coordinate (mm)')
   colormap('hot')
   cx=colorbar;
   cx.Label.String='Hydrostatic Stress(MPa)';
   caxis(HstressLim)
   
%New Strain
subplot(3,2,4)
   imagesc(Xcoord,Ycoord,NHstrain)
   title(strcat(Vnew,' HStrain ',timecount,'\mu s'))
   xlabel('X coordinate (mm)')
   ylabel('Y coordinate (mm)')
   colormap('hot')
   cx=colorbar;
   cx.Label.String='Hydrostatic Strain';
   caxis(HstrainLim)
   
%Stress difference
subplot(3,2,5)
   imagesc(Xcoord,Ycoord,DHstress)
   title(strcat(Vold,'-',Vnew,' HStress ',timecount,'\mu s'))
   xlabel('X coordinate (mm)')
   ylabel('Y coordinate (mm)')
   colormap('hot')
   cx=colorbar;
   cx.Label.String='Hydrostatic Stress(MPa)';
   caxis(HstressDlim)

%Strain difference
subplot(3,2,6)
   imagesc(Xcoord,Ycoord,DHstrain)
   title(strcat(Vold,'-',Vnew,' HStrain ',timecount,'\mu s'))
   xlabel('X coordinate (mm)')
   ylabel('Y coordinate (mm)')
   colormap('hot')
   cx=colorbar;
   cx.Label.String='Hydrostatic Strain';
   caxis(HstrainDlim)
   
saveas(gcf,strcat(pwd,'\',foldername,'\Hydro_frame',framecount,'.fig'));
saveas(gcf,strcat(pwd,'\',foldername,'\Hydro_frame',framecount,'.png'));
      
end
%% Compare Compare shear stresses

fprintf('Plotting Shear Comparison \n')

ShearDiff=Old.ConstStress.xy-New.ConstStress.xy;

ShearDlim=[min(ShearDiff,[],'all'),max(ShearDiff,[],'all')]*10^-6;

MaxShear=max([New.ConstStress.xy,Old.ConstStress.xy],[],'all');
MinShear=min([New.ConstStress.xy,Old.ConstStress.xy],[],'all');
ShearLim=[MinShear,MaxShear]*10^-6;

%make directory
foldername=strcat('Shear',Vold,Vnew);
mkdir(foldername);

figure('units','normalized','outerposition',[0 0 1 1])
for k=1:length(timemic)
 timecount=num2str(timemic(k)); 
 framecount=num2str(k);
 
 NSstress=squeeze(New.ConstStress.xy(:,:,k))*10^-6;
 OSstress=squeeze(Old.ConstStress.xy(:,:,k))*10^-6;
 DSstress=squeeze(ShearDiff(:,:,k))*10^-6;
 
 %Old Shear Stress
 subplot(3,1,1)
 imagesc(Xcoord,Ycoord,OSstress)
   title(strcat(Vold,' ShearStress ',timecount,'\mu s'))
   xlabel('X coordinate (mm)')
   ylabel('Y coordinate (mm)')
   colormap('hot')
   cx=colorbar;
   cx.Label.String='Shear Stress(MPa)';
   caxis(ShearLim)

 %New Shear Stress
 subplot(3,1,2)
 imagesc(Xcoord,Ycoord,NSstress)
   title(strcat(Vnew,' ShearStress ',timecount,'\mu s'))
   xlabel('X coordinate (mm)')
   ylabel('Y coordinate (mm)')
   colormap('hot')
   cx=colorbar;
   cx.Label.String='Shear Stress(MPa)';
   caxis(ShearLim)
   
%Shear difference
 subplot(3,1,3)
 imagesc(Xcoord,Ycoord,DSstress)
   title(strcat(Vold,'-',' ShearStress ',timecount,'\mu s'))
   xlabel('X coordinate (mm)')
   ylabel('Y coordinate (mm)')
   colormap('hot')
   cx=colorbar;
   cx.Label.String='Shear Stress(MPa)';
   caxis(ShearDlim)
   
saveas(gcf,strcat(pwd,'\',foldername,'\Shear_frame',framecount,'.fig'));
saveas(gcf,strcat(pwd,'\',foldername,'\Shear_frame',framecount,'.png'));
end
%% Compare stresses in X direction
fprintf('Plotting XX stress Comparison \n')


XDiff=Old.ConstStress.xx-New.ConstStress.xx;

XDlim=[min(XDiff,[],'all'),max(XDiff,[],'all')]*10^-6;

MaxX=max([New.ConstStress.xx,Old.ConstStress.xx],[],'all');
MinX=min([New.ConstStress.xx,Old.ConstStress.xx],[],'all');
XLim=[MinX,MaxX]*10^-6;

%make directory
foldername=strcat('XX',Vold,Vnew);
mkdir(foldername);

figure('units','normalized','outerposition',[0 0 1 1])
for k=1:length(timemic)
 timecount=num2str(timemic(k)); 
 framecount=num2str(k);
 
 NXstress=squeeze(New.ConstStress.xx(:,:,k))*10^-6;
 OXstress=squeeze(Old.ConstStress.xx(:,:,k))*10^-6;
 DXstress=squeeze(XDiff(:,:,k))*10^-6;
 
 %Old xx Stress
 subplot(3,1,1)
 imagesc(Xcoord,Ycoord,OXstress)
   title(strcat(Vold,' XX Stress ',timecount,'\mu s'))
   xlabel('X coordinate (mm)')
   ylabel('Y coordinate (mm)')
   colormap('hot')
   cx=colorbar;
   cx.Label.String='XX Stress(MPa)';
   caxis(XLim)

 %New xx Stress
 subplot(3,1,2)
 imagesc(Xcoord,Ycoord,NXstress)
   title(strcat(Vnew,' XXStress ',timecount,'\mu s'))
   xlabel('X coordinate (mm)')
   ylabel('Y coordinate (mm)')
   colormap('hot')
   cx=colorbar;
   cx.Label.String='XX Stress(MPa)';
   caxis(XLim)
   
%xx difference
 subplot(3,1,3)
 imagesc(Xcoord,Ycoord,DXstress)
   title(strcat(Vold,'-',' XX Stress ',timecount,'\mu s'))
   xlabel('X coordinate (mm)')
   ylabel('Y coordinate (mm)')
   colormap('hot')
   cx=colorbar;
   cx.Label.String='XX Stress(MPa)';
   caxis(XDlim)
   
saveas(gcf,strcat(pwd,'\',foldername,'\XX_frame',framecount,'.fig'));
saveas(gcf,strcat(pwd,'\',foldername,'\XX_frame',framecount,'.png'));
end

%% Save all data
fprintf('Saving Data \n')
save(strcat(testdeg,'_',Vold,'vs',Vnew,'.mat'))

fprintf('Done \n')