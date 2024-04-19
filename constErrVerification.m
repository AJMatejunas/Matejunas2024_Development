% This code is written to evaluate the error in the constitutive model
    % calculation and produce images for a video of the error propagation
    % with time and space
    
  
    %% initialize
    close all, clear all, clc
    
    %% load data
    datafile=uigetfile('.mat','choose data file');
    propfile=uigetfile('.mat','load material properties file');
    [constfile,constpath]=uigetfile('.mat',...
        'choose file containing constitutive stress');
    
    desig=char(inputdlg('input test designation'));  
    
    fprintf('Loading Data Files \n')
    load(datafile);
    clear MatProps
    load(propfile);
    load(strcat(constpath,'/',constfile),'StressModel');
    testdeg=desig;
    clear desig
    %% extract X and Y vectors
    X_vec=Xq(1,:);
    Y_vec=Yq(:,1);
%     %% calculate constitutive model
%     fprintf('Calculating Constitutive model \n')
%    % get properties into the proper form for inputs into the constitutive
%         % model algorithm
%    MatProps.nu=.26;
%    nu=MatProps.nu;
%    MatProps.Einf=0;
%    MatProps.E0=0;
%    
%    Einf=2.98e9;
%    E1=2.21e9;
% 
%    
%    tau=10e-6;
%    
%    
%    
%    K=MatProps.Ki;
%    G=MatProps.Gi;
%    
%    MatProps.K0=MatProps.Kinf+K;
%    MatProps.G0=MatProps.Ginf+G;
%    
%    StressModel=sfunc_ViscoConstitutiveV5(strain.x,strain.y,strain.s,...
%        time.vec,MatProps,K,G,tau);
%    
%        
%     
    %% calculate error
    fprintf('Calculating Error')
    [stressXXerr,stressXYerr,stressYYerr,...
        XXdiff,YYdiff,Sheardiff]=func_constErrV2(stress,StressModel);
    
  %% convert into proper units
  %time
  timemic=time.vec*10^6; %us
  
  %geometric coordinates
  Xcoord=squeeze(X_vec*1000)';  %mm
  Ycoord=squeeze(Y_vec*1000);   %mm
  
  %% Develop limits for error images
  
  %limits error to errors that are meaningful magnitudes (the goal is to
  %drop this to +- 10%
    ErrLim=[-100,100]; %limit the errors to +-100%
  
  %minimum and maximum strains
  MinXstrain=min(strain.x,[],'all');
  MinXYstrain=min(strain.s,[],'all');
  MaxXstrain=max(strain.x,[],'all');
  MaxXYstrain=max(strain.s,[],'all');
 
  XstrainLim=[MinXstrain,MaxXstrain];
  XYstrainLim=[MinXYstrain,MaxXYstrain];
  
  %minimum and maximum Abaqus Stresses
  MinXstress=min(stress.x,[],'all');
  MinXYstress=min(stress.s,[],'all');
  MaxXstress=max(stress.x,[],'all');
  MaxXYstress=max(stress.s,[],'all');
  
  XstressLim=[MinXstress,MaxXstress]/10^6;
  XYstressLim=[MinXYstress,MaxXYstress]/10^6;
  
  %Constitutive limits
  ConstLim=[min(StressModel.xy,[],'all'),max(StressModel.xy,[],'all')]/10^6;
  
  %acceleration limits
  MinAy=min(accel.y,[],'all');
  MaxAy=max(accel.y,[],'all');
  
  AyLim=[MinAy,MaxAy];
 
  %difference limits
  XXDlim=[min(XXdiff,[],'all'),max(XXdiff,[],'all')]/10^6;
  YYDlim=[min(YYdiff,[],'all'),max(YYdiff,[],'all')]/10^6;
  ShearDlim=[min(Sheardiff,[],'all'),max(Sheardiff,[],'all')]/10^6;
  
% %% Produce error images
%   
% for k=1:length(time.vec)
%     Xerr=squeeze(stressXXerr(:,:,k));
%     Yerr=squeeze(stressYYerr(:,:,k));
%     XYerr=squeeze(stressXYerr(:,:,k));
%     timecount=num2str(timemic(k)); 
%     framecount=num2str(k);
%     
%     Xstress=squeeze(stress.x(:,:,k))/10^6;
%     XYstress=squeeze(stress.s(:,:,k))/10^6;
%     
%     XYstrain=squeeze(strain.s(:,:,k));
%     Ay=squeeze(accel.y(:,:,k));
%     
%     %% Raw images
% %     figure(1)
% %     imagesc(Xcoord,Ycoord,Xerr)
% %     title(strcat('\sigma_{xx} ',timecount,'\mu s'))
% %     xlabel('X coordinate (mm)')
% %     ylabel('Y coordinate (mm)')
% %     cx=colorbar;
% %     cx.Label.String='error (%)';
% %     caxis(ErrLim)
% %     
% %     saveas(gcf,strcat(testdeg,'_Xerror_frame ',framecount,'.fig'))
% %     saveas(gcf,strcat(testdeg,'_Xerror_frame ',framecount,'.png'))
% %     
% %     figure(2)
% %     imagesc(Xcoord,Ycoord,XYerr)
% %     title(strcat('shear ',timecount,'\mu s'))
% %     xlabel('X coordinate (mm)')
% %     ylabel('Y coordinate (mm)')
% %     cx=colorbar;
% %     cx.Label.String='error (%)';
% %     caxis(ErrLim)
% %     
% %     saveas(gcf,strcat(testdeg,'_Shearerror_frame ',framecount,'.fig'))
% %     saveas(gcf,strcat(testdeg,'_Shearerror_frame ',framecount,'.png'))
% %     
% %     figure(3)
% %     imagesc(Xcoord,Ycoord,Yerr)
% %     title(strcat('\sigma_{xx}',timecount,'\mu s'))
% %     xlabel('X coordinate (mm)')
% %     ylabel('Y coordinate (mm)')
% %     cx=colorbar;
% %     cx.Label.String='error (%)';
% %     caxis(ErrLim)
% %     
% %     saveas(gcf,strcat(testdeg,'_Yerror_frame ',framecount,'.fig'))
% %     saveas(gcf,strcat(testdeg,'_Yerror_frame ',framecount,'.png'))
% %     
%     %% Compare error to stress magnitude and stresses
%     figure(4)
%     % Upper left error in XX stress
%     subplot(2,2,2)
%     imagesc(Xcoord,Ycoord,Xerr)
%     title(strcat('\sigma_{xx} error',timecount,'\mu s'))
%     xlabel('X coordinate (mm)')
%     ylabel('Y coordinate (mm)')
%     cx=colorbar;
%     cx.Label.String='error (%)';
%     caxis(ErrLim)
%     
%     % upper right. Magnitude of the stress
%     subplot(2,2,2)
%     imagesc(Xcoord,Ycoord,Xstress)
%     title(strcat('\sigma_{xx} ',timecount,'\mu s'))
%     xlabel('X coordinate (mm)')
%     ylabel('Y coordinate (mm)')
%     cx=colorbar;
%     cx.Label.String='Stress XX (MPa)';
%     caxis(XstressLim)
%     
%     %lower left error in shear stress
%     subplot(2,2,3)
%     imagesc(Xcoord,Ycoord,XYerr)
%     title(strcat('shear Error',timecount,'\mu s'))
%     xlabel('X coordinate (mm)')
%     ylabel('Y coordinate (mm)')
%     cx=colorbar;
%     cx.Label.String='error (%)';
%     caxis(ErrLim)
%     
%     %lower right magnitude of shear stresses
%     % upper right. Magnitude of the stress
%     subplot(2,2,4)
%     imagesc(Xcoord,Ycoord,XYstress)
%     title(strcat('\sigma_{xy} ',timecount,'\mu s'))
%     xlabel('X coordinate (mm)')
%     ylabel('Y coordinate (mm)')
%     cx=colorbar;
%     cx.Label.String='Shear Stress(MPa)';
%     caxis(XYstressLim)
%    
%     
%     saveas(gcf,strcat(testdeg,'_errors_',framecount,'.fig'))
%     saveas(gcf,strcat(testdeg,'_errors_',framecount,'.png'))
%     
%    %% closer look at shear with accelerations and strain magnitudes
%    figure(5)
%    %upper left, shear error
%    subplot(2,2,1)
%     imagesc(Xcoord,Ycoord,abs(XYerr))
%     title(strcat('shear Error',timecount,'\mu s'))
%     xlabel('X coordinate (mm)')
%     ylabel('Y coordinate (mm)')
%     colormap('hot')
%     cx=colorbar;
%     cx.Label.String='error (%)';
%     caxis([0,max(abs(stressXYerr),[],'all')])
%    
%    %upper right acceleration in the Y-direction
%    subplot(2,2,2)
%    imagesc(Xcoord,Ycoord,Ay)
%    title(strcat('\a_{y} ',timecount,'\mu s'))
%    xlabel('X coordinate (mm)')
%    ylabel('Y coordinate (mm)')
%    colormap('hot')
%    cx=colorbar;
%    cx.Label.String='ay (m/s)';
%    caxis(AyLim)
%    
%    %bottom left, shear strain magnitude
%    subplot(2,2,3)
%    imagesc(Xcoord,Ycoord,XYstrain)
%    title(strcat('\epsilon_{xy} ',timecount,'\mu s'))
%    xlabel('X coordinate (mm)')
%    ylabel('Y coordinate (mm)')
%    colormap('hot')
%    cx=colorbar;
%    cx.Label.String='Shear Strain';
%    caxis(XYstrainLim)
%    
%    %bottom right, shear stress magnitude
%    subplot(2,2,4)
%    imagesc(Xcoord,Ycoord,XYstress)
%    title(strcat('\sigma_{xy} ',timecount,'\mu s'))
%    xlabel('X coordinate (mm)')
%    ylabel('Y coordinate (mm)')
%    colormap('hot')
%    cx=colorbar;
%    cx.Label.String='Shear Stress(MPa)';
%    caxis(XYstressLim)
%    
%    saveas(gcf,strcat(testdeg,'_heatPlots_',framecount,'.fig'))
%    saveas(gcf,strcat(testdeg,'_heatPlots_',framecount,'.png'))
% end
% 
% %% plot detailed shear error
% ConstLim=[min(StressModel.xy,[],'all'),max(StressModel.xy,[],'all')]/10^6;
% 
% figure('units','normalized','outerposition',[0 0 1 1])
% for k=1:length(time.vec)
%  absXYerr=abs(squeeze(stressXYerr(:,:,k)));
%  timecount=num2str(timemic(k)); 
%  framecount=num2str(k);
%  XYstress=squeeze(stress.s(:,:,k))/10^6;
%  ConstXY=squeeze(StressModel.xy(:,:,k))/10^6;
%  logerr=log10(absXYerr);
%  logerr(absXYerr==0)=0;
% 
% 
%  
%  % Stress output by abaqus
%  subplot(3,2,1)
%    imagesc(Xcoord,Ycoord,XYstress)
%    title(strcat('ABQ \sigma_{xy} ',timecount,'\mu s'))
%    xlabel('X coordinate (mm)')
%    ylabel('Y coordinate (mm)')
%    colormap('hot')
%    cx=colorbar;
%    cx.Label.String='Shear Stress(MPa)';
%    caxis(XYstressLim)
%    
% % constitutive stress 
% subplot(3,2,2)
%    imagesc(Xcoord,Ycoord,ConstXY)
%    title(strcat('Const \sigma_{xy} ',timecount,'\mu s'))
%    xlabel('X coordinate (mm)')
%    ylabel('Y coordinate (mm)')
%    colormap('hot')
%    cx=colorbar;
%    cx.Label.String='Shear Stress(MPa)';
%    caxis(ConstLim)   
% 
% %absolute value of error
% subplot(3,2,3)
%     imagesc(Xcoord,Ycoord,absXYerr)
%     title(strcat('abs(shear Error)',timecount,'\mu s'))
%     xlabel('X coordinate (mm)')
%     ylabel('Y coordinate (mm)')
%     colormap('hot')
%     cx=colorbar;
%     cx.Label.String='error (%)';
%     caxis([0,max(abs(stressXYerr),[],'all')])
% 
% %absolute value of the error logarithmicly
% subplot(3,2,4)
%     imagesc(Xcoord,Ycoord,logerr)
%     title(strcat('log(abs(shear Error))',timecount,'\mu s'))
%     xlabel('X coordinate (mm)')
%     ylabel('Y coordinate (mm)')
%     colormap('hot')
%     cx=colorbar;
%     cx.Label.String='error (%)';
%     caxis([0,log(max(abs(stressXYerr),[],'all'))])
% 
% %error between 0 and 100%
% subplot(3,2,5)
%     imagesc(Xcoord,Ycoord,absXYerr)
%     title(strcat('abs(shear Error)',timecount,'\mu s'))
%     xlabel('X coordinate (mm)')
%     ylabel('Y coordinate (mm)')
%     colormap('hot')
%     cx=colorbar;
%     cx.Label.String='error (%)';
%     caxis([0,100])
% % error between 0 and 10%
% subplot(3,2,6)
%     imagesc(Xcoord,Ycoord,absXYerr)
%     title(strcat('abs(shear Error)',timecount,'\mu s'))
%     xlabel('X coordinate (mm)')
%     ylabel('Y coordinate (mm)')
%     colormap('hot')
%     cx=colorbar;
%     cx.Label.String='error (%)';
%     caxis([0,10])
% 
% saveas(gcf,strcat(testdeg,'_HeatError_',framecount,'.fig'))
% saveas(gcf,strcat(testdeg,'_HeatError_',framecount,'.png'))
% end
% 
% % %% plot mean and median shear errors
% for m=1:length(time.vec)
%     absXYerr=abs(squeeze(stressXYerr(:,:,m)));
%     XYerr=squeeze(stressXYerr(:,:,m));
%     
%     absmeanErr(m)=mean(absXYerr,'all');
%     absmedianErr(m)=median(absXYerr,'all');
%     
%     meanErr(m)=mean(XYerr,'all');
%     medianErr(m)=median(XYerr,'all');
% end
% figure
% plot(time.vec*10^6,absmeanErr)
% hold on
% plot(time.vec*10^6,absmedianErr)
% hold off
% legend('abs(mean)','abs(median)')
% xlabel('time \mus')
% ylabel('error (%)')
% 
% saveas(gcf,strcat(testdeg,'_absAvgErr.fig'))
% saveas(gcf,strcat(testdeg,'_absAvgErr.png'))
% 
% figure
% plot(time.vec*10^6,absmedianErr)
% xlabel('time \mus')
% ylabel('error (%)')
% legend('abs(median)')
% saveas(gcf,strcat(testdeg,'_absMedianError.fig'))
% saveas(gcf,strcat(testdeg,'_absMedianError.png'))
% 
% figure
% plot(time.vec*10^6,meanErr)
% hold on
% plot(time.vec*10^6,medianErr)
% hold off
% legend('mean','median')
% xlabel('time \mus')
% ylabel('error (%)')
% 
% saveas(gcf,strcat(testdeg,'_AvgErr.fig'))
% saveas(gcf,strcat(testdeg,'_AvgErr.png'))
% 
% figure
% plot(time.vec*10^6,medianErr)
% xlabel('time \mus')
% ylabel('error (%)')
% legend('median')
% saveas(gcf,strcat(testdeg,'_MedianError.fig'))
% saveas(gcf,strcat(testdeg,'_MedianError.png'))
% 
% % Produce images of shear errors with differences
% 
% figure('units','normalized','outerposition',[0 0 1 1])
% for n=1:length(time.vec)
%     Xerr=squeeze(stressXXerr(:,:,n));
%     Yerr=squeeze(stressYYerr(:,:,n));
%     XYerr=squeeze(stressXYerr(:,:,n));
%     timecount=num2str(timemic(n)); 
%     framecount=num2str(n);
%     
%     ConstXX=squeeze(StressModel.xx(:,:,n))/10^6;
%     ConstXY=squeeze(StressModel.xy(:,:,n))/10^6;
%     Xstress=squeeze(stress.x(:,:,n))/10^6;
%     XYstress=squeeze(stress.s(:,:,n))/10^6;
% 
%     DiffX=squeeze(XXdiff(:,:,n))/10^6;
%     DiffXY=squeeze(Sheardiff(:,:,n))/10^6;
%      
% % Shear     
%      
%      
% %   Abaqus stress
%     subplot(2,2,1)
%     imagesc(Xcoord,Ycoord,XYstress)
%     title(strcat('ABQ \sigma_{xy} ',timecount,'\mu s'))
%     xlabel('X coordinate (mm)')
%     ylabel('Y coordinate (mm)')
%     colormap('hot')
%     cx=colorbar;
%     cx.Label.String='Shear Stress(MPa)';
%     caxis(XYstressLim)
%     
%  %   Constitutive Stress
%     subplot(2,2,2)
%     imagesc(Xcoord,Ycoord,ConstXY)
%     title(strcat('Const \sigma_{xy} ',timecount,'\mus'))
%     xlabel('X coordinate (mm)')
%     ylabel('Y coordinate (mm)')
%     colormap('hot')
%     cx=colorbar;
%     cx.Label.String='Shear Stress(MPa)';
%     caxis(ConstLim)   
% 
% %    Differences between abaqus stress and stress calculated with
%  %       constitutive model
%     subplot(2,2,3)
%     imagesc(Xcoord,Ycoord,DiffXY)
%     title(strcat('Const-Abaqus \sigma_{XY} ',timecount,'\mus'))
%     xlabel('X coordinate (mm)')
%     ylabel('Y coordinate (mm)')
%     colormap('hot')
%     cx=colorbar;
%     cx.Label.String='Stress diffence(MPa)';
%     caxis(ShearDlim)
%     
%  % %  Heat Map of Error confined to +-50%\
%    subplot(2,2,4)
%     imagesc(Xcoord,Ycoord,XYerr)
%      title(strcat('shear Error',timecount,'\mu s'))
%      xlabel('X coordinate (mm)')
%      ylabel('Y coordinate (mm)')
%      colormap('hot')
%      cx=colorbar;
%      cx.Label.String='error (%)';
%      caxis([-50,50])
%      
%  saveas(gcf,strcat(testdeg,'_XYdiff_',framecount,'.fig'))
%  saveas(gcf,strcat(testdeg,'_XYdiff_',framecount,'.png'))
% end
% XX stress
figure('units','normalized','outerposition',[0 0 1 1])
for n=1:length(time.vec)

 Xerr=squeeze(stressXXerr(:,:,n));
 Yerr=squeeze(stressYYerr(:,:,n));
 XYerr=squeeze(stressXYerr(:,:,n));
 timecount=num2str(timemic(n)); 
 framecount=num2str(n);
    
 ConstXX=squeeze(StressModel.xx(:,:,n))/10^6;   
 Xstress=squeeze(stress.x(:,:,n))/10^6;
 XYstress=squeeze(stress.s(:,:,n))/10^6;
 DiffX=squeeze(XXdiff(:,:,n))/10^6;
 DiffXY=squeeze(Sheardiff(:,:,n))/10^6; 
     
% %     Abaqus stress
    subplot(2,2,1)
    imagesc(Xcoord,Ycoord,Xstress)
    title(strcat('ABQ \sigma_{xx} ',timecount,'\mus'))
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='Stress(MPa)';
    caxis(XYstressLim)
    
%   % Constitutive Stress
    subplot(2,2,2)
    imagesc(Xcoord,Ycoord,ConstXX)
    title(strcat('Const \sigma_{xx} ',timecount,'\mus'))
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='Stress(MPa)';
    caxis(ConstLim)   

 % %   Differences between abaqus stress and stress calculated with
   %     constitutive model
    subplot(2,2,3)
    imagesc(Xcoord,Ycoord,DiffX)
    title(strcat('Const-Abaqus \sigma_{Xx} ',timecount,'\mus'))
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='Stress diffence(MPa)';
    caxis(XXDlim)
    
% %    Heat Map of Error confined to +-50%\
   subplot(2,2,4)
    imagesc(Xcoord,Ycoord,Xerr)
     title(strcat('\sigma_{xx} Error)',timecount,'\mu s'))
     xlabel('X coordinate (mm)')
     ylabel('Y coordinate (mm)')
     colormap('hot')
     cx=colorbar;
     cx.Label.String='error (%)';
     caxis([-50,50])
     
 saveas(gcf,strcat(testdeg,'_XXdiff_',framecount,'.fig'))
 saveas(gcf,strcat(testdeg,'_XXdiff_',framecount,'.png'))
end
    %% save data
    fprintf('saving data \n')
    save(strcat(testdeg,'_ConstError'),'StressModel','stress','strain',...
        'stressXXerr',...
        'stressXYerr',...
        'stressYYerr',...
        'XXdiff',...
        'YYdiff',...
        'Sheardiff',...
        'time',...
        'X_vec','Y_vec')
%% Save Constitutive model stress only
fprintf('Saving constitutive model stress \n')
save(strcat(testdeg,'_StressModel.mat'),'strain','time','StressModel');



    
 fprintf('done \n')  
    