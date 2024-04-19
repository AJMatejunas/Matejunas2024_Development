function answer = func_Compare3kinFields(data1,data2,data3,...
    Source1,Source2,Source3)
%This code is written to compare the kinematic fields from 3 separate
%sources. It generates comparison plots of Acceleration, Displacement, and
%Strain

%Author: Andrew Matejunas

%Function Input argurments
    %data1- Reference structure of the kinematic fields must contain the
        %following fields
            %accel- accelaration in x and y
            %disp- in-plane displacements in x, y/
            %strain- in-plane strains in xx, yy, and xy
            %X_vec-  vector of X coordinates
            %Y_vec- Vector of Y coordinates
            %time- time data        
     %data2- Other structure of kinematic fields. Must contain the same
        %fields as data1. Time data must be equal        
     %data3- Third structure to be compared
     %Source1 2 and 3- String describing the source of the corresponding
     %data structures
     
%Ourput argument
    %answer- arbitrary output argument. The primary output arguments are
    %image stacks comparing displacements, strains, and accelerations
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define Test Designation
testdeg=char(cell2mat(inputdlg('define test designation')));
%% Get time vector
timemic=data1.time.vec*10^6; %time vector in microseconds


%% Identify kinematic field limits
XDispLim=[min(data1.disp.x,[],'all'),max(data1.disp.x,[],'all')];
YDispLim=[min(data1.disp.y,[],'all'),max(data1.disp.y,[],'all')];

XAccelLim=[min(data1.accel.x,[],'all'),max(data1.accel.x,[],'all')];
YAccelLim=[min(data1.accel.y,[],'all'),max(data1.accel.y,[],'all')];

XXStrainLim=[min(data1.strain.x,[],'all'),max(data1.strain.x,[],'all')];
XYStrainLim=[min(data1.strain.s,[],'all'),max(data1.strain.s,[],'all')];

%% Select save directory
savePath=uigetdir('choose where to save plots');
dispPath=strcat(savePath,'/disp');
accelPath=strcat(savePath,'/accel');
strainPath=strcat(savePath,'/strain');

mkdir(dispPath);
mkdir(strcat(dispPath,'/fig'));
mkdir(strcat(dispPath,'/png'));
mkdir(accelPath);
mkdir(strcat(accelPath,'/fig'));
mkdir(strcat(accelPath,'/png'));
mkdir(strainPath);
mkdir(strcat(strainPath,'/fig'));
mkdir(strcat(strainPath,'/png'));
%% Generate Displacement Plots
figure('units','Centimeters','outerposition',[0 0 54.8 34.4])
for n=1:length(timemic)
   timecount=num2str(timemic(n)); 
   framecount=num2str(n); 
      
   Xdisp1=squeeze(data1.disp.x(:,:,n));
   Ydisp1=squeeze(data1.disp.y(:,:,n));
   Xdisp2=squeeze(data2.disp.x(:,:,n));
   Ydisp2=squeeze(data2.disp.y(:,:,n));
   Xdisp3=squeeze(data3.disp.x(:,:,n));
   Ydisp3=squeeze(data3.disp.y(:,:,n));
   
   
   
   
    %% Source 1 X displacements
    subplot(3,2,1)
    imagesc(data1.X_vec*10^3,data1.Y_vec*10^3,Xdisp1)
    title(strcat(Source1,' $$u_x$$ ',framecount,' $$\mu s$$'),'interpreter','latex')
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='$u_x$ (m)';
    cx.Label.Interpreter='latex';
    caxis(XDispLim)
    
  %% Source1 Y displacements
    subplot(3,2,2)
    imagesc(data1.X_vec*10^3,data1.Y_vec*10^3,Ydisp1)
    title(strcat(Source1,' $$u_y$$ ',framecount,' $$\mu s$$'),'interpreter','latex')
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='$u_y$ (m)';
    cx.Label.Interpreter='latex';
    caxis(YDispLim)
   
  %% Source 2 X displacements
    subplot(3,2,3)
    imagesc(data2.X_vec*10^3,data2.Y_vec*10^3,Xdisp2)
    title(strcat(Source2,' $$ u_x $$ ',framecount, '$$\mu s$$'),'interpreter','latex')
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='$u_x$ (m)';
    cx.Label.Interpreter='latex';
    caxis(XDispLim)
    
  %% Source 2 Y displacements
    subplot(3,2,4)
    imagesc(data2.X_vec*10^3,data2.Y_vec*10^3,Ydisp2)
    title(strcat(Source2,' $$u_y$$ ',framecount,' $$\mu s$$'),'interpreter','latex')
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='$u_y$ (m)';
    cx.Label.Interpreter='latex';
    caxis(YDispLim)
    
%% Source 3 X displacements
    subplot(3,2,5)
    imagesc(data3.X_vec*10^3,data3.Y_vec*10^3,Xdisp3)
    title(strcat(Source3,' $$u_x$$ ',framecount,' $$\mu s$$'),'interpreter','latex')
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='$u_x$ (m)';
    cx.Label.Interpreter='latex';
    caxis(XDispLim)
    
  %% Source 3 Y displacements
    subplot(3,2,6)
    imagesc(data2.X_vec*10^3,data2.Y_vec*10^3,Ydisp3)
    title(strcat(Source3,' $$u_y$$ ',framecount,' $$\mu s$$'),'interpreter','latex')
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='$u_y$ (m)';
    cx.Label.Interpreter='latex';
    caxis(YDispLim)
  
a=findall(gcf,'Type','axes');
h=findall(gcf,'Type','line');
t=findall(gcf,'Type','text');

set(t, ...
    'FontName','Helvetica', ...
    'FontSize',12, ...
    'FontWeight','bold' );

set(a                ,               ...
    'FontName'       , 'Helvetica' , ...
    'FontSize'       , 12     , ...
    'FontWeight'     , 'normal');

    
 saveas(gcf,strcat(dispPath,'/fig/',testdeg,'_disp_',framecount,'.fig'))
 saveas(gcf,strcat(dispPath,'/png/',testdeg,'_disp_',framecount,'.png'))
end

%% Plot Accelerations
figure('units','Centimeters','outerposition',[0 0 54.8 34.4])
for n=1:length(timemic)
   timecount=num2str(timemic(n)); 
   framecount=num2str(n); 
      
   Xaccel1=squeeze(data1.accel.x(:,:,n));
   Yaccel1=squeeze(data1.accel.y(:,:,n));
   Xaccel2=squeeze(data2.accel.x(:,:,n));
   Yaccel2=squeeze(data2.accel.y(:,:,n));
   Xaccel3=squeeze(data3.accel.x(:,:,n));
   Yaccel3=squeeze(data3.accel.y(:,:,n));
   
   
   
   
    %% Source 1 X accelerations
    subplot(3,2,1)
    imagesc(data1.X_vec*10^3,data1.Y_vec*10^3,Xaccel1)
    title(strcat(Source1,' $$a_x$$ ',framecount,' $$\mu s$$'),'interpreter','latex')
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='$a_x$ (m/s^2)';
    cx.Label.Interpreter='latex';
    caxis(XAccelLim)
    
  %% Source1 Y accelerations
    subplot(3,2,2)
    imagesc(data1.X_vec*10^3,data1.Y_vec*10^3,Yaccel1)
    title(strcat(Source1,' $$a_y$$ ',framecount,' $$\mu s$$'),'interpreter','latex')
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='$a_y$ (m/s^2)';
    cx.Label.Interpreter='latex';
    caxis(YAccelLim)
   
  %% Source 2 X accelerations
    subplot(3,2,3)
    imagesc(data2.X_vec*10^3,data2.Y_vec*10^3,Xaccel2)
    title(strcat(Source2,' $$a_x$$ ',framecount,' $$\mu s$$ '),'interpreter','latex')
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='$a_x$ (m/s^2)';
    cx.Label.Interpreter='latex';
    caxis(XAccelLim)
    
  %% Source 2 Y accelerations
    subplot(3,2,4)
    imagesc(data2.X_vec*10^3,data2.Y_vec*10^3,Yaccel2)
    title(strcat(Source2,' $$a_y$$ ',framecount,' $$\mu s$$ '),'interpreter','latex')
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='$a_y$ (m/s^2)';
    cx.Label.Interpreter='latex';
    caxis(YAccelLim)
    
%% Source 3 X accelerations
    subplot(3,2,5)
    imagesc(data3.X_vec*10^3,data3.Y_vec*10^3,Xaccel3)
    title(strcat(Source3,' $$a_x$$ ',framecount,' $$\mu s$$'),'interpreter','latex')
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='$a_x$ (m/s^2)';
    cx.Label.Interpreter='latex';
    caxis(XAccelLim)
    
  %% Source 3 Y accelerations
    subplot(3,2,6)
    imagesc(data2.X_vec*10^3,data2.Y_vec*10^3,Yaccel3)
    title(strcat(Source3,' $$a_y$$ ',framecount,' $$\mu s$$'),'interpreter','latex')
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='$a_y$ (m/s^2)';
    cx.Label.Interpreter='latex';
    caxis(YAccelLim)
  
a=findall(gcf,'Type','axes');
h=findall(gcf,'Type','line');
t=findall(gcf,'Type','text');

set(t, ...
    'FontName','Helvetica', ...
    'FontSize',12, ...
    'FontWeight','bold' );

set(a                ,               ...
    'FontName'       , 'Helvetica' , ...
    'FontSize'       , 12     , ...
    'FontWeight'     , 'normal');

    
 saveas(gcf,strcat(accelPath,'/fig/',testdeg,'_accel_',framecount,'.fig'));
 saveas(gcf,strcat(accelPath,'/png/',testdeg,'_accel_',framecount,'.png'));
end

%% Plot Strains
figure('units','Centimeters','outerposition',[0 0 54.8 34.4])
for n=1:length(timemic)
   timecount=num2str(timemic(n)); 
   framecount=num2str(n); 
      
   Xstrain1=squeeze(data1.strain.x(:,:,n));
   XYstrain1=squeeze(data1.strain.s(:,:,n));
   Xstrain2=squeeze(data2.strain.x(:,:,n));
   XYstrain2=squeeze(data2.strain.s(:,:,n));
   Xstrain3=squeeze(data3.strain.x(:,:,n));
   XYstrain3=squeeze(data3.strain.s(:,:,n));
   
   
   
   
    %% Source 1 X strains
    subplot(3,2,1)
    imagesc(data1.X_vec*10^3,data1.Y_vec*10^3,Xstrain1)
    title(strcat(Source1,' $$\varepsilon{}_{xx}$$ ',framecount,' $$\mu{} s$$'),'interpreter','latex')
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='$\varespsilon{}_{xx}$';
    cx.Label.Interpreter='latex';
    caxis(XXStrainLim)
    
  %% Source1 XY strains
    subplot(3,2,2)
    imagesc(data1.X_vec*10^3,data1.Y_vec*10^3,XYstrain1)
    title(strcat(Source1,' $$\varepsilon{}_{xy}$$ ',framecount,' $$\mu{}s$$'),'interpreter','latex')
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='$\varepsilon{}_{xy}$';
    cx.Label.Interpreter='latex';
    caxis(XYStrainLim)
   
  %% Source 2 X strains
    subplot(3,2,3)
    imagesc(data2.X_vec*10^3,data2.Y_vec*10^3,Xstrain2)
    title(strcat(Source2,' $$\varepsilon{}_{xx}$$ ',framecount,' $$\mu{}s$$'),'interpreter','latex')
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='$\varepsilon{}_{xx}$';
    cx.Label.Interpreter='latex';
    caxis(XXStrainLim)
    
  %% Source 2 XY strains
    subplot(3,2,4)
    imagesc(data2.X_vec*10^3,data2.Y_vec*10^3,XYstrain2)
    title(strcat(Source2,' $$\varepsilon{}_{xy}$$ ',framecount,' $$\mu{}s$$'),'interpreter','latex')
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='$\varepsilon{}_{xy}$';
    cx.Label.Interpreter='latex';
    caxis(XYStrainLim)
    
%% Source 3 X strains
    subplot(3,2,5)
    imagesc(data3.X_vec*10^3,data3.Y_vec*10^3,Xstrain3)
    title(strcat(Source3,' $$\varepsilon{}_{xx}$$ ',framecount,' $$\mu{}s$$'),'interpreter','latex')
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='$\varepsilon{}_{xx}$)';
    cx.Label.Interpreter='latex';
    caxis(XXStrainLim)
    
  %% Source 3 XY strains
    subplot(3,2,6)
    imagesc(data2.X_vec*10^3,data2.Y_vec*10^3,XYstrain3)
    title(strcat(Source3,' $$\varepsilon{}_{xy}$$ ',framecount,' $$\mu{}s$$'),'interpreter','latex')
    xlabel('X coordinate (mm)')
    ylabel('Y coordinate (mm)')
    colormap('hot')
    cx=colorbar;
    cx.Label.String='$\varepsilon{}_{xy}$';
    cx.Label.Interpreter='latex';
    caxis(XYStrainLim)
  
a=findall(gcf,'Type','axes');
h=findall(gcf,'Type','line');
t=findall(gcf,'Type','text');

set(t, ...
    'FontName','Helvetica', ...
    'FontSize',12, ...
    'FontWeight','bold' );

set(a                ,               ...
    'FontName'       , 'Helvetica' , ...
    'FontSize'       , 12     , ...
    'FontWeight'     , 'normal');

    
 saveas(gcf,strcat(strainPath,'/fig/',testdeg,'_strain_',framecount,'.fig'));
 saveas(gcf,strcat(strainPath,'/png/',testdeg,'_strain_',framecount,'.png'));
end

end

