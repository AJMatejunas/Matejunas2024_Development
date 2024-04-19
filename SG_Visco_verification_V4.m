% This code is written to verify viscoelasticity codes using the stress
% gage equation

%Author: Andrew Matejunas

%Version history:
    %V2: Added capability to handle synthetically deformed images without
        %stress information
    %V3: Added shear stress guage functionality
    %v4: Created a single function to calculate all stress gage stresses.
        %Removed option to calculate only normal stress-gage stresses
        %2022-05-31
    
%Change log
    %2022/02/28- Included ProgramVersions data structure to track stress
        %gauge algorithms used (Did not include a version change)


%% initialize
clear all; close all; clc

%% Determine whether analysis is to be performed on pure FE data
    %if the answer is Yes both stress and deformation data will come from
        %the finite element method outputs
    %if the answer is 'No' the stress calculated with the Stress Gage
        %method will be compared to stress calculated with the finite
        %element method
        
quest='Is this pure FE data?';
PureFE=questdlg(quest);

%% load data
prompt='Select file containing deformation information';
[dispfilename,dispPath]=uigetfile('.mat',prompt);
dispfile=strcat(dispPath,'/',dispfilename);
prompt='input test designation';
desig=char(inputdlg(prompt));

fprintf('loading data \n')


switch PureFE
    case 'No'
        prompt='select file containing stress information';
        stressfile=uigetfile('.mat',prompt);
        FE=load(stressfile);
       
       %% Extract spatial coorditnates
        FE.X_vec=FE.Xq(1,:);
        FE.Y_vec=FE.Yq(:,1);
        
       %% calculate average stresses in xx direction
       FE.avgX_stress=mean(FE.stress.x);
       FE.avgX_strain=mean(FE.strain.x);
       FE.strain_surf=squeeze(FE.avgX_strain);
       FE.avg_accel=mean(FE.accel.x);

       %% calculate average stresses and strains in yy direction
       FE.avgY_stress=mean(FE.stress.y);
       FE.avgY_strain=mean(FE.strain.y);

       %% XY direction
       FE.avgXY_stress=squeeze(mean(FE.stress.s));
       FE.avgXY_strain=squeeze(mean(FE.strain.s));
       
       FE.frames=length(FE.time.vec);

       %% FE stress gage
       %allocate memory
        FE.accel_surf=zeros(length(FE.X_vec),FE.frames);
        FE.SG=zeros(length(FE.X_vec),FE.frames);
        FE.stressX_ref=zeros(length(FE.X_vec),FE.frames);
        FE.stressY_ref=zeros(length(FE.X_vec),FE.frames);

        for m=1:length(FE.X_vec)
            for n=1:length(FE.time.vec)
        %surface average acceleration
        FE.accel_surf(m,n)=mean(FE.avg_accel(1,1:m,n));

         %stress guage stress
%         FE.SG(m,n)=FE.material.rho*FE.X_vec(m)*FE.accel_surf(m,n);
        %reference stress from FE data
        FE.stressX_ref(m,n)=FE.avgX_stress(1,m,n);
        FE.stressY_ref(m,n)=FE.avgY_stress(1,m,n);       
            end
        end

        %reshape stress and strain into a vector
        %strain_vec=reshape(strain_surf,[],1);
        FE.strain_vec=reshape(FE.avgX_strain,[],1);

        FE.SG_vec=reshape(FE.SG,[],1);
        FE.stress_vec=reshape(FE.stressX_ref,[],1);

    FE.Full_SG=func_Full_SG(FE.accel,FE.X_vec,FE.time,FE.material.rho);
    FE.SG=FE.Full_SG.x;
    FE.Shear_SG=FE.Full_SG.s;
    
    FE.Shear_Avgx=mean(FE.stress.s);
    
 
    
end

load(dispfile);
%% Define test designation to autosave figures

testdeg=desig;

fprintf('Calculating stress gage \n')
%% Create vectors of spatial coordinates 
switch PureFE
    case 'Yes'
    %in pure FE data X coordinates are stored in Xq    
    X_vec=Xq(1,:);
    Y_vec=Yq(:,1);
    case'No'
    %In processed grid method images x coordinates are stored in pos.x
    X_vec=pos.x;
    Y_vec=pos.y;
end

%% Determine time step and number of frames
dt=time.vec(2)-time.vec(1);
NumFrames=length(time.vec);

%% calculate average strains accelerations in xx direction

avg_accel=mean(accel.x);
avgX_strain=mean(strain.x);
strain_surf=squeeze(avgX_strain);

switch PureFE
    case 'Yes'
    %% Determine average strain and stresses
    %stress in xx direction
    avgX_stress=mean(stress.x);
    %yy direction
    avgY_stress=mean(stress.y);
    avgY_strain=mean(strain.y);

    %xy direction
    avgXY_stress=squeeze(mean(stress.s));
    
    
end
avgXY_strain=squeeze(mean(strain.s));

%% Calculate stress gage stress
%allocate memory
ProgramVersions.NormalSG_alg='SG_visco_verification_V4';
accel_surf=zeros(length(X_vec),NumFrames);


switch PureFE
    case 'Yes'
stressX_ref=zeros(length(X_vec),NumFrames);
stressY_ref=zeros(length(X_vec),NumFrames);
end



for m=1:length(X_vec)
    for n=1:NumFrames
        %surface average acceleration
        accel_surf(m,n)=mean(avg_accel(1,1:m,n));

         %stress guage stress
        %SG(m,n)=material.rho*X_vec(m)*accel_surf(m,n);
        
        switch PureFE
            case 'Yes'
        %reference stress from
        stressX_ref(m,n)=avgX_stress(1,m,n);
        stressY_ref(m,n)=avgY_stress(1,m,n);       
        end
    end
end
Full_SG=func_Full_SG(accel,X_vec,time,material.rho);
    SG=Full_SG.x;
    Shear_SG=Full_SG.s;
    
   
    Shear_SGvec=reshape(Shear_SG,[],1);
    switch PureFE
        case 'Yes'
    Shear_Avgxvec=reshape(avgXY_stress,[],1);
    end
    Shear_strainAV_vec=reshape(avgXY_strain,[],1);
    

%reshape stress and strain into a vector
%strain_vec=reshape(strain_surf,[],1);
strain_vec=reshape(avgX_strain,[],1);

SG_vec=reshape(SG,[],1);
    %%
    switch PureFE
        case 'Yes'
        stress_vec=reshape(stressX_ref,[],1);

        %% Calculate Plane stress
        %sigma_xx-nu*sigma(yy)
        PL_stress=stressX_ref-material.nuxy*stressY_ref;
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
switch PureFE
    case 'Yes'

%% Create comparison plot

figure(1)

%sample 6 locations in the 
    ind1=1;
    ind2=floor(length(X_vec)/5);
    ind3=2*ind2;
    ind4=3*ind2;
    ind5=4*ind2;
    ind6=5*ind2;

    lgd1=num2str(X_vec(1));
    lgd2=num2str(X_vec(ind2));
    lgd3=num2str(X_vec(ind3));
    lgd4=num2str(X_vec(ind4));
    lgd5=num2str(X_vec(ind5));
    lgd6=num2str(X_vec(ind6));

    % plot SG results
    plot(time.vec,SG(ind1,:),'b')
    hold on
    plot(time.vec,SG(ind2,:),'k')
    plot(time.vec,SG(ind3,:),'r')
    plot(time.vec,SG(ind4,:),'g')
    plot(time.vec,SG(ind5,:),'y')
    plot(time.vec,SG(ind6,:),'c')

    %plot ABAQUS stress results
    plot(time.vec,stressX_ref(ind1,:),'--b')
    plot(time.vec,stressX_ref(ind2,:),'--k','handlevisibility','off')
    plot(time.vec,stressX_ref(ind3,:),'--r','handlevisibility','off')
    plot(time.vec,stressX_ref(ind4,:),'--g','handlevisibility','off')
    plot(time.vec,stressX_ref(ind5,:),'--y','handlevisibility','off')
    plot(time.vec,stressX_ref(ind6,:),'--c','handlevisibility','off')
    hold off
    
    xlabel('time (s)')
    ylabel('stress (Pa)')

    legend(strcat(lgd1,' sg'),lgd2,lgd3,lgd4,lgd5,lgd6,strcat(lgd1,' ABAQUS'),'location','southeast')

    saveas(gcf,strcat(testdeg,'_StressTime.fig'))
    saveas(gcf,strcat(testdeg,'_StressTime.png'))
    %% Create stress strain plot


    figure(2)
    scatter(strain_vec,SG_vec)
    hold on
    scatter(strain_vec,stress_vec)
    hold off
    xlabel('strain')
    ylabel('stress (Pa)')
    legend('SG','ABAQUS')

    saveas(gcf,strcat(testdeg,'_StressStrain.fig'))
    saveas(gcf,strcat(testdeg,'_StressStrain.png'))

    %% Plot stress strain at different x locations
    figure(3)
    % plot SG results
    plot(strain_surf(ind1,:),SG(ind1,:),'b')
    hold on
    plot(strain_surf(ind2,:),SG(ind2,:),'k')
    plot(strain_surf(ind3,:),SG(ind3,:),'r')
    plot(strain_surf(ind4,:),SG(ind4,:),'g')
    plot(strain_surf(ind5,:),SG(ind5,:),'y')
    plot(strain_surf(ind6,:),SG(ind6,:),'c')

    %plot ABAQUS stress results
    plot(strain_surf(ind1,:),stressX_ref(ind1,:),'--b')
    plot(strain_surf(ind2,:),stressX_ref(ind2,:),'--k','handlevisibility','off')
    plot(strain_surf(ind3,:),stressX_ref(ind3,:),'--r','handlevisibility','off')
    plot(strain_surf(ind4,:),stressX_ref(ind4,:),'--g','handlevisibility','off')
    plot(strain_surf(ind5,:),stressX_ref(ind5,:),'--y','handlevisibility','off')
    plot(strain_surf(ind6,:),stressX_ref(ind6,:),'--c','handlevisibility','off')
    hold off

    xlabel('strain')
    ylabel('\sigma_{xx} (Pa)')

    legend(strcat(lgd1,' sg'),lgd2,lgd3,lgd4,lgd5,lgd6,strcat(lgd1,' ABAQUS'),'location','southeast')

    saveas(gcf,strcat(testdeg,'_SG.fig'))
    saveas(gcf,strcat(testdeg,'_SG.png'))
    %% plot y stress
    figure(4)
    plot(time.vec,stressY_ref(ind1,:),'--b')
    hold on
    plot(time.vec,stressY_ref(ind2,:),'--k')
    plot(time.vec,stressY_ref(ind3,:),'--r')
    plot(time.vec,stressY_ref(ind4,:),'--g')
    plot(time.vec,stressY_ref(ind5,:),'--y')
    plot(time.vec,stressY_ref(ind6,:),'--c')
    hold off

    ylabel('\sigma_{yy} (Pa)')
    xlabel('time (s)')
    legend(strcat(lgd1,'ABAQUS'),lgd2,lgd3,lgd4,lgd5,lgd6)

    saveas(gcf,strcat(testdeg,'_YStressTime.fig'))
    saveas(gcf,strcat(testdeg,'_YStressTime.png'))

%% Plot plane stress stiffness
%note IBII tests are PLANE STRESS conditions NOT uniaxial!!!!!!

    figure(5)

    plot(strain_surf(ind1,:),PL_stress(ind1,:),'--b')
    hold on
    plot(strain_surf(ind2,:),PL_stress(ind2,:))
    plot(strain_surf(ind3,:),PL_stress(ind3,:))
    plot(strain_surf(ind4,:),PL_stress(ind4,:))
    plot(strain_surf(ind5,:),PL_stress(ind5,:))
    plot(strain_surf(ind6,:),PL_stress(ind6,:))
    hold off

    ylabel('\sigma_{xx}-\nu*\sigma_{yy} (Pa)')
    xlabel('$\varepsilon_{xx} $','interpreter','latex')
    legend(strcat(lgd1,'ABAQUS'),lgd2,lgd3,lgd4,lgd5,lgd6)

    saveas(gcf,strcat(testdeg,'PL_Stress.fig'))
    saveas(gcf,strcat(testdeg,'PL_Stress.png'))

%% Shear Stress Gauge Plots

%% Plot All Shear Stresses
    figure(6)
    scatter(Shear_strainAV_vec,Shear_SGvec)
    hold on
    scatter(Shear_strainAV_vec,Shear_Avgxvec)
    hold off
    legend('SG_{shear}','ABAQUS')
    xlabel('Average strain_{xy}')
    ylabel('Average \sigma_{xy} (Pa)')
    title('All shear stresses')
    
    saveas(gcf,strcat(testdeg,'_XYsgAll.fig'))
    saveas(gcf,strcat(testdeg,'_XYsgAll.png'))
    
    figure(7)
    % plot SG results
    plot(time.vec,Shear_SG(ind1,:),'b')
    hold on
    plot(time.vec,Shear_SG(ind2,:),'k')
    plot(time.vec,Shear_SG(ind3,:),'r')
    plot(time.vec,Shear_SG(ind4,:),'g')
    plot(time.vec,Shear_SG(ind5,:),'y')
    plot(time.vec,Shear_SG(ind6,:),'c')

    %plot ABAQUS stress results
    plot(time.vec,avgXY_stress(ind1,:),'--b')
    plot(time.vec,avgXY_stress(ind2,:),'--k','handlevisibility','off')
    plot(time.vec,avgXY_stress(ind3,:),'--r','handlevisibility','off')
    plot(time.vec,avgXY_stress(ind4,:),'--g','handlevisibility','off')
    plot(time.vec,avgXY_stress(ind5,:),'--y','handlevisibility','off')
    plot(time.vec,avgXY_stress(ind6,:),'--c','handlevisibility','off')
    hold off
    
    xlabel('time (s)')
    ylabel('stress (Pa)')

    legend(strcat(lgd1,' sg'),lgd2,lgd3,lgd4,lgd5,lgd6,strcat(lgd1,' ABAQUS'),'location','southeast')

    saveas(gcf,strcat(testdeg,'_shearStressTime.fig'))
    saveas(gcf,strcat(testdeg,'_shearStressTime.png'))
    
    figure(8)
    % plot SG results
    plot(avgXY_strain(ind1,:),Shear_SG(ind1,:),'b')
    hold on
    plot(avgXY_strain(ind2,:),Shear_SG(ind2,:),'k')
    plot(avgXY_strain(ind3,:),Shear_SG(ind3,:),'r')
    plot(avgXY_strain(ind4,:),Shear_SG(ind4,:),'g')
    plot(avgXY_strain(ind5,:),Shear_SG(ind5,:),'y')
    plot(avgXY_strain(ind6,:),Shear_SG(ind6,:),'c')

    %plot ABAQUS stress results
    plot(avgXY_strain(ind1,:),avgXY_stress(ind1,:),'--b')
    plot(avgXY_strain(ind2,:),avgXY_stress(ind2,:),'--k','handlevisibility','off')
    plot(avgXY_strain(ind3,:),avgXY_stress(ind3,:),'--r','handlevisibility','off')
    plot(avgXY_strain(ind4,:),avgXY_stress(ind4,:),'--g','handlevisibility','off')
    plot(avgXY_strain(ind5,:),avgXY_stress(ind5,:),'--y','handlevisibility','off')
    plot(avgXY_strain(ind6,:),avgXY_stress(ind6,:),'--c','handlevisibility','off')
    hold off

    xlabel('strain')
    ylabel('\sigma_{xy} (Pa)')

    legend(strcat(lgd1,' sg'),lgd2,lgd3,lgd4,lgd5,lgd6,strcat(lgd1,' ABAQUS'),'location','southeast')

    saveas(gcf,strcat(testdeg,'_ShearSG.fig'))
    saveas(gcf,strcat(testdeg,'_ShearSG.png'))
    

    case 'No'
    %% Set up indicies for plot data
    %Data from images
    ind1=1;
    ind2=round(length(X_vec)/5);
    ind3=2*ind2;
    ind4=3*ind2;
    ind5=4*ind2;
    ind6=length(X_vec);
    
    
    lgd1=num2str(X_vec(1));
    lgd2=num2str(X_vec(ind2));
    lgd3=num2str(X_vec(ind3));
    lgd4=num2str(X_vec(ind4));
    lgd5=num2str(X_vec(ind5));
    lgd6=num2str(X_vec(end));
    
    %data from FE file
    FE.ind1=1;
    FE.ind2=round(length(FE.X_vec)/5);
    FE.ind3=2*FE.ind2;
    FE.ind4=3*FE.ind2;
    FE.ind5=4*FE.ind2;
    FE.ind6=length(FE.X_vec);
    
    FE.lgd1=num2str(X_vec(1));
    FE.lgd2=num2str(X_vec(ind2));
    FE.lgd3=num2str(X_vec(ind3));
    FE.lgd4=num2str(X_vec(ind4));
    FE.lgd5=num2str(X_vec(ind5));
    FE.lgd6=num2str(X_vec(end));
    
    %% create scatter plot of all stress and strain
    figure(2)
    scatter(strain_vec,SG_vec)
    hold on
    scatter(FE.strain_vec,FE.stress_vec)
    hold off
    xlabel('strain')
    ylabel('stress_{xx} (Pa)')
    legend('SG (images)','ABAQUS')

    saveas(gcf,strcat(testdeg,'_StressStrain.fig'))
    saveas(gcf,strcat(testdeg,'_StressStrain.png'))
    
     %% Plot stress strain at different x locations
    figure(3)
    % plot SG results from images
    plot(strain_surf(ind1,:),SG(ind1,:),'b')
    hold on
    plot(strain_surf(ind2,:),SG(ind2,:),'k')
    plot(strain_surf(ind3,:),SG(ind3,:),'r')
    plot(strain_surf(ind4,:),SG(ind4,:),'g')
    plot(strain_surf(ind5,:),SG(ind5,:),'y')
    plot(strain_surf(ind6,:),SG(ind6,:),'c')

    %plot ABAQUS stress results
    plot(FE.strain_surf(FE.ind1,:),FE.stressX_ref(FE.ind1,:),'--b')
    plot(FE.strain_surf(FE.ind2,:),FE.stressX_ref(FE.ind2,:),'--k','handlevisibility','off')
    plot(FE.strain_surf(FE.ind3,:),FE.stressX_ref(FE.ind3,:),'--r','handlevisibility','off')
    plot(FE.strain_surf(FE.ind4,:),FE.stressX_ref(FE.ind4,:),'--g','handlevisibility','off')
    plot(FE.strain_surf(FE.ind5,:),FE.stressX_ref(FE.ind5,:),'--y','handlevisibility','off')
    plot(FE.strain_surf(FE.ind6,:),FE.stressX_ref(FE.ind6,:),'--c','handlevisibility','off')
    hold off

    xlabel('strain')
    ylabel('\sigma_{xx} (Pa)')
 
    legend(strcat(lgd1,' sg (images)'),lgd2,lgd3,lgd4,lgd5,lgd6,strcat(lgd1,' FE'),'location','southeast')

    saveas(gcf,strcat(testdeg,'_SG.fig'))
    saveas(gcf,strcat(testdeg,'_SG.png'))
    
    
    %% Plot shear stresses
    
     % plot SG results from images
    plot(avgXY_strain(ind1,:),Shear_SG(ind1,:),'b')
    hold on
    plot(avgXY_strain(ind2,:),Shear_SG(ind2,:),'k')
    plot(avgXY_strain(ind3,:),Shear_SG(ind3,:),'r')
    plot(avgXY_strain(ind4,:),Shear_SG(ind4,:),'g')
    plot(avgXY_strain(ind5,:),Shear_SG(ind5,:),'y')
    plot(avgXY_strain(ind6,:),Shear_SG(ind6,:),'c')

    %plot ABAQUS stress results
    plot(FE.avgXY_strain(FE.ind1,:),FE.avgXY_stress(FE.ind1,:),'--b')
    plot(FE.avgXY_strain(FE.ind2,:),FE.avgXY_stress(FE.ind2,:),'--k','handlevisibility','off')
    plot(FE.avgXY_strain(FE.ind3,:),FE.avgXY_stress(FE.ind3,:),'--r','handlevisibility','off')
    plot(FE.avgXY_strain(FE.ind4,:),FE.avgXY_stress(FE.ind4,:),'--g','handlevisibility','off')
    plot(FE.avgXY_strain(FE.ind5,:),FE.avgXY_stress(FE.ind5,:),'--y','handlevisibility','off')
    plot(FE.avgXY_strain(FE.ind6,:),FE.avgXY_stress(FE.ind6,:),'--c','handlevisibility','off')
    
    xlabel('strain_{xy}')
    ylabel('\sigma_{xy} (Pa)')

    legend(strcat(lgd1,' sg'),lgd2,lgd3,lgd4,lgd5,lgd6,strcat(lgd1,' ABAQUS'),'location','southeast')

    saveas(gcf,strcat(testdeg,'_ShearSG.fig'))
    saveas(gcf,strcat(testdeg,'_ShearSG.png'))
   
end

%% Save data
fprintf('saving data \n')
save(strcat(testdeg,'_SGdata.mat'))

fprintf('Done \n')