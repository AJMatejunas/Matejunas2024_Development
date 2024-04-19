function [X_vec,Y_vec,Full_SG]=...
    func_FEViscoSGCalc(Desig,SaveDir,Xq,Yq,accel,strain,stress,time,...
    material)
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
    %2022-08-25: Converted to a function for only pure FE data to allow for
        %faster calculations


%% Create vectors of spatial coordinates

X_vec=Xq(1,:);
Y_vec=Yq(:,1);

%% Determine time step and number of frames
dt=time.vec(2)-time.vec(1);
NumFrames=length(time.vec);

%% calculate average strains accelerations in xx direction

avg_accel=mean(accel.x);
avgX_strain=mean(strain.x);
strain_surf=squeeze(avgX_strain);

%% Determine average strain and stresses
%stress in xx direction
avgX_stress=mean(stress.x);
%yy direction
avgY_stress=mean(stress.y);
avgY_strain=mean(strain.y);

%xy direction
avgXY_stress=squeeze(mean(stress.s));



avgXY_strain=squeeze(mean(strain.s));

%% Calculate stress gage stress
%allocate memory
ProgramVersions.NormalSG_alg='SG_visco_verification_V4';
accel_surf=zeros(length(X_vec),NumFrames);

stressX_ref=zeros(length(X_vec),NumFrames);
stressY_ref=zeros(length(X_vec),NumFrames);




for m=1:length(X_vec)
    for n=1:NumFrames
        %surface average acceleration
        accel_surf(m,n)=mean(avg_accel(1,1:m,n));

        %stress guage stress
        %SG(m,n)=material.rho*X_vec(m)*accel_surf(m,n);

        %reference stress from
        stressX_ref(m,n)=avgX_stress(1,m,n);
        stressY_ref(m,n)=avgY_stress(1,m,n);
    end
end

Full_SG=func_Full_SG(accel,X_vec,time,material.rho);
    SG=Full_SG.x;
    Shear_SG=Full_SG.s;
    
   
    Shear_SGvec=reshape(Shear_SG,[],1);
    Shear_Avgxvec=reshape(avgXY_stress,[],1);
    Shear_strainAV_vec=reshape(avgXY_strain,[],1);
    

%reshape stress and strain into a vector
%strain_vec=reshape(strain_surf,[],1);
strain_vec=reshape(avgX_strain,[],1);

SG_vec=reshape(SG,[],1);
    %%

        stress_vec=reshape(stressX_ref,[],1);

        %% Calculate Plane stress
        %sigma_xx-nu*sigma(yy)
        PL_stress=stressX_ref-material.nuxy*stressY_ref;

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
SGFigDir=strcat(SaveDir,'\SG');
mkdir(SGFigDir)

%% Create comparison plot

figure(1)

%sample 6 locations in the 
    ind1=1;
    ind2=length(X_vec)/5;
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

    legend(strcat(lgd1,' sg'),lgd2,lgd3,lgd4,lgd5,...
        lgd6,strcat(lgd1,' ABAQUS'),'location','southeast')

    saveas(gcf,strcat(SGFigDir,'\',Desig,'_StressTime.fig'))
    saveas(gcf,strcat(SGFigDir,'\',Desig,'_StressTime.png'))
    %% Create stress strain plot


    figure(2)
    scatter(strain_vec,SG_vec)
    hold on
    scatter(strain_vec,stress_vec)
    hold off
    xlabel('strain')
    ylabel('stress (Pa)')
    legend('SG','ABAQUS')

    saveas(gcf,strcat(SGFigDir,'\',Desig,'_StressStrain.fig'))
    saveas(gcf,strcat(SGFigDir,'\',Desig,'_StressStrain.png'))

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

    saveas(gcf,strcat(SGFigDir,'\',Desig,'_SG.fig'))
    saveas(gcf,strcat(SGFigDir,'\',Desig,'_SG.png'))

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

    saveas(gcf,strcat(SGFigDir,'\',Desig,'_YStressTime.fig'))
    saveas(gcf,strcat(SGFigDir,'\',Desig,'_YStressTime.png'))

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

    saveas(gcf,strcat(SGFigDir,'\',Desig,'PL_Stress.fig'))
    saveas(gcf,strcat(SGFigDir,'\',Desig,'PL_Stress.png'))

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
    
    saveas(gcf,strcat(SGFigDir,'\',Desig,'_XYsgAll.fig'))
    saveas(gcf,strcat(SGFigDir,'\',Desig,'_XYsgAll.png'))
    
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

    saveas(gcf,strcat(SGFigDir,'\',Desig,'_shearStressTime.fig'))
    saveas(gcf,strcat(SGFigDir,'\',Desig,'_shearStressTime.png'))
    
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

    saveas(gcf,strcat(SGFigDir,'\',Desig,'_ShearSG.fig'))
    saveas(gcf,strcat(SGFigDir,'\',Desig,'_ShearSG.png'))
    


%% Save data
fprintf('saving SG data \n')
save(strcat(SaveDir,Desig,'_SGdata.mat'))
close all
fprintf('SG Saved \n')
end