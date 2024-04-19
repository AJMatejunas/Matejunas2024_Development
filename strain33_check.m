%% This script is written to verify that func_ViscoIBIIstrain33 is
%functional

%% initialize
clear all, close all, clc

%% load necessary data
[SGname,SGpath]=uigetfile('.mat','Choose file containing SG data');
SGfile=strcat(SGpath,'/',SGname);

%% load properties
[Propname,Proppath]=uigetfile('.mat',...
    'Choose file containing material properties');
Propfile=strcat(Proppath,'/',Propname);

fprintf('Loading Properties \n')
load(SGfile,'strain','time','testdeg','X_vec');
load(Propfile);

%% Extract necessary properties
G0=MatProps.G0;
K0=MatProps.K0;
Gi=MatProps.Gi;
Ki=MatProps.Ki;

tau=MatProps.tau;
%% Calculate C_visco0
C_Visco0=func_ViscoElasKG(K0,G0);

%% Calculate C_Viscom
C_Viscom=func_ViscoElasKG(Gi,Ki);

%% extract strains
strain11=strain.x;
strain22=strain.y;

%% Calculate strain33
fprintf('Calculating Out of plane strain \n')
strain33=func_ViscoIBIIstrain33_debug(strain11,strain22,... input strains
    time.vec,C_Visco0,C_Viscom,tau);
%% Caclulate averages
strain33Avg=squeeze(mean(strain33));
strain11Avg=squeeze(mean(strain.x));
%% plot results, do they make sense or are they dumb?
fprintf('Plotting results \n')
Max_strain=max(strain11Avg,[],'all');
Min_strain=min(strain11Avg,[],'all');
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


figure(1)
plot(time.vec,strain33Avg(ind1,:),'--b')
hold on
plot(time.vec,strain33Avg(ind2,:),'--k')
plot(time.vec,strain33Avg(ind3,:),'--r')
plot(time.vec,strain33Avg(ind4,:),'--g')
plot(time.vec,strain33Avg(ind5,:),'--y')
plot(time.vec,strain33Avg(ind6,:),'--c')

plot(time.vec,strain11Avg(ind1,:),'-b')
plot(time.vec,strain11Avg(ind2,:),'-k','handlevisibility','off')
plot(time.vec,strain11Avg(ind3,:),'-r','handlevisibility','off')
plot(time.vec,strain11Avg(ind4,:),'-g','handlevisibility','off')
plot(time.vec,strain11Avg(ind5,:),'-y','handlevisibility','off')
plot(time.vec,strain11Avg(ind6,:),'-c','handlevisibility','off')



hold off

xlabel('time (s)')
ylabel('strain')
ylim([Min_strain,Max_strain])

legend(strcat(lgd1,' strain_{33}'),lgd2,lgd3,lgd4,lgd5,lgd6,strcat(lgd1,'strain_{11}'),'location','southeast')

saveas(gcf,strcat(testdeg,'_strain33_11time.fig'))
saveas(gcf,strcat(testdeg,'_strain33_11time.png'))
fprintf('Done \n')

%% plot only strain33
figure(2)
plot(time.vec,strain33Avg(ind1,:),'--b')
hold on
plot(time.vec,strain33Avg(ind2,:),'--k')
plot(time.vec,strain33Avg(ind3,:),'--r')
plot(time.vec,strain33Avg(ind4,:),'--g')
plot(time.vec,strain33Avg(ind5,:),'--y')
plot(time.vec,strain33Avg(ind6,:),'--c')
hold off


xlabel('time (s)')
ylabel('strain_{33}')

legend(strcat(lgd1,' strain_{33}'),lgd2,lgd3,lgd4,lgd5,lgd6)

saveas(gcf,strcat(testdeg,'_strain33time.fig'))
saveas(gcf,strcat(testdeg,'_strain33time.png'))

%% strain33 Center element
figure(3)
X_mid=length(strain11(1,:,1))/2;
Y_mid=length(strain11(:,1,1))/2;

midstrain33=squeeze(strain33(Y_mid,X_mid,:));
plot(time.vec,midstrain33)
legend('strain_{33} at mid point')

xlabel('time')
ylabel('strain_{33}')

saveas(gcf,strcat(testdeg,'_strain33mid.fig'))
saveas(gcf,strcat(testdeg,'_strain33mid.png'))

ylim([Min_strain,Max_strain])

%% Impact edge only
figure(4)
strain33imp=strain33Avg(ind6,:);
strain11imp=strain11Avg(ind6,:);
plot(time.vec,strain33Avg(ind6,:),'--b')
title('Strain_{33} impact edge')
xlabel('time')
ylabel('strain_{33}')

saveas(gcf,strcat(testdeg,'_strain33impact.fig'))
saveas(gcf,strcat(testdeg,'_strain33impact.png'))

