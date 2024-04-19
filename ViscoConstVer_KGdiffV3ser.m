% This code is written to verify my viscoelastic constitutive modelling
% code for a case when k=/=g

%% initialize
clear all; close all; clc


%% Load FEM output data

filename=uigetfile('.mat','choose SG Data');


%% Load material properties structure
filename2=uigetfile('.mat','Choose Material Properties file');

%% Load stress file
stressfile=uigetfile('*.mat','choose file containing FE stress');
%% Define test designation to autosave figures

prompt='input test designation';
deg=char(inputdlg(prompt));

fprintf('Loading data \n')

load(filename);
load(filename2);
load(stressfile,'stress');

testdeg=deg;
clear deg
%% Run constitutive model function
fprintf('evaluating constitutive model \n')
StressModel=func_ViscoConstitutiveV6(strain.x,strain.y,strain.s,...
    time.vec,...
    MatProps,MatProps.Ki,MatProps.Gi,MatProps.tau);

fprintf('generating plots \n')
%% Coordinate vectors
X_vec=pos.x;
Y_vec=pos.y;

%% Calculate average stresses
StressModel.Avxx=squeeze(mean(StressModel.xx));
StressModel.Avyy=squeeze(mean(StressModel.yy));
StressModel.Avxy=squeeze(mean(StressModel.xy));
StressModel.Avhy=squeeze(mean(StressModel.hy));
StressModel.Avhstrain=squeeze(mean(StressModel.hstrain));

stress.avgX=squeeze(mean(stress.x));
stress.avgY=squeeze(mean(stress.y));
stress.avgXY=squeeze(mean(stress.s));

strain.sAvg=squeeze(mean(strain.s));
strain.xAvg=squeeze(mean(strain.x));
strain.yAvg=squeeze(mean(strain.y));
%% Stress Strain in X-direction constitutive vs FEM
strainxx=squeeze(strain.xAvg);
stressxxFE=squeeze(stress.avgX);

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

figure(1)
% plot StressModel.Avxx results
plot(strainxx(ind1,:),StressModel.Avxx(ind1,:),'b')
hold on
plot(strainxx(ind2,:),StressModel.Avxx(ind2,:),'k')
plot(strainxx(ind3,:),StressModel.Avxx(ind3,:),'r')
plot(strainxx(ind4,:),StressModel.Avxx(ind4,:),'g')
plot(strainxx(ind5,:),StressModel.Avxx(ind5,:),'y')
plot(strainxx(ind6,:),StressModel.Avxx(ind6,:),'c')

%plot ABAQUS stress results
plot(strainxx(ind1,:),stressxxFE(ind1,:),'--b')
plot(strainxx(ind2,:),stressxxFE(ind2,:),'--k','handlevisibility','off')
plot(strainxx(ind3,:),stressxxFE(ind3,:),'--r','handlevisibility','off')
plot(strainxx(ind4,:),stressxxFE(ind4,:),'--g','handlevisibility','off')
plot(strainxx(ind5,:),stressxxFE(ind5,:),'--y','handlevisibility','off')
plot(strainxx(ind6,:),stressxxFE(ind6,:),'--c','handlevisibility','off')
hold off

xlabel('strain')
ylabel('\sigma_{xx} (Pa)')

legend(strcat(lgd1,' Constitutive'),lgd2,lgd3,lgd4,lgd5,lgd6,strcat(lgd1,' FE'),'location','southeast')

saveas(gcf,strcat(testdeg,'_constver'))
saveas(gcf,strcat(testdeg,'_constver.png'))

%% Stress_time in the x-direction
figure(2)
plot(time.vec,StressModel.Avxx(ind1,:),'b')
hold on
plot(time.vec,StressModel.Avxx(ind2,:),'k')
plot(time.vec,StressModel.Avxx(ind3,:),'r')
plot(time.vec,StressModel.Avxx(ind4,:),'g')
plot(time.vec,StressModel.Avxx(ind5,:),'y')
plot(time.vec,StressModel.Avxx(ind6,:),'c')

%plot ABAQUS stress results
plot(time.vec,stressxxFE(ind1,:),'--b')
plot(time.vec,stressxxFE(ind2,:),'--k','handlevisibility','off')
plot(time.vec,stressxxFE(ind3,:),'--r','handlevisibility','off')
plot(time.vec,stressxxFE(ind4,:),'--g','handlevisibility','off')
plot(time.vec,stressxxFE(ind5,:),'--y','handlevisibility','off')
plot(time.vec,stressxxFE(ind6,:),'--c','handlevisibility','off')
hold off

xlabel('time (s)')
ylabel('stress (Pa)')

legend(strcat(lgd1,' const model'),lgd2,lgd3,lgd4,lgd5,lgd6,strcat(lgd1,' FE'),'location','southeast')

saveas(gcf,strcat(testdeg,'_ConstTime.fig'))
saveas(gcf,strcat(testdeg,'_ConstTime.png'))

%% Focus on last curve
% plot StressModel.Avxx results
figure(3)
plot(strainxx(ind6,:),StressModel.Avxx(ind6,:),'c')
hold on
plot(strainxx(ind6,:),stressxxFE(ind6,:),'--c')
hold off

xlabel('strain')
ylabel('\sigma_{xx} (Pa)')

legend(strcat(lgd1,' Constitutive'),strcat(lgd1,' FE'),'location','southeast')

saveas(gcf,strcat(testdeg,'_constImp'))
saveas(gcf,strcat(testdeg,'_constImp.png'))

%% Perform error analysis w.r.t time

error=(StressModel.Avxx(ind6,:)-stressxxFE(ind6,:))./stressxxFE(ind6,:)*100;

figure(4)
plot(time.vec,error)
ylabel('error')
xlabel('time (s)')

saveas(gcf,strcat(testdeg,'error'))
saveas(gcf,strcat(testdeg,'error.png'))


%% Stress Strain in Y-direction constitutive vs FEM
strainyy=squeeze(strain.yAvg);
stressyyFE=squeeze(stress.avgY);

figure(5)
% plot StressModel.Avyy results
plot(strainyy(ind1,:),StressModel.Avyy(ind1,:),'b')
hold on
plot(strainyy(ind2,:),StressModel.Avyy(ind2,:),'k')
plot(strainyy(ind3,:),StressModel.Avyy(ind3,:),'r')
plot(strainyy(ind4,:),StressModel.Avyy(ind4,:),'g')
plot(strainyy(ind5,:),StressModel.Avyy(ind5,:),'y')
plot(strainyy(ind6,:),StressModel.Avyy(ind6,:),'c')

%plot ABAQUS stress results
plot(strainyy(ind1,:),stressyyFE(ind1,:),'--b')
plot(strainyy(ind2,:),stressyyFE(ind2,:),'--k','handlevisibility','off')
plot(strainyy(ind3,:),stressyyFE(ind3,:),'--r','handlevisibility','off')
plot(strainyy(ind4,:),stressyyFE(ind4,:),'--g','handlevisibility','off')
plot(strainyy(ind5,:),stressyyFE(ind5,:),'--y','handlevisibility','off')
plot(strainyy(ind6,:),stressyyFE(ind6,:),'--c','handlevisibility','off')
hold off

xlabel('strain')
ylabel('\sigma_{yy} (Pa)')

legend(strcat(lgd1,' Constitutive'),lgd2,lgd3,lgd4,lgd5,lgd6,strcat(lgd1,' FE'),'location','southeast')

saveas(gcf,strcat(testdeg,'_Yconstver'))
saveas(gcf,strcat(testdeg,'_Yconstver.png'))

%% Stress_time in the yy-direction
figure(6)
plot(time.vec,StressModel.Avyy(ind1,:),'b')
hold on
plot(time.vec,StressModel.Avyy(ind2,:),'k')
plot(time.vec,StressModel.Avyy(ind3,:),'r')
plot(time.vec,StressModel.Avyy(ind4,:),'g')
plot(time.vec,StressModel.Avyy(ind5,:),'y')
plot(time.vec,StressModel.Avyy(ind6,:),'c')

%plot ABAQUS stress results
plot(time.vec,stressyyFE(ind1,:),'--b')
plot(time.vec,stressyyFE(ind2,:),'--k','handlevisibility','off')
plot(time.vec,stressyyFE(ind3,:),'--r','handlevisibility','off')
plot(time.vec,stressyyFE(ind4,:),'--g','handlevisibility','off')
plot(time.vec,stressyyFE(ind5,:),'--y','handlevisibility','off')
plot(time.vec,stressyyFE(ind6,:),'--c','handlevisibility','off')
hold off

xlabel('time (s)')
ylabel('\sigma_{yy} (Pa)')

legend(strcat(lgd1,' const model'),lgd2,lgd3,lgd4,lgd5,lgd6,strcat(lgd1,' FE'),'location','southeast')

saveas(gcf,strcat(testdeg,'_YConstTime.fig'))
saveas(gcf,strcat(testdeg,'_YConstTime.png'))

%% Stress_time in the xy-direction
strainxy=squeeze(strain.sAvg);
stressxyFE=squeeze(stress.avgXY);

figure(7)
plot(time.vec,StressModel.Avxy(ind1,:),'b')
hold on
plot(time.vec,StressModel.Avxy(ind2,:),'k')
plot(time.vec,StressModel.Avxy(ind3,:),'r')
plot(time.vec,StressModel.Avxy(ind4,:),'g')
plot(time.vec,StressModel.Avxy(ind5,:),'y')
plot(time.vec,StressModel.Avxy(ind6,:),'c')

%plot ABAQUS stress results
plot(time.vec,stressxyFE(ind1,:),'--b')
plot(time.vec,stressxyFE(ind2,:),'--k','handlevisibility','off')
plot(time.vec,stressxyFE(ind3,:),'--r','handlevisibility','off')
plot(time.vec,stressxyFE(ind4,:),'--g','handlevisibility','off')
plot(time.vec,stressxyFE(ind5,:),'--y','handlevisibility','off')
plot(time.vec,stressxyFE(ind6,:),'--c','handlevisibility','off')
hold off

xlabel('time (s)')
ylabel('\sigma_{xy} (Pa)')

legend(strcat(lgd1,' const model'),lgd2,lgd3,lgd4,lgd5,lgd6,strcat(lgd1,' FE'),'location','southeast')

saveas(gcf,strcat(testdeg,'_ShearConstTime.fig'))
saveas(gcf,strcat(testdeg,'_ShearConstTime.png'))

%% Compare FE and Const model shear results
figure(8)
% plot StressModel.Avxy results
plot(strainxy(ind1,:),StressModel.Avxy(ind1,:),'b')
hold on
plot(strainxy(ind2,:),StressModel.Avxy(ind2,:),'k')
plot(strainxy(ind3,:),StressModel.Avxy(ind3,:),'r')
plot(strainxy(ind4,:),StressModel.Avxy(ind4,:),'g')
plot(strainxy(ind5,:),StressModel.Avxy(ind5,:),'y')
plot(strainxy(ind6,:),StressModel.Avxy(ind6,:),'c')

%plot ABAQUS stress results
plot(strainxy(ind1,:),stressxyFE(ind1,:),'--b')
plot(strainxy(ind2,:),stressxyFE(ind2,:),'--k','handlevisibility','off')
plot(strainxy(ind3,:),stressxyFE(ind3,:),'--r','handlevisibility','off')
plot(strainxy(ind4,:),stressxyFE(ind4,:),'--g','handlevisibility','off')
plot(strainxy(ind5,:),stressxyFE(ind5,:),'--y','handlevisibility','off')
plot(strainxy(ind6,:),stressxyFE(ind6,:),'--c','handlevisibility','off')
hold off

xlabel('strain')
ylabel('\sigma_{xy} (Pa)')

legend(strcat(lgd1,' Constitutive'),lgd2,lgd3,lgd4,lgd5,lgd6,strcat(lgd1,' FE'),'location','southeast')

saveas(gcf,strcat(testdeg,'_Shearconstver'))
saveas(gcf,strcat(testdeg,'_Shearconstver.png'))

%% Plot Hydrostatic Results
P_FE=(stressxxFE+stressyyFE)/3;

figure(9)

% plot StressModel.Avxy results
plot(StressModel.Avhstrain(ind1,:),StressModel.Avhy(ind1,:),'b')
hold on
plot(StressModel.Avhstrain(ind2,:),StressModel.Avhy(ind2,:),'k')
plot(StressModel.Avhstrain(ind3,:),StressModel.Avhy(ind3,:),'r')
plot(StressModel.Avhstrain(ind4,:),StressModel.Avhy(ind4,:),'g')
plot(StressModel.Avhstrain(ind5,:),StressModel.Avhy(ind5,:),'y')
plot(StressModel.Avhstrain(ind6,:),StressModel.Avhy(ind6,:),'c')

%plot ABAQUS stress results
plot(StressModel.Avhstrain(ind1,:),P_FE(ind1,:),'--b')
plot(StressModel.Avhstrain(ind2,:),P_FE(ind2,:),'--k','handlevisibility','off')
plot(StressModel.Avhstrain(ind3,:),P_FE(ind3,:),'--r','handlevisibility','off')
plot(StressModel.Avhstrain(ind4,:),P_FE(ind4,:),'--g','handlevisibility','off')
plot(StressModel.Avhstrain(ind5,:),P_FE(ind5,:),'--y','handlevisibility','off')
plot(StressModel.Avhstrain(ind6,:),P_FE(ind6,:),'--c','handlevisibility','off')
hold off

xlabel('strain')
ylabel('\sigma_{h} (Pa)')

legend(strcat(lgd1,' Constitutive'),lgd2,lgd3,lgd4,lgd5,lgd6,strcat(lgd1,' FE'),'location','southeast')

saveas(gcf,strcat(testdeg,'_hydro'))
saveas(gcf,strcat(testdeg,'_hydro.png'))

%% hydrostatic Stress and Time
figure(10)

% plot StressModel.Avxy results
plot(time.vec,StressModel.Avhy(ind1,:),'b')
hold on
plot(time.vec,StressModel.Avhy(ind2,:),'k')
plot(time.vec,StressModel.Avhy(ind3,:),'r')
plot(time.vec,StressModel.Avhy(ind4,:),'g')
plot(time.vec,StressModel.Avhy(ind5,:),'y')
plot(time.vec,StressModel.Avhy(ind6,:),'c')

%plot ABAQUS stress results
plot(time.vec,P_FE(ind1,:),'--b')
plot(time.vec,P_FE(ind2,:),'--k','handlevisibility','off')
plot(time.vec,P_FE(ind3,:),'--r','handlevisibility','off')
plot(time.vec,P_FE(ind4,:),'--g','handlevisibility','off')
plot(time.vec,P_FE(ind5,:),'--y','handlevisibility','off')
plot(time.vec,P_FE(ind6,:),'--c','handlevisibility','off')
hold off

xlabel('time (s)')
ylabel('\sigma_{h} (Pa)')

legend(strcat(lgd1,' Constitutive'),lgd2,lgd3,lgd4,lgd5,lgd6,strcat(lgd1,' FE'),'location','southeast')

saveas(gcf,strcat(testdeg,'_hydroTime'))
saveas(gcf,strcat(testdeg,'_hydroTime.png'))



% %% Plot all stress and strain data points
% figure()
% 
% scatter(strain.x(:),StressModel.xx(:),'b')
% hold on
% scatter(strain.x(:),stress.x(:),'r')
% hold off
% 
% xlabel('$ \varepsilon_{xx} $','interpreter','latex')
% ylabel('$ \sigma_{xx}$','interpreter','latex')
% legend('Const Model','FE')
% 
% saveas(gcf,strcat(testdeg,'xSSall'))
% saveas(gcf,strcat(testdeg,'xSSall.png'))
% 
% %%Comparison with Stress Gage
% figure()
% plot(strainxx(end,:),StressModel.Avxx(end,:),'--')
% hold on
% plot(strainxx(end,:),SG(end,:))
% 
% xlabel('strain')
% ylabel('\sigma_{xx} (Pa)')
% legend('Constitutive Model','SG','location','southeast')
% 
% saveas(gcf,strcat(testdeg,'Const_SG'))
% saveas(gcf,strcat(testdeg,'Const_SG.png'))
% 
% 
% fprintf('constitutive model evaluations complete \n')

%% Plot Strain33
fprintf('Plotting Strain33 \n')

figure(11)
plot(time.vec,StressModel.ZZAvstrain(ind1,:),'b')
hold on
plot(time.vec,StressModel.ZZAvstrain(ind2,:),'k')
plot(time.vec,StressModel.ZZAvstrain(ind3,:),'r')
plot(time.vec,StressModel.ZZAvstrain(ind4,:),'g')
plot(time.vec,StressModel.ZZAvstrain(ind5,:),'y')
plot(time.vec,StressModel.ZZAvstrain(ind6,:),'c')
hold off
xlabel('Time')
ylabel('Strain_{33}')
legend(strcat(lgd1,' Constitutive'),lgd2,lgd3,lgd4,lgd5,lgd6,'location','southeast')
saveas(gcf,strcat(testdeg,'_strain33.fig'))
saveas(gcf,strcat(testdeg,'_strain33.png'))

%% Save the constitutive model stress information
fprintf('Saving Data')
save(strcat(testdeg,'_constData.mat'),'time','strain','stress','StressModel',...
    'strainxx','stressxxFE','strainyy','stressyyFE',...
    'P_FE','strainxy','stressxyFE');

fprintf('Constitutive Model Verification Complete \n')
