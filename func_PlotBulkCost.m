function [phi_ctau,phi_cG,phi3D] = func_PlotBulkCost(K_vec,G_vec,tau_vec, ...
    G1,tau1,RefPar,...
    strain, time,Full_SG, ...
    CondOpts,SaveDir,Desig)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up input parameters for cost function file
knownParam=zeros(1,4);
knownParam(1)=RefPar.Kinf;
knownParam(2)=RefPar.Ginf;

timevec=time.vec;
%% Transpose G and Tau vectors (to help with parfor loop integration)
G_vec=G_vec';
tau_vec=tau_vec';
%tau_z=reshape(tau_vec,[1,1,length(tau_vec)]);
%% Create parameter vectors for parfor loop integration
Kin_cTau=K_vec.*ones(size(G_vec));
Kvec_ctau=reshape(Kin_cTau,[numel(Kin_cTau),1]);

Kin_cG=K_vec.*ones(size(tau_vec));
Kvec_cG=reshape(Kin_cG,[numel(Kin_cG),1]);

%create G input vector
Gin_cTau=G_vec.*ones(size(K_vec));
Gvec_ctau=reshape(Gin_cTau,[numel(Gin_cTau),1]);

%constant G
tau_in=tau_vec.*ones(size(K_vec));
tau_inVec=reshape(tau_in,[numel(tau_in),1]);

%% Remove stress gage fields from structure
xxSG=Full_SG.x;
ShearSG=Full_SG.s;
%% preallocate memory for phi_ctau
phi_ctau_vec=zeros(size(Kvec_ctau));
constParam=ones(length(Kvec_ctau),3)*tau1;
constParam(:,1)=Kvec_ctau;
constParam(:,2)=Gvec_ctau;

%% Calculate constant tau cost function
fprintf('Calculating cost function for constant tau \n')
parfor k=1:length(Kvec_ctau)
    tempParam=squeeze(constParam(k,:));
     phi_ctau_vec(k)=func_ViscoKGcostSGV5(xxSG,ShearSG,knownParam, ...
         tempParam,strain,timevec,CondOpts);
end
phi_ctau=reshape(phi_ctau_vec,size(Kin_cTau));
phi_ctauLog=log(phi_ctau);
%% Plot constant tau cost function
figure('units','centimeters','InnerPosition',[10,10,18,9])
subplot(1,2,1)
contourf(K_vec*10^-9,G_vec*10^-9,phi_ctau)
cx=colorbar;
cx.Label.String='\phi';
colormap('cool');
title('Constant \tau')
xlabel('K_1 (GPa)')
ylabel('G_1 (GPa)')

subplot(1,2,2)
contourf(K_vec*10^-9,G_vec*10^-9,phi_ctauLog)
cx=colorbar;
cx.Label.String='log(\phi)';
set(gca,'ColorScale','log');
colormap('cool');
title('Constant \tau')
xlabel('K_1 (GPa)')
ylabel('G_1 (GPa)')
%% SaveFigure
SaveName=strcat(SaveDir,'/',Desig,'_KcostConstTau');
saveas(gcf,SaveName,'fig')
saveas(gcf,SaveName,'svg')

%% Calculate Constant G cost function
constParam=ones(length(Kvec_ctau),3)*G1;
constParam(:,1)=Kvec_ctau;
constParam(:,3)=tau_inVec;
phi_cG_vec=zeros(size(Kvec_cG));

fprintf('Calculating cost function for constant G \n')
parfor k=1:length(Kvec_cG)
    tempParam=squeeze(constParam(k,:));
     phi_cG_vec(k)=func_ViscoKGcostSGV5(xxSG,ShearSG,knownParam, ...
         tempParam,strain,timevec,CondOpts);
end
phi_cG=reshape(phi_cG_vec,size(Kin_cG));
phi_cGLog=log(phi_cG);

%% Plot constant G cost function
figure('units','centimeters','InnerPosition',[10,10,18,9])
subplot(1,2,1)
contourf(K_vec*10^-9,tau_vec*10^6,phi_cG)
cx=colorbar;
cx.Label.String='\phi';
colormap('cool');
title('Constant G')
xlabel('K_1 (GPa)')
ylabel('\tau_1 (\mu s)')

subplot(1,2,2)
contourf(K_vec*10^-9,tau_vec*10^6,phi_cGLog)
cx=colorbar;
cx.Label.String='log(\phi)';
set(gca,'ColorScale','log');
colormap('cool');
title('Constant G')
xlabel('K_1 (GPa)')
ylabel('\tau_1 (\mu s)')

SaveName=strcat(SaveDir,'/',Desig,'_KcostConstG');
saveas(gcf,SaveName,'fig')
saveas(gcf,SaveName,'svg')
%3D plot not yet ready to implement
phi3D=false;

end