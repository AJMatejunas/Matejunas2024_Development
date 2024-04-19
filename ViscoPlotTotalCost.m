% This code is written to plot the cost function for simultaneous
% identification of K1, G1, and tau1 in a Prony series formulation of the
% Generalized Maxwell viscoelasticity model

%Author: Andrew Matejunas

%Date Verified:

%Change Log:
    %2021-11-05: Constitutive model algorithm changed to
        %func_ViscoConstitutiveV6 using the Mun2006 stress reconstruction
        %algorithm
    %2022-02-08: Contour plots changed to imageSC for better detail
        %      : Hardcoded K and E exact are now taken directly from
        %          properties file
%NOTE:
    % Code derived from Visco_Shear_SGcost_evaluationV3.m
    

%% Initialize
clear all; close all; clc

%% Find Stress Gage and Parameter data files
SG_filename=uigetfile('.mat','Choose File Containing Stress Gage Data');
param_file=uigetfile('.mat','Chooe File Containing Material Parameters');

%% define test designation
prompt='input test designation';
desig=char(inputdlg(prompt));

%% load files and save test designation
fprintf('Loading data \n')
load(SG_filename);
param=load(param_file);
testdeg=desig;

clear desig

%% define exact material property
MatProps.nu_exact=0.26; 
MatProps.E_exact=2.21E9; %Pa
MatProps.tau_exact=param.MatProps.tau;


MatProps.K_exact=param.MatProps.Ki; %pa
MatProps.G_exact=param.MatProps.Gi;

%% Define guess parameters
fprintf('Defining Search Parameters \n')

E_guess=[.5,.6,.7,.8,.85,.9,.95,1,1.05,1.1,1.15,1.2,1.25,1.3,1.4,1.5]*MatProps.E_exact; %Elastic Modulus
G_guess=[.5,.6,.7,.8,.85,.9,.95,1,1.05,1.1,1.15,1.2,1.25,1.3,1.4,1.5]*MatProps.G_exact;%Shear Modulus
K_guess=[.5,.6,.7,.8,.85,.9,.95,1,1.05,1.1,1.15,1.2,1.25,1.3,1.4,1.5]*MatProps.K_exact;

tau_guess=[.5,.6,.7,.8,.85,.9,.95,1,1.05,1.1,1.15,1.2,1.25,1.3,1.4,1.5]*MatProps.tau_exact;

%save the gueses
save(strcat(testdeg,'_guesses.mat'),'E_guess','tau_guess');

%% Generate vector of known Parameters
MatProps.Kinf=param.MatProps.Kinf;
MatProps.Ginf=param.MatProps.Ginf;
knownParam=[MatProps.Kinf,MatProps.Ginf,MatProps.nu_exact,true];

%% Calculate average strain fields for plotting purposes
AVstrainxx=squeeze(mean(strain.x));
AVstrainyy=squeeze(mean(strain.y));
AVstrainxy=squeeze(mean(strain.s));

%% Generate Cost function with heat maps for varying K and G at fixed Tau
fprintf('Evaluating Cost Function for varying K and G \n')
tau_fixed=MatProps.tau_exact;

phi_KG=zeros(length(E_guess),length(E_guess));

for m=1:length(E_guess)
    parfor n=1:length(E_guess)
        %MatProps.Ei=E_guess(m);
        %MatProps.tau=tau_guess(n);
%% create parameter vector for with input parameters for cost function
    %Note: These are the constParam input vector in func_ViscoShearCostSG
        %for proper order of the components 
    constParam=[K_guess(m),G_guess(n),tau_fixed];

    %THESE LINES NO LONGER NEEDED AS THE MODEL STRESS IS CALCULATED WITHIN
    %THE COST FUNCTION THESE LINES OF CODE WILL BE DELETED ONCE I VERIFY
    %THAT THIS NEW IMPLEMENTATION WORKS
    %         StressModelCF=func_ViscoConstitutiveV4(strain.x,strain.y,strain.s,...
     %       time.vec,MatProps,K_guess(m),G_guess(m),tau_guess(n));
   
  %save(strcat('Modelxx_',num2str(m),'_',num2str(n),'.mat'),'Modelxx')
 
    phi_KG(m,n)=func_ViscoKGcostSGV2(SG,Shear_SG,...
        knownParam,...
        constParam,...
        strain,...
        time.vec);
    %save(strcat('phi_',num2str(m),'_',num2str(n),'_.mat'),'PHI_inst');
    
    %PHI SUM IS NOW CALCULATED WTITHIN THE COST FUNCTION PROGRAM. THIS CAN
    %BE DELETED ONCE I VERIFY THIS CODE WORKS
    %phi_sum(m,n)=sum(PHI_inst,'all');
    
    
    end
end


%% Plot cost function for varying K and G
figure
imagesc(K_guess/MatProps.K_exact,G_guess/MatProps.G_exact,phi_KG')
xlabel('K_{1}/K_{1,ref}')
ylabel('G_{1}/G_{1,ref}')
cx=colorbar;
cx.Label.String='\phi';

saveas(gcf,strcat(testdeg,'KGSGCF_map'))
saveas(gcf,strcat(testdeg,'KGSGCF_map.png'))

%% Plot Logarithm of the contour plot 
logphi_KG=log10(phi_KG);
figure
imagesc(K_guess/MatProps.K_exact,G_guess/MatProps.G_exact,logphi_KG')
xlabel('K_{1}/K_{1,ref}')
ylabel('G_{1}/G_{1,ref}')
cx=colorbar;
cx.Label.String='log(\phi)';


saveas(gcf,strcat(testdeg,'KGlogSGCF_map'))
saveas(gcf,strcat(testdeg,'KGlogSGCF_map.png'))

%% Cost function for fixed G
fprintf('Generating Cost Functions for Fixed G \n')

phi_Ktau=zeros(length(E_guess),length(tau_guess));

G_fixed=MatProps.G_exact;
for m=1:length(E_guess)
    parfor n=1:length(tau_guess)
        %MatProps.Ei=E_guess(m);
        %MatProps.tau=tau_guess(n);
%% create parameter vector for with input parameters for cost function
    %Note: These are the constParam input vector in func_ViscoShearCostSG
        %for proper order of the components 
    constParam=[K_guess(m),G_fixed,tau_guess(n)];

    %THESE LINES NO LONGER NEEDED AS THE MODEL STRESS IS CALCULATED WITHIN
    %THE COST FUNCTION THESE LINES OF CODE WILL BE DELETED ONCE I VERIFY
    %THAT THIS NEW IMPLEMENTATION WORKS
    %         StressModelCF=func_ViscoConstitutiveV4(strain.x,strain.y,strain.s,...
     %       time.vec,MatProps,K_guess(m),G_guess(m),tau_guess(n));
   
  %save(strcat('Modelxx_',num2str(m),'_',num2str(n),'.mat'),'Modelxx')
 
   phi_Ktau(m,n)=func_ViscoKGcostSGV2(SG,Shear_SG,...
        knownParam,...
        constParam,...
        strain,...
        time.vec);
    %save(strcat('phi_',num2str(m),'_',num2str(n),'_.mat'),'PHI_inst');
  
  
    
    end
end

%Plot for fixed G
figure
imagesc(K_guess/MatProps.K_exact,tau_guess/MatProps.tau_exact,phi_Ktau')
xlabel('K_{1}/K_{ref}')
ylabel('\tau_{1}/\tau_{1,ref} (s)')
cx=colorbar;
cx.Label.String='\phi';

saveas(gcf,strcat(testdeg,'KtSGCF_map'))
saveas(gcf,strcat(testdeg,'KtSGCF_map.png'))

% Plot Logarithm of the contour plot 
logphi_Ktau=log10(phi_Ktau);
figure
imagesc(K_guess/MatProps.K_exact,tau_guess/MatProps.tau_exact,logphi_Ktau')
xlabel('K_{1}/K_{ref}')
ylabel('\tau_{1}/\tau_{1,ref} (s)')
zlabel('log(\phi)')
cx=colorbar;
cx.Label.String='log(\phi)';

saveas(gcf,strcat(testdeg,'KtlogSGCF_map'))
saveas(gcf,strcat(testdeg,'KtlogSGCF_map.png'))
 
%% Cost function for fixed K
fprintf('Evaluating Cost Function for Fixed K \n')

phi_Gtau=zeros(length(E_guess),length(tau_guess));

K_fixed=MatProps.K_exact;
for m=1:length(E_guess)
    parfor n=1:length(tau_guess)
        %MatProps.Ei=E_guess(m);
        %MatProps.tau=tau_guess(n);
%% create parameter vector for with input parameters for cost function
    %Note: These are the constParam input vector in func_ViscoShearCostSG
        %for proper order of the components 
    constParam=[K_fixed,G_guess(m),tau_guess(n)];

    %THESE LINES NO LONGER NEEDED AS THE MODEL STRESS IS CALCULATED WITHIN
    %THE COST FUNCTION THESE LINES OF CODE WILL BE DELETED ONCE I VERIFY
    %THAT THIS NEW IMPLEMENTATION WORKS
    %         StressModelCF=func_ViscoConstitutiveV4(strain.x,strain.y,strain.s,...
     %       time.vec,MatProps,K_guess(m),G_guess(m),tau_guess(n));
   
  %save(strcat('Modelxx_',num2str(m),'_',num2str(n),'.mat'),'Modelxx')
 
   phi_Gtau(m,n)=func_ViscoKGcostSGV2(SG,Shear_SG,...
        knownParam,...
        constParam,...
        strain,...
        time.vec);
    %save(strcat('phi_',num2str(m),'_',num2str(n),'_.mat'),'PHI_inst');
  
  
    
    end
end

%Plot for fixed K
figure
imagesc(G_guess/MatProps.G_exact,tau_guess/MatProps.tau_exact,phi_Gtau')
xlabel('G_{1}/G_{ref}')
ylabel('\tau_{1}/\tau_{1,ref} (s)')
colorbar

saveas(gcf,strcat(testdeg,'GtSGCF_map'))
saveas(gcf,strcat(testdeg,'GtSGCF_map.png'))

% Plot Logarithm of the contour plot 
logphi_Gtau=log10(phi_Gtau);
figure
contourf(G_guess/MatProps.G_exact,tau_guess/MatProps.tau_exact,logphi_Gtau')
xlabel('G_{1}/G_{ref}')
ylabel('\tau_{1}/\tau_{1,ref} (s)')
zlabel('log(\phi)')
cx=colorbar;
cx.Label.String='log(\phi)';

saveas(gcf,strcat(testdeg,'GtlogSGCF_map'))
saveas(gcf,strcat(testdeg,'GtlogSGCF_map.png'))
%% Save Cost Function
save(strcat(testdeg,'_KGphiSG.mat'),'phi_KG','phi_Gtau','phi_Ktau');
