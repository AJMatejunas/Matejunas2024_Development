%This code is written to generate evaluation plots for the cost function
%for a viscoelastic material in an IBII test.

% NOTE

%Author: Andrew Matejunas
%Finalized: 1/11/2021
% Rearanged for efficiency 4/5/2021


%% Initialize
clear all; close all; clc

%% Load stress gauge and strain data
SG_filename=uigetfile('.mat','Choose File Containing Stress Gage Data');

%% Load Material Parameters File (these are the real parameters)
param_file=uigetfile('.mat','Chooe File Containing Material Parameters');


%% Get Test designation
prompt='input test designation';
desig=char(inputdlg(prompt));

%% load files
fprintf('loading files \n')
load(param_file);
load(SG_filename);

testdeg=desig;
clear desig

%define exact material property
MatProps.nu_exact=0.26; 
MatProps.E_exact=2.21E9; %Pa
MatProps.tau_exact=input('time constant');

MatProps.K_exact=MatProps.E_exact/(3*(1-2*MatProps.nu_exact)); %pa
MatProps.G_exact=MatProps.E_exact/(2*(1+MatProps.nu_exact));



%% Calculate Stress Using exact parameters
fprintf('evaluating model with reference parameters \n')
 StressModelCFexact=func_ViscoConstitutiveV4(strain.x,strain.y,strain.s,...
        time.vec,MatProps,MatProps.K_exact,MatProps.G_exact,MatProps.tau_exact);

%% Define guess parameters
E_guess=[.5,.6,.7,.8,.85,.9,.95,1,1.05,1.1,1.15,1.2,1.25,1.3,1.4,1.5]*MatProps.E_exact; %Elastic Modulus
G_guess=E_guess/(2*(1+MatProps.nu_exact)); %Shear Modulus
K_guess=E_guess/(3*(1-2*MatProps.nu_exact));

tau_guess=[.5,.6,.7,.8,.85,.9,.95,1,1.05,1.1,1.15,1.2,1.25,1.3,1.4,1.5]*MatProps.tau_exact;

save(strcat(testdeg,'_guesses.mat'),'E_guess','tau_guess');

%% Generate vector of known Parameters
KnownParam=[MatProps.Kinf,MatProps.Ginf,MatProps.nu,true];

%% Calculate average strain fields for plotting purposes
AVstrainxx=squeeze(mean(strain.x));
AVstrainyy=squeeze(mean(strain.y));
AVstrainxy=squeeze(mean(strain.s));

%% Run Constitutive model and Gnerate Cost Function and Generate heat maps
    %of the cost function parfor varying E
    

% phi=zeros(length(SG(:,1)),length(SG(1,:)),length(E_guess));
% parfor k=1:length(E_guess)
%     MatProps.Ei=E_guess(k);
%     StressModel=func_ViscoConstitutiveV4(strain.x,strain.y,strain.s,...
%         time.vec,MatProps);
% %     
%     phi(:,:,k)=func_ViscoCostSG(SG,StressModel);
%     
% figure
% contourf(X_vec,time.vec,phi(:,:,k)')
% title(strcat(testdeg,'-E=',num2str(E_guess(k)/10^9),'-SGCF'))
% xlabel('X coordinate')
% ylabel('time')
% colorbar
% 
% saveas(gcf,strcat(testdeg,'_E_',num2str(k),'_SGCF'))
% saveas(gcf,strcat(testdeg,'_E_',num2str(k),'_SGCF.png'))
% 

%     
% end    

%% Calculate cost function and generate heat maps parfor varying time constant

% parfor k=1:length(tau_guess)
% MatProps.taui=tau_guess(k);
%     StressModel=func_ViscoConstitutiveV4(strain.x,strain.y,strain.s,...
%         time.vec,MatProps);
%     
%     phi(:,:,k)=func_ViscoCostSG(SG,StressModel);
%     
% figure(length(E_guess)+k)
% contourf(X_vec,time.vec,phi(:,:,k)')
% title(strcat(testdeg,'-tau=',num2str(tau_guess(k)),'-SGCF'))
% xlabel('X coordinate')
% ylabel('time')
% colorbar
% 
% saveas(gcf,strcat(testdeg,'_tau_',num2str(k),'_SGCF'))
% saveas(gcf,strcat(testdeg,'_tau_',num2str(k),'_SGCF.png'))
% end
% 
%% Generate heat map parfor total cost function of all initial guess parameter
fprintf('calculating cost functions \n')

phi_sum=zeros(length(E_guess),length(tau_guess));

for m=1:length(E_guess)
    parfor n=1:length(tau_guess)
        %MatProps.Ei=E_guess(m);
        %MatProps.tau=tau_guess(n);
%% create parameter vector for with input parameters for cost function
    %Note: These are the constParam input vector in func_ViscoShearCostSG
        %for proper order of the components 
    constParam=[K_guess(m),G_guess(m),tau_guess(n)];

    %THESE LINES NO LONGER NEEDED AS THE MODEL STRESS IS CALCULATED WITHIN
    %THE COST FUNCTION THESE LINES OF CODE WILL BE DELETED ONCE I VERIFY
    %THAT THIS NEW IMPLEMENTATION WORKS
    %         StressModelCF=func_ViscoConstitutiveV4(strain.x,strain.y,strain.s,...
     %       time.vec,MatProps,K_guess(m),G_guess(m),tau_guess(n));
   
  %save(strcat('Modelxx_',num2str(m),'_',num2str(n),'.mat'),'Modelxx')
 
    PHI_inst=func_ViscoShearCostSG(Shear_SG,...
        KnownParam,...
        constParam,...
        strain,...
        time.vec);
    %save(strcat('phi_',num2str(m),'_',num2str(n),'_.mat'),'PHI_inst');
    
    %PHI SUM IS NOW CALCULATED WTITHIN THE COST FUNCTION PROGRAM. THIS CAN
    %BE DELETED ONCE I VERIFY THIS CODE WORKS
    %phi_sum(m,n)=sum(PHI_inst,'all');
    phi(m,n)=PHI_inst;
    
    end
 end


%% Save Cost Function

fprintf('saving cost function \n')
save(strcat(testdeg,'_phiSG.mat'),'phi');

%% Plot
figure
contourf(G_guess/MatProps.G_exact,tau_guess/MatProps.tau_exact,phi')
xlabel('E_{guess}/E_{exact}')
ylabel('\tau_{guess}/\tau_{exact} (s)')
colorbar

saveas(gcf,strcat(testdeg,'SGCF_map'))
saveas(gcf,strcat(testdeg,'SGCF_map.png'))

%% Plot Logarithm of the contour plot 
logphi=log10(phi);
figure
contourf(E_guess/MatProps.E_exact,tau_guess/MatProps.tau_exact,logphi')
xlabel('G/G_{ref}')
ylabel('\tau/\tau_{ref} (s)')
zlabel('log(\phi)')
colorbar

saveas(gcf,strcat(testdeg,'logSGCF_map'))
saveas(gcf,strcat(testdeg,'logSGCF_map.png'))