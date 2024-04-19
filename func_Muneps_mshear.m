function [eps_12m,eps_13m,eps_23m] = func_Muneps_mshear(strain12,...
     eps_12m0,eps_13m0,eps_23m0,...
    delta_t,tau)

%Author: Andrew Matejunas
%date completed:
%date verified:

%change log
    %2021-10-27: Did not work for multiple tau, removed unnecessary
      %summations in A12, A13, and A23
    %2022-12-05: Performance pass for computational speed
               % Fixed performance analysis comments
               % converted strain 13 and strain23 to a scalar because it
               % will always be zero
               %Removed C_viscom from inputs

%This script is written to compute epsilon_klm from equation 16 in Mun 2006
%'Numerical Computation of Convolution Integral for Linear Viscoelasticity 
%of Asphalt Concrete'. epsilon_klm is computed in the 12, 13, and 23
%directions. 

%Note this is only sufficient in Plane stress conditions for which 
%epsilon_23=epsilon_13=0; For computational efficiency, this script is only
%written to compute eps_klm in shear k=/=l
%function input arguments
    %strain12- array of shear strains in the xy direction at t_(n-1)and t_n 
    %delta_t- time step in seconds
    %Tau- Vector of time constants for each Maxwell elements
    
%function output arguments
    %eps_11m- epsilon_klm in 11 direction from Mun2006
    %eps_13m- epsilon_klm in 13 direction from Mun2006
    %eps_23m- epsilon_klm in 23 direction from Mun2006
    

   
%% Extract number of spacial coordinates
n_coord=length(strain12(:,1));

%% Plane stress conditions define shear strain in 13 and 23 directions
%[strain13,strain23]=deal(zeros(size(strain12)));
strain13=0;
strain23=0;
%% Calculate eps_12m

%intialize eps_klms for computational efficency
[A12,A13,A23]=deal(zeros(n_coord,length(tau)));
for m=1:length(tau)
    %break into smaller constants
A12(:,m)=strain12(:,1)+exp(-delta_t/tau(m))*(eps_12m0(:,m)-strain12(:,1))...
    +(strain12(:,2)-strain12(:,1))/delta_t...
    *(delta_t-tau(m)*(1-exp(-delta_t/tau(m))));
% A13(:,m)=strain13(:,1)+exp(-delta_t/tau(m))*(eps_13m0(:,m)-strain13(:,1))...
%     +(strain13(:,2)-strain13(:,1))/delta_t...
%     *(delta_t-tau(m)*(1-exp(-delta_t/tau(m))));
% A23(:,m)=strain23(:,1)+exp(-delta_t/tau(m))*(eps_23m0(:,m)-strain23(:,1))...
%     +(strain23(:,2)-strain23(:,1))/delta_t...
%     *(delta_t-tau(m)*(1-exp(-delta_t/tau(m))));
A13(:,m)=strain13+exp(-delta_t/tau(m))*(eps_13m0(:,m)-strain13)...
    +(strain13-strain13)/delta_t...
    *(delta_t-tau(m)*(1-exp(-delta_t/tau(m))));
A23(:,m)=strain23+exp(-delta_t/tau(m))*(eps_23m0(:,m)-strain23)...
    +(strain23-strain23)/delta_t...
    *(delta_t-tau(m)*(1-exp(-delta_t/tau(m))));
%calculate eps_klm

end
eps_12m=A12;
eps_13m=A13;
eps_23m=A23;

end

