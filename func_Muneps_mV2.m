function [eps_11m,eps_22m,eps_33m] = func_Muneps_mV2(strain11,...
    strain22,strain33,...
    eps_11m0,eps_22m0,eps_33m0,...
    delta_t,tau)

%Author: Andrew Matejunas
%date completed: 2021-10-04
%date verified:

%change log
    %2021-10-27- Did not work for multiple tau. Deleted summations of A11
    %2022-12-05- Performance pass to improve computational speed
                    %Removed C_viscom from the function 


%This script is written to compute epsilon_klm from equation 16 in Mun 2006
%'Numerical Computatioo of Convolution Integral for Linear Viscoelasticity 
%of Asphalt Concrete'. epsilon_klm is computed in the 11, 22, and 33
%directions. 

%Note this is only sufficient in Plane stress conditions for which 
%epsilon_23=epsilon_13=0; For computational efficiency, this script is only
%written to compute eps_klm in the k=l directions. Shear components will be
%in func_Muneps_mshear.m  

%function input arguments
    %strain11- array of normal strains in the x direction at t_(n-1)and t_n 
    %strainn22- array of normal strains in the y direction
    %strain33- array of normal strains in the z direction
        %Note dimensions of the strain array arrays are
        %[number of spatial coordinates, 2 time increments]
    %delta_t- time step in seconds
    %C_Viscom- Array containing the 4th order elasticity tensor for the
        %spring in each Maxwell element in the material model
    %Tau- Vector of time constants for each Maxwell elements
    
%function output arguments
    %eps_11m- epsilon_klm in 11 direction from Mun2006
    %eps_22m- epsilon_klm in 22 direction from Mun2006
    %eps_33m- epsilon_klm in 33 direction from Mun2006
    

   
%% Extract number of spacial coordinates
n_coord=length(strain11(:,1));


%% Calculate eps_11m

%intialize eps_klms for computational efficency
[A11,A22,A33]=deal(zeros(n_coord,length(tau)));

for m=1:length(tau)
    %break into smaller constants
A11(:,m)=strain11(:,1)+exp(-delta_t/tau(m))*(eps_11m0(:,m)-strain11(:,1))...
    +(strain11(:,2)-strain11(:,1))/delta_t...
    *(delta_t-tau(m)*(1-exp(-delta_t/tau(m))));
A22(:,m)=strain22(:,1)+exp(-delta_t/tau(m))*(eps_22m0(:,m)-strain22(:,1))...
    +(strain22(:,2)-strain22(:,1))/delta_t...
    *(delta_t-tau(m)*(1-exp(-delta_t/tau(m))));
A33(:,m)=strain33(:,1)+exp(-delta_t/tau(m))*(eps_33m0(:,m)-strain33(:,1))...
    +(strain33(:,2)-strain33(:,1))/delta_t...
    *(delta_t-tau(m)*(1-exp(-delta_t/tau(m))));

end
%calculate eps_klm

eps_11m=A11;
eps_22m=A22;
eps_33m=A33;

end

