function [strain33,eps11m,eps22m,eps33m]=func_ViscoIBIIstrain33V2(strain11,strain22,...
    time,C_Visco0,C_Viscom,tau)

%Author: Andrew Matejunas
%Date Completed: 2021/10/04
%Date Verified:2021/10/05

%Current Version:2021/10/05
%Change log
    %2021/10/05- V2 corrected major issue with eps_klm calculation. Do not
        %use elasticity tensor in eq(16) from Mun2006. It is invalid!
    %2022-12-05: Performance pass for improved computation time
                % Added outputs for epsIJm so that it only has to be
                    %calculated one time in calculation of the constitutive
                    %model

%This function is written to calculate out of plane stress in plane stress
%impact on the edge of a viscoelastic plate. The material is characterized
%using the generalized maxwell model, and mathmatical formulations used can
%be found in Mun 2006  'Numerical Computatioo of Convolution Integral for
%Linear Viscoelasticity of Asphalt Concrete'

%Note All equations can be found in the attached PDF.
%"Strain33_MunMethod_equations.pdf"

%Function input arguments
    %strain11- array of normal strains in the x direction 
    %strainn22- array of normal strains in the y direction
    %time- vector containing each time point
    %C_Visco0- Array containing the 4th order instantaneous elasticity
        %tensor for the material
    %C_Viscom- Array containing the 4th order elasticity tensor for the
        %spring in each Maxwell element in the material model
    %Tau- Vector of time constants for each Maxwell elements
    
%Function output argument
   %Srain33- Strain in the zz direction
   %eps11m- elemental viscous strain in the 11 direction
   %eps22m- elemental viscous strain in the 22 direction
   %eps33m- elemental viscous strain in the 33 direction

   
%% Define number of Maxwell elements
Num_mx=length(tau);
   
 %%  Reshape strain arrays from 3 dimensional to 2 dimensional arrays
 %dimensions are [# of spacial coordinates,# of time steps]
 n_coord=length(strain11(:,1,1))*length(strain11(1,:,1)); %number of
    %spatial coordinates
 strain11Time=reshape(strain11,[n_coord,length(time)]);
 strain22Time=reshape(strain22,[n_coord,length(time)]);
 
 %% Calculate vector of Delta_t
 
 Delta_t=zeros(length(time),1);
 
 for k=2:length(time)
     Delta_t(k)=time(k)-time(k-1);
 end

%% Calculate strain33 at 2nd time step
%initialize
[strain33Time,D33]=deal(zeros(size(strain11Time)));

%initialize eps_klm
[eps11m,eps22m,eps33m]=deal(zeros([n_coord,Num_mx,length(time)]));
% eps11m(:,:,1)=0;
% eps22m(:,:,1)=0;
% eps33m(:,:,1)=0;

%initialize known constants
[B11m,R11m,B22m,R22m,R33m,D33m]=deal(zeros(n_coord,Num_mx));
%Separate measureable terms into known constants
for m=1:Num_mx
    %strains in 11 direction
    B11m(:,m)=strain11Time(:,1)...
        +exp(-Delta_t(2)/tau(m)).*(eps11m(:,m,1)-strain11Time(:,1))...
        +(strain11Time(:,2)-strain11Time(:,m,1)).*(1-tau(m)/Delta_t(2)...
        .*(1-exp(-Delta_t(2)/tau(m))));
    R11m(:,m)=C_Viscom(3,3,1,1,m)*B11m(:,m);
    
    %22 direction
   B22m(:,m)=strain22Time(:,1)...
        +exp(-Delta_t(2)/tau(m)).*(eps22m(:,m,1)-strain22Time(:,1))...
        +(strain22Time(:,2)-strain22Time(:,1)).*(1-tau(m)/Delta_t(2)...
        .*(1-exp(-Delta_t(2)/tau(m))));
    R22m(:,m)=C_Viscom(3,3,2,2,m)*B22m(:,m);
    %33 direction
    R33m(:,m)=C_Viscom(3,3,3,3,m)*(strain33Time(:,1)...
       +exp(-Delta_t(2)/tau(m))*(eps33m(:,m,1)-strain33Time(:,1))...
       -strain33Time(:,1)...
       *(1-tau(m)/Delta_t(2)*(1-exp(-Delta_t(2)/tau(m)))));
   %Denominator
   D33m(:,m)=C_Viscom(3,3,3,3,m)*(1-tau(m)/Delta_t(2)...
       *(1-exp(-Delta_t(2)/tau(m))));
end


R11=sum(R11m,2);
R22=sum(R22m,2);
R33=sum(R33m,2);
R=R11+R22+R33;
D33(:,2)=C_Visco0(3,3,3,3)-sum(D33m,2);

strain33Time(:,2)=(R-C_Visco0(3,3,1,1)*strain11Time(:,2)...
    -C_Visco0(3,3,2,2)*strain22Time(:,2))./D33(:,2);





%% Calculate strain33 at further time steps 
for n=2:(length(time)-1)
   %define strain inputs for epsklm function
   strain11input=strain11Time(:,(n-1):n);
   strain22input=strain22Time(:,(n-1):n);
   strain33input=strain33Time(:,(n-1):n);
   
   % caclulate eps_{klm}
   [eps11m(:,:,n),eps22m(:,:,n),eps33m(:,:,n)]=func_Muneps_mV2(strain11input,strain22input,...
       strain33input,... strain inputs for t_{n-1} and t_{n}
   eps11m(:,:,n-1),eps22m(:,:,n-1),eps33m(:,:,n-1),... eps_klm(t_{n-1})
   Delta_t(n),tau);
for m=1:Num_mx
    %strains in 11 direction
    B11m(:,m)=strain11Time(:,n)...
        +exp(-Delta_t(n+1)/tau(m)).*(eps11m(:,m,n)-strain11Time(:,n))...
        +(strain11Time(:,n+1)-strain11Time(:,n)).*(1-tau(m)/Delta_t(n+1)...
        .*(1-exp(-Delta_t(n+1)/tau(m))));
    R11m(:,m)=C_Viscom(3,3,1,1,m)*B11m(:,m);
    
    %22 direction
   B22m(:,m)=strain22Time(:,n)...
        +exp(-Delta_t(n+1)/tau(m)).*(eps22m(:,m,n)-strain22Time(:,n))...
        +(strain22Time(:,n+1)-strain22Time(:,n)).*(1-tau(m)/Delta_t(n+1)...
        .*(1-exp(-Delta_t(n+1)/tau(m))));
    R22m(:,m)=C_Viscom(3,3,2,2,m)*B22m(:,m);
    %33 direction
    R33m(:,m)=C_Viscom(3,3,3,3,m)*(strain33Time(:,n)...
       +exp(-Delta_t(n+1)/tau(m))*(eps33m(:,m,n)-strain33Time(:,n))...
       -strain33Time(:,n)*(1-tau(m)/Delta_t(n+1)...
       *(1-exp(-Delta_t(n+1)/tau(m)))));
   %Denominator
   D33m(:,m)=C_Viscom(3,3,3,3,m)*(1-tau(m)/Delta_t(n+1)...
       *(1-exp(-Delta_t(n+1)/tau(m))));
end        
 
R11(:,n+1)=sum(R11m,2);
R22(:,n+1)=sum(R22m,2);
R33(:,n+1)=sum(R33m,2);
R(:,n+1)=R11(:,n+1)+R22(:,n+1)+R33(:,n+1);
D33(:,n+1)=C_Visco0(3,3,3,3)-sum(D33m,2);

  strain33Time(:,n+1)=(R(:,n+1)-C_Visco0(3,3,1,1)*strain11Time(:,n+1)...
    -C_Visco0(3,3,2,2)*strain22Time(:,n+1))./D33(:,n+1);
end




%% Reshape strain33Time to [# X coordinates, # y_coordinates, time 
strain33=reshape(strain33Time,size(strain11));
end
