function phi=func_ViscoKGcostSGV2(xxSG,ShearSG,knownParam,constParam,strain,time)
%This code is intended to generate a weighted cost function used for
%identification of the pairs for a vscoelastic constituitive model

%Author: Andrew Matejunas (andrew.matejunas@gmail.com)
%Date Verified:

%Change Log:

%% Function input arguments
    %xxSG- Average stress at each x coordinate calculated using the normal
        %stress gage equation
    %ShearSG- Average stress at each coordinated calculated with the shear
        %stress gage equation
    %knownParam- A row vector of known parameters in the following order
        %Kinf- long term bulk modulus
        %Ginf- long term shear modulus
        %nu- Poisson's ration 
        %constnu- true if nu is constant  
     %constParam- vector 
        %[K1,G1,tau1,...,Kn,Gn,taun]; if Kinf&Ginf are known,
        %[Kinf,Ginf,K1,G1,tau1,...Kn,Gn,taun]; if Kinf&Ginf are unknown
%% Function output arguments
    %phi- Weighted cost function

%% Convert Known and Unknown Parameters into inputs for constitutive model

if knownParam(1)~=0
     %% QS moduli are assumed known
       %known Parameters
        MatProps.Kinf=knownParam(1);
        MatProps.Ginf=knownParam(2);

        %unknown Parameters
        constParam=reshape(constParam,[numel(constParam)/3,3]); %a little
        %computationally expensive
    
    
 

        MatProps.Ki=constParam(:,1);
        MatProps.Gi=constParam(:,2);
        MatProps.tau=constParam(:,3);
        MatProps.K0=MatProps.Kinf+sum(MatProps.Ki);
        MatProps.G0=MatProps.Ginf+sum(MatProps.Gi);
    elseif knownParam(1)==0
        MatProps.Kinf=constParam(1);
        MatProps.Ginf=constParam(2);
        
       constParam=constParam(3:end);
       constParam=reshape(constParam,[numel(constParam)/3,3]); %a little
        %computationally expensive
        
        MatProps.Ki=constParam(:,1);
        MatProps.Gi=constParam(:,2);
        MatProps.tau=constParam(:,3);
        MatProps.K0=MatProps.Kinf+sum(MatProps.Ki);
        MatProps.G0=MatProps.Ginf+sum(MatProps.Gi);
end

if knownParam(end)==1
    MatProps.nu=knownParam(3);
else
    MatProps.nu=0;
end

% force a calculation with K and G
MatProps.E0=0;
MatProps.Ei=0;
MatProps.Einf=0;

%% Calculate average stresses for the cost function

StressModel=func_ViscoConstitutiveV6(strain.x,strain.y,strain.s,...
    time,MatProps,0,0,0); %note the zeros are here because I am not working
                          %with exact values here

Avxy=StressModel.Avxy;
Avxx=StressModel.Avxx;

%% Calculate cost function for K
phi_K=(xxSG-Avxx).^2;

%% Calculate cost function for G
phi_G=(ShearSG-Avxy).^2;

%% Calculate total cost function
phi=(phi_K/2+phi_G/2);
phi=sum(phi,'all');

end