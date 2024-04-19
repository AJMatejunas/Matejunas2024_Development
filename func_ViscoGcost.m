function phi = func_ViscoGcost(ShearGage,knownParam,constParam,strain,time)
%This code is intended to generate the cost function used for
%identification of the shear modulus and time constant for a viscoelastic
%constitutive model. This is in an attempt to improve identification by
%identifying shear first, because shear stress is fully decoupled from
%normal strains.

%Author: Andrew Matejunas
%Date Verified:
    
    %Change log:
    

%Function input arguments
    %ShearGage- Average stress at each x coordinate calculated using the
        %ShearStress gage equation
    %knownParam- A row vector of known parameters in the following order
        %Kinf- long term bulk modulus
        %Ginf- long term shear modulus
        %nu- Poisson's ration (
        %constnu- true if nu is constant 
    %constParam- vector 
        %[K1,G1,tau1,...,Kn,Gn,taun]; if Kinf&Ginf are known,
        %[Kinf,Ginf,K1,G1,tau1,...Kn,Gn,taun]; if Kinf&Ginf are unknown
    %strain- structure containing strains with fields
       %s- shear strains. Only shear strains are needed
       %x- normal strains in x direction
       %y- Normal strains in y direction
    %time- vector containing time data points
    
%Function output arguments        
    %Phi- Cost function as a function of X coordinate and time


%% Convert known parameters and unknowkn parameters into proper inputs for
    %constitutive model
      
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

%% Calculate average stresses in the x direction for the cost function
StressModel=func_ViscoConstitutiveV6(zeroes(size(strain.s)),...
    zeros(size(strain.s)),strain.s,...
    time,MatProps,0,0,0); %note the zeros are here because I am not working
                          %with exact values here

Avxy=StressModel.Avxy;

%% Calculate cost function
phi=(ShearGage-Avxy).^2;
phi=sum(phi,'all');

end

