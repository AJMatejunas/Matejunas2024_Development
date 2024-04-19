function phi_K=func_ViscoBulkTauCostSG(xxSG,knownParam,constParam,...
    strain,time,CondOpts)
%This code was created to identify the bulk moduli of a viscoelastic
    %material using a stress gauge approach given previously identified
    %shear moudli and associated time constants.

%Author: Andrew Matejunas (andrew.matejunas@gmail.com)
%Date Verified:

%Change Log:
    %2023/01/16- Adapted the code from func_ViscoKGcostSGV5.
                %Removed shear and time constant calculations from the
                    %cost function.
                %Removed ability to deal with unknown long-term parameters
   %2023/03/15: Added back the ability to calculate time constant and
                    %changed name to func_ViscoBulkTauCostSG from 
                    %func_ViscoBulkCostSG
        
    
%% Function input arguments
    %xxSG- Average stress at each x coordinate calculated using the normal
        %stress gage equation
    %knownParam- A row vector of known parameters in the following order
        %Kinf- long term bulk modulus
        %Ginf- long term shear modulus
        %nu- Poisson's ration 
        %constnu- true if nu is constant
        %Gi- Previously identified shear moduli 
        %tau- Previously identified time constants associated with the
            %shear moduli

     %constParam- vector
        %[K1,tau1,K2,tau2,...,Kn,taun]; if Kinf&Ginf are known,
        %CondOpts- Data Conditioning options. Spatial Downsampling and data
        %sensoring are perfomed before being passed into the function.
        %Temporal downsampling is performed here. Relevant fields
            %TempDS- indicates whether or not temporal downsampling will be
                %performed with options:
                    %true- Temporal downsampling performed
                    %false- No temporal downsampling
            %Tds- Temporal downsampling factor
            
%% Function output arguments
    %phi- Weighted cost function

%% Convert Known and Unknown Parameters into inputs for constitutive model

%% QS moduli are assumed known
%known Parameters
MatProps.Kinf=knownParam(1);
MatProps.Ginf=knownParam(2);

if length(knownParam)==5
    MatProps.Gi=knownParam(5);
elseif length(knownParam)>5
    MatProps.Gi=knownParam(5:end);
end

constParam=reshape(constParam,[numel(constParam)/2,2]);
   
MatProps.Ki=constParam(:,1);
MatProps.tau=constParam(:,2);
MatProps.K0=MatProps.Kinf+sum(MatProps.Ki);
MatProps.G0=MatProps.Ginf+sum(MatProps.Gi);
   

if knownParam(4)==1
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
switch CondOpts.TempDS
    case false %Use every time data point.
        %Avxy=StressModel.Avxy;
        Avxx=StressModel.Avxx;
    
    case true %Downsample temporally
        [Avxx,~,xxSG,~]=func_TempDS(StressModel.Avxx,...
            StressModel.Avxy,xxSG,ShearSG,CondOpts.Tds);
end


%% Calculate cost function for K
phi_K=(xxSG-Avxx).^2;
phi_K=sum(phi_K,'all');

end